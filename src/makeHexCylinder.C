#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "box.h"
#include "l-bfgs.h"
#include "clusterv.h"
#include "spline.h"
#include "util.h"
#include <sys/time.h>

int glob_n_outer = 5;
int glob_n_anneal = 50;
int glob_n_opt = 50;

double del_factor = 0.5;

int do_print = 0;
double dens_factor = 1.0;
double v_penalty = 1000;
int use_box = 0;
//#define USE_BOX
//#define POT_SQR
//#define ELIM
//#define CLEAN
int nfixed = 0;
double r_prot_use = 16.0;
double r_inner = 10;
double particle_inner_r = 2.0;
double last_gmag = 1e10;
double cutoff = 10;
double dsmooth = cutoff/10;
//double point_density = 0.3; //1.0; //0.4;
double PeriodicLengthA = 100; 
double PeriodicLengthB = 100; 
double curvature = 0;
double R = 0;
double **pcen;//[pr][2] = { PeriodicLength/2, PeriodicLength/2 };
int *saved_sense;
int nedge_sense = 0;
int nSenseSpace = 0;
int *nLPts;
int doProtein = 1;
int doShape = BOX_SPHERE;
int *nearp;
int *enable_interactions;
int nLayers = 0;
int nProt   = 0;
static double *cache_r = NULL;
static int cache_np;
double fdf( double *ruv, double *g );
double f( double *ruv );
double anneal( double *ruv, int npoints, double t_init, int do_short_moves );
void minimize( double *ruv, int npoints );
void writeXYZ( double *ruv, double *r, int npoints, int seq );
int pointInsidePerimeter( double *r, double *all );
int cleanBonds( int *bonds, int *nbonds, double *r, int nv, int max_bonds, int op_start  );
void sortBonds( int *bonds, int *nbonds, double *r, int nv, int max_bonds );
int findProblemPoint( int *bonds, int *nbonds, double *r, int nv, int max_bonds, int op_start  );
void MinImageHex( double *dr1, double PBC_A, double PBC_B );
void optimizeStrand( double *theta_points, int nreq, int doProtein );
double triangle_area(double*,double*,double*);
int pointInTriangle( double *pt, double *pt1, double *pt2, double *pt3  );

void getLatticeDelta( double cx, double cy, double *pcen, int nprot, double dhex, int *dx, int *dy, int *pi );
static double *whole_masks = NULL;
int nwhole_masks = 0;

static double *fixed_points = NULL;
static int *is_prot_point = NULL;

static double **protein_perimeter = NULL;
static double **protein_perimeter_seal = NULL;
static double **protein_perimeter_alt = NULL;
static double **protein_perimeter_mid = NULL;
static double **protein_perimeter_tan = NULL;
static double **protein_perimeter_der2 = NULL;

static double **mp_protein_perimeter = NULL;
static double **mp_protein_perimeter_alt = NULL;
static double **mp_protein_perimeter_mid = NULL;
static double **mp_protein_perimeter_tan = NULL;
static double **mp_protein_perimeter_der2 = NULL;



void genuv( double *ruv )
{
	ruv[0] = rand() / (double)RAND_MAX;
	ruv[1] = rand() / (double)RAND_MAX;

	double ruvo = ruv[1];

	if( doShape == BOX_SPHERE )
	{
		ruv[1] = acos(1-2*ruv[1]) / M_PI;
	}
}

double epoint( int p, double *r );
void fdiff_check( double *ruv );
void fdiff_check_particle( double *ruv, int p_outer );

void coord_grads( double u, double v, double dr_du[3], double dr_dv[3] )
{
	if( doShape == BOX_SPHERE )
	{
		dr_du[0] = -R * 2 * M_PI * sin( u * 2 * M_PI ) * sin( v* M_PI );
		dr_du[1] =  R * 2 * M_PI * cos( u * 2 * M_PI ) * sin( v* M_PI );
		dr_du[2] =  0;
		
		dr_dv[0] =  R * M_PI * cos( u * 2 * M_PI ) * cos( v* M_PI );
		dr_dv[1] =  R * M_PI * sin( u * 2 * M_PI ) * cos( v* M_PI );
		dr_dv[2] =  -M_PI * R * sin( v * M_PI );
	}
	else if( doShape == BOX_CYLINDER )
	{
		dr_du[0] = -R * 2 * M_PI * sin( u * 2 * M_PI );
		dr_du[1] =  R * 2 * M_PI * cos( u * 2 * M_PI );
		dr_du[2] =  0;
	
		dr_dv[0] = 0;
		dr_dv[1] = 0;	
		dr_dv[2] =  PeriodicLengthA;
	}
	else if( doShape == BOX_PLANE )
	{
		dr_du[0] = PeriodicLengthA;
		dr_du[1] = 0;
		dr_du[2] = 0;
		
		dr_dv[0] = 0;
		dr_dv[1] = PeriodicLengthB;
		dr_dv[2] = 0;
	}
	else if( doShape == BOX_HEX )
	{
//		double h_vec_2[2] = {  (tLx / (hex_row_dens-1)), 0 };
  //      	double h_vec_1[2] = { -(tLx / (hex_row_dens-1)) * cos( 60.0 * (M_PI/180.0) ), (tLx / (hex_row_dens-1)) * sin( 60.0 * (M_PI/180.0) ) };


		dr_du[0] = PeriodicLengthA;
		dr_du[1] = 0;
		dr_du[2] = 0;
		
		dr_dv[0] = -PeriodicLengthB * cos( 60.0 * (M_PI)/180.0);
		dr_dv[1] =  PeriodicLengthB * sin( 60.0 * (M_PI)/180.0);
		dr_dv[2] = 0;
	}

}

void map_uv_to_r( double u, double v, double r[3] )
{
	if( doShape == BOX_SPHERE )
	{
		r[0] = R * cos( u * 2 * M_PI ) * sin( v * M_PI);
		r[1] = R * sin( u * 2 * M_PI ) * sin( v * M_PI);
		r[2] = R * cos( v * M_PI );		
	}
	else if( doShape == BOX_CYLINDER )
	{
		r[0] = R * cos( u * 2 * M_PI );
		r[1] = R * sin( u * 2 * M_PI );
		r[2] = v * PeriodicLengthA;		

		while( r[2] < 0 ) r[2] += PeriodicLengthA;
		while( r[2] >= PeriodicLengthA ) r[2] -= PeriodicLengthA;
	}
	else if( doShape == BOX_PLANE )
	{
		r[0] = (u) * PeriodicLengthA;
		r[1] = (v) * PeriodicLengthB;
		r[2] = 0;
		
		while( r[0] < -PeriodicLengthA/2 ) r[0] += PeriodicLengthA;
		while( r[0] > PeriodicLengthA/2 ) r[0] -= PeriodicLengthA;
		while( r[1] < -PeriodicLengthB/2 ) r[1] += PeriodicLengthB;
		while( r[1] > PeriodicLengthB/2 ) r[1] -= PeriodicLengthB;
	}
	else if( doShape == BOX_HEX )
	{
		r[0] = (u-0.5) * PeriodicLengthA - (v-0.5) * PeriodicLengthB * cos(60.0*M_PI/180.0);
		r[1] = (v-0.5) * PeriodicLengthB * sin(60.0*M_PI/180.0);
		r[2] = 0;
	
		MinImageHex( r, PeriodicLengthA, PeriodicLengthB );	
	}

	if( ! (r[0] <0 || r[0] >-1) )
	{
		printf("r: %lf %lf %lf u: %lf v: %lf\n", r[0], r[1], r[2], u, v );
	}

}

int ASPECT = 1;

int main( int argc, char **argv )
{
	if( argc < 4 )
	{	
		printf("Syntax: makeHexCylinder.exe meshLen Radius Lc\n");
		return 0;
	}

//	int nlipids = atoi(argv[1]);
	double LC = atof(argv[3]);
	double APL = 65.0;
//	double R = (nlipids*APL)/(2*M_PI*LC);
	double R = atof(argv[2]);
	double meshLen = atof(argv[1]);
	int nlipids = R * 2 * M_PI * LC / APL;

	double Circumf = 2 * M_PI * R;
	double Height  = LC;

	char buffer[4096];
	double *points;

	double hex_spacing = meshLen; //sqrt(Circumf * LC / (nlipids/2) / sqrt(3));
	double special_L = cos(30.0*M_PI/180.0) * 2 * hex_spacing;

	int nx_int = lround(Circumf / special_L/2)*2;
	int ny_int = ASPECT * lround(LC/hex_spacing);



	double LA = nx_int * special_L;
	double LB = ny_int * hex_spacing / ASPECT;

	double zscale = LC/LB;

	double Ly_special = hex_spacing / ASPECT;

	int n_special_points = nx_int * ny_int *2;
	points = (double *)malloc( sizeof(double) * n_special_points * 3 );
	
	PeriodicLengthA = LA;
	PeriodicLengthB = LB;		

	n_special_points = 0;

	R = LA / (2*M_PI);

	for( int ix = 0; ix < nx_int; ix++ )
	{ // x is the cirumferential direction.
		for( int iy = 0; iy < ny_int; iy++ )
		{ // y is the axial direction.
			double x_val0 = ix * special_L                                     - LA/2;
			double x_val1 = ix * special_L +  hex_spacing * cos(30.0*M_PI/180.0)- LA/2; 

			double theta0 = x_val0 / (R);
			double theta1 = x_val1 / (R);
			
			double r[6] = { R * sin(theta0), R*cos(theta0), iy * Ly_special              ,
		                        R * sin(theta1), R*cos(theta1), iy * Ly_special - hex_spacing * sin(30.0*M_PI/180.0)/ASPECT  };	

/*
			double r[6] = { ix * special_L                                     - LA/2, iy * Ly_special              -LB/2, 0.0,
		                        ix * special_L +  hex_spacing * cos(30.0*M_PI/180.0)- LA/2, iy * Ly_special - hex_spacing * sin(30.0*M_PI/180.0) -LB/2, 0.0 };	
*/
			for( int xp = 0; xp < 2; xp++ )
			{	
				int skip = 0;

				points[n_special_points*3+0] = r[xp*3+0];
				points[n_special_points*3+1] = r[xp*3+1];
				points[n_special_points*3+2] = r[xp*3+2];
				n_special_points++;
			}
		}
	}

	int npoints = n_special_points;
	double *r = points;
		
	int max_bonds = 20;

	int *bonds = (int *)malloc( sizeof(int) * npoints * max_bonds );
	int *nbonds = (int*)malloc( sizeof(int) * npoints );

	memset( nbonds, 0, sizeof(int) * npoints );
	double BoxL[2] = { PeriodicLengthA, PeriodicLengthB };
	double *rn_vec = (double *)malloc( sizeof(double) * npoints * 3 );
	for( int x = 0; x < npoints; x++ )
	{
		rn_vec[3*x+0] = (1e-4) * rand() / (double)RAND_MAX;
		rn_vec[3*x+1] = (1e-4) * rand() / (double)RAND_MAX;
		rn_vec[3*x+2] = (1e-4) * rand() / (double)RAND_MAX;

		r[3*x+0] += rn_vec[3*x+0];
		r[3*x+1] += rn_vec[3*x+1];
		r[3*x+2] += rn_vec[3*x+2];
	}
	
	{
		double BoxL[2] = { PeriodicLengthB, 2 * M_PI * R };
		double *r_mod = (double *)malloc( sizeof(double) * 3 * npoints );
		for( int x = 0; x < npoints; x++ )
		{
			r_mod[3*x+0] = r[3*x+2];
			double th = atan2( r[3*x+1], r[3*x+0] );
			r_mod[3*x+1] = R * th;
		}
		get2DVoronoiConnectivity(r_mod, npoints, "unique", bonds, nbonds, max_bonds, BoxL );
		free(r_mod);
	}

	for( int x = 0; x < npoints; x++ )
	{
		r[3*x+0] -= rn_vec[3*x+0];
		r[3*x+1] -= rn_vec[3*x+1];
		r[3*x+2] -= rn_vec[3*x+2];
	}

	freopen("/dev/tty", "w", stdout );
	
	int *tri = (int *)malloc( sizeof(int) * npoints * 3 * 6 );
	int ntri = 0; 

	int *new_bonds = (int *)malloc( sizeof(int) * npoints * max_bonds );
	int *n_new_bonds = (int*)malloc(sizeof(int) * npoints );
	memset( n_new_bonds, 0, sizeof(int) * npoints );


	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < nbonds[i]; jv++ )
		{		
			int consist = 0;
			int j = bonds[i*max_bonds+jv];
		
			for( int kv = 0; kv < nbonds[j]; kv++ )
			{
				if( bonds[j*max_bonds+kv] == i )
					consist = 1;
			}

			if(! consist )
			{
				printf("Inconsistent: %d bonded to %d but not %d to %d.\n", i, j, j, i );
			}
		}
	}
	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < nbonds[i]; jv++ )
		{
			int j = bonds[i*max_bonds+jv];
			if( j <= i ) continue;

			for( int kv = 0; kv < nbonds[j]; kv++ )
			{
				int k = bonds[j*max_bonds+kv];

				if( k <= j ) continue;
				
				int gotit = 0;
			
				for( int iv = 0; iv < nbonds[k]; iv++ )
				{
					int ic = bonds[k*max_bonds+iv];

					if( ic == i ) gotit = 1;
				}		
			
				if( gotit )
				{
					int gb = 0;
		
					int ipairs[6][2] = 
					{
						{i,j}, {i,k},
						{j,i}, {j,k},
						{k,i}, {k,j},
					};

					for( int pairs = 0; pairs< 6; pairs++ )
					{		
						int check_i = ipairs[pairs][0];
						int check_j = ipairs[pairs][1];	
						gb = 0;
						for( int p = 0; p < n_new_bonds[check_i]; p++ )
						{
							if( new_bonds[check_i*max_bonds+p] == check_j )
								gb = 1;
						}
	
						if( !gb )
						{
							new_bonds[check_i*max_bonds+n_new_bonds[check_i]] = check_j;
							n_new_bonds[check_i]++;
						}
					}
					tri[3*ntri+0] = i;
					tri[3*ntri+1] = j;
					tri[3*ntri+2] = k;
					ntri++;
				}
			} 
		}
	}

	// remove points that have no neighbors.

	int npts_use = 0;
	int pt_map[npoints];

	for( int i = 0; i < npoints; i++ )
	{
		pt_map[i] = -1;
		if( n_new_bonds[i] > 0 )
		{
			pt_map[i] = npts_use;
			npts_use++;
		}
/*		else
		{
			printf("Shouldn't happen now that we are reading points.\n");
			exit(1);
		}*/

	}
	

	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < n_new_bonds[i]; jv++ )
		{		
			int consist = 0;
			int j = new_bonds[i*max_bonds+jv];
		
			for( int kv = 0; kv < n_new_bonds[j]; kv++ )
			{
				if( new_bonds[j*max_bonds+kv] == i )
					consist = 1;
			}

			if(! consist )
			{
				printf("Inconsistent: %d bonded to %d but not %d to %d.\n", i, j, j, i );
			}
		}
	}

	int nit = 0;
	int op;

	char fileName[256];
	
	sprintf(fileName,"cylindrical.mesh");

	FILE *theLattice = fopen(fileName,"w");

	fprintf(theLattice, "3d R = %lf\n", R );
	fprintf(theLattice, "%lf 0.0 0.0\n", 3*R);
	fprintf(theLattice, "0.0 %lf 0.0\n", 3*R );
	fprintf(theLattice, "0.0 0.0 %lf\n", PeriodicLengthB*zscale );

	int p = 0;
	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		fprintf(theLattice, "%d %lf %lf %lf %d", pt_map[p], r[3*p+0], r[3*p+1], r[3*p+2]*zscale , n_new_bonds[p] );

		for( int b = 0; b < n_new_bonds[p]; b++ )
			fprintf(theLattice, " %d", pt_map[new_bonds[max_bonds*p+b]] );

		fprintf(theLattice, "\n"); 
	}

	fprintf(theLattice, "ntri %d\n", ntri );

	for( int t = 0 ; t < ntri; t++ )
	{
		fprintf(theLattice, "%d %d %d", pt_map[tri[3*t+0]], pt_map[tri[3*t+1]], pt_map[tri[3*t+2]] );

		// if it's a tri-match..?

		int p1 = tri[3*t+0];
		int p2 = tri[3*t+1];
		int p3 = tri[3*t+2];

		if( p3 < p2 )
		{
			int t = p2;
			p2 = p3;
			p3 = t;
		}
		
		if( p2 < p1 )
		{
			int t = p1;
			p1 = p2;
			p2 = t;
		}
		
		if( p3 < p2 )
		{
			int t = p2;
			p2 = p3;
			p3 = t;
		}

		fprintf(theLattice, "\n");
	}
	
	printf("Writing 'cylindrical.mesh'.\n");
	fclose(theLattice);

#if 0
	for( int write_leaflet = 0; write_leaflet < 2; write_leaflet++ )
	{
		if( write_leaflet >= 0 && !doProtein ) continue;

	char fileName[256];
	if( doShape == BOX_SPHERE )
		sprintf(fileName, "sphere.%d.lattice", npoints );
	else if( doShape == BOX_CYLINDER )
		sprintf(fileName, "cylinder.%d.lattice", npoints );
	else if( doShape == BOX_PLANE || doShape == BOX_HEX )
	{
		
		if( doProtein == 1 )
			sprintf(fileName, "gen.%s.lattice", (write_leaflet==0 ? "surface" : "midplane")  );
		else if( doProtein == 2 )
			sprintf(fileName, "metaII.%s.%d.lattice", (write_leaflet==0 ? "surface" : "midplane"), npoints );
		else
			sprintf(fileName, "plane.%d.lattice", npoints );
	}

	FILE *theLattice = fopen(fileName,"w");

	fprintf(theLattice, "R %lf\n", R );

	if( doShape == BOX_SPHERE )
		fprintf(theLattice, "-1 -1 -1\n");
	else if( doShape == BOX_CYLINDER )
		fprintf(theLattice, "-1 -1 %lf\n", PeriodicLengthA );
	else if( doShape == BOX_PLANE || doShape == BOX_HEX )
	{
		if( doProtein == 0 )
			fprintf(theLattice, "%lf %lf -1\n", PeriodicLengthA, PeriodicLengthA );
		else if( doShape == BOX_PLANE )
		{
			fprintf(theLattice, "%lf 0.0 0.0\n", PeriodicLengthA );
			fprintf(theLattice, "0.0 %lf 0.0\n", PeriodicLengthB );
		}
		else
		{
			fprintf(theLattice, "%lf 0.0 0.0\n", PeriodicLengthA );
			fprintf(theLattice, "%lf %lf 0.0\n", -PeriodicLengthB*cos(60.0*M_PI/180.0), PeriodicLengthB*sin(60.0*M_PI/180.0) );
		}
	}

	p = 0;
	for( int pr = 0; pr < nProt; pr++ )
	{
		for( int tp = 0; tp < nLPts[pr]; tp++, p++)
		{	
			if( pt_map[p] < 0 )
				continue;

			int up = pt_map[p];
	
			if( is_prot_point[p] )
				fprintf(theLattice, "Q%d ", up);
			else
				fprintf(theLattice, "%d ", up );
			if( write_leaflet == 0 )
			{
				fprintf(theLattice, " %lf %lf %lf %d", r[3*p+0], 
								       r[3*p+1], fixed_points[6*p+2], 
					n_new_bonds[p] );
			}
			else
			{
				fprintf(theLattice, " %lf %lf %lf %d",  fixed_points[6*p+3], 	
										     fixed_points[6*p+4], 
											0.0, 
					n_new_bonds[p] );
			}
			for( int b = 0; b < n_new_bonds[p]; b++ )
				fprintf(theLattice, " %d", pt_map[new_bonds[max_bonds*p+b]] );
			fprintf(theLattice, "\n"); 
		}	
	}

	for( ; p < npoints; p++ )
	{
		fprintf(theLattice, "%d %lf %lf %lf %d", pt_map[p], r[3*p+0], r[3*p+1], r[3*p+2] + (write_leaflet == 0  ? avz : 0 ) , n_new_bonds[p] );

		for( int b = 0; b < n_new_bonds[p]; b++ )
			fprintf(theLattice, " %d", pt_map[new_bonds[max_bonds*p+b]] );
		fprintf(theLattice, "\n"); 
	}

	fprintf(theLattice, "ntri %d\n", ntri );

	for( int t = 0 ; t < ntri; t++ )
	{
		fprintf(theLattice, "%d %d %d\n", pt_map[tri[3*t+0]], pt_map[tri[3*t+1]], pt_map[tri[3*t+2]] );
	}

	fprintf(theLattice, "edge_sense %d\n", nedge_sense );

	for( int s = 0; s < nedge_sense; s++ )
		fprintf(theLattice,"%d %d %d\n", pt_map[saved_sense[3*s+0]], pt_map[saved_sense[3*s+1]], pt_map[saved_sense[3*s+2]] ); 

	fclose(theLattice);
	}

#endif
}


void MinImageHex( double *dr1, double PBC_A, double PBC_B)
{
	double PBC_vec[2][3] =
	{
		{PBC_A, 0, 0},
		{-PBC_B*cos(60.0*M_PI/180.0), PBC_B*sin(60.0*M_PI/180.0), 0 }
	};
	int done = 0;


	double tr[3] = { dr1[0], dr1[1], dr1[2] };

	while( !done )
	{
		done = 1;

		double r2 = tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2];

		for( int dx = -1; dx <= 1; dx++ )
		for( int dy = -1; dy <= 1; dy++ )
		{
			if( dx == 0 && dy == 0 ) continue;
			tr[0] += dx * PBC_vec[0][0] + dy * PBC_vec[1][0];
			tr[1] += dx * PBC_vec[0][1] + dy * PBC_vec[1][1];
			tr[2] += dx * PBC_vec[0][2] + dy * PBC_vec[1][2];

			double nr2 = tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2];

			if( nr2 >= r2 )
			{
				tr[0] -= dx * PBC_vec[0][0] + dy * PBC_vec[1][0];
				tr[1] -= dx * PBC_vec[0][1] + dy * PBC_vec[1][1];
				tr[2] -= dx * PBC_vec[0][2] + dy * PBC_vec[1][2];

			}
			else	
			{
				r2 = nr2;
				
				done = 0;
			}
		}
	}

	dr1[0] = tr[0];
	dr1[1] = tr[1];
	dr1[2] = tr[2];
}


double MinImageHexDel( double *dr1, double PBC_len, int *out_dx, int *out_dy)
{
	*out_dx = 0;
	*out_dy = 0;
	double PBC_vec[2][3] =
	{
		{PBC_len, 0, 0},
		{-PBC_len*cos(60.0*M_PI/180.0), PBC_len*sin(60.0*M_PI/180.0), 0 }
	};
	int done = 0;


	double tr[3] = { dr1[0], dr1[1], dr1[2] };
		
	double r2;

	while( !done )
	{
		done = 1;

		r2 = tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2];

		for( int dx = -1; dx <= 1; dx++ )
		for( int dy = -1; dy <= 1; dy++ )
		{
			if( dx == 0 && dy == 0 ) continue;
			tr[0] += dx * PBC_vec[0][0] + dy * PBC_vec[1][0];
			tr[1] += dx * PBC_vec[0][1] + dy * PBC_vec[1][1];
			tr[2] += dx * PBC_vec[0][2] + dy * PBC_vec[1][2];

			double nr2 = tr[0]*tr[0]+tr[1]*tr[1]+tr[2]*tr[2];

			if( nr2 >= r2 )
			{
				tr[0] -= dx * PBC_vec[0][0] + dy * PBC_vec[1][0];
				tr[1] -= dx * PBC_vec[0][1] + dy * PBC_vec[1][1];
				tr[2] -= dx * PBC_vec[0][2] + dy * PBC_vec[1][2];

			}
			else	
			{
				*out_dx += dx;
				*out_dy += dy;

				r2 = nr2;
				
				done = 0;
			}
		}
	}

	return r2;
}

void getLatticeDelta( double cx, double cy, double *pcen, int nprot, double dhex, int *dx, int *dy, int *pi )
{
	double minp = 1e10;

	for( int p = 0; p < nprot; p++ )
	{	
		double del[3] = { cx - pcen[2*p+0], cy - pcen[2*p+1], 0 };	
		int out_dx, out_dy;
		
		if( doShape == BOX_HEX )	
			MinImageHex( del, PeriodicLengthA, PeriodicLengthB );	
		else
		{
			while( del[0] < -PeriodicLengthA/2) del[0] += PeriodicLengthA;
			while( del[0] >  PeriodicLengthA/2) del[0] -= PeriodicLengthA;
			while( del[1] < -PeriodicLengthB/2) del[1] += PeriodicLengthB;
			while( del[1] >  PeriodicLengthB/2) del[1] -= PeriodicLengthB;
		}
		double r2 = del[0]*del[0]+del[1]*del[1];
		double rcheck = MinImageHexDel( del, dhex, &out_dx, &out_dy);

		if( rcheck > 1e-2 )
		{
			printf("what.\n");
			exit(1);
		}

		if( r2 < minp )
		{
			*pi = p;
			minp = r2;
			*dx = out_dx;
			*dy = out_dy;
		}
	}
}

