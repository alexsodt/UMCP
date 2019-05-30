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

#define MAX_EDGES 10
int max_bonds = MAX_EDGES;

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
void optimizeStrand( double *theta_points, int nreq, int doProtein );
double triangle_area(double*,double*,double*);
int pointInTriangle( double *pt, double *pt1, double *pt2, double *pt3  );

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




double epoint( int p, double *r );
void fdiff_check( double *ruv );
void fdiff_check_particle( double *ruv, int p_outer );

void getPlanarMesh( double **rvals, int **edges, int **nedges, int *nv, double Lx, int NXI );
void getCylinderMesh( double **rvals, int **edges, int **nedges, int *nv, int grain, double R, double LC);

void getPath( int start, int *path, int *pl_out, int *ne_arr, int *e_arr, int *interior_tag, int nr );

int main( int argc, char **argv )
{
	if( argc < 5 )
	{	
		printf("Syntax: makePore NX LX neckRadius neckHeight\n");
		return 0;
	}

	int NX = atoi( argv[1] );
	double Lx = atof( argv[2] );
	double *r_pl;
	int *e_pl;
	int *ne_pl;
	int nr_pl;
	int nr_ex = 0;


	printf("Lx: %lf nx: %d\n", Lx, NX );
	getPlanarMesh( &r_pl, &e_pl, &ne_pl, &nr_pl, Lx, NX ); 

	double neckR = atof(argv[3]);
	double neckH = atof(argv[4]);

	int *interior_tag = (int *)malloc( sizeof(int) * nr_pl );

	memset( interior_tag, 0, sizeof(int) * nr_pl );

	for( int v = 0; v < nr_pl; v++ )
	{
		if( r_pl[3*v+0]*r_pl[3*v+0] + r_pl[3*v+1] * r_pl[3*v+1] < neckR*neckR )
			interior_tag[v] = 1; 
		else
			nr_ex++;
	}	
	int start = -1;
	for( int v = 0; v < nr_pl && start == -1; v++ )
	{
		if( interior_tag[v] != 0 ) continue;

		for( int e = 0; e < ne_pl[v]; e++ )
		{
			if( interior_tag[e_pl[MAX_EDGES*v+e]] ) 
				start = v;
		}
	}
	
	int *path = (int *)malloc( sizeof(int) * nr_pl );
	int pl = 0;
	getPath( start, path, &pl, ne_pl, e_pl, interior_tag, nr_pl );

	FILE *xyzFile = fopen("proto.xyz","w");

	for( int v = 0; v < nr_pl; v++ )
	{
		if( !interior_tag[v] )
		{
			fprintf(xyzFile,"C %lf %lf %lf\n", r_pl[3*v+0], r_pl[3*v+1], r_pl[3*v+2] ); 
		}
		else
			fprintf(xyzFile, "O %lf %lf %lf\n", r_pl[3*v+0], r_pl[3*v+1], r_pl[3*v+2] );
	}

	for( int vx = 0; vx < pl; vx++ )
	{
		int v = path[vx];

		fprintf(xyzFile, "N %lf %lf %lf\n", r_pl[3*v+0], r_pl[3*v+1], r_pl[3*v+2] );
	}

	double *r_cyl;
	int *e_cyl;
	int *ne_cyl;
	int nr_cyl;

	int done = 0;
	int *pathc;
	int *cyl_interior_tag;
	int cpl;

	int cur_NX = NX;
	int move_dir = 0;
	double best = 1e10;
	int NXBest = -1;
	int repeat = 0;
	int nr_cyl_ex = 0;

	while( !done && !repeat )
	{
		getCylinderMesh( &r_cyl, &e_cyl, &ne_cyl, &nr_cyl, cur_NX, neckR, neckH * 1.2 );
		
		cyl_interior_tag = (int *)malloc( sizeof(int) * nr_cyl );
		memset( cyl_interior_tag, 0, sizeof(int) * nr_cyl );
		
		for( int ir = 0; ir < nr_cyl; ir++ )
			r_cyl[3*ir+2] -= neckH*0.05;
		
		nr_cyl_ex = 0;	
		for( int v = 0; v < nr_cyl; v++ )
		{
			if( r_cyl[3*v+2] < 0 || r_cyl[3*v+2] > neckH )
				cyl_interior_tag[v] = 1; 
			else
				nr_cyl_ex++;
		}	
		
		pathc = (int *)malloc( sizeof(int) * nr_cyl );

		start = -1;
		for( int v = 0; v < nr_cyl && start == -1; v++ )
		{
			if( cyl_interior_tag[v] != 0 ) continue;
			if( r_cyl[3*v+2] > neckH /2 ) continue;
	
			for( int e = 0; e < ne_cyl[v]; e++ )
			{
				if( cyl_interior_tag[e_cyl[MAX_EDGES*v+e]] ) 
					start = v;
			}
		}

		getPath( start, pathc, &cpl, ne_cyl, e_cyl, cyl_interior_tag, nr_cyl );

		if( fabs(cpl-(double)pl) < best )
		{
			best = fabs(cpl-(double)pl);
			NXBest = cur_NX;
		}

		if( repeat )
		{
			done = 1;
			repeat = 0;
		}
		else if( cpl != pl )
		{
			if( cpl > pl )
			{
				if( move_dir == 1 )
				{
					done = 1;
					repeat = NXBest;		
				}	
				else
					cur_NX--;
				move_dir  = -1;
			}
			else if( move_dir == -1 )
			{
				done = 1;
				repeat = NXBest;			
			}
			else
			{
				cur_NX++;
				move_dir  = 1;
			}
		}
		else
		{
			done = 1;
		}
		if( !done )
		{
			free(pathc);
			free(r_cyl);
			free(e_cyl);
			free(ne_cyl);
			free(cyl_interior_tag);
		}
		if( repeat == cur_NX )
			repeat = 0;
		if( repeat )
			cur_NX = repeat;
	}
		
	start = -1;
	for( int v = 0; v < nr_cyl && start == -1; v++ )
	{
		if( cyl_interior_tag[v] != 0 ) continue;
		if( r_cyl[3*v+2] < neckH /2 ) continue;
	
		for( int e = 0; e < ne_cyl[v]; e++ )
		{
			if( cyl_interior_tag[e_cyl[MAX_EDGES*v+e]] ) 
				start = v;
		}
	}

	int cpl2;
	int *pathc2 = (int *)malloc( sizeof(int) * nr_cyl );
	getPath( start, pathc2, &cpl2, ne_cyl, e_cyl, cyl_interior_tag, nr_cyl );

	printf("PL: %d CPL: %d CPL2: %d\n", pl, cpl, cpl2 );

	for( int v = 0; v < nr_cyl; v++ )
	{
		if( !cyl_interior_tag[v] )
		{
			fprintf(xyzFile,"C %lf %lf %lf\n", r_cyl[3*v+0], r_cyl[3*v+1], r_cyl[3*v+2] ); 
		}
		else
			fprintf(xyzFile,"O %lf %lf %lf\n", r_cyl[3*v+0], r_cyl[3*v+1], r_cyl[3*v+2] );
	}

	for( int vx = 0; vx < cpl; vx++ )
	{
		int v = pathc[vx];

		fprintf(xyzFile,"H %lf %lf %lf\n", r_cyl[3*v+0], r_cyl[3*v+1], r_cyl[3*v+2] );
	}

	fclose(xyzFile);
	
	
	/* connect the first loop */
	int startcyl[2] = {-1,-1};

	int v = path[0];

	for( int layer = 0; layer < 2; layer++ )
	{
		double bestr2 = 1e20;

		for( int x = 0; x < cpl; x++ )
		{
			int v2 = pathc[x];
			if( layer == 1 )
				v2 = pathc2[x];
	
			double dr[3] = { r_cyl[3*v2+0] - r_pl[3*v+0], r_cyl[3*v2+1] - r_pl[3*v+1], r_cyl[3*v2+2] - (r_pl[3*v+2]+layer*neckH) };
			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	
			if( r < bestr2 )
			{
				startcyl[layer] = x;
				bestr2 = r;
			}
		}
		
	
		// which direction to go.
	
		int next3 = startcyl[layer]+3;
		if( next3 > cpl ) next3 -= cpl;
		int prev3 = startcyl[layer]-3;
		if( prev3 < 0 ) prev3 += cpl;

		int v2n,v2p;

		if( layer == 0 )
		{	
			v2n = pathc[next3];
			v2p = pathc[prev3];
		}
		if( layer == 1 )
		{
			v2n = pathc2[next3];
			v2p = pathc2[prev3];

		}

		double rd1[3] = { 
				r_cyl[v2n*3+0] - r_pl[3*(path[3])+0],
				r_cyl[v2n*3+1] - r_pl[3*(path[3])+1],
				r_cyl[v2n*3+2] - (r_pl[3*(path[3])+2]+layer*neckH) };
		
		double rd2[3] = { 
				r_cyl[v2p*3+0] - r_pl[3*(path[3])+0],
				r_cyl[v2p*3+1] - r_pl[3*(path[3])+1],
				r_cyl[v2p*3+2] - (r_pl[3*(path[3])+2]+layer*neckH) };
		
		double rv1 = rd1[0]*rd1[0]+rd1[1]*rd1[1]+rd1[2]*rd1[2];
		double rv2 = rd2[0]*rd2[0]+rd2[1]*rd2[1]+rd2[2]*rd2[2];
	
		if( rv1 < rv2 )
		{
			printf("layer %d Going forward %lf vs %lf\n", layer, rv1, rv2 );
			// forward
		}
		else
		{
			printf("layer %d Going backward %lf vs %lf\n", layer, rv2, rv1 );

			if( layer == 0 )
			{
				// reverse
				int copy[cpl];
				int news = 0;	
				for( int x = 0; x < cpl; x++ )
				{
					copy[cpl-x-1] = pathc[x];
					if( x == startcyl[layer] )
						news = cpl-x-1;
				}
				startcyl[layer] = news;
				memcpy( pathc, copy, sizeof(int) * cpl );
			}
			else
			{
				// reverse
				int copy[cpl2];
				int news = 0;	
				for( int x = 0; x < cpl2; x++ )
				{
					copy[cpl2-x-1] = pathc2[x];
					if( x == startcyl[layer] )
						news = cpl2-x-1;
				}
				startcyl[layer] = news;
				memcpy( pathc2, copy, sizeof(int) * cpl2 );
			}
		}
	}

	int nvtot = 2*nr_ex + nr_cyl_ex;
	double *allr = (double *)malloc( sizeof(double) * 3 * nvtot );
	int    *alle = (int *)malloc( sizeof(int) * nvtot * MAX_EDGES );
	int    *ne = (int *)malloc( sizeof(int) * nvtot );
	memset( ne, 0, sizeof(int) * nvtot );
	int *pl_map = (int *)malloc( sizeof(int) * nr_pl );
	int *pl_map2 = (int *)malloc( sizeof(int) * nr_pl );
	int *cyl_map = (int *)malloc( sizeof(int) * nr_cyl );

	int ii = 0;

	for( int v = 0; v < nr_pl; v++ )
	{
		pl_map[v] = -1;
		pl_map2[v] = -1;
		if( interior_tag[v] ) continue;

 		pl_map[v] = ii;
		pl_map2[v] = ii + nr_ex;
		ii++;
	}
	
	printf("ii: %d nr_ex: %d\n", ii, nr_ex );
	
	ii = 0;

	for( int v = 0; v < nr_cyl; v++ )
	{
		cyl_map[v] = -1;
		if( cyl_interior_tag[v] ) continue;

 		cyl_map[v] = 2*nr_ex+ii;
		ii++;
	}	
	printf("ii: %d nr_cyl_ex: %d\n", ii, nr_cyl_ex );

	int *path2 = (int *)malloc( sizeof(int) * pl );
	memcpy( path2, path, sizeof(int) * pl );

	for( int px = 0; px < pl; px++ )
	{
		path[px]  = pl_map[path[px]];
		path2[px] = pl_map2[path2[px]];
	}

	for( int px = 0; px < cpl; px++ )
		pathc[px] = cyl_map[pathc[px]];
	for( int px = 0; px < cpl2; px++ )
		pathc2[px] = cyl_map[pathc2[px]];
	
	for( int vx = 0; vx < nr_pl; vx++ )
	{
		int v = pl_map[vx];

		if( v >= 0 )
		{	
			allr[3*v+0] = r_pl[3*vx+0];
			allr[3*v+1] = r_pl[3*vx+1];
			allr[3*v+2] = r_pl[3*vx+2];
			
			for( int e = 0; e < ne_pl[vx]; e++ )
			{
				if( pl_map[e_pl[vx*MAX_EDGES+e]] >= 0 )
				{
					alle[v*MAX_EDGES+ne[v]] = pl_map[e_pl[vx*MAX_EDGES+e]];
					ne[v]++; 
				}
			}
		}
	}
	
	for( int vx = 0; vx < nr_pl; vx++ )
	{
		int v = pl_map2[vx];

		if( v >= 0 )
		{	
			allr[3*v+0] = r_pl[3*vx+0];
			allr[3*v+1] = r_pl[3*vx+1];
			allr[3*v+2] = r_pl[3*vx+2] + neckH;
			
			for( int e = 0; e < ne_pl[vx]; e++ )
			{
				if( pl_map2[e_pl[vx*MAX_EDGES+e]] >= 0 )
				{
					alle[v*MAX_EDGES+ne[v]] = pl_map2[e_pl[vx*MAX_EDGES+e]];
					ne[v]++; 
				}
			}
		}
	}
	
	for( int vx = 0; vx < nr_cyl; vx++ )
	{
		int v = cyl_map[vx];

		if( v >= 0 )
		{	

			allr[3*v+0] = r_cyl[3*vx+0];
			allr[3*v+1] = r_cyl[3*vx+1];
			allr[3*v+2] = r_cyl[3*vx+2];
			
			for( int e = 0; e < ne_cyl[vx]; e++ )
			{
				if(  cyl_map[e_cyl[vx*MAX_EDGES+e]] >= 0 )
				{
					alle[v*MAX_EDGES+ne[v]] = cyl_map[e_cyl[vx*MAX_EDGES+e]];
					ne[v]++; 
				}
			}
		}
	}

	// loop over both paths.

	for( int layer = 0; layer < 2; layer++ )
	{
		int *mpath = path;

		if( layer == 1 )
			mpath = path2;

		int p2 = startcyl[layer];
		int clim = cpl;
		if( layer == 1 )
			clim = cpl;
		int *alt_path = pathc;

		if( layer == 1 )
			alt_path = pathc2;

		if( pl <= clim )
		{
			// cylinder path is longer
			for( int p = 0; p < pl; p++, p2++ )
			{
				if( p2 >= clim )
					p2 -= clim;
		
				int v = mpath[p];
				
				int vcp = alt_path[(p2-1<0?p2-1+clim:p2-1)];
				int vc = alt_path[p2];
		 
				alle[v*MAX_EDGES+ne[v]] = vcp; ne[v]++;
				alle[v*MAX_EDGES+ne[v]]  = vc; ne[v]++;
			
				alle[vcp*MAX_EDGES+ne[vcp]] = v; ne[vcp]++;	
				alle[vc*MAX_EDGES+ne[vc]]   = v; ne[vc]++;	
			}
				
			for( int px = pl; px < clim; px++, p2++ )
			{
				int vc = mpath[p2];		
				int v = alt_path[pl-1];
		 
				alle[vc*MAX_EDGES+ne[vc]]  =    v; ne[vc]++;
				alle[v*MAX_EDGES+ne[v]]    =   vc; ne[v]++;	
			}
		}
		else
		{
			// planar path is longer.
	
			for( int p = 0; p < clim; p++, p2++ )
			{
				if( p2 >= clim )
					p2 -= clim;
		
				int vc = alt_path[p2];
				
				int vp = mpath[(p-1<0?p-1+pl:p-1)];
				int v = mpath[p];
		 
				alle[vc*MAX_EDGES+ne[vc]]   = vp; ne[vc]++;
				alle[vc*MAX_EDGES+ne[vc]]  =   v; ne[vc]++;
			
				alle[vp*MAX_EDGES+ne[vp]] = vc; ne[vp]++;	
				alle[v*MAX_EDGES+ne[v]]   = vc; ne[v]++;	
			}

			p2--;
	
			for( int p = clim; p < pl; p++ )
			{
				int vc = alt_path[p2];		
				int v = mpath[p];
		 
				alle[vc*MAX_EDGES+ne[vc]]  =    v; ne[vc]++;
				alle[v*MAX_EDGES+ne[v]]    =   vc; ne[v]++;	
			}
		}
	}	
		
	char fileName[256];
	sprintf(fileName,"fusion.mesh" );

	FILE *theLattice = fopen(fileName,"w");

	fprintf(theLattice, "3d R = %lf\n", R );
	fprintf(theLattice, "%lf 0.0 0.0\n", Lx );
	fprintf(theLattice, "0.0 %lf 0.0\n", Lx );
	fprintf(theLattice, "0.0 0.0 %lf\n", 2*Lx + neckH );

	int p = 0;
	for( ; p < nvtot; p++ )
	{
		if( ne[p] < 5 || ne[p]  > 7 )
			printf("ATTN: Valence %d.\n", ne[p] );
		fprintf(theLattice, "%d %lf %lf %lf %d", p, allr[3*p+0], allr[3*p+1], allr[3*p+2], ne[p] );

		for( int b = 0; b < ne[p]; b++ )
			fprintf(theLattice, " %d", alle[max_bonds*p+b] );

		fprintf(theLattice, "\n"); 
	}

	//fprintf(theLattice, "ntri %d\n", ntri );

//	for( int t = 0 ; t < ntri; t++ )
	{
		//fprintf(theLattice, "%d %d %d", pt_map[tri[3*t+0]], pt_map[tri[3*t+1]], pt_map[tri[3*t+2]] );
		//fprintf(theLattice, "\n");
	}
	
	printf("Writing '%s'.\n", fileName );
	fclose(theLattice);

	
}

void getPlanarMesh( double **rvals, int **edges, int **nedges, int *nv, double Lx, int NXI )
{
	int nx_int = NXI;
	int ny_int;
	PeriodicLengthA = Lx;

	char buffer[4096];
	double *points;

	double special_L = PeriodicLengthA / nx_int;
	double hex_spacing = special_L / ( 2*cos(30.0*M_PI/180.0));
	double pt_spacing = 0.5 * hex_spacing / (sqrt(3.0)/2);
	double Lc = -1;

	double LB = ny_int * hex_spacing;
	ny_int = lround( PeriodicLengthA / hex_spacing);
	double unscaled_LB = ny_int * hex_spacing;
	double scale_B = PeriodicLengthA / unscaled_LB;

	double Ly_special = hex_spacing;
	int n_special_points = nx_int * ny_int *2;
	points = (double *)malloc( sizeof(double) * n_special_points * 3 );
	
	PeriodicLengthB = PeriodicLengthA;	

	n_special_points = 0;

	for( int ix = 0; ix < nx_int; ix++ )
	{
		for( int iy = 0; iy < ny_int; iy++ )
		{
			double r[6] = { ix * special_L                                     - PeriodicLengthA/2, (iy * Ly_special              -PeriodicLengthB/2)*scale_B, 0.0,
		                        ix * special_L +  hex_spacing * cos(30.0*M_PI/180.0)- PeriodicLengthA/2, (iy * Ly_special - hex_spacing * sin(30.0*M_PI/180.0) -PeriodicLengthB/2)*scale_B, 0.0 };	

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
	get2DVoronoiConnectivity( r, npoints, "unique", bonds, nbonds, max_bonds, BoxL );
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

	printf("bond consistency check.\n");

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
	
	printf("bond consistency check.\n");

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

//	fprintf(theLattice, "%lf 0.0 0.0\n", PeriodicLengthA );
//	fprintf(theLattice, "0.0 %lf 0.0\n", PeriodicLengthB );
	int p = 0;

	int npt = 0;
	int net = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		npt++;

		net += n_new_bonds[p];
	}

	(*rvals)  = (double *)malloc( sizeof(double) * 3 * npt );
	(*edges)  = (int *)malloc( sizeof(int) * MAX_EDGES * npt );
	(*nedges) = (int *)malloc( sizeof(int) * npt );
	
	npt = 0;

	p = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		(*rvals)[3*npt+0] = r[3*p+0];
		(*rvals)[3*npt+1] = r[3*p+1];
		(*rvals)[3*npt+2] = r[3*p+2];

		for( int px = 0; px < n_new_bonds[p]; px++ )
			(*edges)[npt*MAX_EDGES+px] = pt_map[new_bonds[max_bonds*p+px]];

		(*nedges)[npt] = n_new_bonds[p];

		npt++;
	}

	*nv = npt;	
	printf("npt: %d\n", npt );
} 
	
void getCylinderMesh( double **rvals, int **edges, int **nedges, int *nv, int grain, double R, double LC)
{ 

	int ny_int = grain;
	double APL = 65.0;
	//double R = (nlipids*APL)/(2*M_PI*LC);

	double Circumf = 2 * M_PI * R;
	double Height  = LC;

	char buffer[4096];
	double *points;

	double hex_spacing = LC / ny_int;
//	double hex_spacing = sqrt(Circumf * LC / (nlipids/2) / sqrt(3));
	double pt_spacing = 0.5 * hex_spacing / (sqrt(3.0)/2);
	double special_L = cos(30.0*M_PI/180.0) * 2 * hex_spacing;

	int nx_int = lround(Circumf / special_L/2)*2;
//	int ny_int = lround(LC/hex_spacing);



	double LA = nx_int * special_L;
	double LB = ny_int * hex_spacing;

	double Ly_special = hex_spacing;

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
		                        R * sin(theta1), R*cos(theta1), iy * Ly_special - hex_spacing * sin(30.0*M_PI/180.0)  };	

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

	int p = 0;

	int npt = 0;
	int net = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		npt++;

		net += n_new_bonds[p];
	}

	(*rvals)  = (double *)malloc( sizeof(double) * 3 * npt );
	(*edges)  = (int *)malloc( sizeof(int) * MAX_EDGES * npt );
	(*nedges) = (int *)malloc( sizeof(int) * npt );
	
	npt = 0;

	p = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		(*rvals)[3*npt+0] = r[3*p+0];
		(*rvals)[3*npt+1] = r[3*p+1];
		(*rvals)[3*npt+2] = r[3*p+2];

		for( int px = 0; px < n_new_bonds[p]; px++ )
			(*edges)[npt*MAX_EDGES+px] = pt_map[new_bonds[max_bonds*p+px]];

		(*nedges)[npt] = n_new_bonds[p];

		npt++;
	}
	
	*nv = npt;	
	printf("npt: %d\n", npt );

	free(points);
	free(bonds);
	free(new_bonds);
	free(n_new_bonds);
	free(tri);
/*
	char fileName[256];
	
	sprintf(fileName,"cylindrical.mesh");

	FILE *theLattice = fopen(fileName,"w");

	fprintf(theLattice, "3d R = %lf\n", R );
	fprintf(theLattice, "%lf 0.0 0.0\n", 4 * R );
	fprintf(theLattice, "0.0 %lf 0.0\n", 4 * R );
	fprintf(theLattice, "0.0 0.0 %lf\n", PeriodicLengthB );

	int p = 0;
	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		fprintf(theLattice, "%d %lf %lf %lf %d", pt_map[p], r[3*p+0], r[3*p+1], r[3*p+2] , n_new_bonds[p] );

		for( int b = 0; b < n_new_bonds[p]; b++ )
			fprintf(theLattice, " %d", pt_map[new_bonds[max_bonds*p+b]] );

		fprintf(theLattice, "\n"); 
	}

	fprintf(theLattice, "ntri %d\n", ntri );
*/
}
	
void getPath( int start, int *path, int *pl_out, int *ne_arr, int *e_arr, int *interior_tag, int nr  )
{
	int pl = 0;

	for( int v = 0; v < nr && start == -1; v++ )
	{
		if( interior_tag[v] ) continue;

		for( int e = 0; e < ne_arr[v] && start == -1; e++ )
		{
			int v2 = e_arr[v*MAX_EDGES+e];

			if( interior_tag[v2] )
			{
				start = v;
			}
		}
	}

	if( start == -1 )
	{
		printf("Couldn't find a place to start.\n");
		exit(1);
	}

	int done = 0;
	
	int *path_tag = (int*)malloc( sizeof(int) * nr );
	memset( path_tag, 0, sizeof(int) * nr );
	path_tag[start] = 1;
	path[pl] = start;
	pl = 1;

	int cur = start;

	while( !done )
	{	
		int ok_stop = 0;

		int v = cur;
		int gotit = -1;

		for( int e = 0; e < ne_arr[v]; e++ )
		{
			int v2 = e_arr[v*MAX_EDGES+e];

			if( interior_tag[v2] )
			{
				for( int e2 = 0; e2 < ne_arr[v2]; e2++ )
				{
					int v3 = e_arr[v2*MAX_EDGES+e2];


					if( v3 == start )
						ok_stop = 1;

					if( path_tag[v3] ) continue;
					if( interior_tag[v3]) continue;
			
					// prospective, v3;

					for( int e3 = 0; e3 < ne_arr[v3]; e3++ )
					{
						if( e_arr[v3*MAX_EDGES+e3] == v )
						{	
							gotit = v3;
						}
					}
				}
			}
		}	
	
		if( gotit >= 0 )	
		{
			path[pl] = gotit;
			path_tag[gotit] = 1;
			pl++;
			cur = gotit;
		}
		else if( ok_stop )
		{
			printf("Successful circuit!! length %d\n", pl);
			done = 1;
		}	
		else
		{
			printf("Failed.\n");
			done = 1;
		}
		if( pl >= nr )		
		{
			printf("Bugged, looping too many times.\n");
			exit(1);
		}		
	}
	*pl_out = pl;
}
