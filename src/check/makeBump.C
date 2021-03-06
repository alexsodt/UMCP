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

const char *program_base = "/Users/sodtaj/git_projects/BD_Membrane";

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
void getSphereMesh( double **rvals, int **edges, int **nedges, int *nv, int ndiv, double R);
void getPath( int start, int *path, int *pl_out, int *ne_arr, int *e_arr, int *interior_tag, int nr );

int main( int argc, char **argv )
{
	if( argc < 4 )
	{	
		printf("Syntax: makeBud NX LX SphereRadius\n");
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

	double bumpR = atof(argv[3]);

	int *interior_tag = (int *)malloc( sizeof(int) * nr_pl );

	memset( interior_tag, 0, sizeof(int) * nr_pl );

	for( int v = 0; v < nr_pl; v++ )
	{
		if( r_pl[3*v+0]*r_pl[3*v+0] + r_pl[3*v+1] * r_pl[3*v+1] < bumpR*bumpR )
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
	// this is the planar path
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

	// now the spherical mesh works with the original path.

	// approximate ntri is two times the number of edges
	int ntri_sph = 12 * 6 / 3; // * pow( 4, nsub );
	// edges per unit area.
	double tri_density = ntri_sph / (4*M_PI*bumpR*bumpR);
	double planar_density = nr_pl * 2 / (Lx*Lx);
	int nsub_try = ceil( log(  planar_density / tri_density ) / log(4.0) );

	double dens1 = pow( 4, nsub_try - 1 ) * 12 * 6 / 3 / (4*M_PI*bumpR*bumpR);
	double dens2 = pow( 4, nsub_try     ) * 12 * 6 / 3 / (4*M_PI*bumpR*bumpR);
	double dens3 = pow( 4, nsub_try + 1 ) * 12 * 6 / 3 / (4*M_PI*bumpR*bumpR);

	printf("pdens %lf dens1 %lf dens2 %lf dens3 %lf\n", planar_density, dens1, dens2, dens3 );

	int nsub = nsub_try -1;
	double sdens = dens1;

	if( fabs(dens2/planar_density-1) < fabs(sdens/planar_density-1) )
	{
		nsub = nsub_try;
		sdens = dens2;
	}
	
	if( fabs(dens3/planar_density-1) < fabs(sdens/planar_density-1) )
	{
		nsub = nsub_try+1;
		sdens = dens3;
	}
	
	printf("selecting nsub %d, density of sphere is %lf times the planar density.\n",
		nsub, sdens / planar_density );

	printf("nsub: %d\n", nsub ); 

	double *r_sph;
	int *e_sph;
	int *ne_sph;
	int nr_sph;
	int nr_sph_ex;
	int *paths;

	getSphereMesh( &r_sph, &e_sph, &ne_sph, &nr_sph, nsub, bumpR );
	int *sph_interior_tag = (int *)malloc( sizeof(int) * nr_sph );

	double zcen = 0;

	int spl = 0;
	int repeat = 0;
	double zbest = zcen;
	double best = 1e10;
	
	int move_dir = 0;
	int done = 0;

	while( ! done )
	{
		memset( sph_interior_tag, 0, sizeof(int) * nr_sph );
		
		for( int ir = 0; ir < nr_sph; ir++ )
			r_sph[3*ir+2] += zcen;
		
		nr_sph_ex = 0;	

		for( int v = 0; v < nr_sph; v++ )
		{
			if( r_sph[3*v+2] < 0 )
				sph_interior_tag[v] = 1; 
			else
				nr_sph_ex++;
//			printf("rsp %lf %lf %lf %s\n", r_sph[3*v+0], r_sph[3*v+1], r_sph[3*v+2], ( sph_interior_tag[v] ? "interior" : "exterior") );
		}	
		

		paths = (int *)malloc( sizeof(int) * nr_sph );

		getPath( -1, paths, &spl, ne_sph, e_sph, sph_interior_tag, nr_sph );

		break;

//		printf("ZCEN: %lf SPL: %d\n", zcen, spl );

		for( int ir = 0; ir < nr_sph; ir++ )
			r_sph[3*ir+2] -= zcen;

		if( fabs(spl-(double)pl) < best )
		{
			best = fabs(spl-(double)pl);
			zbest = zcen;
		}

		if( repeat )
		{
			done = 1;
			repeat = 0;
		}
		else if( spl != pl )
		{
			if( spl > pl )
			{
				if( move_dir == -1 )
				{
					done = 1;
					repeat = 1;		
				}	
				else
					zcen += 1;

				move_dir  = 1;
			}
			else if( move_dir == 1 )
			{
				done = 1;
				repeat = 1;			
			}
			else
			{
				zcen -= 1;
				move_dir  = -1;
			}
		}
		else
		{
			done = 1;
		}

		if( !done )
			free(paths);

		if( repeat )
			zcen = zbest;
	}
		
	for( int ir = 0; ir < nr_sph; ir++ )
		r_sph[3*ir+2] += zcen;

	fclose(xyzFile);
	
	
	/* connect the first loop */
	int startsph = -1;

	{
		double *r_alt = r_pl;
		int alt_start = path[0];	
		int alt_next  = path[3];
		int clim = spl;

		int v = alt_start;

		double bestr2 = 1e20;

		for( int x = 0; x < clim; x++ )
		{
			int v2 = paths[x];
	
			double dr[3] = { r_sph[3*v2+0] - r_alt[3*v+0], r_sph[3*v2+1] - r_alt[3*v+1], r_sph[3*v2+2] - (r_alt[3*v+2]) };
			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r < bestr2 )
			{
				startsph = x;
				bestr2 = r;
			}
		}

		// which direction to go.
	
		int next3 = startsph+3;
		if( next3 >= clim ) next3 -= clim;
		int prev3 = startsph-3;
		if( prev3 < 0 ) prev3 += clim;

		int v2n,v2p;

		v2n = paths[next3];
		v2p = paths[prev3];

		double rd1[3] = { 
				r_sph[v2n*3+0] - r_alt[3*(alt_next)+0],
				r_sph[v2n*3+1] - r_alt[3*(alt_next)+1],
				r_sph[v2n*3+2] - r_alt[3*(alt_next)+2] };
		
		double rd2[3] = { 
				r_sph[v2p*3+0] - r_alt[3*(alt_next)+0],
				r_sph[v2p*3+1] - r_alt[3*(alt_next)+1],
				r_sph[v2p*3+2] - r_alt[3*(alt_next)+2] };
		
		double rv1 = rd1[0]*rd1[0]+rd1[1]*rd1[1]+rd1[2]*rd1[2];
		double rv2 = rd2[0]*rd2[0]+rd2[1]*rd2[1]+rd2[2]*rd2[2];
	
		if( rv1 < rv2 )
		{
			// forward
		}
		else
		{

			// reverse
			int copy[spl];
			int news = 0;	
			for( int x = 0; x < spl; x++ )
			{
				copy[spl-x-1] = paths[x];
				if( x == startsph )
					news = spl-x-1;
			}
			startsph = news;
			memcpy( paths, copy, sizeof(int) * spl );
		}
	}

	int nvtot = nr_ex + nr_sph_ex;
	double *allr = (double *)malloc( sizeof(double) * 3 * nvtot );
	int    *alle = (int *)malloc( sizeof(int) * nvtot * MAX_EDGES );
	int    *ne = (int *)malloc( sizeof(int) * nvtot );
	memset( ne, 0, sizeof(int) * nvtot );
	int *pl_map = (int *)malloc( sizeof(int) * nr_pl );
	int *sph_map = (int *)malloc( sizeof(int) * nr_sph );

	int ii = 0;

	for( int v = 0; v < nr_pl; v++ )
	{
		pl_map[v] = -1;
		if( interior_tag[v] ) continue;

 		pl_map[v] = ii;
		ii++;
	}
	
	printf("ii: %d nr_ex: %d\n", ii, nr_ex );
	
	ii = 0;

	for( int v = 0; v < nr_sph; v++ )
	{
		sph_map[v] = -1;

		if( sph_interior_tag[v] ) continue;

		sph_map[v] = nr_ex + ii;
		ii++;
	}
	
	printf("ii: %d nr_ex: %d\n", ii, nr_sph_ex );
	
	int *path2 = (int *)malloc( sizeof(int) * pl );
	memcpy( path2, path, sizeof(int) * pl );

	for( int px = 0; px < pl; px++ )
		path[px]  = pl_map[path[px]];
	for( int px = 0; px < spl; px++ )
		paths[px] = sph_map[paths[px]];
	
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
	
	for( int vx = 0; vx < nr_sph; vx++ )
	{
		int v = sph_map[vx];

		if( v >= 0 )
		{	
			allr[3*v+0] = r_sph[3*vx+0];
			allr[3*v+1] = r_sph[3*vx+1];
			allr[3*v+2] = r_sph[3*vx+2];
			
			for( int e = 0; e < ne_sph[vx]; e++ )
			{
				if( sph_map[e_sph[vx*MAX_EDGES+e]] >= 0 )
				{
					alle[v*MAX_EDGES+ne[v]] = sph_map[e_sph[vx*MAX_EDGES+e]];
					ne[v]++; 
				}
			}
		}
	}
	
	// loop over both paths.

	{
		int *path1 = path;
		int *path2 = paths;
		int lim1 = pl;
		int lim2 = spl;
		int pstart1 = 0;
		int pstart2 = startsph;

		if( lim2 > lim1 )
		{
			int t = lim1;
			lim1 = lim2;
			lim2 = t;
			int *tp = path1;
			path1 = path2;
			path2 = tp;

			t = pstart1;
			pstart1 = pstart2;
			pstart2 = t;
		}

		int *mpath = path;
		int mlim = pl;


		// planar path is longer.
		int excess = lim1 - lim2;
	
		double excess_cntr = 0;
	
		int p  = pstart1;
		int p2 = pstart2;

	

		for( int pc = 0; pc < lim2; pc++, p2++, p++ )
		{
			if( p >= lim1 )
				p -= lim1;

			if( p2 >= lim2 )
				p2 -= lim2;
	
			int vc = path2[p2];
			
			int vp = path1[(p-1<0?p-1+lim1:p-1)];
			int v = path1[p];
			int vn = path1[(p+1>=lim1?p+1-lim1:p+1)];
					 
//			if( excess_cntr >= 1.0 && ne[vc] <= 4 && excess > 0 )
			if( excess > 0 && ne[vc] <= 4 )
			{
				alle[vc*MAX_EDGES+ne[vc]]   = vp; ne[vc]++;
				alle[vc*MAX_EDGES+ne[vc]]  =   v; ne[vc]++;
				alle[vc*MAX_EDGES+ne[vc]]  =  vn; ne[vc]++;
		
				alle[vp*MAX_EDGES+ne[vp]] = vc; ne[vp]++;	
				alle[v*MAX_EDGES+ne[v]]   = vc; ne[v]++;	
				alle[vn*MAX_EDGES+ne[vn]]   = vc; ne[vn]++;

				printf("Assigning excess to reach valence %d, cntr = %lf\n", ne[vc], excess_cntr );
				p++;	
				excess_cntr -= 1;
				excess -= 1;
			}
			else
			{
				alle[vc*MAX_EDGES+ne[vc]]   = vp; ne[vc]++;
				alle[vc*MAX_EDGES+ne[vc]]  =   v; ne[vc]++;
		
				alle[vp*MAX_EDGES+ne[vp]] = vc; ne[vp]++;	
				alle[v*MAX_EDGES+ne[v]]   = vc; ne[v]++;	
				
			}

			excess_cntr += (lim1-lim2) / (double)lim2;
		}

		// tack on at the end, hopefully we don't get here.

		p2--;
		if( p2 < 0 ) 
			p2 += lim2;

		if( excess > 0 )
		{
			printf("EXCESS: %d\n", excess );
		}	
	
		for( int pc = 0; pc < excess; pc++ ) // lim2; p < lim1; p++ )
		{
			int vc = path2[p2];		
			int v = path1[p];
	
			printf("vc: %d p2: %d\n", vc, p2 );				 
			alle[vc*MAX_EDGES+ne[vc]]  =    v; ne[vc]++;
			alle[v*MAX_EDGES+ne[v]]    =   vc; ne[v]++;	
			
			p += 1;			
			
			if( p >= lim1 )
				p -= lim1;
		}

		printf("Final vertex has valence %d.\n", ne[path2[p2]] );

		printf("lim1: %d\n", lim1 );
		for( int x = 0; x < lim1; x++ )
		{
			printf("ne: %d\n", ne[path1[x]] );
			if( ne[path1[x]] < 5 || ne[path1[x]] > 7 )
			{
					printf("WARNING: valence %d v %d in path linking plane to sphere. PL %d SPL %d\n", ne[path1[x]], path1[x], pl, spl );
			}
			
		}
		printf("lim2: %d\n", lim2 );
		for( int x = 0; x < lim2; x++ )
		{
			printf("ne: %d %d \n", ne[path2[x]], path2[x]  );
			if( ne[path2[x]] < 5 || ne[path2[x]] > 7 )
			{
					printf("WARNING: valence %d v %d in path linking plane to sphere. PL %d CSPL %d\n", ne[path2[x]], path2[x], pl, spl );
			}
		}
	}

	
	char fileName[256];
	sprintf(fileName,"bump.mesh" );

	FILE *theLattice = fopen(fileName,"w");

	fprintf(theLattice, "3d R = %lf\n", R );
	fprintf(theLattice, "%lf 0.0 0.0\n", Lx );
	fprintf(theLattice, "0.0 %lf 0.0\n", Lx );
	fprintf(theLattice, "0.0 0.0 %lf\n", Lx + bumpR );

	int p = 0;
	for( ; p < nvtot; p++ )
	{
		fprintf(theLattice, "%d %lf %lf %lf %d", p, allr[3*p+0], allr[3*p+1], allr[3*p+2], ne[p] );

		if( ne[p] < 5 || ne[p] > 7 )
		{
			printf("WARNING: unusual valence %d.\n", ne[p] );
		}

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



void getSphereMesh( double **r_io, int **edges_io, int **nedges_io, int *nv_io, int ndiv, double R)
{
	// use external calls for this.

	char command[256];
	sprintf(command, "%s/icosahedron %lf > subdiv.mesh", program_base, R);
	system(command);
	
	for( int x = 0; x < ndiv; x++ )
	{
		char command[256];
		sprintf(command, "%s/subdivide subdiv.mesh", program_base );
		system(command);
	}

	// read the mesh.

	FILE *theFile = fopen("subdiv.mesh", "r" );

	if( !theFile )
	{
		printf("Failed to create the sphere mesh.\n");
		exit(1);
	}

	char buffer[4096];

	getLine( theFile, buffer );
	getLine( theFile, buffer );
	getLine( theFile, buffer );
	getLine( theFile, buffer );

	int nvs = 10;
	int nv = 0;
	double *r = (double *)malloc( sizeof(double ) * 3 * nvs );
	int *edges = (int *)malloc( sizeof(int) * MAX_EDGES * nvs );
	int *nedges = (int *)malloc( sizeof(int) * nvs );
	

	while( !feof(theFile) )	
	{
		getLine( theFile, buffer );

		if( !strncasecmp( buffer, "ntri", 4 ) )
			break;

		if( nvs == nv )
		{
			nvs *= 2;
			r = (double *)realloc( r, sizeof(double ) * 3 * nvs );
			edges = (int *)realloc( edges, sizeof(int) * MAX_EDGES * nvs );
			nedges = (int *)realloc( nedges, sizeof(int) * nvs );
		
		}
		
		double x,y,z;
		int ver;
		int edge_buffer[12];
		int nr = sscanf( buffer, "%d %lf %lf %lf %d "
					 "%d %d %d %d %d %d %d %d %d %d %d %d",
					&ver, r+3*nv, r+3*nv+1, r+3*nv+2,
					nedges+nv,
				 	edge_buffer+0,
				 	edge_buffer+1,
				 	edge_buffer+2,
				 	edge_buffer+3,
				 	edge_buffer+4,
				 	edge_buffer+5,
				 	edge_buffer+6,
				 	edge_buffer+7,
				 	edge_buffer+8,
				 	edge_buffer+9,
				 	edge_buffer+10,
				 	edge_buffer+11 );
		if( nr > 5 + MAX_EDGES )
		{
			printf("MAX_EDGES too small to hold mesh.\n");
			exit(1);
		}

		for( int e = 0; e < nedges[nv]; e++ )
			edges[nv*MAX_EDGES+e] = edge_buffer[e];

		nv += 1;
	}

	// subtract out the c.o.m.

	double com[3] = {0,0,0};
	for( int v = 0; v < nv; v++ )
	{
		com[0] += r[3*v+0];
		com[1] += r[3*v+1];
		com[2] += r[3*v+2];
	}

	for( int v = 0; v < nv; v++ )
	{
		r[3*v+0] -= com[0]/nv;
		r[3*v+1] -= com[1]/nv;
		r[3*v+2] -= com[2]/nv;
	}

	fclose(theFile);
	*r_io = r;
	*edges_io = edges;
	*nedges_io = nedges;
	*nv_io = nv;	
}


















