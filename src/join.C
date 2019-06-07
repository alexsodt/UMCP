//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "l-bfgs.h"
#include "m_triangles.h"
#include "parallel.h"
#include "mutil.h"
//#define SKIP_OPT

#define FIX_VOLUME

surface *minSurface = NULL;
double thick = 13.8;
extern int do_print_area_ps;
extern double av_vol_err;
extern double KV_scale;
extern double KA ;
void fdiff_check( double *coords );
int print_energy = 0;
extern double tilt_scale;
surface *upperSurface;
surface *lowerSurface;
volumes *theVolumes; 
 void amoeba(double **p, double y[], int ndim, double ftol,
               double (*funk)(double []), int *nfunk);
void amoebaOpt( double *r, int nv, double *rc, int nrc );
double surface_grad( double *r, double *g );
double surface_e( double *r );
	void curvatureGrad( double *g, surface *surface1, surface *surface2, double *r, double thick );
	double curvatureEnergy( surface *surface1, surface *surface2, double *r, double thick);
	void tiltGrad( double *g, surface *surface1, surface *surface2, double *r );
	void tiltGradExt( double *g, surface *surface1, surface *surface2, double *r );
	double tiltEnergy( surface *surface1, surface *surface2, double *r);
	double tiltEnergyExt( surface *surface1, surface *surface2, double *r, double *cons_energy);
double f(double *r);
double fdf( double *r, double *g);
void writeDualSurface( const char *fileNameXYZ, const char *fileNamePSF, surface *upperSurface, surface *lowerSurface, double *coords, double alpha );

int fix_alpha = 0;

int main( int argc, char **argv )
{
#ifdef PARALLEL
	int ierr;
	ierr = MPI_Init(&argc,&argv);
	if( ierr ) { exit(1); }
	quietParallel();
#endif
#ifdef FIXED_SEED
	srand(1);
#else
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);
#endif

	char buffer[4096];

	if( argc < 7 )
	{
		printf("Syntax: join mesh1 vert1 mesh2 vert2 radius cyl_len\n");
		return -1;
	}
	surface *theSurface1 =(surface *)malloc( sizeof(surface) );
	theSurface1->loadLattice( argv[1], 0. );
	surface *theSurface2 =(surface *)malloc( sizeof(surface) );
	theSurface2->loadLattice( argv[3], 0. );

	theSurface1->generatePlan();
	theSurface2->generatePlan();
	
	double *r1 = (double *)malloc( sizeof(double) * (3 * theSurface1->nv+3) );
	theSurface2->setg0(r1);
	theSurface1->get(r1);
	r1[theSurface1->nv*3] = 1.0;
	r1[theSurface1->nv*3+1] = 1.0;
	r1[theSurface1->nv*3+2] = 1.0;
	
	double *r2 = (double *)malloc( sizeof(double) * (3 * theSurface1->nv+3) );
	theSurface2->setg0(r2);
	theSurface2->get(r2);
	r2[theSurface2->nv*3] = 1.0;
	r2[theSurface2->nv*3+1] = 1.0;
	r2[theSurface2->nv*3+2] = 1.0;

	KA = 1.215;

	int vert1 = atoi(argv[2]);

	if( vert1 >= theSurface1->nv || vert1 < 0 )
	{
		printf("Invalid index for vert1.\n");
		exit(1);

	}
	int vert2 = atoi(argv[4]);

	double join_radius = atof(argv[5]);
	double join_len	   = atof(argv[6]);

	int *path1;
	int pathLen1;

	void getVertPath( surface *theSurface, int vert, double join_radius, int **path_out, int *pathLen_out, double *r);

	getVertPath( theSurface1, vert1, join_radius, &path1, &pathLen1, r1); 

	// now we need a path with the same length, approx, on the opposite mesh.

	int done = 0;

	double scale1 = 0.5;
	double scale2 = 1.5;
		
	int *path2;
	int pathLen2;

	r2[3*theSurface2->nv+0] = scale1;
	r2[3*theSurface2->nv+1] = scale1;
	r2[3*theSurface2->nv+2] = scale1;

	printf("Target is %d.\n", pathLen1 );

	getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
	printf("scale1 path is %d.\n", pathLen2 );

	while( pathLen2 < pathLen1 )
	{
		free(path2);
		scale1 *= 0.95;
		r2[3*theSurface2->nv+0] = scale1;
		r2[3*theSurface2->nv+1] = scale1;
		r2[3*theSurface2->nv+2] = scale1;

		getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
		printf("new scale1 path is %d.\n", pathLen2 );
	}

	free(path2);

	r2[3*theSurface2->nv+0] = scale2;
	r2[3*theSurface2->nv+1] = scale2;
	r2[3*theSurface2->nv+2] = scale2;
	getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
		printf("scale2 path is %d.\n", pathLen2 );

	while( pathLen2 > pathLen1 )
	{
		free(path2);
		scale2 *= 1.05;
		r2[3*theSurface2->nv+0] = scale2;
		r2[3*theSurface2->nv+1] = scale2;
		r2[3*theSurface2->nv+2] = scale2;

		getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );
		printf("new scale2 path is %d.\n", pathLen2 );
	}

	free(path2);

	double trial = (scale1+scale2)/2;
	int niter=0;
	while( !done && niter < 100 )
	{
		r2[3*theSurface2->nv+0] = trial;
		r2[3*theSurface2->nv+1] = trial;
		r2[3*theSurface2->nv+2] = trial;

		getVertPath( theSurface2, vert2, join_radius, &path2, &pathLen2, r2 );

		if( pathLen2 == pathLen1 )
			done = 1;
		else
		{
			free(path2);
				
			printf("trial: %lf scale1: %lf scale2: %lf\n", trial, scale1, scale2  );
			if( pathLen2 > pathLen1 )
				scale1 = trial;
			else
				scale2 = trial;

			trial = (scale1+scale2)/2;
		}	

		niter++;
	}	

	

	return 0;
}
	

void getVertPath( surface *theSurface, int vert, double join_radius, int **path_out, int *pathLen_out, double *r)
{
	// Find an ''edge path'' around object one that matches the radius.
	
	int n_do = 1000;

	int vert_path[n_do];

	double uvcen[2] = { 1.0/3, 1.0/3 };

	double cd1[2], cd2[2], c1,c2;
		
	int f = theSurface->theVertices[vert].faces[0];

	theSurface->c( f, uvcen[0], uvcen[1], r, cd1, cd2, &c1, &c2 );
		
	double drdu[3];
	double drdv[3];

	theSurface->ru( f, uvcen[0], uvcen[1], r, drdu );	
	theSurface->rv( f, uvcen[0], uvcen[1], r, drdv );	

	double vec1[3] = { 
		drdu[0] * cd1[0] + drdv[0] * cd1[1],
		drdu[1] * cd1[0] + drdv[1] * cd1[1],
		drdu[2] * cd1[0] + drdv[2] * cd1[1] };
	
	double vec2[3] = { 
		drdu[0] * cd2[0] + drdv[0] * cd2[1],
		drdu[1] * cd2[0] + drdv[1] * cd2[1],
		drdu[2] * cd2[0] + drdv[2] * cd2[1] };
	double lv1 = normalize(vec1);
	double lv2 = normalize(vec2);

//	printf("%d\n", 2*n_do);
//	printf("circular path\n");

	for( int i = 0; i < n_do; i++ )
	{
		double dx = cos( 2 * M_PI * i /(double)n_do );
		double dy = sin( 2 * M_PI * i /(double)n_do );

		double du = cd1[0] * join_radius * dx/lv1 + cd2[0] * join_radius * dy/lv2;
		double dv = cd1[1] * join_radius * dx/lv1 + cd2[1] * join_radius * dy/lv2;
		
		double u_cen = uvcen[0];
		double v_cen = uvcen[1];

		int f2=f, f_use;
		do 
		{	
			f_use = f2;

			f2 = theSurface->nextFace( f_use, &u_cen, &v_cen, &du, &dv, r );
		} while( f2 != f_use );
		
		double w = 1 - u_cen - v_cen;

		int *cp_use;
		if( f2 < theSurface->nf_faces )
			cp_use = theSurface->theFormulas[f2*theSurface->nf_g_q_p].cp;
		else
			cp_use = theSurface->theIrregularFormulas[(f2-theSurface->nf_faces)*theSurface->nf_irr_pts].cp;

		if( w > u_cen && w > v_cen) // u = small, v = small
			vert_path[i] = cp_use[0];
		else if( u_cen > v_cen ) // u large
			vert_path[i] = cp_use[2];
		else 
			vert_path[i] = cp_use[1];


		double r_eval[3], n_eval[3];

		theSurface->evaluateRNRM( f2, u_cen, v_cen, r_eval, n_eval, r );

//		printf("C %lf %lf %lf # vert %d\n", r[3*vert_path[i]+0], r[3*vert_path[i]+1], r[3*vert_path[i]+2], vert_path[i] );
//		printf("O %lf %lf %lf\n", r_eval[0], r_eval[1], r_eval[2] ); 
	}
	
	int mid_path[n_do];

	mid_path[0] = vert_path[0];
	int nmid = 1;
	int back_trip = 0;
	for( int i = 1; i < n_do; i++ )
	{
		if( vert_path[i] == mid_path[nmid-1] )
			continue; // same as previous point, move to next.

		if( vert_path[i] != mid_path[0] ) back_trip = 1;

		if( back_trip && vert_path[i] == mid_path[0] )
			break;

		mid_path[nmid] = vert_path[i];
		nmid++;	
	}	

	printf("Mid-path, %d points.\n", nmid);	
	
	int final_path[nmid];
	int nfinal = 1;

	final_path[0] = mid_path[0];

	for( int i = 1; i < nmid; i++ )
	{
		// do we add i?

		if( i < nmid-1)
		{
			int next_p = mid_path[i+1];
			int here = final_path[nfinal-1];
			int skip = 0;
			for( int t = 0; t < theSurface->theVertices[here].valence; t++ )
			{
				if( theSurface->theVertices[here].edges[t] == next_p )
					skip = 1;	
			}
		
			if( skip ) continue;
		}
		
		final_path[nfinal] = mid_path[i];
		nfinal++;
	}

	*path_out = (int*)malloc( sizeof(int) * nfinal );
	memcpy( *path_out, final_path, sizeof(int) * nfinal );
	*pathLen_out = nfinal;

	printf("Pathlen: %d\n", nfinal );
}




