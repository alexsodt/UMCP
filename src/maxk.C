#include "interp.h"
#include "uv_map.h"
#include "l-bfgs.h"
#include <math.h>

// a simple utility to find saddle curvature on a surface.

// these are the initial points from which we diverge.
static int src_f, out_f;
static double src_uv[2], out_uv[2];
static double *rsurf = NULL;
static surface *minSurface = NULL;
static double irreg_penalty = 1e5;
double maxk_f( double *parms )
{
	int curf=src_f;

	double uv1[2] = { src_uv[0], src_uv[1] };
	double duv1[2] = { parms[0], parms[1] };
	int f_1,nf=src_f;

	do {
		f_1 = nf;
		nf = minSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf  ); 
	} while( nf != f_1 );

	double k;
	double cvec1[2],cvec2[2],c1,c2;
//	printf("f Eval at %d %.14le %.14le\n", f_1, uv1[0], uv1[1] );
	minSurface->c( f_1, uv1[0], uv1[1], rsurf,  &k, cvec1, cvec2, &c1, &c2 );

	out_f = f_1;
	out_uv[0] = uv1[0];
	out_uv[1] = uv1[1];


	double pot = c1*c2;

	if( f_1 >= minSurface->nf_faces )
	{
		pot += irreg_penalty * pow(1-uv1[0]-uv1[1],2.0);	
	} 

//	printf("f: %.14le --- %d nf_reg: %d %.14le %.14le\n", c1*c2, f_1, minSurface->nf_faces, uv1[0], uv1[1] );
	return pot;	
}

double maxk_fdf( double *parms, double *parms_g )
{
	int curf=src_f;

	double uv1[2] = { src_uv[0], src_uv[1] };
	double duv1[2] = { parms[0], parms[1] };
	double M[4] = { 1,0,0,1};
	int f_1,nf=src_f;
	double null_mom[2]={0,0};
	do {
		f_1 = nf;
		nf = minSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf, null_mom, M ); 
	} while( nf != f_1 );

	double k;

	double cvec1[2],cvec2[2],c1,c2;
	minSurface->c( f_1, uv1[0], uv1[1], rsurf,  &k, cvec1, cvec2, &c1, &c2 );

	double d_c1_duv[2];
	double d_c2_duv[2];

	// these are in terms of the local coordinate system. ``M'' transforms between the original and the local.

//	printf("fdf Eval at %d %.14le %.14le\n", f_1, uv1[0], uv1[1] );
	minSurface->d_c_duv( rsurf, f_1, uv1[0], uv1[1], d_c1_duv, d_c2_duv );

	double pot = c1*c2;

	parms_g[0] = 0;
	parms_g[1] = 0;

	parms_g[0] += M[0*2+0] * d_c1_duv[0] * c2;	
	parms_g[0] += M[1*2+0] * d_c1_duv[1] * c2;	
	
	parms_g[0] += M[0*2+0] * d_c2_duv[0] * c1;	
	parms_g[0] += M[1*2+0] * d_c2_duv[1] * c1;	
	
	parms_g[1] += M[0*2+1] * d_c1_duv[0] * c2;	
	parms_g[1] += M[1*2+1] * d_c1_duv[1] * c2;	
	
	parms_g[1] += M[0*2+1] * d_c2_duv[0] * c1;	
	parms_g[1] += M[1*2+1] * d_c2_duv[1] * c1;	
	
//	printf("fdf: %.14le --- %d nf_reg: %d %.14le %.14le %.14le %.14le\n", c1*c2, f_1, minSurface->nf_faces, uv1[0], uv1[1], parms_g[0], parms_g[1] );

	if( f_1 >= minSurface->nf_faces )
	{
		pot += irreg_penalty * pow(1-uv1[0]-uv1[1],2.0);	

		parms_g[0] -= irreg_penalty * ( 1-uv1[0]-uv1[1]) * M[0*2+0];
		parms_g[0] -= irreg_penalty * ( 1-uv1[0]-uv1[1]) * M[1*2+0];
		parms_g[1] -= irreg_penalty * ( 1-uv1[0]-uv1[1]) * M[0*2+1];
		parms_g[1] -= irreg_penalty * ( 1-uv1[0]-uv1[1]) * M[1*2+1];
	} 

	return pot;	
} 

void max_gauss_c( surface *theSurface, int *f, double *u, double *v, int niter, double *rsurf_in)
{
	disable_random_uv_step();
	minSurface = theSurface;
	rsurf = rsurf_in;

	double *duv = (double *)malloc( sizeof(double) * 2 );
	
	double outk = 0;

	for( int it = 0; it < niter; it++ )
	{
		// where we evaluate the gradient.

		src_f = *f;
		src_uv[0] = *u;
		src_uv[1] = *v;
		duv[0] = 0;
		duv[1] = 0;


		l_bfgs_setup( 4, 2, duv, 0.1, maxk_f, maxk_fdf); 
		
		int nsteps = 200;

		for( int x = 0; x < nsteps; x++ )
		{
			if( ! l_bfgs_iteration( duv ) ) { break; }
		}

		outk = maxk_f( duv );
				

		// do a line minimization just to escape very small wells that the BFGS routines aren't made for.

		double g[2];
		double k0 = maxk_fdf( duv, g );
		

		double deps = 1;	

		double best_f = k0;
		
		*f = out_f;
		*u = out_uv[0];				
		*v = out_uv[1];			
//		printf("FEPS %d %d\n", out_f, minSurface->nf_faces );	

		for( int ieps = 0; ieps < 100; ieps++ )
		{
		 	double trial[2] = {  duv[0] + deps * ieps * g[0],
					     duv[1] + deps * ieps * g[1] };
			double keval = maxk_f(trial);
//			printf("eps: %le eval: %le expec: %le\n", ieps * deps,
//				keval-k0, - deps * ieps * pow(g[0],2.0) - deps * ieps * pow(g[1],2.0)  );  

			if( keval < best_f )
			{
				best_f = keval;
				*f = out_f;
				*u = out_uv[0];				
				*v = out_uv[1];				
			}
		}


		l_bfgs_clear();
	}

	free(duv);

	printf("Final k : %le f: %d uv: %le %le \n", outk, *f, *u, *v );
	enable_random_uv_step();
}
