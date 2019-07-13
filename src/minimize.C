#include "interp.h"
#include "pcomplex.h"
#include <string.h>
#include <math.h>
#include "l-bfgs.h"
#include "p_p.h"
#include "parallel.h"

extern int enable_elastic_interior;
static int min_ncomplex = 0;
static int min_nparams = 0;
static int do_freeze_membrane = 0;
static pcomplex **min_complexes;
static surface *min_surface;
extern double VA,VC;
double surface_f( double *p )
{
#if 0 // pre-simulation change me
#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif

	VA = 0;
	VC = 0;
	double v = min_surface->energy( p, NULL );
	ParallelSum(&VA,1);
	ParallelSum(&VC,1);
	int nparams = 3 * min_surface->nv + 3;

	for( int c = 0; c < min_ncomplex; c++ )
	{
		min_complexes[c]->applyParamsToComplex( p + nparams );
		nparams += min_complexes[c]->nparams();
	}
		
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];

		v += min_complexes[c]->V( min_surface, p );
		v += min_complexes[c]->AttachV( min_surface, p );
	}
	v += Boxed_PP_V( min_surface, p, min_complexes, min_ncomplex );
	
	ParallelSum(&v,1);


#ifdef PARALLEL
	MPI_Bcast( &v, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif
//	printf("v: %le\n", v);
	return v;
#endif
}



double surface_fdf( double *p, double *g)
{
#if 0 // pre-simulation change me
#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif
	int nparams = 3 * min_surface->nv + 3;

	memset( g, 0, sizeof(double) * min_nparams );
	double v = 0;	
	min_surface->grad( p, g );
	
	for( int c = 0; c < min_ncomplex; c++ )
	{
		int np = min_complexes[c]->nparams();
		min_complexes[c]->applyParamsToComplex( p + nparams );
		nparams += np;
	}
	
	for( int c = 0; c < min_ncomplex; c++ )
		memset( min_complexes[c]->save_grad, 0, sizeof(double) * 3 * min_complexes[c]->nsites );

#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];


		int nsites = min_complexes[c]->nsites;

		double rg[3*nsites];
		double ng[3*nsites];

		memset( rg, 0, sizeof(double) * 3 * nsites );
		memset( ng, 0, sizeof(double) * 3 * nsites );
		// derivative wrt positions, normals
		double le = min_complexes[c]->grad( min_surface, p, rg, ng );
			
		v += le;
	
		if( min_complexes[c]->bound )
		{
			// attachment coordinates need gradients.
			
			for( int a = 0; a < min_complexes[c]->nattach; a++ )
			{
				double point_grad[2]={0,0};
				double point_rgrad[3] = { 0,0,0};

				// function updates the gradient of the energy for the surface coordinates.

				min_surface->pointGradient( min_complexes[c]->grad_fs[a], min_complexes[c]->grad_puv[2*a+0], 
										     min_complexes[c]->grad_puv[2*a+1],
										     p, g, point_grad, rg+3*a, ng+3*a );  			
				min_complexes[c]->save_grad[2*a+0] += point_grad[0];		
				min_complexes[c]->save_grad[2*a+1] += point_grad[1];		

				if( !(min_complexes[c]->save_grad[2*a+0] < 0) && !(min_complexes[c]->save_grad[2*a+0] > -1) )
				{
					printf("P GNAN.\n");
					exit(1);
				}
				if( !(min_complexes[c]->save_grad[2*a+1] < 0) && !(min_complexes[c]->save_grad[2*a+1] > -1) )
				{
					printf("P GNAN.\n");
					exit(1);
				}
			}
			
		}
	
		v += min_complexes[c]->AttachG( min_surface, p, g, min_complexes[c]->save_grad );
	}


	nparams = 3*min_surface->nv+3;

	for( int c = 0; c < min_ncomplex; c++ )
	{
		int np = min_complexes[c]->nparams();
		if( min_complexes[c]->bound )
		{
			// attachment coordinates need gradients.
		
			for( int a = 0; a < min_complexes[c]->nattach; a++ )
			{
				int np = min_complexes[c]->nparams();

				double M[4] = { min_complexes[c]->coord_transform[4*a+0],
						min_complexes[c]->coord_transform[4*a+1],
						min_complexes[c]->coord_transform[4*a+2],
						min_complexes[c]->coord_transform[4*a+3] };
				double det = M[0]*M[3]-M[1]*M[2];
				double MINV[4] = { M[3]/det, -M[1]/det, -M[2]/det, M[0]/det };

#if 1 
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0] * M[0];
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+1] * M[2];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+0] * M[1];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1] * M[3];
#else
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1];

#endif

			}
	 	
			int poff = 2 * min_complexes[c]->nattach;
			for( int a = min_complexes[c]->nattach; a < min_complexes[c]->nsites; a++ )
			{
				int ax = a - min_complexes[c]->nattach;
				g[nparams+poff+3*ax+0] += min_complexes[c]->save_grad[3*a+0];
				g[nparams+poff+3*ax+1] += min_complexes[c]->save_grad[3*a+1];
				g[nparams+poff+3*ax+2] += min_complexes[c]->save_grad[3*a+2];
			}			
		}

		nparams += np;
	}	
	
	for( int c = 0; c < min_ncomplex; c++ )
		memset( min_complexes[c]->save_grad, 0, sizeof(double) * 3 * min_complexes[c]->nsites );
		
	v += Boxed_PP_G( min_surface, p, min_complexes, min_ncomplex, g );
	
	nparams = 3 * min_surface->nv + 3;
	for( int c = 0; c < min_ncomplex; c++ )
	{
		int np = min_complexes[c]->nparams();

		for( int a = 0; a < min_complexes[c]->nattach;  a++ )
		{
				double M[4] = { min_complexes[c]->coord_transform[4*a+0],
						min_complexes[c]->coord_transform[4*a+1],
						min_complexes[c]->coord_transform[4*a+2],
						min_complexes[c]->coord_transform[4*a+3] };
				double det = M[0]*M[3]-M[1]*M[2];
				double MINV[4] = { M[3]/det, -M[1]/det, -M[2]/det, M[0]/det };

#if 1 
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0] * M[0];
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+1] * M[2];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+0] * M[1];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1] * M[3];
#else
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1];

#endif
		}
	 	int poff = 2 * min_complexes[c]->nattach;
		for( int a = min_complexes[c]->nattach; a < min_complexes[c]->nsites; a++ )
		{
			int ax = a - min_complexes[c]->nattach;
			g[nparams+poff+3*ax+0] += min_complexes[c]->save_grad[3*a+0];
			g[nparams+poff+3*ax+1] += min_complexes[c]->save_grad[3*a+1];
			g[nparams+poff+3*ax+2] += min_complexes[c]->save_grad[3*a+2];
		}
		nparams += np;
	}
	
	g[3*min_surface->nv+0] = 0;
	g[3*min_surface->nv+1] = 0;
	g[3*min_surface->nv+2] = 0;
	
	if( do_freeze_membrane )
	{
		for( int x = 0; x < 3*min_surface->nv; x++ )
			g[x] = 0;
	}

	v += min_surface->energy( p, NULL );
	ParallelSum( &v, 1 );

	ParallelSum( g, min_nparams );

#if 0 	
	MPI_Barrier(MPI_COMM_WORLD);
	if( par_info.my_id == BASE_TASK )
	{	
		FILE *gradf = fopen("grad.txt","w");
		for( int c = 0; c < min_nparams; c++ )
			fprintf(gradf, "%d %.14le\n", c, g[c] );
		fclose(gradf);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	exit(1);
#endif
	return v;
#endif

}

void full_fd_test( double *p )
{
#if 0 // pre-simulation change me
	double deps = 1e-10;

	double *tp = (double *)malloc( sizeof(double) * min_nparams );
	memset( tp, 0, sizeof(double) * min_nparams );
	double *g = (double *)malloc( sizeof(double) * min_nparams );
	memset( g, 0, sizeof(double) * min_nparams );

	double e0 = surface_fdf( p, g );
	printf("Finite difference test.\n");
	
	memcpy( tp, p, sizeof(double) * min_nparams );
	for( int ttp = 0; ttp < min_nparams; ttp++ )
	{
		double use_eps[2] = { -1e-6, 1e-6 };
		double ens[2];
		for( int ieps = 0; ieps < 2; ieps++ )
		{
			tp[ttp] += use_eps[ieps];
			ens[ieps] = surface_f(tp);
			tp[ttp] -= use_eps[ieps];
		}

		double de = (ens[1]-ens[0])/(use_eps[1]-use_eps[0]);
	
		if( ttp == (3 * min_surface->nv+3) )
			printf("FINISHED MESH FD\n");

		if( fabs(de-g[ttp])/(fabs(1e-20)+fabs(de)+fabs(g[ttp])) < 1e-3 )
			printf("P %d FD %.14le G %.14le OK\n", ttp, de, g[ttp] );
		else
			printf("P %d FD %.14le G %.14le CHECK\n", ttp, de, g[ttp] );
	}

	for( double eps = 0; eps < 20 * deps; eps += deps)
	{
		double expec = 0;
		for( int x = 0; x < min_nparams;x++ )
		{
			tp[x] = p[x] - eps * g[x];
			expec += -g[x] * g[x] * eps;
		}

		double v = surface_f( tp );
		printf("%le %.14le %.14le del %.14le\n", eps, v, e0+expec, (v-(e0+expec)) );
	}

	free(tp);
	free(g);
#endif
}

void fd_test( double *p )
{
#if 0 // pre-simulation change me
	double deps = 1e-10;

	double *tp = (double *)malloc( sizeof(double) * min_nparams );
	memset( tp, 0, sizeof(double) * min_nparams );
	double *g = (double *)malloc( sizeof(double) * min_nparams );
	memset( g, 0, sizeof(double) * min_nparams );

	double e0 = surface_fdf( p, g );
	printf("Finite difference test.\n");
	for( double eps = 0; eps < 20 * deps; eps += deps)
	{
		double expec = 0;
		for( int x = 0; x < min_nparams;x++ )
		{
			tp[x] = p[x] - eps * g[x];
			expec += -g[x] * g[x] * eps;
		}

		double v = surface_f( tp );
		printf("%le %.14le %.14le del %.14le\n", eps, v, e0+expec, (v-(e0+expec)) );
	}

	free(tp);
	free(g);
#endif
}

void surface::minimize( double *r, pcomplex **allComplexes, int ncomplex, int freeze_membrane )
{
#if 0 // pre-simulation change me
	do_freeze_membrane = freeze_membrane;

	int prev_enable = enable_elastic_interior; 

	enable_elastic_interior = 1;

	min_surface = this;
	min_ncomplex = ncomplex;	
	min_complexes = allComplexes;



#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif
	
	int num_params = 3*nv+3;

	for( int c = 0; c < ncomplex; c++ )
		num_params += allComplexes[c]->nparams();
	min_nparams = num_params;
	

	double *p = (double *)malloc( sizeof(double) * num_params );
	memcpy( p, r, sizeof(double) * 3 * (nv+1) );

	int tp = 3*nv+3;

	for( int c = 0; c < ncomplex; c++ )
	{
		allComplexes[c]->getParamsFromComplex( p + tp );
		tp += allComplexes[c]->nparams();
	}	
	
	double *g = (double *)malloc( sizeof(double) * num_params );
	// derivative might be zero (absolutely)
	surface_fdf(p,g);
	double mag_init = 0;
	for( int p = 0; p < num_params; p++ )
		mag_init += g[p]*g[p];
	

	int nsteps = 100;
	
	int use_m = nsteps;
	if( use_m > num_params )
		use_m = num_params;
	
	double e_init = surface_f(p);
	printf("Entering minimize with e_init: %le\n", e_init );
	l_bfgs_setup( use_m, num_params, p, 1.0, surface_f, surface_fdf); 

	if( mag_init > 1e-20 )
	{
		for( int x = 0; x < nsteps; x++ )
		{
			if( ! l_bfgs_iteration( p ) ) { break; }
	
//		if( x %10 == 0 )
//			printf("Sub iteration %d, V: %le\n", x, surface_f(p) );
		}
	}	
	else
	{
		printf("Initial gradient zero.\n");
	}	
	l_bfgs_clear();
//	full_fd_test(p);

	double v =surface_fdf(p,g);
	double e = surface_f(p);

	double rms = 0;
	for( int x = 0; x < num_params; x++ )
		rms += g[x]*g[x];
	rms /= num_params;
	rms = sqrt(rms);

	printf("Minimize: VG: %.14le VV: %.14le VA: %lf VC: %lf grad rms %le\n", v, e, VA, VC, rms );

	memcpy( r, p, sizeof(double) * 3 * (nv+1) );
	
	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->refresh(this, r);
	
	free(p);
	free(g);
	
	enable_elastic_interior = prev_enable;
#endif
}




