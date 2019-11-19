#include "simulation.h"
#include "interp.h"
#include "pcomplex.h"
#include <string.h>
#include <math.h>
#include "l-bfgs.h"
#include "p_p.h"
#include "parallel.h"
#include "fitRho.h"

extern int enable_elastic_interior;
static int min_ncomplex = 0;
static int min_nparams = 0;
static int min_nsurfaceparams = 0;
static int do_freeze_membrane = 0;
static pcomplex **min_complexes;
extern double VA,VC;
Simulation *min_simulation = NULL;
static double cur_rho_thickness[2] = { 15.0, 15.0 };

#define FIT_RHO_ONE_THICKNESS

double surface_f( double *p )
{
	int offset = 3; 
	for( surface_record *sRec = min_simulation->allSurfaces; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;
		p[offset+3*nv+0] = p[0];
		p[offset+3*nv+1] = p[1];
		p[offset+3*nv+2] = p[2];
		
		sRec->r = p+offset;		

//		memcpy( sRec->r, p+offset, sizeof(double) * (3*nv+3) );
		offset += sRec->theSurface->nv*3+3;
	}

#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif

	VA = 0;
	VC = 0;
	double v = 0;
	for( surface_record *sRec = min_simulation->allSurfaces; sRec; sRec = sRec->next )
	{
		v += sRec->theSurface->energy( p+sRec->temp_min_offset, NULL );
	}
	ParallelSum(&VA,1);
	ParallelSum(&VC,1);
	int nparams = min_nsurfaceparams;

	for( int c = 0; c < min_ncomplex; c++ )
	{
		min_complexes[c]->applyParamsToComplex( p + nparams );
		nparams += min_complexes[c]->nparams();
	}
		
	if( fitRho_activated )
	{
#ifdef FIT_RHO_ONE_THICKNESS
		double thick_inner = p[nparams]; nparams++;
		double thick_outer = thick_inner;
#else
		double thick_inner = p[nparams]; nparams++;
		double thick_outer = p[nparams]; nparams++;
#endif
		for( surface_record *sRec = min_simulation->allSurfaces; sRec; sRec = sRec->next )
			v += sRec->theSurface->rhoEnergy( p+sRec->temp_min_offset, min_simulation->PBC_vec, thick_inner, thick_outer );
	}	
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];

		v += min_complexes[c]->V( min_simulation );
		v += min_complexes[c]->AttachV( min_simulation );
	}
	v += Boxed_PP_V( min_simulation );
	
	ParallelSum(&v,1);


#ifdef PARALLEL
	MPI_Bcast( &v, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif
//	printf("v: %le\n", v);
	return v;
}



double surface_fdf( double *p, double *g)
{
#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif

	memset( g, 0, sizeof(double) * min_nparams );
	double v = 0;	
	int offset = 3; 
	
	for( surface_record *sRec = min_simulation->allSurfaces; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;
		p[offset+3*nv+0] = p[0];
		p[offset+3*nv+1] = p[1];
		p[offset+3*nv+2] = p[2];
		sRec->g = g+offset;
		sRec->r = p+offset;		

		sRec->theSurface->grad( p+offset, g+offset );
		offset += sRec->theSurface->nv*3+3;
	}
	int nparams = min_nsurfaceparams;	
	for( int c = 0; c < min_ncomplex; c++ )
	{
		int np = min_complexes[c]->nparams();
		min_complexes[c]->applyParamsToComplex( p + nparams );
		nparams += np;
	}
	
	if( fitRho_activated )
	{
#ifdef FIT_RHO_ONE_THICKNESS
		double thick_inner = p[nparams];
		double *rho_g_i = g+nparams;
		*rho_g_i = 0;

		double thick_outer = p[nparams]; 
		double *rho_g_o = g+nparams;
		*rho_g_o = 0;

		nparams++;
#else
		double thick_inner = p[nparams];
		double *rho_g_i = g+nparams;
		*rho_g_i = 0;
		nparams++;

		double thick_outer = p[nparams]; 
		double *rho_g_o = g+nparams;
		*rho_g_o = 0;
		nparams++;
#endif
		int offset = 3; 

		*rho_g_i = 0;
		*rho_g_o = 0;
	
		for( surface_record *sRec = min_simulation->allSurfaces; sRec; sRec = sRec->next )
		{
			int nv = sRec->theSurface->nv;
			p[offset+3*nv+0] = p[0];
			p[offset+3*nv+1] = p[1];
			p[offset+3*nv+2] = p[2];
			sRec->g = g+offset;
			sRec->r = p+offset;		

			v += sRec->theSurface->rhoGrad( p+offset, g+offset, min_simulation->PBC_vec, thick_inner, thick_outer, rho_g_i, rho_g_o );
			offset += sRec->theSurface->nv*3+3;
		}
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
		double le = min_complexes[c]->grad( min_simulation, rg, ng );
			
		v += le;
	
		if( min_complexes[c]->bound )
		{
			// attachment coordinates need gradients.
			
			for( int a = 0; a < min_complexes[c]->nattach; a++ )
			{
				double point_grad[2]={0,0};
				double point_rgrad[3] = { 0,0,0};
				surface_record *sRec = min_simulation->fetch( min_complexes[c]->sid[a] );
				
				// function updates the gradient of the energy for the surface coordinates.

				sRec->theSurface->pointGradient( min_complexes[c]->grad_fs[a], min_complexes[c]->grad_puv[2*a+0], 
										     min_complexes[c]->grad_puv[2*a+1],
										     sRec->r, sRec->g, point_grad, rg+3*a, ng+3*a );  			
				min_complexes[c]->save_grad[2*a+0] += point_grad[0];		
				min_complexes[c]->save_grad[2*a+1] += point_grad[1];		
			}
			
		}
	
		v += min_complexes[c]->AttachG( min_simulation, min_complexes[c]->save_grad );
	}

	nparams = min_nsurfaceparams;

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
		
	v += Boxed_PP_G( min_simulation);
	
	nparams = min_nsurfaceparams;
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

				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+0] * M[0];
				g[nparams+2*a+0] += min_complexes[c]->save_grad[2*a+1] * M[2];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+0] * M[1];
				g[nparams+2*a+1] += min_complexes[c]->save_grad[2*a+1] * M[3];
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

	// for now zero box dimension gradient.
	
	g[0] = 0;
	g[1] = 0;
	g[2] = 0;

	for( surface_record *sRec = min_simulation->allSurfaces; sRec; sRec = sRec->next )
	{
		g[sRec->temp_min_offset+3*sRec->theSurface->nv+0] = 0;
		g[sRec->temp_min_offset+3*sRec->theSurface->nv+1] = 0;
		g[sRec->temp_min_offset+3*sRec->theSurface->nv+2] = 0;
	}
	
	if( do_freeze_membrane )
	{
		for( int x = 0; x < min_nsurfaceparams; x++ )
			g[x] = 0;
	}

	for( surface_record *sRec = min_simulation->allSurfaces; sRec; sRec = sRec->next )
		v += sRec->theSurface->energy( p+sRec->temp_min_offset, NULL );
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

}

void full_fd_test( double *p )
{
	double deps = 1e-4;

	double *tp = (double *)malloc( sizeof(double) * min_nparams );
	memset( tp, 0, sizeof(double) * min_nparams );
	double *g = (double *)malloc( sizeof(double) * min_nparams );
	memset( g, 0, sizeof(double) * min_nparams );

	double e0 = surface_fdf( p, g );
	printf("Finite difference test.\n");


	memcpy( tp, p, sizeof(double) * min_nparams );
	for( int p = 0; p < min_nparams; p++ )
	{
		double de_pm[2];
	
		for( int im = 0; im < 2; im++ )
		{
			tp[p] += deps * (im == 0 ? 1 : -1);	
		
			double v = surface_f( tp );
			
			de_pm[im] = v;

			tp[p] -= deps * (im == 0 ? 1 : -1);	
		}
		
		printf("parm %d fd_der %.14le g %.14le del %.14le\n", p, (de_pm[0]-de_pm[1])/(2*deps), g[p],  (de_pm[0]-de_pm[1])/(2*deps) - g[p] );
	}
	free(tp);
	free(g);
}

void fd_test( double *p )
{
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
		printf("%le %.14le %.14le del %.14le\n", eps, v-e0, expec, (v-(e0+expec)) );
	}

	free(tp);
	free(g);
}

void Simulation::minimize( int freeze_membrane  )
{
	do_freeze_membrane = freeze_membrane;

	int prev_enable = enable_elastic_interior; 

	enable_elastic_interior = 1;


	min_simulation = this;
	min_ncomplex = ncomplex;	
	min_complexes = allComplexes;

#ifdef PARALLEL
	ParallelSyncComplexes( min_complexes, min_ncomplex );
#endif
	
	int num_params = 3; // alphas

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		num_params += sRec->theSurface->nv*3+3;

	min_nsurfaceparams = num_params;

	for( int c = 0; c < ncomplex; c++ )
		num_params += allComplexes[c]->nparams();

	if( fitRho_activated )
	{
#ifdef FIT_RHO_ONE_THICKNESS
		num_params += 1; // same thickness
#else
		num_params += 2; // inner and outer thickness.
#endif
	}
	min_nparams = num_params;
	
	double *p = (double *)malloc( sizeof(double) * num_params );
	double *g = (double *)malloc( sizeof(double) * num_params );

	p[0] = alpha[0];
	p[1] = alpha[1];
	p[2] = alpha[2];
	int offset = 3;
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;

		memcpy( p+offset, sRec->r, sizeof(double) * 3 * nv );
		sRec->temp_min_offset = offset;
		sRec->temp_r = sRec->r;
		sRec->temp_g = sRec->g;
		sRec->r = p + offset;
		sRec->g = g + offset;

		p[offset+3*nv+0] = alpha[0];
		p[offset+3*nv+1] = alpha[1];
		p[offset+3*nv+2] = alpha[2];

		offset += 3*nv+3;
	}

	int tp = min_nsurfaceparams;

	for( int c = 0; c < ncomplex; c++ )
	{
		allComplexes[c]->getParamsFromComplex( p + tp );
		tp += allComplexes[c]->nparams();
	}	

	int thickness_ptr = 0;

	if( fitRho_activated )
	{
		printf("Current thickness: %lf %lf\n", cur_rho_thickness[0], cur_rho_thickness[1] );
		thickness_ptr = tp;

#ifdef FIT_RHO_ONE_THICKNESS
		p[tp] =cur_rho_thickness[0]; tp++; 
#else
		p[tp] =cur_rho_thickness[0]; tp++; 
		p[tp] =cur_rho_thickness[1]; tp++; 
#endif
	}
	
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

	full_fd_test(p);
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

	double v =surface_fdf(p,g);
	double e = surface_f(p);

	double rms = 0;
	for( int x = 0; x < num_params; x++ )
		rms += g[x]*g[x];
	rms /= num_params;
	rms = sqrt(rms);

	printf("Minimize: VG: %.14le VV: %.14le VA: %lf VC: %lf grad rms %le\n", v, e, VA, VC, rms );

	enable_elastic_interior = prev_enable;

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		int nv = sRec->theSurface->nv;

		memcpy( sRec->temp_r, sRec->r, sizeof(double) * 3 * nv );
		sRec->temp_r[3*nv+0] = p[0];
		sRec->temp_r[3*nv+1] = p[1];
		sRec->temp_r[3*nv+2] = p[2];
		sRec->r = sRec->temp_r;
		sRec->g = sRec->temp_g;
	}
	
	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->refresh(this);

	if( fitRho_activated )
	{
#ifdef FIT_RHO_ONE_THICKNESS
		cur_rho_thickness[0] = p[thickness_ptr];		
		cur_rho_thickness[1] = p[thickness_ptr];		
#else
		cur_rho_thickness[0] = p[thickness_ptr];		
		cur_rho_thickness[1] = p[thickness_ptr+1];		
#endif
	}
	
	free(p);
	free(g);
}

