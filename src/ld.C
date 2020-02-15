#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "mutil.h"
#include "interp.h"
#include "l-bfgs.h"
#include "gsl/gsl_randist.h"
#include "m_triangles.h"
#include "input.h"
#include <assert.h>
#include "srd.h"
#include "units.h"

//#define DEBUG_MOMENTUM_CHANGE
//#define DISABLE_MEMBRANE_MOVE
#define NO_DIFF_R

#ifdef FFTW
#include <fftw3.h>
#endif
int debug_trigger = 0; // wait until all particles have been put in the sim before triggering this check
//#define DEBUG_MODE_ENERGY

#define DIRECT_FT // defined to compute the fourier transform "directly"
//#define TFILE // defined to write  a trajectory

extern double kc;
extern double kg;
extern double KA;
extern double KA4;

double mode_KA = 0; // a ``global'' KA to use when doing modes, because with "z" movements individual elements cannot relax -- this is artificial. 
double dist_nrm = 0;
int do_2p= 0;
surface *sub_surface = NULL;

double pp_grad( double *r, double *gc, int np, double *alpha, double *rads, int * pleaflet);
double particle_particle_energy( double *r, int np_1, int np_2, int *set, int np_set, double rad1, double rad2, double *alpha, double *rads, int *dimer_p, int *pleaflet );
void updateParticleR( int p, int *pfaces, double *puv, double *p_r_m, double *r, surface *sub_surface, int np );

int *global_plist = NULL;

extern double VA, VC, VA4;
double vext = 0;
double max_rad = 0;
int *flag = NULL;

double particle_radius = 0;
double dimerization_radius = 0; // set to zero: never form dimers.
double dimer_eps = -5; // kcal/mol
double search_rad = 0;	
double b;
#define N_BINS 100
#define MAX_R  8.0
double ppdist[N_BINS];
int main( int argc, const char **argv )
{

	char buffer[4096];

	printf("No options are required:\n");
	printf("Syntax: ld [input file] [options ...]\n");
	printf("See input.C for options.\n");

	parameterBlock block;

	int nwarnings = getInput( argv, argc, &block );
	
	srand(block.random_seed);

	/* copy parameters */
	
	int doPlanarTopology = block.planar_topology; // true means the system is planar.	
	int doMonge = block.monge;
	kc = block.kc;
	kg = block.kg;
	KA = block.KA;
	mode_KA = block.mode_KA;
	int nsteps = block.nsteps;
	int nequil = block.nequil;
	int do_srd = block.do_srd;
	int debug = block.debug;
	double temperature = 5.92E-01 * (block.T/300);
	double eta_SI = 8.90*(1e-4); // Joule second per meter cubed)	
	double eta = eta_SI * (1e-30 /* meters per angstrom*/ ) * (1./1000.) * (1/4.184) * (6.022e23);

	printf("eta_SI: %le eta: %le (kcal/angstrom/seconds)\n", eta_SI, eta );
//1.3e-13;
	double dt = block.time_step;

	srd_integrator *srd_i = NULL;
	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( block.meshName , 0. );
	surface *prev_surface = theSurface;
	surface *next_surface;

	int ndiv = 0;
	for( int d = 0; d < ndiv; d++ )
	{
		next_surface = (surface *)malloc( sizeof(surface) );
		int ssize = sizeof(surface);
		char *init = (char *)next_surface;

		next_surface->subdivideSurface( prev_surface );
		prev_surface = next_surface;
	}
	
	sub_surface = prev_surface;
	double Lx = sub_surface->PBC_vec[0][0];
	double Ly = sub_surface->PBC_vec[1][1];
	double Lz = sub_surface->PBC_vec[2][2];

	sub_surface->generatePlan();
	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 );
	double *M[3] = { M5, M6, M7 };
	int mlow = 5;
	int mhigh = 7;
	sub_surface->generateSubdivisionMatrices( M, mlow, mhigh );
	sub_surface->box_system();
	
	int nc = 3 * sub_surface->nv+3;
	int nv = sub_surface->nv;
	double *r = (double *)malloc( sizeof(double) * nc );
	sub_surface->get(r);
	double *g = (double *)malloc( sizeof(double) * nc );
	r[3*nv+0] = 1.0;
	r[3*nv+1] = 1.0;
	r[3*nv+2] = 1.0;
	sub_surface->setg0(r);
	
	double area0;
	double cur_area;
	sub_surface->area(r, -1, &cur_area, &area0 );
	printf("area: %le area0: %le\n", cur_area, area0 );




	double *ro = (double *)malloc( sizeof(double) * nc );
	memcpy( ro, r, sizeof(double) * nc );	

	double T = 300;

	int n = 16;
	
	double V0 = 0;
	double Vtot = 0;
	double Vtot2 = 0;
	double NV = 0;
	double sphere_rad = 1e10;

	if( block.sphere )
	{
		V0 =fabs(sub_surface->volume( r));
		printf("Area: %lf (R: %lf)\n", area0, sqrt(area0/(4*M_PI))  );
		printf("Volume: %lf (R: %lf)\n", V0, pow( V0/(4*M_PI/3.0), 1.0/3.0) );
		sphere_rad = sqrt(area0/(4*M_PI));
	}


	static gsl_rng * rng_x = NULL;
	static const gsl_rng_type *rng_T = gsl_rng_default;
        rng_x = gsl_rng_alloc(rng_T);
        gsl_rng_env_setup();	
	gsl_rng_set( rng_x, block.random_seed );

	double time_step = block.time_step; // one nanosecond.
	double time_step_collision = block.time_step_collision; // one nanosecond.

	// n / cubic angstrom
	double fix_alpha = 0.0;

	if( block.fix_alpha )
		fix_alpha = 1.0;
	
	FILE *tFile = NULL;
	if( block.movie )
	{
		FILE *tpsf = NULL;
		tpsf = fopen("traj.psf","w");
        	sub_surface->writeLimitingSurfacePSF(tpsf);
		fclose(tpsf);
	
		tFile = fopen("traj.xyz","w");
	}
	
	int nfaces = sub_surface->nf_faces + sub_surface->nf_irr_faces;

	// allocate a large array to potentially hold every face that may need to be recomputed if we break up or form a complex.
	int *modified_face_list = (int *)malloc( sizeof(int) * 3 * nfaces );
	int nmod = 0;


	int o_lim = nsteps;
	int o_start = nequil;

	
//#define FDIFF_GRAD

	if( block.loadName )
	{
		FILE *loadFile = fopen(block.loadName, "r");
	
		for( int x = 0; x < sub_surface->nv; x++ )
		{
			getLine(loadFile, buffer );
			if( feof(loadFile) ) break;
			sscanf(buffer, "%lf %lf %lf\n",
				r+3*x+0, r+3*x+1, r+3*x+2 );
		}
		getLine(loadFile, buffer );
		sscanf(buffer, "%lf %lf %lf\n", r+3*nv+0, r+3*nv+1, r+3*nv+2);

		fclose(loadFile);		
	}
	
	double hours = block.hours;
	struct timeval tstart;

	gettimeofday( &tstart, NULL );	

	if( hours > 0 )
	{
		o_lim = 2e9;		
	}
	
	int n_test_pts = 1;
	double *rvals = NULL;

	// Langevin dynamics
	

	double cur_t = 0;

	int done = 0;

	
	force_set *theForceSet = (force_set *)malloc( sizeof(force_set) );

	int plim = 10;
	int ppf = 0;

	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++ )
	{
		double f1 = fi / (double)plim;
		double f2 = fj / (double)plim;

		ppf++;
	}
	
	double mass_per_lipid = 1000 * 760.09 / 6.022e23 / 1000; // POPC kg, wikipedia
	double area_per_lipid = 65.35; // A^2, POPC, interpolated from Kucerka 2011 BBA 2761
	// factor of two for leaflets.
	double mass_per_pt = AMU_PER_KG * 2 * area0 / (sub_surface->nf_faces + sub_surface->nf_irr_faces) / ppf * mass_per_lipid / area_per_lipid; 
	printf("mass per pt: %le\n", mass_per_pt );

	theForceSet->npts = ppf * sub_surface->nt;

	theForceSet->face_set = (int *)malloc( sizeof(int) * theForceSet->npts );
	theForceSet->uv_set = (double *)malloc( sizeof(double) * theForceSet->npts *2 );
	theForceSet->mass = (double *)malloc( sizeof(double) * theForceSet->npts );

	int totp = 0;

	for( int t = 0; t < sub_surface->nt; t++ )
	{
		int f = sub_surface->theTriangles[t].f;
		
		for( int fi = 0; fi <= plim; fi++ )
		for( int fj = 0; fj <= plim-fi; fj++, totp++)
		{
			double f1 = fi / (double)plim;
			double f2 = fj / (double)plim;

			theForceSet->face_set[totp] = f;

			theForceSet->uv_set[2*totp+0] = f1;
			theForceSet->uv_set[2*totp+1] = f2;
			theForceSet->mass[totp] = mass_per_pt;
		}
	}
	
	sub_surface->load_least_squares_fitting( theForceSet );
	
	/*
		Langevin/leapfrog dynamics:

		positions stored for N sites on the surface.
		velocities
		accelerations...

		all are passed through a ``filter'' of the control points:
		the next control point positions are determined by the best fit to 
		dx: = half-step velocities times dt.
	*/
	
	// dynamics in terms of control points as generalized coordinates.	
	// we need partial r partial q
	// this is stored in the force set.
 
	double *xp = (double *)malloc( sizeof(double) * 3 * nv );
	double *vp = (double *)malloc( sizeof(double) * 3 * nv );
	double *ap = (double *)malloc( sizeof(double) * 3 * nv );	

	memset( vp, 0, sizeof(double) * 3 * nv );
	memset( ap, 0, sizeof(double) * 3 * nv );	

	double KE = 0;

	double *effective_mass = (double *)malloc( sizeof(double) * nv * nv );

	sub_surface->getEffectiveMass( theForceSet, effective_mass );
	
	double av_edge_length = 0;
	double n_edge_length = 0;
	for ( int x = 0; x < nv; x++ )
	{
		int val = sub_surface->theVertices[x].valence;

		for( int e = 0; e < val; e++ )
		{
			int j = sub_surface->theVertices[x].edges[e];

			double *epbc = sub_surface->theVertices[x].edge_PBC+3*e;

			double dr[3] = { 
				r[3*x+0] - r[3*j+0] - (epbc[0]*sub_surface->PBC_vec[0][0] + epbc[1] * sub_surface->PBC_vec[1][0] + epbc[2]*sub_surface->PBC_vec[2][0]), 
				r[3*x+1] - r[3*j+1] - (epbc[0]*sub_surface->PBC_vec[0][1] + epbc[1] * sub_surface->PBC_vec[1][1] + epbc[2]*sub_surface->PBC_vec[2][1]), 
				r[3*x+2] - r[3*j+2] - (epbc[0]*sub_surface->PBC_vec[0][2] + epbc[1] * sub_surface->PBC_vec[1][2] + epbc[2]*sub_surface->PBC_vec[2][2]) };
			double r = normalize(dr);

			av_edge_length += r;
			n_edge_length += 1;
		}
	}

	av_edge_length /= n_edge_length;

//#define PERTURB_START

#ifdef PERTURB_START
	if( doPlanarTopology )
	{
		for( int iq = 1; iq < 8; iq++ )
		{
			double q = 2 * M_PI * iq / sub_surface->PBC_vec[0][0];

			double mag = sqrt(temperature / (q*q*q*q*kc*area0));
			for( int x = 0; x < nv; x++ )
			{
				double xv = r[3*x+0];

				r[3*x+2] += mag * sin(xv*q)*2* (rand()/(double)RAND_MAX - 0.5);
			}
		}
	}
	else
	{

		for( int x = 0; x < nv; x++ )
		{
			double rand_mag = 5 * rand() / (double)RAND_MAX;
			r[3*x+0] += rand_mag * 2* (rand()/(double)RAND_MAX - 0.5);
			r[3*x+1] += rand_mag * 2* (rand()/(double)RAND_MAX - 0.5);
			r[3*x+2] += rand_mag * 2* (rand()/(double)RAND_MAX - 0.5);
		}
	}
#endif

	FILE *srdXYZFile = NULL;
	
	if( do_srd )
	{
		srd_i = (srd_integrator *)malloc( sizeof(srd_integrator) );
		
		srd_i->init( av_edge_length, sub_surface->PBC_vec, temperature, eta_SI, time_step_collision, time_step_collision * AKMA_TIME, doPlanarTopology, mass_per_pt, block.srd_M, block.hard_z_boundary ); 
		srd_i->initializeDistances( r, sub_surface, M, mlow, mhigh );
		if( debug )
			srd_i->activateDebugMode();
		srdXYZFile = fopen("srd.xyz", "w" );
	
	}

			
	double V = sub_surface->energy(r,NULL);
	printf("Initial energy: %lf\n", V );
	memset( g, 0, sizeof(double) * nc );
	sub_surface->grad( r, g );

	double *saved_ref_point = (double *)malloc( sizeof(double) * 3 * (nv+1) );	
	memcpy( saved_ref_point, r, sizeof(double) * 3 * (nv+1) );

	double *delta_hull = (double *)malloc( sizeof(double) * block.o_lim );

	double T_save = sub_surface->evaluate_T( vp ) + (do_srd ? srd_i->KE() : 0);
	double *drdt = (double*)malloc( sizeof(double) * 3 * (nv+1) );	

	int o = 0;
	double running_time = 0;
	while( !done )
	{
		memcpy( saved_ref_point, r, sizeof(double) * 3 * (nv+1) );
		if( do_srd )
			srd_i->initializeDistances( r, sub_surface, M, mlow, mhigh );

		double max_tpoint = 0;
		double delta_hull = 0;

		for( int t = 0; t < block.o_lim; t++, cur_t += time_step )
		{	
			double delta_hull = 0;

			for( int v = 0; v < nv; v++ )
			{
				double dr[3] = { r[3*v+0] - saved_ref_point[3*v+0],
						 r[3*v+1] - saved_ref_point[3*v+1],
						 r[3*v+2] - saved_ref_point[3*v+2] };
				double lr = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
				if( lr > delta_hull )
					delta_hull = lr;	
			}

			double use_dhull = max_tpoint + delta_hull;

//			printf("t: %d use_dhull: %le\n", t, use_dhull );

			if( delta_hull > max_tpoint )
				max_tpoint = delta_hull;

			// leapfrog
			
		
			for( int v1 = 0; v1 < nv; v1++ )
			for( int v2 = 0; v2 < nv; v2++ )
			{
				vp[3*v1+0] += (-g[3*v2+0]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step/2;
				vp[3*v1+1] += (-g[3*v2+1]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step/2;
				vp[3*v1+2] += (-g[3*v2+2]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step/2;
			}			
			
			memset( drdt, 0, sizeof(double) * 3 * (nv+1) );
			
			// position/velocity synchronized here.
	
			double T = sub_surface->evaluate_T(  vp );

			if( do_srd && debug )
			{
				double srd_T = srd_i->KE();
				
				printf("T after: %.14le + %.14le = %.14le dKE: %le \n", T, srd_T, T+srd_T, T+srd_T-T_save  );
				T_save = T+srd_T;
			}
		
			for( int v = 0; v < nv; v++ )
			{
				if( do_srd )
				{
					drdt[3*v+0] = vp[3*v+0];
					drdt[3*v+1] = vp[3*v+1];
					drdt[3*v+2] = vp[3*v+2];
				}
				else
				{
					r[3*v+0] += vp[3*v+0] * AKMA_TIME * time_step;
					r[3*v+1] += vp[3*v+1] * AKMA_TIME * time_step;
					r[3*v+2] += vp[3*v+2] * AKMA_TIME * time_step;
				}
			}
			sub_surface->put(r);			
	
			V = sub_surface->energy(r,NULL);
			memset( g, 0, sizeof(double) * nc );

			// dV/dq
			sub_surface->grad( r, g );
			
			int ncol = 0;
	
			sub_surface->rebox_system();

			if( do_srd )
			{
				// reinitialize everything if we are going to collide over and over.
				srd_i->initializeDistances( r, sub_surface, M, mlow, mhigh );
				use_dhull = 0;


				ncol += srd_i->stream_and_collide( r, g, drdt, vp, effective_mass, sub_surface, M, mlow, mhigh, use_dhull, theForceSet, cur_t*AKMA_TIME, time_step * AKMA_TIME, time_step_collision * AKMA_TIME  );
			}
			
			// get particle momentum change.

			double *dv;

			if( debug )
			{
				dv  = (double *)malloc( sizeof(double) * 3 * nv );
				memset( dv, 0, sizeof(double) * 3 * nv );
			}

			for( int v1 = 0; v1 < nv; v1++ )
			for( int v2 = 0; v2 < nv; v2++ )
			{
				vp[3*v1+0] += (-g[3*v2+0]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step/2;
				vp[3*v1+1] += (-g[3*v2+1]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step/2;
				vp[3*v1+2] += (-g[3*v2+2]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step/2;
		
				if( debug )
				{	
					// for debugging..	
					dv[3*v1+0] += (-g[3*v2+0]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step;
					dv[3*v1+1] += (-g[3*v2+1]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step;
					dv[3*v1+2] += (-g[3*v2+2]*effective_mass[v2*nv+v1]) * AKMA_TIME * time_step;
				}
			}		// comes in with particle mass, (momentum) divided by mass (velocity) multiplied by mass ( back to momentum )	

			if( debug )
			{
				// find change in momentum of underlying membrane
		
				double dp_tot[3] = {0,0,0};
				double p_tot[3] = {0,0,0};
				for( int p = 0; p < theForceSet->npts; p++ )
				{
					int nc = theForceSet->frc_ncoef[p];
					for( int c = 0; c < nc; c++ )
					{
						// 13 == maxv, just debugging here.
						int v = theForceSet->frc_coef_list[p*13+c];
		
						// 	
						dp_tot[0] += theForceSet->frc_coef[p*13+c] * dv[3*v+0] * theForceSet->mass[p]; 
						dp_tot[1] += theForceSet->frc_coef[p*13+c] * dv[3*v+1] * theForceSet->mass[p]; 
						dp_tot[2] += theForceSet->frc_coef[p*13+c] * dv[3*v+2] * theForceSet->mass[p]; 
	
						p_tot[0] += theForceSet->frc_coef[p*13+c] * vp[3*v+0] * theForceSet->mass[p]; 
						p_tot[1] += theForceSet->frc_coef[p*13+c] * vp[3*v+1] * theForceSet->mass[p]; 
						p_tot[2] += theForceSet->frc_coef[p*13+c] * vp[3*v+2] * theForceSet->mass[p]; 
					}
				} 
			
				printf("Total momentum: %le %le %le\n", p_tot[0], p_tot[1], p_tot[2] );		
				if( ncol > 0 )
					printf("Timestep membrane momentum change: %le %le %le\n", dp_tot[0], dp_tot[1], dp_tot[2] );
		
				free(dv);
			}	

//			if( do_srd && t % block.srd_collision_freq == 0 )
//				srd_i->collide();
			double dof = 3 * nv;

			double TEMP = 2 * T / dof;

			// in kcal/mol

			TEMP *= (300 / 0.592);

			double srd_ke =0;
			if( do_srd)
				srd_ke= srd_i->KE();

			if( do_srd )
				printf("t: %le ns T: %.8le V: %.8le T+V: %.14le TEMP: %le\n", (cur_t * 1e9), T+srd_ke,V,T+srd_ke+V, TEMP );
			else
				printf("t: %le ns T: %.8le V: %.8le T+V: %.14le TEMP: %le\n", (cur_t * 1e9), T,V,T+V, TEMP );
			fflush(stdout);
		}
		
		if( tFile )
		{
			sub_surface->put(r);
			//sub_surface->writeLimitStructure(tFile);
        		sub_surface->writeLimitingSurface(tFile);
			fflush(tFile);

			if( do_srd )
				srd_i->writeXYZ(srdXYZFile);
		}	
			

		
		if( block.hours > 0 )
		{
			struct timeval tnow;
			gettimeofday( &tnow, NULL );	
				
			double dt = ( (tnow.tv_sec+ tnow.tv_usec*1e-6) / 3600.0 - (tstart.tv_sec + tstart.tv_usec*1e-6) / 3600.0 );
	
			if( dt > hours )
			{
				printf("Currently we have run for %lf hours. Stopping now.\n", dt );
				break;
			}
			if( o % 50 == 0 && hours > 0 )
				printf("Currently we have run for %lf hours, will stop after %lf.\n", dt, hours );
			
		}
		else
		{
			if( o >= nsteps )
				done = 1;
		}

		o++;
	}
	
	if( tFile ) 
		fclose(tFile);
	FILE *saveFile = fopen("file.save", "w");

	for( int x = 0; x < sub_surface->nv; x++ )
	{
		fprintf(saveFile, "%.14le %.14le %.14le\n",
			r[3*x+0], r[3*x+1], r[3*x+2] );
	}
	fprintf(saveFile, "%.14le %.14le %.14le\n", r[3*nv+0], r[3*nv+1], r[3*nv+2] );

	fclose(saveFile);

	clearForceSet(theForceSet);
	free(effective_mass);	
	srd_i->clear();
}



void updateParticleR( int p, int *pfaces, double *puv, double *p_r_m, double *r, surface *sub_surface, int np )
{
	int nv = sub_surface->nv;

	double rp[3];
	double nrm[3];
	sub_surface->evaluateRNRM( pfaces[p], puv[2*p+0], puv[2*p+1], rp, nrm, r);

	p_r_m[3*p+0] = rp[0];			
	p_r_m[3*p+1] = rp[1];			
	p_r_m[3*p+2] = rp[2];		
	
	sub_surface->updateParticle( p_r_m+3*p, p, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
	
	if( do_2p )	
	{
		p_r_m[3*np+3*p+0] = rp[0] + dist_nrm * nrm[0];			
		p_r_m[3*np+3*p+1] = rp[1] + dist_nrm * nrm[1];			
		p_r_m[3*np+3*p+2] = rp[2] + dist_nrm * nrm[2];			
		sub_surface->updateParticle( p_r_m+3*np+3*p, np+p, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
	}
}




