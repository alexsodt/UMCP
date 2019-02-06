#include "local_config.h"
#include "parallel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "mutil.h"
#include "interp.h"
#include "l-bfgs.h"
#include <gsl/gsl_randist.h>
#include "m_triangles.h"
#include "input.h"
#include <assert.h>
#include "srd.h"
#include "units.h"
#include "pcomplex.h"
#include "p_p.h"
#include "npt.h"
//#define MM_METHOD_1
#define BOXED

//#define MONTE_CARLO_HACK

//#define SAVE_RESTARTS
#define NUM_SAVE_BUFFERS 20
//#define DEBUG_MOMENTUM_CHANGE
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
extern double VR;

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
int temp_main( int argc, char **argv );

int main( int argc, char **argv )
{
	return temp_main( argc, argv );
}

int temp_main( int argc, char **argv )
{
#ifdef PARALLEL
	int ierr;
	ierr = MPI_Init(&argc,&argv);
	if( ierr ) { exit(1); }
	quietParallel();
#endif
	printCompilationTime();

	int nprocs = 1;
	int taskid = 0;

	char buffer[4096];

	printf("No options are required:\n");
	printf("Syntax: ld [input file] [options ...]\n");
	printf("See input.C for options.\n");

	parameterBlock block;

	int nwarnings = getInput( (const char **)argv, argc, &block );
	
	srand(block.random_seed);

	/* copy parameters */

	int collect_hk = block.collect_hk;	
	int doPlanarTopology = block.planar_topology; // true means the system is planar.	
	int doMonge = block.monge;
	kc = block.kc;
	kg = block.kg;
	KA = block.KA;
	mode_KA = block.mode_KA;
	int nsteps = block.nsteps;
	int nequil = block.nequil;
	int do_srd = block.do_srd;
	int do_ld = block.do_ld;
	int debug = block.debug;
	double temperature = 5.92E-01 * (block.T/300);
	double eta_SI = block.eta; //8.90*(1e-4); // Joule second per meter cubed)	
	double eta = eta_SI * (1e-30 /* meters per angstrom*/ ) * (1./1000.) * (1/4.184) * (6.022e23);

	printf("eta_SI: %le eta: %le (kcal/angstrom/seconds)\n", eta_SI, eta );
//1.3e-13;
	double dt = block.time_step;

	srd_integrator *srd_i = NULL;
	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( block.meshName , 0. );
	surface *prev_surface = theSurface;
	surface *next_surface;
	

	pcomplex **allComplexes = NULL;

	int ncomplex = 0;


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
	
	int nc = 3 * sub_surface->nv+3;
	int nv = sub_surface->nv;
	double av_edge_length = 0;
	double n_edge_length = 0;
	double *r = (double *)malloc( sizeof(double) * nc );
	sub_surface->get(r);
	double *g = (double *)malloc( sizeof(double) * nc );
	r[3*nv+0] = 1.0;
	r[3*nv+1] = 1.0;
	r[3*nv+2] = 1.0;
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
	
	FILE *inputFile = fopen(argv[1],"r");
	sub_surface->readLipidComposition(inputFile);
	if( inputFile) fclose(inputFile);
	
	sub_surface->setg0(r);

	av_edge_length /= n_edge_length;
	sub_surface->box_system(av_edge_length);
	

	double area0;
	double cur_area;
	sub_surface->area(r, -1, &cur_area, &area0 );
	printf("area: %le area0: %le\n", cur_area, area0 );
	theSurface->loadComplexes( &allComplexes, &ncomplex, &block ); 

	int nmodes = 8;
	int nspacehk = 10;
	int nhk = 0;

	double *hk   = NULL; 
	double *lhq  = NULL; 
	double *lhq2 = NULL; 

	if( collect_hk )
	{
		hk = (double *)malloc( sizeof(double) * nspacehk * nmodes * nmodes );
		lhq  = (double *)malloc( sizeof(double) * nmodes * nmodes * 2 );
		lhq2 = (double *)malloc( sizeof(double) * nmodes * nmodes * 2 );
	}
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
	if( block.movie && par_info.my_id == BASE_TASK )
	{
		FILE *tpsf = NULL;
		char fname[256];
		sprintf(fname, "%s.psf", block.jobName );
		tpsf = fopen(fname,"w");
        	sub_surface->writeLimitingSurfacePSF(tpsf, allComplexes, ncomplex );
		fclose(tpsf);
	
		sprintf(fname, "%s.xyz", block.jobName );
		tFile = fopen(fname,"w");
	}
	
	int nfaces = sub_surface->nf_faces + sub_surface->nf_irr_faces;

	// allocate a large array to potentially hold every face that may need to be recomputed if we break up or form a complex.
	int *modified_face_list = (int *)malloc( sizeof(int) * 3 * nfaces );
	int nmod = 0;


	int o_lim = nsteps;

	
//#define FDIFF_GRAD
/*
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
*/	
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
	
	double MSCALE = 1;//500;
	double mass_per_lipid = MSCALE * 1000 * 760.09 / 6.022e23 / 1000; // POPC kg, wikipedia
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
	
	double *pp_m1   = (double *)malloc( sizeof(double) * 3 * nv );
 	double *qdot_m1 = (double *)malloc( sizeof(double) * 3 * nv );

	double *pp = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	double *next_pp = (double *)malloc( sizeof(double) * 3 * nv );
	double *qdot = (double *)malloc( sizeof(double) * 3 * nv );
	double *qdot0 = (double *)malloc( sizeof(double) * 3 * nv );
	double *qdot_temp = (double *)malloc( sizeof(double) * 3 * nv );
	double *ap = (double *)malloc( sizeof(double) * 3 * nv );	

	memset( pp, 0, sizeof(double) * 3 * nv );
	memset( ap, 0, sizeof(double) * 3 * nv );	
	memset( qdot, 0, sizeof(double) * 3 * nv );	
	memset( qdot0, 0, sizeof(double) * 3 * nv );	
	memset( qdot_temp, 0, sizeof(double) * 3 * nv );	

	double KE = 0;

	double *effective_mass = (double *)malloc( sizeof(double) * nv * nv );
	double *sparse_effective_mass = (double *)malloc( sizeof(double) * nv * nv );
	int *sparse_use = (int *)malloc( sizeof(int) * nv );
	int n_vuse = 0;
	sub_surface->getEffectiveMass( theForceSet, effective_mass );

//#define PERTURB_START

#ifdef PERTURB_START
	if( doPlanarTopology )
	{
		sub_surface->setup_mode_perturb( r, 1, 0, 16, 16, Lx, Ly );
		sub_surface->mode_perturb( r, 100.0, 0 );
		for( int iq = 1; iq <= 1; iq++ )
		{
			double q = 2 * M_PI * iq / sub_surface->PBC_vec[0][0];
	
			double kc_use = 14;
			double mag = 60; //sqrt(temperature / (q*q*q*q*kc_use*area0));

/*			for( int x = 0; x < nv; x++ )
			{
				double xv = r[3*x+0];

				r[3*x+2] += mag * sin(xv*q)*2* (rand()/(double)RAND_MAX - 0.5);
			}
*/
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

	double *saved_ref_point = (double *)malloc( sizeof(double) * 3 * (nv+1) );	
	memcpy( saved_ref_point, r, sizeof(double) * 3 * (nv+1) );

	double *delta_hull = (double *)malloc( sizeof(double) * block.o_lim );


	int o = 0;
	double running_time = 0;
	double gamma_langevin = block.gamma_langevin;
	// report average temperature
	double sum_average_temp = 0;
	double n_temp = 0;
 
	if( debug )
	{
		for( int x = 0; x < nv; x++ )
		{
			double rand_mag = 2500;
//			pp[3*x+0] += rand_mag * 2* (rand()/(double)RAND_MAX - 0.5);
//			pp[3*x+1] += rand_mag * 2* (rand()/(double)RAND_MAX - 0.5);
//			pp[3*x+2] += rand_mag * 2* (rand()/(double)RAND_MAX - 0.5);
		}
	}
	
	if( do_srd )		
		srd_i->clearDKE();

	if( block.loadName )
	{
		printf("Loading %s\n", block.loadName );
		FILE *xyzLoad = fopen(block.loadName,"r");

		if( !xyzLoad )
		{
			printf("Couldn't load save file '%s'.\n", block.loadName );
			exit(1);
		}
		for( int x = 0; x < nv+1; x++ )
		{
			getLine( xyzLoad, buffer );
			int nr = sscanf( buffer, "%lf %lf %lf %lf %lf %lf", r+3*x+0, r+3*x+1, r+3*x+2, pp+3*x+0, pp+3*x+1, pp+3*x+2 );

			if( nr != 3 && nr != 6 )
			{
				printf("LINE %d, %s\n", x+1, buffer );
				printf("ERROR reading load file '%s'.\n", block.loadName );
				exit(1);
			}	
		}

		for( int c = 0; c < ncomplex; c++ )
			allComplexes[c]->loadComplex(xyzLoad,sub_surface,r);
	}
	
	setupParallel( sub_surface, allComplexes, ncomplex );
	SparseMatrix *EFFM;
	sub_surface->getSparseEffectiveMass( theForceSet, sparse_effective_mass, sparse_use, &n_vuse, &EFFM);	
	setupSparseVertexPassing( EFFM, sub_surface->nv );


	if( block.nmin > 0 )
	{
		FILE *mpsf;
		FILE *minFile;

		if( par_info.my_id == BASE_TASK )
		{
			mpsf = fopen("minimize.psf","w");
        		sub_surface->writeLimitingSurfacePSF(mpsf, allComplexes, ncomplex);
			fclose(mpsf);
       
			minFile = fopen("minimize.xyz","w");
		}

		sub_surface->put(r);
		for( int c = 0; c < ncomplex; c++ ) allComplexes[c]->refresh(sub_surface, r );
		if( par_info.my_id == BASE_TASK )
		 	sub_surface->writeLimitingSurface(minFile, allComplexes, ncomplex );
		for( int m = 0; m < block.nmin; m++ )
		{
			sub_surface->minimize( r, allComplexes, ncomplex );
			sub_surface->put(r);
			for( int c = 0; c < ncomplex; c++ ) allComplexes[c]->refresh(sub_surface, r );
			if( par_info.my_id == BASE_TASK )
	 			sub_surface->writeLimitingSurface(minFile, allComplexes, ncomplex );
			if( par_info.my_id == BASE_TASK )
			{
				FILE *minSave = fopen("min.save","w");
			
				for( int x = 0; x < nv+1; x++ )
					fprintf( minSave, "%lf %lf %lf\n", r[3*x+0], r[3*x+1], r[3*x+2] );
	
				for( int c = 0; c < ncomplex; c++ )
					allComplexes[c]->saveComplex(minSave);
				fclose(minSave);
			}
		}
		if( par_info.my_id == BASE_TASK )
			fclose(minFile);
		
	}
		
//	sub_surface->debug_dynamics( r, theForceSet, effective_mass, allComplexes, ncomplex );
	
	if( block.timestep_analysis && taskid == BASE_TASK )
		sub_surface->timestep_analysis( r, theForceSet, effective_mass, allComplexes, ncomplex, dt );

	if( block.nsteps == 0 )
	{
		printf("Requested no dynamics, exitting.\n");
	}
	else
	{

	int switched = 0;

	double PV = 0;

	char *buffer_cycle[NUM_SAVE_BUFFERS];
	for( int x = 0; x < NUM_SAVE_BUFFERS; x++ )
		buffer_cycle[x] = NULL;
	int cur_save = 0;

	int needs_surface_collisions = do_srd;
		
	double step_rate = -1;

	int mom_init = 0;
	double prev_mom[3] = {0,0,0};
	double prev_TV = 0;

	memset( qdot0, 0, sizeof(double) * 3 * nv );
#ifdef MM_METHOD_1
	CartMatVecIncrScale( qdot0, pp, effective_mass, 1.0, nv, r+3*nv );
#elif defined(MM_METHOD_2)
	SparseCartMatVecIncrScale( qdot0, pp, sparse_effective_mass, 1.0, nv, sparse_use, n_vuse, r+3*nv );
#else
	AltSparseCartMatVecIncrScale( qdot0, pp, EFFM, 1.0, nv, sparse_use, n_vuse, r+3*nv );
#endif
	memcpy( qdot, qdot0, sizeof(double) * 3 * nv );
	
	memset( qdot_temp, 0, sizeof(double) * 3 * nv );	
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		allComplexes[c]->compute_qdot( sub_surface, r, effective_mass, qdot0, qdot_temp, pp );			
	}
#ifdef PARALLEL		
	ParallelSum( qdot_temp, 3*nv );
#endif
#ifdef MM_METHOD_1
	CartMatVecIncrScale( qdot, qdot_temp, effective_mass, 1.0, nv, r+3*nv );
#elif defined(MM_METHOD_2)
	SparseCartMatVecIncrScale( qdot, qdot_temp, sparse_effective_mass, 1.0, nv, sparse_use, n_vuse, r+3*nv );
#else
	AltSparseCartMatVecIncrScale( qdot, qdot_temp, EFFM, 1.0, nv, sparse_use, n_vuse, r+3*nv );
#endif
	sub_surface->grad( r, g );


	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		allComplexes[c]->update_dH_dq( sub_surface, r, effective_mass, g, pp, qdot, qdot0 );
	}
#ifdef PARALLEL
	// synchronizes particle positions at this point for computing particle-particle interactions. does not synchronize their gradient.
	if( ncomplex > 0 ) ParallelSyncComplexes( allComplexes, ncomplex );
#endif
			// particle-particle gradient, must be called in this order since save_grad is zero'd in update_dH_dq and added to here.
#ifndef BOXED
	PP_G( theSurface, r, allComplexes, ncomplex, g ); 
#else
	Boxed_PP_G( theSurface, r, allComplexes, ncomplex, g ); 
#endif

#ifdef PARALLEL
	ParallelSum(g,nc);
#endif
	
	memcpy( qdot_m1, qdot, sizeof(double) * 3 * nv );
	memcpy( pp_m1, pp, sizeof(double) * 3 * nv );

#ifdef PARALLEL
	double expec_time[par_info.nprocs];
	expec_time[par_info.my_id] = 0;
	for( int f = 0; f < par_info.nf; f++ )
	{
		if( par_info.faces[f] < sub_surface->nf_faces ) 
			expec_time[par_info.my_id] += sub_surface->nf_g_q_p;	
		else
			expec_time[par_info.my_id] += sub_surface->nf_irr_pts;	
	}
	ParallelGather(expec_time,1);
	double total_work = 0;
	double max_work = 0;
	for( int p =0; p < par_info.nprocs;p++ )
	{
		total_work += expec_time[p];
		if( expec_time[p] > max_work )
			max_work = expec_time[p];
	}
	printf("Max work is %lf times more than optimal.\n", par_info.nprocs * max_work / total_work );

#endif
	for( int f = 0; f < sub_surface->nt; f++ )
		sub_surface->set_g0_from_f(f);

//	sub_surface->debugDeformation( r );

	double p_min = -10000.0;
	double p_max = 10000.0;
	int nbins = 100;
	
	double *DIST0 = (double *)malloc( sizeof(double) * nbins );
	memset( DIST0, 0, sizeof(double) * nbins );


	int global_cntr = 0;
	struct timeval tnow;
	while( !done )
	{

		memcpy( saved_ref_point, r, sizeof(double) * 3 * (nv+1) );
//		if( do_srd )
//			srd_i->initializeDistances( r, sub_surface, M, mlow, mhigh );

		double max_tpoint = 0;
		double delta_hull = 0;


		gettimeofday( &tnow, NULL );	
	
		double time_0 = tnow.tv_sec + (1e-6) * tnow.tv_usec;

		double last_misc_time = 0;
		double last_mesh_time = 0;
		double last_wait_time = 0;
		double last_complex_time = 0;

		for( int t = 0; t < block.o_lim; t++, cur_t += time_step, global_cntr++ )
		{
			// leapfrog

			sub_surface->put(r);			


			/*********** COMPUTE ENERGY ************/	
	
			double T = sub_surface->evaluate_T( qdot, pp, NULL, NULL, r+3*nv ); //, qdot_m1, pp_m1 );
			VR=0;
			V = sub_surface->energy(r,NULL);
			memset( g, 0, sizeof(double) * nc );
		
			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];
				V += allComplexes[c]->V(sub_surface, r );	
				V += allComplexes[c]->AttachV(sub_surface, r );	
			}
			
                        /*********** END COMPUTE ENERGY ************/	

			/*********** SAVE RESTARTS FOR DEBUGGING ***/
#ifdef SAVE_RESTARTS
#ifdef PARALLEL
			if( ncomplex > 0 ) ParallelSyncComplexes( allComplexes, ncomplex );
#endif
			if( buffer_cycle[cur_save] )
				free(buffer_cycle[cur_save]);
			if( par_info.my_id == BASE_TASK )
				sub_surface->saveRestart( buffer_cycle+cur_save, r, pp, allComplexes, ncomplex );
			cur_save++;
			if( cur_save == NUM_SAVE_BUFFERS )
				cur_save = 0;
#endif
			/*********** END SAVE RESTARTS FOR DEBUGGING ***/

			/*********** SAVE PREVIOUS VELOCITIES *************/

			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];
		
				allComplexes[c]->cacheVelocities();
			}
	
			memcpy( qdot_m1, qdot, sizeof(double) * 3*nv);		
			memcpy( pp_m1,   pp,   sizeof(double) * 3*nv);		
	
			int b0 = nbins * (pp[0]- p_min)/(p_max-p_min);
			if( b0 < 0 ) b0 = 0;
			if( b0 >= nbins) b0 = nbins-1;
			DIST0[b0] += 1;

			//printf("b0: %d p: %le\n", b0, pp[0] );
			/********** END SAVE *******************************/


			if( block.lipid_mc_period > 0 )
			{
				PartialSyncVertices(r);
				if( global_cntr % block.lipid_mc_period == 0 )
				{
					sub_surface->local_lipidMCMove( r, allComplexes, ncomplex, time_step, 1.0 / temperature );
				}	 
			}
			
			if( block.npt_mc_period > 0 )
			{
				if( global_cntr % block.npt_mc_period == 0 )
					sub_surface->area_MC_move( r, allComplexes, ncomplex,  1.0 / temperature, qdot, pp, VOLUME_MOVE, &block );
			}
			if( block.cyl_tension_mc_period > 0 )
			{
				if( global_cntr % block.cyl_tension_mc_period == 0 )
					sub_surface->area_MC_move( r, allComplexes, ncomplex,  1.0 / temperature, qdot, pp, CYL_TENSION_MOVE, &block );
			}

			/*********** COMPUTE GRADIENT ******************/

			/* MESH gradient */

			// dV/dq

			gettimeofday(&tnow,NULL);
			double mesh_start = tnow.tv_sec +(1e-6)*tnow.tv_usec;

			sub_surface->grad( r, g );
			
			gettimeofday(&tnow,NULL);
			double mesh_stop = tnow.tv_sec +(1e-6)*tnow.tv_usec;
			
			last_mesh_time += mesh_stop-mesh_start;

			/* BEGIN complex */

			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];

				allComplexes[c]->prepareForGradient();
				allComplexes[c]->setrall(theSurface,r);
			}
			if( ncomplex > 0 ) 
			{
#ifdef PARALLEL
				// synchronizes particle positions at this point for computing particle-particle interactions. does not synchronize their gradient.
				ParallelSyncComplexes( allComplexes, ncomplex );
#endif
#ifndef BOXED
				V += PP_G( theSurface, r, allComplexes, ncomplex, g ); 
#else
				V += Boxed_PP_G( theSurface, r, allComplexes, ncomplex, g ); 
#endif			
				handleElasticCollisions( theSurface, allComplexes, ncomplex,  time_step*AKMA_TIME ); 
			}
			/*********** REPORT ENERGY, CHECK FOR FAILURE ****************/
#ifdef PARALLEL
			gettimeofday( &tnow, NULL );	
			double wait_time_start = tnow.tv_sec + (1e-6) * tnow.tv_usec;
			MPI_Barrier(MPI_COMM_WORLD);
			ParallelSum(&V,1);

			gettimeofday( &tnow, NULL );	
			double wait_time_stop = tnow.tv_sec + (1e-6) * tnow.tv_usec;

			last_wait_time += wait_time_stop - wait_time_start;
#endif

			if( ( V-PV > 1e5) && !(o==0 && t==0) )
			{
				for( int cx = 0; cx < NUM_SAVE_BUFFERS; cx++ )
				{
					if( buffer_cycle[cur_save] )
					{
						char fileName[256];
						sprintf(fileName, "%s_crash%d.save", block.jobName, cx );
						FILE *theFile = fopen(fileName,"w");
						fprintf(theFile, "%s", buffer_cycle[cur_save] );
						fclose(theFile);
					}
					cur_save ++;
		
					if( cur_save == NUM_SAVE_BUFFERS )
						cur_save = 0;
				}
				printf("System unstable at sub-iteration %d: V: %le Previous V: %le, bailing.\n", t, V, PV );
#ifdef PARALLEL 
     	 			 MPI_Barrier(MPI_COMM_WORLD);
      				  MPI_Finalize();
				exit(1);
#endif
			}

			PV = V;

			/*********** END REPORT, CHECK FOR FAILURE ****************/
			memset( qdot_temp, 0, sizeof(double) * 3 * nv );	
			

			/* BEGIN complex */
			gettimeofday(&tnow, NULL);
			double complex_time_start = tnow.tv_sec + (1e-6) * tnow.tv_usec;

			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];
				// currently has the PP gradient.

				int nsites = allComplexes[c]->nsites;
				double save_grad[3*nsites];
				memcpy( save_grad, allComplexes[c]->save_grad, sizeof(double)*3*nsites );		
		

				/*
					Propagation of the particle depends on the details of the surface metric.
					Near irregular vertices the move must be done in tiny pieces.

				*/

				double time_remaining = time_step;

				while( time_remaining > 0 )
				{
					memcpy( allComplexes[c]->save_grad, save_grad, sizeof(double)*3*nsites );

					double dt = allComplexes[c]->update_dH_dq( sub_surface, r, effective_mass, g, pp, qdot, qdot0, time_remaining, time_step );
			
					if( do_ld || o < nequil )
					{					
						allComplexes[c]->applyLangevinFriction( sub_surface, r, dt, gamma_langevin );
						allComplexes[c]->applyLangevinNoise( sub_surface, r, dt,  gamma_langevin, temperature );
					}
	
					allComplexes[c]->propagate_p( sub_surface, r, dt );
					allComplexes[c]->compute_qdot( sub_surface, r, effective_mass, qdot0, qdot_temp, pp, dt/time_step );			
					allComplexes[c]->propagate_q( sub_surface, r, dt );

					time_remaining -= dt;
				}
			}
			gettimeofday(&tnow, NULL);
			double complex_time_stop = tnow.tv_sec + (1e-6) * tnow.tv_usec;

			last_complex_time += complex_time_stop - complex_time_start;

			/* END complex */
			gettimeofday(&tnow,NULL);
			double misc_time_start = tnow.tv_sec + (1e-6)*tnow.tv_usec;
#ifdef PARALLEL
			if( ncomplex == 0 )
			{
				PartialSumVertices(g);
//				ParallelSum(g,nc);
				PartialSyncVertices(qdot_temp);
			}
			else
			{
				ParallelSum(g,nc);
				ParallelSum( qdot_temp, 3*nv );
			}
#endif

			memcpy( next_pp, pp, sizeof(double) * 3 * nv );

#ifdef MONTE_CARLO_HACK
			for( int v = 0; v < nv; v++ )
			{
				double dr[3] = { 
					2*(rand()-0.5)/RAND_MAX,
					2*(rand()-0.5)/RAND_MAX,
					2*(rand()-0.5)/RAND_MAX };
				dr[0] *= MOVE_MAG;
				dr[1] *= MOVE_MAG;
				dr[2] *= MOVE_MAG;
			
				r[3*v+0] += dr[0];		
				r[3*v+1] += dr[1];
				r[3*v+2] += dr[2];

				double E1 = energy( r, NULL );

				double pr = exp(-beta*(E1-E0) );
				double rn = rand()/(double)RAND_MAX;

				if( rn > pr )
				{
					r[3*v+0] -= dr[0];
					r[3*v+1] -= dr[1];
					r[3*v+2] -= dr[2];
				}	
			}
#else
			if( do_ld || o < nequil || (switched) )
			{
#ifdef MM_METHOD_1
				CartMatVecIncrScale( next_pp, pp, effective_mass, -gamma_langevin*AKMA_TIME * time_step, nv, r+3*nv );
#elif defined(MM_METHOD_2)
				SparseCartMatVecIncrScale( next_pp, pp, sparse_effective_mass, -gamma_langevin*AKMA_TIME * time_step, nv, sparse_use, n_vuse, r+3*nv);
#else
				AltSparseCartMatVecIncrScale( next_pp, pp, EFFM, -gamma_langevin*AKMA_TIME * time_step, nv, sparse_use, n_vuse, r+3*nv );
#endif
			}
			for( int v1 = 0; v1 < nv; v1++ )
			{
				next_pp[3*v1+0] += -g[3*v1+0] * AKMA_TIME * time_step;
				next_pp[3*v1+1] += -g[3*v1+1] * AKMA_TIME * time_step;
				next_pp[3*v1+2] += -g[3*v1+2] * AKMA_TIME * time_step;
			}	

			memcpy( pp, next_pp, sizeof(double) * 3 * nv );
			
			if( do_ld || o < nequil )
			{
				for( int v1 = 0; v1 < nv; v1++ )
				{
					double fx = gsl_ran_gaussian(rng_x, sqrt(2*gamma_langevin*temperature*AKMA_TIME*time_step) );
					double fy = gsl_ran_gaussian(rng_x, sqrt(2*gamma_langevin*temperature*AKMA_TIME*time_step) );
					double fz = gsl_ran_gaussian(rng_x, sqrt(2*gamma_langevin*temperature*AKMA_TIME*time_step) );
	
					pp[3*v1+0] += fx;				
					pp[3*v1+1] += fy;				
					pp[3*v1+2] += fz;				
				}
			}
			
			memset( qdot0, 0, sizeof(double) * 3 * nv );
#ifdef MM_METHOD_1
			CartMatVecIncrScale( qdot0, pp, effective_mass, 1.0, nv, r+3*nv );
#elif defined(MM_METHOD_2)
			SparseCartMatVecIncrScale( qdot0, pp, sparse_effective_mass, 1.0, nv, sparse_use, n_vuse, r+3*nv );
#else
			AltSparseCartMatVecIncrScale( qdot0, pp, EFFM, 1.0, nv, sparse_use, n_vuse, r+3*nv );
#endif
			memcpy( qdot, qdot0, sizeof(double) * 3 * nv );
#ifdef MM_METHOD_1
			CartMatVecIncrScale( qdot, qdot_temp, effective_mass, 1.0, nv, r+3*nv );
#elif defined(MM_METHOD_2)
			SparseCartMatVecIncrScale( qdot, qdot_temp, sparse_effective_mass, 1.0, nv, sparse_use, n_vuse, r+3*nv  );
#else
			AltSparseCartMatVecIncrScale( qdot, qdot_temp, EFFM, 1.0, nv, sparse_use, n_vuse, r+3*nv  );
#endif
			for( int v1 = 0; v1 < nv; v1++ )
			{
				r[3*v1+0] += qdot[3*v1+0] * AKMA_TIME * time_step;
				r[3*v1+1] += qdot[3*v1+1] * AKMA_TIME * time_step;
				r[3*v1+2] += qdot[3*v1+2] * AKMA_TIME * time_step;
			}

#ifdef PARALLEL
			PartialSyncVertices(pp);
//			ParallelBroadcast(pp,3*nv);	
#endif

#endif // MONTE_CARLO_HACK

			double PT = 0;
			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];
				PT += allComplexes[c]->T();
			}

#ifdef PARALLEL
			ParallelSum( &PT, 1 );
#endif
			T += PT;
			double dof = 3 * nv;

			for( int c = 0; c < ncomplex; c++ )
				dof += 3 * allComplexes[c]->nsites - allComplexes[c]->nattach;

			double TEMP = 2 * T / dof;

			// in kcal/mol

			TEMP *= (300 / 0.592);

			sum_average_temp += TEMP;
			n_temp += 1;

			fflush(stdout);
			
			if( t == 0 )
			{
				printf("t: %le ns T: %.8le V: %.8le T+V: %.14le TEMP: %le AV_TEMP %le VR: %.3le", (cur_t * 1e9), T,V,T+V, TEMP, sum_average_temp / n_temp,VR );
				if( step_rate > 0 )
					printf(" steps/s: %le", step_rate);
				printf("\n");
				fflush(stdout);
			}
			
			gettimeofday(&tnow,NULL);
			double misc_time_stop = tnow.tv_sec + (1e-6)*tnow.tv_usec;

			last_misc_time += misc_time_stop - misc_time_start;
			switched=0;
		}
		
#ifdef PARALLEL
		double times[par_info.nprocs*4];
		times[par_info.my_id*4+0] = last_mesh_time;
		times[par_info.my_id*4+1] = last_complex_time;
		times[par_info.my_id*4+2] = last_misc_time;
		times[par_info.my_id*4+3] = last_wait_time;

		ParallelGather(times,4);

		double total_mesh = 0, total_complex=0, total_misc=0, total_wait=0;
		
		double max_time = 0;
		double total_time = 0;
		for( int p = 0; p < par_info.nprocs; p++ )
		{	
			double local_time = times[p*4+0] + times[p*4+1] + times[p*4+2] + times[p*4+3];

			if( local_time > max_time ) max_time = local_time;
			total_mesh +=    times[p*4+0];
			total_complex += times[p*4+1];
			total_misc +=    times[p*4+2];
			total_wait +=    times[p*4+3];
		}

		double ineff = 100*max_time / (total_mesh+total_complex+total_misc);
		double factor = par_info.nprocs * (max_time)/(total_mesh+total_complex+total_misc);
		printf("Mesh: %4.2lf (s) Complex: %4.2lf (s) Wait: %.2lf (s) Misc: %4.12lf (s).\n",
			total_mesh, total_complex, total_wait, total_misc );
		printf("Most imbalanced process took %.2lf%% of the time (%lf times more than optimal).\n",
			ineff, factor );
#endif

		gettimeofday( &tnow, NULL );	
		double time_1 = tnow.tv_sec + (1e-6) * tnow.tv_usec;
		
		step_rate = block.o_lim / (time_1-time_0+1e-10);

		if( collect_hk )
		{
			sub_surface->directFT( lhq, lhq2, r, nmodes, nmodes );

			int *q_sorter = (int*)malloc( sizeof(int) * nmodes*nmodes);
			int ind = 0;
			for( int m1 = 0; m1 < nmodes; m1++ )
			for( int m2 = 0; m2 < nmodes; m2++,ind++ )
			{
				q_sorter[ind] = m1*nmodes+m2;
			}
			
			int done = 0;
			while( !done )
			{
				done = 1;
				for( int x = 0; x < nmodes*nmodes-1; x++ )
				{
					if( lhq2[2*q_sorter[x]+0] > lhq2[2*q_sorter[x+1]+0] )
					{
						done = 0;
						int t = q_sorter[x];
						q_sorter[x] = q_sorter[x+1];
						q_sorter[x+1] = t;
					}
				}
			}
			printf("MODES ");
			for( int q = 0; q < nmodes*nmodes; q++)
			{
				printf(" %le", lhq[2*q_sorter[q]+1]/area0);
			}
			printf("\n");
		
			double tsum = 0;	
			printf("MTEMPs ");
			for( int q = 0; q < nmodes*nmodes; q++)
			{
				// temperature of mode.
				//
				double qv = lhq2[2*q_sorter[q]+0];
				if( fabs(qv) < 1e-30 ) { printf(" 0.0 "); continue; }
		
				double expec_hk2 = temperature / ( kc * area0 * qv * qv * qv * qv );
				double val = lhq2[2*q_sorter[q]+1] / area0;
				double rat = val / expec_hk2;

				tsum += rat;

				printf(" %le", rat );
			}
			printf(" av: %lf\n", tsum / (nmodes*nmodes-1) );

			free(q_sorter);
		}
	
		if( ncomplex > 0 ) ParallelSyncComplexes( allComplexes, ncomplex );
		if( tFile && par_info.my_id == BASE_TASK )
		{
			sub_surface->put(r);
			//sub_surface->writeLimitStructure(tFile);
        		sub_surface->writeLimitingSurface(tFile, allComplexes, ncomplex,r+3*nv);
			fflush(tFile);
	
			if( do_srd )
				srd_i->writeXYZ(srdXYZFile);
		}	

		for( int b = 0; b < nbins; b++ )
			printf(" %lf", DIST0[b] );
		printf("\n");

	
		if( block.tachyon && par_info.my_id == BASE_TASK )
		{
			sub_surface->writeTachyon( block.jobName, block.tachyon_res, block.tachyon_interp, 
				r, allComplexes, ncomplex ); 
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

		printf("o: %d\n", o );
		if( o == nequil )
		{
			sum_average_temp = 0;
			n_temp = 0;
			switched = 1;
		}
		
		if( o == block.nve_switch )
		{
			do_ld = 0;
			do_srd = 0;
			sum_average_temp = 0;
			n_temp = 0;
			switched = 1;
		}
	}
	}
	
	if( tFile ) 
		fclose(tFile);

#ifdef PARALLEL
	if( ncomplex > 0 ) ParallelSyncComplexes( allComplexes, ncomplex );
#endif

	if( par_info.my_id == BASE_TASK )
	{
		char fname[256];
	
		sprintf(fname, "%s.save", block.jobName );

		FILE *saveFile = fopen( fname, "w");
	
		for( int x = 0; x < nv; x++ )
			fprintf( saveFile, "%lf %lf %lf %lf %lf %lf\n", r[3*x+0], r[3*x+1], r[3*x+2], pp[3*x+0], pp[3*x+1], pp[3*x+2] );
		
		for( int x = nv; x < nv+1; x++ )
			fprintf( saveFile, "%lf %lf %lf\n", r[3*x+0], r[3*x+1], r[3*x+2]  );
	
		for( int c = 0; c < ncomplex; c++ )
			allComplexes[c]->saveComplex(saveFile);
	
		fclose(saveFile);
	}

/*	FILE *saveFile = fopen("file.save", "w");

	for( int x = 0; x < sub_surface->nv; x++ )
	{
		fprintf(saveFile, "%.14le %.14le %.14le\n",
			r[3*x+0], r[3*x+1], r[3*x+2] );
	}
	fprintf(saveFile, "%.14le %.14le %.14le\n", r[3*nv+0], r[3*nv+1], r[3*nv+2] );

	fclose(saveFile);
*/
	clearForceSet(theForceSet);
	free(effective_mass);	
	if( do_srd )
		srd_i->clear();

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
#endif

	return 0;
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




