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
#define MM_METHOD_1
#define BOXED

//#define FIX_MEMBRANE
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
	double kcal_mol_K = (5.92186663194E-01/298);
	double temperature = kcal_mol_K * (block.T);
	printf("temperature: %le (kcal/mol) %le T\n", temperature, block.T );
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
	printf("mass per mesh vertex: %le\n", AMU_PER_KG * 2 * area0 * mass_per_lipid / area_per_lipid / (sub_surface->nv) );

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
 * for doing generalized coordinates (normal modes here).
 *
 * */
	
	double *scaling_factor = NULL;
	double *gen_transform = NULL;
	int NQ = 0;
	int do_gen_q = 0;
	int doing_spherical_harmonics = 0;
	int doing_planar_harmonics = 0;
	double *output_qvals=NULL;
	double *mass_scaling = NULL;

	if( block.mode_max >= 0 )
	{
		do_gen_q=1;
		if( block.sphere )
		{	
			doing_spherical_harmonics = 1;
			NQ = sub_surface->getSphericalHarmonicModes( r, block.mode_min, block.mode_max, &gen_transform, &output_qvals, &scaling_factor );

	

#if 0 
			doing_spherical_harmonics = 0;
			gen_transform = ( double *)malloc( sizeof(double)*3*nv*3*nv);
			memset( gen_transform, 0, sizeof(double) * 3 * nv * 3 * nv );
	
			for( int v = 0; v < nv*3; v++ )
				gen_transform[v*(3*nv)+v] = 1;
			NQ=3*nv;
#endif
		}
		else
		{
			doing_planar_harmonics = 1;
			NQ = sub_surface->getPlanarHarmonicModes( r, -1, -1, block.mode_min, block.mode_max, &gen_transform, &output_qvals, &scaling_factor );
		}
	}
	else if( block.mode_x >= 0 || block.mode_y >= 0 )
	{
		do_gen_q=1;
		if( block.sphere )
		{
			doing_spherical_harmonics = 1;
			printf("Single spherical harmonic NYI.\n");
			exit(1);
		}
		else
		{
			doing_planar_harmonics = 1;
			int max_l = block.mode_x;
			if( block.mode_y > max_l ) max_l = block.mode_y;

			NQ = sub_surface->getPlanarHarmonicModes( r, block.mode_x, block.mode_y, 0, max_l, &gen_transform, &output_qvals, &scaling_factor );
		}
	}
	
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
	
	double *pp=NULL;
	double *next_pp=NULL;
	double *del_pp=NULL; // for gen_q
	int n_real_q = 3*nv;

	double *nav_Q = NULL;
	double *av_Q = NULL;
	double *av_Q2 = NULL;
	double *av_Q_T = NULL;
	double *nav_Q_T = NULL;

	double *QV = NULL;
	double *Qdot = NULL;
	double *Qdot0 = NULL;
	double *Qdot0_trial = NULL;

	if( do_gen_q )
	{
		n_real_q = NQ;
		pp = (double *)malloc( sizeof(double) * NQ ); 
		memset( pp, 0, sizeof(double) * NQ );
		next_pp = (double *)malloc( sizeof(double) * NQ );
		memset( next_pp, 0, sizeof(double) * NQ );
		del_pp = (double *)malloc( sizeof(double) * NQ );
		QV = (double *)malloc( sizeof(double) * NQ );
		Qdot = (double *)malloc( sizeof(double) * NQ );
		Qdot0 = (double *)malloc( sizeof(double) * NQ );
		Qdot0_trial = (double *)malloc( sizeof(double) * NQ );
		nav_Q = (double *)malloc( sizeof(double) * NQ );
		av_Q = (double *)malloc( sizeof(double) * NQ );
		av_Q2 = (double *)malloc( sizeof(double) * NQ );
		av_Q_T = (double *)malloc( sizeof(double) * NQ );
		nav_Q_T = (double *)malloc( sizeof(double) * NQ );
		memset( nav_Q, 0, sizeof(double) * NQ );
		memset( av_Q, 0, sizeof(double) * NQ );
		memset( av_Q2, 0, sizeof(double) * NQ );
		memset( av_Q_T, 0, sizeof(double) * NQ );
		memset( nav_Q_T, 0, sizeof(double) * NQ );
	}
	else
	{	
		pp = (double *)malloc( sizeof(double) * (3*nv+3) );
		memset( pp, 0, sizeof(double) * (3*nv+3) );
		next_pp = (double *)malloc( sizeof(double) * (3*nv+3) );
		memset( next_pp, 0, sizeof(double) * (3*nv+3) );
		QV = (double *)malloc( sizeof(double) * (3*nv+3) );
		Qdot = (double *)malloc( sizeof(double) * (3*nv+3) );
		Qdot0 = (double *)malloc( sizeof(double) * (3*nv+3) );
		Qdot0_trial = (double *)malloc( sizeof(double) * (3*nv+3) );
	}
	


	double *qdot = (double *)malloc( sizeof(double) * 3 * nv );
	double *qdot0 = (double *)malloc( sizeof(double) * 3 * nv );
	
	
	double *qdot_temp = (double *)malloc( sizeof(double) * 3 * nv );

	double *ap = (double *)malloc( sizeof(double) * 3 * nv );	

	memset( pp, 0, sizeof(double) * NQ );
	memset( ap, 0, sizeof(double) * 3 * nv );	
	memset( qdot, 0, sizeof(double) * 3 * nv );	
	memset( qdot0, 0, sizeof(double) * 3 * nv );	
	memset( qdot_temp, 0, sizeof(double) * 3 * nv );	


	double KE = 0;

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

			
	memset( g, 0, sizeof(double) * nc );

	double *saved_ref_point = (double *)malloc( sizeof(double) * 3 * (nv+1) );	
	memcpy( saved_ref_point, r, sizeof(double) * 3 * (nv+1) );

	double *delta_hull = (double *)malloc( sizeof(double) * block.o_lim );


	int o = 0;
	double running_time = 0;
	double gamma_langevin = block.gamma_langevin;

	if( do_ld )
	{
		double prod =  gamma_langevin * time_step * AKMA_TIME;
	
		printf("gamma_langevin * time_step = %lf\n", gamma_langevin * time_step * AKMA_TIME );

		if( prod > 1 )	
			printf("WARNING: Langevin gamma * time_step greater than one: %lf\n", prod );
	}
	// report average temperature
	double sum_average_temp = 0;
	double n_temp = 0;
	
	if( do_srd )		
		srd_i->clearDKE();


	int debug_seed = -1;

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

			double mom[3];
			int nr = sscanf( buffer, "%lf %lf %lf %lf %lf %lf", r+3*x+0, r+3*x+1, r+3*x+2, mom, mom+1, mom+2 );

			if( nr != 3 && nr != 6 )
			{
				printf("LINE %d, %s\n", x+1, buffer );
				printf("ERROR reading load file '%s'.\n", block.loadName );
				exit(1);
			}	

			if( nr == 6 && ! do_gen_q )
			{
				pp[3*x+0] = mom[0];
				pp[3*x+1] = mom[1];
				pp[3*x+2] = mom[2];
			}
		}

		for( int c = 0; c < ncomplex; c++ )
			allComplexes[c]->loadComplex(xyzLoad,sub_surface,r);

		for( int Q = 0; Q < NQ; Q++ )
		{
			getLine(xyzLoad,buffer);
			if(feof(xyzLoad) ) break;
			double genp;
			int nr = sscanf(buffer, "GQ %lf\n", &genp );
			if( nr == 1 && do_gen_q )
				pp[Q] = genp;
		}

		int nr=	fscanf( xyzLoad, "seed %d", &debug_seed );
		if( nr > 0 )
			printf("Loaded seed %d\n", debug_seed );
	}
	

	setupParallel( sub_surface, allComplexes, ncomplex, ( do_gen_q ? NQ : 0) );
	SparseMatrix *EFFM;
	int max_mat = nv;
	if( NQ > nv )
		max_mat = NQ;
	int *sparse_use = (int *)malloc( sizeof(int) * (NQ > nv ? NQ : nv) );
	int n_vuse = 0;
	sub_surface->getSparseEffectiveMass( theForceSet, sparse_use, &n_vuse, &EFFM, gen_transform, NQ, mass_scaling );	
	setupSparseVertexPassing( EFFM, sub_surface->nv, do_gen_q );



	double V = sub_surface->energy(r,NULL);
#ifdef PARALLEL
	ParallelSum(&V,1);
#endif
	printf("Initial energy: %lf\n", V );
	
	if( !do_gen_q ) // generalized coordinates are simply the control point positions.
		NQ = 3 * nv;
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
			sub_surface->minimize( r, allComplexes, ncomplex, do_gen_q  );
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
		
	
//	if( block.timestep_analysis && taskid == BASE_TASK )
//		sub_surface->timestep_analysis( r, theForceSet, effective_mass, allComplexes, ncomplex, dt );

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
	if( do_gen_q )
	{
		memset( Qdot0, 0, sizeof(double) * NQ );
		GenQMatVecIncrScale( Qdot0, pp, EFFM, 1.0 );
		MatVec( gen_transform, Qdot0, qdot0, NQ, 3*nv ); 
	}
	else
		AltSparseCartMatVecIncrScale( qdot0, pp, EFFM, 1.0, r+3*nv );
	
	memcpy( qdot, qdot0, sizeof(double) * 3 * nv );
	
	memset( qdot_temp, 0, sizeof(double) * 3 * nv );	
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		allComplexes[c]->compute_qdot( sub_surface, r, qdot0, qdot_temp );			
	}
#ifdef PARALLEL		
	ParallelSum( qdot_temp, 3*nv );
#endif
	if( ! do_gen_q )
		AltSparseCartMatVecIncrScale( qdot, qdot_temp, EFFM, 1.0, r+3*nv );
	sub_surface->grad( r, g );


	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		allComplexes[c]->update_dH_dq( sub_surface, r, g, qdot, qdot0 );
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
	
	
	FILE *alphaFile = NULL;

	if( block.write_alpha_period > 0 )
	{
		char fileName[256];
		sprintf(fileName, "%s.alpha", block.jobName );
		
		alphaFile = fopen(fileName,"w");
	}



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
			
		if( debug && par_info.my_id == BASE_TASK )
		{
			char fileName[256];
			sprintf(fileName, "debug_%s.save", block.jobName );
			FILE *debugSave = fopen(fileName,"w");
		
			for( int x = 0; x < nv+1; x++ )
				fprintf( debugSave, "%lf %lf %lf\n", r[3*x+0], r[3*x+1], r[3*x+2] );
	
			for( int c = 0; c < ncomplex; c++ )
				allComplexes[c]->saveComplex(debugSave);
			int new_seed = rand();
			fprintf(debugSave, "seed %d\n", new_seed );
			srand(new_seed);
			fclose(debugSave);
		}
		else if( debug_seed >= 0 )
		{
			srand(debug_seed);
			debug_seed = -1;
		}

		for( int t = 0; t < block.o_lim; t++, cur_t += time_step, global_cntr++ )
		{
			// leapfrog

			sub_surface->put(r);			

			/*********** COMPUTE ENERGY ************/	
	
			VR=0;
			V = sub_surface->energy(r,NULL);
			double VMEM = V; 
			memset( g, 0, sizeof(double) * nc );
		
			double VP = 0;
			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];
				VP += allComplexes[c]->V(sub_surface, r );	
				VP += allComplexes[c]->AttachV(sub_surface, r );	
			}

			V += VP;
			
                        /*********** END COMPUTE ENERGY ************/	

			/*********** SAVE RESTARTS FOR DEBUGGING ***/
#ifdef SAVE_RESTARTS
#ifdef PARALLEL
			if( ncomplex > 0 ) ParallelSyncComplexes( allComplexes, ncomplex );
#endif
			if( buffer_cycle[cur_save] )
				free(buffer_cycle[cur_save]);
			if( par_info.my_id == BASE_TASK )
				sub_surface->saveRestart( buffer_cycle+cur_save, r, pp, allComplexes, ncomplex, NQ );
			cur_save++;
			if( cur_save == NUM_SAVE_BUFFERS )
				cur_save = 0;
#endif
			/*********** END SAVE RESTARTS FOR DEBUGGING ***/

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

			// LEAPFROG: gives p-dot at t, we have q(t), p(t-eps/2)
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
//				handleElasticCollisions( theSurface, r, allComplexes, ncomplex,  time_step*AKMA_TIME ); 
			}
			/*********** REPORT ENERGY, CHECK FOR FAILURE ****************/
#ifdef PARALLEL
			gettimeofday( &tnow, NULL );	
			double wait_time_start = tnow.tv_sec + (1e-6) * tnow.tv_usec;
			MPI_Barrier(MPI_COMM_WORLD);
			ParallelSum(&V,1);
			ParallelSum(&VMEM,1);
			ParallelSum(&VP,1);

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
			
			double PT = 0;

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

					double dt = allComplexes[c]->update_dH_dq( sub_surface, r, g, qdot, qdot0, time_remaining, time_step );
			
					if( do_ld || o < nequil )
					{					
						allComplexes[c]->applyLangevinFriction( sub_surface, r, dt, gamma_langevin );
						allComplexes[c]->applyLangevinNoise( sub_surface, r, dt,  gamma_langevin, temperature );
					}
	
					allComplexes[c]->propagate_p( sub_surface, r, dt/2 );
					allComplexes[c]->compute_qdot( sub_surface, r, qdot0, qdot_temp, dt/time_step );			
			
					// close enough.
					PT += allComplexes[c]->T(sub_surface,r)  * (dt/time_step);

					allComplexes[c]->propagate_p( sub_surface, r, dt/2 );
					allComplexes[c]->compute_qdot( sub_surface, r, qdot0, qdot_temp,  dt/time_step );			

					allComplexes[c]->propagate_surface_q( sub_surface, r, dt );

					time_remaining -= dt;
				}
			}
			
			// has special routines for handling elastic collisions.
			propagateSolutionParticles( sub_surface, r, allComplexes, ncomplex, time_step * AKMA_TIME );

			gettimeofday(&tnow, NULL);
			double complex_time_stop = tnow.tv_sec + (1e-6) * tnow.tv_usec;

			last_complex_time += complex_time_stop - complex_time_start;

			/* END complex */

			/* BEGIN SRD */
			
			if( do_srd )
			{
				sub_surface->rebox_system();
				double use_dhull = 0;

				double ncol = srd_i->stream_and_collide( r, g, qdot, EFFM, sub_surface, M, mlow, mhigh, use_dhull, theForceSet, cur_t*AKMA_TIME, time_step * AKMA_TIME, time_step_collision * AKMA_TIME  );
			}


			/* END SRD */

			gettimeofday(&tnow,NULL);
			double misc_time_start = tnow.tv_sec + (1e-6)*tnow.tv_usec;
#ifdef PARALLEL
			if( ncomplex == 0 )
			{
				if( do_gen_q )
				{

				}
				else
				{
					PartialSumVertices(g);
					PartialSyncVertices(qdot_temp);
				}
			}
			else
			{
				if( do_gen_q )
				{
				}
				else
				{
					ParallelSum(g,nc);
					ParallelSum( qdot_temp, 3*nv );
				}
			}
#endif

			memcpy( next_pp, pp, sizeof(double) * NQ );
			if( !block.disable_mesh )
			{
				if( do_ld || o < nequil || (switched) )
				{
					if( do_gen_q )
						GenQMatVecIncrScale( next_pp, pp, EFFM, -gamma_langevin*AKMA_TIME*time_step );
					else
						AltSparseCartMatVecIncrScale( next_pp, pp, EFFM, -gamma_langevin*AKMA_TIME * time_step, r+3*nv );
				}
			}

			// LEAPFROG: increment p by 1/2 eps, we have q(t), p(t), report properties for this state (perform Monte Carlo?)

			if( !block.disable_mesh )
			{

				if( do_gen_q )
				{
					memset(del_pp,0,sizeof(double)*NQ);
					
					for( int Q = 0; Q < NQ; Q++ )
					for( int v1 = 0; v1 < nv; v1++ )
					{
						del_pp[Q] += gen_transform[Q*3*nv+v1*3+0] * -g[3*v1+0] * AKMA_TIME * time_step/2;
						del_pp[Q] += gen_transform[Q*3*nv+v1*3+1] * -g[3*v1+1] * AKMA_TIME * time_step/2;
						del_pp[Q] += gen_transform[Q*3*nv+v1*3+2] * -g[3*v1+2] * AKMA_TIME * time_step/2;
					}

#ifdef PARALLEL
					ParallelSum( del_pp, NQ );
					ParallelBroadcast( del_pp, NQ );
#endif
					for( int Q = 0; Q < NQ; Q++ )
						next_pp[Q] += del_pp[Q];
				}
				else
				{
					for( int v1 = 0; v1 < nv; v1++ )
					{
						next_pp[3*v1+0] += -g[3*v1+0] * AKMA_TIME * time_step/2;
						next_pp[3*v1+1] += -g[3*v1+1] * AKMA_TIME * time_step/2;
						next_pp[3*v1+2] += -g[3*v1+2] * AKMA_TIME * time_step/2;
					}
				}	
			}
			memcpy( pp, next_pp, sizeof(double) * NQ );
			memset( Qdot0_trial, 0, sizeof(double) * NQ );

			if( do_gen_q )
				GenQMatVecIncrScale( Qdot0_trial, pp, EFFM, 1.0 );
			else
				AltSparseCartMatVecIncrScale( Qdot0_trial, pp, EFFM, 1.0, r+3*nv );

			double T = 0;				
			for( int Q = 0; Q < NQ; Q++ )
			{
				if( o >= nequil && do_gen_q ) { av_Q_T[Q] += pp[Q]*Qdot0_trial[Q] * 0.5; 
								nav_Q_T[Q] += 1;
//								printf("av_Q_T[%d] %le\n", Q, av_Q_T[Q]/nav_Q[Q]);
				}
				T += pp[Q] * Qdot0_trial[Q] * 0.5;
			}

			
			double srd_T = 0;
			double srd_dof = 0;

			if( do_srd )
			{
				double *vp = srd_i->vp;

				srd_dof = 3 * srd_i->np;
				for( int p = 0; p < srd_i->np; p++ )
					srd_T += 0.5 * srd_i->mass * ( vp[3*p+0]*vp[3*p+0]+vp[3*p+1]*vp[3*p+1]+vp[3*p+2]*vp[3*p+2] ); 
			}


#ifdef PARALLEL
//			if( ! do_gen_q )
//				ParallelSum(&T,1);
#endif
				
			/*****************************
 * 				Right now, we are on the *mesh* canonical ensemble with {q},{p}
 * 			 *****************************/


			if( !block.disable_mesh )
			{
				if( do_ld || o < nequil )
				{
					for( int Q = 0; Q < NQ; Q++ )
						next_pp[Q] += gsl_ran_gaussian(rng_x, sqrt(2*gamma_langevin*temperature*AKMA_TIME*time_step));

#ifdef PARALLEL
					ParallelBroadcast(next_pp,NQ);
#endif
				}

				// LEAPFROG: increment p by 1/2 eps, we have q(t), p(t+eps/2)
				if( do_gen_q )
				{
					for( int Q = 0; Q < NQ; Q++ )
						next_pp[Q] += del_pp[Q];
				}
				else
				{
					for( int v1 = 0; v1 < nv; v1++ )
					{
						next_pp[3*v1+0] += -g[3*v1+0] * AKMA_TIME * time_step/2;
						next_pp[3*v1+1] += -g[3*v1+1] * AKMA_TIME * time_step/2;
						next_pp[3*v1+2] += -g[3*v1+2] * AKMA_TIME * time_step/2;
					}	
				}
			}
			memcpy( pp, next_pp, sizeof(double) * NQ );
			memset( qdot0, 0, sizeof(double) * 3 * nv );

			if( do_gen_q )
			{	
				memset( Qdot0, 0, sizeof(double) * NQ );
				GenQMatVecIncrScale( Qdot0, pp, EFFM, 1.0 );
				MatVec( gen_transform, Qdot0, qdot0, NQ, 3*nv ); 
			}
			else
				AltSparseCartMatVecIncrScale( qdot0, pp, EFFM, 1.0, r+3*nv );

			memcpy( qdot, qdot0, sizeof(double) * 3 * nv );

			if( !do_gen_q )
				AltSparseCartMatVecIncrScale( qdot, qdot_temp, EFFM, 1.0, r+3*nv  );

			// LEAPFROG: increment q by eps, we have q(t+eps), p(t+eps/2)
			for( int v1 = 0; v1 < nv; v1++ )
			{
				r[3*v1+0] += qdot[3*v1+0] * AKMA_TIME * time_step;
				r[3*v1+1] += qdot[3*v1+1] * AKMA_TIME * time_step;
				r[3*v1+2] += qdot[3*v1+2] * AKMA_TIME * time_step;
			}
				

			if( do_gen_q )
			{

				for( int Q = 0; Q < NQ; Q++ )
				{
					QV[Q] += Qdot0[Q] * AKMA_TIME * time_step;
					if( o >= nequil )
					{
						av_Q[Q] += QV[Q];
						av_Q2[Q] += QV[Q]*QV[Q];
						nav_Q[Q] += 1;
					}
				}
			}
#ifdef PARALLEL
			if( do_gen_q )
			{

			}
			else
				PartialSyncVertices(pp);
#endif


#ifdef PARALLEL
			ParallelSum( &PT, 1 );
#endif
			double dof = NQ;
			if( !do_gen_q )
				dof = 3*nv-3;
			double mem_T = 2*T / dof;
			
			T += PT;
			T += srd_T;
			dof += srd_dof;
			for( int c = 0; c < ncomplex; c++ )
#ifdef DISABLE_ON_MEMBRANE_T
				dof += 3 * (allComplexes[c]->nsites - allComplexes[c]->nattach);
#else
				dof += 3 * allComplexes[c]->nsites - allComplexes[c]->nattach;
#endif
			double TEMP = 2 * T / dof;

			// in kcal/mol

			TEMP /= kcal_mol_K;
			mem_T /= kcal_mol_K;
			sum_average_temp += TEMP;
			n_temp += 1;

			fflush(stdout);
			
			if( t == 0 )
			{
				printf("t: %le ns o: %d T: %.8le V: %.8le T+V: %.14le TEMP: %le MEM_TEMP: %le AV_TEMP %le VR: %.3le VMEM: %le VP: %le", (cur_t * 1e9), o, T, V,T+V, TEMP, mem_T, sum_average_temp / n_temp,VR, VMEM, VP );
				if( step_rate > 0 )
					printf(" steps/s: %le", step_rate );
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
#ifdef PRINT_TIMING
		printf("Mesh: %4.2lf (s) Complex: %4.2lf (s) Wait: %.2lf (s) Misc: %4.12lf (s).\n",
			total_mesh, total_complex, total_wait, total_misc );
		printf("Most imbalanced process took %.2lf%% of the time (%lf times more than optimal).\n",
			ineff, factor );
#endif
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
	
#ifdef PARALLEL
		if( ncomplex > 0 ) ParallelSyncComplexes( allComplexes, ncomplex );
#endif
		if( tFile && par_info.my_id == BASE_TASK )
		{
			sub_surface->put(r);
			//sub_surface->writeLimitStructure(tFile);
        		sub_surface->writeLimitingSurface(tFile, allComplexes, ncomplex,r+3*nv);
			fflush(tFile);
	
			if( do_srd )
				srd_i->writeXYZ(srdXYZFile);
		}	

		if( block.write_alpha_period > 0 && alphaFile )
		{
			fprintf(alphaFile, "%le %le %le\n", r[3*nv+0], r[3*nv+1], r[3*nv+2] );
			fflush(alphaFile);
		}

	
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
		{
			if( do_gen_q )
				fprintf( saveFile, "%lf %lf %lf\n", r[3*x+0], r[3*x+1], r[3*x+2]  );
			else
				fprintf( saveFile, "%lf %lf %lf %lf %lf %lf\n", r[3*x+0], r[3*x+1], r[3*x+2], pp[3*x+0], pp[3*x+1], pp[3*x+2] );
		}
		for( int x = nv; x < nv+1; x++ )
			fprintf( saveFile, "%lf %lf %lf\n", r[3*x+0], r[3*x+1], r[3*x+2]  );
	
		for( int c = 0; c < ncomplex; c++ )
			allComplexes[c]->saveComplex(saveFile);
	
		for( int Q = 0; Q < NQ; Q++ )
			fprintf(saveFile, "GQ %.14le\n", pp[Q] );

		fclose(saveFile);
	}

	if( doing_spherical_harmonics || doing_planar_harmonics )
	{
		printf("------ Harmonic general variables ------\n");
		for( int Q = 0; Q < NQ; Q++ )
		{
			av_Q2[Q] *= scaling_factor[Q] * scaling_factor[Q];
			av_Q[Q] *= scaling_factor[Q];
			av_Q2[Q] /= nav_Q[Q];	
			av_Q[Q] /= nav_Q[Q];

			double kc_app = -1;
			double expected = 0;
			if( doing_spherical_harmonics )
			{	
				double l = output_qvals[Q];
				expected = temperature / ( kc * (l+2)*(l-1)*l*(l+1) ); 
			}
			else
			{	
				double q = output_qvals[Q];
				expected = temperature / ( kc * area0 *q*q*q*q); 
			}

			double var = av_Q2[Q] - av_Q[Q] * av_Q[Q];

			kc_app = kc * (expected / var);

			printf("Mode_index %d", Q );
			if( doing_planar_harmonics )
				printf(" q %le", output_qvals[Q] );
			else if( doing_spherical_harmonics )
				printf(" l %d", lround(output_qvals[Q]) ); 
			printf("<h> %le <h^2> %le k_c_apparent %le\n",
					av_Q[Q], av_Q2[Q], kc_app );
		}
		printf("------ done\n");
		
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




