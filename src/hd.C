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
#include "random_global.h"
#include "sans.h"
#include "globals.h"
#ifdef USE_CUDA
#include "local_cuda.h"
#endif
#include "simulation.h"
#include "meshCollisionLib.h"
#include "rd.h"

#define MM_METHOD_1
#define BOXED


int global_delete_this = 0;

#define OLD_LANGEVIN

//#define MONTE_CARLO_HACK

//#define SAVE_RESTARTS
#define NUM_SAVE_BUFFERS 50
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

double pp_grad( double *r, double *gc, int np, double *alpha, double *rads, int * pleaflet);
double particle_particle_energy( double *r, int np_1, int np_2, int *set, int np_set, double rad1, double rad2, double *alpha, double *rads, int *dimer_p, int *pleaflet );
void updateParticleR( int p, int *pfaces, double *puv, double *p_r_m, double *r, surface *theSurface, int np );

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

#ifdef USE_CUDA
	showCUDAStats();
#endif


	int nprocs = 1;
	int taskid = 0;

	char buffer[4096];

	printf("No options are required:\n");
	printf("Syntax: ld [input file] [options ...]\n");
	printf("See input.C for options.\n");

	parameterBlock block;

	int nwarnings = getInput( (const char **)argv, argc, &block );

	srand(block.random_seed);
	if( !rng_x ) init_random(block.random_seed);
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
	int on_surface = block.on_surface;
	int do_bd_membrane = block.do_bd_membrane;
	int do_bd_particles = block.do_bd_particles;
	int do_rd = block.do_rd;
	int debug = block.debug;
	global_debug_mode = block.debug;
	double kcal_mol_K = (5.92186663194E-01/298);
	kT = (5.92186663194E-01/298) * (block.T);
	double temperature = kcal_mol_K * (block.T);
	printf("temperature: %le (kcal/mol) %le T\n", temperature, block.T );
	double eta_SI = block.eta; //8.90*(1e-4); // Joule second per meter cubed)	
	double eta = eta_SI * (1e-30 /* meters per angstrom*/ ) * (1./1000.) * (1/4.184) * (6.022e23);

	printf("eta_SI: %le eta: %le (kcal/angstrom/seconds)\n", eta_SI, eta );
//1.3e-13;
	double dt = block.time_step;

	Simulation *theSimulation = (Simulation *)malloc( sizeof(Simulation) );
	theSimulation->visualization_cache = 0;

	srd_integrator *srd_i = NULL;
	theSimulation->allSurfaces = NULL;	

	surface *theSurface1 =(surface *)malloc( sizeof(surface) );
	double *r1 = NULL;
	theSurface1->loadLattice( block.meshName , 0. );

	for( int t = 0; t < theSurface1->nv; t++ )
	{
		theSurface1->theVertices[t].r[0] += block.shift[0];
		theSurface1->theVertices[t].r[1] += block.shift[1];
		theSurface1->theVertices[t].r[2] += block.shift[2];
	}



	double Lx = theSurface1->PBC_vec[0][0];
	double Ly = theSurface1->PBC_vec[1][1];
	double Lz = theSurface1->PBC_vec[2][2];

	theSurface1->generatePlan();

	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 );
	double *M[3] = { M5, M6, M7 };
	int mlow = 5;
	int mhigh = 7;
	theSurface1->generateSubdivisionMatrices( M, mlow, mhigh );

	// parameters from SOPC, table 1 from Rand/Parsegian BBA 988/1989/351-376, orig Ref Rand/Fuller/Parsegian/Rau 1988 v27 20 7711-7722, Table II
	double collision_alpha = 2.11*2;
	double collision_v0    = 2.95972e-7; // kcal/mol/A, per A^4
//	double collision_v0    = 1e-3;
	surface_record *newSurface = (struct surface_record*)malloc( sizeof(surface_record) );
	newSurface->theSurface = theSurface1;
	newSurface->id = 0;
	newSurface->next = theSimulation->allSurfaces;
	theSimulation->allSurfaces = newSurface;

	if( block.meshName2 )
	{
		// load and propagate a second mesh. just a hack for now until we simulate lots of objects.
	
		surface *addSurface = (surface *)malloc( sizeof(surface) );
		addSurface->loadLattice( block.meshName2, 0. );
	
		for( int t = 0; t < addSurface->nv; t++ )
		{
			addSurface->theVertices[t].r[0] += block.shift[0];
			addSurface->theVertices[t].r[1] += block.shift[1];
			addSurface->theVertices[t].r[2] += block.shift[2];
		}
	
		surface_record *newSurface = (struct surface_record*)malloc( sizeof(surface_record) );
		newSurface->theSurface = addSurface;
		newSurface->id = 1;
		newSurface->next = theSimulation->allSurfaces;
		theSimulation->allSurfaces = newSurface;
	
	}

	double av_edge_length = 0;
	double n_edge_length = 0;

	theSimulation->alpha[0] = 1.0;
	theSimulation->alpha[1] = 1.0;
	theSimulation->alpha[2] = 1.0;

	for( int x = 0; x < 3; x++ )
	for( int y = 0; y < 3; y++ )
		theSimulation->PBC_vec[x][y] = theSurface1->PBC_vec[x][y];


	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		surface *useSurface =sRec->theSurface;
		useSurface->surface_id = sRec->id;

		int nc = 3 * useSurface->nv+3;
		int nv = useSurface->nv;
		sRec->r = (double *)malloc( sizeof(double) * nc );
		sRec->nc = nc;
		useSurface->get(sRec->r);

		if( useSurface == theSurface1 )
			r1 = sRec->r;

		sRec->g = (double *)malloc( sizeof(double) * nc );
		memset( sRec->g, 0, sizeof(double) * nc );
		double *r = sRec->r;

		r[3*nv+0] = 1.0;
		r[3*nv+1] = 1.0;
		r[3*nv+2] = 1.0;

		for ( int x = 0; x < nv; x++ )
		{
			int val = useSurface->theVertices[x].valence;
	
			for( int e = 0; e < val; e++ )
			{
				int j = useSurface->theVertices[x].edges[e];
	
				double *epbc = useSurface->theVertices[x].edge_PBC+3*e;
	
				double dr[3] = { 
					r[3*x+0] - r[3*j+0] - (epbc[0]*useSurface->PBC_vec[0][0] + epbc[1] * useSurface->PBC_vec[1][0] + epbc[2]*useSurface->PBC_vec[2][0]), 
					r[3*x+1] - r[3*j+1] - (epbc[0]*useSurface->PBC_vec[0][1] + epbc[1] * useSurface->PBC_vec[1][1] + epbc[2]*useSurface->PBC_vec[2][1]), 
					r[3*x+2] - r[3*j+2] - (epbc[0]*useSurface->PBC_vec[0][2] + epbc[1] * useSurface->PBC_vec[1][2] + epbc[2]*useSurface->PBC_vec[2][2]) };
				double r = normalize(dr);
	
				av_edge_length += r;
				n_edge_length += 1;
			}
		}
	
		FILE *inputFile = fopen(argv[1],"r");
		useSurface->readLipidComposition(inputFile);
		if( inputFile) fclose(inputFile);
		useSurface->setg0(sRec->r);
	
		double area0;
		double cur_area;
		useSurface->area(r, -1, &cur_area, &area0 );
		printf("Surface %d area: %le area0: %le\n", sRec->id, cur_area, area0 );
	}
	
	av_edge_length /= n_edge_length;
	theSurface1->box_system(av_edge_length);
	
	// for now we are putting all complexes on the first surface.
	
	theSurface1->loadComplexes( &(theSimulation->allComplexes), &(theSimulation->ncomplex), &block ); 
	theSimulation->ncomplexSpace = theSimulation->ncomplex;

	// for now only collect kc info from first surface.
	// INITIALIZE REACTION DIFFUSION
	RD *rd = NULL;

	if(do_rd)
	{
		rd = (RD *)malloc( sizeof(RD) );
		rd->init(theSimulation, dt, &block );
		if(debug)
			printf("debug RD 1st ncomplexes: %d\n", theSimulation->ncomplex);

		int nsites_total = 0;
		for( int c = 0; c < theSimulation->ncomplex; c++ )
			nsites_total += theSimulation->allComplexes[c]->nsites;

		theSimulation->nsites_at_psfwrite = nsites_total;
		theSimulation->visualization_cache = 3 * nsites_total;
	}

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

	double T = 300;

	int n = 16;
	
	double Vtot = 0;
	double Vtot2 = 0;
	double NV = 0;
	double area0=0;

	if( block.sphere )
	{
		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		{
			double V0 = 0;
			V0 =fabs(theSurface1->volume( sRec->r));
//			printf("Area: %lf (R: %lf)\n", area0, sqrt(area0/(4*M_PI))  );
			printf("Volume: %lf (R: %lf)\n", V0, pow( V0/(4*M_PI/3.0), 1.0/3.0) );
		}
	}

	double *qvals = NULL;

	initSANS( &block, &qvals, &(block.nq) );
	double A2dz2_sampled = 0;
	double *B_hist = NULL;
	double sans_max_r = 0;
	double sans_bin_width = 0.25;
	int nsans_bins = 0;

	if( block.s_q )
	{
		sans_max_r = theSurface1->PBC_vec[0][0];
		if( sans_max_r < theSurface1->PBC_vec[1][1] )
			sans_max_r = theSurface1->PBC_vec[1][1];
		if( sans_max_r < theSurface1->PBC_vec[2][2] )
			sans_max_r = theSurface1->PBC_vec[2][2];
		nsans_bins = sans_max_r / sans_bin_width;
		B_hist = (double *)malloc( sizeof(double) * nsans_bins );
	}


//	static gsl_rng * rng_x = NULL;
//	static const gsl_rng_type *rng_T = gsl_rng_default;
  //      rng_x = gsl_rng_alloc(rng_T);
 //       gsl_rng_env_setup();	
//	gsl_rng_set( rng_x, block.random_seed );

	double time_step = block.time_step; // one nanosecond.
	double time_step_collision = block.time_step_collision; // one nanosecond.
	int nsurfaces = 0;
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		nsurfaces++;

	FILE *tpsf = NULL;
	char fname[256];
	sprintf(fname, "%s.psf", block.jobName );
	tpsf = fopen(fname,"w");
	sprintf(fname, "%s.xyz", block.jobName );
	FILE *tFile = fopen(fname,"w");
	theSimulation->writeLimitingSurfacePSF(tpsf);
	fclose(tpsf);
	int o_lim = nsteps;

	double hours = block.hours;
	struct timeval tstart;

	gettimeofday( &tstart, NULL );	

	if( hours > 0 )
		o_lim = 2e9;		
	
	// Langevin dynamics

	double cur_t = 0;

	int done = 0;

	// BLOCK FOR SETTING UP POINT EVALULATIONS ON SURFACE
	int plim = 10;
	int ppf = 0;

	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++ )
	{
		double f1 = fi / (double)plim;
		double f2 = fj / (double)plim;

		ppf++;
	}
	

	double min_mass_per_pt = -1;

	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		surface *useSurface = sRec->theSurface;

		int nc = 3 * useSurface->nv+3;
		int nv = useSurface->nv;

		sRec->theForceSet = (force_set *)malloc( sizeof(force_set) );
		force_set *theForceSet = sRec->theForceSet;

		theForceSet->npts = ppf * useSurface->nt;
		theForceSet->face_set = (int *)malloc( sizeof(int) * theForceSet->npts );
		theForceSet->uv_set = (double *)malloc( sizeof(double) * theForceSet->npts *2 );
		theForceSet->mass = (double *)malloc( sizeof(double) * theForceSet->npts );

		double cur_area,area0;
		useSurface->area(sRec->r, -1, &cur_area, &area0 );
	
		double MSCALE = 1;//500;
		double mass_per_lipid = MSCALE * 1000 * 760.09 / 6.022e23 / 1000; // POPC kg, wikipedia
		double area_per_lipid = 65.35; // A^2, POPC, interpolated from Kucerka 2011 BBA 2761
		// factor of two for leaflets.
		double mass_per_pt = AMU_PER_KG * 2 * area0 / (useSurface->nf_faces + useSurface->nf_irr_faces) / ppf * mass_per_lipid / area_per_lipid; 
		printf("mass per mesh vertex: %le\n", AMU_PER_KG * 2 * area0 * mass_per_lipid / area_per_lipid / (useSurface->nv) );

		if( mass_per_pt < min_mass_per_pt || min_mass_per_pt < 0 )
			min_mass_per_pt = mass_per_pt;

		int totp = 0;

		for( int t = 0; t < useSurface->nt; t++ )
		{
			int f = useSurface->theTriangles[t].f;
			
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
	
		useSurface->load_least_squares_fitting( theForceSet );
	
	// END BLOCK FOR SETTING UP POINT EVALULATIONS ON SURFACE

/*
 * for doing generalized coordinates (normal modes here).
 *
 * */
	
		sRec->scaling_factor = NULL;
		sRec->gen_transform = NULL;
		sRec->NQ = 0;
		sRec->do_gen_q = 0;
		sRec->doing_spherical_harmonics = 0;
		sRec->doing_planar_harmonics = 0;
		sRec->output_qvals = NULL;	

		if( block.mode_max >= 0 )
		{
			sRec->do_gen_q=1;
			if( block.sphere )
			{	
				sRec->doing_spherical_harmonics = 1;
				sRec->NQ = useSurface->getSphericalHarmonicModes( sRec->r, block.mode_min, block.mode_max, &sRec->gen_transform, &sRec->output_qvals, &sRec->scaling_factor );
			}
			else
			{
				sRec->doing_planar_harmonics = 1;
				sRec->NQ = useSurface->getPlanarHarmonicModes( sRec->r, -1, -1, block.mode_min, block.mode_max, &sRec->gen_transform, &sRec->output_qvals, &sRec->scaling_factor );
			}
		}
		else if( block.mode_x >= 0 || block.mode_y >= 0 )
		{
			sRec->do_gen_q=1;
			if( block.sphere )
			{
				sRec->doing_spherical_harmonics = 1;
				printf("Single spherical harmonic NYI.\n");
				exit(1);
			}
			else
			{
				sRec->doing_planar_harmonics = 1;
				int max_l = block.mode_x;
				if( block.mode_y > max_l ) max_l = block.mode_y;
	
				sRec->NQ = useSurface->getPlanarHarmonicModes( sRec->r, block.mode_x, block.mode_y, 0, max_l, &sRec->gen_transform, &sRec->output_qvals, &sRec->scaling_factor );
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

		sRec->pp = NULL;
		sRec->next_pp = NULL;
		sRec->del_pp = NULL;
		sRec->n_real_q = 3*nv;
		sRec->nav_Q = NULL;
		sRec->av_Q = NULL;
		sRec->av_Q2 = NULL;
		sRec->av_Q_T = NULL;
		sRec->nav_Q_T = NULL;
		sRec->QV = NULL;
		sRec->Qdot = NULL;
		sRec->Qdot0 = NULL;
		sRec->Qdot0_trial = NULL;		

		if( sRec->do_gen_q )
		{
			sRec->n_real_q = sRec->NQ;
			sRec->pp = (double *)malloc( sizeof(double) * sRec->NQ ); 
			memset( sRec->pp, 0, sizeof(double) * sRec->NQ );
			sRec->next_pp = (double *)malloc( sizeof(double) * sRec->NQ );
			memset( sRec->next_pp, 0, sizeof(double) * sRec->NQ );
			sRec->del_pp = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->QV = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->Qdot = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->Qdot0 = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->Qdot0_trial = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->nav_Q = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->av_Q = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->av_Q2 = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->av_Q_T = (double *)malloc( sizeof(double) * sRec->NQ );
			sRec->nav_Q_T = (double *)malloc( sizeof(double) * sRec->NQ );
			memset( sRec->nav_Q, 0, sizeof(double) * sRec->NQ );
			memset( sRec->av_Q, 0, sizeof(double) * sRec->NQ );
			memset( sRec->av_Q2, 0, sizeof(double) * sRec->NQ );
			memset( sRec->av_Q_T, 0, sizeof(double) * sRec->NQ );
			memset( sRec->nav_Q_T, 0, sizeof(double) * sRec->NQ );
		}
		else
		{
			sRec->pp = (double *)malloc( sizeof(double) * (3*nv+3) );
			memset( sRec->pp, 0, sizeof(double) * (3*nv+3) );
			sRec->next_pp = (double *)malloc( sizeof(double) * (3*nv+3) );
			memset( sRec->next_pp, 0, sizeof(double) * (3*nv+3) );
	
			sRec->QV = (double *)malloc( sizeof(double) * (3*nv+3) );
			sRec->Qdot = (double *)malloc( sizeof(double) * (3*nv+3) );
			sRec->Qdot0 = (double *)malloc( sizeof(double) * (3*nv+3) );
			sRec->Qdot0_trial = (double *)malloc( sizeof(double) * (3*nv+3) );
		}
		
		sRec->qdot = (double *)malloc( sizeof(double) * (3*nv+3) );
		sRec->qdot0 = (double *)malloc( sizeof(double) * (3*nv+3) );
		sRec->qdot_temp = (double *)malloc( sizeof(double) * (3*nv+3) );
	
		memset( sRec->pp, 0, sizeof(double) * sRec->NQ );
		memset( sRec->qdot, 0, sizeof(double) * 3 * nv );	
		memset( sRec->qdot0, 0, sizeof(double) * 3 * nv );	
		memset( sRec->qdot_temp, 0, sizeof(double) * 3 * nv );	
	}
	
	double KE = 0;

	FILE *srdXYZFile = NULL;
	
	if( do_srd )
	{
		srd_i = (srd_integrator *)malloc( sizeof(srd_integrator) );
		
		srd_i->init( av_edge_length, theSurface1->PBC_vec, temperature, eta_SI, time_step_collision, time_step_collision * AKMA_TIME, doPlanarTopology, min_mass_per_pt, block.srd_M, block.hard_z_boundary ); 
		srd_i->initializeDistances( theSimulation, M, mlow, mhigh );
		srdXYZFile = fopen("srd.xyz", "w" );
	
	}
			

	double *delta_hull = (double *)malloc( sizeof(double) * block.o_lim );

	int o = 0;
	double running_time = 0;
	double gamma_langevin = block.gamma_langevin;

	if( do_ld || do_bd_membrane || do_bd_particles )
	{
#ifndef OLD_LANGEVIN
		double prod =  gamma_langevin * time_step * AKMA_TIME;	
		printf("gamma_langevin * time_step = %lf\n", gamma_langevin * time_step * AKMA_TIME );
		if( prod > 1 )	
			printf("WARNING: Langevin gamma * time_step greater than one: %lf\n", prod );
#endif
	}
	// report average temperature
	double sum_average_temp = 0;
	double n_temp = 0;
	
	if( do_srd )		
		srd_i->clearDKE();


	int use_seed = block.random_seed;
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

		theSimulation->loadRestart(xyzLoad,&debug_seed);

		if( debug_seed > 0 )	
			use_seed = debug_seed;
#ifdef OLD_RESTART
		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		{
			double *r = sRec->r;
			int nv = sRec->theSurface->nv;
			double *pp = sRec->pp;

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
	
				if( nr == 6 && ! sRec->do_gen_q )
				{
					pp[3*x+0] = mom[0];
					pp[3*x+1] = mom[1];
					pp[3*x+2] = mom[2];
				}
			}
		}
		for( int c = 0; c < theSimulation->ncomplex; c++ )
			theSimulation->allComplexes[c]->loadComplex(xyzLoad,theSimulation,theSurface1->surface_id );
		
		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		{
			double *pp = sRec->pp;
			int NQ = sRec->NQ;

			for( int Q = 0; Q < NQ; Q++ )
			{
				getLine(xyzLoad,buffer);
				if(feof(xyzLoad) ) break;
				double genp;
				int nr = sscanf(buffer, "GQ %lf\n", &genp );
				if( nr == 1 && sRec->do_gen_q )
					pp[Q] = genp;
			}
		}
		int nr = fscanf( xyzLoad, "seed %d", &debug_seed );

		if( nr > 0 )
		{
			printf("Loaded seed %d\n", debug_seed );
			use_seed = debug_seed;
		}
#endif
	}

	
	setupParallel( theSimulation ); 

	if( nsteps > 0 )
	{	
		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		{
			surface *useSurface = sRec->theSurface;
			int nv = useSurface->nv;
			int NQ = sRec->NQ;
	
			sRec->EFFM = NULL;
			sRec->MMat = NULL;
	
			int max_mat = nv;
	
			if( NQ > nv )
				max_mat = NQ;
			int *sparse_use = (int *)malloc( sizeof(int) * (NQ > nv ? NQ : nv) );
			int n_vuse = 0;
			double *mass_scaling = NULL;
		
			useSurface->getSparseEffectiveMass( sRec->theForceSet, sparse_use, &n_vuse, &sRec->EFFM, sRec->gen_transform, NQ, mass_scaling );	
			if( do_bd_membrane )
				useSurface->getSparseRoot( sRec->theForceSet, &sRec->MMat ); 
	
			free(sparse_use);
		}
			
		setupSparseVertexPassing( theSimulation );
	}

	double V = 0;
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		V += sRec->theSurface->energy( sRec->r,NULL);
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		V += theSimulation->allComplexes[c]->V(theSimulation);	
		V += theSimulation->allComplexes[c]->AttachV(theSimulation);	
	}
	V += Boxed_PP_V( theSimulation ); 
#ifdef PARALLEL
	ParallelSum(&V,1);
#endif
	printf("Initial energy: %lf\n", V );

	int any_do_gen_q = 0;
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		if( !sRec->do_gen_q ) // generalized coordinates are simply the control point positions.
			sRec->NQ = 3 * sRec->theSurface->nv;
		else
			any_do_gen_q = 1;	
	}


	if( block.nmin > 0 )
	{
		FILE *mpsf;
		FILE *minFile;

		if( par_info.my_id == BASE_TASK )
		{
			mpsf = fopen("minimize.psf","w");
        		theSimulation->writeLimitingSurfacePSF(mpsf);
			fclose(mpsf);
       
			minFile = fopen("minimize.xyz","w");
		}

		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			sRec->theSurface->put(sRec->r);
		for( int c = 0; c < theSimulation->ncomplex; c++ ) theSimulation->allComplexes[c]->refresh(theSimulation );
		if( par_info.my_id == BASE_TASK )
		 	theSimulation->writeLimitingSurface(minFile );
		for( int m = 0; m < block.nmin; m++ )
		{
			theSimulation->minimize(any_do_gen_q);
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
				sRec->theSurface->put(sRec->r);
			for( int c = 0; c < theSimulation->ncomplex; c++ ) theSimulation->allComplexes[c]->refresh(theSimulation);
			if( par_info.my_id == BASE_TASK )
	 			theSimulation->writeLimitingSurface(minFile);
		}
		if( par_info.my_id == BASE_TASK )
		{
			FILE *minSave = fopen("min.save","w");
			theSimulation->saveRestart(minSave,-1);
			fclose(minFile);
		}
		
	}
	
//	if( block.timestep_analysis && taskid == BASE_TASK )
//		theSurface->timestep_analysis( r, theForceSet, effective_mass, allComplexes, theSimulation->ncomplex, dt );

	double	navc = 0;
	double *avc = NULL;
	if( block.record_curvature )
	{
		avc = (double *)malloc( sizeof(double) * theSimulation->ncomplex );
		memset( avc, 0, sizeof(double) * theSimulation->ncomplex );
	}

	if( nsteps == 0 )
	{

#ifdef MULTISURFACE_TACHYON
		if( block.tachyon && par_info.my_id == BASE_TASK )
		{
			printf("Writing tachyon object file.\n");
			theSurface->writeTachyon( block.jobName, block.tachyon_res, 1, 
				r, theSimulation->allComplexes, theSimulation->ncomplex, &block, srd_i, block.tachyon_tri_center ); 
		}
#endif
		printf("Requested no dynamics, exiting.\n");
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

	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		memset( sRec->qdot0, 0, sizeof(double) * 3 * sRec->theSurface->nv );
		if( sRec->do_gen_q )
		{
			memset( sRec->Qdot0, 0, sizeof(double) * sRec->NQ );
			GenQMatVecIncrScale( sRec->Qdot0, sRec->pp, sRec->EFFM, 1.0  );
			MatVec( sRec->gen_transform, sRec->Qdot0, sRec->qdot0, sRec->NQ, 3*sRec->theSurface->nv ); 
		}
		else
			AltSparseCartMatVecIncrScale( sRec->qdot0, sRec->pp, sRec->EFFM, 1.0, theSimulation->alpha, sRec->id );
		
		memcpy( sRec->qdot, sRec->qdot0, sizeof(double) * 3 * sRec->theSurface->nv );
		memset( sRec->qdot_temp, 0, sizeof(double) * 3 * sRec->theSurface->nv );
	}
	
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		theSimulation->allComplexes[c]->compute_qdot( theSimulation );			
	}
#ifdef PARALLEL
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		ParallelSum( sRec->qdot_temp, 3*sRec->theSurface->nv );
#endif
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		if( ! sRec->do_gen_q )
			AltSparseCartMatVecIncrScale( sRec->qdot, sRec->qdot_temp, sRec->EFFM, 1.0, theSimulation->alpha, sRec->id );
		sRec->theSurface->grad( sRec->r, sRec->g );
	}


	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		// NEEDS TO BE UPDATED FOR WHICH SURFACE THE PARTICLE IS ON.
		theSimulation->allComplexes[c]->update_dH_dq( theSimulation );
	}
#ifdef PARALLEL
	// synchronizes particle positions at this point for computing particle-particle interactions. does not synchronize their gradient.
	if( theSimulation->ncomplex > 0 ) ParallelSyncComplexes( theSimulation->allComplexes, theSimulation->ncomplex );
#endif
			// particle-particle gradient, must be called in this order since save_grad is zero'd in update_dH_dq and added to here.
#ifndef BOXED
	PP_G( theSimulation  ); 
#else
	Boxed_PP_G( theSimulation  ); 
#endif

#ifdef PARALLEL
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		ParallelSum(sRec->g,sRec->nc);
#endif
#ifdef PARALLEL
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		double expec_time[par_info.nprocs];
		expec_time[par_info.my_id] = 0;
		for( int f = 0; f < par_info.nf[sRec->id]; f++ )
		{
			if( par_info.faces[sRec->id][f] < sRec->theSurface->nf_faces ) 
				expec_time[par_info.my_id] += sRec->theSurface->nf_g_q_p;	
			else
				expec_time[par_info.my_id] += sRec->theSurface->nf_irr_pts;	
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
	}
#endif
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	for( int f = 0; f < sRec->theSurface->nt; f++ )
		sRec->theSurface->set_g0_from_f(f);

//	theSurface->debugDeformation( r );

	double p_min = -10000.0;
	double p_max = 10000.0;
	int nbins = 100;
	
	
	if( do_bd_particles )
	{
		for( int p = 0; p < theSimulation->ncomplex; p++ )
			theSimulation->allComplexes[p]->activateBrownianDynamics();
	}
	FILE *alphaFile = NULL;

	if( block.write_alpha_period > 0 )
	{
		char fileName[256];
		sprintf(fileName, "%s.alpha", block.jobName );
		
		alphaFile = fopen(fileName,"w");
	}

#ifndef OLD_LANGEVIN
	double *noise_vec = (double *)malloc( sizeof(double) * 3 * nv );
#endif

	srand(use_seed);
	my_gsl_reseed( use_seed );

	int global_cntr = 0;
	struct timeval tnow;


	while( !done )
	{

//		if( do_srd )
//			srd_i->initializeDistances( r, theSurface, M, mlow, mhigh );

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
#ifdef SAVE_RESTARTS
			int save_seed = rand();
			srand(save_seed);
			my_gsl_reseed(save_seed);
#endif
			if( debug_seed >= 0 )
			{
				srand(debug_seed);
				my_gsl_reseed(debug_seed);
				debug_seed = -1;
			}

			// leapfrog

			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
				sRec->theSurface->put(sRec->r);			

			/*********** COMPUTE ENERGY ************/	
			V=0;
			VR=0;
			double VMEM = 0;
			double VP = 0;
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			{
				V += sRec->theSurface->energy(sRec->r,NULL);

				VMEM += V; 
				memset( sRec->g, 0, sizeof(double) * sRec->nc );
			}

			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];
				VP += theSimulation->allComplexes[c]->V(theSimulation);	
				VP += theSimulation->allComplexes[c]->AttachV(theSimulation);	
			}


			V += VP;
			
                        /*********** END COMPUTE ENERGY ************/	
                                                
			/*********** SAVE RESTARTS FOR DEBUGGING ***/
#if 0 // PRESIMULATION, requires updating
#ifdef SAVE_RESTARTS
#ifdef PARALLEL
			if( theSimulation->ncomplex > 0 ) ParallelSyncComplexes( allComplexes, theSimulation->ncomplex );
#endif
			if( buffer_cycle[cur_save] )
				free(buffer_cycle[cur_save]);
			if( par_info.my_id == BASE_TASK )
				theSurface->saveRestart( buffer_cycle+cur_save, r, pp, allComplexes, theSimulation->ncomplex, (do_gen_q ? NQ : 0), save_seed );
			cur_save++;
			if( cur_save == NUM_SAVE_BUFFERS )
				cur_save = 0;
#endif
#endif
			/*********** END SAVE RESTARTS FOR DEBUGGING ***/

			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			{
				if( block.lipid_mc_period > 0 )
				{
					PartialSyncVertices(sRec->r,sRec->id);
					if( global_cntr % block.lipid_mc_period == 0 )
					{
						sRec->theSurface->local_lipidMCMove( sRec->r, theSimulation->allComplexes, theSimulation->ncomplex, time_step, 1.0 / temperature );
					}
				}	 
			
			}
			if( block.cyl_tension_mc_period > 0 )
			{
				if( global_cntr % block.cyl_tension_mc_period == 0 )
					theSimulation->area_MC_move( 1.0 / temperature, CYL_TENSION_MOVE, &block );
			}
				
			if( block.npt_mc_period > 0 )
			{
				if( global_cntr % block.npt_mc_period == 0 )
					theSimulation->area_MC_move( 1.0 / temperature, VOLUME_MOVE, &block );
			}

			/*********** COMPUTE GRADIENT ******************/

			/* MESH gradient */

			// dV/dq

			gettimeofday(&tnow,NULL);
			double mesh_start = tnow.tv_sec +(1e-6)*tnow.tv_usec;

			// LEAPFROG: gives p-dot at t, we have q(t), p(t-eps/2)
			
			
			// surface-surface collision

			for( surface_record *sRec1 = theSimulation->allSurfaces; sRec1; sRec1 = sRec1->next )
			for( surface_record *sRec2 = theSimulation->allSurfaces; sRec2; sRec2 = sRec2->next )
			{
				if( sRec2->id <= sRec1->id )
					continue;
				surfaceSurfaceCollisionForces( 
					sRec1->theSurface,
					sRec2->theSurface,
					sRec1->g,		
					sRec2->g,
					collision_alpha,
					collision_v0,
					M, mlow, mhigh );
			}
			
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
				sRec->theSurface->grad( sRec->r, sRec->g );

			gettimeofday(&tnow,NULL);
			double mesh_stop = tnow.tv_sec +(1e-6)*tnow.tv_usec;
			
			last_mesh_time += mesh_stop-mesh_start;

			/* BEGIN complex */

			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];

				theSimulation->allComplexes[c]->prepareForGradient();
				theSimulation->allComplexes[c]->setrall(theSimulation);
			}
			if( theSimulation->ncomplex > 0 ) 
			{
#ifdef PARALLEL
				// synchronizes particle positions at this point for computing particle-particle interactions. does not synchronize their gradient.
				ParallelSyncComplexes( theSimulation->allComplexes, theSimulation->ncomplex );
#endif
#ifndef BOXED
				V += PP_G( theSimulation ); 
#else
				V += Boxed_PP_G( theSimulation ); 
#endif			
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

			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
				memset( sRec->qdot_temp, 0, sizeof(double) * 3 * sRec->theSurface->nv );	

			/* BEGIN complex */
			gettimeofday(&tnow, NULL);
			double complex_time_start = tnow.tv_sec + (1e-6) * tnow.tv_usec;
			
			double PT = 0;

			for( int cx = 0; cx < par_info.nc; cx++ )
			{
				int c = par_info.complexes[cx];
				// currently has the PP gradient.

				int nsites = theSimulation->allComplexes[c]->nsites;
				double save_grad[3*nsites];
				memcpy( save_grad, theSimulation->allComplexes[c]->save_grad, sizeof(double)*3*nsites );		
		

				/*
					Propagation of the particle depends on the details of the surface metric.
					Near irregular vertices the move must be done in tiny pieces.

				*/

				double time_remaining = time_step;

				while( time_remaining > 0 )
				{
					memcpy( theSimulation->allComplexes[c]->save_grad, save_grad, sizeof(double)*3*nsites );

					double dt = theSimulation->allComplexes[c]->update_dH_dq( theSimulation, time_remaining, time_step );
			
					if(  do_ld || o < nequil )
					{					
						theSimulation->allComplexes[c]->applyLangevinFriction( theSimulation, dt, gamma_langevin );
						theSimulation->allComplexes[c]->applyLangevinNoise( theSimulation, dt,  gamma_langevin, temperature );
					}

					if( !theSimulation->allComplexes[c]->do_bd )
					{
						theSimulation->allComplexes[c]->propagate_p( theSimulation, dt/2 );
						theSimulation->allComplexes[c]->compute_qdot( theSimulation, dt/time_step );			
			
						// close enough.
						PT += theSimulation->allComplexes[c]->T(theSimulation)  * (dt/time_step);
	
						theSimulation->allComplexes[c]->propagate_p( theSimulation, dt/2 );
						theSimulation->allComplexes[c]->compute_qdot( theSimulation,  dt/time_step );			
					}
						
					theSimulation->allComplexes[c]->propagate_surface_q( theSimulation, dt );

					time_remaining -= dt;
				}
			}

			// has special routines for handling elastic collisions.
			propagateSolutionParticles( theSimulation, time_step );

			gettimeofday(&tnow, NULL);
			double complex_time_stop = tnow.tv_sec + (1e-6) * tnow.tv_usec;

			last_complex_time += complex_time_stop - complex_time_start;

			/* END complex */

			/* BEGIN SRD */
			
			if( do_srd && ! block.fix_membrane )
			{
				for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
				sRec->theSurface->rebox_system();
				double use_dhull = 0;

				double ncol = srd_i->stream_and_collide( theSimulation, M, mlow, mhigh, use_dhull, cur_t*AKMA_TIME, time_step * AKMA_TIME, time_step_collision * AKMA_TIME  );
			}


			/* END SRD */

			gettimeofday(&tnow,NULL);
			double misc_time_start = tnow.tv_sec + (1e-6)*tnow.tv_usec;
			double T = 0;				
			double srd_T = 0;
			double srd_dof = 0;
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			{
				int nv = sRec->theSurface->nv;
#ifdef PARALLEL
				if( theSimulation->ncomplex == 0 )
				{
					if( sRec->do_gen_q )
					{

					}
					else
					{
						PartialSumVertices(sRec->g,sRec->id);
						PartialSyncVertices(sRec->qdot_temp,sRec->id);
					}
				}
				else
				{
					if( sRec->do_gen_q )
					{
					}
					else
					{
						ParallelSum(sRec->g,sRec->nc);
						ParallelSum( sRec->qdot_temp, 3*sRec->theSurface->nv);
					}
				}
#endif

				memcpy( sRec->next_pp, sRec->pp, sizeof(double) * sRec->NQ );
				if( !block.disable_mesh )
				{
					if( do_bd_membrane || do_ld || (switched) )
					{
						if( sRec->do_gen_q )
							GenQMatVecIncrScale( sRec->next_pp, sRec->pp, sRec->EFFM, -gamma_langevin*AKMA_TIME*time_step );
						else
						{
#ifdef OLD_LANGEVIN
							AltSparseCartMatVecIncrScale( sRec->next_pp, sRec->pp, sRec->EFFM, -gamma_langevin*AKMA_TIME*time_step, theSimulation->alpha, sRec->id );
#else
							for( int xv = 0; xv < 3*sRec->theSurface->nv; xv++ ) // here, gamma_langevin has units per time.
								sRec->next_pp[xv] -= gamma_langevin*AKMA_TIME * time_step * sRec->next_pp[xv];
#endif
						}
						
					}
				}

				// LEAPFROG: increment p by 1/2 eps, we have q(t), p(t), report properties for this state (perform Monte Carlo?)

				if( !block.disable_mesh )
				{

					if( sRec->do_gen_q )
					{
						memset(sRec->del_pp,0,sizeof(double)*sRec->NQ);
						
						for( int Q = 0; Q < sRec->NQ; Q++ )
						for( int v1 = 0; v1 <sRec->theSurface->nv; v1++ )
						{
							sRec->del_pp[Q] += sRec->gen_transform[Q*3*nv+v1*3+0] * -sRec->g[3*v1+0] * AKMA_TIME * time_step/2;
							sRec->del_pp[Q] += sRec->gen_transform[Q*3*nv+v1*3+1] * -sRec->g[3*v1+1] * AKMA_TIME * time_step/2;
							sRec->del_pp[Q] += sRec->gen_transform[Q*3*nv+v1*3+2] * -sRec->g[3*v1+2] * AKMA_TIME * time_step/2;
						}

#ifdef PARALLEL
						ParallelSum( sRec->del_pp, sRec->NQ );
						ParallelBroadcast( sRec->del_pp, sRec->NQ );
#endif
						for( int Q = 0; Q < sRec->NQ; Q++ )
							sRec->next_pp[Q] += sRec->del_pp[Q];
					}
					else
					{
						for( int v1 = 0; v1 < nv; v1++ )
						{
							sRec->next_pp[3*v1+0] += -sRec->g[3*v1+0] * AKMA_TIME * time_step/2;
							sRec->next_pp[3*v1+1] += -sRec->g[3*v1+1] * AKMA_TIME * time_step/2;
							sRec->next_pp[3*v1+2] += -sRec->g[3*v1+2] * AKMA_TIME * time_step/2;
						}
					}	
				}
				memcpy( sRec->pp, sRec->next_pp, sizeof(double) * sRec->NQ );
				memset( sRec->Qdot0_trial, 0, sizeof(double) * sRec->NQ );
				

				if( sRec->do_gen_q )
					GenQMatVecIncrScale( sRec->Qdot0_trial, sRec->pp, sRec->EFFM, 1.0 );
				else
					AltSparseCartMatVecIncrScale( sRec->Qdot0_trial, sRec->pp, sRec->EFFM, 1.0, theSimulation->alpha, sRec->id );
					
				for( int Q = 0; Q < sRec->NQ; Q++ )
				{
					if( o >= nequil && sRec->do_gen_q ) { sRec->av_Q_T[Q] += sRec->pp[Q]*sRec->Qdot0_trial[Q] * 0.5; 
									sRec->nav_Q_T[Q] += 1;
//								printf("av_Q_T[%d] %le\n", Q, av_Q_T[Q]/nav_Q[Q]);
					}
					T += sRec->pp[Q] * sRec->Qdot0_trial[Q] * 0.5;
				}

				

				if( do_srd )
				{
					double *vp = srd_i->vp;

					srd_dof = 3 * srd_i->np;
					for( int p = 0; p < srd_i->np; p++ )
						srd_T += 0.5 * srd_i->mass * ( vp[3*p+0]*vp[3*p+0]+vp[3*p+1]*vp[3*p+1]+vp[3*p+2]*vp[3*p+2] ); 
				}
			}


#ifdef PARALLEL
//			if( ! do_gen_q )
//				ParallelSum(&T,1);
#endif
				
			/*****************************
 * 				Right now, we are on the *mesh* canonical ensemble with {q},{p}
 * 			 *****************************/


			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			{
				if( !block.disable_mesh )
				{
					if( do_bd_membrane || do_ld )
					{
#ifdef OLD_LANGEVIN
						for( int Q = 0; Q < sRec->NQ; Q++ )
							sRec->next_pp[Q] += gsl_ran_gaussian(rng_x, sqrt(2*gamma_langevin*temperature*AKMA_TIME*time_step));
#else
						if( sRec->do_gen_q )
						{
							for( int Q = 0; Q < sRec->NQ; Q++ )
								sRec->next_pp[Q] += gsl_ran_gaussian(rng_x, sqrt(2*gamma_langevin*temperature*AKMA_TIME*time_step));
						}
						else
						{

							for( int Q = 0; Q < 3*sRec->theSurface->nv; Q++ )
								noise_vec[Q] = gsl_ran_gaussian(rng_x, sqrt(2*gamma_langevin*temperature*AKMA_TIME*time_step));

							AltSparseCartMatVecIncrScale( sRec->next_pp, noise_vec, sRec->MMat, 1.0, theSimulation->alpha, sRec->id  );
							
						}
#endif
#ifdef PARALLEL
						ParallelBroadcast(sRec->next_pp,sRec->NQ );
#endif
					}

					// LEAPFROG: increment p by 1/2 eps, we have q(t), p(t+eps/2)
					if( sRec->do_gen_q )
					{
						for( int Q = 0; Q < sRec->NQ; Q++ )
							sRec->next_pp[Q] += sRec->del_pp[Q];
					}
					else
					{
						for( int v1 = 0; v1 < sRec->theSurface->nv; v1++ )
						{
							sRec->next_pp[3*v1+0] += -sRec->g[3*v1+0] * AKMA_TIME * time_step/2;
							sRec->next_pp[3*v1+1] += -sRec->g[3*v1+1] * AKMA_TIME * time_step/2;
							sRec->next_pp[3*v1+2] += -sRec->g[3*v1+2] * AKMA_TIME * time_step/2;
						}	
					}
				}
				memcpy( sRec->pp, sRec->next_pp, sizeof(double) * sRec->NQ );
				memset( sRec->qdot0, 0, sizeof(double) * 3 * sRec->theSurface->nv );

				if( sRec->do_gen_q )
				{	
					memset( sRec->Qdot0, 0, sizeof(double) * sRec->NQ );
					GenQMatVecIncrScale( sRec->Qdot0, sRec->pp, sRec->EFFM, 1.0 );
					MatVec( sRec->gen_transform, sRec->Qdot0, sRec->qdot0, sRec->NQ, 3*sRec->theSurface->nv ); 
				}
				else
					AltSparseCartMatVecIncrScale( sRec->qdot0, sRec->pp, sRec->EFFM, 1.0, theSimulation->alpha, sRec->id );

				memcpy( sRec->qdot, sRec->qdot0, sizeof(double) * 3 * sRec->theSurface->nv );

				if( !sRec->do_gen_q )
					AltSparseCartMatVecIncrScale( sRec->qdot, sRec->qdot_temp, sRec->EFFM, 1.0, theSimulation->alpha, sRec->id  );

				// LEAPFROG: increment q by eps, we have q(t+eps), p(t+eps/2)

				for( int v1 = 0; v1 < sRec->theSurface->nv; v1++ )
				{
					sRec->r[3*v1+0] += sRec->qdot[3*v1+0] * AKMA_TIME * time_step;
					sRec->r[3*v1+1] += sRec->qdot[3*v1+1] * AKMA_TIME * time_step;
					sRec->r[3*v1+2] += sRec->qdot[3*v1+2] * AKMA_TIME * time_step;
				}

				if( sRec->do_gen_q )
				{

					for( int Q = 0; Q < sRec->NQ; Q++ )
					{
						sRec->QV[Q] += sRec->Qdot0[Q] * AKMA_TIME * time_step;
						if( o >= nequil )
						{
							sRec->av_Q[Q] += sRec->QV[Q];
							sRec->av_Q2[Q] += sRec->QV[Q]*sRec->QV[Q];
							sRec->nav_Q[Q] += 1;
						}
					}
				}
#ifdef PARALLEL
				if( sRec->do_gen_q )
				{

				}
				else
					PartialSyncVertices(sRec->pp, sRec->id);
#endif
			}
			// ************* DO REACTION DIFFUSION

			if(do_rd)
			{	
				// get_tracked
				if(debug)
					printf("debug RD 2nd ncomplexes: %d\n", theSimulation->ncomplex);

				rd->do_rd(theSimulation); 

				if(debug)
				{
					for(int p = 0; p < theSimulation->ncomplex; p++)
					{
						printf("debug RD 3rd id: %d ntracked: %d\n", p, rd->tracked[p]->ntracked);
					}
				}
				//run RD
			}

			// ************* END REACTION DIFFUSION
#ifdef PARALLEL
			ParallelSum( &PT, 1 );
#endif
			double dof = 0;

			
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			{
				double ldof = 0;
				ldof = sRec->NQ;
				if( !sRec->do_gen_q )
					ldof = 3*sRec->theSurface->nv-3;
				dof += ldof;
			}
			double mem_T = 2*T / dof;
			
			T += PT;
			T += srd_T;
			dof += srd_dof;
			for( int c = 0; c < theSimulation->ncomplex; c++ )
			{
				if( ! theSimulation->allComplexes[c]->do_bd ) 
				{	
#ifdef DISABLE_ON_MEMBRANE_T
					dof += 3 * (theSimulation->allComplexes[c]->nsites - theSimulation->allComplexes[c]->nattach);
#else
					dof += 3 * theSimulation->allComplexes[c]->nsites - theSimulation->allComplexes[c]->nattach;
#endif
				}
			}
			double TEMP = 2 * T / dof;

			// in kcal/mol

			TEMP /= kcal_mol_K;
			mem_T /= kcal_mol_K;
			sum_average_temp += TEMP;
			n_temp += 1;

			fflush(stdout);
			
			if( t == 0 )
			{
				printf("t: %le ns o: %d T: %.8le V: %.12le T+V: %.14le TEMP: %le MEM_TEMP: %le AV_TEMP %le VR: %.3le VMEM: %le VP: %le", (cur_t * 1e9), o, T, V, T+V, TEMP, mem_T, sum_average_temp / n_temp,VR, VMEM, VP );

				global_delete_this = o;
				if( step_rate > 0 )
					printf(" steps/s: %le", step_rate );
				printf("\n");
				fflush(stdout);
			}
			
			gettimeofday(&tnow,NULL);
			double misc_time_stop = tnow.tv_sec + (1e-6)*tnow.tv_usec;

			last_misc_time += misc_time_stop - misc_time_start;
			switched=0;

			if( block.record_curvature && o >= nequil  )
			{
				for( int c = 0; c < theSimulation->ncomplex; c++ )
					avc[c] += theSimulation->allComplexes[c]->local_curvature( theSimulation );
				navc+=1;
			}

			if( block.s_q && global_cntr % block.s_q_period == 0 && o >= nequil)
			{
				// Update SANS B-histogram.
				theSimulation->sample_B_hist( B_hist, &A2dz2_sampled, SANS_SAMPLE_NRM, 100000, sans_max_r, nsans_bins, block.shape_correction );  
			}

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
		
		if( block.lipid_mc_period > 0 )
		{
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
				sRec->theSurface->measureLipidCurvature(sRec->r, o < nequil );
		}
#if 0 // pre-simulation, change me!
		if( collect_hk )
		{
			theSurface->directFT( lhq, lhq2, r, nmodes, nmodes );

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
#endif
	
#ifdef PARALLEL
		if( theSimulation->ncomplex > 0 ) ParallelSyncComplexes( theSimulation->allComplexes, theSimulation->ncomplex );
#endif

		if( tFile && par_info.my_id == BASE_TASK )
		{
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
				sRec->theSurface->put(sRec->r);
			//theSurface->writeLimitStructure(tFile);
        		theSimulation->writeLimitingSurface(tFile);
			fflush(tFile);
	
			if( do_srd )
				srd_i->writeXYZ(srdXYZFile);
		}	

		if( block.write_alpha_period > 0 && alphaFile )
		{
			fprintf(alphaFile, "%le %le %le\n", theSimulation->alpha[0], theSimulation->alpha[1], theSimulation->alpha[2] );
			fflush(alphaFile);
		}

	
#if 0 // pre-simulation, change me!
		if( block.tachyon && par_info.my_id == BASE_TASK )
		{
			theSurface->writeTachyon( block.jobName, block.tachyon_res, block.tachyon_interp, 
				r, allComplexes,theSimulation->ncomplex, &block, srd_i, block.tachyon_tri_center ); 
		}
#endif			

		
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
			do_bd_membrane = 0;
			do_bd_particles = 0;
			do_srd = 0;
			sum_average_temp = 0;
			n_temp = 0;
			switched = 1;
		}

		if( block.s_q && o >= block.nequil)
		{
			// update SANS-intensity on disk.
				
			char fileName[256];
			sprintf(fileName, "%s.sq", block.jobName );

			writeSq( fileName, B_hist, A2dz2_sampled, sans_max_r, nsans_bins, block.q_min, block.q_max, block.nq, qvals );
			
		}

	}

#ifndef OLD_LANGEVIN
	free(noise_vec);
#endif
	}

#ifdef PARALLEL
	if( theSimulation->ncomplex > 0 ) ParallelSyncComplexes( theSimulation->allComplexes, theSimulation->ncomplex );
#endif

#if 1 // pre-simulation change me


	/*********** COMPUTE ENERGY at SAVE ************/	

	// limit context.
	{
		V=0;
		VR=0;
		double VMEM = 0;
		double VP = 0;
		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			V += sRec->theSurface->energy(sRec->r,NULL);
	
		for( int cx = 0; cx < par_info.nc; cx++ )
		{
			int c = par_info.complexes[cx];
			VP += theSimulation->allComplexes[c]->V(theSimulation);	
			VP += theSimulation->allComplexes[c]->AttachV(theSimulation);	
		}
	
		V += VP;
		V += Boxed_PP_V( theSimulation ); 
		ParallelSum(&V,1);
	
		printf("Energy at save: %.14le\n", V );
		
		if( par_info.my_id == BASE_TASK )
		{
			char fname[256];
		
			sprintf(fname, "%s.save", block.jobName );
	
			FILE *saveFile = fopen( fname, "w");
			theSimulation->saveRestart(saveFile,-1);
			fclose(saveFile);
		}
	}
#endif

	if( block.record_curvature )
	{
		for( int c = 0; c < theSimulation->ncomplex; c++ )
		{
			char *typeName = NULL;
			theSimulation->allComplexes[c]->print_type(&typeName);
			printf("Complex %d type %s average curvature %lf\n", c, typeName, avc[c]/(navc+1e-11) );
			if( typeName) free(typeName);
		}
	}

#if 0 // pre-simulation change me!
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
				printf(" l %d", (int)lround(output_qvals[Q]) ); 
			printf("<h> %le <h^2> %le k_c_apparent %le\n",
					av_Q[Q], av_Q2[Q], kc_app );
		}
		printf("------ done\n");
		
	}
#endif
	if( block.create_all_atom )
	{
		printf("WARNING: creating all atom from the first surface.\n");
		theSimulation->allSurfaces->theSurface->createAllAtom(  &block );
	}
/*	FILE *saveFile = fopen("file.save", "w");

	for( int x = 0; x < theSurface->nv; x++ )
	{
		fprintf(saveFile, "%.14le %.14le %.14le\n",
			r[3*x+0], r[3*x+1], r[3*x+2] );
	}
	fprintf(saveFile, "%.14le %.14le %.14le\n", r[3*nv+0], r[3*nv+1], r[3*nv+2] );

	fclose(saveFile);
*/
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		clearForceSet(sRec->theForceSet);
	if( do_srd )
		srd_i->clear();

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
#endif

	return 0;
}



void updateParticleR( int p, int *pfaces, double *puv, double *p_r_m, double *r, surface *theSurface, int np )
{
	int nv = theSurface->nv;

	double rp[3];
	double nrm[3];
	theSurface->evaluateRNRM( pfaces[p], puv[2*p+0], puv[2*p+1], rp, nrm, r);

	p_r_m[3*p+0] = rp[0];			
	p_r_m[3*p+1] = rp[1];			
	p_r_m[3*p+2] = rp[2];		
	
	theSurface->updateParticle( p_r_m+3*p, p, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
	
	if( do_2p )	
	{
		p_r_m[3*np+3*p+0] = rp[0] + dist_nrm * nrm[0];			
		p_r_m[3*np+3*p+1] = rp[1] + dist_nrm * nrm[1];			
		p_r_m[3*np+3*p+2] = rp[2] + dist_nrm * nrm[2];			
		theSurface->updateParticle( p_r_m+3*np+3*p, np+p, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
	}
}




