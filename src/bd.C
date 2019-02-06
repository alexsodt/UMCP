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
	printf("Syntax: bd [input file] [options ...]\n");
	printf("See input.C for options.\n");

	parameterBlock block;

	int nwarnings = getInput( argv, argc, &block );
	
	srand(block.random_seed);


	/* copy parameters */
		
	int doMonge = block.monge;
	kc = block.kc;
	kg = block.kg;
	KA = block.KA;
	mode_KA = block.mode_KA;
	double particle_density = block.rho;
	double particle_footprint = block.footprint;
	double particle_c0 = block.c0;
	int nsteps = block.nsteps;
	int nequil = block.nequil;
	double dimer_c0 = block.dimer_c0;
	dist_nrm =block.dist_nrm;
	particle_radius  = block.radius1; 
	double particle_radius2 = block.radius2;

	if( particle_radius2 > 1e-7 )
	{
		printf("Second particle enabled.\n");
		do_2p = 1;	
	}

	if( max_rad < particle_radius ) max_rad = particle_radius;
	if( max_rad < particle_radius2 ) max_rad = particle_radius2;

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

	int nc = 3 * sub_surface->nv+3;
	int nv = sub_surface->nv;

	double *r = (double *)malloc( sizeof(double) * nc );
	double *g = (double *)malloc( sizeof(double) * nc );
	sub_surface->get(r);
	

	r[3*nv+0] = 1.0;
	r[3*nv+1] = 1.0;
	r[3*nv+2] = 1.0;
	sub_surface->setg0(r);


/* The surface has been defined at this point */

	for( int t = 0; t < sub_surface->nt; t++ )
	{
		// I stored information about the triangles and information I used to compute properties (like curvature) on the surfaces defined by the triangles in different structures.
		// here is the triangle:
		triangle *tri = sub_surface->theTriangles+t;

		// f is the index into my structures for computing properties.
		int f = tri->f;

		if( f >= sub_surface->nf_faces )
		{
			f -= sub_surface->nf_faces;
			int formulas_per_face = sub_surface->nf_irr_pts;
			int *indices = sub_surface->theIrregularFormulas[f*formulas_per_face].cp;
			int np = sub_surface->theIrregularFormulas[f*formulas_per_face].ncoor;
			int base_vertex = sub_surface->theIrregularFormulas[f*formulas_per_face].vertex;
			
			printf("Triangle %d has irregular base vertex %d. Its points are:", t, base_vertex );
			for( int i = 0; i < np; i++ )
				printf(" %d", indices[i] );
			printf("\n");
		}
		else
		{
			int formulas_per_face = sub_surface->nf_g_q_p;
			int *indices = sub_surface->theFormulas[f*formulas_per_face].cp;
			int np = sub_surface->theFormulas[f*formulas_per_face].ncoor;
			int base_vertex = sub_surface->theFormulas[f*formulas_per_face].vertex;
			
			printf("Triangle %d is regular with base vertex %d. Its points are:", t, base_vertex );
			for( int i = 0; i < np; i++ )
				printf(" %d", indices[i] );
			printf("\n");
		}			
		 
	}











	double *ro = (double *)malloc( sizeof(double) * nc );
	memcpy( ro, r, sizeof(double) * nc );	

	double T = 300;

	double nacc = 0;
	double nrej = 0;
	double macc = 0;
	double mrej = 0;
	
	double n_p_acc = 0;
	double n_p_rej = 0;

	double n_area_acc = 0;
	double n_area_rej = 0;
	
	double n_mode_acc = 0;
	double n_mode_rej = 0;

	double area0;
	double cur_area;
	sub_surface->area(r, -1, &cur_area, &area0 );
	printf("area: %le area0: %le\n", cur_area, area0 );
	int n = 16;
	
	int mode_x = block.mode_x;
	int mode_y = block.mode_y;
	int do_mode = 0;	
	double mode_mag = 1;
	double temperature = 5.92E-01;
	double q2 = 0;

	temperature *= (block.T/298.0);

	if( block.mode_max >= 0 )
	{
		if( block.sphere )
		{	
			// currently only works for sphere
			double u2 = temperature / ( kc * (block.mode_max +2) * (block.mode_max-1) * (block.mode_max*(block.mode_max+1)) );
			mode_mag = sqrt(fabs(u2));
			do_mode = 1;
		}
		else
		{
			double qmax = block.mode_max * 2 * M_PI / sqrt(area0);
			double u2 = temperature / ( area0 * kc * qmax*qmax*qmax*qmax );
			printf("planar u2: %le\n", u2 );
			mode_mag = sqrt(fabs(u2));
			do_mode = 1;
		}
	}
	else if( block.sphere && mode_x >= 0 )
	{
		if( mode_y < 0 && mode_y < -mode_x )
		{
			printf("Spherical harmonic m %d cannot be less than %d.\n", mode_y, -mode_x );
			exit(1);
		}	
		if( mode_y >= 0 && mode_y > mode_x )
		{
			printf("Spherical harmonic m %d cannot be greater than %d.\n", mode_y, mode_x );
			exit(1);
		}	
		// q: 1 / L
		
		do_mode = 1;

		double u2 = temperature / ( kc * (mode_x +2) * (mode_x-1) * (mode_x*(mode_x+1)) );
		mode_mag = sqrt(fabs(u2));
		printf("mode magnitude: %lf\n", mode_mag );
	}
	else if( mode_x > 0 || mode_y >0 )
	{
		if( !(mode_x >= 0 && mode_x < n && mode_y >=0 && mode_y < n) || (mode_x == 0 && mode_y == 0) )
		{	
			printf("Invalid modes: %d and %d must be between 0 and %d, and one must be > 0.\n", mode_x, mode_y, n-1 );
			exit(1);
		}
		
		// q: 1 / L
		
		double qx = mode_x * 2 * M_PI / Lx;
		double qy = mode_y * 2 * M_PI / Ly;
		
		if( mode_x > n/2 ) qx = -(n-mode_x) * 2 * M_PI / Lx;
		if( mode_y > n/2 ) qy = -(n-mode_y) * 2 * M_PI / Ly;

		q2 = qx*qx+qy*qy;

		mode_mag = 2*sqrt( temperature / ( q2 * q2 ) / Lx / Ly / kc ); 
		do_mode = 1;

		if( block.sphere )
		{
			double u2 = temperature / ( kc * (mode_x +2) * (mode_x-1) * (mode_x*(mode_x+1)) );
			mode_mag = sqrt(fabs(u2));
		}
		printf("mode magnitude: %lf\n", mode_mag );
	}

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

	dimerization_radius = block.dimer_radius;
	dimer_eps = block.dimer_eps;

	if( dimer_eps >= 0 )
		dimerization_radius = particle_radius;

	search_rad = particle_radius;
	if( dimerization_radius > search_rad )
		search_rad = dimerization_radius;
	b = particle_radius / pow(2., 1.0/6.0);


	// kinetic MC parameters.
	double time_step = 1/mode_mag*1e-6; //1e-7; // 100 nanonseconds: Dr. Sapp: please set according to good membrane step size (mode_mag)
	

	int n_bd_modes = 0;
	
	double *qvals = NULL;
	
	if( block.sphere )
	{	
		if( block.mode_max >= 2 )
		{
			n_bd_modes = sub_surface->setup_spherical_set( ro, block.mode_min, block.mode_max, &qvals );

			for( int x = 0; x < n_bd_modes; x++ )
				qvals[x] /= sphere_rad;
		}
		else if( block.mode_x >= 0 )
		{
			sub_surface->setup_spherical_perturb( ro, mode_x, mode_y );
			qvals = (double *)malloc( sizeof(double) * 1 );
			qvals[0] = mode_x / sphere_rad;
			n_bd_modes=1;
		}
	}
	else
	{
		if( block.mode_max >= 1 )
		{
			n_bd_modes = sub_surface->setup_planar_set( ro, block.mode_min, block.mode_max, &qvals );
			printf("n_bd_modes: %d\n", n_bd_modes );
		}
		else if( block.mode_x >= 0 )
		{
			sub_surface->setup_mode_perturb( ro, mode_x, mode_y, 16, 16, Lx, Ly );
			qvals    = (double *)malloc( sizeof(double) * 2 );
			qvals[0] = sqrt(mode_x * mode_x+mode_y*mode_y) * 2 * M_PI / Lx;
			qvals[1] = sqrt(mode_x * mode_x+mode_y*mode_y) * 2 * M_PI / Lx;
			n_bd_modes=2;
		}
	}
	

	if( n_bd_modes == 0 )
	{
		printf("This code needs some modes.\n");
		exit(1);
	}


	double *mode_t = (double *)malloc( sizeof(double) * n_bd_modes );
	memset( mode_t, 0, sizeof(double) * n_bd_modes );
	
	double eta = 1.3e-13;	
	
	static gsl_rng * rng_x = NULL;
	static const gsl_rng_type *rng_T = gsl_rng_default;
        rng_x = gsl_rng_alloc(rng_T);
        gsl_rng_env_setup();	
	gsl_rng_set( rng_x, block.random_seed );

	time_step = block.time_step; // one nanosecond.

	double particle_kinetic_mc_move = 2 * sqrt( block.diffc * time_step ); 
	printf("time step: %.10e (seconds) particle_move: %lf\n", time_step, particle_kinetic_mc_move);
	
	// n / cubic angstrom
	double fix_alpha = 0.0;

	if( block.fix_alpha )
		fix_alpha = 1.0;
	
	int np = lround( particle_density * area0 );
	printf("np: %d\n", np  );
	printf("Area per vertex: %lf\n", area0 / sub_surface->nv );
	printf("Using %d particles\n", np );

	global_plist = (int *)malloc( sizeof(int) * (np + (do_2p ? np : 0 ) ) );

	double min_box = 3.0;

	if( particle_radius > min_box )
		min_box = particle_radius;

	sub_surface->setupBoxing( NULL, NULL, np + (do_2p ? np : 0 ) , min_box, 1 );

	int *pfaces = (int *)malloc( sizeof(int) * np );
	int *pleaflet = (int *)malloc( sizeof(int) * np*2 );

	double *puv = (double *)malloc( sizeof(double) * np * 2 );
	double *p_r_m = (double *)malloc( sizeof(double) * 3 * (2*np) );
	double *p_r_n = (double *)malloc( sizeof(double) * 3 * (2*np) );
	double *p_r_off = p_r_m + np*3;
	double *rads = (double *)malloc( sizeof(double) * 2*np);
	double *place_rads = (double*)malloc( sizeof(double) * 2 * np );

	FILE *tFile = NULL;
	if( block.movie )
	{
		FILE *tpsf = NULL;
		tpsf = fopen("traj.psf","w");
		sub_surface->writeLimitPSF(tpsf,np, dist_nrm);
		fclose(tpsf);
	
		tFile = fopen("traj.xyz","w");
	}
	
	// loop through the particles and assure that none are interacting.

	int p = 0;

	int max_particle_placement = 10000 * np;
	int particle_iterations = 0;
	
	// the index of the other pair if part of a dimer.
	int *dimer_p = (int *)malloc( sizeof(int) * np );

	double rad_search = dimerization_radius;
	if( particle_radius > rad_search )
		rad_search = particle_radius;			
	for( int p = 0; p < np; p++ )
	{
		place_rads[p] = rad_search;
		place_rads[p+np] = particle_radius2;
	}
	
	flag = (int *)malloc( sizeof(int) * 2 * np );
	memset( flag, 0, sizeof(int) * 2 * np );


	for( int p = 0; p < np; p++ )
	{
		int f;

		double u;
		double v;

		int mode_qualify = 0;
		double rp[3];
		double nrm[3];
		sub_surface->randomPointOnSurface( &f, &u, &v );



		double rn = rand()/(double)RAND_MAX;
		//pleaflet[p] = (rand() % 2 ? 1 : -1);
		pleaflet[p] = (rn < block.leaflet_fraction ? 1 : -1);
		pleaflet[p+np] = pleaflet[p];
		pfaces[p]  = f;
		puv[2*p]   = u;
		puv[2*p+1] = v;
		

		sub_surface->evaluateRNRM( f, u, v, rp, nrm, r);

		p_r_m[3*p+0] = rp[0];			
		p_r_m[3*p+1] = rp[1];			
		p_r_m[3*p+2] = rp[2];			
		
		p_r_n[3*p+0] = nrm[0];			
		p_r_n[3*p+1] = nrm[1];			
		p_r_n[3*p+2] = nrm[2];			


		rads[p] = particle_radius;	

		if( do_2p )	
		{
			p_r_off[3*p+0] = rp[0] + dist_nrm * nrm[0];			
			p_r_off[3*p+1] = rp[1] + dist_nrm * nrm[1];			
			p_r_off[3*p+2] = rp[2] + dist_nrm * nrm[2];			
			rads[p+np] = particle_radius2;	
		}
			
		sub_surface->addParticle( p_r_m+3*p, p, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
		if( do_2p ) sub_surface->addParticle( p_r_off+3*p, p + np, r[3*nv+0], r[3*nv+1], r[3*nv+2]  ); 
	}
	
	particle_particle_energy( p_r_m, np, (do_2p ? np : 0), NULL, 0,  particle_radius, particle_radius2, r+3*nv, rads, dimer_p, pleaflet ); 
	
	while( p < np && particle_iterations < max_particle_placement)
	{
		double test_e = particle_particle_energy( p_r_m, np, (do_2p ? np : 0), &p, 1, particle_radius, particle_radius2, r+3*nv, rads, dimer_p, pleaflet ); 

		if( test_e > 1e-3 || dimer_p[p] != -1 )
		{
			int f;

			double u;
			double v;
			double rp[3];
			double nrm[3];

			sub_surface->randomPointOnSurface( &f, &u, &v );

			pleaflet[p] = (rand() % 2 ? 1 : -1);
			pleaflet[p+np] = pleaflet[p];
			pfaces[p]  = f;
			puv[2*p]   = u;
			puv[2*p+1] = v;

			sub_surface->evaluateRNRM( f, u, v, rp, nrm, r);

			p_r_m[3*p+0] = rp[0];			
			p_r_m[3*p+1] = rp[1];			
			p_r_m[3*p+2] = rp[2];			
			
			p_r_n[3*p+0] = nrm[0];			
			p_r_n[3*p+1] = nrm[1];			
			p_r_n[3*p+2] = nrm[2];			

			rads[p] = particle_radius;	
	
			if( do_2p )	
			{
				p_r_off[3*p+0] = rp[0] + dist_nrm * nrm[0];			
				p_r_off[3*p+1] = rp[1] + dist_nrm * nrm[1];			
				p_r_off[3*p+2] = rp[2] + dist_nrm * nrm[2];			
				rads[p+np] = particle_radius2;	
			}
						
			updateParticleR( p, pfaces, puv, p_r_m, r, sub_surface, np ); 
		}
		else	
		{
			sub_surface->addParticleToFace( pfaces[p], p, particle_c0 * pleaflet[p], particle_footprint );
			p++;
		}
	}


	debug_trigger = 1;

	if( p != np )
	{
		printf("Failed to place particles on the surface without clashes. Is the density very very high?\n");
		exit(1);
	}

	// We can now assume that *nothing* exists as a complex.
	
	// we will save the state of dimers beforehand.
	int *saved_dimer_p = (int *)malloc( sizeof(int) * np );
	int *dimer_list    = (int *)malloc( sizeof(int) * np );
	int *dimer_spot    = (int *)malloc( sizeof(int) * np );

	for( int p = 0; p < np; p++ )
		dimer_spot[p] = -1;

	int ndimers = 0;

	double *c0_p = (double *)malloc( sizeof(double) * np );

	for( int p = 0; p < np; p++ )
	{
		dimer_p[p] = -1;
		c0_p[p] = particle_c0 * pleaflet[p];
	}

	int nfaces = sub_surface->nf_faces + sub_surface->nf_irr_faces;

	// allocate a large array to potentially hold every face that may need to be recomputed if we break up or form a complex.
	int *modified_face_list = (int *)malloc( sizeof(int) * 3 * nfaces );
	int nmod = 0;

	double *avh = (double *)malloc( sizeof(double) * n * n );
	memset( avh, 0, sizeof(double) * n * n );
	double *lh_big = (double *)malloc( sizeof(double) * n * n*2*2 );
	double *lh = (double *)malloc( sizeof(double) * n * n *2 );
	double *lhq = (double *)malloc( sizeof(double) * n * n*2 );
	double *avq = (double *)malloc( sizeof(double) * n * n*2 );
	memset( avq, 0, sizeof(double) * n * n*2 );
	double ncnt = 0;
	

	double *temp_avq = (double *)malloc( sizeof(double) * n * n*2 );
	memset( temp_avq, 0, sizeof(double) * n * n*2 );
	double temp_ncnt = 0;

	int nbins = 100;
	double histogram[nbins];
	memset( histogram, 0, sizeof(double) * nbins );


	int o_lim = nsteps;
	int o_start = nequil;

	double point_density = (sub_surface->nv) / (Lx*Ly);


	int *plist = (int *)malloc( sizeof(int) * np );
	
	double expec_dimer = 0;
	double nexpec_dimer = 0;
	double expec_e = 0;
	double nexpec_e = 0;
	double *av_sph_h2 = (double *)malloc( sizeof(double) * n_bd_modes );
	memset( av_sph_h2, 0, sizeof(double) * n_bd_modes );
	double *nav_sph_h2 = (double *)malloc( sizeof(double) * n_bd_modes );
	memset( nav_sph_h2, 0, sizeof(double) * n_bd_modes );
	double nav_h2  = 0;
	double av_h2 = 0;
	double mode_s = 0;
	double mode_c = 0;

	
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

		for( int p = 0; p < np; p++ )
		{
			getLine(loadFile, buffer );
			if( feof(loadFile) ) break;
			sscanf(buffer,"%d %lf %lf", pfaces+p, puv+2*p, puv+2*p+1 ); 
	
		}
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


	int *f_test = (int *)malloc (sizeof(int) *  n_test_pts );
	double *uv_test = (double *)malloc( sizeof(double) * n_test_pts*2 );

	for( int px = 0; px < n_test_pts; px++ )
	{
		int f;

		double u;
		double v;

		sub_surface->randomPointOnSurface( &f, &u, &v );

		f_test[px] = f;
		uv_test[2*px+0] = u;
		uv_test[2*px+1] = v;
	}
	

//	int nspace = block.o_lim * 10;
	int nspace = 1;
	int nrm_nspace = 1;
	int rho_nspace = 1;
	rvals = (double *)malloc( sizeof(double) * nspace * n_test_pts );
#ifndef NO_DIFF_R
	double *diff_r = (double *)malloc( sizeof(double) * np * 3 * nspace );	
#endif


	int n_q = block.nq;
	double q_min = block.q_min;
	double q_max = block.q_max;
	
	double *uq = NULL; // monitor one mode's decorrelation time to check.
	double *Sq = NULL;
	double *Sq_t = NULL;

	double *av_full = NULL;
	double *av_ncorr_full = NULL;
	double *av_corr_fn = NULL;
	double *av_ncorr_fn = NULL;
	double *nrm_t = NULL;


	int n_corr_av = 0;	
	int n_ncorr_av = 0;	
	int n_rho_corr_av = 0;	
	double weight_local = 0;

	int l_nmax = 0;

	uq = (double *)malloc( sizeof(double) * nspace );
	int nsq = 0;
	int counter_reset = -1;

	int n_ncorr_pts = 10;
	int *ncorr_f = (int *)malloc( sizeof(int) * n_ncorr_pts );
	double *ncorr_uv = (double *)malloc( sizeof(double) * n_ncorr_pts * 2 );
	
	double *av_rho_corr_fn = NULL;
	double *av_rho_corr_full = NULL;

	for( int p = 0; p < n_ncorr_pts; p++ )
		sub_surface->randomPointOnSurface( ncorr_f+p, ncorr_uv+2*p+0, ncorr_uv+2*p+1 );
	int n_modes_tracked = 0;
	
	if( block.track_rho )
	{
		n_modes_tracked += 1;
	}

	double n_rho_corr_fn = 0;
	double *inst_rho_corr_fn = NULL;
	double *tracked_modes = NULL;

	if( block.s_q || block.nse )
	{	
		Sq = (double *)malloc( sizeof(double) * n_q * 2);
		memset( Sq, 0, sizeof(double) * n_q * 2 );
		if( block.nse )
		{
			Sq_t = (double *)malloc( sizeof(double) * n_q * 2 * nspace );
		}
	
		double n_max = block.max_time / ( time_step * block.kinetic_corr_period );
		counter_reset = lround(n_max);

		if( counter_reset < block.o_lim / block.kinetic_corr_period )
		{
			counter_reset = block.o_lim / block.kinetic_corr_period;

			printf("Resetting counter max to minimum time %le (%d).\n", counter_reset * time_step * block.kinetic_corr_period, counter_reset );
				
		}

		if( n_q * n_max * 8 / (1024*1024) > 1000 )
		{
			printf("WARNING: %lf GB of memory dedicated to saving the correlation function in core.\n",
				n_max * 8 / (1024*1024*1024) );
		}
	
		if( n_max > pow(2,31) )
		{
			printf("The array is too big.\n");
			exit(1);
		}


		av_full = (double*) malloc( sizeof(double) * counter_reset * n_q );
		memset( av_full, 0, sizeof(double) * counter_reset );
 
		av_corr_fn = (double *)malloc( sizeof(double) * counter_reset * n_q);
		memset( av_corr_fn, 0, sizeof(double) * counter_reset * n_q );
	}
	else if( block.ncorr )
	{
		double n_max = block.max_time / ( time_step * block.kinetic_corr_period );
		counter_reset = lround(n_max);

		if( counter_reset < block.o_lim / block.kinetic_corr_period )
		{
			counter_reset = block.o_lim / block.kinetic_corr_period;

			printf("Resetting counter max to minimum time %le (%d).\n", counter_reset * time_step * block.kinetic_corr_period, counter_reset );
				
		}
		if( n_ncorr_pts * n_max * 8 / (1024*1024) > 1000 )
		{
			printf("WARNING: %lf GB of memory dedicated to saving the correlation function in core.\n",
				n_max * 8 / (1024*1024*1024) );
		}
	
		if( n_max > pow(2,31) )
		{
			printf("The array is too big.\n");
			exit(1);
		}
		nrm_t = (double *)malloc( sizeof(double) * n_ncorr_pts * 3 * nrm_nspace );

		av_ncorr_full = (double*) malloc( sizeof(double) * counter_reset*3  );
		memset( av_ncorr_full, 0, sizeof(double) * counter_reset*3  );
 
		av_ncorr_fn = (double *)malloc( sizeof(double) * counter_reset*3);
		memset( av_ncorr_fn, 0, sizeof(double) * counter_reset*3 );
	}
	else if( block.track_rho)
	{
		double n_max = block.max_time / ( time_step * block.kinetic_corr_period );
		counter_reset = lround(n_max);

		if( counter_reset < block.o_lim / block.kinetic_corr_period )
		{
			counter_reset = block.o_lim / block.kinetic_corr_period;

			printf("Resetting counter max to minimum time %le (%d).\n", counter_reset * time_step * block.kinetic_corr_period, counter_reset );
				
		}
		if( n_modes_tracked * n_max * 8 / (1024*1024) > 1000 )
		{
			printf("WARNING: %lf GB of memory dedicated to saving the correlation function in core.\n",
				n_max * 8 / (1024*1024*1024) );
		}
	
		if( n_max > pow(2,31) )
		{
			printf("The array is too big.\n");
			exit(1);
		}
		tracked_modes = (double *)malloc( sizeof(double) * n_modes_tracked * 4 * counter_reset );

		av_rho_corr_full = (double*) malloc( sizeof(double) * counter_reset*4*n_modes_tracked  );
		memset( av_rho_corr_full, 0, sizeof(double) * counter_reset*4*n_modes_tracked );
 
		av_rho_corr_fn = (double *)malloc( sizeof(double) * counter_reset*4*n_modes_tracked);
		memset( av_rho_corr_fn, 0, sizeof(double) * counter_reset*4*n_modes_tracked );
	}

	double *av_t2 = (double *)malloc( sizeof(double) * n_bd_modes );
	memset( av_t2, 0, sizeof(double) * n_bd_modes );
	double nav = 0;
	double *D = (double *)malloc( sizeof(double) * n_bd_modes );
	double *fc = (double *)malloc( sizeof(double) * n_bd_modes );
	int m_cntr = 0;

	double running_fcoef = 0;
	double running_fcoef_2 = 0;
	double n_running = 0;

	double p_l_val = -1;
	for( int m = 0; m < n_bd_modes; m++ )
	{
		double q = qvals[m];

		double quote_solvent_scaling = 1;//1.0/3.0;

		if( block.sphere )
		{
			double l = qvals[m] * sphere_rad;
		
			double Zl = (2*l+1)*(2*l*l+2*l-1) / ( l*(l+1)*(l+2)*(l-1) );


			double TimeScale = eta * sphere_rad * sphere_rad * sphere_rad * Zl / (kc * l*(l+1));
			fc[m] = kc*(l+2)*(l-1)*l*(l+1);
			D[m] = quote_solvent_scaling * temperature / (fc[m]*TimeScale); 
			double time_constant = D[m] * fc[m] / temperature;
		
			if( fabs(p_l_val -qvals[m] *  sphere_rad ) > 1e-4 )
				printf("Mode %d has relaxation time %le nanoseconds. Particle relaxation time %le nanoseconds\n", m,
					(1.0/time_constant)*(1e9), 1.0/(block.diffc * q * q) * (1e9) );
			p_l_val = l;
		}
		else
		{
			double partial_factor = 2;
			fc[m] = kc*Lx*Ly*(q*q*q*q)/partial_factor;
			D[m] = quote_solvent_scaling * partial_factor*temperature / (4 * eta * q * Lx * Ly); 
			double time_constant = D[m] * fc[m] / temperature;
				printf("Mode %d has relaxation time %le nanoseconds. Particle relaxation time %le nanoseconds\n", m,
					(1.0/time_constant)*(1e9), 1.0/(block.diffc * q * q) * (1e9) );
		}

	}

	// brownian dynamics
	
	double *pgrad = (double *)malloc( sizeof(double) * 2 * np + 1 ); // particle uv gradient.
	double *pcgrad = (double *)malloc( sizeof(double) * 3 * np + 1 ); // particle cartesian gradient from any external forces.
	double *mgrad = (double *)malloc( sizeof(double) * n_bd_modes );	
	
	int do_grad = 1;

	if( np == 0 ) //use analytical grad. disable this to compare particle simulations
		do_grad = 0;	

	memset( ppdist, 0, sizeof(double) * N_BINS );

	for( int p = 0; p < np; p++ )
		updateParticleR( p, pfaces, puv, p_r_m, r, sub_surface, np ); 
	double *mode_t_p = (double *)malloc( sizeof(double) * n_bd_modes );
	double dt = block.time_step;
	int corr_step = 0;
	int ncorr_step = 0;
	int rho_corr_step = 0;

	char code = 'A';


	int done = 0;
	int o = 0;

	double cur_t = 0;
					
	double i_fcoef_c = 0;
	double i_fcoef_s = 0;
	for( int i = 0; i < np; i++ )
	{
		i_fcoef_c += 2*cos( p_r_m[3*i+0] * 2*M_PI * block.mode_x / Lx ) / area0; 
		i_fcoef_s += 2*sin( p_r_m[3*i+0] * 2*M_PI * block.mode_x / Lx ) / area0; 
	}

	for( int m = 0; m < n_bd_modes; m++ )
	{
		mode_t_p[m] = mode_t[m];
		sub_surface->mode_perturb( r, mode_t[m], m );
	}

#ifdef DEBUG_MODE_DENSITY_COUPLING
	double min_uq = 0, min_uq_e = 1e100;
	double cur_pert = mode_t[0];
	for( double uq = -10; uq <= 10; uq += 0.1 )
	{
		sub_surface->mode_perturb( r, uq - cur_pert, 0 );
		double en = sub_surface->energy(r,puv,-1,NULL,NULL,doMonge);
		if( en < min_uq_e ) 
		{
			min_uq_e = en;
			min_uq = uq;
		}
		printf("uq: %le en: %le ifcoef_c: %le ifcoef_s: %le\n", uq, en, i_fcoef_c, i_fcoef_s  );
		cur_pert = uq;
	}

	printf("i_f/s_coef: %le %le min_uq: %le\n", i_fcoef_c, i_fcoef_s, min_uq );
	mode_t[0] = cur_pert;
	mode_t_p[0] = mode_t[0];
#endif

	for( int iq = 0; iq < n_modes_tracked; iq++ )
	{
		double q = qvals[iq];
		double D = block.diffc;

		double expec_t = 1.0 / (D * q * q);
		double force_k = (temperature * area0 / particle_density / 2);
		double var = temperature / (force_k);
		 
		printf("Tracked density mode %d q %le diffusion constant %lf expected relaxation time %le expected variance %le\n",
			iq, q, D, expec_t, var );	
		var = 2*temperature / (kc * q*q*q*q * area0);

		double partial_factor = 2;
		double fc = kc*area0*(q*q*q*q)/partial_factor;
		D = partial_factor*temperature / (4 * eta * q * Lx * Ly); 
		double time_constant = D * fc / temperature;

		printf("Tracked membrane mode %d q %le expected relaxation time %le expected variance %le\n",
			iq, q, 1.0 / time_constant, var);	
	}

	while( !done )
	{
//		printf("Current time step %le nanoseconds.\n", corr_step * time_step * block.kinetic_corr_period *(1e9) );
		double membrane_en = sub_surface->energy(r,puv,-1,NULL,NULL,doMonge);
		double particle_en = 0;

		if( !block.non_interacting )
			particle_en = pp_grad( p_r_m, pcgrad, np, r+3*nv, rads, pleaflet ); 
		printf("t: %le ns Vtot: %le Vmem: %le Vpp: %le\n", cur_t*(1e9), membrane_en + particle_en, membrane_en, particle_en );

		fflush(stdout);
		for( int t = 0; t < block.o_lim; t++, cur_t += time_step )
		{
			if( do_grad )
			{
				memcpy( mode_t_p, mode_t, sizeof(double) * n_bd_modes );
				memset( g, 0, sizeof(double) * nc );
				memset(mgrad, 0, sizeof(double) * n_bd_modes );
				memset(pgrad, 0, sizeof(double) * 2 * np );

				// pgrad has the derivative of the particle energy with respect to its movement on the membrane.
				sub_surface->grad( r, g, puv, pgrad );
				double ec = sub_surface->energy( r, puv, -1, NULL, NULL, doMonge );
				// particle particle interactions.
				
				// dV/dR

				// brute force to debug.
				
				memset( pcgrad, 0, sizeof(double ) * 3 * np );

				double ep = 0;
				if( block.non_interacting == 0 )
				{
					ep = pp_grad( p_r_m, pcgrad, np, r+3*nv, rads, pleaflet ); 

					if( ep > 1e2 )
						printf("Ec: %.14le Ep: %.14le\n", ec, ep );

					for( int p = 0; p < np; p++ )
					{
						double de_dnrm[3] = { 0, 0, 0 };
	
						sub_surface->pointGradient( pfaces[p], puv[2*p+0], puv[2*p+1], r, g, pgrad+2*p, pcgrad+3*p, de_dnrm );  					
					}
				}
				fflush(stdout);


				if( do_grad )
					sub_surface->project_grad_to_modes( r, g, mgrad );		
			}		

			for( int mo = 0; mo < n_bd_modes; mo++ )
			{
				double dmo_dt = 0;

				if( do_grad )
				{
					dmo_dt = - mgrad[mo]; 
					if( np == 0 )
						printf(" grad/analytical: %le gcode: %.14le gpaper: %.14le\n", (dmo_dt / (-fc[mo]*mode_t[mo]+1e-20) ), dmo_dt, -fc[mo]*mode_t[mo] );
				}
				else
				{
					dmo_dt = - fc[mo] * mode_t[mo];
				}
				dmo_dt = dmo_dt * D[mo] / temperature + sqrt( 2 * D[mo] / dt ) * gsl_ran_gaussian(rng_x, 1.0 );
				mode_t[mo] +=	dmo_dt * dt;
				av_t2[mo] += mode_t[mo] * mode_t[mo];
			}

			if( do_grad )
			{
				for( int m = 0; m < n_bd_modes; m++ )
					sub_surface->mode_perturb( r, mode_t[m]-mode_t_p[m], m );
			}
			
			for( int p = 0; p < np; p++ )
			{
				int f=pfaces[p];
				double u = puv[2*p+0], v = puv[2*p+1];
				double frc_duv[2] = { -pgrad[2*p+0], -pgrad[2*p+1]};

				
				frc_duv[0] *= block.diffc / temperature;
				frc_duv[1] *= block.diffc / temperature;

				double sigma = sqrt(2.0 * block.diffc);

//				double pe_before = sub_surface->penergy( f, r, puv, 0, 0 ); 
				int of = f;
				double fstep[3]={0,0,0};
				double ouv[2] = {u,v};
				sub_surface->localMove( &f, &u, &v, sigma, r, frc_duv, dt, fstep, 10); 	


				pfaces[p] = f;
				puv[2*p+0] = u;
				puv[2*p+1] = v;
				
				double old_p[3] = { p_r_m[3*p+0], p_r_m[3*p+1], p_r_m[3*p+2] };

				updateParticleR( p, pfaces, puv, p_r_m, r, sub_surface, np ); 

				double pdr[3] = { p_r_m[3*p+0] - old_p[0], p_r_m[3*p+1] - old_p[1], p_r_m[3*p+2] - old_p[2] };
				double dr_squared = pdr[0]*pdr[0]+pdr[1]*pdr[1]+pdr[2]*pdr[2];


				if( f != of )
				{
					sub_surface->removeParticleFromFace( of, p, c0_p[p], particle_footprint );
					sub_surface->addParticleToFace( pfaces[p], p, c0_p[p], particle_footprint );
				}	
			}

			nav++;
			if( t % block.kinetic_corr_period == 0 )
			{
				if( !do_grad )
				{
					memcpy( r, ro, sizeof(double) * (3 * nv + 3) );
	
					for( int m = 0; m < n_bd_modes; m++ )
						sub_surface->mode_perturb( r, mode_t[m], m );
				}

				if( corr_step == nspace )
				{
					nspace *= 2;
					uq = (double *)realloc( uq, sizeof(double) * nspace );
					rvals = (double *)realloc( rvals, sizeof(double) * n_test_pts * nspace );
#ifndef NO_DIFF_R
					diff_r = (double *)realloc( diff_r, sizeof(double) * 3 * np * nspace );
#endif
					if( block.nse )
						Sq_t = (double *)realloc( Sq_t, sizeof(double) * 2 * n_q * nspace );
				}		
				
				if( rho_corr_step == rho_nspace )
				{
					rho_nspace *= 2;
					if( block.track_rho )
						tracked_modes = (double *)realloc( tracked_modes, sizeof(double) * n_modes_tracked * rho_nspace * 4 );	
				}
					
				if( ncorr_step == nrm_nspace )
				{
					nrm_nspace *= 2;
					if( block.ncorr )
						nrm_t = (double *)realloc( nrm_t, sizeof(double) * n_ncorr_pts * nrm_nspace *3 );	
				}
				if( n_bd_modes > 0 )
					uq[corr_step] = mode_t[0];
				else
					uq[corr_step] = 0;


			

				for( int px = 0; px < n_test_pts; px++ )
				{
					double rv[3];
					double nrm[3];
					sub_surface->evaluateRNRM(f_test[px], uv_test[2*px+0], uv_test[2*px+1], rv, nrm, r );

					if( block.sphere )
						rvals[corr_step*n_test_pts+px] = sqrt(rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2] );
					else
						rvals[corr_step*n_test_pts+px] = rv[2];

				}

#ifndef NO_DIFF_R
				for( int p = 0 ; p < np; p++ )
				{	
					double r[3] = { p_r_m[3*p+0], p_r_m[3*p+1], p_r_m[3*p+2] };
		
					if( corr_step > 0 )
					{
						while( r[0] - diff_r[3*p+0] < -Lx/2 ) r[0] += Lx;
						while( r[0] - diff_r[3*p+0] >  Lx/2 ) r[0] -= Lx;
							
						while( r[1] - diff_r[3*p+1] < -Ly/2 ) r[1] += Ly;
						while( r[1] - diff_r[3*p+1] >  Ly/2 ) r[1] -= Ly;
					
						while( r[2] - diff_r[3*p+2] < -Lz/2 ) r[2] += Lz;
						while( r[2] - diff_r[3*p+2] >  Lz/2 ) r[2] -= Lz;
					}	
						
					diff_r[corr_step*3*np+p*3+0] = r[0];
					diff_r[corr_step*3*np+p*3+1] = r[1];
					diff_r[corr_step*3*np+p*3+2] = r[2];
				}
#endif
				if( block.track_rho ) 
				{
					double fcoef_c = 0;
					double fcoef_s = 0;

					for( int i = 0; i < np; i++ )
					{
						fcoef_c += 2*cos( p_r_m[3*i+0] * 2*M_PI * block.mode_x / Lx ) / area0; 
						fcoef_s += 2*sin( p_r_m[3*i+0] * 2*M_PI * block.mode_x / Lx ) / area0; 
					}

					tracked_modes[rho_corr_step*4*n_modes_tracked+0] = fcoef_c;
					tracked_modes[rho_corr_step*4*n_modes_tracked+1] = fcoef_s;
					
					tracked_modes[rho_corr_step*4*n_modes_tracked+2] = mode_t[0];
					tracked_modes[rho_corr_step*4*n_modes_tracked+3] = mode_t[1];

					running_fcoef += fcoef_c;
					running_fcoef_2 += fcoef_c * fcoef_c;
					n_running += 1;

					double var = (running_fcoef_2/n_running) - pow( running_fcoef/n_running,2.0);

//					printf("rho fcoef %le %lf %lf running variance %le\n", cur_t, fcoef_c, fcoef_s, var );
				}

				ncorr_step += 1;	
				corr_step += 1;	
				rho_corr_step += 1;
			}
		}


		
		if( block.track_rho )
		{
/****************/
			fftw_complex *cin = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * corr_step );
			fftw_complex *cou = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * corr_step );
			fftw_plan theForwardPlan = fftw_plan_dft_1d( corr_step, cin, cou, FFTW_FORWARD, FFTW_ESTIMATE );
			fftw_plan theBackwardPlan = fftw_plan_dft_1d( corr_step, cin, cou, FFTW_BACKWARD, FFTW_ESTIMATE );

			double q_test = 0.03;

			double *expec_exp = (double *)malloc( sizeof(double) * corr_step  * 2);
			double *expec_gauss = (double *)malloc( sizeof(double) * corr_step );
			double *expec_fn = (double *)malloc( sizeof(double) * corr_step );
			memset( expec_exp,   0, sizeof(double) * corr_step*2 );
			memset( expec_gauss, 0, sizeof(double) * corr_step );
			memset( expec_fn, 0, sizeof(double) * corr_step );

			int do_eiq = 0;

			for( int t = 0; t < n_test_pts; t++ )
			{
				if( do_eiq )
				{
					for( int t1 = 0; t1 < corr_step; t1++ )
					for( int dt = 0; dt < corr_step-t1; dt++ )
					{
						expec_exp[2*dt+0] += cos( q_test * (rvals[(t1+dt)*n_test_pts+t] - rvals[t1*n_test_pts+t]) );
						expec_exp[2*dt+1] += sin( q_test * (rvals[(t1+dt)*n_test_pts+t] - rvals[t1*n_test_pts+t]) );
	      				
						expec_gauss[dt]   += pow(           rvals[(t1+dt)*n_test_pts+t] - rvals[t1*n_test_pts+t],2.0);
						expec_fn[dt] += 1;
					}
				}

				double av_r2 = 0, av_r = 0;
				for( int i = 0; i < corr_step; i++ )
				{
					av_r2 += rvals[i*n_test_pts+t] * rvals[i*n_test_pts+t];
					av_r +=  rvals[i*n_test_pts+t];
				}

				av_r /= corr_step;
				av_r2 /= corr_step;

				for( int i = 0; i < corr_step; i++)
				{
					cin[i][0] = rvals[i*n_test_pts+t];
					cin[i][1] = 0;
				}
			
				fftw_execute(theForwardPlan);
				double average = 0;
				for( int i = 0; i < corr_step; i++ )
				{
					double r = cou[i][0];
					double c = cou[i][1];
			
					cin[i][0] = r * r + c * c;
					cin[i][1] = 0;
				}
			
				fftw_execute(theBackwardPlan);

			}
		}



		if( o % 30 == 0 && block.mode_max > 0)
		{
			printf("kc:");
			int mo = 0;
			for( int l = block.mode_min; l <= block.mode_max; l++ )
			{
				for( int m = -l; m <= l; m++, mo++)
				{
					printf(" %lf", temperature / ( (av_t2[mo]/nav) * (l+2)*(l-1)*l*(l+1) ) );
				}
			}	
			printf("\n");
		}
		if( o % 10 == 0 && block.nse)
		{
			if( code == 'A' ) code = 'B';
			else if( code == 'B' ) code = 'A';
			printf("Writing MSD file.\n");
			char fileName[256];
			sprintf(fileName, "%s_msd_%c.txt", block.jobName, code );
			FILE *msdFile = fopen(fileName,"w");
			double *r2_hist = (double *)malloc( sizeof(double) * corr_step );
			double *fft_r2_hist = (double *)malloc( sizeof(double) * corr_step );
			double *n_r2_hist = (double *)malloc( sizeof(double) * corr_step );
			memset( r2_hist, 0, sizeof(double) * corr_step );
			memset( fft_r2_hist, 0, sizeof(double) * corr_step );
			memset( n_r2_hist, 0, sizeof(double) * corr_step );

/****************/
			fftw_complex *cin = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * corr_step );
			fftw_complex *cou = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * corr_step );
			fftw_plan theForwardPlan = fftw_plan_dft_1d( corr_step, cin, cou, FFTW_FORWARD, FFTW_ESTIMATE );
			fftw_plan theBackwardPlan = fftw_plan_dft_1d( corr_step, cin, cou, FFTW_BACKWARD, FFTW_ESTIMATE );

			double q_test = 0.03;

			double *expec_exp = (double *)malloc( sizeof(double) * corr_step  * 2);
			double *expec_gauss = (double *)malloc( sizeof(double) * corr_step );
			double *expec_fn = (double *)malloc( sizeof(double) * corr_step );
			memset( expec_exp,   0, sizeof(double) * corr_step*2 );
			memset( expec_gauss, 0, sizeof(double) * corr_step );
			memset( expec_fn, 0, sizeof(double) * corr_step );

			int do_eiq = 0;

			for( int t = 0; t < n_test_pts; t++ )
			{
				if( do_eiq )
				{
					for( int t1 = 0; t1 < corr_step; t1++ )
					for( int dt = 0; dt < corr_step-t1; dt++ )
					{
						expec_exp[2*dt+0] += cos( q_test * (rvals[(t1+dt)*n_test_pts+t] - rvals[t1*n_test_pts+t]) );
						expec_exp[2*dt+1] += sin( q_test * (rvals[(t1+dt)*n_test_pts+t] - rvals[t1*n_test_pts+t]) );
	      				
						expec_gauss[dt]   += pow(           rvals[(t1+dt)*n_test_pts+t] - rvals[t1*n_test_pts+t],2.0);
						expec_fn[dt] += 1;
					}
				}

				double av_r2 = 0, av_r = 0;
				for( int i = 0; i < corr_step; i++ )
				{
					av_r2 += rvals[i*n_test_pts+t] * rvals[i*n_test_pts+t];
					av_r +=  rvals[i*n_test_pts+t];
				}

				av_r /= corr_step;
				av_r2 /= corr_step;

				for( int i = 0; i < corr_step; i++)
				{
					cin[i][0] = rvals[i*n_test_pts+t];
					cin[i][1] = 0;
				}
			
				fftw_execute(theForwardPlan);
				double average = 0;
				for( int i = 0; i < corr_step; i++ )
				{
					double r = cou[i][0];
					double c = cou[i][1];
			
					cin[i][0] = r * r + c * c;
					cin[i][1] = 0;
				}
			
				fftw_execute(theBackwardPlan);

				for( int i = 0; i < corr_step; i++ )
				{
					fft_r2_hist[i] += (-2 * cou[i][0]/(double)corr_step/(double)corr_step + 2 * cou[0][0]/(double)corr_step/(double)corr_step) / n_test_pts ;
				}
			}
/****************/

			if( do_eiq )
			{
				for( int dt = 0; dt < corr_step; dt++ )
				{
					printf("monkey %le %le %le %le %le\n", dt * time_step * block.kinetic_corr_period, 
						expec_exp[2*dt+0]/expec_fn[dt],  
						expec_exp[2*dt+1]/expec_fn[dt],  
						exp( - q_test*q_test/2.0 * expec_gauss[dt] / expec_fn[dt] ),
						exp( - q_test*q_test/2.0 * fft_r2_hist[dt])  );
				}
			}

			for( int dt =0; dt < corr_step; dt++ )	
			{
			
				fprintf(msdFile, "%le %le\n", dt *time_step * block.kinetic_corr_period, fft_r2_hist[dt] ); 
			}
			free(r2_hist);
			free(fft_r2_hist);
			free(n_r2_hist);
			fclose(msdFile);


		
		}

		if( block.track_rho )
		{
			char fileName[256];

			double *all_corr = (double *)malloc( sizeof(double)  * rho_corr_step * n_modes_tracked * 4 );
			
			fftw_complex *cin = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * rho_corr_step );
			fftw_complex *cou = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * rho_corr_step );
			fftw_plan theForwardPlan = fftw_plan_dft_1d( rho_corr_step, cin, cou, FFTW_FORWARD, FFTW_ESTIMATE );
			fftw_plan theBackwardPlan = fftw_plan_dft_1d( rho_corr_step, cin, cou, FFTW_BACKWARD, FFTW_ESTIMATE );

	
			sprintf(fileName, "%s_temp_mode.txt", block.jobName );
			FILE *tempModes = fopen(fileName,"w");
			sprintf(fileName, "%s_tracked_modes.txt", block.jobName );
			FILE *trackedModes = fopen(fileName,"w");
			
			for( int i = 0; i < rho_corr_step; i++)
			{
				for( int ix = 0; ix < 4; ix++ )
				{
					for( int iq = 0; iq < n_modes_tracked; iq++ )
					{	
						fprintf(tempModes, " %le", tracked_modes[i*n_modes_tracked*4+iq*4+ix] );
					}
				}
				fprintf(tempModes,"\n");
			}
			fclose(tempModes);

			for( int ix = 0; ix < 4; ix++ )
			{
				for( int iq = 0; iq < n_modes_tracked; iq++ )
				{	
					for( int i = 0; i < rho_corr_step; i++)
					{
						cin[i][0] = tracked_modes[i*n_modes_tracked*4+iq*4+ix];
						cin[i][1] = 0;
					}
			
					fftw_execute(theForwardPlan);
					double average = 0;
					for( int i = 0; i < rho_corr_step; i++ )
					{
						double r = cou[i][0];
						double c = cou[i][1];
			
						cin[i][0] = r * r + c * c;
						cin[i][1] = 0;
			
						average += r*r+c*c;
					}
	
					average /= rho_corr_step;	
					average /= rho_corr_step;
			
					fftw_execute(theBackwardPlan);
	
					for( int i = 0; i < rho_corr_step; i++ )
					{
						all_corr[iq*rho_corr_step*4+ix*rho_corr_step+i] = cou[i][0];
						all_corr[iq*rho_corr_step*4+ix*rho_corr_step+i] /= rho_corr_step;
						all_corr[iq*rho_corr_step*4+ix*rho_corr_step+i] /= rho_corr_step;
					//	printf("%d %d %d %le\n", ix, iq, i,  all_corr[iq*rho_corr_step*+ix*rho_corr_step+i] );
					}
				}
			}

			fprintf(trackedModes,"t c s uq_c, uq_s");
			fprintf(trackedModes,"\n");

			double local_weight = (double)rho_corr_step / (double)counter_reset;

		/// print variance, I was debugging this. seems good.			
	//		printf("Var: ");
			for( int it = 0; it < counter_reset; it++ )
			{
				fprintf(trackedModes,"%le", it * time_step * block.kinetic_corr_period ); 
					
				if( it >= rho_corr_step && n_rho_corr_av == 0 )
					break;
	
				for( int iq = 0; iq < n_modes_tracked; iq++ )
				for( int ix = 0; ix < 4; ix++ )
				{
					double val_tot = 0;

					if( it < rho_corr_step )
					{
	
						val_tot += local_weight * (all_corr[(iq*4+ix)*rho_corr_step+it]) ;
					}	
					else
					{

					}
						
					val_tot += av_rho_corr_fn[(iq*4+ix)];
					val_tot /= ( n_rho_corr_av + local_weight);
// print variance
//					if( it == 0 )
//						printf(" %le", val_tot );

					fprintf(trackedModes," %le", val_tot );
				}

		// print variance
		//		if( it == 0 )
		//			printf("\n" );
				fprintf(trackedModes,"\n");
			}
					
	
			if( rho_corr_step >= counter_reset )
			{
				for( int iq = 0; iq < n_modes_tracked; iq++ )
				for( int ix = 0; ix < 4; ix++ )
				{
					double val0 = all_corr[(iq*4+ix)*rho_corr_step+0];

					av_rho_corr_full[0] += val0;

					for( int i = 0; i < rho_corr_step; i++ )
						av_rho_corr_fn[(iq*4+ix)*rho_corr_step+i] += all_corr[(iq*4+ix)*rho_corr_step+i];
					
				}

				printf("Resetting correlation counter.\n");
				n_rho_corr_av += 1;
				rho_corr_step = 0;
			}

			fclose(trackedModes );

			free(all_corr);
			fftw_destroy_plan(theForwardPlan);
			fftw_destroy_plan(theBackwardPlan);
			fftw_free(cin);
			fftw_free(cou);
		}
		
		if( block.ncorr )
		{
			char fileName[256];

			double *all_corr = (double *)malloc( sizeof(double)  * ncorr_step * n_ncorr_pts * 3 );
			
			fftw_complex *cin = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * ncorr_step );
			fftw_complex *cou = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * ncorr_step );
			fftw_plan theForwardPlan = fftw_plan_dft_1d( ncorr_step, cin, cou, FFTW_FORWARD, FFTW_ESTIMATE );
			fftw_plan theBackwardPlan = fftw_plan_dft_1d( ncorr_step, cin, cou, FFTW_BACKWARD, FFTW_ESTIMATE );

	
			sprintf(fileName, "%s_ncorr.txt", block.jobName );
			FILE *nrm_times = fopen(fileName,"w");
		


			for( int ix = 0; ix < 3; ix++ )
			{
				for( int iq = 0; iq < n_ncorr_pts; iq++ )
				{	
					for( int i = 0; i < ncorr_step; i++)
					{
						cin[i][0] = nrm_t[i*n_ncorr_pts*3+iq*3+ix];
						cin[i][1] = 0;
						//printf("%d %lf\n", i, nrm_t[i*n_ncorr_pts+iq] );
					}
			
					fftw_execute(theForwardPlan);
					double average = 0;
					for( int i = 0; i < ncorr_step; i++ )
					{
						double r = cou[i][0];
						double c = cou[i][1];
			
						cin[i][0] = r * r + c * c;
						cin[i][1] = 0;
			
						average += r*r+c*c;
					}
	
					average /= ncorr_step;	
					average /= ncorr_step;
			
					fftw_execute(theBackwardPlan);
	
					for( int i = 0; i < ncorr_step; i++ )
					{
						all_corr[iq*ncorr_step*3+i*3+ix] = cou[i][0];
						all_corr[iq*ncorr_step*3+i*3+ix] /= ncorr_step;
						all_corr[iq*ncorr_step*3+i*3+ix] /= ncorr_step;
					}
				}
			}

			fprintf(nrm_times,"t xx yy zz ");
			fprintf(nrm_times,"\n");

			double local_weight = (double)ncorr_step / (double)counter_reset;
			for( int it = 0; it < counter_reset; it++ )
			{
				fprintf(nrm_times,"%le", it * time_step * block.kinetic_corr_period ); 
					
				if( it >= ncorr_step && n_ncorr_av == 0 )
					break;
	
				for( int ix = 0; ix < 3; ix++ )
				{
					double val_tot = 0;
					if( it < ncorr_step )
					{
	
						for( int iq = 0; iq < n_ncorr_pts; iq++ )
							val_tot += local_weight * (all_corr[(iq*ncorr_step+it)*3+ix]) / n_ncorr_pts;
					}	
					else
					{

					}
						
					val_tot += av_ncorr_fn[it*3+ix];
					val_tot /= ( n_ncorr_av + local_weight);

					fprintf(nrm_times," %le", val_tot );
				}

				fprintf(nrm_times,"\n");
			}

					
	
			if( ncorr_step >= counter_reset )
			{
				for( int iq = 0; iq < n_ncorr_pts; iq++ )
				{
					double val0 = all_corr[iq*ncorr_step+0];

					av_ncorr_full[0] += val0/n_ncorr_pts;

					for( int i = 0; i < ncorr_step; i++ )
					for( int ix = 0; ix < 3; ix++ )
						av_ncorr_fn[i*3+ix] += all_corr[(iq*ncorr_step+i)*3+ix] / n_ncorr_pts;
					
				}

				printf("Resetting correlation counter.\n");
				n_ncorr_av += 1;
				ncorr_step = 0;
			}

			fclose(nrm_times);

			free(all_corr);
			fftw_destroy_plan(theForwardPlan);
			fftw_destroy_plan(theBackwardPlan);
			fftw_free(cin);
			fftw_free(cou);
		}



		if( tFile )
		{
			sub_surface->put(r);
			sub_surface->writeLimitStructure(tFile,pfaces,puv,np,dist_nrm);
			fflush(tFile);
		}
		if( o % 2 == 0 )
		{
			char fileName[256];

			double *all_corr = (double *)malloc( sizeof(double) * n_q * corr_step );
			
			fftw_complex *cin = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * corr_step );
			fftw_complex *cou = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * corr_step );
			fftw_plan theForwardPlan = fftw_plan_dft_1d( corr_step, cin, cou, FFTW_FORWARD, FFTW_ESTIMATE );
			fftw_plan theBackwardPlan = fftw_plan_dft_1d( corr_step, cin, cou, FFTW_BACKWARD, FFTW_ESTIMATE );

			if(block.nse )
			{
				printf("Writing correlation file.\n");
	
				sprintf(fileName, "%s_sq_corr.txt", block.jobName );
				FILE *Sq_times = fopen(fileName,"w");
			
				
				double *av_sq_r = (double *)malloc( sizeof(double) * n_q );
				double *av_sq_c = (double *)malloc( sizeof(double) * n_q );
				
				for( int iq = 0; iq < n_q; iq++ )
				{
					av_sq_r[iq] = 0;
					av_sq_c[iq] = 0;
				
					for( int i = 0; i < corr_step; i++)
					{
						av_sq_r[iq] += Sq_t[i*2*n_q+2*iq+0];
						av_sq_c[iq] += Sq_t[i*2*n_q+2*iq+1];
					}
			
					av_sq_r[iq] /= (double)corr_step;
					av_sq_c[iq] /= (double)corr_step;
				}
			
				for( int iq = 0; iq < n_q; iq++ )
				{	
					for( int i = 0; i < corr_step; i++)
					{
						cin[i][0] = Sq_t[i*2*n_q+2*iq+0];// - av_sq_r[iq];
						cin[i][1] = Sq_t[i*2*n_q+2*iq+1];// - av_sq_c[iq];
				//		printf("1 %d %le %le\n", i, cin[i][0], cin[i][1] );
					}
			
					fftw_execute(theForwardPlan);
					double average = 0;
					for( int i = 0; i < corr_step; i++ )
					{
						double r = cou[i][0];
						double c = cou[i][1];
			
						cin[i][0] = r * r + c * c;
						cin[i][1] = 0;
				//		printf("2 %d %le %le\n", i, r, c );
			
						average += r*r+c*c;
					}
					average /= corr_step;	
					average /= corr_step;
				//	printf("Average: %le\n", average );
			
					fftw_execute(theBackwardPlan);
					for( int i = 0; i < corr_step; i++ )
					{
				//		printf("3 %d %le %le\n", i, cou[i][0], cou[i][1] );
						all_corr[iq*corr_step+i] = cou[i][0];
						all_corr[iq*corr_step+i] /= corr_step;
						all_corr[iq*corr_step+i] /= corr_step;
					}
				}
		
	


				fprintf(Sq_times,"t");
				if( n_q > 1 )
				{
					for( int iq = 0; iq < n_q; iq++ )
						fprintf(Sq_times, " %lf", q_min + iq*(q_max-q_min)/(n_q-1) );
				}
				else
				{
					fprintf(Sq_times, " %lf", (q_max+q_min)/2 );
				}
				fprintf(Sq_times,"\n");

				double local_weight = (double)corr_step / (double)counter_reset;

				for( int iq = 0; iq < n_q; iq++ )
				{
					double val0 = all_corr[iq*corr_step+0];
//					fprintf(Sq_times, " %lf", val0);
					fprintf(Sq_times, " %lf", (av_full[iq]+local_weight*val0) / (n_corr_av+local_weight)  );
				}
				fprintf(Sq_times,"\n");

				for( int it = 0; it < corr_step; it++ )
				{
					fprintf(Sq_times,"%le", it * time_step * block.kinetic_corr_period ); 
					for( int iq = 0; iq < n_q; iq++ )
					{
						double val0 = all_corr[iq*corr_step+0];
//						fprintf(Sq_times," %le", all_corr[iq*corr_step+it]/val0 );
						fprintf(Sq_times," %le", (av_corr_fn[it*n_q+iq]+local_weight*all_corr[iq*corr_step+it]/val0) / (n_corr_av+local_weight) );
					}
				//	for( int iq = 0; iq < n_q; iq++ )
				//		fprintf(Sq_times," %le %le", Sq_t[it*2*n_q+2*iq+0], 
				//					     Sq_t[it*2*n_q+2*iq+1] );
					fprintf(Sq_times,"\n");
					
				}
				
				if( corr_step >= counter_reset )
				{
					for( int iq = 0; iq < n_q; iq++ )
					{
						double val0 = all_corr[iq*corr_step+0];
	
						av_full[iq] += val0;
	
						for( int i = 0; i < corr_step; i++ )
							av_corr_fn[i*n_q+iq] += all_corr[iq*corr_step+i] / val0;
						
					}

					printf("Resetting correlation counter.\n");
					n_corr_av += 1;
					corr_step = 0;
				}

				free(av_sq_r);
				free(av_sq_c);
				fclose(Sq_times);
			}

			

			free(all_corr);
			fftw_destroy_plan(theForwardPlan);
			fftw_destroy_plan(theBackwardPlan);
			fftw_free(cin);
			fftw_free(cou);
		}
				
		o++;
		
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

	for( int p = 0; p < np; p++ )
		fprintf(saveFile, "%d %lf %lf\n", pfaces[p], puv[2*p+0], puv[2*p+1] );
	fclose(saveFile);


}

//#define DEBUG_NEAR 

double particle_particle_energy( double *r, int np_1, int np_2, int *set, int nset,  double rad1, double rad2, double *alpha, double *rads, int *dimer_p, int * pleaflet)
{
	double La = sub_surface->PBC_vec[0][0]*alpha[0];
	double Lb = sub_surface->PBC_vec[1][1]*alpha[1];
	double Lc = sub_surface->PBC_vec[2][2]*alpha[2];
	double eps = 1e4;
	double en = 0;

	if( set  )
	  {
	    int nnear = 0;
	    for( int t = 0; t < nset; t++ )
	      {
		flag[set[t]] = 1;

		if( np_2 > 0 )
		  flag[set[t]+np_1] = 1;
			
		// form new pairs with everything.	
		if( dimer_p[set[t]] != -1 ) // undimerize.
			dimer_p[dimer_p[set[t]]] = -1;
		dimer_p[set[t]] = -1;
	      }

	    for( int t = 0; t < nset; t++ )
	      {
		int select = set[t];

		int nnear = sub_surface->getNear( r+select *3, select, r, search_rad*2,    alpha[0], alpha[1], alpha[2], global_plist, rads );

		for( int x = 0; x < nnear; x++ )
		  {
			if( pleaflet[select] != pleaflet[global_plist[x]] ) continue;

		    double dr[3] = {
		      r[select*3+0] - r[global_plist[x]*3+0],
		      r[select*3+1] - r[global_plist[x]*3+1],
		      r[select*3+2] - r[global_plist[x]*3+2] };
		    while( dr[0] < -La/2 ) dr[0] += La;
		    while( dr[0] >  La/2 ) dr[0] -= La;
		    while( dr[1] < -Lb/2 ) dr[1] += Lb;
		    while( dr[1] >  Lb/2 ) dr[1] -= Lb;
		    while( dr[2] < -Lc/2 ) dr[2] += Lc;
		    while( dr[2] >  Lc/2 ) dr[2] -= Lc;

		    double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( rval < rads[select] + rads[global_plist[x]] )
			{
				if( flag[global_plist[x]] )
				  en += 0.5 * eps;
				else
				  en += eps;
			}
			if( global_plist[x] >= np_1 ) continue;

			if( rval < 2*dimerization_radius )
			{ // if the new particle is part of a dimer already, DO NOT SET. here, we set:
			 	if( dimer_p[global_plist[x]] == -1 && dimer_p[select] == -1 )
				{
					dimer_p[select] = global_plist[x];
					dimer_p[global_plist[x]] = select;
				}
				else if( dimer_p[select] != global_plist[x] )
					en += eps;
			}

			if( rval < 2*dimerization_radius )
			{
				if( flag[global_plist[x]] )
				  en += 0.5 * dimer_eps;
				else
				  en += dimer_eps;
			}

		  }
	
#ifdef DEBUG_NEAR
			static int dbg_cntr = 0;
			int alt_list[np_1+np_2];
			int nalt = 0;
	
			for( int p = 0; p < np_1+np_2; p++ )
			{
				if( p == select ) continue;
				if( pleaflet[p] != pleaflet[select] ) continue;
				double dr[3] = {
					 r[select*3+0] - r[3*p+0],
					 r[select*3+1] - r[3*p+1],
					 r[select*3+2] - r[3*p+2] };
	
				while( dr[0] < -La/2 ) dr[0] += La;
				while( dr[0] >  La/2 ) dr[0] -= La;
				while( dr[1] < -Lb/2 ) dr[1] += Lb;
				while( dr[1] >  Lb/2 ) dr[1] -= Lb;
				while( dr[2] < -Lc/2 ) dr[2] += Lc;
				while( dr[2] >  Lc/2 ) dr[2] -= Lc;
	
				double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	
				if( rval < (rads[select]+rads[p]) )
				{
					alt_list[nalt] = p;
					nalt++;
				}	
			}
			if( nalt != nnear )
			{
				printf("ERROR #\n");
				exit(1);
			}
		
			for( int x = 0; x < nnear; x++ )
			{
				int got = 0;
	
				for( int p = 0; p < nalt; p++ )
				{
					if( global_plist[x] == alt_list[p] )
						got = 1;
				}
	
				if(!got)
				{
					printf("FIND ERROR\n");
					exit(1);
				}
			}
#endif	
			int nnear_2 = 0;
		
	
			if( np_2 > 0 )
			{
				nnear_2    = sub_surface->getNear( r+ 3*np_1 + select *3, select+np_1, r, 2*search_rad, 	alpha[0], alpha[1], alpha[2], global_plist, rads );	
			
				for( int x = 0; x < nnear_2; x++ )
				{
					if( pleaflet[np_1+select] != pleaflet[global_plist[x]] ) continue;
					double dr[3] = {
					r[(np_1+select)*3+0] - r[global_plist[x]*3+0],
		      			r[(np_1+select)*3+1] - r[global_plist[x]*3+1],
					r[(np_1+select)*3+2] - r[global_plist[x]*3+2] };
					    while( dr[0] < -La/2 ) dr[0] += La;
					    while( dr[0] >  La/2 ) dr[0] -= La;
					    while( dr[1] < -Lb/2 ) dr[1] += Lb;
					    while( dr[1] >  Lb/2 ) dr[1] -= Lb;
					    while( dr[2] < -Lc/2 ) dr[2] += Lc;
					    while( dr[2] >  Lc/2 ) dr[2] -= Lc;
	
					    double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

					if( rval < rads[np_1+select] + rads[global_plist[x]] )
					{
						if( flag[global_plist[x]] )
							en += 0.5 * eps;
						else
							en += eps;
					}
				}
			}			
#ifdef DEBUG_NEAR
			dbg_cntr++;
			nalt = 0;
			for( int p = 0; p < np_1+np_2; p++ )
			{
				if( p == select+np_1 ) continue;
				if( pleaflet[p] != pleaflet[select+np_1] ) continue;
				double dr[3] = {
					 r[(np_1+select)*3+0] - r[3*p+0],
					 r[(np_1+select)*3+1] - r[3*p+1],
					 r[(np_1+select)*3+2] - r[3*p+2] };
	
				while( dr[0] < -La/2 ) dr[0] += La;
				while( dr[0] >  La/2 ) dr[0] -= La;
				while( dr[1] < -Lb/2 ) dr[1] += Lb;
				while( dr[1] >  Lb/2 ) dr[1] -= Lb;
				while( dr[2] < -Lc/2 ) dr[2] += Lc;
				while( dr[2] >  Lc/2 ) dr[2] -= Lc;
	
				double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	
				if( rval < (rads[np_1+select]+rads[p]) )
				{
					alt_list[nalt] = p;
					nalt++;
				}	
			}
			if( nalt != nnear_2 )
			{
				printf("ERROR #\n");
				exit(1);
			}
		
			for( int x = 0; x < nnear_2; x++ )
			{
				int got = 0;
	
				for( int p = 0; p < nalt; p++ )
				{
					if( global_plist[x] == alt_list[p] )
						got = 1;
				}
	
				if(!got)
				{
					printf("FIND ERROR\n");
					exit(1);
				}
			}
#endif		
		}
		for( int t = 0; t < nset; t++ )
		{
			flag[set[t]] = 0;

			if( np_2 > 0 )
				flag[set[t]+np_1] = 0;
		}

//		printf("PPE: %lf\n", nnear * eps );
		return en;	
	}

	for( int p = 0; p < np_1; p++ )
		dimer_p[p] = -1;

	for( int p = 0; p < np_1; p++ )
	{
		int nnear = sub_surface->getNear( r+p *3, p, r, 2*search_rad, 	alpha[0], alpha[1], alpha[2], global_plist, rads);	
		for( int x = 0; x < nnear; x++ )
		{			
			if( pleaflet[p] != pleaflet[global_plist[x]] ) continue;
			double dr[3] = {
				r[p*3+0] - r[3*global_plist[x]+0],
				r[p*3+1] - r[3*global_plist[x]+1],
				r[p*3+2] - r[3*global_plist[x]+2] };

			while( dr[0] < -La/2 ) dr[0] += La;
			while( dr[0] >  La/2 ) dr[0] -= La;
			while( dr[1] < -Lb/2 ) dr[1] += Lb;
			while( dr[1] >  Lb/2 ) dr[1] -= Lb;
			while( dr[2] < -Lc/2 ) dr[2] += Lc;
			while( dr[2] >  Lc/2 ) dr[2] -= Lc;
	
			double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( rval < (rads[p]+rads[global_plist[x]]) )
			{
				int p1 = global_plist[x];
				en += 0.5 * eps;
			}
			// only accept the membrane-bound particle for dimerization.
			if( global_plist[x] > np_1 )
				continue;
			
			if( rval < 2*dimerization_radius )
				  en += 0.5*dimer_eps;

			if( rval < 2*dimerization_radius )
			{
				// they have dimerized.
				if( dimer_p[global_plist[x]] == -1 && dimer_p[p] == -1 )
				{
					dimer_p[p] = global_plist[x];
					dimer_p[global_plist[x]] = p;
				}
				else if( dimer_p[p] != global_plist[x] )
				{
					en += 0.5*eps;
				}
			}
			
		}
	}	
	for( int p = 0; p < np_2; p++ )
	{
		int nnear_2 = sub_surface->getNear( r+ 3*np_1 + p *3, p+np_1, r, 2*max_rad, 	alpha[0], alpha[1], alpha[2], global_plist, rads  );	
		
		for( int x = 0; x < nnear_2; x++ )
		{			
			if( pleaflet[np_1+p] != pleaflet[global_plist[x]] ) continue;
			double dr[3] = {
				r[(np_1+p)*3+0] - r[3*global_plist[x]+0],
				r[(np_1+p)*3+1] - r[3*global_plist[x]+1],
				r[(np_1+p)*3+2] - r[3*global_plist[x]+2] };

			while( dr[0] < -La/2 ) dr[0] += La;
			while( dr[0] >  La/2 ) dr[0] -= La;
			while( dr[1] < -Lb/2 ) dr[1] += Lb;
			while( dr[1] >  Lb/2 ) dr[1] -= Lb;
			while( dr[2] < -Lc/2 ) dr[2] += Lc;
			while( dr[2] >  Lc/2 ) dr[2] -= Lc;
				
	
			double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( rval < rads[np_1+p] + rads[global_plist[x]] )
			{
				  en += 0.5 * eps;
			}
		}
	}

//#define DEBUG_WHOLE	
#ifdef DEBUG_WHOLE
	if( debug_trigger  )
	{
	double dbg_en = 0;
	for( int p = 0; p < np_1; p++ )
		dimer_p[p] = -1;
	for( int p1 = 0; p1 < np_1+np_2; p1++ )
	{
		for( int p = 0; p < np_1+np_2; p++ )
		{
			if( p == p1 ) continue;
			if( pleaflet[p] != pleaflet[p1] ) continue;
			double dr[3] = {
				 r[(p1)*3+0] - r[3*p+0],
				 r[(p1)*3+1] - r[3*p+1],
				 r[(p1)*3+2] - r[3*p+2] };
	
			while( dr[0] < -La/2 ) dr[0] += La;
			while( dr[0] >  La/2 ) dr[0] -= La;
			while( dr[1] < -Lb/2 ) dr[1] += Lb;
			while( dr[1] >  Lb/2 ) dr[1] -= Lb;
			while( dr[2] < -Lc/2 ) dr[2] += Lc;
			while( dr[2] >  Lc/2 ) dr[2] -= Lc;
	
			double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
			
			if( rval < (rads[p1]+rads[p]) )
			{
				dbg_en += 0.5 * eps;
			}
			if( p1 >= np_1 || p >= np_1 )
				continue;
			if( rval < 2*dimerization_radius )
			{ // if the new particle is part of a dimer already, DO NOT SET. here, we set:
			 	if( dimer_p[p1] == -1 && dimer_p[p] == -1 )
				{
					dimer_p[p] = p1;
					dimer_p[p1] = p;
				}
				else if( dimer_p[p1] != p )
				{
					dbg_en += 0.5*eps;
				}
			}

			if( rval < 2*dimerization_radius )
			{
				  dbg_en += 0.5 * dimer_eps;
			}
		}
	}
	if( fabs(dbg_en-en) > 1e-7 )
	{
		printf("DBG energy check problem %lf vs %lf\n", en, dbg_en);
		exit(1);
	}
	}
#endif
	return en;
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

double pp_grad( double *r, double *gc, int np, double *alpha, double *rads, int * pleaflet)
{
	double La = sub_surface->PBC_vec[0][0]*alpha[0];
	double Lb = sub_surface->PBC_vec[1][1]*alpha[1];
	double Lc = sub_surface->PBC_vec[2][2]*alpha[2];
	double en = 0;



	// cartesian gradient.
	memset( gc, 0, sizeof(double) * np * 3 );

//#define DEBUG_BOXING

#ifdef DEBUG_BOXING
	double epass[2] = {0,0};
	for( int pass = 0; pass < 2; pass++ )
	{
		en = 0;
#endif
	for( int p = 0; p < np; p++ )
	{
		int nnear = sub_surface->getNear( r+p*3, p, r, 2*search_rad, 	alpha[0], alpha[1], alpha[2], global_plist, rads);	
#ifdef DEBUG_BOXING
		if( pass == 1 )
		{
			nnear = np;

			for( int xp =0; xp< np; xp++ )
				global_plist[xp] = xp;
		}
#endif

		for( int x = 0; x < nnear; x++ )
		{			
			if( pleaflet[p] != pleaflet[global_plist[x]] ) continue;
			int p2 = global_plist[x];

			if( p >= p2 ) continue;

			double dr[3] = {
				r[p*3+0] - r[3*global_plist[x]+0],
				r[p*3+1] - r[3*global_plist[x]+1],
				r[p*3+2] - r[3*global_plist[x]+2] };

			while( dr[0] < -La/2 ) dr[0] += La;
			while( dr[0] >  La/2 ) dr[0] -= La;
			while( dr[1] < -Lb/2 ) dr[1] += Lb;
			while( dr[1] >  Lb/2 ) dr[1] -= Lb;
			while( dr[2] < -Lc/2 ) dr[2] += Lc;
			while( dr[2] >  Lc/2 ) dr[2] -= Lc;
	
			double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( rval < 2 * particle_radius )
			{
				int b = rval / (particle_radius) * N_BINS;
				if( b < 0 ) b = 0;
				if( b < N_BINS )
					ppdist[b] += 1;

			//	printf("particles %d and %d have distance %le\n", p, p2, rval );
			}

			if( rval > dimerization_radius  )
				continue;

			if( rval < particle_radius )		
			{
#ifdef WCA		
				double l = (b/rval);
				double l2 = l*l;
				double l6 = l2*l2*l2;
				double l12 = l6*l6;
	
				en += 1 + 4 * ( l12 - l6 ) + dimer_eps;	
				double dvdr = -48 * l12 / rval + 24 * l6 / rval;
	

				// trying the negative to debug.
				gc[3*p+0] += dvdr * (dr[0]/rval);	
				gc[3*p+1] += dvdr * (dr[1]/rval);	
				gc[3*p+2] += dvdr * (dr[2]/rval);	
				
				gc[3*p2+0] -= dvdr * (dr[0]/rval);	
				gc[3*p2+1] -= dvdr * (dr[1]/rval);	
				gc[3*p2+2] -= dvdr * (dr[2]/rval);	
#else
				double kval = 3;
				en += (kval/2) * (rval - particle_radius) * (rval-particle_radius) + dimer_eps;
				double dvdr = kval * (rval - particle_radius);
//				en += 1 + 4 * ( l12 - l6 ) + dimer_eps;	
//				double dvdr = -48 * l12 / rval + 24 * l6 / rval;

				// trying the negative to debug.
				gc[3*p+0] += dvdr * (dr[0]/rval);	
				gc[3*p+1] += dvdr * (dr[1]/rval);	
				gc[3*p+2] += dvdr * (dr[2]/rval);	
				
				gc[3*p2+0] -= dvdr * (dr[0]/rval);	
				gc[3*p2+1] -= dvdr * (dr[1]/rval);	
				gc[3*p2+2] -= dvdr * (dr[2]/rval);	


#endif
			}
			else
			{
				double wc = dimerization_radius-particle_radius;
				double c = cos( M_PI * (rval - particle_radius) / (2 * wc ));
				double s = sin( M_PI * (rval - particle_radius) / (    wc ));
				en += dimer_eps * c*c;
				double dvdr = dimer_eps * M_PI / (2*wc);
	
				gc[3*p+0] += dvdr * (dr[0]/rval);	
				gc[3*p+1] += dvdr * (dr[1]/rval);	
				gc[3*p+2] += dvdr * (dr[2]/rval);	
				
				gc[3*p2+0] -= dvdr * (dr[0]/rval);	
				gc[3*p2+1] -= dvdr * (dr[1]/rval);	
				gc[3*p2+2] -= dvdr * (dr[2]/rval);	
			}
		}
	}
#ifdef DEBUG_BOXING
		epass[pass]  = en;
	}
	double de = epass[0]-epass[1];
	if( fabs(de) > 1e-8 ) 
	{
		printf("ERROR.\n");
		printf("Boxing energy: %.14le Full energy: %.14le\n", epass[0], epass[1] );
		
				
		int tp = 218; 				
		int nnear = sub_surface->getNear( r+tp*3, tp, r, 2*search_rad, 	alpha[0], alpha[1], alpha[2], global_plist, rads);	

		exit(1);
	}
	else if( epass[0] > 1e2 )
	{
		printf("BIG BOXING ENERGY.\n");
	}
#endif	
//	printf("En: %le\n", en );
	return en;
}



