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
#include "m_triangles.h"
#include "input.h"
#include <assert.h>
#include "random_global.h"
#include "fpr_subroutines/fpr.h"
#include "fpr_subroutines/Faddeeva.hh"

#define NUMERICAL_INTEGRATION
//#define DO_HISTOGRAM
//#define MOVE_MEMBRANE

double rtol = 1E-10;
double kr_on = 1e11;
double kr_off = 1e4;
double binding_radius = 0;
double sp_ap_displacement = 10.0;

#define STATE_MASK			(1+2+4)
#define STATE_REACTING			0
#define STATE_IN_COMPLEX		1	
#define STATE_FREE			2
#define REACTED_BIT			(1<<7)

typedef struct tracked_pair
{
	int surface_id;	
	int aqueous_id;

	double prev_norm;
	double ps_prev;
	double sep;
	double prev_sep;
	double weight;
	double on_rate;
	double off_rate;
	int flag;
	int state;
	struct tracked_pair *next;
} tracked_pair;

extern double kc;
extern double kg;
extern double KA;
extern double KA4;
double defaultR = 250;
void updateParticleR( int p, int *pfaces, double *puv, double *p_r_m, double *r, surface *sub_surface, int np );

double particle_radius = 0;
double b;
int *global_plist = NULL;
double *global_rads = NULL;
int get_pair_list( surface *sub_surface, double *r, double *aqueous_r, int **pairs, int *npairsSpace, int surface_np, int aqueous_np, double cutoff );

int main( int argc, const char **argv )
{
	char buffer[4096];

	printf("No options are required:\n");
	printf("Syntax: rd_vesicle [input file] [options ...]\n");
	printf("See input.C for options.\n");

	parameterBlock block;
	
	setDefaults(&block);	
	
	const char *defaultMesh = "sphere.mesh";
	free(block.meshName);
	block.meshName = (char *)malloc( sizeof(char) * (strlen(defaultMesh)+1) );
	sprintf(block.meshName, "%s", defaultMesh );

	block.sphere = 1;
	block.mode_max = 6;
//	block.fix_alpha = 1;
	block.sigma = sp_ap_displacement;

	int nwarnings = getInput( argv, argc, &block );

	sp_ap_displacement = block.sigma;
	kr_on = block.k_on;
	kr_off = block.k_off;

	printf("On-rate: %le Off-rate: %le\n", kr_on, kr_off );

//	srand(1);	
	srand(block.random_seed);

	/* copy parameters */
		
	kc = block.kc;
	kg = block.kg;
	KA = block.KA;


	int N_BINS_PHIST = 100;
	double *sum_prob_distance = (double *)malloc( sizeof(double) * N_BINS_PHIST );
	double *n_prob_distance = (double *)malloc( sizeof(double) * N_BINS_PHIST );
	memset( sum_prob_distance, 0, sizeof(double) * N_BINS_PHIST );
	memset( n_prob_distance, 0, sizeof(double) * N_BINS_PHIST );

	double particle_density = block.rho;
	double particle_footprint = block.footprint;
	double particle_c0 = block.c0;
	int nsteps = block.nsteps;
	int nequil = block.nequil;
	particle_radius  = block.radius1; 

	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( block.meshName , 0. );
	surface *sub_surface = theSurface;

	double Lx = sub_surface->PBC_vec[0][0];
	double Ly = sub_surface->PBC_vec[1][1];
	double Lz = sub_surface->PBC_vec[2][2];

	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 );

	double *M[3] = { M5, M6, M7 };
	int mlow = 5;
	int mhigh = 7;
	sub_surface->generateSubdivisionMatrices( M, mlow, mhigh );
	sub_surface->generatePlan();
	sub_surface->box_system();
	
	int nc = 3 * sub_surface->nv+3;
	int nv = sub_surface->nv;
	double *r = (double *)malloc( sizeof(double) * nc );
	sub_surface->get(r);
	r[3*nv+0] = 1.0;
	r[3*nv+1] = 1.0;
	r[3*nv+2] = 1.0;
	sub_surface->setg0(r);
	double *g = (double *)malloc( sizeof(double) * nc );
	
	double area0;
	double cur_area;
	sub_surface->area(r, -1, &cur_area, &area0 );
	printf("area: %le area0: %le\n", cur_area, area0 );
	
	force_set *theForceSet = (force_set *)malloc( sizeof(force_set) );

	int plim = 5;
	int ppf = 0;
	for( int fi = 0; fi <= plim; fi++ )
	for( int fj = 0; fj <= plim-fi; fj++ )
	{
		double f1 = fi / (double)plim;
		double f2 = fj / (double)plim;

		ppf++;
	}

	theForceSet->npts = ppf * sub_surface->nt;
	theForceSet->face_set = (int *)malloc( sizeof(int) * theForceSet->npts );
	theForceSet->uv_set = (double *)malloc( sizeof(double) * theForceSet->npts *2 );
	theForceSet->mass = (double *)malloc( sizeof(double) * theForceSet->npts );

	int totp = 0;
	
	double mass_per_lipid = 760.09 / 6.022e23 / 1000; // POPC kg, wikipedia
	double area_per_lipid = 65.35; // A^2, POPC, interpolated from Kucerka 2011 BBA 2761
	// factor of two for leaflets.
	double mass_per_pt = 2 * area0 / (sub_surface->nf_faces + sub_surface->nf_irr_faces) / ppf * mass_per_lipid / area_per_lipid; 
	// in kg


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

#ifdef TEST
	sub_surface->load_least_squares_fitting( theForceSet );
	sub_surface->test_force_set( theForceSet );

	exit(1);	
#endif
	double LA = sub_surface->PBC_vec[0][0];
	double LB = sub_surface->PBC_vec[1][1];
	double LC = sub_surface->PBC_vec[2][2];


	



	double *ro = (double *)malloc( sizeof(double) * nc );
	memcpy( ro, r, sizeof(double) * nc );	

	double T = 300;

	
	double temperature = 5.92E-01;

	double V0 = 0;
	
	double sphere_rad = 1.0;

	if( block.sphere )
	{	
		V0 =fabs(sub_surface->volume( r));
		sphere_rad = sqrt(area0/(4*M_PI));
			printf("Area: %lf (R: %lf)\n", area0, sqrt(area0/(4*M_PI))  );
			printf("Volume: %lf (R: %lf)\n", V0, pow( V0/(4*M_PI/3.0), 1.0/3.0) );
	}
	else
	{
		double tvec[3];
		cross( sub_surface->PBC_vec[0], sub_surface->PBC_vec[1], tvec );

		V0 = 
			sub_surface->PBC_vec[2][0] * tvec[0] + 
			sub_surface->PBC_vec[2][1] * tvec[1] + 
			sub_surface->PBC_vec[2][2] * tvec[2];
	}
	b = particle_radius / pow(2., 1.0/6.0);

	char fileName[1024];
#ifdef DO_HISTOGRAM
	sprintf(fileName, "%s_histo.txt", block.jobName );
	FILE *histo_file = fopen(fileName,"w");
#endif

	int n_bd_modes = 0;
	
	double *qvals;
	n_bd_modes = sub_surface->setup_spherical_set( ro, block.mode_min, block.mode_max, &qvals );

	double *mode_t = (double *)malloc( sizeof(double) * n_bd_modes );
	memset( mode_t, 0, sizeof(double) * n_bd_modes );
	
	double eta = 1.3e-13;	

	if( !rng_x ) init_random(block.random_seed);

	double time_step = block.time_step; // one nanosecond.

	// n / cubic angstrom
	double fix_alpha = 1.0;
	

	printf("V0: %le conc %le\n", V0, block.concentration );
	int aqueous_np = lround( block.concentration * V0 );
	printf("area0: %le particle_density: %le\n", area0, particle_density );
	int surface_np = lround( particle_density * area0 );
	int surface_np_meanfield = surface_np;
	int *mean_field_activated = NULL;

	if( block.mean_field )
	{
		printf("Mean field surface np: %d Explicit aqueous particles: %d\n", surface_np, aqueous_np  );
		surface_np = aqueous_np; // a placeholder for each aqueous particle.	
		mean_field_activated = (int *)malloc( sizeof(int) * surface_np );
		memset( mean_field_activated, 0, sizeof(int) * surface_np );
	}
	else
		printf("surface np: %d aqueous: %d\n", surface_np, aqueous_np  );

	global_plist = (int *)malloc( sizeof(int) * (surface_np + aqueous_np) );
	global_rads = (double *)malloc( sizeof(double) * (surface_np + aqueous_np) );
	
	for( int x = 0; x < surface_np+aqueous_np; x++ )
		global_rads[x] = 1.0;

	double min_box = 3.0;

	if( particle_radius > min_box )
		min_box = particle_radius;
	
	binding_radius = sp_ap_displacement;
	double Rmax = binding_radius + 3 * sqrt(2 * block.aqueous_diffc * time_step); 

	sub_surface->setupBoxing( NULL, NULL, surface_np+aqueous_np , Rmax, 1 );

	printf("SURFACE NP %d AQUEOUS NP %d\n", surface_np, aqueous_np );

	int *pfaces = (int *)malloc( sizeof(int) * surface_np );
	int *pleaflet = (int *)malloc( sizeof(int) * surface_np );

	double *puv = (double *)malloc( sizeof(double) * surface_np * 2 );
	double *surface_p_r_m = (double *)malloc( sizeof(double) * 3 * (surface_np+aqueous_np) );
	double *surface_p_r_n = (double *)malloc( sizeof(double) * 3 * (surface_np) );
	double *aqueous_r = surface_p_r_m + 3 * surface_np;
	double *aqueous_r_last = (double *)malloc( sizeof(double) * 3 * aqueous_np );
	double *distance = (double *)malloc( sizeof(double) * aqueous_np );

	int *surface_status = (int *)malloc( sizeof(int) * surface_np );
	int *aqueous_status = (int *)malloc( sizeof(int) * aqueous_np );


	tracked_pair **paired_complex = (tracked_pair **)malloc( sizeof(tracked_pair *) * aqueous_np );

	tracked_pair *tracked_pairs = NULL;


	FILE *tFile = NULL;
	FILE *pFile = NULL;

	if( block.movie )
	{
		FILE *tpsf = NULL;
		tpsf = fopen("traj.psf","w");
        	sub_surface->writeLimitingSurfacePSF(tpsf);
		fclose(tpsf);
	
		tFile = fopen("traj.xyz","w");
		pFile = fopen("pts.xyz","w");
	}
	
	// loop through the particles and assure that none are interacting.



	double *c0_p = (double *)malloc( sizeof(double) * surface_np );

	for( int p = 0; p < surface_np; p++ )
		c0_p[p] = particle_c0 * pleaflet[p];

	int nfaces = sub_surface->nf_faces + sub_surface->nf_irr_faces;

	// allocate a large array to potentially hold every face that may need to be recomputed if we break up or form a complex.
	int *modified_face_list = (int *)malloc( sizeof(int) * 3 * nfaces );
	int nmod = 0;


	int o_lim = nsteps;
	int o_start = nequil;

	int *plist = (int *)malloc( sizeof(int) * surface_np );
	

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

		for( int p = 0; p < surface_np; p++ )
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
		o_lim = 2e9;		
	
	int n_test_pts = 1;
	double *rvals = NULL;

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
			{
				printf("Mode %d has relaxation time %le nanoseconds. Particle relaxation time %le nanoseconds\n", m,
					(1.0/time_constant)*(1e9), 1.0/(block.diffc * q * q) * (1e9) );
			}
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
	
	double *pgrad = (double *)malloc( sizeof(double) * (2 * surface_np + 1) ); // particle uv gradient.
	double *pcgrad = (double *)malloc( sizeof(double) * (3 * surface_np + 1) ); // particle cartesian gradient from any external forces.
	double *mgrad = (double *)malloc( sizeof(double) * n_bd_modes );	
					
	memset( g, 0, sizeof(double) * nc );
	memset(mgrad, 0, sizeof(double) * n_bd_modes );
	memset(pgrad, 0, sizeof(double) * (2 * surface_np+1) );
	memset(pcgrad, 0, sizeof(double) * (3 * surface_np+1) );
	
	int move_membrane = 1;
#ifndef MOVE_MEMBRANE
	move_membrane = 0;
#endif
	int do_grad = 1;

	if( surface_np == 0 ) //use analytical grad. disable this to compare particle simulations
		do_grad = 0;	

	double *mode_t_p = (double *)malloc( sizeof(double) * n_bd_modes );
	double dt = block.time_step;
	int corr_step = 0;
	int ncorr_step = 0;
	int rho_corr_step = 0;

	char code = 'A';


	int done = 0;
	int o = 0;

					
	for( int m = 0; m < n_bd_modes; m++ )
	{
		mode_t_p[m] = mode_t[m];
		sub_surface->mode_perturb( r, mode_t[m], m );
	}

	double membrane_fudge = 0;
	double nmembrane_fudge = 0;

	printf("binding_radius: %le Rmax: %le\n", binding_radius, Rmax );
	int *pair_list = (int *)malloc( sizeof(int) * surface_np * 2 );
	int npairsSpace = surface_np;
	int npairs = 0;

	double tot_on = 0;
	double tot_off = 0;

	int *pfaces_cache = (int *)malloc( sizeof(int) * surface_np );
	double *puv_cache = (double *)malloc( sizeof(double) * surface_np*2 );
	double *aqueous_r_cache = (double *)malloc( sizeof(double) * aqueous_np*3 );

	int ncomplex_size = 2000;
	double *ncomplex = (double *)malloc( sizeof(double) * ncomplex_size );
	memset( ncomplex, 0, sizeof(double) * ncomplex_size);
	

	double t_per_complex = time_step; 

	int nruns = block.nruns;

	for( int run = 0; run < nruns; run++ )
	{
		int ncind = 0;

		char fileName[256];
		sprintf(fileName, "%s_ncomplex.txt", block.jobName );
		FILE *n_complex_file = fopen(fileName,"w");

		tracked_pair *next =NULL;
		for( tracked_pair *apair = tracked_pairs; apair; apair = next )
		{
			next = apair->next;
			free(apair);
		}	
		tracked_pairs = NULL;
		

		for( int sp = 0; sp < surface_np; sp++ )
			surface_status[sp] = STATE_FREE;
		for( int ap = 0; ap < aqueous_np; ap++ )
		{
			aqueous_status[ap] = STATE_FREE;
			if( block.mean_field )
				mean_field_activated[ap] = 0;
		}

		for( int p = 0; p < surface_np; p++ )
		{
			if( run != 0 )
				sub_surface->removeParticleFromFace( pfaces[p], p, c0_p[p], particle_footprint );

			int f;
	
			double u;
			double v;
	
			int mode_qualify = 0;
			double rp[3];
			double nrm[3];
			sub_surface->randomPointOnSurface( &f, &u, &v );
	
			pleaflet[p] = 1;
			pfaces[p]  = f;
			puv[2*p]   = u;
			puv[2*p+1] = v;
			
	
			sub_surface->evaluateRNRM( f, u, v, rp, nrm, r);
	
			surface_p_r_m[3*p+0] = rp[0];			
			surface_p_r_m[3*p+1] = rp[1];			
			surface_p_r_m[3*p+2] = rp[2];			
			
			surface_p_r_n[3*p+0] = nrm[0];			
			surface_p_r_n[3*p+1] = nrm[1];			
			surface_p_r_n[3*p+2] = nrm[2];			
	
	
			if( run == 0 )	
				sub_surface->addParticle( surface_p_r_m+3*p, p, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
		}
		
		
		for( int p = 0; p < surface_np; p++ )
			sub_surface->addParticleToFace( pfaces[p], p, particle_c0 * pleaflet[p], particle_footprint );
		// initialize 
	
		for( int p = 0; p < surface_np; p++ )
			updateParticleR( p, pfaces, puv, surface_p_r_m, r, sub_surface, surface_np ); 
	
		for (int i = 0; i < aqueous_np; i++) {
			double radius;
			radius = defaultR;

			int outside = 1;

			while( outside )
			{
				aqueous_r_last[3*i+0] = aqueous_r[3*i+0] = -LA/2 + LA * rand() / (double)RAND_MAX;	
				aqueous_r_last[3*i+1] = aqueous_r[3*i+1] = -LB/2 + LB * rand() / (double)RAND_MAX;	
				aqueous_r_last[3*i+2] = aqueous_r[3*i+2] = -LC/2 + LC * rand() / (double)RAND_MAX;	
			
		
				int f;
				double u,v;
				
				if( sub_surface->withinBoxedSurface( aqueous_r+3*i, &f, &u, &v, M, mlow, mhigh, distance+i, -1 ) )
				{
					outside = 0; 
			
					distance[i] -= 1;
	
					if( distance[i] < 1e-3 )
						distance[i] = 1e-3;
				}
				
			}
			
			if( run == 0 )
				sub_surface->addParticle( aqueous_r+3*i, surface_np+i, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
			else
				sub_surface->updateParticle( aqueous_r+3*i, surface_np+i, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
		}

		double cur_t = 0;
		int run_done = 0;
		while( !run_done )
		{
	//		printf("Current time step %le nanoseconds.\n", corr_step * time_step * block.kinetic_corr_period *(1e9) );
			double membrane_en = sub_surface->energy(r,puv,-1,NULL,NULL,0);
			double particle_en = 0;
					
			printf("t: %le ns Vtot: %le Vmem: %le Vpp: %le\n", cur_t*(1e9), membrane_en + particle_en, membrane_en, particle_en );
	
			fflush(stdout);
			for( int t = 0; t < block.o_lim; t++, cur_t += time_step )
			{
				for( int sp = 0; sp < surface_np; sp++ )
					surface_status[sp] -= (surface_status[sp] & REACTED_BIT);
				for( int ap = 0; ap < aqueous_np; ap++ )
					aqueous_status[ap] -= (aqueous_status[ap] & REACTED_BIT);
	
				double Dtot = ( (2.0/3.0) * block.diffc + block.aqueous_diffc);
				double kdiff = 4 * M_PI * Dtot * binding_radius;
				double kact = kr_on;
				double fact = 1.0 + kact / kdiff;
				double alpha = fact * sqrt(Dtot) / binding_radius;
				
				if( block.mean_field )
				{
					for( int ap = 0; ap < aqueous_np; ap++ )
					{
						// are we within range of the membrane?
						// distance[ap] is a lower bound on the minimum distance to the membrane.
	
	
						if( mean_field_activated[ap] )
						{
							double p = 1 - exp( -kr_off * time_step );
						
							if( (double)rand() / (double)RAND_MAX < p )
							{
								aqueous_status[ap] = STATE_REACTING + REACTED_BIT;	
								mean_field_activated[ap] = 0;
								tot_off += 1;
							}		
						}
						else if( distance[ap] < Rmax )
						{
							// alt, find all positions on membrane.
#ifdef NUMERICAL_INTEGRATION							
							int *f_list = NULL;
							double *puv_list = NULL;
							double *areas = NULL;
							int npts = 0;

							// returns a list of triangles that are near the particle.
							sub_surface->assembleNearList( aqueous_r+3*ap, &f_list, &puv_list, &areas, &npts, M, mlow, mhigh, Rmax, 4 );
							if( npts > 0 )
							{
							double *probs = (double *)malloc( sizeof(double) * (npts+1) );

							double running_prob = 0;

							probs[0] = 0;

							for( int dA = 0; dA < npts; dA++ )
							{

								double uv0[3] = { puv_list[6*dA+0], puv_list[6*dA+1],0 };
								double uv1[3] = { puv_list[6*dA+2], puv_list[6*dA+3],0 };
								double uv2[3] = { puv_list[6*dA+4], puv_list[6*dA+5],0 };

								double dudv = triangle_area( uv0, uv1, uv2 );

								double midp[2] = { 
									(uv0[0]+uv1[0]+uv2[0])/3.0,
									(uv0[1]+uv1[1]+uv2[1])/3.0 };
								// the surface metric at the middle of the triangle.
								double g = sub_surface->g(f_list[dA], 
									midp[0],
									midp[1], r );

								
					
								double rmid[3], nmid[3];

								sub_surface->evaluateRNRM( f_list[dA], midp[0], midp[1], rmid, nmid, r );

								double dr[3] = { aqueous_r[3*ap+0] - rmid[0],
										 aqueous_r[3*ap+1] - rmid[1],
										 aqueous_r[3*ap+2] - rmid[2] };
								double sep = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
								// the main integrand.	
								double val = particle_density * passocF(sep, time_step, Dtot, binding_radius, alpha, kact/(kact+kdiff));

								probs[dA+1] = running_prob + val * g * dudv; // zeroeth order integration.
								running_prob = probs[dA+1];

							}
						

#endif

							// evaluate the real distance to the membrane.
	
							int sf;
							double su,sv;
							sub_surface->nearPointOnBoxedSurface( aqueous_r+3*ap, &sf, &su, &sv, M, mlow, mhigh, distance+ap, -1 );
	
							double c_tot = sub_surface->c(sf, su, sv, r );
	
							// updated distance may be greater.
							if( distance[ap] > Rmax )
								continue;
	
							double probvec1 = 0;
	
							double dz = distance[ap];
	
#ifdef EXPR_1	
							if( dz < binding_radius )
							{
								probvec1 = -1 + exp( time_step * (kact+kdiff) * (kact+kdiff) / (16*Dtot*M_PI*M_PI*binding_radius*binding_radius*binding_radius*binding_radius) ) * erfc( time_step * (kact+kdiff) / (4*M_PI * sqrt(Dtot*time_step) * binding_radius * binding_radius) );
								probvec1 *= 2 * Dtot * M_PI * sqrt(M_PI) * binding_radius*binding_radius;
								probvec1 += sqrt(time_step*Dtot) * (kact+kdiff);
								probvec1 = probvec1 * 4 * particle_density * kact * sqrt(M_PI) * binding_radius / (kact+kdiff) / (kact+kdiff);
							}
							else
							{
								probvec1 = sqrt(Dtot)*exp(alpha*(alpha*sqrt(Dtot)*time_step-binding_radius+dz)/sqrt(Dtot))*erfc((2*alpha*sqrt(Dtot)*time_step-binding_radius+dz)/(2*sqrt(Dtot*time_step)))/alpha;
								probvec1 = probvec1 + (binding_radius-dz-sqrt(Dtot)/alpha)*erfc((dz-binding_radius)/(2*sqrt(Dtot*time_step)));
								probvec1 = probvec1 + 2*sqrt(Dtot*time_step)*exp(-(binding_radius-dz)*(binding_radius-dz)/(4*Dtot*time_step))/sqrt(M_PI);
								probvec1 = 2*M_PI*binding_radius*kact/(kact+kdiff)*probvec1*particle_density;
							}
#else
							// the upper indefinite limit.
							probvec1 = binding_radius - sqrt(Dtot)/alpha;
							
							if( dz < binding_radius )
							{
								probvec1 += 2 * sqrt(Dtot * time_step) / sqrt(M_PI);	
								probvec1 -= binding_radius;
								probvec1 += sqrt(Dtot) * Faddeeva::erfcx( sqrt(time_step) * alpha ) / alpha;
							}
							else
							{ // use zero as lower limit.
								probvec1 += 2 * exp( -(dz-binding_radius) * (dz-binding_radius) / (4*Dtot*time_step)) * sqrt(Dtot*time_step) / sqrt(M_PI);
								probvec1 -= dz;
								probvec1 += (sqrt(Dtot) + alpha * (dz-binding_radius)) * erf((dz-binding_radius)/(2*sqrt(Dtot*time_step)) ) / alpha;
								double erfc_arg =  (dz+2*sqrt(Dtot) * time_step * alpha-binding_radius)/(2*sqrt(Dtot*time_step));
								probvec1    += sqrt(Dtot)*exp(alpha*(dz+sqrt(Dtot)*time_step*alpha-binding_radius)/sqrt(Dtot) - erfc_arg * erfc_arg) * Faddeeva::erfcx( erfc_arg ) / alpha;
							}

								
							// outer scaling.
							probvec1 *= 2*M_PI*particle_density * kact * binding_radius / (kact+ kdiff);
							
		
#endif
	
	//						printf("testing: %le distance: %le\n", probvec1, distance[ap] );
	
	#ifdef DO_HISTOGRAM
							int dbin = dz * N_BINS_PHIST / Rmax;
							if( dbin >= N_BINS_PHIST )
								dbin = N_BINS_PHIST-1;
							sum_prob_distance[dbin] += probvec1;
	//						n_prob_distance[dbin] += 1; 
	#endif						
	
							double rn = rand() / (double)RAND_MAX;


#ifdef NUMERICAL_INTEGRATION
							if( rn < running_prob )
							{
								int lower_p = 0;
								int upper_p = npts;

								if( running_prob > 1 )
									rn *= running_prob;

								int done = 0;

								while(!done)
								{
									int p_face_trial = (lower_p+upper_p)/2;

									if( rn < probs[p_face_trial] )
										upper_p = p_face_trial;
									else
										lower_p = p_face_trial;


									if( upper_p == lower_p + 1 )
										done = 1;
								}
								if( lower_p >= npts )
									lower_p = npts-1;

								int sf = f_list[lower_p];
								double xpuv[2] = { 
									(puv_list[6*lower_p+0]+puv_list[6*lower_p+2]+puv_list[6*lower_p+4])/3,
									(puv_list[6*lower_p+1]+puv_list[6*lower_p+3]+puv_list[6*lower_p+5])/3 };
								double su = xpuv[0];							
								double sv = xpuv[1];							
#else
							if( rn < probvec1 )
							{

#endif
								// move the aqueous particle into ``complex'' position.
								
					
								if( sf != pfaces[ap] )
								{	
									sub_surface->removeParticleFromFace( pfaces[ap], ap, c0_p[ap], particle_footprint );
									sub_surface->addParticleToFace( sf, ap, c0_p[ap], particle_footprint );
								}
								pfaces[ap] = sf;
								puv[2*ap+0] = su;
								puv[2*ap+1] = sv;
	
								sub_surface->evaluateRNRM( sf, su, sv, surface_p_r_m+3*ap, surface_p_r_n+3*ap, r );
	
								aqueous_r[3*ap+0] = surface_p_r_m[3*ap+0] + surface_p_r_n[3*ap+0] * sp_ap_displacement;
								aqueous_r[3*ap+1] = surface_p_r_m[3*ap+1] + surface_p_r_n[3*ap+1] * sp_ap_displacement;
								aqueous_r[3*ap+2] = surface_p_r_m[3*ap+2] + surface_p_r_n[3*ap+2] * sp_ap_displacement;
	
								mean_field_activated[ap] = 1;
								aqueous_status[ap] = STATE_IN_COMPLEX + REACTED_BIT;
	
								tot_on += 1;
							}
							else
								aqueous_status[ap] = STATE_REACTING;
							
#ifdef NUMERICAL_INTEGRATION
								free(probs);
							}
							free(f_list);
							free(puv_list);
							free(areas);
#endif
						}	
					}
				}
				else
				{
					for( tracked_pair *apair = tracked_pairs; apair; apair = apair->next )
					{
						int sp = apair->surface_id;
						int ap = apair->aqueous_id;
		
						if( apair->state== STATE_REACTING || apair->state == STATE_FREE )
						{
							if( (surface_status[sp]&STATE_MASK) == STATE_IN_COMPLEX ||
							    (aqueous_status[ap]&STATE_MASK) == STATE_IN_COMPLEX )		
								continue;
		
							double curnorm = 1.0;
							double p0_ratio = 1.0;
		
							double probvec1 = passocF(apair->sep, time_step, Dtot, binding_radius, alpha, kact/(kact+kdiff));
							double currnorm = 1.0;
				
							if( apair->state == STATE_REACTING ) 
							{	// this pair were previously in the reacting zone.
			
								double p0_ratio = pirr_pfree_ratio_psF(apair->sep, apair->prev_sep, time_step, Dtot, binding_radius, alpha, apair->ps_prev, rtol);
					
								currnorm = apair->prev_norm * p0_ratio;
	//							printf("Modifying ratio! %le * %le = %le\n", apair->prev_norm, p0_ratio, currnorm );
		
							}
							
	#ifdef DO_HISTOGRAM
							int sf;
							double su,sv;
							sub_surface->nearPointOnBoxedSurface( aqueous_r+3*ap, &sf, &su, &sv, M, mlow, mhigh, distance+ap, -1 );
							int dbin = distance[ap] * N_BINS_PHIST / Rmax;
							if( dbin >= N_BINS_PHIST )
								dbin = N_BINS_PHIST-1;
							sum_prob_distance[dbin] += currnorm * probvec1;
	#endif						
							double p = currnorm * probvec1;
					
							if( rand() / (double)RAND_MAX < p )
							{
		//						printf("Reacted! p: %le \n", p);
								apair->state = STATE_IN_COMPLEX;
			
								// move the aqueous particle into ``complex'' position.
			
								aqueous_r[3*ap+0] = surface_p_r_m[3*sp+0] + surface_p_r_n[3*sp+0] * sp_ap_displacement;
								aqueous_r[3*ap+1] = surface_p_r_m[3*sp+1] + surface_p_r_n[3*sp+1] * sp_ap_displacement;
								aqueous_r[3*ap+2] = surface_p_r_m[3*sp+2] + surface_p_r_n[3*sp+2] * sp_ap_displacement;
		
								tot_on += 1;
								surface_status[sp] = STATE_IN_COMPLEX + REACTED_BIT;
								aqueous_status[ap] = STATE_IN_COMPLEX + REACTED_BIT;
							}
							else
							{
								apair->prev_norm = currnorm;
								apair->ps_prev = 1 - probvec1 * currnorm;
								apair->state = STATE_REACTING;
								
								surface_status[sp] = STATE_REACTING;
								aqueous_status[ap] = STATE_REACTING;
							}
						}
						else if( apair->state == STATE_IN_COMPLEX )
						{
							double p = 1 - exp( -kr_off * time_step );
						
							if( (double)rand() / (double)RAND_MAX < p )
							{
								apair->state = STATE_REACTING;
		//						printf("Dissociated! p: %le \n", p );
								surface_status[sp] = STATE_REACTING + REACTED_BIT;	
								aqueous_status[ap] = STATE_REACTING + REACTED_BIT;	
								tot_off += 1;
							}		
						}
					}
				}
	
	#ifdef DO_HISTOGRAM
				for( int p = 0; p < aqueous_np; p++ )
				{
					if( distance[p] < Rmax )
					{
						int dbin = distance[p] * N_BINS_PHIST / Rmax;
						if( dbin >= N_BINS_PHIST )
							dbin = N_BINS_PHIST-1;
						n_prob_distance[dbin] += 1;	
					}
				}
	#endif
				for( int x = 0; x < nv; x++ )
				{
					double dr[3] = { r[3*x+0] - ro[3*x+0], 
						         r[3*x+1] - ro[3*x+1],
						         r[3*x+2] - ro[3*x+2] };
	
					membrane_fudge += dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
					nmembrane_fudge += 1;
				}
	
				if( move_membrane )
				{
					if( do_grad )
					{
						memcpy( mode_t_p, mode_t, sizeof(double) * n_bd_modes );
						memset( g, 0, sizeof(double) * nc );
						memset(mgrad, 0, sizeof(double) * n_bd_modes );
						memset(pgrad, 0, sizeof(double) * (2 * surface_np+1) );
		
						// pgrad has the derivative of the particle energy with respect to its movement on the membrane.
						sub_surface->grad( r, g, puv, pgrad );
						double ec = sub_surface->energy( r, puv, -1, NULL, NULL, 0 );
						// particle particle interactions.
						
						// dV/dR
		
						// brute force to debug.
						
						memset( pcgrad, 0, sizeof(double ) * 3 * surface_np );
		
						double ep = 0;
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
							if( surface_np == 0 )
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
				}
	
				double *vertex_data = NULL;
				int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
				int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
				sub_surface->buildFaceData( &vertex_data, ptr_to_data, nump );
				
				int okay_surface_flag[surface_np];
				int okay_aqueous_flag[aqueous_np];
				memset( okay_surface_flag, 0, sizeof(int) * surface_np );
				memset( okay_aqueous_flag, 0, sizeof(int) * aqueous_np );
	
	
				memcpy( pfaces_cache, pfaces, sizeof(int) * surface_np );
				memcpy( puv_cache, puv, sizeof(double) * 2 * surface_np );
				memcpy( aqueous_r_cache, aqueous_r, sizeof(double) * 3 * aqueous_np );
	
				int done = 0;
	
				// loop moves.
	
				int iters = 0;
				while ( ! done ) 
				{
					if( iters > 100 )
					{
						printf("Terrible collision problem.\n");
						exit(1);
					}
	
					done = 1;
	
					for( int p = 0; p < surface_np; p++ )
					{
						if( okay_surface_flag[p] ) continue;
		
						if( surface_status[p] & REACTED_BIT ) continue;
						if( block.mean_field && !mean_field_activated[p] )
							continue; 
		
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
						
						double old_p[3] = { surface_p_r_m[3*p+0], surface_p_r_m[3*p+1], surface_p_r_m[3*p+2] };
		
						updateParticleR( p, pfaces, puv, surface_p_r_m, r, sub_surface, surface_np ); 
		
						double pdr[3] = { surface_p_r_m[3*p+0] - old_p[0], surface_p_r_m[3*p+1] - old_p[1], surface_p_r_m[3*p+2] - old_p[2] };
						double dr_squared = pdr[0]*pdr[0]+pdr[1]*pdr[1]+pdr[2]*pdr[2];
		
						if( f != of )
						{
							sub_surface->removeParticleFromFace( of, p, c0_p[p], particle_footprint );
							sub_surface->addParticleToFace( pfaces[p], p, c0_p[p], particle_footprint );
						}	
		
						// am I carrying an aqueous particle with me?		
						
						if( block.mean_field || (surface_status[p]&STATE_MASK) == STATE_IN_COMPLEX )
						{
							int ap = -1;
		
							if( block.mean_field )
								ap = p;
							else
							{ 
								for( tracked_pair *apair = tracked_pairs; apair; apair = apair->next )
								{
									if( apair->surface_id == p && apair->state == STATE_IN_COMPLEX )
										ap = apair->aqueous_id;
								}
							}
		
							if( ap == -1 )
							{
								printf("MAJOR PROBLEM here.\n");
								exit(1);
							}	
							aqueous_r[3*ap+0] = surface_p_r_m[3*p+0] + surface_p_r_n[3*p+0] * sp_ap_displacement;	
							aqueous_r[3*ap+1] = surface_p_r_m[3*p+1] + surface_p_r_n[3*p+1] * sp_ap_displacement;	
							aqueous_r[3*ap+2] = surface_p_r_m[3*p+2] + surface_p_r_n[3*p+2] * sp_ap_displacement;	
								
							sub_surface->updateParticle( aqueous_r+3*ap, surface_np+ap, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
						}
						
						okay_surface_flag[p] = 1;
					}
		
					for( int p = 0; p < aqueous_np; p++ )
					{
						if( okay_aqueous_flag[p] ) continue;
		
						if( aqueous_status[p] & REACTED_BIT || aqueous_status[p] == STATE_IN_COMPLEX ) { continue; }
		
						int done = 0;
		
						while( !done )
						{
							done = 1;
		
							double dx = gsl_ran_gaussian(rng_x,  sqrt(2 * block.aqueous_diffc * time_step) );
							double dy = gsl_ran_gaussian(rng_x,  sqrt(2 * block.aqueous_diffc * time_step) );
							double dz = gsl_ran_gaussian(rng_x,  sqrt(2 * block.aqueous_diffc * time_step) );
		
							aqueous_r[3*p+0] += dx;
							aqueous_r[3*p+1] += dy;
							aqueous_r[3*p+2] += dz;
		
							double dmove[3] = { 
								aqueous_r[3*p+0] - aqueous_r_last[3*p+0],
								aqueous_r[3*p+1] - aqueous_r_last[3*p+1],
								aqueous_r[3*p+2] - aqueous_r_last[3*p+2] };
							double r_move = sqrt(dmove[0]*dmove[0]+dmove[1]*dmove[1]+dmove[2]*dmove[2]);
							if( r_move > distance[p] )
							{
								int f;
								double u,v;
								double dist;
			
								double rad = sub_surface->returnRadius( aqueous_r+3*p, &f, &u, &v,  M, mlow, mhigh,  r_move, distance[p],  -1,  defaultR, vertex_data, ptr_to_data, nump, 1 );
							
								if( rad < 0 )
								{
								//	printf("Outside, redoing.\n");
									aqueous_r[3*p+0] -= dx;
									aqueous_r[3*p+1] -= dy;
									aqueous_r[3*p+2] -= dz;
									done = 0;
								}
								else
								{
							//		printf("rad: %lf pt: %le %le %le\n", rad, aqueous_r[3*p+0], aqueous_r[3*p+1], aqueous_r[3*p+2] );
									aqueous_r_last[3*p+0] = aqueous_r[3*p+0];
									aqueous_r_last[3*p+1] = aqueous_r[3*p+1];
									aqueous_r_last[3*p+2] = aqueous_r[3*p+2];
									distance[p] = rad;
								}
							}
		
							if( done )
								sub_surface->updateParticle( aqueous_r+3*p, surface_np+p, r[3*nv+0], r[3*nv+1], r[3*nv+2] );
						}
		
						while( aqueous_r[3*p+0] < -LA/2 ) { aqueous_r[3*p+0] += LA; aqueous_r_last[3*p+0] += LA; }
						while( aqueous_r[3*p+0] > LA/2 )  { aqueous_r[3*p+0] -= LA; aqueous_r_last[3*p+0] -= LA; }
						while( aqueous_r[3*p+1] < -LB/2 ) { aqueous_r[3*p+1] += LB; aqueous_r_last[3*p+1] += LB; }
						while( aqueous_r[3*p+1] > LB/2 )  { aqueous_r[3*p+1] -= LB; aqueous_r_last[3*p+1] -= LB; }
						while( aqueous_r[3*p+2] < -LC/2 ) { aqueous_r[3*p+2] += LC; aqueous_r_last[3*p+2] += LC; }
						while( aqueous_r[3*p+2] > LC/2 )  { aqueous_r[3*p+2] -= LC; aqueous_r_last[3*p+2] -= LC; }
		
						okay_aqueous_flag[p] = 1;
					}
	
					// check for collisions.
					//
	
					if( !block.mean_field )
					{
						int npairs_collide = get_pair_list( sub_surface, surface_p_r_m, aqueous_r, &pair_list, &npairsSpace, surface_np, aqueous_np, binding_radius );
						int ncol = 0;
							
						for( int x = 0; x < npairs_collide; x++ )
						{
							int sp = pair_list[2*x+0];
							int ap = pair_list[2*x+1];
							if( (surface_status[sp]&STATE_MASK) == STATE_IN_COMPLEX )
							{
								int apx = -1;
		
								for( tracked_pair *apair = tracked_pairs; apair; apair = apair->next )
								{
									if( apair->surface_id == sp && apair->state == STATE_IN_COMPLEX )
										apx = apair->aqueous_id;
								}
	
								if( apx == ap ) continue;
							}	
		
							if( okay_surface_flag[pair_list[2*x+0]] )
							{ 
								okay_surface_flag[pair_list[2*x+0]] = 0; 
						
								if( pfaces[pair_list[2*x+0]] != pfaces_cache[pair_list[2*x+0]] )
								{
									sub_surface->removeParticleFromFace( pfaces[pair_list[2*x+0]], sp, c0_p[sp], particle_footprint );
									sub_surface->addParticleToFace( pfaces_cache[pair_list[2*x+0]], sp, c0_p[sp], particle_footprint );
								}	
								pfaces[pair_list[2*x+0]] = pfaces_cache[pair_list[2*x+0]];
								puv[pair_list[2*x+0]*2+0] = puv_cache[pair_list[2*x+0]*2+0];
								puv[pair_list[2*x+0]*2+1] = puv_cache[pair_list[2*x+0]*2+1];
								ncol++;
								done = 0; 
							}
							if( okay_aqueous_flag[pair_list[2*x+1]] )
							{
								okay_aqueous_flag[pair_list[2*x+1]] = 0;
								aqueous_r[pair_list[2*x+1]*3+0] = aqueous_r_cache[pair_list[2*x+1]*3+0];
								aqueous_r[pair_list[2*x+1]*3+1] = aqueous_r_cache[pair_list[2*x+1]*3+1];
								aqueous_r[pair_list[2*x+1]*3+2] = aqueous_r_cache[pair_list[2*x+1]*3+2];
								sub_surface->updateParticle( aqueous_r+3*pair_list[2*x+1], surface_np+pair_list[2*x+1], r[3*nv+0], r[3*nv+1], r[3*nv+2] );
								done=0;
							}
						}
					}
	
					iters++;
				}
	
				if( block.mean_field )
				{
				}
				else
				{
					for( tracked_pair *apair = tracked_pairs; apair; apair = apair->next )
						apair->flag = 0;
	
					int npairs = get_pair_list( sub_surface, surface_p_r_m, aqueous_r, &pair_list, &npairsSpace, surface_np, aqueous_np, Rmax );
			
					for( int px = 0; px < npairs; px++ )
					{
						tracked_pair *gotit = NULL;
						for( tracked_pair *apair = tracked_pairs; apair; apair = apair->next )
						{
							if( apair->surface_id == pair_list[2*px+0] && apair->aqueous_id == pair_list[2*px+1] )
							{
								gotit = apair;
								apair->flag = 1;
								break;
							}
						}
		
						if( !gotit )
						{		
							tracked_pair *apair = (tracked_pair *)malloc( sizeof(tracked_pair) );
							apair->surface_id = pair_list[2*px+0];
							apair->aqueous_id = pair_list[2*px+1];
							apair->prev_norm = 1.0;
							apair->ps_prev   = 1.0;
							apair->weight = 1.0;
							apair->state = STATE_FREE;
							apair->flag = 2; // new
							apair->sep = -1;
							apair->prev_sep = -1;
							apair->next = tracked_pairs;
							tracked_pairs = apair;
		
							gotit = apair;
						}
					
						int sp = pair_list[2*px+0];
						int ap = pair_list[2*px+1];
		
						double dr[3] = {
						surface_p_r_m[sp*3+0] - aqueous_r[3*ap+0],
						surface_p_r_m[sp*3+1] - aqueous_r[3*ap+1],
						surface_p_r_m[sp*3+2] - aqueous_r[3*ap+2] };
		
						while( dr[0] < -LA/2 ) dr[0] += LA;
						while( dr[0] >  LA/2 ) dr[0] -= LA;
						while( dr[1] < -LB/2 ) dr[1] += LB;
						while( dr[1] >  LB/2 ) dr[1] -= LB;
						while( dr[2] < -LC/2 ) dr[2] += LC;
						while( dr[2] >  LC/2 ) dr[2] -= LC;
		
						double sep = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
		
						gotit->prev_sep = gotit->sep;
						gotit->sep = sep;
					}
	#define DEBUG_PAIRS
					int ncnt = 0;
	#ifdef DEBUG_PAIRS
					int nfound = 0;
					for( int sp = 0; sp < surface_np; sp++ )
					{
						for( int ap = 0; ap < aqueous_np; ap++ )
						{
							double dr[3] = { 
								aqueous_r[3*ap+0] - surface_p_r_m[3*sp+0],
								aqueous_r[3*ap+1] - surface_p_r_m[3*sp+1],
								aqueous_r[3*ap+2] - surface_p_r_m[3*sp+2] };
		
							while( dr[0] < -LA/2 ) dr[0] += LA;
							while( dr[1] < -LB/2 ) dr[1] += LB;
							while( dr[2] < -LC/2 ) dr[2] += LC;
							
							while( dr[0] >  LA/2 ) dr[0] -= LA;
							while( dr[1] >  LB/2 ) dr[1] -= LB;
							while( dr[2] >  LC/2 ) dr[2] -= LC;
		
							double rl = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	
							// 1070 2965	
							if( sp == 1070 && ap == 2965 )
							{
	//							printf("D 1070-2965: %le\n", rl );
							}
		
							if( rl < Rmax )
							{
								nfound++;
								int pf = 0, pf_noflag = 0;
								for( tracked_pair *apair = tracked_pairs; apair; apair = apair->next )
								{
									if( apair->surface_id == sp && apair->aqueous_id == ap  )
									{
										pf_noflag = 1;
	//									if( sp == 1070 && ap == 2965 )
	//									{
	//									printf("1070/2965 apair: %p flag: %d\n", apair, apair->flag );
	//									}
									}
									if( apair->surface_id == sp && apair->aqueous_id == ap && apair->flag )
									{
										pf = 1;
									}
	
								}
			
								if( !pf )	
								{
									if( pf_noflag ) printf("NOFLAG.\n");
									printf("Pair not found rl: %lf/%lf sp: %d ap: %d.\n", rl, Rmax, sp, ap  );
									exit(1);
								}		
							}	
						}
					}
	#endif
		
					tracked_pair *prev = NULL;
					tracked_pair *next = NULL;
					for( tracked_pair *apair = tracked_pairs; apair; apair = next)
					{
						next = apair->next;
						if( apair->flag == 0 )
						{
							if( apair->state == STATE_IN_COMPLEX )
							{
								int sp = apair->surface_id;
								int ap = apair->aqueous_id;
								double dr[3] = { 
									aqueous_r[3*ap+0] - surface_p_r_m[3*sp+0],
									aqueous_r[3*ap+1] - surface_p_r_m[3*sp+1],
									aqueous_r[3*ap+2] - surface_p_r_m[3*sp+2] };
								double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
								printf("DESTROYING COMPLEX STATE r: %le %d %d .\n", r, sp, ap);
								exit(1);
							}
							if( prev )
								prev->next = apair->next;
							else
								tracked_pairs = apair->next;
		
							free(apair);
						}
						else
							prev = apair;	
					}
					
					for( tracked_pair *apair = tracked_pairs; apair; apair = apair->next )
						ncnt++;	
	
				}	
	
	
				nav++;
	
				int npaired = 0;
			
				if( block.mean_field )
				{
					for( int ap = 0; ap < aqueous_np; ap++ )
					{
						if( mean_field_activated[ap] )
							npaired++;	
					}
		
				}
				else
				{
					for( int s = 0; s < surface_np; s++ )
					{
						if( (surface_status[s]&STATE_MASK) == STATE_IN_COMPLEX )
							npaired++;
					}
				}

				if( ncind < ncomplex_size )	
				{
					ncomplex[ncind] += npaired;
					ncind++;
				}
				fprintf(n_complex_file, "%le %d %lf %lf\n", cur_t*(1e9), npaired, tot_on, tot_off );
//				printf("ncind: %d ncompl %le\n", ncind, npaired );
				if( ncind >= ncomplex_size && nruns > 1 )
				{
					printf("RUN DONE stopping.\n");
					run_done = 1;
				}
				free(vertex_data);
				free(ptr_to_data);
				free(nump);
			}
			fflush(n_complex_file);
	
			int npaired=0;
	
			if( block.mean_field )
			{
				for( int ap = 0; ap < aqueous_np; ap++ )
				{
					if( mean_field_activated[ap] )
						npaired++;	
				}
				printf("%d/%d aqueous particles are in complex [Mean Field].\n", npaired, aqueous_np );
	
			}
			else
			{
				for( int s = 0; s < surface_np; s++ )
				{
					if( (surface_status[s]&STATE_MASK) == STATE_IN_COMPLEX )
						npaired++;
				}
				printf("%d/%d surface particles are in complex.\n", npaired, surface_np );
			}		
	#ifdef DO_HISTOGRAM
			for( int b = 0; b < N_BINS_PHIST; b++ )	
				fprintf( histo_file, " %le", sum_prob_distance[b]/(1e-15+n_prob_distance[b]) );
			fprintf(histo_file, "\n"); 
			fflush(histo_file);
	#endif						
	
			if( tFile )
			{
				sub_surface->put(r);
	        		sub_surface->writeLimitingSurface(tFile);
				if( !block.mean_field )
					fprintf(pFile, "%d\n", surface_np + aqueous_np );
				else
					fprintf(pFile, "%d\n",  aqueous_np );
				fprintf(pFile, "surface bound points (C) aqueous points (O)\n");
				if( !block.mean_field )
				{
					for( int p = 0; p < surface_np; p++ )
						fprintf(pFile, "C %lf %lf %lf\n", surface_p_r_m[3*p+0], surface_p_r_m[3*p+1], surface_p_r_m[3*p+2] );
				}
				for( int p = 0; p < aqueous_np; p++ )
					fprintf(pFile, "O %lf %lf %lf\n", aqueous_r[3*p+0], aqueous_r[3*p+1], aqueous_r[3*p+2] );
	 
				fflush(tFile);
				fflush(pFile);
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
					run_done = 1;
					break;
				}
				if( o % 50 == 0 && hours > 0 )
					printf("Currently we have run for %lf hours, will stop after %lf.\n", dt, hours );
				
			}
			else
			{
				if( o >= nsteps )
					run_done = 1;
			}
			
	//		printf("membrane fudge: %le\n", membrane_fudge / nmembrane_fudge );
		}

		fclose(n_complex_file);

		sprintf(fileName, "%s_multi_complex.txt", block.jobName );
		FILE *multi_complex_file = fopen(fileName,"w");
		for( int c = 0; c < ncind; c++ )
			fprintf(multi_complex_file, "%le %le\n", c*t_per_complex, ncomplex[c] / (run+1) );
		fflush( multi_complex_file);
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

	for( int p = 0; p < surface_np; p++ )
		fprintf(saveFile, "%d %lf %lf\n", pfaces[p], puv[2*p+0], puv[2*p+1] );
	fclose(saveFile);

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
	
}


int get_pair_list( surface *sub_surface, double *surface_r, double *aqueous_r, int **pairs, int *npairsSpace, int surface_np, int aqueous_np, double cutoff )
{
	int *out_list = *pairs;

	int npairs = 0;

	double La = sub_surface->PBC_vec[0][0];
	double Lb = sub_surface->PBC_vec[1][1];
	double Lc = sub_surface->PBC_vec[2][2];

	for( int p = 0; p < surface_np; p++ )
	{
		int nnear = sub_surface->getNear( surface_r+p*3, p, surface_r, cutoff, 1., 1., 1., global_plist, global_rads);	

		for( int x = 0; x < nnear; x++ )
		{			
			int p2 = global_plist[x];

			if( p2 < surface_np ) continue;

			int p2x = p2 - surface_np;

			double dr[3] = {
				surface_r[p*3+0] - aqueous_r[3*p2x+0],
				surface_r[p*3+1] - aqueous_r[3*p2x+1],
				surface_r[p*3+2] - aqueous_r[3*p2x+2] };

			while( dr[0] < -La/2 ) dr[0] += La;
			while( dr[0] >  La/2 ) dr[0] -= La;
			while( dr[1] < -Lb/2 ) dr[1] += Lb;
			while( dr[1] >  Lb/2 ) dr[1] -= Lb;
			while( dr[2] < -Lc/2 ) dr[2] += Lc;
			while( dr[2] >  Lc/2 ) dr[2] -= Lc;
	
			double rval = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( rval < cutoff)
			{
				if( npairs == *npairsSpace )
				{
					*npairsSpace *= 2;
		
					*pairs = (int *)realloc( *pairs, sizeof(int) * 2 * *npairsSpace );
					out_list = *pairs;
				}

				out_list[2*npairs+0] = p;
				out_list[2*npairs+1] = p2x;
			
				npairs++;
			}
		}
	}

	return npairs;
}
	

