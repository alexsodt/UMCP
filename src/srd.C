/* SRD algorithm */
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
#include "alignSet.h"
#include "srd.h"
#include "units.h"
//#define SPECIFIC_PARTICLE
#define TAG_MECHANISM
#define MEMBRANE_LEVEL_TOLERANCE (0.0)
static double FUDGE_DISTANCE = 1.0;
static double SAFETY_FACTOR  = 2.0;

static gsl_rng * r_gen = NULL;
static double running_temp = 0;
static int nrunning = 0;
//#define DEBUG_MOMENTUM_CHANGE
// void rotateArbitrary( double *xyz_rot, double *axis, double *origin, int nat, double val );

int ipw2( int y )
{
	double val = 1.0;
	for( int ty = 0; ty < y; ty++ )
		val *= 2;
	return val;
}

void srd_integrator::init( double av_edge_length, double PBC_vec_in[3][3], double temp_in /* kcal/mol */, double eta_in /* SI */, double dt_si, double dt_in /* seconds */, int doPlanarTopology, double mass_in, double srd_M, int hard_z_boundary_in)
{
	hard_z_boundary = hard_z_boundary_in;
	debug_mode = 0;
	// the mass of a membrane unit that receives the impulse force.
	membrane_point_mass = mass_in; 

	double binw = av_edge_length;
	double M = srd_M; //(double)np_in / (double)(grain_x * grain_y * grain_z); 
	planar = doPlanarTopology;
	// compute the mass to get the correct eta.
	// collision frequency is 20 timesteps.
	for( int c1 = 0; c1 < 3; c1++ )
	for( int c2 = 0; c2 < 3; c2++ )
		PBC_vec[c1][c2] = PBC_vec_in[c1][c2];
	double sim_vol = PBC_vec[0][0]*PBC_vec[1][1]*PBC_vec[2][2]; 
	grain_x = PBC_vec[0][0] / binw;
	grain_y = PBC_vec[1][1] / binw;
	grain_z = PBC_vec[2][2] / binw;
	int np_in =  lround(M * grain_x*grain_y*grain_z);	
	printf("Solvent particles: %d\n", np_in );
	M = (double)np_in / (double)(grain_x*grain_y*grain_z);
	double tau_col = dt_si;
	
	double temperature_joules = temp_in * 4.184 * 1000 / (6.022e23);
	double alpha = rotation_angle = M_PI/2;
	
	// solve for approximate bin width that gives us the target mass.

	// M = np_in / (grain_cubed)
	// grain_cubed == vol / (a^3)
	// M = np_in * a^3 / vol

	// I guess this depends only on density (for a simple gas)	
	a = binw * (1e-10); // meters

	double cell_vol = (1e-30)*PBC_vec[0][0] *PBC_vec[1][1] * PBC_vec[2][2] / (grain_x * grain_y * grain_z);
	printf("M: %lf\n", M );
	double eta_streaming = temperature_joules * tau_col * M / cell_vol / ( 2 ) * ( (5 * M)/((M-1+exp(-M))*(2-cos(alpha)-cos(2*alpha)))-1);
	double alt_streaming = temperature_joules * tau_col * M / cell_vol * (2*M+3) / ( 6 * M-1);
	double eta_collision_per_mass = (1.0/(18.0*a*tau_col)) * ( M - 1 + exp(-M)) * (1- cos(rotation_angle));
	
	mass= (eta_in - eta_streaming)/eta_collision_per_mass * AMU_PER_KG;
	double eta_collisional = eta_collision_per_mass * mass / AMU_PER_KG;

	printf("With # particles requested:\n");
	printf("eta in: %le eta_streaming: %le mass: %le\n", eta_in, eta_streaming, mass );
	printf("total_viscosity %le\n", eta_streaming + eta_collisional );
	
	//exit(1);


	a *= (1e10); // back into angstroms.
	
	
	
	int nbins = grain_x * grain_y * grain_z ;
	bins = (bin *)malloc( sizeof(bin) * nbins );

	for( int x = 0; x < nbins; x++ )
	{
		bins[x].np = 0;
		bins[x].npSpace = 5;
		bins[x].list = (int *)malloc( sizeof(int) * bins[x].npSpace );
	}

	np = np_in;
	temp = temp_in;
	eta = eta_in;

	if( ! r_gen ) {
		r_gen = gsl_rng_alloc (gsl_rng_taus);
		gsl_rng_set( r_gen, rand() );
	}
	rp = (double *)malloc( sizeof(double) * 3 * np_in );
	rp_prev = (double *)malloc( sizeof(double) * 3 * np_in );
	rp_dist = (double *)malloc( sizeof(double) * 3 * np_in );
	distance = (double *)malloc( sizeof(double) * np_in );
	vp = (double *)malloc( sizeof(double) * 3 * np_in );
	inside_outside = (int *)malloc( sizeof(int) * np_in );
	topological_tag = (double *)malloc( sizeof(double) * np_in );
	alt_topological_tag = (double *)malloc( sizeof(double) * np_in );
	memset( topological_tag, 0, sizeof(double) * np_in );
	memset( alt_topological_tag, 0, sizeof(double) * np_in );
	last_known_tag = (double *)malloc( sizeof(double) * np_in );
	alt_last_known_tag = (double *)malloc( sizeof(double) * np_in );
	memset( last_known_tag, 0, sizeof(double) * np_in );
	memset( alt_last_known_tag, 0, sizeof(double) * np_in );
	did_collide = (int *)malloc( sizeof(int) * np_in );
	for( int x = 0; x < np_in; x++ )
		inside_outside[x] = POINT_UNKNOWN; // UNKNOWN! (for now we will cross this bridge when we come to it).

	double DOF = 3.0;


	double maxv = 0;

	for( int x = 0; x < 100; x++ )
	{
		double ke = gsl_ran_chisq( r_gen, DOF );
		double vel = sqrt(ke);
		// draw random orientation
		double dx,dy,dz;
		gsl_ran_dir_3d( r_gen, &dx,&dy,&dz); 

		if( vel > maxv )
			maxv = vel;
	}

#ifdef DIST	
	int nbins = 100;
	double hist[nbins];
	memset(hist,0,sizeof(double)*nbins);
	
	for( int x = 0; x < 100000; x++ )
	{
		double ke = gsl_ran_chisq( r_gen, DOF );
		double vel = sqrt(ke);
		// draw random orientation
		double dx,dy,dz;
		gsl_ran_dir_3d( r_gen, &dx,&dy,&dz); 

		int vbin = nbins * vel / maxv;

		if( vbin >= nbins ) vbin = nbins-1;

		hist[vbin]+=1;
	}
	
	for( int x = 0; x < nbins; x++ )
	{
		printf("%lf %lf\n", (x+0.5) * maxv / nbins, hist[x] );
	}
#endif

	double NDOF = 3*np;
	double sum_ke = 0;

	for( int x = 0; x < np; x++ )
	{
		rp[3*x+0] = -PBC_vec[0][0]/2 + PBC_vec[0][0] * (double)rand() / (double)RAND_MAX;
		rp[3*x+1] = -PBC_vec[1][1]/2 +PBC_vec[1][1] * (double)rand() / (double)RAND_MAX;
		rp[3*x+2] = -PBC_vec[2][2]/2 +PBC_vec[2][2] * (double)rand() / (double)RAND_MAX;
		

		// drawn from PDF: 1 / sqrt(M_PI/2) * x^2 exp(-x/2) 
		double ke = gsl_ran_chisq( r_gen, DOF );
		double vel = sqrt(ke * temp / mass );
		// draw random orientation
		double dx,dy,dz;
		gsl_ran_dir_3d( r_gen, &dx,&dy,&dz); 

		vp[3*x+0] = vel * dx;
		vp[3*x+1] = vel * dy;
		vp[3*x+2] = vel * dz;

		sum_ke += mass * (vp[3*x+0]*vp[3*x+0] + vp[3*x+1]*vp[3*x+1] + vp[3*x+2]*vp[3*x+2]) / 2;  

		distance[x] = 0;
	}

	memcpy( rp_prev, rp, sizeof(double) * 3 * np );
	memcpy( rp_dist, rp, sizeof(double) * 3 * np );

	double average_velocity = sqrt((temp*3/mass));
	double time_free_path_grain = a / average_velocity;		
	printf("SRD temperature: %le\n", sum_ke / NDOF / 0.592 );
	srdCollisionEstimator = average_velocity * dt_si * AKMA_TIME * 2 * SAFETY_FACTOR;
	if( srdCollisionEstimator > a )
		srdCollisionEstimator = a;
	printf("srdCollisionEstimator: %le\n", srdCollisionEstimator );

	{	
		double vol = PBC_vec[0][0] * PBC_vec[1][1] * PBC_vec[2][2];
		double rho = np / vol;
	
		// lowest level: approx 5 per.
		double low_vol = 5 / rho;
		double approx_len = pow( low_vol, 1.0 / 3.0 );
	
		grain_x_p = 1 + log( PBC_vec[0][0] / approx_len ) / log(2);
		grain_y_p = 1 + log( PBC_vec[1][1] / approx_len ) / log(2);
		grain_z_p = 1 + log( PBC_vec[2][2] / approx_len ) / log(2);
		
		id_grain_x = (int) pow( 2, grain_x_p );
		id_grain_y = (int) pow( 2, grain_x_p );
		id_grain_z = (int) pow( 2, grain_x_p );
	
		// find lowest power, this will limit our log search
		int min_p = grain_x_p;
		if( grain_y_p < min_p ) min_p = grain_y_p;
		if( grain_z_p < min_p ) min_p = grain_z_p;
	
		int nbins_id = id_grain_x * id_grain_y * id_grain_z;
	
		bins_id = (bin *)malloc( sizeof(bin) * nbins_id );
		for( int x = 0; x < nbins_id; x++ )
		{
			bins_id[x].np = 0;
			bins_id[x].npSpace = 5;
			bins_id[x].list = (int *)malloc( sizeof(int) * bins_id[x].npSpace );
		}
		bin_level = (int *)malloc( sizeof(int) * nbins_id );

	}
}

void srd_integrator::stream( double dt )
{
	memcpy( rp_prev, rp, sizeof(double) * 3 * np );

	for( int xp = 0; xp < np; xp++ )
	{
		rp[3*xp+0] += vp[3*xp+0] * dt;
		rp[3*xp+1] += vp[3*xp+1] * dt;
		rp[3*xp+2] += vp[3*xp+2] * dt;

#ifdef SPECIFIC_PARTICLE
//		if( xp == 2267 )
//		printf("2267 advancing from %lf %lf %lf to %lf %lf %lf vp: %le %le %le\n", 
//				rp_prev[3*xp+0], rp_prev[3*xp+1], rp_prev[3*xp+2],
//				rp[3*xp+0], rp[3*xp+1], rp[3*xp+2], vp[3*xp+0], vp[3*xp+1], vp[3*xp+2] );
#endif
		while( rp[3*xp+0] < -PBC_vec[0][0]/2 ) { rp[3*xp+0] += PBC_vec[0][0]; rp_prev[3*xp+0] += PBC_vec[0][0]; }
		while( rp[3*xp+1] < -PBC_vec[1][1]/2 ) { rp[3*xp+1] += PBC_vec[1][1]; rp_prev[3*xp+1] += PBC_vec[1][1]; }
		while( rp[3*xp+0] > PBC_vec[0][0]/2 ) { rp[3*xp+0] -= PBC_vec[0][0]; rp_prev[3*xp+0] -= PBC_vec[0][0]; }
		while( rp[3*xp+1] > PBC_vec[1][1]/2 ) { rp[3*xp+1] -= PBC_vec[1][1]; rp_prev[3*xp+1] -= PBC_vec[1][1]; }

		if( hard_z_boundary )
		{
			double pdel = rp[3*xp+2] - rp_prev[3*xp+2];
			if( rp[3*xp+2] < -PBC_vec[2][2]/2 )
			{
				double del =  -PBC_vec[2][2]/2 - rp[3*xp+2];

				rp[3*xp+2] = -PBC_vec[2][2]/2 + del;

				vp[3*xp+2] *= -1;
			}
			else if( rp[3*xp+2] > PBC_vec[2][2]/2 )
			{
				double del =  PBC_vec[2][2]/2 - rp[3*xp+2];

				rp[3*xp+2] = PBC_vec[2][2]/2 + del;

				vp[3*xp+2] *= -1;
			}

			rp_prev[3*xp+2] = rp[3*xp+2] - pdel;
		}
		else
		{
			while( rp[3*xp+2] < -PBC_vec[2][2]/2 ) { rp[3*xp+2] += PBC_vec[2][2]; rp_prev[3*xp+2] += PBC_vec[2][2]; }
			while( rp[3*xp+2] > PBC_vec[2][2]/2 ) { rp[3*xp+2] -= PBC_vec[2][2]; rp_prev[3*xp+2] -= PBC_vec[2][2]; }
		}
	}	
}

void srd_integrator::collide( void )
{

	double KE_init = KE();

	double offset_dx = PBC_vec[0][0] * ( rand() / (double)RAND_MAX -0.5) / grain_x; 
	double offset_dy = PBC_vec[1][1] * ( rand() / (double)RAND_MAX -0.5) / grain_y; 
	double offset_dz = PBC_vec[2][2] * ( rand() / (double)RAND_MAX -0.5) / grain_z; 

	int nbins = grain_x * grain_y * grain_z;

	for( int x = 0; x < nbins;x ++ )
		bins[x].np = 0;

	for( int p = 0; p < np; p++ )
	{
		double fx = rp[3*p+0] / PBC_vec[0][0];
		while( fx < 0 ) fx += 1.0;
		while( fx >= 1.0 ) fx -= 1.0;
		int bx = fx * grain_x;
		if( bx == grain_x ) bx--;
		
		double fy = rp[3*p+1] / PBC_vec[1][1];
		while( fy < 0 ) fy += 1.0;
		while( fy >= 1.0 ) fy -= 1.0;
		int by = fy * grain_y;
		if( by == grain_y ) by--;
		
		double fz = rp[3*p+2] / PBC_vec[2][2];
		while( fz < 0 ) fz += 1.0;
		while( fz >= 1.0 ) fz -= 1.0;
		int bz = fz * grain_z;
		if( bz == grain_z ) bz--;

		int bin = (bx * grain_y + by ) * grain_z + bz;

		if( bins[bin].npSpace == bins[bin].np )
		{
			bins[bin].npSpace *= 2;
			bins[bin].list = (int *)realloc( bins[bin].list, sizeof(int) * bins[bin].npSpace );
		}

		bins[bin].list[bins[bin].np] = p;
		bins[bin].np++;
	}

	double zero[3] = {0,0,0};

	double eps = 0.1;
		

	double inst_temp = 0;


	for( int b = 0; b < nbins; b++ )
	{
		double S = (1+eps);

		if( rand() % 2 == 0 )
			S = 1.0 / (1+eps);

		if( bins[b].np < 1 ) continue;

		double vav[3] = {0,0,0};


		for( int px = 0; px < bins[b].np; px++ )
		{
			int p = bins[b].list[px];

			vav[0] += vp[3*p+0]; 
			vav[1] += vp[3*p+1]; 
			vav[2] += vp[3*p+2]; 
		}

		vav[0] /= bins[b].np;
		vav[1] /= bins[b].np;
		vav[2] /= bins[b].np;

		// draw random orientation
		double axis[3];
		gsl_ran_dir_3d( r_gen, axis, axis+1, axis+2 ); 

		double RELKE = 0;

		for( int px = 0; px < bins[b].np; px++ )
		{
			int p = bins[b].list[px];
			
			double dvel[3] = { 
				vp[3*p+0] - vav[0],
				vp[3*p+1] - vav[1],
				vp[3*p+2] - vav[2] };

			RELKE += dvel[0]*dvel[0]+dvel[1]*dvel[1]+dvel[2]*dvel[2];
			
			rotateArbitrary( dvel, axis, zero, 1, rotation_angle );

			vp[3*p+0] = vav[0] + dvel[0];
			vp[3*p+1] = vav[1] + dvel[1]; 
			vp[3*p+2] = vav[2] + dvel[2]; 
				
			inst_temp += 0.5 * mass * (vp[3*p+0] * vp[3*p+0] + vp[3*p+1] * vp[3*p+1] + vp[3*p+2] * vp[3*p+2]);	
		}
			
		if( !debug_mode && bins[b].np > 1 ) 
		{
			double A = pow( S, 3 * (bins[b].np - 1) ) * exp( - mass * ( S*S-1) / (2 * temp) * RELKE );
	
			double rn = rand()/(double)RAND_MAX;
	
			if( rn < A )
			{
				for( int px = 0; px < bins[b].np; px++ )
				{
					int p = bins[b].list[px];
					
					double dvel[3] = { 
						vp[3*p+0] - vav[0],
						vp[3*p+1] - vav[1],
						vp[3*p+2] - vav[2] };
		
					vp[3*p+0] = vav[0] + dvel[0]*S;
					vp[3*p+1] = vav[1] + dvel[1]*S; 
					vp[3*p+2] = vav[2] + dvel[2]*S; 
				}
			}
		}

	}

		
	inst_temp /= (3 * (np));

	running_temp += inst_temp;
	nrunning += 1;

//	if( nrunning % 100 == 0 )		
		printf("temp: %le av %le target: %le\n", inst_temp / 0.592, running_temp / (nrunning) / 0.592, temp );
	
	double KE_out = KE();

	dKE_out_thermostat += KE_out - KE_init;

//	printf("dKE_CHECK SRD dKE: %le dMEM: 0\n", KE_out-KE_init );
}

double srd_integrator::KE( void )
{
	double l_KE = 0;
	for( int p = 0; p < np; p++ )
	{
		l_KE += 0.5 * vp[3*p+0] * vp[3*p+0] * mass;
		l_KE += 0.5 * vp[3*p+1] * vp[3*p+1] * mass;
		l_KE += 0.5 * vp[3*p+2] * vp[3*p+2] * mass;
	}

	return l_KE;
}
/*
void srd_integrator::initializeDistances_works( double *r, surface *theSurface, double **M, int mlow, int mhigh )
{
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );

#if 1 
	int nbins = grain_x * grain_y * grain_z;

	for( int x = 0; x < nbins;x ++ )
		bins[x].np = 0;

	for( int p = 0; p < np; p++ )
	{
		double fx = rp[3*p+0] / PBC_vec[0][0];
		while( fx < 0 ) fx += 1.0;
		while( fx >= 1.0 ) fx -= 1.0;
		int bx = fx * grain_x;
		if( bx == grain_x ) bx--;
		
		double fy = rp[3*p+1] / PBC_vec[1][1];
		while( fy < 0 ) fy += 1.0;
		while( fy >= 1.0 ) fy -= 1.0;
		int by = fy * grain_y;
		if( by == grain_y ) by--;
		
		double fz = rp[3*p+2] / PBC_vec[2][2];
		while( fz < 0 ) fz += 1.0;
		while( fz >= 1.0 ) fz -= 1.0;
		int bz = fz * grain_z;
		if( bz == grain_z ) bz--;

		int bin = (bx * grain_y + by ) * grain_z + bz;

		if( bins[bin].npSpace == bins[bin].np )
		{
			bins[bin].npSpace *= 2;
			bins[bin].list = (int *)realloc( bins[bin].list, sizeof(int) * bins[bin].npSpace );
		}

		bins[bin].list[bins[bin].np] = p;
		bins[bin].np++;
	}

	for( int bx = 0; bx < grain_x; bx++ )
	for( int by = 0; by < grain_y; by++ )
	for( int bz = 0; bz < grain_z; bz++ )
	{
		int b = (bx * grain_y + by ) * grain_z + bz;
		double cen[3] = {  
			(bx+0.5) * PBC_vec[0][0] / grain_x,
			(by+0.5) * PBC_vec[1][1] / grain_y,
			(bz+0.5) * PBC_vec[2][2] / grain_z };

		if( bins[b].np == 0 ) continue;

		double r_max = 0;
		double dr[bins[b].np];

		for( int px = 0; px < bins[b].np; px++ )
		{
			int p = bins[b].list[px];

			double dr[3] = { rp[3*p+0] - cen[0],
					 rp[3*p+1] - cen[1],
					 rp[3*p+2] - cen[2] };

			while( dr[0] < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0];
			while( dr[1] < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1];
			while( dr[2] < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2];
			while( dr[0] >  PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0];
			while( dr[1] >  PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1];
			while( dr[2] >  PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];

			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r > r_max )
				r_max = r;

			dr[px] = r;
					
		}
		
		int f;
		double u,v;
		printf("Checking %d particles near %lf %lf %lf ", bins[b].np, cen[0], cen[1], cen[2] );
		if( ! theSurface->withinRadius( cen, &f, &u, &v,  M, mlow, mhigh, -1,  r_max + srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
		{
			printf("... all are clear.\n");
			for( int px = 0; px < bins[b].np; px++ )
			{
				int p = bins[b].list[px];
				distance[p] = (r_max - dr[px]) + srdCollisionEstimator;

				rp_dist[3*p+0] = rp[3*p+0];
				rp_dist[3*p+1] = rp[3*p+1];
				rp_dist[3*p+2] = rp[3*p+2];				
			}
		}	
		else
		{
			printf("...checking them all...");
			double nok = 0, nbad=0;
			for( int px = 0; px < bins[b].np; px++ )
			{
				int p = bins[b].list[px];
			
				if( ! theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
				{
					nok += 1;
					distance[p] = srdCollisionEstimator;

					rp_dist[3*p+0] = rp[3*p+0];
					rp_dist[3*p+1] = rp[3*p+1];
					rp_dist[3*p+2] = rp[3*p+2];				
				}
				else
				{
					nbad += 1;
					theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, -1 );
					rp_dist[3*p+0] = rp[3*p+0];
					rp_dist[3*p+1] = rp[3*p+1];
					rp_dist[3*p+2] = rp[3*p+2];

					double rpt[3], rnrm[3];

					theSurface->evaluateRNRM( f, u, v, rpt, rnrm, r );

					double dr[3] = { rp[3*p+0] - rpt[0],
						 rp[3*p+1] - rpt[1],
						 rp[3*p+2] - rpt[2] };
	
					double len = normalize(dr);
				
					double dp = dr[0] * rnrm[0] + dr[1] * rnrm[1] + dr[2] * rnrm[2];

//			printf("init len: %le dp: %le dist: %le\n", len, dp, distance[p] );

					if( !planar )
					{
						if( dp > 0 )
							inside_outside[p] = POINT_INSIDE;
						else
							inside_outside[p] = POINT_OUTSIDE;
					}
				}
				printf(" %lf%% were clear.\n", 100 * nok / (nok+nbad) );
			}
		}
	}

#else
	for( int p = 0; p < np; p++ )
	{	
		int f;
		double u,v;	
		if( ! theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
		{
#ifdef CHECK_ID
			theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, srdCollisionEstimator );
			
			double rpt[3], rnrm[3];

			theSurface->evaluateRNRM( f, u, v, rpt, rnrm, r );

			double dr[3] = { rp[3*p+0] - rpt[0],
					 rp[3*p+1] - rpt[1],
					 rp[3*p+2] - rpt[2] };

			double len = normalize(dr);

			if( len < srdCollisionEstimator )
			{
				printf("ID FAILURE!! %lf %lf\n", len, srdCollisionEstimator );	
		
				 theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  srdCollisionEstimator, vertex_data, ptr_to_data, nump );
			}
#endif
			distance[p] = srdCollisionEstimator;
			rp_dist[3*p+0] = rp[3*p+0];
			rp_dist[3*p+1] = rp[3*p+1];
			rp_dist[3*p+2] = rp[3*p+2];
		}
		else if( planar )
		{
			distance[p] = 0;
			rp_dist[3*p+0] = rp[3*p+0];
			rp_dist[3*p+1] = rp[3*p+1];
			rp_dist[3*p+2] = rp[3*p+2];
		}
		else
		{
			theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, -1 );
			rp_dist[3*p+0] = rp[3*p+0];
			rp_dist[3*p+1] = rp[3*p+1];
			rp_dist[3*p+2] = rp[3*p+2];

			double rpt[3], rnrm[3];

			theSurface->evaluateRNRM( f, u, v, rpt, rnrm, r );

			double dr[3] = { rp[3*p+0] - rpt[0],
					 rp[3*p+1] - rpt[1],
					 rp[3*p+2] - rpt[2] };

			double len = normalize(dr);
			
			double dp = dr[0] * rnrm[0] + dr[1] * rnrm[1] + dr[2] * rnrm[2];

//			printf("init len: %le dp: %le dist: %le\n", len, dp, distance[p] );

			if( !planar )
			{
				if( dp > 0 )
					inside_outside[p] = POINT_INSIDE;
				else
					inside_outside[p] = POINT_OUTSIDE;
			}
		}	
	}
#endif	
	free(vertex_data);
	free(ptr_to_data);
	free(nump);
}
*/
void srd_integrator::initializeDistances( Simulation *theSimulation,  double **M, int mlow, int mhigh )
{
#if 0 // pre-simulation, change me!
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );
		
	// find lowest power, this will limit our log search
	int min_p = grain_x_p;
	if( grain_y_p < min_p ) min_p = grain_y_p;
	if( grain_z_p < min_p ) min_p = grain_z_p;
	
	int nbins_id = id_grain_x * id_grain_y * id_grain_z;

	for( int x = 0; x < nbins_id; x++ )
	{
		bin_level[x] = -1; // unchecked.
		bins_id[x].np = 0;
	}
	for( int p = 0; p < np; p++ )
	{
		double fx = rp[3*p+0] / PBC_vec[0][0];
		while( fx < 0 ) fx += 1.0;
		while( fx >= 1.0 ) fx -= 1.0;
		int bx = fx * grain_x;
		if( bx == grain_x ) bx--;
		
		double fy = rp[3*p+1] / PBC_vec[1][1];
		while( fy < 0 ) fy += 1.0;
		while( fy >= 1.0 ) fy -= 1.0;
		int by = fy * grain_y;
		if( by == grain_y ) by--;
		
		double fz = rp[3*p+2] / PBC_vec[2][2];
		while( fz < 0 ) fz += 1.0;
		while( fz >= 1.0 ) fz -= 1.0;
		int bz = fz * grain_z;
		if( bz == grain_z ) bz--;

		int bin = (bx * id_grain_y + by ) * id_grain_z + bz;

		if( bins_id[bin].npSpace == bins_id[bin].np )
		{
			bins_id[bin].npSpace *= 2;
			bins_id[bin].list = (int *)realloc( bins_id[bin].list, sizeof(int) * bins_id[bin].npSpace );
		}

		bins_id[bin].list[bins_id[bin].np] = p;
		bins_id[bin].np++;
	}

	int *upper_level = NULL;

	int nradius = 0;
	int ntiny = 0;
	int nexplicit = 0;

	// this is the power search.
	for( int level = min_p-2; level >= 0; level-- )
	{
		int pfac = ipw2(level);

		int gx = id_grain_x / pfac; 
		int gy = id_grain_y / pfac; 
		int gz = id_grain_z / pfac; 

		int *our_level = (int *)malloc( sizeof(int) * gx * gy * gz );
		memset( our_level,0,sizeof(int) * gx*gy*gz);
		for( int bx = 0; bx < gx; bx++ )
		for( int by = 0; by < gy; by++ )
		for( int bz = 0; bz < gz; bz++ )
		{
			// do we need to check this or was it cleared at a higher level?

			if( upper_level )
			{
				int tbx = bx/2;
				int tby = by/2;
				int tbz = bz/2;

				if( upper_level[tbx*(gy/2)*(gz/2)+tby*(gz/2)+tbz] )
				{
//					printf("%d %d %d cleared at a higher level!\n", bx, by, bz );
					our_level[bx*gy*gz+by*gz+bz] = 1;
					continue; 
				}
			}

			double cen[3] = { (bx+0.5) * PBC_vec[0][0] / gx,
					  (by+0.5) * PBC_vec[1][1] / gy,
					  (bz+0.5) * PBC_vec[2][2] / gz };

			double LX = PBC_vec[0][0] / gx;
			double LY = PBC_vec[1][1] / gy;
			double LZ = PBC_vec[2][2] / gz;

		
			double maxr = sqrt( (LX*LX + LY*LY + LZ*LZ) / 4 );
	
			int f;
			double u,v;

//			printf("LEVEL %d cen %le %le %le\n", level, cen[0], cen[1], cen[2] );
			nradius++;
			if( ! theSurface->withinRadius( cen, &f, &u, &v,  M, mlow, mhigh, -1, maxr + srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
			{
				// cleared at this level.
				our_level[bx*gy*gz+by*gz+bz] = 1;
				int ncleared = 0;
				for( int tbx = bx * pfac; tbx < (bx+1)*pfac; tbx++ )
				for( int tby = by * pfac; tby < (by+1)*pfac; tby++ )
				for( int tbz = bz * pfac; tbz < (bz+1)*pfac; tbz++ )
				{
					int tbin = tbx*id_grain_y*id_grain_z+tby*id_grain_z+tbz;
					for( int bx = 0; bx < bins_id[tbin].np; bx++ )
					{
						int p = bins_id[tbin].list[bx];
						double dr[3] = { rp[3*p+0] - cen[0], rp[3*p+1] - cen[1], rp[3*p+2] - cen[2] };
						while( dr[0] < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0];
						while( dr[1] < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1];
						while( dr[2] < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2];
						while( dr[0] > PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0];
						while( dr[1] > PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1];
						while( dr[2] > PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];
						double r_cen = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
						distance[p] = (maxr - r_cen) + srdCollisionEstimator;
						ncleared++;
						rp_dist[3*p+0] = rp[3*p+0];
						rp_dist[3*p+1] = rp[3*p+1];
						rp_dist[3*p+2] = rp[3*p+2];					
					}
				}
//				printf("CLEARED! with %d particles.\n", ncleared );
			}
			else if( level == 0 ) 
			{
//				printf("failed to clear at lowest level.\n");
				int bin = bx*gy*gz+by*gz+bz;
				for( int px = 0; px < bins_id[bin].np; px++ )
				{
					int p = bins_id[bin].list[px];

					ntiny++;
					if( ! theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1, maxr + srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
					{
						rp_dist[3*p+0] = rp[3*p+0];
						rp_dist[3*p+1] = rp[3*p+1];
						rp_dist[3*p+2] = rp[3*p+2];
						distance[p] = srdCollisionEstimator;
					}
					else
					{
						nexplicit++;
						theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, -1 );
						rp_dist[3*p+0] = rp[3*p+0];
						rp_dist[3*p+1] = rp[3*p+1];
						rp_dist[3*p+2] = rp[3*p+2];
	
						double rpt[3], rnrm[3];
	
						theSurface->evaluateRNRM( f, u, v, rpt, rnrm, r );
	
						double dr[3] = { rp[3*p+0] - rpt[0],
							 rp[3*p+1] - rpt[1],
							 rp[3*p+2] - rpt[2] };
		
						double len = normalize(dr);
					
						double dp = dr[0] * rnrm[0] + dr[1] * rnrm[1] + dr[2] * rnrm[2];

//			printf("init len: %le dp: %le dist: %le\n", len, dp, distance[p] );

						if( !planar )
						{
							if( dp > 0 )
								inside_outside[p] = POINT_INSIDE;
							else
								inside_outside[p] = POINT_OUTSIDE;
						}
					}
				}
			}
			else
			{
//				printf("...failed to clear.\n");
			}
		}

		if( upper_level )
			free(upper_level);
		upper_level = our_level;
	}  

	printf("Ran %d radius checks.\n", nradius );
	printf("Ran %d tiny radius checks.\n", ntiny );
	printf("Ran explicit check on %d/%d particles.\n", nexplicit, np );
	

/*
	for( int bx = 0; bx < id_grain_x; bx++ )
	for( int by = 0; by < id_grain_y; by++ )
	for( int bz = 0; bz < id_grain_z; bz++ )
	{
		int b = (bx * grain_y + by ) * grain_z + bz;

		

		double cen[3] = {  
			(bx+0.5) * PBC_vec[0][0] / grain_x,
			(by+0.5) * PBC_vec[1][1] / grain_y,
			(bz+0.5) * PBC_vec[2][2] / grain_z };

		if( bins[b].np == 0 ) continue;

		double r_max = 0;
		double dr[bins[b].np];

		for( int px = 0; px < bins[b].np; px++ )
		{
			int p = bins[b].list[px];

			double dr[3] = { rp[3*p+0] - cen[0],
					 rp[3*p+1] - cen[1],
					 rp[3*p+2] - cen[2] };

			while( dr[0] < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0];
			while( dr[1] < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1];
			while( dr[2] < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2];
			while( dr[0] >  PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0];
			while( dr[1] >  PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1];
			while( dr[2] >  PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];

			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r > r_max )
				r_max = r;

			dr[px] = r;
					
		}
		
		int f;
		double u,v;
		printf("FINAL LEVEL Checking %d particles near %lf %lf %lf ", bins[b].np, cen[0], cen[1], cen[2] );
		if( ! theSurface->withinRadius( cen, &f, &u, &v,  M, mlow, mhigh, -1,  r_max + srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
		{
			printf("... all are clear.\n");
			for( int px = 0; px < bins[b].np; px++ )
			{
				int p = bins[b].list[px];
				distance[p] = (r_max - dr[px]) + srdCollisionEstimator;

				rp_dist[3*p+0] = rp[3*p+0];
				rp_dist[3*p+1] = rp[3*p+1];
				rp_dist[3*p+2] = rp[3*p+2];				
			}
		}	
		else
		{
			printf("...checking them all...");
			double nok = 0, nbad=0;
			for( int px = 0; px < bins[b].np; px++ )
			{
				int p = bins[b].list[px];
			
				if( ! theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
				{
					nok += 1;
					distance[p] = srdCollisionEstimator;

					rp_dist[3*p+0] = rp[3*p+0];
					rp_dist[3*p+1] = rp[3*p+1];
					rp_dist[3*p+2] = rp[3*p+2];				
				}
				else
				{
					nbad += 1;
					theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, -1 );
					rp_dist[3*p+0] = rp[3*p+0];
					rp_dist[3*p+1] = rp[3*p+1];
					rp_dist[3*p+2] = rp[3*p+2];

					double rpt[3], rnrm[3];

					theSurface->evaluateRNRM( f, u, v, rpt, rnrm, r );

					double dr[3] = { rp[3*p+0] - rpt[0],
						 rp[3*p+1] - rpt[1],
						 rp[3*p+2] - rpt[2] };
	
					double len = normalize(dr);
				
					double dp = dr[0] * rnrm[0] + dr[1] * rnrm[1] + dr[2] * rnrm[2];

//			printf("init len: %le dp: %le dist: %le\n", len, dp, distance[p] );

					if( !planar )
					{
						if( dp > 0 )
							inside_outside[p] = POINT_INSIDE;
						else
							inside_outside[p] = POINT_OUTSIDE;
					}
				}
				printf(" %lf%% were clear.\n", 100 * nok / (nok+nbad) );
			}
		}
	}
*/
	free(vertex_data);
	free(ptr_to_data);
	free(nump);
#endif
}

int srd_integrator::stream_and_collide( Simulation *theSimulation, double **M, int mlow, int mhigh, double del_hull,  double running_time, double integrate_time, double time_step_collision )
{
	int ncol = 0;
#if 0  // pre-simulation change me!!


	double force_factor = time_step_collision / integrate_time;

	double average_velocity = sqrt((temp*3/mass));
	double time_free_path_grain = a / average_velocity;		
	// in AKMA time.		

	double dt = time_free_path_grain;

	if( dt > time_step_collision )
		dt = time_step_collision;

	double use_running_time = running_time;
	if( use_running_time/time_step_collision > 1e9 )
		use_running_time -= time_step_collision * floor(use_running_time/time_step_collision);

	int icol = use_running_time / time_step_collision;

//	printf("Free path time: %le collision time: %le running_time: %le\n",
//		time_free_path_grain, time_step_collision, running_time );
		
	setMembraneMinMax(r,theSurface);

	int iter = 0;
	for( double t = use_running_time; t < use_running_time + integrate_time - 1e-14; t += dt, iter++)
	{
		double use_dt = dt;

		if( t + dt > use_running_time+integrate_time )
			use_dt = use_running_time+integrate_time - t;
		
		double path_expec = SAFETY_FACTOR * average_velocity * use_dt;
//		printf("path expec: %le\n", path_expec );
		
		// refreshes the tags on particles that have moved into range of the membrane.		
		tagParticlesForCollision( r, theSurface, path_expec+del_hull, M, mlow, mhigh );
		stream(use_dt);


		int ncol_local = 0;				

		ncol_local = resolveCollision( r, g, vmem, effm, path_expec+del_hull, theSurface, M, mlow, mhigh, theForceSet, use_dt, force_factor );
		ncol += ncol_local;

		if( icol != (int)use_running_time/time_step_collision )
		{
			icol = (int)use_running_time/time_step_collision;
			
			collide();
		}
	}

#endif
	return ncol;

}

int srd_integrator::collide_with_membrane( double *r, double *g, surface *theSurface, double **M, int mlow, int mhigh, double del_hull, force_set *theForceSet, double time_step, double force_factor, double *vmem, SparseMatrix *effm )
{
	memset( did_collide, 0, sizeof(int) * np );
	if( planar )
		return collide_with_membrane_planar( r, g, theSurface, M, mlow, mhigh, del_hull, theForceSet, time_step, force_factor, vmem, effm );
	else
		return collide_with_membrane_inside_outside( r, g, theSurface, M, mlow, mhigh, del_hull, theForceSet, time_step, force_factor,vmem, effm );
}

static int col_cntr = 0; // debugging only 

int srd_integrator::collide_with_membrane_planar( double *r, double *g, surface *theSurface, double **M, int mlow, int mhigh, double del_hull, force_set *theForceSet, double time_step, double force_factor, double *vmem, SparseMatrix *effm  )
{
	// this is the vector of forces projected onto the control points ( the generalized coordinates )
	// d v / d xi * dxi / dcontrol_pt. The change in control points giving the proper change in velocities must be solved at the end
	// the solution is a linear problem using the control point inverse matrix.
 
	double *force_vector = (double *)malloc( sizeof(double) * theSurface->nv * 3 );
	memset( force_vector, 0, sizeof(double) * theSurface->nv * 3 );
	double total_dp[3] = { 0,0,0};
	
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );

	int n_col = 0;

	int n_col_err = 0;
	for( int p = 0; p < np; p++ )
	{
		double dmove[3] = { 
			rp[3*p+0] - rp_dist[3*p+0],	
			rp[3*p+1] - rp_dist[3*p+1],	
			rp[3*p+2] - rp_dist[3*p+2] };

		while( dmove[0] < -PBC_vec[0][0]/2 ) dmove[0] += PBC_vec[0][0];	
		while( dmove[1] < -PBC_vec[1][1]/2 ) dmove[1] += PBC_vec[1][1];	
		while( dmove[2] < -PBC_vec[2][2]/2 ) dmove[2] += PBC_vec[2][2];
	
		while( dmove[0] >  PBC_vec[0][0]/2 ) dmove[0] -= PBC_vec[0][0];	
		while( dmove[1] >  PBC_vec[1][1]/2 ) dmove[1] -= PBC_vec[1][1];	
		while( dmove[2] >  PBC_vec[2][2]/2 ) dmove[2] -= PBC_vec[2][2];

		double r_move = sqrt(dmove[0]*dmove[0]+dmove[1]*dmove[1]+dmove[2]*dmove[2]);

//		printf("r_move: %le distance: %le\n", r_move, distance[p] );

		if( r_move > distance[p] - del_hull  )
		{
			int f;
			double u,v;
			double dist;

			double rad = -1;

			// = theSurface->returnRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh,  r_move, distance[p],  -1,  srdCollisionEstimator, inside_outside[p] );

#ifdef WITHIN_CHECK
			if( ! theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
				rad = srdCollisionEstimator;		
	
			rp_dist[3*p+0] = rp[3*p+0];		
			rp_dist[3*p+1] = rp[3*p+1];		
			rp_dist[3*p+2] = rp[3*p+2];		

			if( rad >= 0 )	
				distance[p] = rad;
			else
				distance[p] =  0;
#endif

			int col_f;
			double col_u, col_v;

			int collide = theSurface->linearCollisionPoint( rp_prev+3*p, rp+3*p, &col_f, &col_u, &col_v, M, mlow, mhigh );

			if( collide )
			{
				if( p == 1436 )
				{
					printf("check this.\n");
				}

				did_collide[p] = 1;
				double rcol[3], nrm[3];
		
				theSurface->evaluateRNRM( col_f, col_u, col_v, rcol, nrm, r );
				
				double vp_proj = (vp[3*p+0]) * nrm[0] + (vp[3*p+1]) * nrm[1] + (vp[3*p+2]) * nrm[2];
	
				double mem_v[3];
				theSurface->velocityAtPoint( col_f, col_u, col_v, vp, mem_v ); 

/*
				The force will be directed along the membrane normal.
				We are solving for the magnitude of the force, x.
				The constraints are conservation of momentum and conservation of kinetic energy.

				compute the first and second derivatives of the KE.

				the force is the change in momentum of the particle divide by the timestep.			
*/
				// the ratio of the magnitude of the force to the change in momentum of the particle.

				double df_dp = 1.0;//force_factor;
				double dKE_mem_dx = 0;
				double d2KE_mem_dx2 = 0;

				double dKE_p_dx = vp_proj; // mass weighted 
				double d2KE_p_dx2 = 1.0 / mass;

				theSurface->dKE_dx_and2( theForceSet, effm, vmem, nrm, col_f, col_u, col_v, &dKE_mem_dx, &d2KE_mem_dx2 );
//				theSurface->debug_dKE_dx_and2( theForceSet, effm, vmem, nrm, col_f, col_u, col_v );
				dKE_mem_dx *= df_dp;
				d2KE_mem_dx2 *= df_dp * df_dp;

				double mag = - 2* (dKE_p_dx + dKE_mem_dx) / (d2KE_p_dx2 + d2KE_mem_dx2); 


				double dKE_p = dKE_p_dx * mag + d2KE_p_dx2 * 0.5 * mag*mag;
				double dKE_mem = -dKE_mem_dx * mag + d2KE_mem_dx2 * 0.5 * mag*mag;
		
				dKE_out_collision += dKE_mem;

	
				// modify the position and velocity of the colliding particle.

				double proj_out_move[3] = {rp[3*p+0] - rp_prev[3*p+0],
							   rp[3*p+1] - rp_prev[3*p+1], 
							   rp[3*p+2] - rp_prev[3*p+2] };
				//change in momentum
				double rp_proj =  proj_out_move[0] * nrm[0] + proj_out_move[1] * nrm[1] + proj_out_move[2] * nrm[2];

				/*
					vp_proj is the velocity projected along the normal. the solution to
					constant momentum and kinetic energy is

					dvp_nrm = - 2 * membrane_point_mass * vnrm / (membrane_point_mass + mass)
				        dvm_nrm = - 2 * mp * vnrm / ( membrane_point_mass + mass)

				*/

				double sign_change = -1.0;

				rp[3*p+0] = rcol[0] - rp_proj * nrm[0] / 2;
				rp[3*p+1] = rcol[1] - rp_proj * nrm[1] / 2;
				rp[3*p+2] = rcol[2] - rp_proj * nrm[2] / 2;

				//change in momentum

				double dp[3] = {0,0,0};

				dp[0] = mag * nrm[0];
				dp[1] = mag * nrm[1];
				dp[2] = mag * nrm[2];

				double ke_before = 0.5 * mass * (vp[3*p+0]*vp[3*p+0]+vp[3*p+1]*vp[3*p+1]+vp[3*p+2]*vp[3*p+2]);


				vp[3*p+0] += dp[0] / mass;
				vp[3*p+1] += dp[1] / mass;
				vp[3*p+2] += dp[2] / mass;

				double ke_after = 0.5 * mass * (vp[3*p+0]*vp[3*p+0]+vp[3*p+1]*vp[3*p+1]+vp[3*p+2]*vp[3*p+2]);


				// apply force to the membrane.

				int vert,edge;
				// getting the vertex and edge of a face is a convenient way to get the local coordinate system.
				if( col_f < theSurface->nf_faces )
				{
					vert = theSurface->theFormulas[col_f*theSurface->nf_g_q_p].vertex;
					edge = theSurface->theFormulas[col_f*theSurface->nf_g_q_p].edge;
				}
				else	
				{
					vert = theSurface->theIrregularFormulas[(col_f-theSurface->nf_faces)*theSurface->nf_irr_pts].vertex;
					edge = theSurface->theIrregularFormulas[(col_f-theSurface->nf_faces)*theSurface->nf_irr_pts].edge;
				}

				int j = theSurface->theVertices[vert].edges[edge];
				int ep1 = edge+1;
				if( ep1 >= theSurface->theVertices[vert].valence )
					ep1 -= theSurface->theVertices[vert].valence;
				int k = theSurface->theVertices[vert].edges[ep1];

				// I want a particular momentum spread on the face.
				// I can do a least-squares fit using the inverse operator.
				// this appears to be equivalent to the effective-mass matrix inverse, so I'm not sure at this point.

				theSurface->applyForceAtPoint( col_f, col_u, col_v, dp, force_vector, theForceSet );

				// at this point, dp is in the direction of the particle's momentum/velocity change. we will change in the opposite direction	
				for( int v = 0; v < theSurface->nv; v++ )
				{
					force_vector[3*v+0] *= -1;
					force_vector[3*v+1] *= -1;
					force_vector[3*v+2] *= -1;
				}

				total_dp[0] += dp[0];
				total_dp[1] += dp[1];
				total_dp[2] += dp[2];
				n_col++;
				col_cntr++;
			}

		}
	}		

	int nv = theSurface->nv;

	for( int v1 = 0; v1 < nv; v1++ )
	{
		// force_factor is the ratio of the membrane time step (over which the force is integrated) to the SRD time step

		// gradient is the negative of the force, so a minus sign here.
		g[3*v1+0] -= force_factor * force_vector[3*v1+0] / time_step;		
		g[3*v1+1] -= force_factor * force_vector[3*v1+1] / time_step;	
		g[3*v1+2] -= force_factor * force_vector[3*v1+2] / time_step;		
	}

	free(force_vector);
	free(vertex_data);
	free(ptr_to_data);
	free(nump);

	return n_col;
}

int srd_integrator::collide_with_membrane_inside_outside( double *r, double *g, surface *theSurface, double **M, int mlow, int mhigh, double del_hull, force_set *theForceSet, double time_step, double force_factor, double *vmem, SparseMatrix *effm )
{

	printf("This needs to be updated to reflect the proper collision.\n");
	exit(1);

	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );
	
	// this is the vector of forces projected onto the control points ( the generalized coordinates )
	// d v / d xi * dxi / dcontrol_pt. The change in control points giving the proper change in velocities must be solved at the end
	// the solution is a linear problem using the control point inverse matrix.
 
	double *force_vector = (double *)malloc( sizeof(double) * theSurface->nv * 3 );
	memset( force_vector, 0, sizeof(double) * theSurface->nv * 3 );
	double total_dp[3] = { 0,0,0};

	int n_col = 0;

	int n_col_err = 0;
	for( int p = 0; p < np; p++ )
	{
		double dmove[3] = { 
			rp[3*p+0] - rp_dist[3*p+0],	
			rp[3*p+1] - rp_dist[3*p+1],	
			rp[3*p+2] - rp_dist[3*p+2] };

		while( dmove[0] < -PBC_vec[0][0]/2 ) dmove[0] += PBC_vec[0][0];	
		while( dmove[1] < -PBC_vec[1][1]/2 ) dmove[1] += PBC_vec[1][1];	
		while( dmove[2] < -PBC_vec[2][2]/2 ) dmove[2] += PBC_vec[2][2];
	
		while( dmove[0] >  PBC_vec[0][0]/2 ) dmove[0] -= PBC_vec[0][0];	
		while( dmove[1] >  PBC_vec[1][1]/2 ) dmove[1] -= PBC_vec[1][1];	
		while( dmove[2] >  PBC_vec[2][2]/2 ) dmove[2] -= PBC_vec[2][2];

		double r_move = sqrt(dmove[0]*dmove[0]+dmove[1]*dmove[1]+dmove[2]*dmove[2]);

//		printf("r_move: %le distance: %le\n", r_move, distance[p] );

		if( inside_outside[p] == POINT_UNKNOWN && distance[p] < srdCollisionEstimator )
		{
			int f;
			double u, v;

			theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, srdCollisionEstimator );
			rp_dist[3*p+0] = rp[3*p+0];
			rp_dist[3*p+1] = rp[3*p+1];
			rp_dist[3*p+2] = rp[3*p+2];

			double rpt[3], rnrm[3];

			theSurface->evaluateRNRM( f, u, v, rpt, rnrm, r );

			double dr[3] = { rp[3*p+0] - rpt[0],
					 rp[3*p+1] - rpt[1],
					 rp[3*p+2] - rpt[2] };


			while( dr[0] < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0];	
			while( dr[1] < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1];	
			while( dr[2] < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2];
	
			while( dr[0] >  PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0];	
			while( dr[1] >  PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1];	
			while( dr[2] >  PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];
			double len = normalize(dr);
			
			double dp = dr[0] * rnrm[0] + dr[1] * rnrm[1] + dr[2] * rnrm[2];

//			printf("len: %le dp: %le dist: %le\n", len, dp, distance[p] );

			if( dp > 0 )
				inside_outside[p] = POINT_INSIDE;
			else
				inside_outside[p] = POINT_OUTSIDE;	
		}


		if( r_move > distance[p] - del_hull  )
		{
			int f;
			double u,v;
			double dist;

			double rad = 0;

			rad = theSurface->returnRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh,  r_move, distance[p],  -1,  srdCollisionEstimator, vertex_data, ptr_to_data, nump, inside_outside[p] );
			if( rad < FUDGE_DISTANCE )
			{
				// The particle has collided with the surface.
				// we need to find a point on the membrane, reflect the particle,
				// and apply a collision force to the membrane with the correct momentum transfer.

				int col_f;
				double col_u, col_v;

				double endp[3]   = { rp[3*p+0], rp[3*p+1], rp[3*p+2] };
				double startp[3] = { rp_prev[3*p+0], rp_prev[3*p+1], rp_prev[3*p+2] };

				double dr[3] = { endp[0] - startp[0], endp[1] - startp[1], endp[2] - startp[2] };

				double len = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

				if( len < MEMBRANE_LEVEL_TOLERANCE )
				{
					double midp[3] = { (startp[0] +endp[0])/2, (startp[1]+endp[1])/2, (startp[2]+endp[2])/2 };
					normalize(dr);
					startp[0] = midp[0] - dr[0] * MEMBRANE_LEVEL_TOLERANCE/2;
					startp[1] = midp[1] - dr[1] * MEMBRANE_LEVEL_TOLERANCE/2;
					startp[2] = midp[2] - dr[2] * MEMBRANE_LEVEL_TOLERANCE/2;
					endp[0] = midp[0] + dr[0] * MEMBRANE_LEVEL_TOLERANCE/2;
					endp[1] = midp[1] + dr[1] * MEMBRANE_LEVEL_TOLERANCE/2;
					endp[2] = midp[2] + dr[2] * MEMBRANE_LEVEL_TOLERANCE/2;
				}

				int collide = theSurface->linearCollisionPoint( startp, endp, &col_f, &col_u, &col_v, M, mlow, mhigh );

				if( collide )
				{
//					printf("%d COLLISION !!!!!!\n", p);	
					double rcol[3], nrm[3];
			
					theSurface->evaluateRNRM( col_f, col_u, col_v, rcol, nrm, r );
		
					// modify the position and velocity of the colliding particle.
	
					double proj_out_move[3] = {rp[3*p+0] - rp_prev[3*p+0],
								   rp[3*p+1] - rp_prev[3*p+1], 
								   rp[3*p+2] - rp_prev[3*p+2] };
					//change in momentum
					double rp_proj =  proj_out_move[0] * nrm[0] + proj_out_move[1] * nrm[1] + proj_out_move[2] * nrm[2];
					double vp_proj = vp[3*p+0] * nrm[0] + vp[3*p+1] * nrm[1] + vp[3*p+2] * nrm[2];
					
					double sign_change = -1.0;
	
//					if( inside_outside[p] == POINT_INSIDE )
//						sign_change = -1;

					// MAKE SURE IT STAYS INSIDE!!!
	
					if( inside_outside[p] == POINT_INSIDE )
					{
						if( rp_proj > 0 )
						{ // a displacement defect.
							rp[3*p+0] = rcol[0] + fabs(rp_proj) * nrm[0] / 2;
							rp[3*p+1] = rcol[1] + fabs(rp_proj) * nrm[1] / 2;
							rp[3*p+2] = rcol[2] + fabs(rp_proj) * nrm[2] / 2;
						}
						else
						{
							rp[3*p+0] = rcol[0] - rp_proj * nrm[0] / 2;
							rp[3*p+1] = rcol[1] - rp_proj * nrm[1] / 2;
							rp[3*p+2] = rcol[2] - rp_proj * nrm[2] / 2;
						}
					}
					else
					{
						if( rp_proj < 0 )
						{ // I should find out more about these defects but they do seem rare...
							rp[3*p+0] = rcol[0] - fabs(rp_proj) * nrm[0] / 2;
							rp[3*p+1] = rcol[1] - fabs(rp_proj) * nrm[1] / 2;
							rp[3*p+2] = rcol[2] - fabs(rp_proj) * nrm[2] / 2;
						}
						else
						{
							rp[3*p+0] = rcol[0] - rp_proj * nrm[0] / 2;
							rp[3*p+1] = rcol[1] - rp_proj * nrm[1] / 2;
							rp[3*p+2] = rcol[2] - rp_proj * nrm[2] / 2;
						}
					}

					//change in momentum
	
					double dp[3] = {0,0,0};
					
					dp[0] = sign_change*2 *vp_proj * nrm[0] * mass;
					dp[1] = sign_change*2 *vp_proj * nrm[1] * mass;
					dp[2] = sign_change*2 *vp_proj * nrm[2] * mass;
	
					vp[3*p+0] += sign_change*2*vp_proj * nrm[0];
					vp[3*p+1] += sign_change*2*vp_proj * nrm[1];
					vp[3*p+2] += sign_change*2*vp_proj * nrm[2];
				
					// apply force to the membrane.
				
					int vert,edge;
					// getting the vertex and edge of a face is a convenient way to get the local coordinate system.
					if( f < theSurface->nf_faces )
					{
						vert = theSurface->theFormulas[f*theSurface->nf_g_q_p].vertex;
						edge = theSurface->theFormulas[f*theSurface->nf_g_q_p].edge;
					}
					else	
					{
						vert = theSurface->theIrregularFormulas[(f-theSurface->nf_faces)*theSurface->nf_irr_pts].vertex;
						edge = theSurface->theIrregularFormulas[(f-theSurface->nf_faces)*theSurface->nf_irr_pts].edge;
					}

					int j = theSurface->theVertices[vert].edges[edge];
					int ep1 = edge+1;
					if( ep1 >= theSurface->theVertices[vert].valence )
						ep1 -= theSurface->theVertices[vert].valence;
					int k = theSurface->theVertices[vert].edges[ep1];

					// I want a particular momentum spread on the face.
					// I can do a least-squares fit using the inverse operator.
					// this appears to be equivalent to the effective-mass matrix inverse, so I'm not sure at this point.
		
					theSurface->applyForceAtPoint( col_f, col_u, col_v, dp, force_vector, theForceSet );
					total_dp[0] += dp[0];
					total_dp[1] += dp[1];
					total_dp[2] += dp[2];
					n_col++;
				}

				distance[p] = 0;
			}
			else
				distance[p] = rad;

			rp_dist[3*p+0] = rp[3*p+0];		
			rp_dist[3*p+1] = rp[3*p+1];		
			rp_dist[3*p+2] = rp[3*p+2];		
		}
	}		

	int nv = theSurface->nv;

	for( int v1 = 0; v1 < nv; v1++ )
	{
		g[3*v1+0] += force_factor * force_vector[3*v1+0] / time_step;		
		g[3*v1+1] += force_factor * force_vector[3*v1+1] / time_step;	
		g[3*v1+2] += force_factor * force_vector[3*v1+2] / time_step;		
	}

#ifdef DEBUG_MOMENTUM_CHANGE
	if( n_col > 0 )
		printf("Collision total solvent momentum change: %le %le %le\n", total_dp[0], total_dp[1], total_dp[2] ); 
#endif
	free(force_vector);

	free(vertex_data);
	free(ptr_to_data);
	free(nump);

	return n_col;
}

void srd_integrator::writeXYZ( FILE *theXYZ )
{	
	if( theXYZ )
	{
		fprintf(theXYZ, "%d\n", np );
		fprintf(theXYZ, "SRD virtual solvent particles\n");
		for( int p = 0; p < np; p++ )
		{
			if( inside_outside[p] == POINT_INSIDE )
				fprintf(theXYZ, "C %lf %lf %lf\n", rp[3*p+0], rp[3*p+1], rp[3*p+2] );
			else if( inside_outside[p] == POINT_OUTSIDE )
				fprintf(theXYZ, "O %lf %lf %lf\n", rp[3*p+0], rp[3*p+1], rp[3*p+2] );
			else
				fprintf(theXYZ, "N %lf %lf %lf\n", rp[3*p+0], rp[3*p+1], rp[3*p+2] );
		}
		fflush(theXYZ);
	}
}	

void srd_integrator::clear( void )
{
	int nbins = grain_x * grain_y * grain_z ;

	for( int x = 0; x < nbins; x++ )
		free(bins[x].list);
	free( bins);

	free(rp);
	free(rp_prev);
	free(rp_dist);
	free(distance);
	free(vp);
	free(inside_outside);
}

//hi
/*
void srd_integrator::tagParticlesForCollision( double * r, surface * theSurface, double delta_hull_collision, double **M, int mlow, int mhigh )
{
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );

	double search_radius = 50.0;
	//double search_radius = delta_hull_collision;

	double LA = PBC_vec[0][0];
	double LB = PBC_vec[1][1];
	double LC = PBC_vec[2][2];

	int nnear = 0;
	for( int p = 0; p < np; p++ )
	{
		// only re-tag particles if their tag is zero

//		if( fabs(topological_tag[p]) < 1e-14 )
		{
			int f;
			double u, v;
	
			int exclude = 0;

			for( int c = 0; c < 3; c++ )
			{
				if( rp[3*p+c] < membrane_min[c] - delta_hull_collision ) exclude = 1;
				if( rp[3*p+c] > membrane_max[c] + delta_hull_collision ) exclude = 1;
			}
		
			if( exclude ) continue;
	
			int did_col = theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1, delta_hull_collision, vertex_data, ptr_to_data, nump ); 
				
			if( did_col )
			{
				nnear += 1;
				theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, delta_hull_collision );
				rp_dist[3*p+0] = rp[3*p+0];
				rp_dist[3*p+1] = rp[3*p+1];
				rp_dist[3*p+2] = rp[3*p+2];

				double rnrm[3], nrm[3];

				theSurface->evaluateRNRM(f, u, v, rnrm, nrm, r ); 

				double dr[3] = { rp[3*p+0] - rnrm[0],
						 rp[3*p+1] - rnrm[1],
						 rp[3*p+2] - rnrm[2] };

				while( dr[0] < -LA/2 ) dr[0] += LA;
				while( dr[0] >  LA/2 ) dr[0] -= LA;
				while( dr[1] < -LB/2 ) dr[1] += LB;
				while( dr[1] >  LB/2 ) dr[1] -= LB;
				while( dr[2] < -LC/2 ) dr[2] += LC;
				while( dr[2] >  LC/2 ) dr[2] -= LC;

				double dp = dr[0] * nrm[0] + dr[1] * nrm[1] + dr[2] * nrm[2];
				
				if( dp < 0 )
					topological_tag[p] = dp;
				else
					topological_tag[p] = dp;

				if( last_known_tag[p] * topological_tag[p] < 0 )
				{
					printf("ABSOLUTELY MISSED COLLISION %d tag1 %lf tag2 %lf\n", p, last_known_tag[p], topological_tag[p]);
				}

				last_known_tag[p] = dp;
			}
			else
				last_known_tag[p] = 0;
		} 
//		else
//			last_known_tag[p] = 0;
	}

	free(vertex_data);
	free(ptr_to_data);
	free(nump);
}
*/

void srd_integrator::tagParticlesForCollision( double * r, surface * theSurface, double delta_hull_collision, double **M, int mlow, int mhigh )
{
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );
		
	// find lowest power, this will limit our log search
	int min_p = grain_x_p;
	if( grain_y_p < min_p ) min_p = grain_y_p;
	if( grain_z_p < min_p ) min_p = grain_z_p;
	
	int nbins_id = id_grain_x * id_grain_y * id_grain_z;

	for( int x = 0; x < nbins_id; x++ )
	{
		bin_level[x] = -1; // unchecked.
		bins_id[x].np = 0;
	}
	for( int p = 0; p < np; p++ )
	{
		double fx = rp[3*p+0] / PBC_vec[0][0];
		while( fx < 0 ) fx += 1.0;
		while( fx >= 1.0 ) fx -= 1.0;
		int bx = fx * id_grain_x;
		if( bx == id_grain_x ) bx--;
		
		double fy = rp[3*p+1] / PBC_vec[1][1];
		while( fy < 0 ) fy += 1.0;
		while( fy >= 1.0 ) fy -= 1.0;
		int by = fy * id_grain_y;
		if( by == id_grain_y ) by--;
		
		double fz = rp[3*p+2] / PBC_vec[2][2];
		while( fz < 0 ) fz += 1.0;
		while( fz >= 1.0 ) fz -= 1.0;
		int bz = fz * id_grain_z;
		if( bz == id_grain_z ) bz--;

		int bin = (bx * id_grain_y + by ) * id_grain_z + bz;

		if( bins_id[bin].npSpace == bins_id[bin].np )
		{
			bins_id[bin].npSpace *= 2;
			bins_id[bin].list = (int *)realloc( bins_id[bin].list, sizeof(int) * bins_id[bin].npSpace );
		}

		bins_id[bin].list[bins_id[bin].np] = p;
		bins_id[bin].np++;
	}

	int *upper_level = NULL;

	int nradius = 0;
	int ntiny = 0;
	int nexplicit = 0;

	double LA = PBC_vec[0][0];
	double LB = PBC_vec[1][1];
	double LC = PBC_vec[2][2];

	// this is the power search.
	for( int level = min_p-2; level >= 0; level-- )
	{
		int pfac = ipw2(level);

		int gx = id_grain_x / pfac; 
		int gy = id_grain_y / pfac; 
		int gz = id_grain_z / pfac; 

		int *our_level = (int *)malloc( sizeof(int) * gx * gy * gz );
		memset( our_level,0,sizeof(int) * gx*gy*gz);
		for( int bx = 0; bx < gx; bx++ )
		for( int by = 0; by < gy; by++ )
		for( int bz = 0; bz < gz; bz++ )
		{
			// do we need to check this or was it cleared at a higher level?

			if( upper_level )
			{
				int tbx = bx/2;
				int tby = by/2;
				int tbz = bz/2;

				if( upper_level[tbx*(gy/2)*(gz/2)+tby*(gz/2)+tbz] )
				{
//					printf("%d %d %d cleared at a higher level!\n", bx, by, bz );
					our_level[bx*gy*gz+by*gz+bz] = 1;
					continue; 
				}
			}

			double cen[3] = { (bx+0.5) * PBC_vec[0][0] / gx,
					  (by+0.5) * PBC_vec[1][1] / gy,
					  (bz+0.5) * PBC_vec[2][2] / gz };

			double LX = PBC_vec[0][0] / gx;
			double LY = PBC_vec[1][1] / gy;
			double LZ = PBC_vec[2][2] / gz;

		
			double maxr = sqrt( (LX*LX + LY*LY + LZ*LZ) / 4 );
	
			int f;
			double u,v;

//			printf("LEVEL %d cen %le %le %le\n", level, cen[0], cen[1], cen[2] );
			nradius++;
			if( ! theSurface->withinRadius( cen, &f, &u, &v,  M, mlow, mhigh, -1, maxr + srdCollisionEstimator + delta_hull_collision, vertex_data, ptr_to_data, nump ) )
			{
				// cleared at this level.
				our_level[bx*gy*gz+by*gz+bz] = 1;
				int ncleared = 0;
				for( int tbx = bx * pfac; tbx < (bx+1)*pfac; tbx++ )
				for( int tby = by * pfac; tby < (by+1)*pfac; tby++ )
				for( int tbz = bz * pfac; tbz < (bz+1)*pfac; tbz++ )
				{
					int tbin = tbx*id_grain_y*id_grain_z+tby*id_grain_z+tbz;
					for( int bx = 0; bx < bins_id[tbin].np; bx++ )
					{
						int p = bins_id[tbin].list[bx];
						double dr[3] = { rp[3*p+0] - cen[0], rp[3*p+1] - cen[1], rp[3*p+2] - cen[2] };
						while( dr[0] < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0];
						while( dr[1] < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1];
						while( dr[2] < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2];
						while( dr[0] > PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0];
						while( dr[1] > PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1];
						while( dr[2] > PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];
						double r_cen = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
						distance[p] = (maxr - r_cen) + srdCollisionEstimator + delta_hull_collision;
						ncleared++;
						rp_dist[3*p+0] = rp[3*p+0];
						rp_dist[3*p+1] = rp[3*p+1];
						rp_dist[3*p+2] = rp[3*p+2];					
					}
				}
//				printf("CLEARED! with %d particles.\n", ncleared );
			}
			else if( level == 0 ) 
			{
//				printf("failed to clear at lowest level.\n");
				int bin = bx*gy*gz+by*gz+bz;
				for( int px = 0; px < bins_id[bin].np; px++ )
				{
					int p = bins_id[bin].list[px];

					ntiny++;
					if( ! theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  srdCollisionEstimator + delta_hull_collision, vertex_data, ptr_to_data, nump ) )
					{
						rp_dist[3*p+0] = rp[3*p+0];
						rp_dist[3*p+1] = rp[3*p+1];
						rp_dist[3*p+2] = rp[3*p+2];
						distance[p] = srdCollisionEstimator;
					}
					else
					{
						nexplicit++;
						theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, -1 );

						rp_dist[3*p+0] = rp[3*p+0];
						rp_dist[3*p+1] = rp[3*p+1];
						rp_dist[3*p+2] = rp[3*p+2];
		
						double rnrm[3], nrm[3];
		
						theSurface->evaluateRNRM(f, u, v, rnrm, nrm, r ); 
		
						double dr[3] = { rp[3*p+0] - rnrm[0],
								 rp[3*p+1] - rnrm[1],
								 rp[3*p+2] - rnrm[2] };
		
						while( dr[0] < -LA/2 ) dr[0] += LA;
						while( dr[0] >  LA/2 ) dr[0] -= LA;
						while( dr[1] < -LB/2 ) dr[1] += LB;
						while( dr[1] >  LB/2 ) dr[1] -= LB;
						while( dr[2] < -LC/2 ) dr[2] += LC;
						while( dr[2] >  LC/2 ) dr[2] -= LC;
		
						double dp = dr[0] * nrm[0] + dr[1] * nrm[1] + dr[2] * nrm[2];
						
						if( dp < 0 )
							topological_tag[p] = dp;
						else
							topological_tag[p] = dp;
		
						if( last_known_tag[p] * topological_tag[p] < 0 )
						{
//							printf("ABSOLUTELY MISSED COLLISION %d tag1 %lf tag2 %lf\n", p, last_known_tag[p], topological_tag[p]);
						}
		
						last_known_tag[p] = dp;
					}
				}
			}
			else
			{
//				printf("...failed to clear.\n");
			}
		}

		if( upper_level )
			free(upper_level);
		upper_level = our_level;
	}  

	printf("tag Ran %d radius checks.\n", nradius );
	printf("tag Ran %d tiny radius checks.\n", ntiny );
	printf("tag Ran explicit check on %d/%d particles.\n", nexplicit, np );
	

/*
	for( int bx = 0; bx < id_grain_x; bx++ )
	for( int by = 0; by < id_grain_y; by++ )
	for( int bz = 0; bz < id_grain_z; bz++ )
	{
		int b = (bx * grain_y + by ) * grain_z + bz;

		

		double cen[3] = {  
			(bx+0.5) * PBC_vec[0][0] / grain_x,
			(by+0.5) * PBC_vec[1][1] / grain_y,
			(bz+0.5) * PBC_vec[2][2] / grain_z };

		if( bins[b].np == 0 ) continue;

		double r_max = 0;
		double dr[bins[b].np];

		for( int px = 0; px < bins[b].np; px++ )
		{
			int p = bins[b].list[px];

			double dr[3] = { rp[3*p+0] - cen[0],
					 rp[3*p+1] - cen[1],
					 rp[3*p+2] - cen[2] };

			while( dr[0] < -PBC_vec[0][0]/2 ) dr[0] += PBC_vec[0][0];
			while( dr[1] < -PBC_vec[1][1]/2 ) dr[1] += PBC_vec[1][1];
			while( dr[2] < -PBC_vec[2][2]/2 ) dr[2] += PBC_vec[2][2];
			while( dr[0] >  PBC_vec[0][0]/2 ) dr[0] -= PBC_vec[0][0];
			while( dr[1] >  PBC_vec[1][1]/2 ) dr[1] -= PBC_vec[1][1];
			while( dr[2] >  PBC_vec[2][2]/2 ) dr[2] -= PBC_vec[2][2];

			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r > r_max )
				r_max = r;

			dr[px] = r;
					
		}
		
		int f;
		double u,v;
		printf("FINAL LEVEL Checking %d particles near %lf %lf %lf ", bins[b].np, cen[0], cen[1], cen[2] );
		if( ! theSurface->withinRadius( cen, &f, &u, &v,  M, mlow, mhigh, -1,  r_max + srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
		{
			printf("... all are clear.\n");
			for( int px = 0; px < bins[b].np; px++ )
			{
				int p = bins[b].list[px];
				distance[p] = (r_max - dr[px]) + srdCollisionEstimator;

				rp_dist[3*p+0] = rp[3*p+0];
				rp_dist[3*p+1] = rp[3*p+1];
				rp_dist[3*p+2] = rp[3*p+2];				
			}
		}	
		else
		{
			printf("...checking them all...");
			double nok = 0, nbad=0;
			for( int px = 0; px < bins[b].np; px++ )
			{
				int p = bins[b].list[px];
			
				if( ! theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  srdCollisionEstimator, vertex_data, ptr_to_data, nump ) )
				{
					nok += 1;
					distance[p] = srdCollisionEstimator;

					rp_dist[3*p+0] = rp[3*p+0];
					rp_dist[3*p+1] = rp[3*p+1];
					rp_dist[3*p+2] = rp[3*p+2];				
				}
				else
				{
					nbad += 1;
					theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, -1 );
					rp_dist[3*p+0] = rp[3*p+0];
					rp_dist[3*p+1] = rp[3*p+1];
					rp_dist[3*p+2] = rp[3*p+2];

					double rpt[3], rnrm[3];

					theSurface->evaluateRNRM( f, u, v, rpt, rnrm, r );

					double dr[3] = { rp[3*p+0] - rpt[0],
						 rp[3*p+1] - rpt[1],
						 rp[3*p+2] - rpt[2] };
	
					double len = normalize(dr);
				
					double dp = dr[0] * rnrm[0] + dr[1] * rnrm[1] + dr[2] * rnrm[2];

//			printf("init len: %le dp: %le dist: %le\n", len, dp, distance[p] );

					if( !planar )
					{
						if( dp > 0 )
							inside_outside[p] = POINT_INSIDE;
						else
							inside_outside[p] = POINT_OUTSIDE;
					}
				}
				printf(" %lf%% were clear.\n", 100 * nok / (nok+nbad) );
			}
		}
	}
*/
	free(vertex_data);
	free(ptr_to_data);
	free(nump);
}

void srd_integrator::altTag( double * r, surface * theSurface, double delta_hull_collision, double **M, int mlow, int mhigh )
{
	double LA = PBC_vec[0][0];
	double LB = PBC_vec[1][1];
	double LC = PBC_vec[2][2];
	
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );

	double search_radius = 25;
	for( int p = 0; p < np; p++ )
	{
		// only re-tag particles if their tag is zero.

//		if( fabs(alt_topological_tag[p]) < 1e-14 )
		{
			int f;
			double u, v;
		
			int did_col = theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1, search_radius, vertex_data, ptr_to_data, nump ); 
				
			if( did_col )
			{
				theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, search_radius );
				rp_dist[3*p+0] = rp[3*p+0];
				rp_dist[3*p+1] = rp[3*p+1];
				rp_dist[3*p+2] = rp[3*p+2];

				double rnrm[3], nrm[3];

				theSurface->evaluateRNRM(f, u, v, rnrm, nrm, r ); 

				double dr[3] = { rp[3*p+0] - rnrm[0],
						 rp[3*p+1] - rnrm[1],
						 rp[3*p+2] - rnrm[2] };

				while( dr[0] < -LA/2 ) dr[0] += LA;
				while( dr[0] >  LA/2 ) dr[0] -= LA;
				while( dr[1] < -LB/2 ) dr[1] += LB;
				while( dr[1] >  LB/2 ) dr[1] -= LB;
				while( dr[2] < -LC/2 ) dr[2] += LC;
				while( dr[2] >  LC/2 ) dr[2] -= LC;
				double dp = dr[0] * nrm[0] + dr[1] * nrm[1] + dr[2] * nrm[2];
				
				alt_topological_tag[p] = dp;

//				printf("alt_tag_pair: %lf %lf\n", alt_topological_tag[p], alt_last_known_tag[p] );

				alt_last_known_tag[p] = dp;
			}
			else
				alt_last_known_tag[p] = 0;
		} 
//		else
//			last_known_tag[p] = 0;
	}
	free(vertex_data);
	free(ptr_to_data);
	free(nump);
}

void srd_integrator::checkResolve( double * r, surface * theSurface, double delta_hull_collision, double **M, int mlow, int mhigh )
{
	double LA = PBC_vec[0][0];
	double LB = PBC_vec[1][1];
	double LC = PBC_vec[2][2];
	
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );

	double search_radius = 25;
	for( int p = 0; p < np; p++ )
	{
		// only re-tag particles if their tag is zero.

//		if( fabs(alt_topological_tag[p]) < 1e-14 )
		{
			int f;
			double u, v;
		
			int did_col = theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1, search_radius, vertex_data, ptr_to_data, nump ); 
				
			if( did_col )
			{
				theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, search_radius );
				rp_dist[3*p+0] = rp[3*p+0];
				rp_dist[3*p+1] = rp[3*p+1];
				rp_dist[3*p+2] = rp[3*p+2];

				double rnrm[3], nrm[3];

				theSurface->evaluateRNRM(f, u, v, rnrm, nrm, r ); 

				double dr[3] = { rp[3*p+0] - rnrm[0],
						 rp[3*p+1] - rnrm[1],
						 rp[3*p+2] - rnrm[2] };

				while( dr[0] < -LA/2 ) dr[0] += LA;
				while( dr[0] >  LA/2 ) dr[0] -= LA;
				while( dr[1] < -LB/2 ) dr[1] += LB;
				while( dr[1] >  LB/2 ) dr[1] -= LB;
				while( dr[2] < -LC/2 ) dr[2] += LC;
				while( dr[2] >  LC/2 ) dr[2] -= LC;
				double dp = dr[0] * nrm[0] + dr[1] * nrm[1] + dr[2] * nrm[2];
				
				
				if( dp * alt_topological_tag[p] < 0 )
				{
//					printf("ABSOLUTELY MISSED COLLISION %d tag1 %lf tag2 %lf\n", p, dp, alt_topological_tag[p] );
				}

				alt_last_known_tag[p] = dp;
			}
			else
				alt_last_known_tag[p] = 0;
		} 
	}
	free(vertex_data);
	free(ptr_to_data);
	free(nump);
}

int srd_integrator::resolveCollision( double *r, double *g, double *qdot, SparseMatrix *effm, double delta_hull_collision, surface *theSurface, double **M, int mlow, int mhigh, force_set *theForceSet, double time_step, double force_factor )
{
	double *vertex_data = NULL;
	int *ptr_to_data = (int *)malloc( sizeof(int) * theSurface->nt );
	int *nump = (int *)malloc( sizeof(int) * theSurface->nt );
	theSurface->buildFaceData( &vertex_data, ptr_to_data, nump );

	double LA = PBC_vec[0][0];
	double LB = PBC_vec[1][1];
	double LC = PBC_vec[2][2];
	double search_radius = 20.0;
	//double search_radius = delta_hull_collision;

	double *dp_vector = (double *)malloc( sizeof(double) * theSurface->nv * 3 );
	memset( dp_vector, 0, sizeof(double) * theSurface->nv * 3 );

	double total_dp[3]={0,0,0};

	int n_col = 0;

	for( int p = 0; p < np; p++ )
	{
		if( fabs(topological_tag[p]) < 1e-14 ) continue;
	

		int f;
		double u, v;
	
		int did_col = theSurface->withinRadius( rp+3*p, &f, &u, &v,  M, mlow, mhigh, -1,  delta_hull_collision, vertex_data, ptr_to_data, nump ); 
			
		if( did_col )
		{
			theSurface->nearPointOnBoxedSurface( rp+3*p, &f, &u, &v, M, mlow, mhigh, distance+p, delta_hull_collision );
			rp_dist[3*p+0] = rp[3*p+0];
			rp_dist[3*p+1] = rp[3*p+1];
			rp_dist[3*p+2] = rp[3*p+2];
			double rnrm[3], nrm[3];

			theSurface->evaluateRNRM(f, u, v, rnrm, nrm, r ); 

			double dr[3] = { rp[3*p+0] - rnrm[0],
					 rp[3*p+1] - rnrm[1],
					 rp[3*p+2] - rnrm[2] };

			while( dr[0] < -LA/2 ) dr[0] += LA;
			while( dr[0] >  LA/2 ) dr[0] -= LA;
			while( dr[1] < -LB/2 ) dr[1] += LB;
			while( dr[1] >  LB/2 ) dr[1] -= LB;
			while( dr[2] < -LC/2 ) dr[2] += LC;
			while( dr[2] >  LC/2 ) dr[2] -= LC;
			double dp = dr[0] * nrm[0] + dr[1] * nrm[1] + dr[2] * nrm[2];
			

			
			
			int do_collide = 0;

			if( dp < 0 && topological_tag[p] > 0 )
			{

				do_collide = 1;	

//				printf("SWEEP COLLISION! (-)\n");
			}
			else if( dp > 0 && topological_tag[p] < 0 )
			{
				do_collide = 1;	

//				printf("SWEEP COLLISION! (+)\n");
			}

/*			if( do_collide && did_collide[p] )
			{
				printf("COLLISION, but topology changed.\n");
			}
			else if( !do_collide && did_collide[p] )
			{
			}

			if( do_collide && !did_collide[p] )
*/

#ifdef SPECIFIC_PARTICLE
			if( p == 2267 )				printf("2267 last %lf cur %lf tag %lf %s z %lf\n", last_known_tag[p], dp, topological_tag[p], ( do_collide ? "COLLIDE" : "NOPE" ), rp[3*p+2] );
#endif

			if( do_collide )
			{
				// DOCUMENTATION MARKER: SRD/Mesh collision.

				int col_f;	
				double col_u, col_v;
			
				theSurface->nearPointOnBoxedSurface( rp+3*p, &col_f, &col_u, &col_v, M, mlow, mhigh, distance+p, search_radius );

				double col_w = 1-col_u-col_v;

				double dA = (col_u-1)*(col_u-1) + col_v*col_v+col_w*col_w;
				double dB = (col_v-1)*(col_v-1) + col_u*col_u+col_w*col_w;
				double dC = (col_w-1)*(col_w-1) + col_u*col_u+col_v*col_v;

				double r_collision[3];
				double n_collision[3];
				
				int use_code = 0; // use the main face vertex if dC is smallest.

				if( dA < dC && dA < dB )
					use_code = 1; // du
				else if( dB < dC )
					use_code = 2;

				// get n_collision.
				theSurface->evaluateRNRM(col_f, col_u, col_v, r_collision, n_collision, r ); 

				int tri, near_vertex;
	
				if( col_f < theSurface->nf_faces )
					near_vertex = theSurface->theFormulas[col_f*theSurface->nf_g_q_p].cp[use_code];
				else	
					near_vertex = theSurface->theIrregularFormulas[(col_f-theSurface->nf_faces)*theSurface->nf_irr_pts].cp[use_code];

				double M_ii = effm->diagonal_element[near_vertex]; 
	
				double m_srd_p_ncol = vp[3*p+0] * n_collision[0] +
						      vp[3*p+1] * n_collision[1] +  
						      vp[3*p+2] * n_collision[2];

				double m_mesh_p_ncol = qdot[3*near_vertex+0] * n_collision[0] +
						       qdot[3*near_vertex+1] * n_collision[1] +  
						       qdot[3*near_vertex+2] * n_collision[2]; 

				double alpha_i = 2 * ( m_srd_p_ncol - m_mesh_p_ncol) / ( 1.0/mass + M_ii);
				double alpha_SRD = -alpha_i;

				
				double rp_proj = (vp[3*p+0] * n_collision[0] + vp[3*p+1] * n_collision[1] + vp[3*p+2] * n_collision[2]) * time_step;			
				/*
					vp_proj is the velocity projected along the normal. the solution to
					constant momentum and kinetic energy is

					dvp_nrm = - 2 * membrane_point_mass * vnrm / (membrane_point_mass + mass)
				        dvm_nrm = - 2 * mp * vnrm / ( membrane_point_mass + mass)

				*/
	
				// move the SRD particle ... could compute the real time of collision I guess...

				rp[3*p+0] = r_collision[0] - rp_proj * n_collision[0] / 2;
				rp[3*p+1] = r_collision[1] - rp_proj * n_collision[1] / 2;
				rp[3*p+2] = r_collision[2] - rp_proj * n_collision[2] / 2;

				double dr[3] = { rp[3*p+0] - r_collision[0],
						 rp[3*p+1] - r_collision[1],
						 rp[3*p+2] - r_collision[2] };

				while( dr[0] < -LA/2 ) dr[0] += LA;
				while( dr[0] >  LA/2 ) dr[0] -= LA;
				while( dr[1] < -LB/2 ) dr[1] += LB;
				while( dr[1] >  LB/2 ) dr[1] -= LB;
				while( dr[2] < -LC/2 ) dr[2] += LC;
				while( dr[2] >  LC/2 ) dr[2] -= LC;

				double dpv = dr[0] * nrm[0] + dr[1] * nrm[1] + dr[2] * nrm[2];

				topological_tag[p] = dpv; 
				last_known_tag[p] = topological_tag[p];

				//change in momentum

				double dp[3] = {0,0,0};
						
				dp[0] = alpha_SRD * n_collision[0];
				dp[1] = alpha_SRD * n_collision[1];
				dp[2] = alpha_SRD * n_collision[2];

				vp[3*p+0] += dp[0] / mass;
				vp[3*p+1] += dp[1] / mass;
				vp[3*p+2] += dp[2] / mass;
			
				dp_vector[3*near_vertex+0] = alpha_i * n_collision[0];	
				dp_vector[3*near_vertex+1] = alpha_i * n_collision[1];	
				dp_vector[3*near_vertex+2] = alpha_i * n_collision[2];	
	
				n_col++;
			}
		}
	}


	int nv = theSurface->nv;

	for( int v1 = 0; v1 < nv; v1++ )
	{
		g[3*v1+0] -= force_factor * dp_vector[3*v1+0] / time_step;		
		g[3*v1+1] -= force_factor * dp_vector[3*v1+1] / time_step;	
		g[3*v1+2] -= force_factor * dp_vector[3*v1+2] / time_step;		
	}

#ifdef DEBUG_MOMENTUM_CHANGE
	if( n_col > 0 )
		printf("Collision total solvent momentum change: %le %le %le\n", total_dp[0], total_dp[1], total_dp[2] ); 
#endif
	free(dp_vector);
#ifdef SPECIFIC_PARTICLE
	printf("tag 2267 %lf %lf\n", topological_tag[2267], last_known_tag[2267] );
#endif
	free(vertex_data);
	free(ptr_to_data);
	free(nump);

	return n_col;
}

void srd_integrator::subPos( double *dr )
{
	for( int p = 0; p < np; p++ )
	{
		rp[3*p+0] -= dr[0];
		rp[3*p+1] -= dr[1];
		rp[3*p+2] -= dr[2];
		
		rp_prev[3*p+0] -= dr[0];
		rp_prev[3*p+1] -= dr[1];
		rp_prev[3*p+2] -= dr[2];
	}
}

void srd_integrator::setMembraneMinMax( double *r, surface *theSurface )
{
	for( int c = 0; c < 3; c++ )
	{
		membrane_min[c] = 1e100;
		membrane_max[c] = -1e100;
	}

	for( int v = 0; v < theSurface->nv; v++ )
	for( int c = 0; c < 3; c++ )
	{
		if( r[3*v+c] > membrane_max[c] ) membrane_max[c] = r[3*v+c]; 
		if( r[3*v+c] < membrane_min[c] ) membrane_min[c] = r[3*v+c]; 
	}

	for( int c = 0; c < 3; c++ )
	{
		if( membrane_max[c] - membrane_min[c] > 0.7 * PBC_vec[c][c] )
		{
		// don't try to use this dimension to exclude collisions.
			membrane_max[c] = 1e100;
			membrane_min[c] = -1e100;
		}
	}
}
