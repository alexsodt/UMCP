#include "pcomplex.h"
#include "interp.h"
#include "util.h"
#include "mutil.h"
#include "units.h"
#include "gsl_random_globals.h"
#include <math.h>
#include <string.h>
#include "parallel.h"
#include "p_p.h"
#include <typeinfo>
#include <ctype.h>
#include "globals.h"
#include "rd.h"
#include "M_matrix.h"

//#define DISABLE_POINT_GRADIENT_DEBUG

extern double kc;
static double default_particle_area = 65;
double lipid_DC = 1e10; //Angstrom^2/s
double solution_DC = 1e10; //Angstrom^2/s

// For doing Newtonian/Langevin dynamics:
static double default_mass = 20000;
// For doing brownian dynamics:
static double default_dc   = (1e-6)*(1e8)*(1e8); // cm^2/s to Angstroms^2/s

//#define DISABLE_ATTACH
#define PART_A
#define PART_B
//#define TURN_OFF_TERM2
//#define TURN_OFF_TERM4
//#define DO_FIRST_ORDER
//
void pcomplex::base_init( void )
{
	do_bd = 0;
	debug = DEBUG_OFF;
	nwatchers = 0;
	disabled = 0;
	alive = 1;
}

void pcomplex::copyParentParameters( pcomplex *parent )
{
	do_bd = parent->do_bd;
	is_inside = parent->is_inside;
}

/* general routines */

void pcomplex::alloc( void )
{
	mass = (double *)malloc( sizeof(double) * nsites );
	for( int s = 0; s < nsites; s++ )
		mass[s] = default_mass;
	rall = (double *)malloc( sizeof(double) * 3 * nsites );
	memset( rall, 0, sizeof(double) * 3 * nsites );

	sid = (int*)malloc( sizeof(int) * nsites );
	for( int s = 0; s < nsites; s++ )
		sid[s] = -1;

	rd_timestep_disabled = (int *)malloc( sizeof(int) * nsites );
	memset( rd_timestep_disabled, 0, sizeof(int) * nsites );
	stype = (int*)malloc( sizeof(int) * nsites );

	fs = (int *)malloc( sizeof(int) * nsites );
	puv = (double *)malloc( sizeof(double) * 2 * nsites );

	p = (double *)malloc( sizeof(double) * nsites * 3 );
	qdot = (double *)malloc( sizeof(double) * nsites * 3 );
	memset( p, 0, sizeof(double) * 3 * nsites );
	memset( qdot, 0, sizeof(double) * 3 * nsites );

	DC = (double *)malloc( sizeof(double) * nsites );
	for( int p = 0; p < nsites; p++ )		
		DC[p] = lipid_DC;


	grad_fs = (int *)malloc( sizeof(int) * nsites );

	grad_puv = (double *)malloc( sizeof(double) * 2 * nsites );
	save_grad = (double *)malloc( sizeof(double) * 3 * nsites );
	memset( save_grad, 0, sizeof(double) * 3 * nsites );

	/* caches */	
	cache_grad = (double *)malloc( sizeof(double) * 3 * nsites );
	cache_rall = (double *)malloc( sizeof(double) * 3 * nsites );
	cache_f = (int *)malloc( sizeof(int) * nsites );
	cache_puv = (double *)malloc( sizeof(double) * 2 * nsites);

	sigma = (double *)malloc( sizeof(double) * nsites );
	memset( sigma, 0, sizeof(double) * nsites );

	att_eps = (double *)malloc( sizeof(double) * nsites );
	memset( att_eps, 0, sizeof(double) * nsites );
	
	att_sigma = (double *)malloc( sizeof(double) * nsites );
	memset( att_sigma, 0, sizeof(double) * nsites );
	
	p_area = (double *)malloc( sizeof(double) * nsites );
	p_c0 = ( double *)malloc( sizeof(double) * nsites );

	for( int s = 0; s < nsites; s++ )
	{
		p_area[s] = default_particle_area;
		p_c0[s] = 0;
	}

	coord_transform = (double *)malloc( sizeof(double) * 4 * nattach );

	PBC_ext = (double *)malloc( sizeof(double) * 3 * nsites );
	memset( PBC_ext, 0, sizeof(double) * 3 * nsites );
	last_pos = (double *)malloc( sizeof(double) * 3 * nsites );
	memset( last_pos, 0, sizeof(double) * 3 * nsites );

}

void pcomplex::print_type( char **outp )
{
	// override this if you want a better name printed.
	int len = strlen(typeid(*this).name());
	char *temp = (char *)malloc( sizeof(char) * (len+1) );

	sprintf(temp, "%s", typeid(*this).name());

	char *p = temp;
	while( *p && !isalpha(*p) ) p += 1;

	*outp = (char *)malloc( sizeof(char) * (len+1) );
	sprintf(*outp, "%s", p );
	
	free(temp);
}

void pcomplex::destroy( void )
{
}

void pcomplex::activateBrownianDynamics( void )
{
	do_bd = 1;
}



void pcomplex::move_inside( void )
{
	is_inside = 1;
}

void pcomplex::move_outside( void )
{
	is_inside = 0;
}

/*

Propagate p, qdot, q, etc.

*/

void pcomplex::debug_dPinv_dq( surface * theSurface, double *rsurf  )
{
	int s = 0;	
			
	double coefs[3*(MAX_VALENCE+6)];
	int coef_list[MAX_VALENCE+6];
	int ncoef=0;
	theSurface->get_pt_dcoeffs( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], coefs, coef_list, &ncoef ); 
	double *dcoef = coefs + ncoef;
	
	double *dPinv_dq = (double *)malloc( sizeof(double) * 4 * ncoef );

	int f = grad_fs[s];
	double u = grad_puv[2*s+0];
	double v = grad_puv[2*s+1];
			
	double drdu[3];
	double drdv[3];

	theSurface->ru( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1],
			rsurf, drdu );
	theSurface->rv( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1],
			rsurf, drdv );
	
	double P_uv[4] = {
		0,	0,
		0,	0
	};
	
	theSurface->fetchPuv( f, u, v, P_uv, rsurf );
	
	P_uv[0] *= mass[s];
	P_uv[1] *= mass[s];
	P_uv[2] *= mass[s];
	P_uv[3] *= mass[s];

	double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
	double Pinv0[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };

	double EPS_VAL = 0.1;
	
	double der_Pinv[4*ncoef*3];

	memset( der_Pinv, 0, sizeof(double) * 3 * ncoef * 4 );
	
	for( int c = 0; c < ncoef; c++ )
	for( int u = 0; u < 2; u++ )
	for( int v = 0; v < 2; v++ )
	{
		double der_P_uu_c_x =  mass[s] *( drdu[0] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[0]);
		double der_P_uu_c_y =  mass[s] *( drdu[1] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[1]);
		double der_P_uu_c_z =  mass[s] *( drdu[2] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[2]);
		
		double der_P_uv_c_x =  mass[s] *( drdu[0] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[0]);
		double der_P_uv_c_y =  mass[s] *( drdu[1] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[1]);
		double der_P_uv_c_z =  mass[s] *( drdu[2] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[2]);
		
		double der_P_vu_c_x =  mass[s] *( drdu[0] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[0]);
		double der_P_vu_c_y =  mass[s] *( drdu[1] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[1]);
		double der_P_vu_c_z =  mass[s] *( drdu[2] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[2]);
	
		double der_P_vv_c_x =  mass[s] *( drdv[0] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[0]);
		double der_P_vv_c_y =  mass[s] *( drdv[1] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[1]);
		double der_P_vv_c_z =  mass[s] *( drdv[2] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[2]);

		der_Pinv[c*4*3+u*2+v] -= Pinv0[u*2+0] * der_P_uu_c_x * Pinv0[0*2+v];
		der_Pinv[c*4*3+u*2+v] -= Pinv0[u*2+0] * der_P_uv_c_x * Pinv0[1*2+v];
		der_Pinv[c*4*3+u*2+v] -= Pinv0[u*2+1] * der_P_vu_c_x * Pinv0[0*2+v];
		der_Pinv[c*4*3+u*2+v] -= Pinv0[u*2+1] * der_P_vv_c_x * Pinv0[1*2+v];
	}

	for( int c = 0; c < ncoef; c++ )
	{
		printf("COEF %d\n", c );
		for( int ieps = 0; ieps < 20; ieps++ )
		{
			rsurf[3*coef_list[c]+0] += ieps * EPS_VAL;
			
			double P_uv[4] = {
				0,	0,
				0,	0
			};
	
			theSurface->fetchPuv( f, u, v, P_uv, rsurf );
	
			P_uv[0] *= mass[s];
			P_uv[1] *= mass[s];
			P_uv[2] *= mass[s];
			P_uv[3] *= mass[s];

			double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
			double Pinv[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };

			printf("INST PUV: %.14le %.14le %.14le %.14le EXTRAP %.14le %.14le %.14le %.14le\n",
				Pinv[0], Pinv[1], Pinv[2], Pinv[3], 
				Pinv0[0] + der_Pinv[c*4*3+0] * ieps*EPS_VAL,
				Pinv0[1] + der_Pinv[c*4*3+1] * ieps*EPS_VAL,
				Pinv0[2] + der_Pinv[c*4*3+2] * ieps*EPS_VAL,
				Pinv0[3] + der_Pinv[c*4*3+3] * ieps*EPS_VAL
		 );

			rsurf[3*coef_list[c]+0] -= ieps * EPS_VAL;
		}
	}
}


void pcomplex::prepareForGradient( void )
{
//	for( int s = 0; s < nsites; s++ )
//	{
//		save_grad[2*s+0] = 0;
//		save_grad[2*s+1] = 0;
//	}
	memset( save_grad, 0, sizeof(double) * 3 * nsites );
}

double pcomplex::update_dH_dq( Simulation *theSimulation, double time_step, double timestep_total)
{
	if( disabled ) return time_step;

	double *alphas = theSimulation->alpha;

	// select the minimum timestep we can do.

	double use_dt = time_step;
	
	int is_irreg = 0;
	for( int s = 0; s < nattach; s++ )
	{
		surface *theSurface;
		double *rsurf;
		theSimulation->fetch(sid[s],&theSurface,&rsurf);

		if( grad_fs[s] >= theSurface->nf_faces )
			is_irreg = 1;
	}

	// frac mult is used for distributing the mesh gradients throughout propagation..	
	double frac_mult=1.0;

	if( is_irreg && time_step > 0 && ! do_bd )
	{
		double fraction = 1.0;

		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = NULL;
			surface *theSurface =NULL;
			for( surface_record *tRec = theSimulation->allSurfaces; tRec; tRec = tRec->next )
			{
				if( tRec->id == sid[s] )
				{
					sRec = tRec;	
					theSurface = sRec->theSurface;
					break;
				}
			}

			double *rsurf = sRec->r;
			double *mesh_grad = sRec->g;
			double *mesh_qdot = sRec->qdot;
			double *mesh_qdot0 = sRec->qdot0;

			int f = grad_fs[s];
	
			double u = grad_puv[2*s+0];
			double v = grad_puv[2*s+1];

#ifdef USE_KIN_E	
			double P_uv[4] = {
				0,	0,
				0,	0
			};
	
			theSurface->fetchPuv( f, u, v, P_uv, rsurf );
		
			P_uv[0] *= mass[s];
			P_uv[1] *= mass[s];
			P_uv[2] *= mass[s];
			P_uv[3] *= mass[s];
	
			double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
			double Pinv[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };
	
			double qdot0[2] = { Pinv[0] * p[2*s+0] + Pinv[1] * p[2*s+1],
					    Pinv[2] * p[2*s+0] + Pinv[3] * p[2*s+1] };
			double *d_P_duv;
	
			int nc = theSurface->fetchdP_duv( f, u, v, &d_P_duv, rsurf );
			for( int t = 0; t < 8; t++ )
				d_P_duv[t] *= mass[s];
		
			// delta-p =  changes by dt qdot (dP) qdot 
	
			//save_grad[2*s+0] -=  0.5 * qdot0[0] * d_P_duv[0*2+0] * qdot0[0]; 
		
			// expected dKE:
	
			double inter1[2]={0,0};
			double inter2[2]={0,0};

			// du.

			inter1[0] += d_P_duv[0] * qdot0[0];
			inter1[0] += d_P_duv[1] * qdot0[1];
			
			inter1[1] += d_P_duv[2] * qdot0[0];
			inter1[1] += d_P_duv[3] * qdot0[1];	

			double dPu_dt = inter1[0] * qdot0[0] + inter1[1] * qdot0[1];

			inter2[0] += d_P_duv[0] * qdot0[0];
			inter2[0] += d_P_duv[1] * qdot0[1];
			
			inter2[1] += d_P_duv[2] * qdot0[0];
			inter2[1] += d_P_duv[3] * qdot0[1];	
			double dPv_dt = inter2[0] * qdot0[0] + inter2[1] * qdot0[1];
	
			double dKE = 0.25 * (dPu_dt * qdot0[0] + dPu_dt * qdot0[1])*time_step*AKMA_TIME;
			double tolerance = 0.1;
			double fr = tolerance * fabs(0.592/dKE); // just use 300K for now.

			if( fr < fraction )
			{
	//			printf("f u v %d %le %le dKE_exp %.14le fraction %le\n", f, u, v, dKE, fr ); 
				fraction = fr;
			}

			free(d_P_duv);
#else

			double P_uv[4] = {
				0,	0,
				0,	0
			};
	
			theSurface->fetchPuv( f, u, v, P_uv, rsurf );
		
			P_uv[0] *= mass[s];
			P_uv[1] *= mass[s];
			P_uv[2] *= mass[s];
			P_uv[3] *= mass[s];
	
			double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
			double Pinv[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };
	
			double qdot0[2] = { Pinv[0] * p[2*s+0] + Pinv[1] * p[2*s+1],
					    Pinv[2] * p[2*s+0] + Pinv[3] * p[2*s+1] };

			// from propagating qdot:
			double expected_dq[2] = { AKMA_TIME * time_step * qdot0[0], AKMA_TIME * time_step * qdot[1] };

			

			// d(Minv)_{uv}_u

			double *d_P_duv;
	
			int nc = theSurface->fetchdP_duv( f, u, v, &d_P_duv, rsurf );
			for( int t = 0; t < 8; t++ )
				d_P_duv[t] *= mass[s];
	

			double ax[2] = { d_P_duv[0], d_P_duv[4] };
			double bx[2] = { d_P_duv[1], d_P_duv[5] };
			double cx[2] = { d_P_duv[2], d_P_duv[6] };
			double dx[2] = { d_P_duv[3], d_P_duv[7] };

			double a=P_uv[0];
			double b=P_uv[1];
			double c=P_uv[2];
			double d=P_uv[3];

			double dM0_det[2] = { ax[0]*d-bx[0]*c-cx[0]*b+dx[0]*a,
					      ax[1]*d-bx[1]*c-cx[1]*b+dx[1]*a };
			dM0_det[0] /= (b*c-a*d);
			dM0_det[1] /= (b*c-a*d);		
			double fr = fabs(dM0_det[0] * expected_dq[0] + dM0_det[1] * expected_dq[1]);
			double tolerance = 0.005;

			if( tolerance/fr < fraction )
			{
//				printf("f u v %d %le %le d_metric: %le setting fraction\n", f, u, v, fr, tolerance/fr ); 
				fraction = tolerance/fr;
			}

			free(d_P_duv);
#endif
		}

		time_step *= fraction;
		frac_mult = time_step / timestep_total;
//		printf("Fraction: %le Frac_mult: %le\n", fraction, frac_mult );
	}
	

	

	if( nattach > 0  )
	{	

		double surfacer_g[3*nattach];
		double surfacen_g[3*nattach];

		memset( surfacer_g, 0, sizeof(double)*3*nattach);
		memset( surfacen_g, 0, sizeof(double)*3*nattach);

		double pg[2*nattach];
		memset( pg, 0, sizeof(double) * 2 * nattach );
		grad( theSimulation, surfacer_g, surfacen_g );
		if( debug == DEBUG_OFF || debug == DEBUG_NO_T )
			AttachG( theSimulation, pg ); 
	
		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = theSimulation->fetch(sid[s]);
			surface *theSurface = sRec->theSurface;

			int nv = theSurface->nv;
			double *rsurf = sRec->r;
			double *mesh_grad = sRec->g;
			double *mesh_qdot = sRec->qdot;
			double *mesh_qdot0 = sRec->qdot0;

			double coefs[3*(MAX_VALENCE+6)];
			int coef_list[MAX_VALENCE+6];
			int ncoef=0;
			theSurface->get_pt_dcoeffs( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], coefs, coef_list, &ncoef ); 
			double *dcoef = coefs + ncoef;

			// the attachment gradient.
			double uv_grad[2] = { pg[2*s+0],pg[2*s+1] };

#ifndef DISABLE_POINT_GRADIENT_DEBUG
			if( debug == DEBUG_OFF || debug == DEBUG_NO_T )
			theSurface->pointGradient( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1],
				rsurf, mesh_grad, uv_grad, surfacer_g+3*s, surfacen_g+3*s, frac_mult ); 
#endif
			// the contribution from dV/du and dV/dv
			save_grad[2*s+0] += uv_grad[0];
			save_grad[2*s+1] += uv_grad[1];

			if( debug == DEBUG_NO_T )
				continue;

			double drdu[3]={0,0,0};
			double drdv[3]={0,0,0};

			theSurface->ru( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1],
					rsurf, drdu );
			theSurface->rv( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1],
					rsurf, drdv );

			double ruu[3] = {0,0,0};
			double ruv[3] = {0,0,0};	
			double rvv[3] = {0,0,0};

			theSurface->r2der( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], rsurf, ruu, ruv, rvv );

		// gradient is uv_grad.

			// compute instantaneous inverse mass matrix and its u/v/q derivatives that show up as effective forces.
			int f = grad_fs[s];
			double u = grad_puv[2*s+0];
			double v = grad_puv[2*s+1];
	
			double P_uv[4] = {
				0,	0,
				0,	0
			};

		
			theSurface->fetchPuv( f, u, v, P_uv, rsurf );
		
			P_uv[0] *= mass[s];
			P_uv[1] *= mass[s];
			P_uv[2] *= mass[s];
			P_uv[3] *= mass[s];


			double *C_iu;
		
			int nc = theSurface->fetchCiu( f, u, v, &C_iu, rsurf );
	
			for( int t = 0; t < 3 * ncoef * 2; t++ )
				C_iu[t] *= mass[s];	

				
			double *d_P_duv;
				
			// these derivatives include the mesh derivatives wrt u/v
			nc = theSurface->fetchdP_duv( f, u, v, &d_P_duv, rsurf );
			for( int t = 0; t < 8; t++ )
				d_P_duv[t] *= mass[s];

			for( int c = 0; c < ncoef; c++ )
			for( int uv = 0; uv < 4; uv++ )
			for( int xyz = 0; xyz < 3; xyz++ )
				d_P_duv[8+c*4*3+uv*3+xyz] *= mass[s];
 
			// SCALE d_P_duv BY THE MASS?
		

			// compute d inverse d u and dinverse dv.

			double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
			double Pinv[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };

			double qdot0[2] = { Pinv[0] * p[2*s+0] + Pinv[1] * p[2*s+1],
					    Pinv[2] * p[2*s+0] + Pinv[3] * p[2*s+1] };

			double d_inv_P_uv_du[4], T[4], d_inv_P_uv_dv[4];	

			Mul22( Pinv, d_P_duv, T );
			Mul22( T, Pinv, d_inv_P_uv_du );  	
			Mul22( Pinv, d_P_duv+4, T );
			Mul22( T, Pinv, d_inv_P_uv_dv );  	

			// Pinv and Mij are the block inverses of the particle and mesh, respectively.
			// Ciu is the coupling element.
			//
			// Ciu is stored as [ ncoords * (6) + 2 (u/v dot) x 3 (cartesian) 
			//
			// tM[i,v] = Ciu[i u] P[u,v]

			double *C_Pinv = (double *)malloc( sizeof(double) * 3 * ncoef * 2 );
			memset( C_Pinv, 0, sizeof(double) * 3 * ncoef * 2 );

			// right multiply by Pinv.
			for( int c = 0; c < ncoef; c++ )
			for( int i = 0; i < 2; i++ )
			{
				for( int j = 0; j < 2; j++ )
				{
					C_Pinv[c*3*2+i*3+0] += C_iu[c*3*2+j*3+0] * Pinv[i*2+j]; 
					C_Pinv[c*3*2+i*3+1] += C_iu[c*3*2+j*3+1] * Pinv[i*2+j]; 
					C_Pinv[c*3*2+i*3+2] += C_iu[c*3*2+j*3+2] * Pinv[i*2+j]; 
				}
			}
			
			
			// COMPUTE qdot
	
			// Pcouple is now the piece of the matrix inverse coupling the mesh and particle at first order.	

			// to propagate p we need u/v/q derivatives of the kinetic energy.
	
			// the contribution from dT/du and dT/dv, where T is the mesh kinetic energy.
			
			// d( 1/2 pi Mij^{-1} pj )	

			// this is the zero-order piece. It will diverge near an irregular vertex.

			// TERM 1A, dmu
			save_grad[2*s+0] -=  0.5 * qdot0[0] * d_P_duv[0*2+0] * qdot0[0]; 
			save_grad[2*s+0] -=  0.5 * qdot0[0] * d_P_duv[0*2+1] * qdot0[1]; 
			save_grad[2*s+0] -=  0.5 * qdot0[1] * d_P_duv[1*2+0] * qdot0[0]; 
			save_grad[2*s+0] -=  0.5 * qdot0[1] * d_P_duv[1*2+1] * qdot0[1]; 
			
			// TERM 1A, dnu
			save_grad[2*s+1] -=  0.5 * qdot0[0] * d_P_duv[4+0*2+0] * qdot0[0]; 
			save_grad[2*s+1] -=  0.5 * qdot0[0] * d_P_duv[4+0*2+1] * qdot0[1]; 
			save_grad[2*s+1] -=  0.5 * qdot0[1] * d_P_duv[4+1*2+0] * qdot0[0]; 
			save_grad[2*s+1] -=  0.5 * qdot0[1] * d_P_duv[4+1*2+1] * qdot0[1]; 
			
			// derivative of the particle-particle kinetic energy wrt the mesh coordinates (kinetic energy of bare particle not underlying mesh)

			// TERM 1B, d q
			for( int c1 = 0; c1 < ncoef; c1++ )
			for( int uv2 = 0; uv2 < 2; uv2++ )
			{	// factor of two cancels one half p dot qdot
				mesh_grad[3*coef_list[c1]+0] -=  (mass[s] * qdot0[0] * drdu[0] * dcoef[(uv2)*ncoef+c1] * qdot0[uv2])*frac_mult;
				mesh_grad[3*coef_list[c1]+1] -=  (mass[s] * qdot0[0] * drdu[1] * dcoef[(uv2)*ncoef+c1] * qdot0[uv2])*frac_mult;
				mesh_grad[3*coef_list[c1]+2] -=  (mass[s] * qdot0[0] * drdu[2] * dcoef[(uv2)*ncoef+c1] * qdot0[uv2])*frac_mult;
                                                                                                                 
				mesh_grad[3*coef_list[c1]+0] -=  (mass[s] * qdot0[1] * drdv[0] * dcoef[(uv2)*ncoef+c1] * qdot0[uv2])*frac_mult;
				mesh_grad[3*coef_list[c1]+1] -=  (mass[s] * qdot0[1] * drdv[1] * dcoef[(uv2)*ncoef+c1] * qdot0[uv2])*frac_mult;
				mesh_grad[3*coef_list[c1]+2] -=  (mass[s] * qdot0[1] * drdv[2] * dcoef[(uv2)*ncoef+c1] * qdot0[uv2])*frac_mult;
			}
			// now the first order piece. d ( - Mij^{-1} C Mij^{-1} ) / du		

			// first, the derivative of the zero-order inverse matrix. only the particle flavor changes with u/v

			// can't use qdots on this one because I am taking a derivative of an essential piece of qdot.

#ifdef DO_FIRST_ORDER
			double cross_factor = 2.0;

#ifndef TURN_OFF_TERM2
			// TERM 2A
			// See compute_qdot: TERM 2: this is - Pinv_{mu nu} C{nu i} Pinv_{i j} qdot j 	
			for( int uv1 = 0; uv1 < 2; uv1++ )
			for( int uv2 = 0; uv2 < 2; uv2++ )
			{
				for( int c = 0; c < ncoef; c++ )
				{
					// (-) signs cancel
					save_grad[2*s+0] += cross_factor * 0.5*qdot0[uv1] * d_P_duv[uv1*2+uv2] * C_Pinv[c*3*2+uv2*3+0] * mesh_qdot0[3*coef_list[c]+0]; 
					save_grad[2*s+0] += cross_factor * 0.5*qdot0[uv1] * d_P_duv[uv1*2+uv2] * C_Pinv[c*3*2+uv2*3+1] * mesh_qdot0[3*coef_list[c]+1]; 
					save_grad[2*s+0] += cross_factor * 0.5*qdot0[uv1] * d_P_duv[uv1*2+uv2] * C_Pinv[c*3*2+uv2*3+2] * mesh_qdot0[3*coef_list[c]+2]; 
					
					save_grad[2*s+1] += cross_factor * 0.5*qdot0[uv1] * d_P_duv[4+uv1*2+uv2] * C_Pinv[c*3*2+uv2*3+0] * mesh_qdot0[3*coef_list[c]+0]; 
					save_grad[2*s+1] += cross_factor * 0.5*qdot0[uv1] * d_P_duv[4+uv1*2+uv2] * C_Pinv[c*3*2+uv2*3+1] * mesh_qdot0[3*coef_list[c]+1]; 
					save_grad[2*s+1] += cross_factor * 0.5*qdot0[uv1] * d_P_duv[4+uv1*2+uv2] * C_Pinv[c*3*2+uv2*3+2] * mesh_qdot0[3*coef_list[c]+2]; 
				}
			}
			// second, the derivative of the mesh-particle coupling element itself.
				
			// TERM 2B, derivative of inverse term wrt q

			double f1 = 1;
			double f2 = 1;

			for( int c    = 0; c   < ncoef;  c++ )
			for( int c2    = 0; c2   < ncoef;  c2++ )
			for( int xyz = 0; xyz < 3; xyz++ )
			{		// (-) signs cancel

				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[xyz] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[xyz]) * C_Pinv[c2*3*2+0*3+0] * mesh_qdot0[coef_list[c2]*3+0])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[xyz] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[xyz]) * C_Pinv[c2*3*2+0*3+1] * mesh_qdot0[coef_list[c2]*3+1])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[xyz] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[xyz]) * C_Pinv[c2*3*2+0*3+2] * mesh_qdot0[coef_list[c2]*3+2])*frac_mult; 
				                                                                                                                                         
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[xyz] * dcoef[0*ncoef+c] + dcoef[1*ncoef+c] * drdu[xyz]) * C_Pinv[c2*3*2+0*3+0] * mesh_qdot0[coef_list[c2]*3+0])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[xyz] * dcoef[0*ncoef+c] + dcoef[1*ncoef+c] * drdu[xyz]) * C_Pinv[c2*3*2+0*3+1] * mesh_qdot0[coef_list[c2]*3+1])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[xyz] * dcoef[0*ncoef+c] + dcoef[1*ncoef+c] * drdu[xyz]) * C_Pinv[c2*3*2+0*3+2] * mesh_qdot0[coef_list[c2]*3+2])*frac_mult; 
				                                                                                                                                         
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[xyz] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[xyz]) * C_Pinv[c2*3*2+1*3+0] * mesh_qdot0[coef_list[c2]*3+0])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[xyz] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[xyz]) * C_Pinv[c2*3*2+1*3+1] * mesh_qdot0[coef_list[c2]*3+1])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[xyz] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[xyz]) * C_Pinv[c2*3*2+1*3+2] * mesh_qdot0[coef_list[c2]*3+2])*frac_mult; 
				                                                                                                                                         
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[xyz] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[xyz]) * C_Pinv[c2*3*2+1*3+0] * mesh_qdot0[coef_list[c2]*3+0])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[xyz] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[xyz]) * C_Pinv[c2*3*2+1*3+1] * mesh_qdot0[coef_list[c2]*3+1])*frac_mult; 
				mesh_grad[3*coef_list[c]+xyz] += (f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[xyz] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[xyz]) * C_Pinv[c2*3*2+1*3+2] * mesh_qdot0[coef_list[c2]*3+2])*frac_mult; 

/*				mesh_grad[3*coef_list[c]+1] += f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[1] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[1]) * C_Pinv[c2*3*2+0*3+1] * mesh_qdot0[coef_list[c2]*3+1]; 
				mesh_grad[3*coef_list[c]+2] += f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[2] * dcoef[0*ncoef+c] + dcoef[0*ncoef+c] * drdu[2]) * C_Pinv[c2*3*2+0*3+2] * mesh_qdot0[coef_list[c2]*3+2]; 
				                                                                                                          
				mesh_grad[3*coef_list[c]+0] += f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[0] * dcoef[0*ncoef+c] + dcoef[1*ncoef+c] * drdu[0]) * C_Pinv[c2*3*2+0*3+0] * mesh_qdot0[coef_list[c2]*3+0]; 
				mesh_grad[3*coef_list[c]+1] += f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[1] * dcoef[0*ncoef+c] + dcoef[1*ncoef+c] * drdu[1]) * C_Pinv[c2*3*2+0*3+1] * mesh_qdot0[coef_list[c2]*3+1]; 
				mesh_grad[3*coef_list[c]+2] += f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[2] * dcoef[0*ncoef+c] + dcoef[1*ncoef+c] * drdu[2]) * C_Pinv[c2*3*2+0*3+2] * mesh_qdot0[coef_list[c2]*3+2]; 
				                                                                                                          
				mesh_grad[3*coef_list[c]+0] += f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[0] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[0]) * C_Pinv[c2*3*2+1*3+0] * mesh_qdot0[coef_list[c2]*3+0]; 
				mesh_grad[3*coef_list[c]+1] += f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[1] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[1]) * C_Pinv[c2*3*2+1*3+1] * mesh_qdot0[coef_list[c2]*3+1]; 
				mesh_grad[3*coef_list[c]+2] += f1*cross_factor * 0.5 * qdot0[0] * mass[s] *( drdu[2] * dcoef[1*ncoef+c] + dcoef[0*ncoef+c] * drdv[2]) * C_Pinv[c2*3*2+1*3+2] * mesh_qdot0[coef_list[c2]*3+2]; 
				                                                                                                          
				mesh_grad[3*coef_list[c]+0] += f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[0] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[0]) * C_Pinv[c2*3*2+1*3+0] * mesh_qdot0[coef_list[c2]*3+0]; 
				mesh_grad[3*coef_list[c]+1] += f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[1] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[1]) * C_Pinv[c2*3*2+1*3+1] * mesh_qdot0[coef_list[c2]*3+1]; 
				mesh_grad[3*coef_list[c]+2] += f1*cross_factor * 0.5 * qdot0[1] * mass[s] *( drdv[2] * dcoef[1*ncoef+c] + dcoef[1*ncoef+c] * drdv[2]) * C_Pinv[c2*3*2+1*3+2] * mesh_qdot0[coef_list[c2]*3+2]; 
*/
			}

			// TERM 2Ci, derivative of 0.5 * p_u udot
	
			// udot is modified at first order by the mesh kinetic energy. C is the coupling element.
			// this is the derivative with respect to drdu[0], wrt u
			for( int ci = 0; ci < ncoef; ci++ )
			{
				// p[0] += -( Pinv_{uv} C_{iu} Pinv_{ij} qdot_j
				
				save_grad[2*s+0] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+0] * coefs[ci] * ruu[0] * qdot0[0]; 
				save_grad[2*s+0] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+1] * coefs[ci] * ruu[1] * qdot0[0]; 
				save_grad[2*s+0] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+2] * coefs[ci] * ruu[2] * qdot0[0]; 
				
				save_grad[2*s+1] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+0] * coefs[ci] * ruv[0] * qdot0[0]; 
				save_grad[2*s+1] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+1] * coefs[ci] * ruv[1] * qdot0[0]; 
				save_grad[2*s+1] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+2] * coefs[ci] * ruv[2] * qdot0[0]; 
				
				save_grad[2*s+0] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+0] * coefs[ci] * ruv[0] * qdot0[1]; 
				save_grad[2*s+0] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+1] * coefs[ci] * ruv[1] * qdot0[1]; 
				save_grad[2*s+0] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+2] * coefs[ci] * ruv[2] * qdot0[1]; 
				
				save_grad[2*s+1] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+0] * coefs[ci] * rvv[0] * qdot0[1]; 
				save_grad[2*s+1] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+1] * coefs[ci] * rvv[1] * qdot0[1]; 
				save_grad[2*s+1] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+2] * coefs[ci] * rvv[2] * qdot0[1]; 
			}
			
			// TERM 2Cii
			// now, d/du -> dr/dqi  : requires coefficients of drdu/v.
			for( int ci = 0; ci < ncoef; ci++ )
			for( int uv = 0; uv < 2; uv++ )
			{
				save_grad[2*s+uv] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+0] * dcoef[(uv)*ncoef+ci] * drdu[0] * qdot0[0]; 
				save_grad[2*s+uv] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+1] * dcoef[(uv)*ncoef+ci] * drdu[1] * qdot0[0]; 
				save_grad[2*s+uv] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+2] * dcoef[(uv)*ncoef+ci] * drdu[2] * qdot0[0]; 
				                          
				save_grad[2*s+uv] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+0] * dcoef[(uv)*ncoef+ci] * drdv[0] * qdot0[1]; 
				save_grad[2*s+uv] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+1] * dcoef[(uv)*ncoef+ci] * drdv[1] * qdot0[1]; 
				save_grad[2*s+uv] -= cross_factor * 0.5 * mass[s] * mesh_qdot0[coef_list[ci]*3+2] * dcoef[(uv)*ncoef+ci] * drdv[2] * qdot0[1]; 
			}
	
			// TERM 2D
			// derivative of coupling element wrt control point. from KE, which has dotted pairs velocities, only one mesh point is left in the coupling matrix.

			for( int ck = 0; ck < ncoef; ck++ )
			for( int ci = 0; ci < ncoef; ci++ )
			for( int uv = 0; uv < 2; uv++ )
			{
				// DOES THIS NEED A FACTOR OF 1/2 ??? 
				mesh_grad[3*coef_list[ck]+0] -= (f2*0.5 * cross_factor * mass[s] * mesh_qdot0[coef_list[ci]*3+0] * coefs[ci] * dcoef[(uv)*ncoef+ck] * qdot0[uv])*frac_mult;	
				mesh_grad[3*coef_list[ck]+1] -= (f2*0.5 * cross_factor * mass[s] * mesh_qdot0[coef_list[ci]*3+1] * coefs[ci] * dcoef[(uv)*ncoef+ck] * qdot0[uv])*frac_mult;	
				mesh_grad[3*coef_list[ck]+2] -= (f2*0.5 * cross_factor * mass[s] * mesh_qdot0[coef_list[ci]*3+2] * coefs[ci] * dcoef[(uv)*ncoef+ck] * qdot0[uv])*frac_mult;	
			}

#endif
			// TERM 4
			// derivative of the mesh-mesh coupling wrt u.


#ifndef TURN_OFF_TERM4
			for( int uv = 0; uv < 2; uv++ )
			for( int ci = 0; ci < ncoef; ci++ )
			for( int cj = 0; cj < ncoef; cj++ )
			{	// factor of two is included.
				save_grad[2*s+uv] -= mass[s] * mesh_qdot0[coef_list[ci]*3+0] * mesh_qdot0[coef_list[cj]*3+0] * dcoef[(uv)*ncoef+ci] * coefs[cj];
				save_grad[2*s+uv] -= mass[s] * mesh_qdot0[coef_list[ci]*3+1] * mesh_qdot0[coef_list[cj]*3+1] * dcoef[(uv)*ncoef+ci] * coefs[cj];
				save_grad[2*s+uv] -= mass[s] * mesh_qdot0[coef_list[ci]*3+2] * mesh_qdot0[coef_list[cj]*3+2] * dcoef[(uv)*ncoef+ci] * coefs[cj];
			}
#endif

#endif
			
			// that's it for derivatives wrt u/v


			free(C_Pinv);
			free(d_P_duv);
			free(C_iu);
		}
	}

	return time_step;
}

void pcomplex::compute_qdot( Simulation *theSimulation,  double frac_mult )
{
	double *alphas = theSimulation->alpha;

	if( nattach > 0  )
	{	
	
		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = NULL;
			surface *theSurface =NULL;
			for( surface_record *tRec = theSimulation->allSurfaces; tRec; tRec = tRec->next )
			{
				if( tRec->id == sid[s] )
				{
					sRec = tRec;	
					theSurface = sRec->theSurface;
					break;
				}
			}

			int nv = theSurface->nv;
			double *rsurf = sRec->r;
			double *mesh_grad = sRec->g;
			double *mesh_qdot = sRec->qdot;
			double *mesh_qdot0 = sRec->qdot0;

			double zero_order[2] = { 0,0 };
			double first_order[2] = {0,0};

			double coefs[MAX_VALENCE+6];
			int coef_list[MAX_VALENCE+6];
			int ncoef=0;
			theSurface->get_pt_coeffs( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], coefs, coef_list, &ncoef ); 

			double drdu[3];
			double drdv[3];

			theSurface->ru( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1],
					rsurf, drdu );
			theSurface->rv( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1],
					rsurf, drdv );

			// gradient is uv_grad.

			// compute instantaneous inverse mass matrix and its u/v/q derivatives that show up as effective forces.
			int f = grad_fs[s];
			double u = grad_puv[2*s+0];
			double v = grad_puv[2*s+1];
	
			double P_uv[4] = {
				0,	0,
				0,	0
			};
	
			theSurface->fetchPuv( f, u, v, P_uv, rsurf );
		
			P_uv[0] *= mass[s];
			P_uv[1] *= mass[s];
			P_uv[2] *= mass[s];
			P_uv[3] *= mass[s];

			double *C_iu;
		
			int nc = theSurface->fetchCiu( f, u, v, &C_iu, rsurf );
	
			for( int t = 0; t < 3 * ncoef * 2; t++ )
				C_iu[t] *= mass[s];	

			// compute d inverse d u and dinverse dv.

			double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
			double Pinv[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };
			
			double qdot0[2] = { Pinv[0] * p[2*s+0] + Pinv[1] * p[2*s+1],
					    Pinv[2] * p[2*s+0] + Pinv[3] * p[2*s+1] };

			// Pinv and Mij are the block inverses of the particle and mesh, respectively.
			// Ciu is the coupling element.
			//
			// Ciu is stored as [ ncoords * (6) + 2 (u/v dot) x 3 (cartesian) 
			//
			// tM[i,v] = Ciu[i u] P[u,v]

			double *C_Pinv = (double *)malloc( sizeof(double) * 3 * ncoef * 2 );
			memset( C_Pinv, 0, sizeof(double) * 3 * ncoef * 2 );

			// right multiply by Pinv.
			for( int c = 0; c < ncoef; c++ )
			for( int i = 0; i < 2; i++ )
			{
				for( int j = 0; j < 2; j++ )
				{
					C_Pinv[c*3*2+i*3+0] += C_iu[c*3*2+j*3+0] * Pinv[i*2+j]; 
					C_Pinv[c*3*2+i*3+1] += C_iu[c*3*2+j*3+1] * Pinv[i*2+j]; 
					C_Pinv[c*3*2+i*3+2] += C_iu[c*3*2+j*3+2] * Pinv[i*2+j]; 
				}
			}

			// COMPUTE qdot
	
			// Pcouple is now the piece of the matrix inverse coupling the mesh and particle at first order.	
			
			// the uncoupled component.
			qdot[2*s+0] = qdot0[0];
			qdot[2*s+1] = qdot0[1];

			zero_order[0] += qdot0[0];
			zero_order[1] += qdot0[1];
			// TERM 1: it has one derivative: qdot_nu Pinv (dP/dX) Pinv qdot_mu
			// computed as p (dP/dX) p
			
			// term has both mu/nu 1A and qi 1B derivatives (see d_H_d_q)
			

#ifdef DO_FIRST_ORDER
			// TERM 2: this is - Pinv_{mu nu} C{nu i} Pinv_{i j} qdot j 	

			// Pinv_{mu nu} has a derivative wrt mu/nu, ihis is 2A
			// Pinv_{q} has a derivative wrt mu/nu, ihis is 2B

			// C{nu i} has a derivative wrt mu/nu, this is 2Ci and 2Cii (there are two terms to take derivatives of).
			// C{nu i} has a derivative wrt q, this is 2D
		
#ifndef TURN_OFF_TERM2

			for( int uv = 0; uv < 2; uv++ )
			for( int uv2 = 0; uv2 < 2; uv2++ )
			for( int c = 0; c < ncoef; c++ )
			{	
				qdot[2*s+uv] -= Pinv[uv*2+uv2] * C_iu[c*3*2+uv2*3+0] * mesh_qdot0[3*coef_list[c]+0];
				qdot[2*s+uv] -= Pinv[uv*2+uv2] * C_iu[c*3*2+uv2*3+1] * mesh_qdot0[3*coef_list[c]+1];
				qdot[2*s+uv] -= Pinv[uv*2+uv2] * C_iu[c*3*2+uv2*3+2] * mesh_qdot0[3*coef_list[c]+2];
			
				first_order[uv] -= Pinv[uv*2+uv2] * C_iu[c*3*2+uv2*3+0] * mesh_qdot0[3*coef_list[c]+0];
				first_order[uv] -= Pinv[uv*2+uv2] * C_iu[c*3*2+uv2*3+1] * mesh_qdot0[3*coef_list[c]+1];
				first_order[uv] -= Pinv[uv*2+uv2] * C_iu[c*3*2+uv2*3+2] * mesh_qdot0[3*coef_list[c]+2];
			}	

			// compute contribution to mesh_qdot.
			// This term has derivatives of the coupling element wrt u,v
			// TERM 3Ai, derivative of C_iu wrt u,v part 1
			// TERM 3Aii, derivative of C_iu wrt u,v part 2
			// TERM 3B, derivative of C_iu wrt mesh.
			// TERM 3C, derivative of inverse wrt uv
			// TERM 3D, derivative of inverse wrt mesh.

			for( int c = 0; c < ncoef; c++ )
			for( int uv = 0; uv < 2; uv++ )
			{
				mesh_qdot_temp[3*coef_list[c]+0] -= C_iu[c*3*2+uv*3+0] * qdot0[uv] * frac_mult; 
				mesh_qdot_temp[3*coef_list[c]+1] -= C_iu[c*3*2+uv*3+1] * qdot0[uv] * frac_mult; 
				mesh_qdot_temp[3*coef_list[c]+2] -= C_iu[c*3*2+uv*3+2] * qdot0[uv] * frac_mult; 
			}
#endif
			// mesh/mesh coupling. 
			// TERM 4, derivative of coefs wrt u and v. (This has no control point derivatives)



#ifndef TURN_OFF_TERM4
			for( int c = 0; c < ncoef; c++ )
			for( int c2 = 0; c2 < ncoef; c2++ )
			{
				mesh_qdot_temp[3*coef_list[c]+0] -= mass[s] * coefs[c] * coefs[c2] * mesh_qdot0[3*coef_list[c2]+0]*frac_mult;
				mesh_qdot_temp[3*coef_list[c]+1] -= mass[s] * coefs[c] * coefs[c2] * mesh_qdot0[3*coef_list[c2]+1]*frac_mult;
				mesh_qdot_temp[3*coef_list[c]+2] -= mass[s] * coefs[c] * coefs[c2] * mesh_qdot0[3*coef_list[c2]+2]*frac_mult;
			}	
#endif

#if 0
			for( int c = 0; c < nv; c++ )
			{
				mesh_qdot[3*c+0] += Pcouple[c*6+0*3+0] * p[2*s+0];
				mesh_qdot[3*c+1] += Pcouple[c*6+0*3+1] * p[2*s+0];
				mesh_qdot[3*c+2] += Pcouple[c*6+0*3+2] * p[2*s+0];
				
				mesh_qdot[3*c+0] += Pcouple[c*6+1*3+0] * p[2*s+1];
				mesh_qdot[3*c+1] += Pcouple[c*6+1*3+1] * p[2*s+1];
				mesh_qdot[3*c+2] += Pcouple[c*6+1*3+2] * p[2*s+1];
			}
#endif
#endif

//			if( s == 0 )
//				printf("p: %le %le qdot: %le %le KE: %le\n", p[2*s+0], p[2*s+1], qdot[2*s+0], qdot[2*s+1], 0.5 * (p[2*s+0]*qdot[2*s+0] + p[2*s+1]*qdot[2*s+1]) );

//			printf("ZO %le %le FO %le %le\n", zero_order[0], zero_order[1], first_order[0], first_order[1] );

			free(C_Pinv);
			free(C_iu);
#if 0
			free(Pcouple);
#endif
		}
	}

	for( int s = nattach; s < nsites; s++ )
	{
		qdot[3*s+0] = p[3*s+0] / mass[s];
		qdot[3*s+1] = p[3*s+1] / mass[s];
		qdot[3*s+2] = p[3*s+2] / mass[s];
	}
}

void pcomplex::propagate_surface_q( Simulation *theSimulation,  double dt )
{
	if( disabled )
		return;

	static double av_dstep = 0;
	static double av_dstep2 = 0;
	static double n_dstep = 0;

	static double av_root_g = 0;
	static double av_root_g2 = 0;
	static double n_av_root_g = 0;
	// move a step.

	for( int s = 0; s < nattach; s++ )
	{
		int okay = 1;
		int check_iter = 0;
 
		surface_record *sRec = NULL;
		surface *theSurface =NULL;
		for( surface_record *tRec = theSimulation->allSurfaces; tRec; tRec = tRec->next )
		{
			if( tRec->id == sid[s] )
			{
				sRec = tRec;	
				theSurface = sRec->theSurface;
				break;
			}
		}

		double *rsurf = sRec->r;
		double *mesh_grad = sRec->g;
		double *mesh_qdot = sRec->qdot;
		double *mesh_qdot0 = sRec->qdot0;
		int last_f = fs[s];
		double last_metric = theSurface->g( fs[s], puv[2*s+0], puv[2*s+1], rsurf );
		double duv[2];

		if( do_bd && !rd_timestep_disabled[s] )
		{
			int was_irreg = fs[s] >= theSurface->nf_faces;
			double curp[3] = { rall[0], rall[1], rall[2] };
			double gmat[4];
			double g_u; // derivative of |g|, not root wrt u
			double g_v; //                             wrt v
			
	
			// gval is the sqrt metric
			double root_g = theSurface->metric( fs[s], puv[2*s+0], puv[2*s+1], rsurf, gmat, &g_u, &g_v );

#ifdef TRACK_UNUSUAL
			if( n_av_root_g > 10 )
			{
				double var = av_root_g2/n_av_root_g - pow(av_root_g/n_av_root_g,2);
				double stdd = sqrt(var);

				if( fabs( (root_g - av_root_g/n_av_root_g)/stdd ) > 3 )
				{
					printf("Unusual root g %le av %le stdd %le\n",
						root_g, av_root_g/n_av_root_g, stdd );
				}
			}

			av_root_g2 += root_g*root_g;
			av_root_g += root_g;
			n_av_root_g += 1;
#endif
			double val = sqrt(gmat[0]*gmat[0]+4*gmat[1]*gmat[1]-2*gmat[0]*gmat[3]+gmat[3]*gmat[3]);

			double gamma_neg_half[4] = { 
				(1.0/val) * ((-gmat[0]+gmat[3]+val)/sqrt(2*(gmat[0]+gmat[3]-val))-(-gmat[0]+gmat[3]-val)/sqrt(2*(gmat[0]+gmat[3]+val))),   	
				(1.0/val) * ( -sqrt(2)*gmat[1]/sqrt(gmat[0]+gmat[3]-val) + sqrt(2)*gmat[1]/sqrt(gmat[0]+gmat[3]+val) ),
				(1.0/val) * ( -sqrt(2)*gmat[1]/sqrt(gmat[0]+gmat[3]-val) + sqrt(2)*gmat[1]/sqrt(gmat[0]+gmat[3]+val) ),
				(1.0/val) * (( gmat[0]-gmat[3]+val)/sqrt(2*(gmat[0]+gmat[3]-val))+(-gmat[0]+gmat[3]+val)/sqrt(2*(gmat[0]+gmat[3]+val))),   
				};

			double gamma_inv[4] = {
				gmat[0]/(root_g*root_g), -gmat[1]/(root_g*root_g), -gmat[1]/(root_g*root_g), gmat[3]/(root_g*root_g)
			};

			double ru[3],rv[3];
			theSurface->ru( fs[s], puv[2*s+0], puv[2*s+1], rsurf, ru );
			theSurface->rv( fs[s], puv[2*s+0], puv[2*s+1], rsurf, rv );
//			printf("ru: %.14lf %.14lf %.14lf rv: %.14lf %.14lf %.14lf\n", ru[0], ru[1], ru[2], rv[0], rv[1], rv[2] );

			// g_u + means metric increases in u direction, move particle there.
			double pre_drift[2] = { 
					save_grad[2*s+0]/kT + (-g_u/(2*root_g*root_g)), // 1/g done with sqrt
					save_grad[2*s+1]/kT + (-g_v/(2*root_g*root_g))
						};

			double tu = puv[2*s+0]; 
			double tv = puv[2*s+1]; 
			int fu = fs[s];

#ifndef PRINT_PTS
			int ntests = 0;
#else
			int ntests = 500;
			printf("%d\n", ntests+1);
			printf("test\n");
			printf("O %le %le %le\n", rall[0], rall[1], rall[2] );	
#endif
			check_iter=0;
			do {
				fs[s] = fu;
				puv[2*s+0] = tu;
				puv[2*s+1] = tv;

				okay = 1;

				double noise[2] = { 
					gsl_ran_gaussian(r_gen_global, 1 ),
					gsl_ran_gaussian(r_gen_global, 1 )
					};
	
				double corr_noise[2] = { gamma_neg_half[0] * noise[0] + gamma_neg_half[1] * noise[1],
							 gamma_neg_half[2] * noise[0] + gamma_neg_half[3] * noise[1] }; 
			
				// negative sign means move in the direction of force not grad.
				duv[0] = -(gamma_inv[0] * pre_drift[0] + gamma_inv[1] * pre_drift[1]) * DC[s] * dt + sqrt(2 * DC[s] * dt) * corr_noise[0];
				duv[1] = -(gamma_inv[2] * pre_drift[0] + gamma_inv[3] * pre_drift[1]) * DC[s] * dt + sqrt(2 * DC[s] * dt) * corr_noise[1];	
				
				puv[2*s+0] += duv[0];
				puv[2*s+1] += duv[1];				

				/*
				if( theSimulation->rd && stype && stype[s] >= 0 && check_iter < 100 )
				{
					refresh(theSimulation);
					// this particle is configured for reaction/diffusion.
					// we must block it from entering the reaction zone for all reactions for which it can participate.
		
					if( theSimulation->rd->check_RD_blocked( theSimulation, my_id, s ) )
					{
	//						printf("BLOCKED RD MOVE.\n");
						okay = 0; 
					}
					check_iter++;

					if( check_iter > 50 )
					{
					}
				}*/
				refresh(theSimulation);

//#define DEBUG_DIFFUSION		
#ifdef DEBUG_DIFFUSION
		if( fabs(dt-1e-12) < 1e-11)
		{
			double newp[3] = { rall[0], rall[1], rall[2] };
			double dr[3] = {newp[0]-curp[0],newp[1]-curp[1],newp[2]-curp[2] };
			printf("%d DEL %le dx: %le dy: %le\n", was_irreg, dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2], dr[0], dr[1] );
		}
#endif		

			
#ifdef PRINT_PTS
				printf("C %le %le %le\n", rall[0], rall[1], rall[2] );
#endif
				check_iter++;
			} while ( !okay || check_iter < ntests);
		}
		else
		{
			duv[0] = qdot[2*s+0] * dt * AKMA_TIME;
			duv[1] = qdot[2*s+1] * dt * AKMA_TIME;
		
			puv[2*s+0] += duv[0];
			puv[2*s+1] += duv[1];
		}



	
		if( check_iter >= 100 )
			printf("totally blocked.\n");	

		//double new_metric = theSurface->g( fs[s], puv[2*s+0], puv[2*s+1], rsurf );
		//double fr = fabs(1-new_metric/last_metric);
	}
	
	if( bound )	
	{
		refresh(theSimulation);
		// 

		if( puv[0] < 0 || puv[1] < 0 || puv[0]+puv[1] > 1 )
			printf("after propagate and setrall, uv: %le %le\n",puv[0], puv[1] );
	}
	else
	{

	} 
#if 0
#define VERIFY_BLOCKED_CHECK
#ifdef VERIFY_BLOCKED_CHECK
	int nblocked = 0;
	int I_was_blocked = 0;
	for( int c = 0; c < theSimulation->ncomplex; c++ )
	{
		for( int s = 0; s < theSimulation->allComplexes[c]->nsites; s++ )
		{
			if( theSimulation->rd && theSimulation->allComplexes[c]->stype[s] >= 0 )
			{
				if( theSimulation->rd->check_RD_blocked( theSimulation, c, s ) ) {
					nblocked++;
	
					if( c == my_id )
						I_was_blocked=1;
				}
			}
		}
	}
	if( nblocked > 0 )
	{
		printf("after propagate NBLOCKED: %d", nblocked );
		if( I_was_blocked) printf(" including me.");
		printf("\n");
	}
#endif	
#endif	
}

void pcomplex::propagate_p( Simulation *theSimulation, double dt )
{
	if( disabled ) return;

		// move a step.

#if 0 // presim
	double T_before = T(theSurface,rsurf);
#endif 
	for( int s = 0; s < nattach; s++ )
	{
		p[2*s+0] -= save_grad[2*s+0] * dt * AKMA_TIME;
		p[2*s+1] -= save_grad[2*s+1] * dt * AKMA_TIME;
	}
	
	for( int s = nattach; s < nsites; s++ )
	{
		p[3*s+0] -= save_grad[3*s+0] * dt * AKMA_TIME;
		p[3*s+1] -= save_grad[3*s+1] * dt * AKMA_TIME;
		p[3*s+2] -= save_grad[3*s+2] * dt * AKMA_TIME;
	}

#if 0 //presim	
	double T_after = T(theSurface,rsurf);

	if( fabs(T_after-T_before) > 0.592  )
	{
		printf("dKE: %le f: %d u: %le v: %le\n", T_after-T_before, fs[0], puv[0], puv[1] );
	}
#endif
	if( bound )	
	{
		// 
	}
	else
	{

	} 
}

void pcomplex::pfree( void )
{
	free(DC);
	free(rall);
	free(fs);
	free(puv);
	free(p);
	free(qdot);
	free(mass);
	free(grad_fs);
	free(grad_puv);
}

void pcomplex::init( double *r )
{
	base_init();

	nsites = 1;
	nattach = 0;
	
	alloc();

	rall[0] = r[0];	
	rall[1] = r[1];	
	rall[2] = r[2];	

	mass[0] = default_mass;

	bound = 0;
}

void pcomplex::init( surface *theSurface, double *rsurf, int f, double u, double v )
{
	base_init();

	nsites = 1;
	nattach = 1;

	alloc();

	grad_fs[0] = fs[0] = f;
	grad_puv[0] = puv[0] = u;
	grad_puv[1] = puv[1] = v;
	
	mass[0] = default_mass;



	bound = 1;
}

void pcomplex::bind( int f, double u, double v )
{
	bound = 1;

	fs[0] = f;
	puv[0] = u;
	puv[1] = v;
}

void pcomplex::unbind( void )
{
	bound = 0;
}

void pcomplex::applyParamsToComplex( double *p )
{
	if( bound )
	{
		for( int s = 0; s < nattach; s++ )
		{
			puv[2*s+0]=p[2*s+0];
			puv[2*s+1]=p[2*s+1];
		}

		for( int t = nattach; t < nsites; t++ )
		{
			rall[3*t+0]=p[2*nattach+3*(t-nattach)+0];
			rall[3*t+1]=p[2*nattach+3*(t-nattach)+1];
			rall[3*t+2]=p[2*nattach+3*(t-nattach)+2];
		}
	}
	else
	{
		for( int t = 0; t < nsites; t++ )
		{
			rall[3*t+0]=p[3*(t)+0];
			rall[3*t+1]=p[3*(t)+1];
			rall[3*t+2]=p[3*(t)+2];
		}
	}
}

void pcomplex::getParamsFromComplex( double *p )
{
	if( bound )
	{
		for( int s = 0; s < nattach; s++ )
		{
			p[2*s+0] = puv[2*s+0];
			p[2*s+1] = puv[2*s+1];
		}

		for( int t = nattach; t < nsites; t++ )
		{
			p[2*nattach+3*(t-nattach)+0] = rall[3*t+0];
			p[2*nattach+3*(t-nattach)+1] = rall[3*t+1];
			p[2*nattach+3*(t-nattach)+2] = rall[3*t+2];
		}
	}
	else
	{
		for( int t = 0; t < nsites; t++ )
		{
			p[3*(t)+0] = rall[3*t+0];
			p[3*(t)+1] = rall[3*t+1];
			p[3*(t)+2] = rall[3*t+2];
		}
	}
}

int pcomplex::nparams( void )
{
	if( bound )
		return 2*nattach + 3 * (nsites-nattach);
	else
		return 3 * nsites;
}

void pcomplex::refresh( Simulation *theSimulation )
{
	if( disabled ) return;

	double *alphas = theSimulation->alpha;
	double center_on[3] = { 0,0,0};
	int update_attach = 1;

	if( nattach < nsites )
	{
		update_attach = 0;
		center_on[0] = rall[3*nattach+0];
		center_on[1] = rall[3*nattach+1];
		center_on[2] = rall[3*nattach+2];
	}

	for( int s = 0; s < nattach; s++ )
	{
		surface_record *sRec = theSimulation->fetch(sid[s]);

		surface *theSurface = sRec->theSurface;
		double *rsurf = sRec->r;

		double nrm[3];

		int o_f = fs[s];
		int f_1 = fs[s], nf = fs[s];
		double uv1[2] = { 0.33, 0.33 };
		double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };

		double momentum_uv[2] = { qdot[2*s+0], qdot[2*s+1] }; 

		if( puv[2*s+0] <= 0 || puv[2*s+1] <= 0 || puv[2*s+0] + puv[2*s+1] >= 1.0 )
		{
			double ro[3];
			theSurface->evaluateRNRM( fs[s], 0.33, 0.33, ro, nrm, rsurf );  

			do {
				f_1 = nf;
				nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf, momentum_uv ); 
			} while( nf != f_1 );

			uv1[0] += duv1[0];		
			uv1[1] += duv1[1];		
			
			double KE_before = qdot[2*s+0] * p[2*s+0] + qdot[2*s+1] * p[2*s+1];

			qdot[2*s+0] = momentum_uv[0];
			qdot[2*s+1] = momentum_uv[1];

		
			fs[s] = grad_fs[s] = f_1;
		
			puv[2*s+0] = grad_puv[2*s+0] = uv1[0];
			puv[2*s+1] = grad_puv[2*s+1] = uv1[1];
			
			theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], rall+3*s, nrm, rsurf );  
				
			double dr[3] = { rall[3*s+0] - ro[0], rall[3*s+1] - ro[1], rall[3*s+2] - ro[2] };
			double del[3];
			MinImage3D( dr, theSurface->PBC_vec, del, rsurf+theSurface->nv*3 );

			PBC_ext[3*s+0] += del[0];
			PBC_ext[3*s+1] += del[1];
			PBC_ext[3*s+2] += del[2];

			// re-compute p.
#ifdef DO_FIRST_ORDER
//			printf("THIS MUST BE IMPLEMENTED.\n");
//			exit(1);
#endif	
			double P_uv[4] = {
				0,	0,
				0,	0
			};
	
			theSurface->fetchPuv( fs[s], puv[2*s+0], puv[2*s+1], P_uv, rsurf );
	
			P_uv[0] *= mass[s];
			P_uv[1] *= mass[s];
			P_uv[2] *= mass[s];
			P_uv[3] *= mass[s];

			double ppre[2] = { p[2*s+0], p[2*s+1] };
			p[2*s+0] = P_uv[0] * qdot[2*s+0] + P_uv[1] * qdot[2*s+1];
			p[2*s+1] = P_uv[2] * qdot[2*s+0] + P_uv[3] * qdot[2*s+1];
			
			double KE_after = qdot[2*s+0] * p[2*s+0] + qdot[2*s+1] * p[2*s+1];
		
//			printf("p before %le %le p after %le %le KE before: %le KE after: %le\n", ppre[0], ppre[1], p[2*s+0], p[2*s+1], KE_before, KE_after );
		}
		else
			theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], rall+3*s, nrm, rsurf );  

/*	
		double dr[3] = { rall[3*s+0] - last_pos[3*s+0],
				 rall[3*s+1] - last_pos[3*s+1],
				 rall[3*s+2] - last_pos[3*s+2] };
		double del[3] = {0,0,0};

		MinImage3D( dr, theSurface->PBC_vec, del, rsurf+theSurface->nv*3 );

		PBC_ext[3*s+0] = del[0];
		PBC_ext[3*s+1] = del[1];
		PBC_ext[3*s+2] = del[2];
*/		
		rall[3*s+0] += theSimulation->PBC_vec[0][0] * PBC_ext[3*s+0] * alphas[0] + theSimulation->PBC_vec[1][0] * PBC_ext[3*s+1] * alphas[0] + theSimulation->PBC_vec[2][0] * PBC_ext[3*s+2] * alphas[0];
		rall[3*s+1] += theSimulation->PBC_vec[0][1] * PBC_ext[3*s+0] * alphas[1] + theSimulation->PBC_vec[1][1] * PBC_ext[3*s+1] * alphas[1] + theSimulation->PBC_vec[2][1] * PBC_ext[3*s+2] * alphas[1];
		rall[3*s+2] += theSimulation->PBC_vec[0][2] * PBC_ext[3*s+0] * alphas[2] + theSimulation->PBC_vec[1][2] * PBC_ext[3*s+1] * alphas[2] + theSimulation->PBC_vec[2][2] * PBC_ext[3*s+2] * alphas[2];


	}

	if( puv[0] < 0 || puv[1] < 0 || puv[0]+puv[1] > 1.0 )
	{
		printf("DOMAIN ERROR.\n");
	}

	memcpy( last_pos, rall, sizeof(double) * nsites*3);
}

pcomplex *loadComplex( const char *name )
{
	pcomplex *the_complex = NULL;

	if( !strcasecmp( name, "NBAR" ) )
		the_complex = new NBAR;
	else if( !strcasecmp( name, "dimer" ) )
		the_complex = new dimer;
	else if( !strcasecmp( name, "MAB" ) )
		the_complex = new MAB;
	else if( !strcasecmp( name, "crowder" ) )
		the_complex = new elasticCrowder;
	else if( !strcasecmp( name, "simpleParticle" ) )
		the_complex = new simpleParticle;
	else if( !strcasecmp( name, "simpleBound" ) )
		the_complex = new simpleBound;
	else if( !strcasecmp( name, "simpleLipid" ) )
		the_complex = new simpleLipid;
	else if( !strcasecmp( name, "simpleDimer" ) )
		the_complex = new simpleDimer;
	else
	{
		printf("Unknown complex name '%s'.\n", name );
		exit(1);
	}

	the_complex->complex_name = (char *)malloc( sizeof(char) * (1 + strlen(name) ) );
	strcpy( the_complex->complex_name, name );

	return the_complex;
}


double pcomplex::V( Simulation *theSimulation )
{
	return 0;
}

double pcomplex::grad( Simulation *theSimulation, double *surface_g, double *particle_g )
{
	memset( particle_g, 0, sizeof(double) * 3 * nattach );
	memset( surface_g, 0, sizeof(double) * 3 * nattach );

	for( int s = 0; s < nattach; s++ )
	{
		grad_fs[s] = fs[s];
		grad_puv[2*s+0] = puv[2*s+0];
		grad_puv[2*s+1] = puv[2*s+1];
			
		coord_transform[4*s+0] = 1;
		coord_transform[4*s+1] = 0;
		coord_transform[4*s+2] = 0;
		coord_transform[4*s+3] = 1;	
	}	

	// return the energy.
	return 0;
}

void surface::loadComplexes( pcomplex ***allComplexes, int *ncomplex, parameterBlock *block )
{
	int nspace = 1;
	(*allComplexes) = (pcomplex **)malloc( sizeof(pcomplex *) );

	double **M;
	int mlow = 5;
	int mhigh = 7;
	getM( &M, &mlow, &mhigh );
	double *rsurf = (double *)malloc( sizeof(double) * ( 3*nv+3) );
	get(rsurf);

	rsurf[3*nv+0] = 1.0;
	rsurf[3*nv+1] = 1.0;
	rsurf[3*nv+2] = 1.0;

	double *alphas = rsurf+3*nv;

	double the_area, area0;
	
	area( rsurf, -1, &the_area, &area0 ); 

	double vol_inside = volume(rsurf);
	double total_volume = cellVolume();	
	double vol_outside = total_volume - vol_inside;

	*ncomplex = 0;

	for( complex_record *rec = block->complex_records; rec; rec = rec->next )
	{
		int num_inside_sol = 0;
		int num_outside_sol = 0;
		int num_anywhere_sol = 0;
		
		int num_inside_surf = 0;
		int num_outside_surf = 0;
		int num_anywhere_surf = 0;

		if( rec->inside_outside == 1 )
		{
			num_outside_sol = rec->nsolution;
			num_outside_surf = rec->nbound;
			num_outside_sol += rec->concentration * vol_outside;
			num_outside_surf += rec->coverage * the_area;
		}
		else if( rec->inside_outside == -1 )
		{
			num_inside_sol = rec->nsolution;
			num_inside_surf = rec->nbound;
			num_inside_sol += rec->concentration * vol_inside;
			num_inside_surf += rec->coverage * the_area;
		}
		else
		{
			num_anywhere_sol = rec->nsolution;
			num_anywhere_surf = rec->nbound;
			num_anywhere_sol += rec->concentration * total_volume;
			num_anywhere_surf += rec->coverage * the_area;
		}

		int nsol_tot = num_inside_sol + num_outside_sol + num_anywhere_sol;
	
		for( int p = 0; p < nsol_tot; p++ )
		{
			int io_check = 0;
		
			int test = 0; // anywhere

			if( p < num_inside_sol )
				test = -1;
			else if( p < num_inside_sol + num_outside_sol )
				test = 1;
			else
				test = 0;

			double tp[3];

			while( !io_check )
			{
				double f1 = rand() / (double)RAND_MAX;
				double f2 = rand() / (double)RAND_MAX;
				double f3 = rand() / (double)RAND_MAX;

				tp[0] = f1 * PBC_vec[0][0] * alphas[0] + f2 * PBC_vec[1][0] * alphas[0] + f3 * PBC_vec[2][0] * alphas[0];
				tp[1] = f1 * PBC_vec[0][1] * alphas[1] + f2 * PBC_vec[1][1] * alphas[1] + f3 * PBC_vec[2][1] * alphas[1];
				tp[2] = f1 * PBC_vec[0][2] * alphas[2] + f2 * PBC_vec[1][2] * alphas[2] + f3 * PBC_vec[2][2] * alphas[2];
				
				if( test != 0 )
				{
					int fcol;
					double ucol,vcol;
					double dummy;

					if( withinBoxedSurface( tp, &fcol, &ucol, &vcol, M, mlow, mhigh, &dummy ) )
					{
						if( test == -1 )
							io_check = 1;
					}
					else if( test == 1 )
						io_check = 1;
				}
				else
					io_check = 1; 
			}

			pcomplex *prot = loadComplex( rec->name );
			prot->loadParams(block);
			prot->init(tp); 
		
			if( *ncomplex  == nspace )
			{
				nspace *= 2;
				*allComplexes = (pcomplex **)realloc( *allComplexes, sizeof(pcomplex *) * nspace );
			}	
		
			(*allComplexes)[*ncomplex] = prot;
			prot->my_id = *ncomplex;
			(*ncomplex)++;
		}
		
		int nsurf_tot = num_inside_surf + num_outside_surf + num_anywhere_surf;

		for( int p = 0; p < nsurf_tot; p++ )
		{
			int f;
			double u,v;
		
			randomPointOnSurface( &f, &u, &v );			

			pcomplex *prot = loadComplex( rec->name );

			prot->loadParams(block);
			if( p < num_inside_surf || (p > num_inside_surf + num_outside_surf && rand() % 2 == 0) )
				prot->move_inside();
			else
				prot->move_outside();

			prot->init(this, rsurf, f,u,v); 				 

			for( int t = 0; t < prot->nattach; t++ )
				prot->sid[t] = surface_id; 

			if( *ncomplex  == nspace )
			{
				nspace *= 2;
				*allComplexes = (pcomplex **)realloc( *allComplexes, sizeof(pcomplex *) * nspace );
			}	

			(*allComplexes)[*ncomplex] = prot;
			prot->my_id = *ncomplex;

			(*ncomplex)++;
		}

	}	

	free(rsurf);
}


void pcomplex::fd_grad_debug(surface *theSurface, double *rsurf )
{
	double eps = 5e-4;
	
	int nfree = nsites-nattach;
 
	double g_attach[3*nattach];
	double g_free[3*nfree];

	for( int c = 0; c < nattach; c++ )
	{
		
	}

	for( int c = 0; c < nfree; c++ )
	{
	}
}

void pcomplex::loadParams( parameterBlock *block )
{
	
	
}

void simpleLipid::loadParams( parameterBlock *block )
{
	c0_val = block->c0;
}

void simpleDimer::loadParams( parameterBlock *block )
{
	c0_val = block->dimer_c0;
}

void simpleBound::loadParams( parameterBlock *block )
{
	bound_sigma = block->bound_sigma;
}



void simpleLipid::init( surface *theSurface, double *rsurf, int f, double u, double v )
{
	pcomplex::init( theSurface, rsurf, f, u, v );

	p_c0[0] = c0_val;
}

void simpleDimer::init( surface *theSurface, double *rsurf, int f, double u, double v )
{
	simpleLipid::init( theSurface, rsurf, f, u, v );

	p_c0[0] = c0_val;
	p_area[0] = 2 * default_particle_area;
}

void simpleBound::init( surface *theSurface, double *rsurf, int f, double u, double v )
{
	base_init();

	nsites = 2;
	nattach = 1;

	alloc();

	grad_fs[0] = fs[0] = f;
	grad_puv[0] = puv[0] = u;
	grad_puv[1] = puv[1] = v;
	
	mass[0] = default_mass;
	mass[1] = default_mass;

	bound = 1;
	
	double rp[3];
	double nrm[3];

	DC[0] = lipid_DC;
	DC[1] = solution_DC;
	
	grad_puv[0] = puv[0] = u;
	grad_puv[1] = puv[1] = v;
	grad_fs[0] = fs[0] = f;

	theSurface->evaluateRNRM( f, u, v, rp, nrm, rsurf);
	
	rall[0] = rp[0];
	rall[1] = rp[1];
	rall[2] = rp[2];
	
	// ``solution'' particle.

	rall[3] = rp[0] + bound_sigma * nrm[0];
	rall[4] = rp[1] + bound_sigma * nrm[1];
	rall[5] = rp[2] + bound_sigma * nrm[2];
}


void pcomplex::evaluate_momentum( surface *theSurface, double *rsurf, double *pout )
{
	for( int s = 0; s < nattach; s++ )
	{
		// du dt
		double ru[3], rv[3];

		theSurface->ru( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], rsurf, ru );
		theSurface->rv( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], rsurf, rv );

		pout[0] += mass[s] * ( ru[0] * qdot[2*s+0] + rv[0] * qdot[2*s+1]);
		pout[1] += mass[s] * ( ru[1] * qdot[2*s+0] + rv[1] * qdot[2*s+1]);
		pout[2] += mass[s] * ( ru[2] * qdot[2*s+0] + rv[2] * qdot[2*s+1]);
	}
}

double pcomplex::T( Simulation *theSimulation, int subp)
{
	double T = 0;

	if( do_bd ) return 0;

	for( int s = 0; s < nattach; s++ )
	{
		if( subp != -1 && s!= subp) continue;
		surface_record *sRec = theSimulation->fetch(sid[s]);
		surface *theSurface = sRec->theSurface;
		double *rsurf = sRec->r;

		int f = fs[s];
		double u = puv[2*s+0];
		double v = puv[2*s+1];
		double P_uv[4] = {0,0,0};

		theSurface->fetchPuv( f, u, v, P_uv, rsurf );
		
		P_uv[0] *= mass[s];
		P_uv[1] *= mass[s];
		P_uv[2] *= mass[s];
		P_uv[3] *= mass[s];
	
		double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
		double Pinv[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };
	
		double qdot0[2] = { Pinv[0] * p[2*s+0] + Pinv[1] * p[2*s+1],
				    Pinv[2] * p[2*s+0] + Pinv[3] * p[2*s+1] };

		// T = 0.5 * p * qdot
		T += 0.5 * (p[2*s+0]) * (qdot0[0]);
		T += 0.5 * (p[2*s+1]) * (qdot0[1]);
	}

	for( int s = nattach; s < nsites; s++ )
	{
		if( subp != -1 && s!= subp) continue;
		T += 0.5 * (p[3*s+0]) * (qdot[3*s+0]);
		T += 0.5 * (p[3*s+1]) * (qdot[3*s+1]);
		T += 0.5 * (p[3*s+2]) * (qdot[3*s+2]);
	}

	return T;
}

void pcomplex::putBonds( int *bond_list )
{
}

#if 0 // pre sim
void pcomplex::debug_dynamics( surface *theSurface,
		      double *rsurf,
		
		      double *mesh_qdot,
		      double *mesh_qdot0,

		      double *dV_dmesh, 
		      double *dT_dmesh,
		      double *dV_dcomplex,
		      double *dT_dcomplex )
{	
	// this retrieves dV_dq and dT_dq

	compute_qdot( theSurface, rsurf, mesh_qdot0, mesh_qdot );	
	debug = DEBUG_NO_V;
	for( int s = 0; s < nattach; s++ )
	{
		save_grad[2*s+0] = 0;
		save_grad[2*s+1] = 0;
	}
	update_dH_dq( theSurface, rsurf, dT_dmesh, mesh_qdot, mesh_qdot0 ); 
	for( int s = 0; s < nattach; s++ )
	{
		dT_dcomplex[2*s+0] += save_grad[2*s+0];
		dT_dcomplex[2*s+1] += save_grad[2*s+1];
	}
	debug = DEBUG_NO_T;
	
	for( int s = 0; s < nattach; s++ )
	{
		save_grad[2*s+0] = 0;
		save_grad[2*s+1] = 0;
	}
	update_dH_dq( theSurface, rsurf, dV_dmesh, mesh_qdot, mesh_qdot0 ); 
	for( int s = 0; s < nattach; s++ )
	{
		dV_dcomplex[2*s+0] += save_grad[2*s+0];
		dV_dcomplex[2*s+1] += save_grad[2*s+1];
	}
}
#endif

void pcomplex::applyLangevinFriction( Simulation *theSimulation, double dt, double gamma )
{

	for( int s = 0; s < nattach; s++ )
	{
		double ppre[2] = { p[2*s+0], p[2*s+1] };
		p[2*s+0] -= qdot[2*s+0] * gamma * AKMA_TIME * dt;
		p[2*s+1] -= qdot[2*s+1] * gamma * AKMA_TIME * dt;
	}
		
	for( int s = nattach; s < nsites; s++ )
	{
		p[3*s+0] -= qdot[3*s+0] * gamma * AKMA_TIME * dt;
		p[3*s+1] -= qdot[3*s+1] * gamma * AKMA_TIME * dt;
		p[3*s+2] -= qdot[3*s+2] * gamma * AKMA_TIME * dt;
	}
}

void pcomplex::applyLangevinNoise( Simulation *theSimulation,  double dt, double gamma, double temperature )
{

	for( int s = 0; s < nattach; s++ )
	{
		double fx = gsl_ran_gaussian(r_gen_global, sqrt(2*gamma*temperature*AKMA_TIME*dt) );
		double fy = gsl_ran_gaussian(r_gen_global, sqrt(2*gamma*temperature*AKMA_TIME*dt) );

		p[2*s+0] += fx;
		p[2*s+1] += fy;
	}
	
	for( int s = nattach; s < nsites; s++ )
	{
		double fx = gsl_ran_gaussian(r_gen_global, sqrt(2*gamma*temperature*AKMA_TIME*dt) );
		double fy = gsl_ran_gaussian(r_gen_global, sqrt(2*gamma*temperature*AKMA_TIME*dt) );
		double fz = gsl_ran_gaussian(r_gen_global, sqrt(2*gamma*temperature*AKMA_TIME*dt) );

		p[3*s+0] += fx;
		p[3*s+1] += fy;
		p[3*s+2] += fz;
	}
}

void pcomplex::loadComplex( FILE *theFile, Simulation *theSimulation, int load_to )
{
	char buffer[4096];

	getLine( theFile, buffer );

	int nr = sscanf( buffer, "%d", &bound );

	if( bound )
	{
	}
	else
	{
	}

	for( int s = 0; s < nattach; s++ )
	{
		getLine( theFile, buffer );
		int nr = sscanf( buffer, "%d %lf %lf %lf %lf %lf %lf %lf", fs+s, puv+2*s+0, puv+2*s+1, PBC_ext+3*s, PBC_ext+3*s+1, PBC_ext+3*s+2, p+2*s+0, p+2*s+1 );

		sid[s] = load_to;
		grad_fs[s] = fs[s];
		grad_puv[2*s+0] = puv[2*s+0];
		grad_puv[2*s+1] = puv[2*s+1];

		if( nr != 6 && nr != 8 )
		{
			printf("LINE %s\n", buffer );
			printf("ERROR reading complex.\n");
			exit(1);				
		}

		if( nr == 6 )
		{
			p[2*s+0] = 0;
			p[2*s+1] = 0;
		}
	}
	
	for( int s = nattach; s < nsites; s++ )
	{
		getLine( theFile, buffer );
		int nr = sscanf( buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", rall+3*s+0, rall+3*s+1, rall+3*s+2, PBC_ext+3*s, PBC_ext+3*s+1, PBC_ext+3*s+2, p+3*s+0, p+3*s+1, p+3*s+2 );

		if( nr != 9 )
		{
			printf("LINE %s\n", buffer );
			printf("ERROR reading complex.\n");
			exit(1);				
		}
	}

	refresh(theSimulation);
}

int pcomplex::saveComplex( char *write_to, int *bytes_written, int max_write)
{
	int ptr = 0;

	char *tbuf = (char *)malloc( sizeof(char) * 10000 );

	sprintf( tbuf, "%d\n", bound );

	if( strlen(tbuf) >= max_write )
		return 1;
	strcpy( write_to+ptr, tbuf );
	ptr += strlen(tbuf);	


	for( int s = 0; s < nattach; s++ )
	{
		sprintf(tbuf, "%d %lf %lf %lf %lf %lf %lf %lf\n", fs[s], puv[2*s+0], puv[2*s+1], PBC_ext[3*s+0], PBC_ext[3*s+1], PBC_ext[3*s+2], p[2*s+0], p[2*s+1] );
		if( ptr + strlen(tbuf) >= max_write )
			return 1;
		strcpy( write_to+ptr, tbuf );
		ptr += strlen(tbuf);		
	}
	for( int s = nattach; s < nsites; s++ )
	{
		sprintf(tbuf, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", rall[3*s+0], rall[3*s+1], rall[3*s+2], PBC_ext[3*s], PBC_ext[3*s+1], PBC_ext[3*s+2], p[3*s+0], p[3*s+1], p[3*s+2] );
		
		if( ptr + strlen(tbuf) >= max_write )
			return 1;
		strcpy( write_to+ptr, tbuf );
		ptr += strlen(tbuf);		
	}
	free(tbuf);

	*bytes_written = ptr;

	return 0;
}

void pcomplex::saveComplex( FILE *theFile )
{

	char *buf_write = (char *)malloc( sizeof(char) * 20000 );

	int len = 0;
	int maxWrite = 20000;

	saveComplex( buf_write, &len, maxWrite );

	fprintf( theFile, "%s", buf_write );
		

	free( buf_write );
}


// the surface energy of attachments.
double pcomplex::AttachV( Simulation *theSimulation )
{	
	if( disabled ) return 0;

	double pot = 0;

#ifndef DISABLE_ATTACH
	for( int s = 0; s < nattach; s++ )
	{
		struct surface_record *sRec = theSimulation->fetch(sid[s]);
		double k;
		double curv = sRec->theSurface->c( fs[s],puv[2*s+0], puv[2*s+1], sRec->r, &k);


		pot += 0.5 * kc * p_area[s] * ( curv - p_c0[s]) * (curv - p_c0[s]); 
	}
#endif

	return pot;
}

double pcomplex::AttachG( Simulation *theSimulation,  double *pg )
{
	if( disabled ) return 0;

	double pot = 0;

#ifndef DISABLE_ATTACH
	for( int s = 0; s < nattach; s++ )
	{
		struct surface_record *sRec = theSimulation->fetch(sid[s]);
		sRec->theSurface->particle_H_grad( sRec->r, sRec->g, grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], p_area[s], p_c0[s], pg + 2 *s );
		double k;
		double curv = sRec->theSurface->c( grad_fs[s], grad_puv[2*s+0], grad_puv[2*s+1], sRec->r, &k );

		pot += 0.5 * kc * p_area[s] * ( curv - p_c0[s]) * (curv - p_c0[s]); 
	}
#endif

	return pot;
	
}

double pcomplex::local_curvature( Simulation *theSimulation )
{
	double av = 0;

	for( int s = 0; s < nattach; s++ )
	{
		struct surface_record *sRec = theSimulation->fetch(sid[s]);
		double k;
		double curv = sRec->theSurface->c( fs[s], puv[2*s+0], puv[2*s+1], sRec->r, &k );

		av += curv/nattach;
	}
	return av;
}

void pcomplex::setrall( Simulation *theSimulation  )
{
	double *alphas = theSimulation->alpha ;
	if( bound )
	{
		double njunk[3];
		for( int s = 0; s < nattach; s++ )
		{
			double *rsurf;
			surface *theSurface;
			theSimulation->fetch( sid[s], &theSurface, &rsurf );

			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
			
			if( puv[2*s+0] <= 0 || puv[2*s+1] <= 0 || puv[2*s+0]+puv[2*s+1] >= 1 )
			{
				double ro[3];
				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], ro, njunk, rsurf );  

				do {
					f_1 = nf;
					nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
				} while( nf != f_1 );
	
				uv1[0] += duv1[0];		
				uv1[1] += duv1[1];		
			
				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rall+3*s, njunk, rsurf );  

				double dr[3] = { rall[3*s+0] - ro[0], rall[3*s+1] - ro[1], rall[3*s+2] - ro[2] };
				double del[3];
				MinImage3D( dr, theSurface->PBC_vec, del, rsurf+theSurface->nv*3 );

				rall[3*s+0] = ro[0] + dr[0];
				rall[3*s+1] = ro[1] + dr[1];
				rall[3*s+2] = ro[2] + dr[2];
			}
			else
			{
				theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], rall+3*s, njunk, rsurf );
			}			
				
			rall[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] * alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] * alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] * alphas[0];
			rall[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] * alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] * alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] * alphas[1];
			rall[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] * alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] * alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] * alphas[2];

		}
	}
}

void pcomplex::loadCoords( surface *theSurface, double *rsurf, double *r, double *n )
{
	double *alphas = rsurf+3*theSurface->nv;
	if( bound )
	{
		for( int s = 0; s < nattach; s++ )
		{
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
			
			if( puv[2*s+0] <= 0 || puv[2*s+1] <= 0 || puv[2*s+0]+puv[2*s+1] >= 1 )
			{
				double ro[3];
				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], ro, n+3*s, rsurf );  

				do {
					f_1 = nf;
					nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
				} while( nf != f_1 );
	
				uv1[0] += duv1[0];		
				uv1[1] += duv1[1];		
			
				grad_fs[s] = f_1;
				grad_puv[2*s+0] = uv1[0];
				grad_puv[2*s+1] = uv1[1];

				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], r+3*s, n+3*s, rsurf );  

				double dr[3] = { r[3*s+0] - ro[0], r[3*s+1] - ro[1], r[3*s+2] - ro[2] };
				double del[3];
				MinImage3D( dr, theSurface->PBC_vec, del, rsurf+theSurface->nv*3 );

				r[3*s+0] = ro[0] + dr[0];
				r[3*s+1] = ro[1] + dr[1];
				r[3*s+2] = ro[2] + dr[2];
			}
			else
			{
				theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], r+3*s, n+3*s, rsurf );
				
				grad_fs[s] = fs[s];
				grad_puv[2*s+0] = puv[2*s+0];
				grad_puv[2*s+1] = puv[2*s+1];
			
			}			
				
			r[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] * alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] *alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] * alphas[0];
			r[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] * alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] *alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] * alphas[1];
			r[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] * alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] *alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] * alphas[2];

		}
	}
	else
	{
		memcpy(r , rall, sizeof(double) * nsites );
		memset(n, 0, sizeof(double) * nsites ); 
	}
}

int pcomplex::packLenF( void )
{
	// double positions to include gradient.
	// puv, grad_puv, save_Grad (3*nattach*2)
	// 2 * membrane embedded + 3 * solution
	// 3 * nsites for PBC ext
	// conjugate moment, for now 2 * nattach...

	return   2 * nattach*2 
			+ 3 * nsites * 5 // qdot, p, grad
		; 
}

int pcomplex::packLenI( void )
{
	// just the list of faces.

	return nattach*2;
}

void pcomplex::packageI( int *io )
{
	int io_cntr =0;
	for( int c = 0; c < nattach; c++ )
	{
		io[io_cntr] = fs[c];      io_cntr++;
		io[io_cntr] = grad_fs[c]; io_cntr++;
	}
}

void pcomplex::unpackI( int *io )
{
	int io_cntr = 0;
	for( int c = 0; c < nattach; c++ )
	{
		fs[c]      = io[io_cntr]; io_cntr++;	
		grad_fs[c] = io[io_cntr]; io_cntr++;	
	}
}

void pcomplex::packageF( double *io )
{
	int io_cntr = 0;
	
	for( int c = 0; c < nattach; c++ )
	{
		io[io_cntr] = puv[2*c+0]; io_cntr++;
		io[io_cntr] = puv[2*c+1]; io_cntr++;
		io[io_cntr] = grad_puv[2*c+0]; io_cntr++;
		io[io_cntr] = grad_puv[2*c+1]; io_cntr++;
	}
	
	for( int c = 0; c < nsites; c++ )
	{
		io[io_cntr] = save_grad[3*c+0]; io_cntr++;
		io[io_cntr] = save_grad[3*c+1]; io_cntr++;
		io[io_cntr] = save_grad[3*c+2]; io_cntr++;
	}

	
	for( int c = 0; c < nsites; c++ )
	{
		io[io_cntr] = qdot[3*c+0]; io_cntr++;
		io[io_cntr] = qdot[3*c+1]; io_cntr++;
		io[io_cntr] = qdot[3*c+2]; io_cntr++;
	}
	
	for( int c = 0; c < nsites; c++ )
	{
		io[io_cntr] = p[3*c+0]; io_cntr++;
		io[io_cntr] = p[3*c+1]; io_cntr++;
		io[io_cntr] = p[3*c+2]; io_cntr++;
	}

	// position

	for( int s = 0; s < nsites; s++ )
	{
		io[io_cntr] = rall[3*s+0]; io_cntr++;
		io[io_cntr] = rall[3*s+1]; io_cntr++;
		io[io_cntr] = rall[3*s+2]; io_cntr++;
	}
	
	for( int s = 0; s < nsites; s++ )
	{
		io[io_cntr] = PBC_ext[3*s+0]; io_cntr++;
		io[io_cntr] = PBC_ext[3*s+1]; io_cntr++;
		io[io_cntr] = PBC_ext[3*s+2]; io_cntr++;
	}
	
}

void pcomplex::unpackF( double *io )
{
	int io_cntr = 0;
	
	for( int c = 0; c < nattach; c++ )
	{
		puv[2*c+0]=io[io_cntr]; io_cntr++;
		puv[2*c+1]=io[io_cntr]; io_cntr++;
		grad_puv[2*c+0]=io[io_cntr]; io_cntr++;
		grad_puv[2*c+1]=io[io_cntr]; io_cntr++;
	}
	
	for( int c = 0; c < nsites; c++ )
	{
		save_grad[3*c+0]=io[io_cntr]; io_cntr++;
		save_grad[3*c+1]=io[io_cntr]; io_cntr++;
		save_grad[3*c+2]=io[io_cntr]; io_cntr++;
	}

	// dummy transfer of solution gradient.
	
	for( int c = 0; c < nsites; c++ )
	{
		qdot[3*c+0] = io[io_cntr]; io_cntr++;
		qdot[3*c+1] = io[io_cntr]; io_cntr++;
		qdot[3*c+2] = io[io_cntr]; io_cntr++;
	}
	
	for( int c = 0; c < nsites; c++ )
	{
		p[3*c+0] = io[io_cntr]; io_cntr++;
		p[3*c+1] = io[io_cntr]; io_cntr++;
		p[3*c+2] = io[io_cntr]; io_cntr++;
	}

	// position

	for( int s = 0; s < nsites; s++ )
	{
		rall[3*s+0] = io[io_cntr]; io_cntr++;
		rall[3*s+1] = io[io_cntr]; io_cntr++;
		rall[3*s+2] = io[io_cntr]; io_cntr++;
	}
	
	for( int s = 0; s < nsites; s++ )
	{
		PBC_ext[3*s+0]=io[io_cntr]; io_cntr++;
		PBC_ext[3*s+1]=io[io_cntr]; io_cntr++;
		PBC_ext[3*s+2]=io[io_cntr]; io_cntr++;
	}
}

void propagateSolutionParticles( Simulation *theSimulation, double dt )
{
	pcomplex **allComplexes = theSimulation->allComplexes;
	int ncomplex = theSimulation->ncomplex;
	double ncol = 0;
	static int icntr=0;
	double *alphas = theSimulation->alpha;
	int done = 0;

	for( int c = 0; c < ncomplex; c++ )
	{
		if( allComplexes[c]->disabled ) continue;

		if( allComplexes[c]->do_bd  )
		{
			for( int s = allComplexes[c]->nattach; s < allComplexes[c]->nsites; s++ )
			{ // DC is in seconds, same with dt.  AKMA_TIME converts.
				allComplexes[c]->qdot[3*s+0] = (-allComplexes[c]->save_grad[3*s+0] * allComplexes[c]->DC[s]/kT + sqrt(2 * allComplexes[c]->DC[s]/dt ) * gsl_ran_gaussian(r_gen_global, 1 )) / AKMA_TIME;
				allComplexes[c]->qdot[3*s+1] = (-allComplexes[c]->save_grad[3*s+1] * allComplexes[c]->DC[s]/kT + sqrt(2 * allComplexes[c]->DC[s]/dt ) * gsl_ran_gaussian(r_gen_global, 1 )) / AKMA_TIME;
				allComplexes[c]->qdot[3*s+2] = (-allComplexes[c]->save_grad[3*s+2] * allComplexes[c]->DC[s]/kT + sqrt(2 * allComplexes[c]->DC[s]/dt ) * gsl_ran_gaussian(r_gen_global, 1 )) / AKMA_TIME;
			}
		}
	}
	
	// propagate them all forward, then backtrack as we see fit.
	dt *= AKMA_TIME;
	double propagated = dt;
		
	for( int c = 0; c < ncomplex; c++ )
	{
		if( allComplexes[c]->disabled ) continue;

		for( int s = allComplexes[c]->nattach; s < allComplexes[c]->nsites; s++ )
		{
			allComplexes[c]->rall[3*s+0] += allComplexes[c]->qdot[3*s+0] * dt;
			allComplexes[c]->rall[3*s+1] += allComplexes[c]->qdot[3*s+1] * dt;
			allComplexes[c]->rall[3*s+2] += allComplexes[c]->qdot[3*s+2] * dt;
		}
	}
	double to_propagate = 0;

//	printf("STEP %lf\n", dt );
	while( !done )
	{
		int col1=-1, p1=-1;
		int col2=-1, p2=-1;
			
		ParallelSyncComplexes( allComplexes, ncomplex );

		double toc = timePrecedingElasticCollision( theSimulation, allComplexes, ncomplex, propagated, &col1, &p1, &col2, &p2  );
		
	
		if( col1 >= 0 )
		{
//			printf("Collision between %d and %d, rewinding by %lf.\n", col1, col2, toc );
			// backtrack toc amount.
		
//			for( int cx = 0; cx < par_info.nc; cx++ )
			for( int c = 0; c < ncomplex; c++ )
			{
				if( allComplexes[c]->disabled ) continue;
				for( int s = allComplexes[c]->nattach; s < allComplexes[c]->nsites; s++ )
				{
					allComplexes[c]->rall[3*s+0] -= allComplexes[c]->qdot[3*s+0] * toc;
					allComplexes[c]->rall[3*s+1] -= allComplexes[c]->qdot[3*s+1] * toc;
					allComplexes[c]->rall[3*s+2] -= allComplexes[c]->qdot[3*s+2] * toc;
				}
			}
	
			// collide!!

			int c1 = col1;
			int c2 = col2;

			double r1[3] = { 
				allComplexes[c1]->rall[3*p1+0],
				allComplexes[c1]->rall[3*p1+1],
				allComplexes[c1]->rall[3*p1+2] };
			
			double v1[3] = { 
				allComplexes[c1]->qdot[3*p1+0],
				allComplexes[c1]->qdot[3*p1+1],
				allComplexes[c1]->qdot[3*p1+2] };
		
			double r2[3] = { 
				allComplexes[c2]->rall[3*p2+0],
				allComplexes[c2]->rall[3*p2+1],
				allComplexes[c2]->rall[3*p2+2] };
			
			double v2[3] = { 
				allComplexes[c2]->qdot[3*p2+0],
				allComplexes[c2]->qdot[3*p2+1],
				allComplexes[c2]->qdot[3*p2+2] };
	
			double pos1_at_c[3] = { r1[0],
						r1[1],
						r1[2] };
		
			double pos2_at_c[3] = { r2[0],
						r2[1],
						r2[2] };

			double dv[3] = { v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]};
	
			double m1 = allComplexes[c1]->mass[p1];
			double m2 = allComplexes[c2]->mass[p2];
	
			// direction of dv.
			double dr_at_c[3] = { 
					pos1_at_c[0] - pos2_at_c[0],
					pos1_at_c[1] - pos2_at_c[1],
					pos1_at_c[2] - pos2_at_c[2] };
			double len = length3(dr_at_c);
			if( len > 1.1 * (allComplexes[c1]->sigma[p1] + allComplexes[c2]->sigma[p2]) )
				theSimulation->wrapPBC( dr_at_c, alphas );

			normalize(dr_at_c);
			// solve for dKE
		
			double lambda2 = (2 * m1 * (dr_at_c[0] * dv[0] +dr_at_c[1]*dv[1] + dr_at_c[2] * dv[2])/(m1+m2));
			double lambda1 = -lambda2 * (m2/m1);
	
	
			double v1_after[3] = { v1[0] + lambda1 * dr_at_c[0],
					      v1[1] + lambda1 * dr_at_c[1],			
					      v1[2] + lambda1 * dr_at_c[2] };			
			double v2_after[3] = { v2[0] + lambda2 * dr_at_c[0],
					      v2[1] + lambda2 * dr_at_c[1],			
					      v2[2] + lambda2 * dr_at_c[2] };			
		
			double KE_before = 0.5 * m1 * (v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]) +0.5 * m2 * (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
			double KE_after = 0.5 * m1 * (v1_after[0]*v1_after[0]+v1_after[1]*v1_after[1]+v1_after[2]*v1_after[2]) +0.5 * m2 * (v2_after[0]*v2_after[0]+v2_after[1]*v2_after[1]+v2_after[2]*v2_after[2]);
	
	//		printf("COLLISION %d %d vafter: %le %le %le and %le %le %le KEBA %le/%le\n", c1, c2, 
	//				v1_after[0], v1_after[1], v1_after[2],	
	//				v2_after[0], v2_after[1], v2_after[2], KE_before, KE_after );	
	
			allComplexes[c1]->qdot[3*p1+0] = v1_after[0];
			allComplexes[c1]->qdot[3*p1+1] = v1_after[1];
			allComplexes[c1]->qdot[3*p1+2] = v1_after[2];

			allComplexes[c1]->p[3*p1+0] = allComplexes[c1]->qdot[3*p1+0] * m1;
			allComplexes[c1]->p[3*p1+1] = allComplexes[c1]->qdot[3*p1+1] * m1;
			allComplexes[c1]->p[3*p1+2] = allComplexes[c1]->qdot[3*p1+2] * m1;
			
			allComplexes[c2]->qdot[3*p2+0] = v2_after[0];
			allComplexes[c2]->qdot[3*p2+1] = v2_after[1];
			allComplexes[c2]->qdot[3*p2+2] = v2_after[2];
			
			allComplexes[c2]->p[3*p2+0] = allComplexes[c2]->qdot[3*p2+0] * m2;
			allComplexes[c2]->p[3*p2+1] = allComplexes[c2]->qdot[3*p2+1] * m2;
			allComplexes[c2]->p[3*p2+2] = allComplexes[c2]->qdot[3*p2+2] * m2;
		
	
			for( int c = 0; c < ncomplex; c++ )
			{
				for( int s = allComplexes[c]->nattach; s < allComplexes[c]->nsites; s++ )
				{
					allComplexes[c]->rall[3*s+0] += allComplexes[c]->qdot[3*s+0] * toc;
					allComplexes[c]->rall[3*s+1] += allComplexes[c]->qdot[3*s+1] * toc;
					allComplexes[c]->rall[3*s+2] += allComplexes[c]->qdot[3*s+2] * toc;
				}
			}
	
			ncol++;			
			propagated = toc;
		}
		else
			done = 1;	
	}

	icntr++;

	if( global_debug_mode && ncol > 0 )
		printf("NCollisions: %lf\n", ncol );

//	if( icntr % 100 == 0 )
//		printf("NCollisions/check %le\n", ncol / icntr );
}

void pcomplex::watch( void )
{
	nwatchers++;
}

void pcomplex::forget( void )
{
	nwatchers--;
}

void pcomplex::cache(void)
{
	memcpy( cache_grad, save_grad, sizeof(double) * 3 * nsites );
	memcpy( cache_f, fs, sizeof(int) * nattach );
	memcpy( cache_puv, puv, sizeof(double) * 2*nattach );
	memcpy( cache_rall, rall, sizeof(double) * 3 * nsites );
}

void pcomplex::uncache(void)
{
	memcpy( save_grad, cache_grad, sizeof(double) * 3 * nsites );
	memcpy(        fs, cache_f,    sizeof(int) * nattach );
	memcpy(       puv, cache_puv,  sizeof(double) * 2*nattach );
	memcpy(      rall, cache_rall, sizeof(double) * 3 * nsites );
}


