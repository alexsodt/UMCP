#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interp.h"
#include "pcomplex.h"
#include "parallel.h"
#include "p_p.h"
#include "npt.h"
extern double VC,AVC;
void surface::area_MC_move( double *r, pcomplex **allComplexes, int ncomplex, double beta, double *qdot, double *pp, int move_type, parameterBlock *block )
{
	static double nacc = 0;
	static double nrej = 0;
	static int icntr = 0;
	static double VMAG = 0.0003;
	static double TMAG = 0.0003;
	
#ifdef PARALLEL
	ParallelSyncComplexes( allComplexes, ncomplex );
#endif
	double E0 = energy( r, NULL );

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];

		E0 += allComplexes[c]->V(this, r );
		E0 += allComplexes[c]->AttachV( this, r );
	}

	E0 += PP_V( this, r, allComplexes, ncomplex );
	
	ParallelSum(&E0,1);
					
	double npt_KE0 = evaluate_T( qdot, pp, NULL, NULL, r+3*nv );
	double PT = 0;
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		PT += allComplexes[c]->T();
	}

#ifdef PARALLEL
	ParallelSum( &PT, 1 );
#endif
	npt_KE0 += PT;

	double ERestraint = 0;
	if( block->alpha_restraint_k )
	{
		if( block->alpha_restraint_x > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( log(r[3*nv+0]) - block->alpha_restraint_x, 2.0); 
		if( block->alpha_restraint_y > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( log(r[3*nv+1]) - block->alpha_restraint_y, 2.0); 
		if( block->alpha_restraint_z > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( log(r[3*nv+2]) - block->alpha_restraint_z, 2.0); 
	}

	E0 += ERestraint;
	
#ifdef PARALLEL
	MPI_Bcast( &E0, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif
	
	double nElastic0 = nElasticCollisions( this, allComplexes, ncomplex ); 

#ifdef PARALLEL
	ParallelSum(&nElastic0,1);
#endif
	double rn = 2 * ((double)rand() / (double)RAND_MAX-0.5);

	if( move_type == VOLUME_MOVE )
		rn *= VMAG;
	else if( move_type == CYL_TENSION_MOVE )
		rn *= TMAG;

	double dAlpha = exp( - rn );

#ifdef PARALLEL
	MPI_Bcast( &dAlpha, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif

//#define T_DEBUG
#ifdef T_DEBUG
	int nalpha = 100;
	for( int ix = 0; ix < nalpha; ix++ )
	{
		r[3*nv+0] *= exp(0.001);
		r[3*nv+1] *= exp(0.001);
		r[3*nv+2] *= exp(-0.001);

		double area0;
		double cur_area;
		area(r, -1, &cur_area, &area0 );
		//printf("area: %le area0: %le\n", cur_area, area0 );
		double E = energy( r, NULL );

		printf("%le %.14le VC: %.14le AVC: %le area: %le\n", r[3*nv+0], E, VC, AVC, cur_area ); 
	}
	exit(1);
#endif
	double scale_fac[3]={1,1,1};	
	if( move_type == VOLUME_MOVE )
	{
		scale_fac[0] = dAlpha;
		scale_fac[1] = dAlpha;
		scale_fac[2] = dAlpha;

	}
	else if( move_type == CYL_TENSION_MOVE )
	{
		scale_fac[0] = 1.0/dAlpha;
		scale_fac[1] = 1.0/dAlpha;
		scale_fac[2] = dAlpha;
	}
	
	r[3*nv+0] *= scale_fac[0];
	r[3*nv+1] *= scale_fac[1];
	r[3*nv+2] *= scale_fac[2];

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		for( int p = allComplexes[c]->nattach; p < allComplexes[c]->nsites; p++ )
		{
			allComplexes[c]->rall[3*p+0] *= scale_fac[0];
			allComplexes[c]->rall[3*p+1] *= scale_fac[1];
			allComplexes[c]->rall[3*p+2] *= scale_fac[2];
		}
	}
	double E1 = energy( r, NULL );

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];

		E1 += allComplexes[c]->V(this, r );
		E1 += allComplexes[c]->AttachV( this, r );
	}

	E1 += PP_V( this, r, allComplexes, ncomplex );
	
	ParallelSum(&E1,1);
	
	double npt_KE1 = evaluate_T( qdot, pp, NULL, NULL, r+3*nv );
	PT = 0;
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		PT += allComplexes[c]->T();
	}

#ifdef PARALLEL
	ParallelSum( &PT, 1 );
#endif
	npt_KE1 += PT;

	ERestraint = 0;
	if( block->alpha_restraint_k )
	{
		if( block->alpha_restraint_x > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( log(r[3*nv+0]) - block->alpha_restraint_x, 2.0); 
		if( block->alpha_restraint_y > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( log(r[3*nv+1]) - block->alpha_restraint_y, 2.0); 
		if( block->alpha_restraint_z > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( log(r[3*nv+2]) - block->alpha_restraint_z, 2.0); 
	}
	E1 += ERestraint;
	
	// probably should ignore particles that are colliding before move, or both before and after -- those particles wouldn't be colliding anyway.
	double nElastic1 = nElasticCollisions( this, allComplexes, ncomplex ); 

#ifdef PARALLEL
	ParallelSum(&nElastic1,1);
#endif

	// Change of kinetic energy due to scaling.

	//printf("alpha0: %le E0: %le alpha1: %le E1: %le, dE: %.14le\n",
//		r[3*nv+0]/dAlpha, E0, r[3*nv+0], E1, E1-E0 );
	

//	printf("dAlpha: %le KE1: %le KE0: %le dKE: %le\n", dAlpha-1, npt_KE1, npt_KE0, npt_KE1-npt_KE0 );

	printf("dE: %.14le dKE: %.14le dAlpha: %.14le cur_alpha_x: %le\n", E1-E0, npt_KE1-npt_KE0, dAlpha, r[3*nv+0]*dAlpha );

	double pr = exp(-beta*( (E1-E0) + (npt_KE1-npt_KE0) ));

	if( nElastic1 > nElastic0 + 0.5)
		pr = 0;
	else if( nElastic0 > nElastic1 + 0.5)
		pr = 1;

	if( rand()/(double)RAND_MAX < pr )
	{
		nacc++;
	}
	else
	{
		r[3*nv+0] /= scale_fac[0];
		r[3*nv+1] /= scale_fac[1];
		r[3*nv+2] /= scale_fac[2];
	
		for( int cx = 0; cx < par_info.nc; cx++ )
		{
			int c = par_info.complexes[cx];
			for( int p = allComplexes[c]->nattach; p < allComplexes[c]->nsites; p++ )
			{
				allComplexes[c]->rall[3*p+0] /= scale_fac[0];
				allComplexes[c]->rall[3*p+1] /= scale_fac[1];
				allComplexes[c]->rall[3*p+2] /= scale_fac[2];
			}
		}
		nrej++;
	}
	
	icntr++;

	if( icntr == 100 )
	{
		double fr = (nacc)/(nacc+nrej);

		if( fr < 0.3 )		
		{
			if( move_type == VOLUME_MOVE )
				VMAG *= 0.8;	
			else
				TMAG *= 0.8;
		}
		else if( fr > 0.7 )
		{
			if( move_type == VOLUME_MOVE )
				VMAG /= 0.8;	
			else
				TMAG /= 0.8;
		}

		double area0, o_area;

		area( r, -1, &area0, &o_area );

		printf("FR: %le VMAG: %le alpha: %le area %le\n", fr, VMAG, r[3*nv+0], area0 );
		nacc = 0;
		nrej = 0;
		icntr = 0;
	}
#ifdef PARALLEL
	ParallelBroadcast(r+3*nv,3);	
#endif

//	if( block->alpha_restraint_k > 0 )
	if( move_type == CYL_TENSION_MOVE )
		printf("ALPHAS: %.14le %.14le %.14le\n", r[3*nv], r[3*nv+1], r[3*nv+2] );

}
