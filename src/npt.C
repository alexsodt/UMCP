#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interp.h"
#include "pcomplex.h"
#include "parallel.h"
#include "p_p.h"
#include "npt.h"
extern double VC,AVC;
void Simulation::area_MC_move( double beta, int move_type, parameterBlock *block )
{
	static double nacc = 0;
	static double nrej = 0;
	static int icntr = 0;
	static double VMAG = 0.002;
	static double TMAG = 0.002;
	static double N_REXPANSIONS_PREVENTED= 0; // dAlpha > 0
	static double N_RCONTRACTIONS_PREVENTED = 0; // dAlpha < 0
	static double NCONTRACT_COL_FRCD = 0; // dAlpha > 0
	static double NEXPAND_COL_FRCD = 0; // dAlpha < 0

#ifdef PARALLEL
	ParallelSyncComplexes( allComplexes, ncomplex );
#endif
	double E0 = 0;
	double MEME0 = 0;
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		MEME0 += sRec->theSurface->energy(sRec->r,NULL);
	E0 += MEME0;
	ParallelSum(&MEME0,1);		
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];

		double dE = allComplexes[c]->V(this );
		E0 += dE;
		E0 += allComplexes[c]->AttachV( this );
	}

	E0 += Boxed_PP_V( this );
	
	ParallelSum(&E0,1);

	double d_complex = 0;
			
	double npt_KE0 = 0;
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		npt_KE0 += sRec->theSurface->evaluate_T( sRec->qdot, sRec->pp, NULL, NULL, alpha );
	double PT = 0;

	double n_bd_sol_p = 0;

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];
		for( int s = allComplexes[c]->nattach; s < allComplexes[c]->nsites; s++ )
			PT += allComplexes[c]->T(this,s);
		if( allComplexes[c]->do_bd )
			n_bd_sol_p += allComplexes[c]->nsites - allComplexes[c]->nattach;
	}

#ifdef PARALLEL
	ParallelSum( &PT, 1 );
	ParallelSum( &n_bd_sol_p, 1 );
#endif
	double npt_KE0_mesh = npt_KE0;
	double npt_KE0_p    = PT;
	npt_KE0 += PT;

	double ERestraint = 0;
	if( block->alpha_restraint_k )
	{
		if( block->alpha_restraint_x > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( alpha[0] - block->alpha_restraint_x, 2.0); 
		if( block->alpha_restraint_y > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( alpha[1] - block->alpha_restraint_y, 2.0); 
		if( block->alpha_restraint_z > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( alpha[2] - block->alpha_restraint_z, 2.0); 
	}

	E0 += ERestraint;
	
#ifdef PARALLEL
	MPI_Bcast( &E0, 1, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif
	
	double nElastic0 = nElasticCollisions( this, allComplexes, ncomplex  ); 

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
	
	alpha[0] *= scale_fac[0];
	alpha[1] *= scale_fac[1];
	alpha[2] *= scale_fac[2];
	
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		sRec->r[3*sRec->theSurface->nv+0] = alpha[0];
		sRec->r[3*sRec->theSurface->nv+1] = alpha[1];
		sRec->r[3*sRec->theSurface->nv+2] = alpha[2];
	}

//	for( int cx = 0; cx < par_info.nc; cx++ )
//	{

	for( int c = 0; c < ncomplex; c++ )
	{
//	int c = par_info.complexes[cx];
		for( int p = allComplexes[c]->nattach; p < allComplexes[c]->nsites; p++ )
		{
			allComplexes[c]->rall[3*p+0] *= scale_fac[0];
			allComplexes[c]->rall[3*p+1] *= scale_fac[1];
			allComplexes[c]->rall[3*p+2] *= scale_fac[2];
			
			allComplexes[c]->p[3*p+0] /= scale_fac[0];
			allComplexes[c]->p[3*p+1] /= scale_fac[1];
			allComplexes[c]->p[3*p+2] /= scale_fac[2];
			
			allComplexes[c]->qdot[3*p+0] /= scale_fac[0];
			allComplexes[c]->qdot[3*p+1] /= scale_fac[1];
			allComplexes[c]->qdot[3*p+2] /= scale_fac[2];
		}
	}
	double E1=0;
	double MEME1 = 0;
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		MEME1+= sRec->theSurface->energy(sRec->r,NULL);
	E1 += MEME1;
	ParallelSum(&MEME1,1);		

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c = par_info.complexes[cx];

		double dE = allComplexes[c]->V(this );
		E1 += dE;
		E1 += allComplexes[c]->AttachV( this );
	}

	E1 += Boxed_PP_V( this );
	
	ParallelSum(&E1,1);



	// reciprocal relationship for scaling factor through effective mass matrix which is inverse of position-position dot product.
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	for( int x = 0; x < sRec->theSurface->nv; x++ )
	{
		sRec->qdot[3*x+0] /= scale_fac[0]*scale_fac[0];
		sRec->qdot[3*x+1] /= scale_fac[1]*scale_fac[1];
		sRec->qdot[3*x+2] /= scale_fac[2]*scale_fac[2];
	}
	
	double npt_KE1 = 0;
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		npt_KE1 += sRec->theSurface->evaluate_T( sRec->qdot, sRec->pp, NULL, NULL, alpha );
	PT = 0;
	for( int cx = 0; cx < par_info.nc; cx++ )
	{

		int c = par_info.complexes[cx];
		for( int s = allComplexes[c]->nattach; s < allComplexes[c]->nsites; s++ )
			PT += allComplexes[c]->T(this,s);
	}

#ifdef PARALLEL
	ParallelSum( &PT, 1 );
#endif
	double npt_KE1_mesh = npt_KE1;
	double npt_KE1_p    = PT;
	npt_KE1 += PT;

	ERestraint = 0;
	if( block->alpha_restraint_k )
	{
		if( block->alpha_restraint_x > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( alpha[0] - block->alpha_restraint_x, 2.0); 
		if( block->alpha_restraint_y > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( alpha[1] - block->alpha_restraint_y, 2.0); 
		if( block->alpha_restraint_z > 0 ) ERestraint += 0.5 * block->alpha_restraint_k * pow( alpha[2] - block->alpha_restraint_z, 2.0); 
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

//	printf("dE: %.14le dKE: %.14le dAlpha: %.14le cur_alpha_x: %le\n", E1-E0, npt_KE1-npt_KE0, dAlpha, r[3*nv+0]*dAlpha );

	double d_virtual_KE_BD = 0;

	// each cartesian KE has distribution exp(-m v^2 /2 kT)

	d_virtual_KE_BD += n_bd_sol_p * (-1.0/beta)*log(scale_fac[0]); 
	d_virtual_KE_BD += n_bd_sol_p * (-1.0/beta)*log(scale_fac[1]); 
	d_virtual_KE_BD += n_bd_sol_p * (-1.0/beta)*log(scale_fac[2]); 

	double pr = exp(-beta*( (E1-E0) + (npt_KE1-npt_KE0) + d_virtual_KE_BD ));
	if( nElastic1 > nElastic0 + 0.5)
	{
		if( dAlpha > 1 ) // contraction prevented ... R -> R / dalpha, prevented ... cylinder ``expands''
			N_RCONTRACTIONS_PREVENTED += 1;
		else // expand, prevented
			N_REXPANSIONS_PREVENTED += 1; // expansion prevented ... R -> R / dalpha, prevented ... cylinder ``collapses''
		pr = 0;
//		printf("NELASTIC1: %lf NELASTIC0: %lf pr: %le\n", nElastic1, nElastic0, pr );
	}
	else if( nElastic0 > nElastic1 + 0.5)
	{
		pr = 1;
		if( dAlpha > 1 ) // contract, forced
			NEXPAND_COL_FRCD += 1;
		else // expand, forced
			NCONTRACT_COL_FRCD += 1;
//		printf("NELASTIC1: %lf NELASTIC0: %lf pr: %le\n", nElastic1, nElastic0, pr );
	}

//	printf("NELASTIC1: %lf NELASTIC0: %lf pr: %le dalpha: %.14le dE: %.14le dMEM: %.14le cur: %.14le\n", nElastic1, nElastic0, pr, dAlpha, E1-E0, MEME1-MEME0, r[3*nv+0]/scale_fac[0] );

	int decision = 0;
	if( rand() / (double)RAND_MAX < pr )
		decision = 1;
#ifdef PARALLEL
	ParallelBroadcast(&decision,1);	
#endif
	if( decision )
	{
#if 0
		printf("ACCEPTED dKE_mesh: %le dKE_p: %le dV_mesh: %le dV_p: %le dAlpha: %le\n",
			npt_KE1_mesh-npt_KE0_mesh,	
			npt_KE1_p - npt_KE0_p + d_virtual_KE_BD,
			MEME1-MEME0,
			E1-E0-(MEME1-MEME0),  dAlpha );		
#endif
		nacc++;
	}
	else
	{
		alpha[0] /= scale_fac[0];
		alpha[1] /= scale_fac[1];
		alpha[2] /= scale_fac[2];
	
		for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		{
			sRec->r[3*sRec->theSurface->nv+0] = alpha[0];
			sRec->r[3*sRec->theSurface->nv+1] = alpha[1];
			sRec->r[3*sRec->theSurface->nv+2] = alpha[2];
		}

		for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		{
			for( int x = 0; x < sRec->theSurface->nv; x++ )
			{
				sRec->qdot[3*x+0] *= scale_fac[0]*scale_fac[0];
				sRec->qdot[3*x+1] *= scale_fac[1]*scale_fac[1];
				sRec->qdot[3*x+2] *= scale_fac[2]*scale_fac[2];
			}
		}

		for( int c = 0; c < ncomplex; c++ )
		{
//			int c = par_info.complexes[cx];
			for( int p = allComplexes[c]->nattach; p < allComplexes[c]->nsites; p++ )
			{
				allComplexes[c]->rall[3*p+0] /= scale_fac[0];
				allComplexes[c]->rall[3*p+1] /= scale_fac[1];
				allComplexes[c]->rall[3*p+2] /= scale_fac[2];
			
				allComplexes[c]->p[3*p+0] *= scale_fac[0];
				allComplexes[c]->p[3*p+1] *= scale_fac[1];
				allComplexes[c]->p[3*p+2] *= scale_fac[2];
				
				allComplexes[c]->qdot[3*p+0] *= scale_fac[0];
				allComplexes[c]->qdot[3*p+1] *= scale_fac[1];
				allComplexes[c]->qdot[3*p+2] *= scale_fac[2];
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

//		area( r, -1, &area0, &o_area );

		printf("N_RADIAL_EXPANSIONS_PREVENTED %lf N_RADIAL_CONTRACTIONS_PREVENTED %lf NEXP,FRCD %lf NCOL,FRCD %lf\n",
			N_REXPANSIONS_PREVENTED,N_RCONTRACTIONS_PREVENTED, NEXPAND_COL_FRCD, NCONTRACT_COL_FRCD );
//		printf("FR: %le VMAG: %le alpha: %le area %le\n", fr, VMAG, r[3*nv+0], area0 );
		nacc = 0;
		nrej = 0;
		icntr = 0;
	}
#ifdef PARALLEL
	ParallelBroadcast(alpha,3);	
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		sRec->r[3*sRec->theSurface->nv+0] = alpha[0];
		sRec->r[3*sRec->theSurface->nv+1] = alpha[1];
		sRec->r[3*sRec->theSurface->nv+2] = alpha[2];
	}
#endif

//	if( block->alpha_restraint_k > 0 )
//	if( move_type == CYL_TENSION_MOVE )
//		printf("ALPHAS: %.14le %.14le %.14le\n", r[3*nv], r[3*nv+1], r[3*nv+2] );

}
