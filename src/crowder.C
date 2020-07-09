#include "pcomplex.h"
#include <math.h>
#include "mutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// defaults
extern double lipid_DC;
extern double solution_DC;

void elasticCrowder::init( Simulation *theSimulation, surface *theSurface, double *rsurf, int f, double u, double v )
{
	base_init();

	nsites = 2;
	nattach = 1;

	alloc();

	sigma[0] = att_sigma_inner;
	sigma[1] = hard_r;
	
	att_sigma[0] = att_sigma_val;
	att_sigma[1] = 0;
		
	att_eps[0] = att_eps_val;
	att_eps[1] = 0;
	

	bound = 1;

	double rp[3];
	double nrm[3];

	mass[0] = crowder_mass;
	mass[1] = crowder_mass;

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

	rall[3] = rp[0] + d * nrm[0];
	rall[4] = rp[1] + d * nrm[1];
	rall[5] = rp[2] + d * nrm[2];
	
} 

void elasticCrowder::putBonds( int *list )
{
	list[0] = 0;
	list[1]  =1;
}

void elasticCrowder::loadParams( parameterBlock *block )
{
	hard_r = block->crowder_r;	
	d = fabs(block->crowder_d);
	bond_k = fabs(block->crowder_bond_k);
	att_eps_val = block->crowder_attraction;
	att_sigma_val = block->crowder_attraction_r;
	att_sigma_inner = block->crowder_attraction_r_rep;
	crowder_mass = block->crowder_mass;
}

void elasticCrowder::move_inside( void )
{
	d = -fabs(d);
}

void elasticCrowder::move_outside( void )
{
	d = fabs(d);
}

double elasticCrowder::V( Simulation *theSimulation )
{
	double *alphas =theSimulation->alpha;
	double r[6];
	double n[6];

	if( bound )
	{
		// evaluate the real-space coordinates and normals based on the membrane surface coordinates.
		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = theSimulation->fetch(sid[s]);
			surface *theSurface = sRec->theSurface;
			double *rsurf = sRec->r;

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
				MinImage3D( dr, theSurface->PBC_vec, del, rsurf+3*theSurface->nv );

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
				
			r[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] * alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] * alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] * alphas[0];
			r[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] * alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] * alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] * alphas[1];
			r[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] * alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] * alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] * alphas[2];

		}

		// solution particles
		memcpy( r+3, rall+3, sizeof(double) * 3 );
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 6 );
		memset(n, 0, sizeof(double) * 6 ); 
	}
	double dr[3] = { r[3] - (r[0] + d * n[0]),
			  r[4] - (r[1] + d * n[1]),
			  r[5] - (r[2] + d * n[2]) };


	double bond_en = 0.5 * bond_k * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]); 

	double pot = bond_en;

	return pot;
}

double elasticCrowder::grad( Simulation *theSimulation,  double *surfacer_g, double *surfacen_g )
{
	double *alphas = theSimulation->alpha;
	double r[6];
	double n[6];

	if( bound )
	{
		// evaluate the real-space coordinates and normals based on the membrane surface coordinates.
		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = theSimulation->fetch(sid[s]);
			surface *theSurface = sRec->theSurface;
			double *rsurf = sRec->r;
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
			
			coord_transform[4*s+0] = 1;
			coord_transform[4*s+1] = 0;
			coord_transform[4*s+2] = 0;
			coord_transform[4*s+3] = 1;	
			double null_mom[3] = {0,0,0};

			if( puv[2*s+0] <= 0 || puv[2*s+1] <= 0 || puv[2*s+0]+puv[2*s+1] >= 1 )
			{
				double ro[3];
				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], ro, n+3*s, rsurf );  

				do {
					f_1 = nf;
					nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf, null_mom, coord_transform+4*s ); 
				} while( nf != f_1 );
	
				uv1[0] += duv1[0];		
				uv1[1] += duv1[1];		
			
				grad_fs[s] = f_1;
				grad_puv[2*s+0] = uv1[0];
				grad_puv[2*s+1] = uv1[1];

				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], r+3*s, n+3*s, rsurf );  

				double dr[3] = { r[3*s+0] - ro[0], r[3*s+1] - ro[1], r[3*s+2] - ro[2] };
				double del[3];
				MinImage3D( dr, theSurface->PBC_vec, del, rsurf+3*theSurface->nv );

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
				
			r[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] * alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] * alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] * alphas[0];
			r[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] * alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] * alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] * alphas[1];
			r[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] * alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] * alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] * alphas[2];

		}

		// solution particles
		memcpy( r+3, rall+3, sizeof(double) * 3 );
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 6 );
		memset(n, 0, sizeof(double) * 6 ); 
	}
	

	double dr[3] = { r[3] - (r[0] + d * n[0]),
			  r[4] - (r[1] + d * n[1]),
			  r[5] - (r[2] + d * n[2]) };


	double bond_en = 0.5 * bond_k * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]); 
	double pot = bond_en;

	surfacer_g[0] += bond_k * ( - dr[0] );
	surfacer_g[1] += bond_k * ( - dr[1] );
	surfacer_g[2] += bond_k * ( - dr[2] );

	surfacen_g[0] += bond_k * d * (- dr[0] );	
	surfacen_g[1] += bond_k * d * (- dr[1] );	
	surfacen_g[2] += bond_k * d * (- dr[2] );	
	
	save_grad[3] -= bond_k * ( - dr[0] );
	save_grad[4] -= bond_k * ( - dr[1] );
	save_grad[5] -= bond_k * ( - dr[2] );

	return pot;
}
