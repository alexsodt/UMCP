#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mutil.h"
#include "gsl_random_globals.h"
#include "alignSet.h"
#define SQRT_SAFETY (1e-7)

#define EPS_SMALL (1e-14)
#define WORKING

//SEKVGGAEGTKLDDDFKEMERKVDVTSRAVMEIMTKTIEYLQPNPASRAKLSMINTMSKIRGQEKGPGYPQAEALLAEAMLKFGRELGDDCNFGPALGEVGEAMRELSEVKDSLDIEVKQNFIDPLQNLHDKDLREIQHHLKKLEGRRLDFDYKKKRQGKIPDEELRQALEKFDESKEIAESSMFNLLEMDIEQVSQLSALVQAQLEYHKQAVQILQQVTVRLEERIRQASS
static double monomer_MW = 10*26492.14; //amu 
static double NBAR_P_RADIUS = 25.0;

// initialize the BAR domain in solution.

void NBAR::init( double *r )
{
	base_init();

	nsites = 3;
	nattach = 3;

	alloc();

	sigma[0] = NBAR_P_RADIUS;
	sigma[1] = NBAR_P_RADIUS;
	sigma[2] = NBAR_P_RADIUS;

	mass[0] = monomer_MW * 2 / 3.0;
	mass[1] = monomer_MW * 2 / 3.0;
	mass[2] = monomer_MW * 2 / 3.0;
	bound = 0;

	rall[3] = r[0];
	rall[4] = r[1];
	rall[5] = r[2];

	// random orientation.

	check_random_init();
		
	// draw random orientation
	double dx,dy,dz;
	gsl_ran_dir_3d( r_gen_global, &dx,&dy,&dz); 
	double dx2,dy2,dz2;
	gsl_ran_dir_3d( r_gen_global, &dx2,&dy2,&dz2); 
	
	double axis[3] ={ dx, dy, dz };
	double rand[3] ={ dx2, dy2, dz2 };

	rall[0] = rall[3] + dx * bond_length;
	rall[1] = rall[4] + dy * bond_length;
	rall[2] = rall[5] + dz * bond_length;
	
	rall[6] = rall[3] + dx * bond_length;
	rall[7] = rall[4] + dy * bond_length;
	rall[8] = rall[5] + dz * bond_length;

	double perp[3];

	cross( axis, rand, perp );
	normalize(perp);

	rotateArbitrary( rall, perp, rall+3, 1, theta_0/2 ); 
	rotateArbitrary( rall+6, perp, rall+3, 1, -theta_0/2 ); 


	memcpy( last_pos, rall, sizeof(double) * 9 );
	memset( PBC_ext, 0, sizeof(double) * 9 );
}

// initialize the BAR domain on the membrane.

void NBAR::init( surface *theSurface, double *rsurf, int f, double u, double v )
{
	base_init();

	nsites = 3;
	nattach = 3;

	alloc();

	bound = 1;
	mass[0] = monomer_MW * 2 / 3.0;
	mass[1] = monomer_MW * 2 / 3.0;
	mass[2] = monomer_MW * 2 / 3.0;
	
	sigma[0] = NBAR_P_RADIUS;
	sigma[1] = NBAR_P_RADIUS;
	sigma[2] = NBAR_P_RADIUS;

	// place the NBAR here, randomly oriented.

	double rp[3];
	double nrm[3];
	theSurface->evaluateRNRM( f, u, v, rp, nrm, rsurf);
	
	rall[3] = rp[0];
	rall[4] = rp[1];
	rall[5] = rp[2];

	puv[2] = u;
	puv[3] = v;
	fs[1] = f;
	
	double dr_u[3];
	theSurface->ru( f, u, v, rsurf, dr_u );
	double drdu = normalize(dr_u);
	double dr_v[3];
	theSurface->rv( f, u, v, rsurf, dr_v );
	double drdv = normalize(dr_v);

	double druv[2] = {1,0};

	double expec[3] = { 
			druv[0] * drdu * dr_u[0] + druv[1] * drdv * dr_v[0],
			druv[0] * drdu * dr_u[1] + druv[1] * drdv * dr_v[1],
			druv[0] * drdu * dr_u[2] + druv[1] * drdv * dr_v[2] };
	double len = normalize(expec);

	// start on our face at the center, choose directions for ends of BAR.
	int eoff[2] = {0,2};
	int dirs[2] = { 1, -1 };
	for( int ends = 0; ends < 2; ends++ )
	{
		int f_1 = f, nf = f;
		double uv1[2] = { u, v };
		double duv1[2] = { dirs[ends] * bond_length * druv[0] / len, dirs[ends] * bond_length * druv[1] / len };
	
		do {
			f_1 = nf;
			nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
		} while( nf != f_1 );

		uv1[0] += duv1[0];		
		uv1[1] += duv1[1];		

		puv[0+2*eoff[ends]]   = uv1[0];
		puv[1+2*eoff[ends]]   = uv1[1];
		fs[eoff[ends]]    = f_1;	

		double rp[3];
		double nrm[3];
		theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rp, nrm, rsurf);
	
		rall[0+3*eoff[ends]] = rp[0];	
		rall[1+3*eoff[ends]] = rp[1];	
		rall[2+3*eoff[ends]] = rp[2];	
	}	

	memcpy( grad_fs, fs, sizeof(int) * 3 );
	memcpy( grad_puv, puv, sizeof(double) * 6 );
	
	double dr1[3] = { rall[0] - rall[3], rall[1] - rall[4], rall[2] - rall[5] };
	theSurface->wrapPBC(dr1,rsurf+3*theSurface->nv);
	rall[0] = rall[3] + dr1[0];
	rall[1] = rall[4] + dr1[1];
	rall[2] = rall[5] + dr1[2];
	double dr2[3] = { rall[6] - rall[3], rall[7] - rall[4], rall[8] - rall[5] };
	theSurface->wrapPBC(dr2,rsurf+3*theSurface->nv);
	rall[6] = rall[3] + dr2[0];
	rall[7] = rall[4] + dr2[1];
	rall[8] = rall[5] + dr2[2];

	memcpy( last_pos, rall, sizeof(double) * 9 );
	memset( PBC_ext, 0, sizeof(double) * 9 );
}

// custom orient procedure.

void NBAR::orient( surface *theSurface, double *rsurf )
{
	/* 
		For flat sections of the membrane (presumably, a planar surface):
		orient the BAR so it points in the major PBC direction 

		For surfaces with anisotropic curvature, orient the BAR domain along max positive curvature.
	*/	

}


void NBAR::bind( int f, double u, double v)
{
	bound = 1;
}

void NBAR::unbind( void )
{
	bound = 0;
}


void NBAR::loadParams( parameterBlock *block )
{
	bond_length = block->bar_bond_length;
	k_phi = block->bar_phi_k;
	k_theta = block->bar_theta_k;
	bond_k = block->bar_bond_k;
	theta_0 = block->bar_theta_0;
	phi_0 = block->bar_phi_0;		
}

int NBAR::getNBonds( void )
{
	return 2;
}

void NBAR::putBonds( int *bond_list )
{
	bond_list[0] = 0;
	bond_list[1] = 1;

	bond_list[2] = 1;
	bond_list[3] = 2;
}

double NBAR::V( surface *theSurface, double *rsurf )
{
	double *alphas = rsurf+3*theSurface->nv;
	double pot = 0;

	double r[9];
	double n[9];

	if( bound )
	{
		// evaluate the real-space coordinates and normals based on the membrane surface coordinates.
#if 0
		for( int s = 0; s < nattach; s++ )
		{
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
	
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
		}
#else		
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
				
			r[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] *alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] *alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] *alphas[0];
			r[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] *alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] *alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] *alphas[1];
			r[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] *alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] *alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] *alphas[2];

		}
#endif
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 9 );
		memset(n, 0, sizeof(double) * 9 ); 
	}
	double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	
	double dr2[3] = { r[6] - r[3], r[7] - r[4], r[8] - r[5] };	

	//lengths of bonds and normal	
	double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));
	double dr2_length = sqrt((dr2[0]*dr2[0]) + (dr2[1]*dr2[1]) + (dr2[2]*dr2[2]));

	double bond_en = 0.5 * bond_k * ( 
			(dr1_length-bond_length)*(dr1_length-bond_length) + 
			(dr2_length-bond_length)*(dr2_length-bond_length) );

	pot += bond_en;

	//Calculates dot product between dr1 and dr2
	double dp_dr1_dr2 = dr1[0]*dr2[0] + dr1[1]*dr2[1] + dr1[2]*dr2[2];


	if( bound )
	{	
		//Calculates dot product between dr and normal
		double nrm_length = sqrt((n[3]*n[3]) + (n[4]*n[4]) + (n[5]*n[5]));
		double dp_dr1_nrm = dr1[0]*n[3] + dr1[1]*n[4] + dr1[2]*n[5];
		double dp_dr2_nrm = dr2[0]*n[3] + dr2[1]*n[4] + dr2[2]*n[5];

	   //Angles between dr1, dr2 and normal
		double theta_dr1_rad = acos((1-SQRT_SAFETY)*dp_dr1_nrm/(dr1_length*nrm_length));
		double theta_dr1_deg = theta_dr1_rad*180/M_PI;
		double theta_dr2_rad = acos((1-SQRT_SAFETY)*dp_dr2_nrm/(dr2_length*nrm_length));
	    double theta_dr2_deg = theta_dr2_rad*180/M_PI;
	
		//arrays hold xyx positions for tangent vectors
		double dr1_t [3] = {dr1[0] - dp_dr1_nrm * n[3], 
			dr1[1] - dp_dr1_nrm * n[4], 
			dr1[2] - dp_dr1_nrm * n[5]};
	
		double dr2_t[3] = {dr2[0] - dp_dr2_nrm * n[3],
			dr2[1] - dp_dr2_nrm * n[4],
			dr2[2] - dp_dr2_nrm * n[5]};
	
		//Calculates angle between dr1 and dr2
		double theta_rad = acos( (1-SQRT_SAFETY)*(dp_dr1_dr2)/(dr1_length*dr2_length));
		double theta_deg = theta_rad*180.0/M_PI; //angle in degrees
	
		
		//Dot product between tangent vectors
		double dp_dr1t_dr2t = dr1_t[0]*dr2_t[0] + dr1_t[1]*dr2_t[1] + dr1_t[2]*dr2_t[2];
	
		//lengths of tangent vectors
		double dr1_t_length = sqrt(dr1_t[0]*dr1_t[0] + dr1_t[1]*dr1_t[1] + dr1_t[2]*dr1_t[2]);
		double dr2_t_length = sqrt(dr2_t[0]*dr2_t[0] + dr2_t[1]*dr2_t[1] + dr2_t[2]*dr2_t[2]);
	
		//Calculates angle between dr1_t and dr2_t
		double phi_rad = acos((1-SQRT_SAFETY)*(dp_dr1t_dr2t)/(dr1_t_length*dr2_t_length));
		double phi_deg = phi_rad*180.0/M_PI;
	
		double theta_dr1_en = 0.5 *  k_theta * ( (theta_dr1_deg - theta_0)*(theta_dr1_deg - theta_0) ); //add theta energy to bond energy
		double theta_dr2_en = 0.5 *  k_theta * ( (theta_dr2_deg - theta_0)*(theta_dr2_deg - theta_0) );
	
		double phi_en = 0.5 * k_phi * (180.0/M_PI) * (180.0/M_PI) * pow( (dp_dr1t_dr2t)/(dr1_t_length*dr2_t_length)+1, 2.0 );

		pot += theta_dr1_en + theta_dr2_en + phi_en;
	}
	else
	{
		double theta = acos( dp_dr1_dr2 );

		double theta_en = 0.5 * k_theta * (theta - theta_0) * (theta - theta_0);

		pot += theta_en;
	}

	return pot;
}

// gets derivative of internal energy relative to position (surfacer_g) and the normal (surfacen_g).

double NBAR::grad(surface *theSurface, double *rsurf, double *surfacer_g, double *surfacen_g )
{
	double *alphas = rsurf+3*theSurface->nv;
	double pot = 0;

	double r[9];
	double n[9];

	if( bound )
	{
#if OLD_METHOD
		for( int s = 0; s < nattach; s++ )
		{
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
	
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
		}
#else
		for( int s = 0; s < nattach; s++ )
		{
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
			
			double null_mom[2] = {0,0};

			coord_transform[4*s+0] = 1;
			coord_transform[4*s+1] = 0;
			coord_transform[4*s+2] = 0;
			coord_transform[4*s+3] = 1;

			if( puv[2*s+0] <= 0 || puv[2*s+1] <= 0 || puv[2*s+0]+puv[2*s+1] >= 1 )
			{
				double ro[3];
				theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], ro, n+3*s, rsurf );  

				do {
					f_1 = nf;
					nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf, null_mom, coord_transform+4*s  ); 
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
				grad_fs[s] = fs[s];
				grad_puv[2*s+0] = puv[2*s+0];
				grad_puv[2*s+1] = puv[2*s+1];
				theSurface->evaluateRNRM( fs[s], puv[2*s+0], puv[2*s+1], r+3*s, n+3*s, rsurf );
	
			}			
				
			r[3*s+0] += theSurface->PBC_vec[0][0] * PBC_ext[3*s+0] * alphas[0] + theSurface->PBC_vec[1][0] * PBC_ext[3*s+1] * alphas[0] + theSurface->PBC_vec[2][0] * PBC_ext[3*s+2] * alphas[0];
			r[3*s+1] += theSurface->PBC_vec[0][1] * PBC_ext[3*s+0] * alphas[1] + theSurface->PBC_vec[1][1] * PBC_ext[3*s+1] * alphas[1] + theSurface->PBC_vec[2][1] * PBC_ext[3*s+2] * alphas[1];
			r[3*s+2] += theSurface->PBC_vec[0][2] * PBC_ext[3*s+0] * alphas[2] + theSurface->PBC_vec[1][2] * PBC_ext[3*s+1] * alphas[2] + theSurface->PBC_vec[2][2] * PBC_ext[3*s+2] * alphas[2];

		}
#endif		
		double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	
		double dr2[3] = { r[6] - r[3], r[7] - r[4], r[8] - r[5] };	
	
			// nx ny nz n[3] n[4] n[5]
		
		//lengths of bonds and normal	
		double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));
		double dr2_length = sqrt((dr2[0]*dr2[0]) + (dr2[1]*dr2[1]) + (dr2[2]*dr2[2]));
		double nrm_length = sqrt((n[3]*n[3]) + (n[4]*n[4]) + (n[5]*n[5]));
	
		double bond_en = 0.5 * bond_k * ( 
				(dr1_length-bond_length)*(dr1_length-bond_length) + 
				(dr2_length-bond_length)*(dr2_length-bond_length) );
		//Calculates dot product between dr1 and dr2
		double dp_dr1_dr2 = dr1[0]*dr2[0] + dr1[1]*dr2[1] + dr1[2]*dr2[2];
		
		//Calculates dot product between dr and normal
		double dp_dr1_nrm = dr1[0]*n[3] + dr1[1]*n[4] + dr1[2]*n[5];
		double dp_dr2_nrm = dr2[0]*n[3] + dr2[1]*n[4] + dr2[2]*n[5];
	
	   //Angles between dr1, dr2 and normal
		double theta_dr1_rad = acos((1-SQRT_SAFETY)*dp_dr1_nrm/(dr1_length*nrm_length));
		double theta_dr1_deg = theta_dr1_rad*180/M_PI;
		double theta_dr2_rad = acos((1-SQRT_SAFETY)*dp_dr2_nrm/(dr2_length*nrm_length));
	    double theta_dr2_deg = theta_dr2_rad*180/M_PI;
	
		//arrays hold xyx positions for tangent vectors
		double dr1_t [3] = {dr1[0] - dp_dr1_nrm * n[3], 
			dr1[1] - dp_dr1_nrm * n[4], 
			dr1[2] - dp_dr1_nrm * n[5]};
	
		double dr2_t[3] = {dr2[0] - dp_dr2_nrm * n[3],
			dr2[1] - dp_dr2_nrm * n[4],
			dr2[2] - dp_dr2_nrm * n[5]};
	
		//Calculates angle between dr1 and dr2
		double theta_rad = acos( (1-SQRT_SAFETY)*(dp_dr1_dr2)/(dr1_length*dr2_length));
		double theta_deg = theta_rad*180.0/M_PI; //angle in degrees
		
		//Dot product between tangent vectors
		double dp_dr1t_dr2t = dr1_t[0]*dr2_t[0] + dr1_t[1]*dr2_t[1] + dr1_t[2]*dr2_t[2];
	
		//lengths of tangent vectors
		double dr1_t_length = sqrt(dr1_t[0]*dr1_t[0] + dr1_t[1]*dr1_t[1] + dr1_t[2]*dr1_t[2]);
		double dr2_t_length = sqrt(dr2_t[0]*dr2_t[0] + dr2_t[1]*dr2_t[1] + dr2_t[2]*dr2_t[2]);
	
		//Calculates angle between dr1_t and dr2_t
		double phi_rad = acos((1-SQRT_SAFETY)*(dp_dr1t_dr2t)/(dr1_t_length*dr2_t_length));
		double phi_deg = phi_rad*180.0/M_PI;
	
		// FIRST part of gradient, theta.
	
		// the derivative of the energy with respect to theta_dr1_deg.
	
		double d_e_d_theta_dr1 = k_theta * ( theta_dr1_deg - theta_0);
		double d_e_d_theta_dr2 = k_theta * ( theta_dr2_deg - theta_0);
	
		double d_theta_dr1_deg_d_r1x = (180.0/M_PI) * (-n[3] * dr1_length*dr1_length + dp_dr1_nrm * dr1[0] ) / ( nrm_length * sqrt(  EPS_SMALL + 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r1y = (180.0/M_PI) * (-n[4] * dr1_length*dr1_length + dp_dr1_nrm * dr1[1] ) / ( nrm_length * sqrt(  EPS_SMALL + 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r1z = (180.0/M_PI) * (-n[5] * dr1_length*dr1_length + dp_dr1_nrm * dr1[2] ) / ( nrm_length * sqrt(  EPS_SMALL + 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		
		double d_theta_dr1_deg_d_r2x = -(180.0/M_PI) * (-n[3] * dr1_length*dr1_length + dp_dr1_nrm * dr1[0] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r2y = -(180.0/M_PI) * (-n[4] * dr1_length*dr1_length + dp_dr1_nrm * dr1[1] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r2z = -(180.0/M_PI) * (-n[5] * dr1_length*dr1_length + dp_dr1_nrm * dr1[2] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
	
		double d_theta_dr2_deg_d_r3x = (180.0/M_PI) * (-n[3] * dr2_length*dr2_length + dp_dr2_nrm * dr2[0] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length ); 
		double d_theta_dr2_deg_d_r3y = (180.0/M_PI) * (-n[4] * dr2_length*dr2_length + dp_dr2_nrm * dr2[1] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		double d_theta_dr2_deg_d_r3z = (180.0/M_PI) * (-n[5] * dr2_length*dr2_length + dp_dr2_nrm * dr2[2] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		
		double d_theta_dr2_deg_d_r2x = -(180.0/M_PI) * (-n[3] * dr2_length*dr2_length + dp_dr2_nrm * dr2[0] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		double d_theta_dr2_deg_d_r2y = -(180.0/M_PI) * (-n[4] * dr2_length*dr2_length + dp_dr2_nrm * dr2[1] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		double d_theta_dr2_deg_d_r2z = -(180.0/M_PI) * (-n[5] * dr2_length*dr2_length + dp_dr2_nrm * dr2[2] ) / ( nrm_length * sqrt( EPS_SMALL +  1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
	
		double d_theta_dr1_deg_d_nrmx = (180.0/M_PI) * (dp_dr1_nrm * n[3] + nrm_length*nrm_length * (-dr1[0])) / (nrm_length*nrm_length*nrm_length * sqrt(  EPS_SMALL + dr1_length * dr1_length - dp_dr1_nrm * dp_dr1_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr1_deg_d_nrmy = (180.0/M_PI) * (dp_dr1_nrm * n[4] + nrm_length*nrm_length * (-dr1[1])) / (nrm_length*nrm_length*nrm_length * sqrt(  EPS_SMALL + dr1_length * dr1_length - dp_dr1_nrm * dp_dr1_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr1_deg_d_nrmz = (180.0/M_PI) * (dp_dr1_nrm * n[5] + nrm_length*nrm_length * (-dr1[2])) / (nrm_length*nrm_length*nrm_length * sqrt(  EPS_SMALL + dr1_length * dr1_length - dp_dr1_nrm * dp_dr1_nrm / nrm_length / nrm_length ) ); 
		
		double d_theta_dr2_deg_d_nrmx = (180.0/M_PI) * (dp_dr2_nrm * n[3] + nrm_length*nrm_length * (-dr2[0])) / (nrm_length*nrm_length*nrm_length * sqrt( EPS_SMALL + dr2_length * dr2_length - dp_dr2_nrm * dp_dr2_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr2_deg_d_nrmy = (180.0/M_PI) * (dp_dr2_nrm * n[4] + nrm_length*nrm_length * (-dr2[1])) / (nrm_length*nrm_length*nrm_length * sqrt( EPS_SMALL +  dr2_length * dr2_length - dp_dr2_nrm * dp_dr2_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr2_deg_d_nrmz = (180.0/M_PI) * (dp_dr2_nrm * n[5] + nrm_length*nrm_length * (-dr2[2])) / (nrm_length*nrm_length*nrm_length * sqrt( EPS_SMALL +  dr2_length * dr2_length - dp_dr2_nrm * dp_dr2_nrm / nrm_length / nrm_length ) ); 
	
		double g_r[9] = { 0,0,0,0,0,0,0,0,0};
		double g_n[9] = { 0,0,0,0,0,0,0,0,0};
	
		g_r[0] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r1x;		
		g_r[1] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r1y;		
		g_r[2] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r1z;		
		
		g_r[3] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r2x;		
		g_r[4] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r2y;		
		g_r[5] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r2z;		
		
		g_r[3] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r2x;		
		g_r[4] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r2y;		
		g_r[5] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r2z;		
		
		g_r[6] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r3x;		
		g_r[7] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r3y;		
		g_r[8] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r3z;		
	
		g_n[3] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_nrmx;
		g_n[4] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_nrmy;
		g_n[5] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_nrmz;
		
		g_n[3] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_nrmx;
		g_n[4] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_nrmy;
		g_n[5] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_nrmz;
	
		// END OF DERIVATIVE OF THETA ENERGIES.
	
		// DERIVATIVE of PHI ENERGY.
	
		// convert to derivative of phi in radians.
		double d_e_d_phi = k_phi * (phi_deg - phi_0) * (180.0/M_PI);
	
		double nrmx = n[3];
		double nrmy = n[4];
		double nrmz = n[5];
	
		double dr1x = dr1[0];
		double dr1y = dr1[1];
		double dr1z = dr1[2];
		
		double dr2x = dr2[0];
		double dr2y = dr2[1];
		double dr2z = dr2[2];
	
		double barPhiK = k_phi * (180.0/M_PI) * (180.0/M_PI);
	
		double dPR1N = dr1[0] * nrmx + dr1[1] * nrmy + dr1[2] * nrmz; 
		double dPR2N = dr2[0] * nrmx + dr2[1] * nrmy + dr2[2] * nrmz; 
	
		double tLen12 =  (dr1[0]-dPR1N*nrmx)*  (dr1[0]-dPR1N*nrmx) +  (dr1[1]-dPR1N*nrmy)*  (dr1[1]-dPR1N*nrmy) + (dr1[2]-dPR1N*nrmz)*  (dr1[2]-dPR1N*nrmz); 
		double tLen22 =  (dr2[0]-dPR2N*nrmx)*  (dr2[0]-dPR2N*nrmx) +  (dr2[1]-dPR2N*nrmy)*  (dr2[1]-dPR2N*nrmy) + (dr2[2]-dPR2N*nrmz)*  (dr2[2]-dPR2N*nrmz); 
	
		double den_dr1x =(barPhiK*(-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr1x*(-1 + Power(nrmx,2)) + 
	            2*nrmx*(-(dr1y*nrmy) - dr1z*nrmz + 
	               dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr2x - dr2x*Power(nrmx,2) + 
	          nrmx*(-(dr2y*nrmy) - dr2z*nrmz + 
	             dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12)*
	     (1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22))))/
	   (2.*Power(tLen12,1.5)*sqrt(tLen22));
	double den_dr1y =(barPhiK*(-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr1y*(-1 + Power(nrmy,2)) + 
	            2*nrmy*(-(dr1x*nrmx) - dr1z*nrmz + 
	               dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr2y - dr2y*Power(nrmy,2) + 
	          nrmy*(-(dr2x*nrmx) - dr2z*nrmz + 
	             dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12)*
	     (1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22))))/
	   (2.*Power(tLen12,1.5)*sqrt(tLen22));
	double den_dr1z =(barPhiK*(-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr1z*(-1 + Power(nrmz,2)) + 
	            2*nrmz*(-(dr1x*nrmx) - dr1y*nrmy + 
	               dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr2z - dr2z*Power(nrmz,2) + 
	          nrmz*(-(dr2x*nrmx) - dr2y*nrmy + 
	             dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12)*
	     (1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22))))/
	   (2.*Power(tLen12,1.5)*sqrt(tLen22));
	
	double den_dr2x = (barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr2x*(-1 + Power(nrmx,2)) + 
	            2*nrmx*(-(dr2y*nrmy) - dr2z*nrmz + 
	               dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr1x - dr1x*Power(nrmx,2) + 
	          nrmx*(-(dr1y*nrmy) - dr1z*nrmz + 
	             dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen22))/
	   (2.*sqrt(tLen12)*Power(tLen22,1.5));
	double den_dr2y =(barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr2y*(-1 + Power(nrmy,2)) + 
	            2*nrmy*(-(dr2x*nrmx) - dr2z*nrmz + 
	               dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr1y - dr1y*Power(nrmy,2) + 
	          nrmy*(-(dr1x*nrmx) - dr1z*nrmz + 
	             dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen22))/
	   (2.*sqrt(tLen12)*Power(tLen22,1.5));
	double den_dr2z =(barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr2z*(-1 + Power(nrmz,2)) + 
	            2*nrmz*(-(dr2x*nrmx) - dr2y*nrmy + 
	               dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr1z - dr1z*Power(nrmz,2) + 
	          nrmz*(-(dr1x*nrmx) - dr1y*nrmy + 
	             dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen22))/
	   (2.*sqrt(tLen12)*Power(tLen22,1.5));
	
	double den_nrmx = (barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (2*dr2x*nrmy*(-dr2y + dPR2N*nrmy) + 2*dr2x*nrmz*(-dr2z + dPR2N*nrmz) - 
	            2*(dr2x - dPR2N*nrmx)*(2*dr2x*nrmx + dr2y*nrmy + dr2z*nrmz))*tLen12) - 
	       ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	        (2*dr1x*nrmy*(-dr1y + dPR1N*nrmy) + 2*dr1x*nrmz*(-dr1z + dPR1N*nrmz) - 
	          2*(dr1x - dPR1N*nrmx)*(2*dr1x*nrmx + dr1y*nrmy + dr1z*nrmz))*tLen22 + 
	       2*(-((2*dr2x - dPR2N*nrmx)*(dr1y*nrmy + dr1z*nrmz)) + 
	          dr1x*(-4*dr2x*nrmx - 2*(dr2y*nrmy + dr2z*nrmz) + 
	             dPR2N*(2*Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))) + 
	          dPR1N*(nrmx*(dr2y*nrmy + dr2z*nrmz) + 
	             dr2x*(2*Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12*tLen22))/
	   (2.*Power(tLen12,1.5)*Power(tLen22,1.5));
	double den_nrmy = (barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (2*dr2y*nrmx*(-dr2x + dPR2N*nrmx) + 2*dr2y*nrmz*(-dr2z + dPR2N*nrmz) - 
	            2*(dr2y - dPR2N*nrmy)*(dr2x*nrmx + 2*dr2y*nrmy + dr2z*nrmz))*tLen12) - 
	       ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	        (2*dr1y*nrmx*(-dr1x + dPR1N*nrmx) + 2*dr1y*nrmz*(-dr1z + dPR1N*nrmz) - 
	          2*(dr1y - dPR1N*nrmy)*(dr1x*nrmx + 2*dr1y*nrmy + dr1z*nrmz))*tLen22 + 
	       2*(-2*dr1x*dr2y*nrmx + dPR1N*dr2y*Power(nrmx,2) + dPR2N*dr1x*nrmx*nrmy + 
	          dPR1N*dr2x*nrmx*nrmy + 2*dPR1N*dr2y*Power(nrmy,2) - 2*dr1z*dr2y*nrmz + 
	          dPR2N*dr1z*nrmy*nrmz + dPR1N*dr2z*nrmy*nrmz + dPR1N*dr2y*Power(nrmz,2) + 
	          dr1y*(-2*dr2x*nrmx - 2*(2*dr2y*nrmy + dr2z*nrmz) + 
	             dPR2N*(Power(nrmx,2) + 2*Power(nrmy,2) + Power(nrmz,2))))*tLen12*tLen22))/
	   (2.*Power(tLen12,1.5)*Power(tLen22,1.5));
	double den_nrmz =(barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (2*dr2z*nrmx*(-dr2x + dPR2N*nrmx) + 2*dr2z*nrmy*(-dr2y + dPR2N*nrmy) - 
	            2*(dr2z - dPR2N*nrmz)*(dr2x*nrmx + dr2y*nrmy + 2*dr2z*nrmz))*tLen12) - 
	       ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	        (2*dr1z*nrmx*(-dr1x + dPR1N*nrmx) + 2*dr1z*nrmy*(-dr1y + dPR1N*nrmy) - 
	          2*(dr1z - dPR1N*nrmz)*(dr1x*nrmx + dr1y*nrmy + 2*dr1z*nrmz))*tLen22 + 
	       2*(-2*dr1x*dr2z*nrmx + dPR1N*dr2z*Power(nrmx,2) - 2*dr1y*dr2z*nrmy + 
	          dPR1N*dr2z*Power(nrmy,2) + dPR2N*dr1x*nrmx*nrmz + dPR1N*dr2x*nrmx*nrmz + 
	          dPR2N*dr1y*nrmy*nrmz + dPR1N*dr2y*nrmy*nrmz + 2*dPR1N*dr2z*Power(nrmz,2) + 
	          dr1z*(-2*dr2x*nrmx - 2*(dr2y*nrmy + 2*dr2z*nrmz) + 
	             dPR2N*(Power(nrmx,2) + Power(nrmy,2) + 2*Power(nrmz,2))))*tLen12*tLen22))/
	   (2.*Power(tLen12,1.5)*Power(tLen22,1.5));
	
	
		double den_d_r1x = den_dr1x; 
		double den_d_r1y = den_dr1y; 
		double den_d_r1z = den_dr1z; 
	
		double den_d_r2x = -(den_dr1x+den_dr2x);
		double den_d_r2y = -(den_dr1y+den_dr2y);
		double den_d_r2z = -(den_dr1z+den_dr2z);
		
		double den_d_r3x = den_dr2x; 
		double den_d_r3y = den_dr2y; 
		double den_d_r3z = den_dr2z; 
	
	
		g_r[0] += den_d_r1x;	
		g_r[1] += den_d_r1y;	
		g_r[2] += den_d_r1z;	
		
		g_r[3] += den_d_r2x;	
		g_r[4] += den_d_r2y;	
		g_r[5] += den_d_r2z;	
		
		g_r[6] += den_d_r3x;	
		g_r[7] += den_d_r3y;	
		g_r[8] += den_d_r3z;	
	
	
		g_n[3] += den_nrmx;
		g_n[4] += den_nrmy;
		g_n[5] += den_nrmz;
	
	
		g_r[0] += bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		g_r[1] += bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		g_r[2] += bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
		
		g_r[3] -= bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		g_r[4] -= bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		g_r[5] -= bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
		
		g_r[3] -= bond_k * ( dr2_length - bond_length) * (dr2[0]/dr2_length);
		g_r[4] -= bond_k * ( dr2_length - bond_length) * (dr2[1]/dr2_length);
		g_r[5] -= bond_k * ( dr2_length - bond_length) * (dr2[2]/dr2_length);
		
		g_r[6] += bond_k * ( dr2_length - bond_length) * (dr2[0]/dr2_length);
		g_r[7] += bond_k * ( dr2_length - bond_length) * (dr2[1]/dr2_length);
		g_r[8] += bond_k * ( dr2_length - bond_length) * (dr2[2]/dr2_length);
		
		double theta_dr1_en = 0.5 *  k_theta * ( (theta_dr1_deg - theta_0)*(theta_dr1_deg - theta_0) ); //add theta energy to bond energy
		double theta_dr2_en = 0.5 *  k_theta * ( (theta_dr2_deg - theta_0)*(theta_dr2_deg - theta_0) );
		double phi_en = 0.5 * k_phi * (180.0/M_PI) * (180.0/M_PI) * pow( (dp_dr1t_dr2t)/(dr1_t_length*dr2_t_length)+1, 2.0 );
	
	
		for( int x = 0; x < 9; x++ ) 
		{
			surfacer_g[x] += g_r[x];
			surfacen_g[x] += g_n[x];
		}
	
		pot = bond_en + theta_dr1_en + theta_dr2_en + phi_en;
		double gmag = g_r[0]*g_r[0]+g_r[1]*g_r[1]+g_r[2]*g_r[2]+g_n[0]*g_n[0]+g_n[1]*g_n[1]+g_n[2]*g_n[2];
		gmag = sqrt(gmag);
//		printf("Grad: %.2le GMAG: %le dr1_len %.2le dr2_len %.2le theta_dr1_deg: %.2le theta_dr2_deg: %.2le dp_dr1t_dr2t: %.2le dr1: %.2le %.2le %.2le dr2: %.2le %.2le %.2le\n", pot, gmag, dr1_length, dr2_length, theta_dr1_deg, theta_dr2_deg, dp_dr1t_dr2t/(dr1_t_length*dr2_t_length), dr1[0], dr1[1], dr1[2], dr2[0], dr2[1], dr2[2] );
	}
	else
	{
		memcpy( r, rall, sizeof(double) * 9 );
 	
		double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	
		double dr2[3] = { r[6] - r[3], r[7] - r[4], r[8] - r[5] };	
	
		theSurface->wrapPBC( dr1, rsurf+theSurface->nv*3 );
		theSurface->wrapPBC( dr2, rsurf+theSurface->nv*3 );
	
			// nx ny nz n[3] n[4] n[5]
		
		//lengths of bonds and normal	
		double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));
		double dr2_length = sqrt((dr2[0]*dr2[0]) + (dr2[1]*dr2[1]) + (dr2[2]*dr2[2]));
		double nrm_length = sqrt((n[3]*n[3]) + (n[4]*n[4]) + (n[5]*n[5]));
	
		double bond_en = 0.5 * bond_k * ( 
				(dr1_length-bond_length)*(dr1_length-bond_length) + 
				(dr2_length-bond_length)*(dr2_length-bond_length) );
		//Calculates dot product between dr1 and dr2
		double dp_dr1_dr2 = dr1[0]*dr2[0] + dr1[1]*dr2[1] + dr1[2]*dr2[2];
		
		//Calculates dot product between dr and normal
		double dp_dr1_nrm = dr1[0]*n[3] + dr1[1]*n[4] + dr1[2]*n[5];
		double dp_dr2_nrm = dr2[0]*n[3] + dr2[1]*n[4] + dr2[2]*n[5];
	
	   //Angles between dr1, dr2 and normal
		double theta_dr1_rad = acos((1-SQRT_SAFETY)*dp_dr1_nrm/(dr1_length*nrm_length));
		double theta_dr1_deg = theta_dr1_rad*180/M_PI;
		double theta_dr2_rad = acos((1-SQRT_SAFETY)*dp_dr2_nrm/(dr2_length*nrm_length));
	    double theta_dr2_deg = theta_dr2_rad*180/M_PI;
	
		//arrays hold xyx positions for tangent vectors
		double dr1_t [3] = {dr1[0] - dp_dr1_nrm * n[3], 
			dr1[1] - dp_dr1_nrm * n[4], 
			dr1[2] - dp_dr1_nrm * n[5]};
	
		double dr2_t[3] = {dr2[0] - dp_dr2_nrm * n[3],
			dr2[1] - dp_dr2_nrm * n[4],
			dr2[2] - dp_dr2_nrm * n[5]};
	
		//Calculates angle between dr1 and dr2
		double theta_rad = acos( (1-SQRT_SAFETY)*(dp_dr1_dr2)/(dr1_length*dr2_length));
		double theta_deg = theta_rad*180.0/M_PI; //angle in degrees
		
		//Dot product between tangent vectors
		double dp_dr1t_dr2t = dr1_t[0]*dr2_t[0] + dr1_t[1]*dr2_t[1] + dr1_t[2]*dr2_t[2];
	
		//lengths of tangent vectors
		double dr1_t_length = sqrt(dr1_t[0]*dr1_t[0] + dr1_t[1]*dr1_t[1] + dr1_t[2]*dr1_t[2]);
		double dr2_t_length = sqrt(dr2_t[0]*dr2_t[0] + dr2_t[1]*dr2_t[1] + dr2_t[2]*dr2_t[2]);
	
		//Calculates angle between dr1_t and dr2_t
		double phi_rad = acos((1-SQRT_SAFETY)*(dp_dr1t_dr2t)/(dr1_t_length*dr2_t_length));
		double phi_deg = phi_rad*180.0/M_PI;
	
		// FIRST part of gradient, theta.
	
		// the derivative of the energy with respect to theta_dr1_deg.
	
		double d_e_d_theta_dr1 = k_theta * ( theta_dr1_deg - theta_0);
		double d_e_d_theta_dr2 = k_theta * ( theta_dr2_deg - theta_0);
	
		double d_theta_dr1_deg_d_r1x = (180.0/M_PI) * (-n[3] * dr1_length*dr1_length + dp_dr1_nrm * dr1[0] ) / ( nrm_length * sqrt( 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r1y = (180.0/M_PI) * (-n[4] * dr1_length*dr1_length + dp_dr1_nrm * dr1[1] ) / ( nrm_length * sqrt( 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r1z = (180.0/M_PI) * (-n[5] * dr1_length*dr1_length + dp_dr1_nrm * dr1[2] ) / ( nrm_length * sqrt( 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		
		double d_theta_dr1_deg_d_r2x = -(180.0/M_PI) * (-n[3] * dr1_length*dr1_length + dp_dr1_nrm * dr1[0] ) / ( nrm_length * sqrt( 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r2y = -(180.0/M_PI) * (-n[4] * dr1_length*dr1_length + dp_dr1_nrm * dr1[1] ) / ( nrm_length * sqrt( 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
		double d_theta_dr1_deg_d_r2z = -(180.0/M_PI) * (-n[5] * dr1_length*dr1_length + dp_dr1_nrm * dr1[2] ) / ( nrm_length * sqrt( 1 - (dp_dr1_nrm * dp_dr1_nrm ) / ( nrm_length*nrm_length*dr1_length*dr1_length) )*dr1_length*dr1_length*dr1_length ); 
	
		double d_theta_dr2_deg_d_r3x = (180.0/M_PI) * (-n[3] * dr2_length*dr2_length + dp_dr2_nrm * dr2[0] ) / ( nrm_length * sqrt( 1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length ); 
		double d_theta_dr2_deg_d_r3y = (180.0/M_PI) * (-n[4] * dr2_length*dr2_length + dp_dr2_nrm * dr2[1] ) / ( nrm_length * sqrt( 1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		double d_theta_dr2_deg_d_r3z = (180.0/M_PI) * (-n[5] * dr2_length*dr2_length + dp_dr2_nrm * dr2[2] ) / ( nrm_length * sqrt( 1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		
		double d_theta_dr2_deg_d_r2x = -(180.0/M_PI) * (-n[3] * dr2_length*dr2_length + dp_dr2_nrm * dr2[0] ) / ( nrm_length * sqrt( 1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		double d_theta_dr2_deg_d_r2y = -(180.0/M_PI) * (-n[4] * dr2_length*dr2_length + dp_dr2_nrm * dr2[1] ) / ( nrm_length * sqrt( 1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
		double d_theta_dr2_deg_d_r2z = -(180.0/M_PI) * (-n[5] * dr2_length*dr2_length + dp_dr2_nrm * dr2[2] ) / ( nrm_length * sqrt( 1 - (dp_dr2_nrm * dp_dr2_nrm ) / ( nrm_length*nrm_length*dr2_length*dr2_length) )*dr2_length*dr2_length*dr2_length  ); 
	
		double d_theta_dr1_deg_d_nrmx = (180.0/M_PI) * (dp_dr1_nrm * n[3] + nrm_length*nrm_length * (-dr1[0])) / (nrm_length*nrm_length*nrm_length * sqrt( dr1_length * dr1_length - dp_dr1_nrm * dp_dr1_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr1_deg_d_nrmy = (180.0/M_PI) * (dp_dr1_nrm * n[4] + nrm_length*nrm_length * (-dr1[1])) / (nrm_length*nrm_length*nrm_length * sqrt( dr1_length * dr1_length - dp_dr1_nrm * dp_dr1_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr1_deg_d_nrmz = (180.0/M_PI) * (dp_dr1_nrm * n[5] + nrm_length*nrm_length * (-dr1[2])) / (nrm_length*nrm_length*nrm_length * sqrt( dr1_length * dr1_length - dp_dr1_nrm * dp_dr1_nrm / nrm_length / nrm_length ) ); 
		
		double d_theta_dr2_deg_d_nrmx = (180.0/M_PI) * (dp_dr2_nrm * n[3] + nrm_length*nrm_length * (-dr2[0])) / (nrm_length*nrm_length*nrm_length * sqrt( dr2_length * dr2_length - dp_dr2_nrm * dp_dr2_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr2_deg_d_nrmy = (180.0/M_PI) * (dp_dr2_nrm * n[4] + nrm_length*nrm_length * (-dr2[1])) / (nrm_length*nrm_length*nrm_length * sqrt( dr2_length * dr2_length - dp_dr2_nrm * dp_dr2_nrm / nrm_length / nrm_length ) ); 
		double d_theta_dr2_deg_d_nrmz = (180.0/M_PI) * (dp_dr2_nrm * n[5] + nrm_length*nrm_length * (-dr2[2])) / (nrm_length*nrm_length*nrm_length * sqrt( dr2_length * dr2_length - dp_dr2_nrm * dp_dr2_nrm / nrm_length / nrm_length ) ); 
	
		double g_r[9] = { 0,0,0,0,0,0,0,0,0};
		double g_n[9] = { 0,0,0,0,0,0,0,0,0};
	
		g_r[0] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r1x;		
		g_r[1] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r1y;		
		g_r[2] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r1z;		
		
		g_r[3] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r2x;		
		g_r[4] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r2y;		
		g_r[5] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_r2z;		
		
		g_r[3] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r2x;		
		g_r[4] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r2y;		
		g_r[5] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r2z;		
		
		g_r[6] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r3x;		
		g_r[7] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r3y;		
		g_r[8] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_r3z;		
	
		g_n[3] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_nrmx;
		g_n[4] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_nrmy;
		g_n[5] += d_e_d_theta_dr1 * d_theta_dr1_deg_d_nrmz;
		
		g_n[3] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_nrmx;
		g_n[4] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_nrmy;
		g_n[5] += d_e_d_theta_dr2 * d_theta_dr2_deg_d_nrmz;
	
		// END OF DERIVATIVE OF THETA ENERGIES.
	
		// DERIVATIVE of PHI ENERGY.
	
		// convert to derivative of phi in radians.
		double d_e_d_phi = k_phi * (phi_deg - phi_0) * (180.0/M_PI);
	
		double nrmx = n[3];
		double nrmy = n[4];
		double nrmz = n[5];
	
		double dr1x = dr1[0];
		double dr1y = dr1[1];
		double dr1z = dr1[2];
		
		double dr2x = dr2[0];
		double dr2y = dr2[1];
		double dr2z = dr2[2];
	
		double barPhiK = k_phi * (180.0/M_PI) * (180.0/M_PI);
	
		double dPR1N = dr1[0] * nrmx + dr1[1] * nrmy + dr1[2] * nrmz; 
		double dPR2N = dr2[0] * nrmx + dr2[1] * nrmy + dr2[2] * nrmz; 
	
		double tLen12 =  (dr1[0]-dPR1N*nrmx)*  (dr1[0]-dPR1N*nrmx) +  (dr1[1]-dPR1N*nrmy)*  (dr1[1]-dPR1N*nrmy) + (dr1[2]-dPR1N*nrmz)*  (dr1[2]-dPR1N*nrmz); 
		double tLen22 =  (dr2[0]-dPR2N*nrmx)*  (dr2[0]-dPR2N*nrmx) +  (dr2[1]-dPR2N*nrmy)*  (dr2[1]-dPR2N*nrmy) + (dr2[2]-dPR2N*nrmz)*  (dr2[2]-dPR2N*nrmz); 
	
		double den_dr1x =(barPhiK*(-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr1x*(-1 + Power(nrmx,2)) + 
	            2*nrmx*(-(dr1y*nrmy) - dr1z*nrmz + 
	               dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr2x - dr2x*Power(nrmx,2) + 
	          nrmx*(-(dr2y*nrmy) - dr2z*nrmz + 
	             dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12)*
	     (1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22))))/
	   (2.*Power(tLen12,1.5)*sqrt(tLen22));
	double den_dr1y =(barPhiK*(-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr1y*(-1 + Power(nrmy,2)) + 
	            2*nrmy*(-(dr1x*nrmx) - dr1z*nrmz + 
	               dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr2y - dr2y*Power(nrmy,2) + 
	          nrmy*(-(dr2x*nrmx) - dr2z*nrmz + 
	             dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12)*
	     (1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22))))/
	   (2.*Power(tLen12,1.5)*sqrt(tLen22));
	double den_dr1z =(barPhiK*(-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr1z*(-1 + Power(nrmz,2)) + 
	            2*nrmz*(-(dr1x*nrmx) - dr1y*nrmy + 
	               dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr2z - dr2z*Power(nrmz,2) + 
	          nrmz*(-(dr2x*nrmx) - dr2y*nrmy + 
	             dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12)*
	     (1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22))))/
	   (2.*Power(tLen12,1.5)*sqrt(tLen22));
	
	double den_dr2x = (barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr2x*(-1 + Power(nrmx,2)) + 
	            2*nrmx*(-(dr2y*nrmy) - dr2z*nrmz + 
	               dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr1x - dr1x*Power(nrmx,2) + 
	          nrmx*(-(dr1y*nrmy) - dr1z*nrmz + 
	             dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen22))/
	   (2.*sqrt(tLen12)*Power(tLen22,1.5));
	double den_dr2y =(barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr2y*(-1 + Power(nrmy,2)) + 
	            2*nrmy*(-(dr2x*nrmx) - dr2z*nrmz + 
	               dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr1y - dr1y*Power(nrmy,2) + 
	          nrmy*(-(dr1x*nrmx) - dr1z*nrmz + 
	             dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen22))/
	   (2.*sqrt(tLen12)*Power(tLen22,1.5));
	double den_dr2z =(barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (-2*dr2z*(-1 + Power(nrmz,2)) + 
	            2*nrmz*(-(dr2x*nrmx) - dr2y*nrmy + 
	               dPR2N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))) + 
	       2*(dr1z - dr1z*Power(nrmz,2) + 
	          nrmz*(-(dr1x*nrmx) - dr1y*nrmy + 
	             dPR1N*(-1 + Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen22))/
	   (2.*sqrt(tLen12)*Power(tLen22,1.5));
	
	double den_nrmx = (barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (2*dr2x*nrmy*(-dr2y + dPR2N*nrmy) + 2*dr2x*nrmz*(-dr2z + dPR2N*nrmz) - 
	            2*(dr2x - dPR2N*nrmx)*(2*dr2x*nrmx + dr2y*nrmy + dr2z*nrmz))*tLen12) - 
	       ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	        (2*dr1x*nrmy*(-dr1y + dPR1N*nrmy) + 2*dr1x*nrmz*(-dr1z + dPR1N*nrmz) - 
	          2*(dr1x - dPR1N*nrmx)*(2*dr1x*nrmx + dr1y*nrmy + dr1z*nrmz))*tLen22 + 
	       2*(-((2*dr2x - dPR2N*nrmx)*(dr1y*nrmy + dr1z*nrmz)) + 
	          dr1x*(-4*dr2x*nrmx - 2*(dr2y*nrmy + dr2z*nrmz) + 
	             dPR2N*(2*Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))) + 
	          dPR1N*(nrmx*(dr2y*nrmy + dr2z*nrmz) + 
	             dr2x*(2*Power(nrmx,2) + Power(nrmy,2) + Power(nrmz,2))))*tLen12*tLen22))/
	   (2.*Power(tLen12,1.5)*Power(tLen22,1.5));
	double den_nrmy = (barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (2*dr2y*nrmx*(-dr2x + dPR2N*nrmx) + 2*dr2y*nrmz*(-dr2z + dPR2N*nrmz) - 
	            2*(dr2y - dPR2N*nrmy)*(dr2x*nrmx + 2*dr2y*nrmy + dr2z*nrmz))*tLen12) - 
	       ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	        (2*dr1y*nrmx*(-dr1x + dPR1N*nrmx) + 2*dr1y*nrmz*(-dr1z + dPR1N*nrmz) - 
	          2*(dr1y - dPR1N*nrmy)*(dr1x*nrmx + 2*dr1y*nrmy + dr1z*nrmz))*tLen22 + 
	       2*(-2*dr1x*dr2y*nrmx + dPR1N*dr2y*Power(nrmx,2) + dPR2N*dr1x*nrmx*nrmy + 
	          dPR1N*dr2x*nrmx*nrmy + 2*dPR1N*dr2y*Power(nrmy,2) - 2*dr1z*dr2y*nrmz + 
	          dPR2N*dr1z*nrmy*nrmz + dPR1N*dr2z*nrmy*nrmz + dPR1N*dr2y*Power(nrmz,2) + 
	          dr1y*(-2*dr2x*nrmx - 2*(2*dr2y*nrmy + dr2z*nrmz) + 
	             dPR2N*(Power(nrmx,2) + 2*Power(nrmy,2) + Power(nrmz,2))))*tLen12*tLen22))/
	   (2.*Power(tLen12,1.5)*Power(tLen22,1.5));
	double den_nrmz =(barPhiK*(1 + ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))/(sqrt(tLen12)*sqrt(tLen22)))*
	     (-(((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	            (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	            (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	          (2*dr2z*nrmx*(-dr2x + dPR2N*nrmx) + 2*dr2z*nrmy*(-dr2y + dPR2N*nrmy) - 
	            2*(dr2z - dPR2N*nrmz)*(dr2x*nrmx + dr2y*nrmy + 2*dr2z*nrmz))*tLen12) - 
	       ((dr1x - dPR1N*nrmx)*(dr2x - dPR2N*nrmx) + 
	          (dr1y - dPR1N*nrmy)*(dr2y - dPR2N*nrmy) + 
	          (dr1z - dPR1N*nrmz)*(dr2z - dPR2N*nrmz))*
	        (2*dr1z*nrmx*(-dr1x + dPR1N*nrmx) + 2*dr1z*nrmy*(-dr1y + dPR1N*nrmy) - 
	          2*(dr1z - dPR1N*nrmz)*(dr1x*nrmx + dr1y*nrmy + 2*dr1z*nrmz))*tLen22 + 
	       2*(-2*dr1x*dr2z*nrmx + dPR1N*dr2z*Power(nrmx,2) - 2*dr1y*dr2z*nrmy + 
	          dPR1N*dr2z*Power(nrmy,2) + dPR2N*dr1x*nrmx*nrmz + dPR1N*dr2x*nrmx*nrmz + 
	          dPR2N*dr1y*nrmy*nrmz + dPR1N*dr2y*nrmy*nrmz + 2*dPR1N*dr2z*Power(nrmz,2) + 
	          dr1z*(-2*dr2x*nrmx - 2*(dr2y*nrmy + 2*dr2z*nrmz) + 
	             dPR2N*(Power(nrmx,2) + Power(nrmy,2) + 2*Power(nrmz,2))))*tLen12*tLen22))/
	   (2.*Power(tLen12,1.5)*Power(tLen22,1.5));
	
	
		double den_d_r1x = den_dr1x; 
		double den_d_r1y = den_dr1y; 
		double den_d_r1z = den_dr1z; 
	
		double den_d_r2x = -(den_dr1x+den_dr2x);
		double den_d_r2y = -(den_dr1y+den_dr2y);
		double den_d_r2z = -(den_dr1z+den_dr2z);
		
		double den_d_r3x = den_dr2x; 
		double den_d_r3y = den_dr2y; 
		double den_d_r3z = den_dr2z; 
	
	
		g_r[0] += den_d_r1x;	
		g_r[1] += den_d_r1y;	
		g_r[2] += den_d_r1z;	
		
		g_r[3] += den_d_r2x;	
		g_r[4] += den_d_r2y;	
		g_r[5] += den_d_r2z;	
		
		g_r[6] += den_d_r3x;	
		g_r[7] += den_d_r3y;	
		g_r[8] += den_d_r3z;	
	
	
		g_n[3] += den_nrmx;
		g_n[4] += den_nrmy;
		g_n[5] += den_nrmz;
	
	
		g_r[0] += bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		g_r[1] += bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		g_r[2] += bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
		
		g_r[3] -= bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		g_r[4] -= bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		g_r[5] -= bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
		
		g_r[3] -= bond_k * ( dr2_length - bond_length) * (dr2[0]/dr2_length);
		g_r[4] -= bond_k * ( dr2_length - bond_length) * (dr2[1]/dr2_length);
		g_r[5] -= bond_k * ( dr2_length - bond_length) * (dr2[2]/dr2_length);
		
		g_r[6] += bond_k * ( dr2_length - bond_length) * (dr2[0]/dr2_length);
		g_r[7] += bond_k * ( dr2_length - bond_length) * (dr2[1]/dr2_length);
		g_r[8] += bond_k * ( dr2_length - bond_length) * (dr2[2]/dr2_length);
		
		double theta_dr1_en = 0.5 *  k_theta * ( (theta_dr1_deg - theta_0)*(theta_dr1_deg - theta_0) ); //add theta energy to bond energy
		double theta_dr2_en = 0.5 *  k_theta * ( (theta_dr2_deg - theta_0)*(theta_dr2_deg - theta_0) );
		double phi_en = 0.5 * k_phi * (180.0/M_PI) * (180.0/M_PI) * pow( (dp_dr1t_dr2t)/(dr1_t_length*dr2_t_length)+1, 2.0 );
	
	
		for( int x = 0; x < 9; x++ ) 
		{
			surfacer_g[x] += g_r[x];
			surfacen_g[x] += g_n[x];
		}
	
		pot = bond_en + theta_dr1_en + theta_dr2_en + phi_en;
		
	}

	return pot;
}

void NBAR::move_inside( void )
{
	pcomplex::move_inside();

	theta_0 = fabs(M_PI-theta_0);
	
}

void NBAR::move_outside( void )
{
	pcomplex::move_outside();

	theta_0 = fabs(theta_0);
	
}
