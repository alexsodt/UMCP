#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mutil.h"
#include "gsl_random_globals.h"
#include "alignSet.h"
#define SQRT_SAFETY (1e-7)

static double monomer_MW = 10000; //amu 
static double monomer_radius = 40;
// initialize the BAR domain in solution.

void dimer::init( double *r )
{
	base_init();

	nsites = 2;
	nattach = 0;

	alloc();

	sigma[0] = monomer_radius;
	sigma[1] = monomer_radius;

	mass[0] = monomer_MW;
	mass[1] = monomer_MW;
	bound = 0;


	// random orientation.

	check_random_init();
		
	// draw random orientation
	double dx,dy,dz;
	gsl_ran_dir_3d( r_gen_global, &dx,&dy,&dz); 
	double dx2,dy2,dz2;
	gsl_ran_dir_3d( r_gen_global, &dx2,&dy2,&dz2); 
	
	double axis[3] ={ dx, dy, dz };
	double rand[3] ={ dx2, dy2, dz2 };

	rall[0] = r[0] + dx * bond_length;
	rall[1] = r[1] + dy * bond_length;
	rall[2] = r[2] + dz * bond_length;
	
	rall[3] = r[0] - dx * bond_length;
	rall[4] = r[1] - dy * bond_length;
	rall[5] = r[2] - dz * bond_length;
}

// initialize the BAR domain on the membrane.

void dimer::init( surface *theSurface, double *rsurf, int f, double u, double v )
{
	base_init();

	nsites = 2;
	nattach = 2;

	alloc();

	bound = 1;
	mass[0] = monomer_MW;
	mass[1] = monomer_MW;
	
	sigma[0] = monomer_radius;
	sigma[1] = monomer_radius;


	double rp[3];
	double nrm[3];
	theSurface->evaluateRNRM( f, u, v, rp, nrm, rsurf);
	
	puv[0] = u;
	puv[1] = v;
	fs[0] = f;
	rall[0] = rp[0];	
	rall[1] = rp[1];	
	rall[2] = rp[2];	

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
	int f_1 = f, nf = f;
	double uv1[2] = { u, v };
	double duv1[2] = { bond_length * druv[0] / len, bond_length * druv[1] / len };
	
	do {
		f_1 = nf;
		nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
	} while( nf != f_1 );

	uv1[0] += duv1[0];		
	uv1[1] += duv1[1];		

	puv[2]   = uv1[0];
	puv[3]   = uv1[1];
	fs[1]    = f_1;	

	theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], rp, nrm, rsurf);

	rall[3] = rp[0];	
	rall[4] = rp[1];	
	rall[5] = rp[2];	

	memcpy( grad_fs, fs, sizeof(int) * 2 );
	memcpy( grad_puv, puv, sizeof(double) * 4 );
}


double dimer::V( Simulation *theSimulation )
{
	double pot = 0;

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
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 6 );
		memset(n, 0, sizeof(double) * 6 ); 
	}
	double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	

	theSimulation->wrapPBC( dr1, theSimulation->alpha );

	//lengths of bonds and normal	
	double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));

	double bond_en = 0.5 * bond_k * ( 
			(dr1_length-bond_length)*(dr1_length-bond_length) ); 

	pot += bond_en;

	return pot;
}

// gets derivative of internal energy relative to position (surfacer_g) and the normal (surfacen_g).

double dimer::grad(Simulation *theSimulation, double *surfacer_g, double *surfacen_g )
{
	double pot = 0;

	double r[6];
	double n[6];

	if( bound )
	{
		for( int s = 0; s < nattach; s++ )
		{
			surface_record *sRec = theSimulation->fetch(sid[s]);
			surface *theSurface = sRec->theSurface;
			double *rsurf = sRec->r;
			int f_1 = fs[s], nf = fs[s];
			double uv1[2] = { 0.33, 0.33 };
			double duv1[2] = { puv[2*s+0]-uv1[0], puv[2*s+1]-uv1[1] };
			
			double null_mom[2] = {0,0};

			coord_transform[4*s+0] = 1;
			coord_transform[4*s+1] = 0;
			coord_transform[4*s+2] = 0;
			coord_transform[4*s+3] = 1;
	
			do {
				f_1 = nf;
				nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf, null_mom, coord_transform+4*s); 
			} while( nf != f_1 );

			uv1[0] += duv1[0];		
			uv1[1] += duv1[1];		
			

			grad_fs[s] = f_1;
			grad_puv[2*s+0] = uv1[0];
			grad_puv[2*s+1] = uv1[1];

			theSurface->evaluateRNRM( f_1, uv1[0], uv1[1], r+3*s, n+3*s, rsurf );  
		}
		
		double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	
	
		theSimulation->wrapPBC( dr1, theSimulation->alpha );
	
			// nx ny nz n[3] n[4] n[5]
		
		//lengths of bonds and normal	
		double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));
	
		double bond_en = 0.5 * bond_k * ( 
				(dr1_length-bond_length)*(dr1_length-bond_length) ); 
	
		surfacer_g[0] += bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		surfacer_g[1] += bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		surfacer_g[2] += bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
		
		surfacer_g[3] -= bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		surfacer_g[4] -= bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		surfacer_g[5] -= bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
	
		pot = bond_en;
	}
	else
	{
		memcpy( r, rall, sizeof(double) * 6 );
 	
		double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	
		
		theSimulation->wrapPBC( dr1, theSimulation->alpha );
	
			// nx ny nz n[3] n[4] n[5]
		
		//lengths of bonds and normal	
		double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));

		double g_r[6];

		g_r[0] += bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		g_r[1] += bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		g_r[2] += bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
		
		g_r[3] -= bond_k * ( dr1_length - bond_length) * (dr1[0]/dr1_length);
		g_r[4] -= bond_k * ( dr1_length - bond_length) * (dr1[1]/dr1_length);
		g_r[5] -= bond_k * ( dr1_length - bond_length) * (dr1[2]/dr1_length);
		
		double bond_en = 0.5 * bond_k * ( 
				(dr1_length-bond_length)*(dr1_length-bond_length) ); 
	
		pot = bond_en;
		
	}

	return pot;
} 

void dimer::loadParams( parameterBlock *block )
{
	bond_length = block->bar_bond_length;
	bond_k = block->default_bond_k;
}

int dimer::getNBonds( void )
{
	return 1;
}

void dimer::putBonds( int *bond_list )
{
	bond_list[0] = 0;
	bond_list[1] = 1;
}

void dimer::bind( int f, double u, double v )
{
	bound = 1;
}
void dimer::unbind( void )
{
	bound = 0;
}


/*

MAB

*/

double MAB::V(  Simulation *theSimulation  )
{
	double v = 0;

	v += dimer::V( theSimulation  );

	// additional potential for their normals.
	
	double r[6];
	double n[6];

	surface_record *sRec = theSimulation->fetch(sid[0]);
	surface *theSurface = sRec->theSurface;
	double *rsurf = sRec->r;

	loadCoords( theSurface, rsurf, r, n );


	double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	

	double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));

	dr1[0] /= dr1_length;
	dr1[1] /= dr1_length;
	dr1[2] /= dr1_length;

	double dp1 = n[0] * dr1[0] + n[1] * dr1[1] + n[2] * dr1[2];
	double dp2 = n[3] * dr1[0] + n[4] * dr1[1] + n[5] * dr1[2];
	
	double dth = (180.0/M_PI)*acos(dp1)-(180.0/M_PI)*acos(dp2);

	v += (k_theta/2) * (dth-dtheta_0) * (dth - dtheta_0);

	return v;
}

double MAB::grad( Simulation *theSimulation, double *surface_g, double *normal_g)
{
	double v = 0;
	surface_record *sRec = theSimulation->fetch(sid[0]);
	surface *theSurface = sRec->theSurface;
	double *rsurf = sRec->r;

	v += dimer::grad( theSimulation, surface_g, normal_g );
	
	double r[6];
	double n[6];

	loadCoords( theSurface, rsurf, r, n );

	double dr1[3] = { r[0] - r[3], r[1] - r[4], r[2] - r[5] };	

	double dr1_length = sqrt((dr1[0]*dr1[0]) + (dr1[1]*dr1[1]) + (dr1[2]*dr1[2]));

	dr1[0] /= dr1_length;
	dr1[1] /= dr1_length;
	dr1[2] /= dr1_length;


	double dp1 = n[0] * dr1[0] + n[1] * dr1[1] + n[2] * dr1[2];
	double dp2 = n[3] * dr1[0] + n[4] * dr1[1] + n[5] * dr1[2];

	double dth = (180.0/M_PI)*acos(dp1)-(180.0/M_PI)*acos(dp2);

	double d_dp1_d_nx = dr1[0];
	double d_dp1_d_ny = dr1[1];
	double d_dp1_d_nz = dr1[2];
	
	double d_dp2_d_nx = dr1[0];
	double d_dp2_d_ny = dr1[1];
	double d_dp2_d_nz = dr1[2];

	double dv_d_dp1 = -k_theta  * (180.0/M_PI)* (dth-dtheta_0 ) / sqrt(1-dp1*dp1);
	double dv_d_dp2 =  k_theta  * (180.0/M_PI)* (dth-dtheta_0 ) / sqrt(1-dp2*dp2);

	v += (k_theta/2) * (dth-dtheta_0) * (dth - dtheta_0);
	// derivatives of the normals:

	normal_g[0] += d_dp1_d_nx * dv_d_dp1;
	normal_g[1] += d_dp1_d_ny * dv_d_dp1;
	normal_g[2] += d_dp1_d_nz * dv_d_dp1;
	
	normal_g[3] += d_dp2_d_nx * dv_d_dp2;
	normal_g[4] += d_dp2_d_ny * dv_d_dp2;
	normal_g[5] += d_dp2_d_nz * dv_d_dp2;

	// derivatives of the positions:

	double d_dp1_dr1x = n[0] * (dr1[1]*dr1[1] + dr1[2]*dr1[2]) / dr1_length;
	double d_dp1_dr1y = n[1] * (dr1[0]*dr1[0] + dr1[2]*dr1[2]) / dr1_length;
	double d_dp1_dr1z = n[2] * (dr1[0]*dr1[0] + dr1[1]*dr1[1]) / dr1_length;

	d_dp1_dr1x -= n[1] * dr1[0]*dr1[1] / dr1_length; 
	d_dp1_dr1x -= n[2] * dr1[0]*dr1[2] / dr1_length; 

	d_dp1_dr1y -= n[0] * dr1[1]*dr1[0] / dr1_length; 
	d_dp1_dr1y -= n[2] * dr1[1]*dr1[2] / dr1_length; 

	d_dp1_dr1z -= n[0] * dr1[2]*dr1[0] / dr1_length; 
	d_dp1_dr1z -= n[1] * dr1[2]*dr1[1] / dr1_length; 

	double d_dp2_dr1x = -n[3] * (dr1[1]*dr1[1] + dr1[2]*dr1[2]) / dr1_length;
	double d_dp2_dr1y = -n[4] * (dr1[0]*dr1[0] + dr1[2]*dr1[2]) / dr1_length;
	double d_dp2_dr1z = -n[5] * (dr1[0]*dr1[0] + dr1[1]*dr1[1]) / dr1_length;

	d_dp2_dr1x += n[4] * dr1[0]*dr1[1] / dr1_length; 
	d_dp2_dr1x += n[5] * dr1[0]*dr1[2] / dr1_length; 

	d_dp2_dr1y += n[3] * dr1[1]*dr1[0] / dr1_length; 
	d_dp2_dr1y += n[5] * dr1[1]*dr1[2] / dr1_length; 

	d_dp2_dr1z += n[3] * dr1[2]*dr1[0] / dr1_length; 
	d_dp2_dr1z += n[4] * dr1[2]*dr1[1] / dr1_length; 

	surface_g[0] += d_dp1_dr1x * dv_d_dp1 - d_dp2_dr1x * dv_d_dp2; 
	surface_g[1] += d_dp1_dr1y * dv_d_dp1 - d_dp2_dr1y * dv_d_dp2; 
	surface_g[2] += d_dp1_dr1z * dv_d_dp1 - d_dp2_dr1z * dv_d_dp2; 
	
	surface_g[3] += d_dp2_dr1x * dv_d_dp2 - d_dp1_dr1x * dv_d_dp1; 
	surface_g[4] += d_dp2_dr1y * dv_d_dp2 - d_dp1_dr1y * dv_d_dp1; 
	surface_g[5] += d_dp2_dr1z * dv_d_dp2 - d_dp1_dr1z * dv_d_dp1; 

	return v;
}

void MAB::loadParams( parameterBlock *block )
{
	dimer::loadParams( block );

	k_theta = block->mab_k_theta;
	dtheta_0 = block->mab_d_theta;
	bond_k = block->mab_bond_k;
	bond_length = block->mab_bond_length;
}

void MAB::move_inside( void )
{
	pcomplex::move_inside();

	dtheta_0 = -fabs(dtheta_0);
}

void MAB::move_outside( void )
{
	pcomplex::move_outside();

	dtheta_0 = +fabs(dtheta_0);
}


