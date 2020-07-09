#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mutil.h"
#include "gsl_random_globals.h"
#include "alignSet.h"
#include "pdbFetch.h"

#define SQRT_SAFETY (1e-7)

#define EPS_SMALL (1e-14)
#define WORKING

static double scale_k = 0.01;
	
static double attach_bond_k = 10.0 * scale_k;
static double inter_bond_k = 10.0 * scale_k;
static double dihedral_k = 1000.0 * scale_k;
static double angle_k = 3000.0 * scale_k;
static int nbonds=5;
static int nangles=6;
static int ndihe=1;


static int bonds[5][2] = {
	{0,4},
	{1,4},
	{2,5},
	{3,5},
	{4,5}
};	

static int angles[6][3] = {
	{0,4,1},
	{0,4,5},
	{1,4,5},
	{2,5,3},
	{2,5,4},
	{3,5,4}
};

static int dihedrals[1][4] = {
	3,5,4,1
};


static double monomer_MW = 10*26492.14; //amu 
static double C2_P_RADIUS = 25.0;
static double attach_p_radius = 5.0;
static double bond_length_attach_long = 45.0;
static double bond_length_attach_short = 25.0;
static double bond_length_inter = 30.0;

static double bond_k[5][2] =
{
	{ bond_length_attach_long, attach_bond_k },	
	{ bond_length_attach_short, attach_bond_k },	
	{ bond_length_attach_long, attach_bond_k },	
	{ bond_length_attach_short, attach_bond_k },	
	{ bond_length_inter, inter_bond_k }
};

void syt7::init( double *r )
{
	printf("Currently Syt7 can only be initialized on the membrane.\n");
	exit(1);
}

// initialize the BAR domain on the membrane.

void syt7::init( Simulation *theSimulation, surface *theSurface, double *rsurf, int f, double u, double v )
{
	// assume for now this is one of the points on the membrane neck.

	base_init();

	nsites = 6;
	nattach = 4;

	alloc();

	double vdw_r = 5.0;
	sigma[0] = vdw_r;
	sigma[1] = vdw_r;
	sigma[2] = vdw_r;
	sigma[3] = vdw_r;
	sigma[4] = vdw_r;
	sigma[5] = vdw_r;

	mass[0] = monomer_MW*100;
	mass[1] = monomer_MW*100;
	mass[2] = monomer_MW*100;
	mass[3] = monomer_MW*100;
	mass[4] = monomer_MW*100;
	mass[5] = monomer_MW*100;

	bound = 1;
	
	double rpt_attach1[3], nrm_attach1[3];
	theSurface->evaluateRNRM( f, u, v, rpt_attach1, nrm_attach1, rsurf );

	double rp[3];
	double nrm[3];
	theSurface->evaluateRNRM( f, u, v, rp, nrm, rsurf);

	sid[0] = theSurface->surface_id;
	sid[1] = theSurface->surface_id;
	sid[2] = theSurface->surface_id;
	sid[3] = theSurface->surface_id;
	
	//
	// This is the neck-attachment point of the first protein.
	// goes into site-position "0"

	puv[0] = u;
	puv[1] = v;
	fs[0] = f;

	// find the other points on the membrane.
	// assume we are on a saddle, in which case we will move along the positive curvature direction to place the second point.


	double dr_u[3];
	theSurface->ru( f, u, v, rsurf, dr_u );
	double drdu = normalize(dr_u);
	double dr_v[3];
	theSurface->rv( f, u, v, rsurf, dr_v );
	double drdv = normalize(dr_v);

	double c_vec_1[3];
	double c_vec_2[3];
	double c_val1, c_val2;
	double k;
	theSurface->c( f, u, v, rsurf, &k, c_vec_1, c_vec_2, &c_val1, &c_val2 ); 

	//
	// Find the neck-attachment point of the second protein.
	// goes into site-position "2"

	double d_move_1[2] = { c_vec_1[0], c_vec_1[1] };
	double d_move_2[2] = { c_vec_2[0], c_vec_2[1] };

	if( c_val2 > c_val1 )
	{
		d_move_1[0] = c_vec_2[0];
		d_move_1[1] = c_vec_2[1];
		d_move_2[0] = c_vec_1[0];
		d_move_2[1] = c_vec_1[1];
	} 
	
	double d_move_expec[3] = {0,0,0};
	
	d_move_expec[0] += d_move_1[0] * dr_u[0]*drdu;
	d_move_expec[1] += d_move_1[0] * dr_u[1]*drdu;
	d_move_expec[2] += d_move_1[0] * dr_u[2]*drdu;
	
	d_move_expec[0] += d_move_1[1] * dr_v[0]*drdv;
	d_move_expec[1] += d_move_1[1] * dr_v[1]*drdv;
	d_move_expec[2] += d_move_1[1] * dr_v[2]*drdv;
	
	// save this to see if we need to reorient the other normal.
	double save_lat[3];
	memcpy( save_lat, d_move_expec, sizeof(double) * 3 );

	double len = normalize(d_move_expec);
	double r_move = bond_length_inter;

	d_move_1[0] *= r_move/len;
	d_move_1[1] *= r_move/len;

	int f_1 = f, nf = f;
	double uv1[2] = { u, v };
	double duv1[2] = { d_move_1[0], d_move_1[1] }; 
	
	do {
		f_1 = nf;
		nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
	} while( nf != f_1 );
	
	// the neck-attachment point of the second protein.
	puv[4] = uv1[0];
	puv[5] = uv1[1];
	fs[2] = f_1;

	//
	// find the upper-leaflet attachment point of the first protein
	//	

	uv1[0] = u;
	uv1[1] = v;

	duv1[0] = d_move_2[0];
	duv1[1] = d_move_2[1];

	
	
	d_move_expec[0] = duv1[0] * dr_u[0]*drdu;
	d_move_expec[1] = duv1[0] * dr_u[1]*drdu;
	d_move_expec[2] = duv1[0] * dr_u[2]*drdu;
	
	d_move_expec[0] += duv1[1] * dr_v[0]*drdv;
	d_move_expec[1] += duv1[1] * dr_v[1]*drdv;
	d_move_expec[2] += duv1[1] * dr_v[2]*drdv;
	
	double test_v[3];

	cross(d_move_expec, save_lat, test_v );

	double dp = test_v[0] * nrm_attach1[0] + test_v[1] * nrm_attach1[1] + test_v[2] * nrm_attach1[2];

	if( dp > 0 )
	{
		duv1[0] *= -1;
		duv1[1] *= -1;
		d_move_expec[0] *= -1;
		d_move_expec[1] *= -1;
		d_move_expec[2] *= -1;
	}

	len = normalize(d_move_expec);
	r_move = bond_length_attach_short;
	// save the attachment displacement: we want the second-protein attachment displacement to be in the opposite direction (approximately).
	double attach_displacement_1[3] = { d_move_expec[0], d_move_expec[1], d_move_expec[2] };

	duv1[0] *= (r_move/len);
	duv1[1] *= (r_move/len);
	
	f_1 = f;
	nf = f;
	
	do {
		f_1 = nf;
		nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
	} while( nf != f_1 );
	
	// the upper-leaflet attachment point of the first protein.
	puv[2] = uv1[0];
	puv[3] = uv1[1];
	fs[1] = f_1;

	//
	// find the lower-leaflet attachment point of the second protein
	//	

	uv1[0] = puv[4];
	uv1[1] = puv[5];
	f_1 = fs[2];
	nf = fs[2];
	
	double dr_u_2[3];
	theSurface->ru( f_1, uv1[0], uv1[1], rsurf, dr_u_2 );
	double drdu2 = normalize(dr_u_2);
	double dr_v_2[3];
	theSurface->rv( f_1, uv1[0],uv1[1], rsurf, dr_v_2 );
	double drdv2 = normalize(dr_v_2);

	double c_vec_1_2[3];
	double c_vec_2_2[3];
	double c_val1_2, c_val2_2;
	double k_2;
	theSurface->c( f_1, uv1[0], uv1[1], rsurf, &k_2, c_vec_1_2, c_vec_2_2, &c_val1_2, &c_val2_2 ); 

	// choose the negative curvature direction.

	duv1[0] = c_vec_1_2[0];	
	duv1[1] = c_vec_1_2[1];	

	if( c_val1_2 > c_val2_2 )
	{
		duv1[0] = c_vec_2_2[0];	
		duv1[1] = c_vec_2_2[1];	
	}
	
	d_move_expec[0] = duv1[0] * dr_u_2[0]*drdu2;
	d_move_expec[1] = duv1[0] * dr_u_2[1]*drdu2;
	d_move_expec[2] = duv1[0] * dr_u_2[2]*drdu2;
	
	d_move_expec[0] += duv1[1] * dr_v_2[0]*drdv2;
	d_move_expec[1] += duv1[1] * dr_v_2[1]*drdv2;
	d_move_expec[2] += duv1[1] * dr_v_2[2]*drdv2;

	len = normalize(d_move_expec);
	double dpx =
		 d_move_expec[0] * attach_displacement_1[0] + 
		 d_move_expec[1] * attach_displacement_1[1] +
		 d_move_expec[2] * attach_displacement_1[2];

	if( dpx > 0 )
	{
		// wrong direction, send the other way.

		duv1[0] *= -1;
		duv1[1] *= -1;
	} 

	r_move = bond_length_attach_short;
	
	duv1[0] *= (r_move/len);
	duv1[1] *= (r_move/len);
	
	do {
		f_1 = nf;
		nf = theSurface->nextFace( f_1, uv1+0, uv1+1, duv1+0, duv1+1, rsurf ); 
	} while( nf != f_1 );

	puv[6] = uv1[0];
	puv[7] = uv1[1];
	fs[3] = f_1;

	// now for the ``aqueous'' sites.  they will be projected off the normal of the attachment sites.

	double io_sign = (is_inside ? 1 : -1 );

	rall[4*3+0] = rpt_attach1[0] + nrm_attach1[0] * bond_length_attach_long * io_sign;	
	rall[4*3+1] = rpt_attach1[1] + nrm_attach1[1] * bond_length_attach_long * io_sign;	
	rall[4*3+2] = rpt_attach1[2] + nrm_attach1[2] * bond_length_attach_long * io_sign;	
	
	double rpt_attach2[3], nrm_attach2[3];
	theSurface->evaluateRNRM( fs[2], puv[4], puv[5], rpt_attach2, nrm_attach2, rsurf );

	rall[5*3+0] = rpt_attach2[0] + nrm_attach2[0] * bond_length_attach_long * io_sign;	
	rall[5*3+1] = rpt_attach2[1] + nrm_attach2[1] * bond_length_attach_long * io_sign;	
	rall[5*3+2] = rpt_attach2[2] + nrm_attach2[2] * bond_length_attach_long * io_sign;	

	memcpy( grad_fs, fs, sizeof(int) * 4 );
	memcpy( grad_puv, puv, sizeof(double) * 4*2 );

	memset( PBC_ext, 0, sizeof(double) * 3*6 );
	
	setrall(theSimulation);
	
	// wrap them all to one atom

	for( int a = 1; a < 6; a++ )
	{
		double dr_r[3] = { rall[3*a+0] - rall[0], rall[3*a+1] - rall[1], rall[3*a+2] - rall[2] };

		double put[3];

		MinImage3D( dr_r, theSimulation->PBC_vec, put, theSimulation->alpha );

		PBC_ext[3*a+0] = put[0];
		PBC_ext[3*a+1] = put[1];
		PBC_ext[3*a+2] = put[2];
	}

}

// custom orient procedure.



void syt7::bind( int f, double u, double v)
{
	bound = 1;
}

void syt7::unbind( void )
{
	bound = 0;
}


void syt7::loadParams( parameterBlock *block )
{
}


int syt7::getNBonds( void )
{
	return 5;
}

void syt7::putBonds( int *bond_list )
{
	bond_list[0] = 0;
	bond_list[1] = 4;
	
	bond_list[2] = 1;
	bond_list[3] = 4;
	
	bond_list[4] = 2;
	bond_list[5] = 5;
	
	bond_list[6] = 3;
	bond_list[7] = 5;

	bond_list[8] = 4;
	bond_list[9] = 5;
}

double syt7::V( Simulation *theSimulation  )
{
	double *alphas = theSimulation->alpha;

	double r[3*nsites];
	double n[3*nsites];

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
		
		memcpy(r+4*3 , rall+4*3, sizeof(double) * 3*2 );
	}
	else
	{
		memcpy(r , rall, sizeof(double) * 3*nsites );
		memset(n, 0, sizeof(double) * 3*nsites ); 
	}

	double pot = 0;

	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = bonds[b][0];
		int a2 = bonds[b][1];

		double r0 = bond_k[b][0];
		double k = bond_k[b][1];

		double dr[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double ln = normalize(dr);

		pot += 0.5 * k * (ln - r0) * (ln-r0);	
	}
	
	for( int a = 0; a < nangles; a++ )
	{
		int a1 = angles[a][0];
		int a2 = angles[a][1];
		int a3 = angles[a][2];

		double dr1[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr2[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		double ln1 = normalize(dr1);
		double ln2 = normalize(dr2);

		double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];
		double nrm_dp = dp;
		pot += 0.5 * angle_k *nrm_dp*nrm_dp;	
	}

	for( int d = 0; d < ndihe; d++ )
	{
		int a1 = dihedrals[d][0];
		int a2 = dihedrals[d][1];
		int a3 = dihedrals[d][2];
		int a4 = dihedrals[d][3];
		
		double dr1A[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr1B[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		
		double dr2A[3] = { r[3*a2+0] - r[3*a3+0],
				 r[3*a2+1] - r[3*a3+1],
				 r[3*a2+2] - r[3*a3+2] };
		double dr2B[3] = { r[3*a4+0] - r[3*a3+0],
				 r[3*a4+1] - r[3*a3+1],
				 r[3*a4+2] - r[3*a3+2] };

		double n1[3], n2[3];

		cross( dr1A,dr1B,n1);
		cross( dr2A,dr2B,n2);
		normalize(n1);
		normalize(n2);
		double dp = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
		double theta = acos(dp);

		pot += 0.5 * dihedral_k * (theta - M_PI) * (theta - M_PI);
	}	

	return pot;
}

// gets derivative of internal energy relative to position (surfacer_g) and the normal (surfacen_g).

double syt7::grad( Simulation *theSimulation,  double *surfacer_g, double *surfacen_g )
{
	double *alphas = theSimulation->alpha;
	double r[3*nsites];
	double n[3*nsites];

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
		
		memcpy(r+4*3 , rall+4*3, sizeof(double) * 3*2 );
	}

	double pot = 0;

	for( int b = 0; b < nbonds; b++ )
	{
		int a1 = bonds[b][0];
		int a2 = bonds[b][1];

		double r0 = bond_k[b][0];
		double k = bond_k[b][1];

		double dr[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double ln = normalize(dr);

		pot += 0.5 * k * (ln - r0) * (ln-r0);	

		double f = k * (ln-r0); // dV / dx: 
		
		surfacer_g[3*a1+0] += dr[0] * f;
		surfacer_g[3*a1+1] += dr[1] * f;
		surfacer_g[3*a1+2] += dr[2] * f;

		surfacer_g[3*a2+0] -= dr[0] * f;
		surfacer_g[3*a2+1] -= dr[1] * f;
		surfacer_g[3*a2+2] -= dr[2] * f;
	}
	
	for( int a = 0; a < nangles; a++ )
	{
		int a1 = angles[a][0];
		int a2 = angles[a][1];
		int a3 = angles[a][2];

		double dr1[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr2[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		double ln1 = normalize(dr1);
		double ln2 = normalize(dr2);

		double dp = dr1[0]*dr2[0]+dr1[1]*dr2[1]+dr1[2]*dr2[2];

		// numerator derivative

		surfacer_g[3*a1+0] += angle_k * dp * dr2[0] / ln1;
		surfacer_g[3*a1+1] += angle_k * dp * dr2[1] / ln1;
		surfacer_g[3*a1+2] += angle_k * dp * dr2[2] / ln1;
		
		surfacer_g[3*a2+0] += angle_k * dp * (-dr2[0]) / ln1;
		surfacer_g[3*a2+1] += angle_k * dp * (-dr2[1]) / ln1;
		surfacer_g[3*a2+2] += angle_k * dp * (-dr2[2]) / ln1;
		
		surfacer_g[3*a2+0] += angle_k * dp * (-dr1[0]) / ln2;
		surfacer_g[3*a2+1] += angle_k * dp * (-dr1[1]) / ln2;
		surfacer_g[3*a2+2] += angle_k * dp * (-dr1[2]) / ln2;
		
		surfacer_g[3*a3+0] += angle_k * dp * (dr1[0]) / ln2;
		surfacer_g[3*a3+1] += angle_k * dp * (dr1[1]) / ln2;
		surfacer_g[3*a3+2] += angle_k * dp * (dr1[2]) / ln2;

		// denominator derivative.
		
		double ln1_2 = ln1*ln1;
		double ln2_2 = ln2*ln2;

		surfacer_g[3*a1+0] += -angle_k * dp * dp * dr1[0] / ln1;
		surfacer_g[3*a1+1] += -angle_k * dp * dp * dr1[1] / ln1;
		surfacer_g[3*a1+2] += -angle_k * dp * dp * dr1[2] / ln1;
		
		surfacer_g[3*a2+0] += angle_k * dp * dp * dr1[0] / ln1;
		surfacer_g[3*a2+1] += angle_k * dp * dp * dr1[1] / ln1;
		surfacer_g[3*a2+2] += angle_k * dp * dp * dr1[2] / ln1;
		
		surfacer_g[3*a2+0] += angle_k * dp * dp * dr2[0] / ln2;
		surfacer_g[3*a2+1] += angle_k * dp * dp * dr2[1] / ln2;
		surfacer_g[3*a2+2] += angle_k * dp * dp * dr2[2] / ln2;
		
		surfacer_g[3*a3+0] += -angle_k * dp * dp * dr2[0] / ln2;
		surfacer_g[3*a3+1] += -angle_k * dp * dp * dr2[1] / ln2;
		surfacer_g[3*a3+2] += -angle_k * dp * dp * dr2[2] / ln2;

		pot += 0.5 * angle_k * dp*dp;	
	}

	for( int d = 0; d < ndihe; d++ )
	{
		int a1 = dihedrals[d][0];
		int a2 = dihedrals[d][1];
		int a3 = dihedrals[d][2];
		int a4 = dihedrals[d][3];
		
		double dr1A[3] = { r[3*a1+0] - r[3*a2+0],
				 r[3*a1+1] - r[3*a2+1],
				 r[3*a1+2] - r[3*a2+2] };
		double dr1B[3] = { r[3*a3+0] - r[3*a2+0],
				 r[3*a3+1] - r[3*a2+1],
				 r[3*a3+2] - r[3*a2+2] };
		
		double dr2A[3] = { r[3*a2+0] - r[3*a3+0],
				 r[3*a2+1] - r[3*a3+1],
				 r[3*a2+2] - r[3*a3+2] };
		double dr2B[3] = { r[3*a4+0] - r[3*a3+0],
				 r[3*a4+1] - r[3*a3+1],
				 r[3*a4+2] - r[3*a3+2] };

		double n1[3], n2[3];

		cross( dr1A,dr1B,n1);
		cross( dr2A,dr2B,n2);
		normalize(n1);
		normalize(n2);
		double dp = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
		double theta = acos(dp);

		// derivative of nrm wrt the three vectors.
		// slow: cartesian component of the normal
		// mid: which vector (r1,r2,r3)
		// fast: cartesian component of that vector
		double dnrm1[27];
		double dnrm2[27];

		// order here is mixed up, sorry.
		normal_cp_der( r+3*a2, r+3*a1, r+3*a3, dnrm1 );
		normal_cp_der( r+3*a3, r+3*a2, r+3*a4, dnrm2 );

		// everything prior to d(dp)/dr1
		double base_der = dihedral_k * (theta - M_PI) * (-1 / sqrt(1-dp*dp));
	
		// everything else is d(dp)/dx_i, etc.
		for( int pc = 0; pc < 3; pc++ )
		for( int nc = 0; nc < 3; nc++ )
		{
			surfacer_g[3*a1+pc] += base_der * n2[nc] * dnrm1[nc*9+1*3+pc];
			surfacer_g[3*a2+pc] += base_der * n2[nc] * dnrm1[nc*9+0*3+pc];
			surfacer_g[3*a3+pc] += base_der * n2[nc] * dnrm1[nc*9+2*3+pc];
			
			surfacer_g[3*a2+pc] += base_der * n1[nc] * dnrm2[nc*9+1*3+pc];
			surfacer_g[3*a3+pc] += base_der * n1[nc] * dnrm2[nc*9+0*3+pc];
			surfacer_g[3*a4+pc] += base_der * n1[nc] * dnrm2[nc*9+2*3+pc];
		}
	
		pot += 0.5 * dihedral_k * (theta - M_PI) * (theta - M_PI);
	}	

	return pot;
}

void syt7::move_inside( void )
{
	pcomplex::move_inside();
}

void syt7::move_outside( void )
{
	pcomplex::move_outside();
}

void syt7::writeStructure( Simulation *theSimulation, struct atom_rec **at_out, int *nat_out )
{
	struct atom_rec *C2A = NULL;
	struct atom_rec *C2B = NULL;

	int nC2A=0;
	int nC2B=0;
	int nNTH=0;	

	pdbFetch( &C2A, &nC2A, "syt7", "C2A" );
	pdbFetch( &C2B, &nC2B, "syt7", "C2B" );
#ifdef DO_NTH	
	struct atom_rec *NTH = NULL;
	pdbFetch( &NTH, &nNTH, "syt7", "NTH" );
#endif
	double io_sign = (is_inside ? 1 : -1 );

	//
	// alignments: rotate and translate the protein to align it to the sub-sites.
	//

	//	
	//	FIRST: C2A. match virtual sites of model to specific atoms in the pdb, then do alignment.
	//
	double site_offset = 20.0;
	double virtual_site[3] = { rall[4*3+0] + 0.3 * ( rall[3]-rall[4*3+0]),
				   rall[4*3+1] + 0.3 * ( rall[4]-rall[4*3+1]),
				   rall[4*3+2] + 0.3 * ( rall[5]-rall[4*3+2]) };
	double vec_1[3];
	double virtual_site_neck[3];
	
	vec_1[0]= rall[0] - rall[4*3+0];
	vec_1[1]= rall[1] - rall[4*3+1];
	vec_1[2]= rall[2] - rall[4*3+2];

	normalize(vec_1);

	virtual_site_neck[0] = rall[0*3+0] - vec_1[0] * site_offset;
	virtual_site_neck[1] = rall[0*3+1] - vec_1[1] * site_offset; 
	virtual_site_neck[2] = rall[0*3+2] - vec_1[2] * site_offset; 

	// the three "C2A" sites:
	double C2A_align[9] = { 
				virtual_site_neck[0], virtual_site_neck[1], virtual_site_neck[2],
				virtual_site[0], virtual_site[1], virtual_site[2],
				rall[4*3+0], rall[4*3+1], rall[4*3+2] };

	int set_C2A[3] = {0,1,2};
	int set_C2A_pdb[3] = {-1,-1,-1};
	int res_C2A_pdb[3] = { 166, 185, 154 };
	const char *at_C2A_pdb[3] = { "CA", "CA", "CA" };
	double *pcopy = (double *)malloc( sizeof(double) * 3 * nC2A );

	for( int t = 0; t < nC2A; t++ )
	{
		// copy the coordinates
		pcopy[3*t+0] = C2A[t].x;
		pcopy[3*t+1] = C2A[t].y;
		pcopy[3*t+2] = C2A[t].z;
		// look for the sites
		for( int x = 0; x < 3; x++ )
		{
			if( res_C2A_pdb[x] == C2A[t].res && !strcasecmp( C2A[t].atname, at_C2A_pdb[x] ) )
				set_C2A_pdb[x] = t;
		}
	}	
	
	if( set_C2A_pdb[0] < 0 || set_C2A_pdb[1] < 0 || set_C2A_pdb[2] < 0 )
	{
		printf("Failed to find site connections for SYT7.\n");
		exit(1);
	}
	
	alignStructuresOnAtomSet( C2A_align, set_C2A, pcopy, set_C2A_pdb, 3, nC2A ); 
	
	for( int t = 0; t < nC2A; t++ )
	{
		// copy the coordinates out
		C2A[t].x=pcopy[3*t+0];
		C2A[t].y=pcopy[3*t+1];
		C2A[t].z=pcopy[3*t+2];
	}

	//
	//	NEXT: C2B
	//

	virtual_site[0] = rall[5*3+0] + 0.3 * (rall[9]  - rall[5*3+0]);
	virtual_site[1] = rall[5*3+1] + 0.3 * (rall[10] - rall[5*3+1]);
	virtual_site[2] = rall[5*3+2] + 0.3 * (rall[11] - rall[5*3+2]);
	
	vec_1[0]= rall[2*3+0] - rall[5*3+0];
	vec_1[1]= rall[2*3+1] - rall[5*3+1];
	vec_1[2]= rall[2*3+2] - rall[5*3+2];

	normalize(vec_1);

	virtual_site_neck[0] = rall[2*3+0] - vec_1[0] * site_offset;
	virtual_site_neck[1] = rall[2*3+1] - vec_1[1] * site_offset; 
	virtual_site_neck[2] = rall[2*3+2] - vec_1[2] * site_offset; 

	
	// the three "C2A" sites:
	double C2B_align[9] = { virtual_site_neck[0], virtual_site_neck[1], virtual_site_neck[2],
				virtual_site[0], virtual_site[1], virtual_site[2],
				rall[5*3+0], rall[5*3+1], rall[5*3+2] };

	int set_C2B[3] = {0,1,2};
	int set_C2B_pdb[3] = {-1,-1,-1};
	int res_C2B_pdb[3] = { 297, 319, 389 };
	const char *at_C2B_pdb[3] = { "CA", "CA", "CA" };

	pcopy = (double *)realloc( pcopy, sizeof(double) * 3 * nC2B );

	for( int t = 0; t < nC2B; t++ )
	{
		// copy the coordinates
		pcopy[3*t+0] = C2B[t].x;
		pcopy[3*t+1] = C2B[t].y;
		pcopy[3*t+2] = C2B[t].z;
		// look for the sites
		for( int x = 0; x < 3; x++ )
		{
			if( res_C2B_pdb[x] == C2B[t].res && !strcasecmp( C2B[t].atname, at_C2B_pdb[x] ) )
				set_C2B_pdb[x] = t;
		}
	}	
	
	if( set_C2B_pdb[0] < 0 || set_C2B_pdb[1] < 0 || set_C2B_pdb[2] < 0 )
	{
		printf("Failed to find site connections for SYT7.\n");
		exit(1);
	}
	
	alignStructuresOnAtomSet( C2B_align, set_C2B, pcopy, set_C2B_pdb, 3, nC2B ); 
	
	for( int t = 0; t < nC2B; t++ )
	{
		// copy the coordinates out
		C2B[t].x=pcopy[3*t+0];
		C2B[t].y=pcopy[3*t+1];
		C2B[t].z=pcopy[3*t+2];
	}

#ifdef DO_NTH
	//
	// Finally: the N-terminal transmembrane helix.
	//
	//

	// get the normal.

	surface *theSurface;
	double *rsurf;
	theSimulation->fetch(sid[1],&theSurface,&rsurf);

	double rpt[3], rnrm[3];
	theSurface->evaluateRNRM( fs[1], puv[1*2+0], puv[1*2+1], rpt, rnrm, rsurf );	
	
	rnrm[0] *= io_sign;
	rnrm[1] *= io_sign;
	rnrm[2] *= io_sign;

	virtual_site[0] = rall[3] + rnrm[0] * 12.0; 
	virtual_site[1] = rall[4] + rnrm[1] * 12.0;
	virtual_site[2] = rall[5] + rnrm[2] * 12.0;
	
	// another arbitrary site we could configure later if it matters?
	
	double arb[3] = { rand(), rand(), rand() };
	normalize(arb);
	double perp[3];
	cross( rnrm, arb, perp );
	normalize(perp);
	
	double virtual_site_2[3] = { virtual_site[0] + perp[0] *2,
				     virtual_site[1] + perp[1] *2,
				     virtual_site[2] + perp[2] *2 };

	double NTH_align[9] = { rall[3], rall[4], rall[5],
				virtual_site[0], virtual_site[1], virtual_site[2],
				virtual_site_2[0], virtual_site_2[1], virtual_site_2[2],
				 };

	int set_NTH[3] = {0,1,2};
	int set_NTH_pdb[3] = {-1,-1,-1};
	int res_NTH_pdb[3] = { 10, 22, 23 };
	const char *at_NTH_pdb[3] = { "CA", "CA", "CA" };

	pcopy = (double *)realloc( pcopy, sizeof(double) * 3 * nNTH );

	for( int t = 0; t < nNTH; t++ )
	{
		// copy the coordinates
		pcopy[3*t+0] = NTH[t].x;
		pcopy[3*t+1] = NTH[t].y;
		pcopy[3*t+2] = NTH[t].z;
		// look for the sites
		for( int x = 0; x < 3; x++ )
		{
			if( res_NTH_pdb[x] == NTH[t].res && !strcasecmp( NTH[t].atname, at_NTH_pdb[x] ) )
				set_NTH_pdb[x] = t;
		}
	}	
	
	if( set_NTH_pdb[0] < 0 || set_NTH_pdb[1] < 0 || set_NTH_pdb[2] < 0 )
	{
		printf("Failed to find site connections for SYT7.\n");
		exit(1);
	}
	
	alignStructuresOnAtomSet( NTH_align, set_NTH, pcopy, set_NTH_pdb, 3, nNTH ); 
	
	for( int t = 0; t < nNTH; t++ )
	{
		// copy the coordinates out
		NTH[t].x=pcopy[3*t+0];
		NTH[t].y=pcopy[3*t+1];
		NTH[t].z=pcopy[3*t+2];
	}
#endif

	(*at_out) = (struct atom_rec *)realloc( *at_out, sizeof(struct atom_rec) * (*nat_out + nC2A + nC2B + nNTH ) );

	for( int a = 0; a < nC2A; a++ )
	{
		(*at_out)[*nat_out] = C2A[a];
		(*nat_out) += 1;
	}
	
	for( int a = 0; a < nC2B; a++ )
	{
		(*at_out)[*nat_out] = C2B[a];
		(*nat_out) += 1;
	}

#ifdef DO_NTH	
	for( int a = 0; a < nNTH; a++ )
	{
		(*at_out)[*nat_out] = NTH[a];
		(*nat_out) += 1;
	}
#endif	
	free(pcopy);	
}


