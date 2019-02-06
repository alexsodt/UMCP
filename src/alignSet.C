// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "local_config.h"
#include "lapack_we_use.h"


struct quat
{
	public:
	double w;
	double x,y,z;

	quat() {};

	quat( double w_in, double *vec )
	{
		w = w_in;
		x = vec[0];
		y = vec[1];
		z = vec[2];
	}

};

	quat quat_invert( quat a )
	{
		quat t_quat;

		t_quat.w = a.w;
		t_quat.x = -a.x;
		t_quat.y = -a.y;
		t_quat.z = -a.z;

		return t_quat;
	}

	quat quat_mult( quat a, quat b )
	{
		quat t_quat;

		t_quat.w = ( a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z);
		t_quat.x = ( a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y);
		t_quat.y = ( a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z);
		t_quat.z = ( a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x);

		return t_quat;
	}

	quat quat_imult( quat a, quat b )
	{
		quat t_quat;

		t_quat.w = ( a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z);
		t_quat.x = ( a.w * b.x - a.x * b.w - a.y * b.z + a.z * b.y);
		t_quat.y = ( a.w * b.y - a.y * b.w - a.z * b.x + a.x * b.z);
		t_quat.z = ( a.w * b.z - a.z * b.w - a.x * b.y + a.y * b.x);

		return t_quat;
	}


double alignStructuresOnAtomSet( double *xyz1_in, 
			int *atoms_align1,
			double *xyz2_in,
			int *atoms_align2,
			int nat,
			int nat_tot )
{	
	double *xyz2 = (double *)malloc( sizeof(double) * 3 * nat );
	double *xyz1 = (double *)malloc( sizeof(double) * 3 * nat );
	
	for( int i = 0; i < nat; i++ )
	{
		xyz1[3*i+0] = xyz1_in[3*atoms_align1[i]+0];
		xyz1[3*i+1] = xyz1_in[3*atoms_align1[i]+1];
		xyz1[3*i+2] = xyz1_in[3*atoms_align1[i]+2];
		
		xyz2[3*i+0] = xyz2_in[3*atoms_align2[i]+0];
		xyz2[3*i+1] = xyz2_in[3*atoms_align2[i]+1];
		xyz2[3*i+2] = xyz2_in[3*atoms_align2[i]+2];
	}

	// subtract off center of mass of each xyz.

	double cm1[3] = {0,0,0};

	for( int a = 0; a < nat; a++ )
	{
		cm1[0] += xyz1[3*a+0];
		cm1[1] += xyz1[3*a+1];
		cm1[2] += xyz1[3*a+2];
	}

	cm1[0] /= nat;
	cm1[1] /= nat;
	cm1[2] /= nat;

	for( int a = 0; a < nat; a++ )
	{
		xyz1[3*a+0] -= cm1[0];
		xyz1[3*a+1] -= cm1[1];
		xyz1[3*a+2] -= cm1[2];
	}

	double cm2[3] = {0,0,0};

	for( int a = 0; a < nat; a++ )
	{
		cm2[0] += xyz2[3*a+0];
		cm2[1] += xyz2[3*a+1];
		cm2[2] += xyz2[3*a+2];
	}

	cm2[0] /= nat;
	cm2[1] /= nat;
	cm2[2] /= nat;

	for( int a = 0; a < nat; a++ )
	{
		xyz2[3*a+0] -= cm2[0];
		xyz2[3*a+1] -= cm2[1];
		xyz2[3*a+2] -= cm2[2];
	}

/*
	for( int a = 0; a < nat; a++ )
		printf("C %lf %lf %lf\n", xyz1[3*a+0], xyz1[3*a+1], xyz1[3*a+2] );

	for( int a = 0; a < nat; a++ )
		printf("C %lf %lf %lf\n", xyz2[3*a+0], xyz2[3*a+1], xyz2[3*a+2] );

	exit(1);
*/
	double dcm[3] = { cm2[0]-cm1[0], cm2[1]-cm1[1], cm2[2]-cm1[2] };

	double S_xx = 0, S_xy = 0, S_xz = 0;	
	double S_yx = 0, S_yy = 0, S_yz = 0;	
	double S_zx = 0, S_zy = 0, S_zz = 0;	

	double G_A=0, G_B=0;

	for( int a = 0; a < nat; a++ )
	{
		G_A += xyz1[3*a+0] * xyz1[3*a+0] + xyz1[3*a+1] * xyz1[3*a+1] + xyz1[3*a+2] * xyz1[3*a+2];
		G_B += xyz2[3*a+0] * xyz2[3*a+0] + xyz2[3*a+1] * xyz2[3*a+1] + xyz2[3*a+2] * xyz2[3*a+2];

		S_xx += xyz1[3*a+0] * xyz2[3*a+0];
		S_xy += xyz1[3*a+0] * xyz2[3*a+1];
		S_xz += xyz1[3*a+0] * xyz2[3*a+2];

		S_yx += xyz1[3*a+1] * xyz2[3*a+0];
		S_yy += xyz1[3*a+1] * xyz2[3*a+1];
		S_yz += xyz1[3*a+1] * xyz2[3*a+2];

		S_zx += xyz1[3*a+2] * xyz2[3*a+0];
		S_zy += xyz1[3*a+2] * xyz2[3*a+1];
		S_zz += xyz1[3*a+2] * xyz2[3*a+2];
	}

	double D = pow(S_xy * S_xy + S_xz * S_xz - S_yx * S_yx - S_zx * S_zx, 2.0);
	
	double E = ( - S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + S_yz * S_yz + S_zy * S_zy - 2 * ( S_yy * S_zz - S_yz * S_zy)) *
			( - S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + S_yz * S_yz + S_zy * S_zy + 2 * ( S_yy * S_zz - S_yz * S_zy));
	
	double F = ( -(S_xz + S_zx)*(S_yz - S_zy) + (S_xy - S_yx)*(S_xx - S_yy - S_zz)) *
			( -(S_xz - S_zx)*(S_yz + S_zy) + (S_xy - S_yx)*(S_xx - S_yy + S_zz));	

	double G = ( -(S_xz + S_zx)*(S_yz + S_zy) - (S_xy + S_yx)*(S_xx + S_yy - S_zz)) *
			( -(S_xz - S_zx)*(S_yz - S_zy) - (S_xy + S_yx)*(S_xx + S_yy + S_zz));

	double H = ((S_xy + S_yx) * (S_yz + S_zy) + (S_xz + S_zx)*(S_xx - S_yy + S_zz)) *
			(-(S_xy - S_yx) * (S_yz - S_zy) + (S_xz + S_zx)*(S_xx + S_yy + S_zz));

	double I = ((S_xy + S_yx)*(S_yz - S_zy) + (S_xz - S_zx) * (S_xx - S_yy - S_zz)) *
			(-(S_xy - S_yx)*(S_yz + S_zy) + (S_xz - S_zx) * (S_xx + S_yy - S_zz));

	double C2 = -2 *(S_xx * S_xx + S_xy * S_xy + S_xz * S_xz + S_yx * S_yx + S_yy * S_yy + S_yz * S_yz + S_zx * S_zx + S_zy * S_zy + S_zz * S_zz );
	
	double C1 = 8 * ( S_xx * S_yz * S_zy + S_yy * S_zx * S_xz + S_zz * S_xy * S_yx)
			- 8 *( S_xx * S_yy * S_zz + S_yz * S_zx * S_xy + S_zy * S_yx * S_xz);
	double C0 = D + E + F +G + H + I;

	double lam = 0.5 * (G_A + G_B);

	double lam_old;

//	printf("D: %lf E: %lf F: %lf G: %lf H: %lf I: %lf C0: %.12le C1: %.12le C2: %.12le\n",
//		D,E,F,G,H,I,C0,C1,C2 );

	int niter = 0;

	do {
		lam_old = lam;

		double P_lam = C0 + lam * (C1 + lam * ( C2 + lam * lam));
		double P_lam_d = C1 + lam * ( 2 * C2 + 4 * lam * lam );

		lam = lam - P_lam / P_lam_d; 	
		niter++;
		if( niter > 1000 )
		{
			printf("niter: %d %.12le %.12le.\n", niter, lam, lam_old );	
			exit(1);
		}
	} while( fabs(lam -lam_old) > 1e-7 );
	

//	printf("G_A + G_B: %.12le\n", G_A + G_B );

//	printf("RMSD: %lf\n", (G_A + G_B - 2 * lam)/nat );


	double *K = (double *)malloc( sizeof(double) * 4 * 4 );

	K[0*4+0] = (S_xx + S_yy + S_zz);	K[0*4+1] = (S_yz - S_zy);	K[0*4+2] = (S_zx - S_xz);	K[0*4+3] = (S_xy - S_yx);
	K[1*4+0] = (S_yz - S_zy);		K[1*4+1] = (S_xx-S_yy-S_zz);	K[1*4+2] = (S_xy + S_yx);	K[1*4+3] = (S_xz + S_xz);
	K[2*4+0] = (S_zx - S_xz);		K[2*4+1] = (S_xy + S_yx);	K[2*4+2] = (-S_xx+S_yy-S_zz);	K[2*4+3] = (S_yz + S_zy);
	K[3*4+0] = (S_xy - S_yx);		K[3*4+1] = (S_zx+S_xz);		K[3*4+2] = (S_yz + S_zy);	K[3*4+3] = (-S_xx - S_yy + S_zz);	

	char jobz = 'V';
	char uplo = 'U';
	int order =4;
	double *ev = (double *)malloc( sizeof(double) * 4 );	
	double owork = 0;
	int lwork = -1;	
	int info = 0;

	dsyev(&jobz,&uplo,&order,K,&order,ev,&owork,&lwork,&info);
	double *work = (double *)malloc( sizeof(int)*(int)lround(owork));

	lwork = lround(owork);
	dsyev(&jobz,&uplo,&order,K,&order,ev,work,&lwork,&info);
	double qscale =  K[3*4+0];
	double q_vec[3] = { K[3*4+1], K[3*4+2], K[3*4+3] };
	
	quat rot_quat( qscale, q_vec);
	quat rot_quat_i;

	rot_quat_i = quat_invert(rot_quat);

	for( int i = 0; i < nat_tot; i++ )
	{
		double ovec[3] = {xyz2_in[3*i+0] - cm2[0], xyz2_in[3*i+1] - cm2[1], xyz2_in[3*i+2] - cm2[2]};

		double nvec[3] = {0,0,0};
		
		quat o_quat( 0, ovec );

		quat prod_quat = quat_imult( rot_quat, quat_mult( o_quat, rot_quat) );

		xyz2_in[3*i+0] = prod_quat.x + cm1[0];
		xyz2_in[3*i+1] = prod_quat.y + cm1[1];
		xyz2_in[3*i+2] = prod_quat.z + cm1[2];
	}

	return (G_A + G_B - 2 * lam)/nat;
}

void displacePlanes( double *xyz_to_displace, int nat, double fac )
{
	double d1[3] = { 
			xyz_to_displace[3*0+0] - xyz_to_displace[3*1+0],
			xyz_to_displace[3*0+1] - xyz_to_displace[3*1+1],
			xyz_to_displace[3*0+2] - xyz_to_displace[3*1+2] };
	
	double d2[3] = { 
			xyz_to_displace[3*2+0] - xyz_to_displace[3*1+0],
			xyz_to_displace[3*2+1] - xyz_to_displace[3*1+1],
			xyz_to_displace[3*2+2] - xyz_to_displace[3*1+2] };
	
	double c1[3] = {  (d1[1]*d2[2] - d1[2] * d2[1]),
			 -(d1[0]*d2[2] - d1[2] * d2[0]),
			  (d1[0]*d2[1] - d1[1] * d2[0]) };

	double lc = sqrt(c1[0]*c1[0]+c1[1]*c1[1]+c1[2]*c1[2]);
		
	c1[0] /= lc;
	c1[1] /= lc;
	c1[2] /= lc;
		
	for( int i = 0; i < nat; i++ )
	{	
		xyz_to_displace[3*i+0] += c1[0] * fac;
		xyz_to_displace[3*i+1] += c1[1] * fac;
		xyz_to_displace[3*i+2] += c1[2] * fac;
	}
}


void rotatePlanar( double *xyz_rot, double pt[3], int nat, double val )
{
	double d1[3] = { 
			xyz_rot[3*0+0] - xyz_rot[3*1+0],
			xyz_rot[3*0+1] - xyz_rot[3*1+1],
			xyz_rot[3*0+2] - xyz_rot[3*1+2] };
	
	double d2[3] = { 
			xyz_rot[3*2+0] - xyz_rot[3*1+0],
			xyz_rot[3*2+1] - xyz_rot[3*1+1],
			xyz_rot[3*2+2] - xyz_rot[3*1+2] };
	
	double c1[3] = {  (d1[1]*d2[2] - d1[2] * d2[1]),
			 -(d1[0]*d2[2] - d1[2] * d2[0]),
			  (d1[0]*d2[1] - d1[1] * d2[0]) };

	double lc = sqrt(c1[0]*c1[0]+c1[1]*c1[1]+c1[2]*c1[2]);
		
	c1[0] /= lc;
	c1[1] /= lc;
	c1[2] /= lc;
	
	double qscale =  val;
	double q_vec[3] = { c1[0], c1[1], c1[2] };
	
	quat rot_quat( qscale, q_vec);
	quat rot_quat_i;

	rot_quat_i = quat_invert(rot_quat);

	for( int i = 0; i < nat; i++ )
	{	
		xyz_rot[3*i+0] -= pt[0];
		xyz_rot[3*i+1] -= pt[1];
		xyz_rot[3*i+2] -= pt[2];

		double uvcross[3] = { (q_vec[1] * xyz_rot[3*i+2] - q_vec[2] * xyz_rot[3*i+1]),
				     -(q_vec[0] * xyz_rot[3*i+2] - q_vec[2] * xyz_rot[3*i+0]),
				      (q_vec[0] * xyz_rot[3*i+1] - q_vec[1] * xyz_rot[3*i+0]) };
		double udu = q_vec[0] * q_vec[0] + q_vec[1] * q_vec[1] + q_vec[2] * q_vec[2];
		double udv = q_vec[0] * xyz_rot[3*i+0] + q_vec[1] * xyz_rot[3*i+1] + q_vec[2] * xyz_rot[3*i+2];
		double cha = cos( val / (2) );		
		double sha = sin( val / (2) );		

		double rvec[3] = {
					xyz_rot[3*i+0] * (cha*cha - sha*sha) + uvcross[0] * (2 * sha * cha ) + q_vec[0] * (udv) * 2 *sha*sha,
					xyz_rot[3*i+1] * (cha*cha - sha*sha) + uvcross[1] * (2 * sha * cha ) + q_vec[1] * (udv) * 2 *sha*sha,
					xyz_rot[3*i+2] * (cha*cha - sha*sha) + uvcross[2] * (2 * sha * cha ) + q_vec[2] * (udv) * 2 *sha*sha };

	
/*		quat o_quat( 0, xyz_rot+3*i );

		quat prod_quat = quat_imult( rot_quat, quat_mult( o_quat, rot_quat) );
*/
		xyz_rot[3*i+0] = rvec[0] + pt[0];
		xyz_rot[3*i+1] = rvec[1] + pt[1];
		xyz_rot[3*i+2] = rvec[2] + pt[2];
	}
}

void flipCoords( double *xyz, int nat, double *porig, double *n )
{
	double val = porig[0]*n[0] + porig[1] * n[1] + porig[2] * n[2];

	double nn = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

	for( int a = 0; a < nat; a++ )
	{
		double px[3] = {xyz[3*a+0], xyz[3*a+1], xyz[3*a+2] };
		
		double np = n[0] * px[0] + n[1] * px[1] + n[2] * px[2];

		double lam = -2*(val - np) / (nn);
		
		double xx[3] = { (1.0/2.0) * (-lam * n[0] + 2 * px[0]),
				 (1.0/2.0) * (-lam * n[1] + 2 * px[1]), 
				 (1.0/2.0) * (-lam * n[2] + 2 * px[2]) };

		xyz[3*a+0] = xx[0] - (px[0]-xx[0]); 
		xyz[3*a+1] = xx[1] - (px[1]-xx[1]); 
		xyz[3*a+2] = xx[2] - (px[2]-xx[2]); 
	}
}

void rotateAxial( double *xyz_rot, int at1, int at2, int nat, double val )
{
	double c1[3] = { xyz_rot[3*at1+0] - xyz_rot[3*at2+0],
			 xyz_rot[3*at1+1] - xyz_rot[3*at2+1],
			 xyz_rot[3*at1+2] - xyz_rot[3*at2+2] };

	double lc = sqrt(c1[0]*c1[0]+c1[1]*c1[1]+c1[2]*c1[2]);
		
	c1[0] /= lc;
	c1[1] /= lc;
	c1[2] /= lc;
	
	double qscale =  val;
	double q_vec[3] = { c1[0], c1[1], c1[2] };
	double pt[3] = { xyz_rot[3*at2+0], xyz_rot[3*at2+1], xyz_rot[3*at2+2] };	

	quat rot_quat( qscale, q_vec);
	quat rot_quat_i;

	rot_quat_i = quat_invert(rot_quat);

	for( int i = 0; i < nat; i++ )
	{	
		xyz_rot[3*i+0] -= pt[0];
		xyz_rot[3*i+1] -= pt[1];
		xyz_rot[3*i+2] -= pt[2];

		double uvcross[3] = { (q_vec[1] * xyz_rot[3*i+2] - q_vec[2] * xyz_rot[3*i+1]),
				     -(q_vec[0] * xyz_rot[3*i+2] - q_vec[2] * xyz_rot[3*i+0]),
				      (q_vec[0] * xyz_rot[3*i+1] - q_vec[1] * xyz_rot[3*i+0]) };
		double udu = q_vec[0] * q_vec[0] + q_vec[1] * q_vec[1] + q_vec[2] * q_vec[2];
		double udv = q_vec[0] * xyz_rot[3*i+0] + q_vec[1] * xyz_rot[3*i+1] + q_vec[2] * xyz_rot[3*i+2];
		double cha = cos( val / (2) );		
		double sha = sin( val / (2) );		

		double rvec[3] = {
					xyz_rot[3*i+0] * (cha*cha - sha*sha) + uvcross[0] * (2 * sha * cha ) + q_vec[0] * (udv) * 2 *sha*sha,
					xyz_rot[3*i+1] * (cha*cha - sha*sha) + uvcross[1] * (2 * sha * cha ) + q_vec[1] * (udv) * 2 *sha*sha,
					xyz_rot[3*i+2] * (cha*cha - sha*sha) + uvcross[2] * (2 * sha * cha ) + q_vec[2] * (udv) * 2 *sha*sha };

	
/*		quat o_quat( 0, xyz_rot+3*i );

		quat prod_quat = quat_imult( rot_quat, quat_mult( o_quat, rot_quat) );
*/
		xyz_rot[3*i+0] = rvec[0] + pt[0];
		xyz_rot[3*i+1] = rvec[1] + pt[1];
		xyz_rot[3*i+2] = rvec[2] + pt[2];
	}
}






void rotateArbitrary( double *xyz_rot, double *axis, double *origin, int nat, double val )
{
	double c1[3] = { axis[0], axis[1], axis[2] };

	double lc = sqrt(c1[0]*c1[0]+c1[1]*c1[1]+c1[2]*c1[2]);
		
	c1[0] /= lc;
	c1[1] /= lc;
	c1[2] /= lc;
	
	double qscale =  val;
	double q_vec[3] = { c1[0], c1[1], c1[2] };
	double pt[3] = { origin[0], origin[1], origin[2] };	

	quat rot_quat( qscale, q_vec);
	quat rot_quat_i;

	rot_quat_i = quat_invert(rot_quat);
		
	double cha = cos( val / (2) );		
	double sha = sin( val / (2) );		

	for( int i = 0; i < nat; i++ )
	{	
		xyz_rot[3*i+0] -= pt[0];
		xyz_rot[3*i+1] -= pt[1];
		xyz_rot[3*i+2] -= pt[2];

		double uvcross[3] = { (q_vec[1] * xyz_rot[3*i+2] - q_vec[2] * xyz_rot[3*i+1]),
				     -(q_vec[0] * xyz_rot[3*i+2] - q_vec[2] * xyz_rot[3*i+0]),
				      (q_vec[0] * xyz_rot[3*i+1] - q_vec[1] * xyz_rot[3*i+0]) };
		double udu = q_vec[0] * q_vec[0] + q_vec[1] * q_vec[1] + q_vec[2] * q_vec[2];
		double udv = q_vec[0] * xyz_rot[3*i+0] + q_vec[1] * xyz_rot[3*i+1] + q_vec[2] * xyz_rot[3*i+2];

		double rvec[3] = {
					xyz_rot[3*i+0] * (cha*cha - sha*sha) + uvcross[0] * (2 * sha * cha ) + q_vec[0] * (udv) * 2 *sha*sha,
					xyz_rot[3*i+1] * (cha*cha - sha*sha) + uvcross[1] * (2 * sha * cha ) + q_vec[1] * (udv) * 2 *sha*sha,
					xyz_rot[3*i+2] * (cha*cha - sha*sha) + uvcross[2] * (2 * sha * cha ) + q_vec[2] * (udv) * 2 *sha*sha };

	
/*		quat o_quat( 0, xyz_rot+3*i );

		quat prod_quat = quat_imult( rot_quat, quat_mult( o_quat, rot_quat) );
*/
		xyz_rot[3*i+0] = rvec[0] + pt[0];
		xyz_rot[3*i+1] = rvec[1] + pt[1];
		xyz_rot[3*i+2] = rvec[2] + pt[2];
	}
}

double getChi2( double *xyz1, 
			double *xyz2,
			int nat )
{	

	// subtract off center of mass of each xyz.

	double cm1[3] = {0,0,0};

	for( int a = 0; a < nat; a++ )
	{
		cm1[0] += xyz1[3*a+0];
		cm1[1] += xyz1[3*a+1];
		cm1[2] += xyz1[3*a+2];
	}

	cm1[0] /= nat;
	cm1[1] /= nat;
	cm1[2] /= nat;

	for( int a = 0; a < nat; a++ )
	{
		xyz1[3*a+0] -= cm1[0];
		xyz1[3*a+1] -= cm1[1];
		xyz1[3*a+2] -= cm1[2];
	}

	double cm2[3] = {0,0,0};

	for( int a = 0; a < nat; a++ )
	{
		cm2[0] += xyz2[3*a+0];
		cm2[1] += xyz2[3*a+1];
		cm2[2] += xyz2[3*a+2];
	}

	cm2[0] /= nat;
	cm2[1] /= nat;
	cm2[2] /= nat;

	for( int a = 0; a < nat; a++ )
	{
		xyz2[3*a+0] -= cm2[0];
		xyz2[3*a+1] -= cm2[1];
		xyz2[3*a+2] -= cm2[2];
	}

	double dcm[3] = { cm2[0]-cm1[0], cm2[1]-cm1[1], cm2[2]-cm1[2] };

	double S_xx = 0, S_xy = 0, S_xz = 0;	
	double S_yx = 0, S_yy = 0, S_yz = 0;	
	double S_zx = 0, S_zy = 0, S_zz = 0;	

	double G_A=0, G_B=0;

	for( int a = 0; a < nat; a++ )
	{
		G_A += xyz1[3*a+0] * xyz1[3*a+0] + xyz1[3*a+1] * xyz1[3*a+1] + xyz1[3*a+2] * xyz1[3*a+2];
		G_B += xyz2[3*a+0] * xyz2[3*a+0] + xyz2[3*a+1] * xyz2[3*a+1] + xyz2[3*a+2] * xyz2[3*a+2];

		S_xx += xyz1[3*a+0] * xyz2[3*a+0];
		S_xy += xyz1[3*a+0] * xyz2[3*a+1];
		S_xz += xyz1[3*a+0] * xyz2[3*a+2];

		S_yx += xyz1[3*a+1] * xyz2[3*a+0];
		S_yy += xyz1[3*a+1] * xyz2[3*a+1];
		S_yz += xyz1[3*a+1] * xyz2[3*a+2];

		S_zx += xyz1[3*a+2] * xyz2[3*a+0];
		S_zy += xyz1[3*a+2] * xyz2[3*a+1];
		S_zz += xyz1[3*a+2] * xyz2[3*a+2];
	}

	double D = pow(S_xy * S_xy + S_xz * S_xz - S_yx * S_yx - S_zx * S_zx, 2.0);
	
	double E = ( - S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + S_yz * S_yz + S_zy * S_zy - 2 * ( S_yy * S_zz - S_yz * S_zy)) *
			( - S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + S_yz * S_yz + S_zy * S_zy + 2 * ( S_yy * S_zz - S_yz * S_zy));
	
	double F = ( -(S_xz + S_zx)*(S_yz - S_zy) + (S_xy - S_yx)*(S_xx - S_yy - S_zz)) *
			( -(S_xz - S_zx)*(S_yz + S_zy) + (S_xy - S_yx)*(S_xx - S_yy + S_zz));	

	double G = ( -(S_xz + S_zx)*(S_yz + S_zy) - (S_xy + S_yx)*(S_xx + S_yy - S_zz)) *
			( -(S_xz - S_zx)*(S_yz - S_zy) - (S_xy + S_yx)*(S_xx + S_yy + S_zz));

	double H = ((S_xy + S_yx) * (S_yz + S_zy) + (S_xz + S_zx)*(S_xx - S_yy + S_zz)) *
			(-(S_xy - S_yx) * (S_yz - S_zy) + (S_xz + S_zx)*(S_xx + S_yy + S_zz));

	double I = ((S_xy + S_yx)*(S_yz - S_zy) + (S_xz - S_zx) * (S_xx - S_yy - S_zz)) *
			(-(S_xy - S_yx)*(S_yz + S_zy) + (S_xz - S_zx) * (S_xx + S_yy - S_zz));

	double C2 = -2 *(S_xx * S_xx + S_xy * S_xy + S_xz * S_xz + S_yx * S_yx + S_yy * S_yy + S_yz * S_yz + S_zx * S_zx + S_zy * S_zy + S_zz * S_zz );
	
	double C1 = 8 * ( S_xx * S_yz * S_zy + S_yy * S_zx * S_xz + S_zz * S_xy * S_yx)
			- 8 *( S_xx * S_yy * S_zz + S_yz * S_zx * S_xy + S_zy * S_yx * S_xz);
	double C0 = D + E + F +G + H + I;

	double lam = 0.5 * (G_A + G_B);

	double lam_old;

//	printf("D: %lf E: %lf F: %lf G: %lf H: %lf I: %lf C0: %.12le C1: %.12le C2: %.12le\n",
//		D,E,F,G,H,I,C0,C1,C2 );

	int niter = 0;

	do {
		lam_old = lam;

		double P_lam = C0 + lam * (C1 + lam * ( C2 + lam * lam));
		double P_lam_d = C1 + lam * ( 2 * C2 + 4 * lam * lam );

		lam = lam - P_lam / P_lam_d; 	
		niter++;
		if( niter > 1000 )
		{
			printf("niter: %d %.12le %.12le.\n", niter, lam, lam_old );	
			exit(1);
		}
	} while( fabs(lam -lam_old) > 1e-7 );
	
	return (G_A + G_B - 2 * lam)/nat;
}

double alignStructuresOnAtomSetXY( double *xyz1_in, 
			int *atoms_align1,
			double *xyz2_in,
			int *atoms_align2,
			int nat,
			int nat_tot )
{	
	double *xyz2 = (double *)malloc( sizeof(double) * 3 * nat );
	double *xyz1 = (double *)malloc( sizeof(double) * 3 * nat );
	
	for( int i = 0; i < nat; i++ )
	{
		xyz1[3*i+0] = xyz1_in[3*atoms_align1[i]+0];
		xyz1[3*i+1] = xyz1_in[3*atoms_align1[i]+1];
		xyz1[3*i+2] = xyz1_in[3*atoms_align1[i]+2];
		
		xyz2[3*i+0] = xyz2_in[3*atoms_align2[i]+0];
		xyz2[3*i+1] = xyz2_in[3*atoms_align2[i]+1];
		xyz2[3*i+2] = xyz2_in[3*atoms_align2[i]+2];
	}

	// subtract off center of mass of each xyz.

	double cm1[3] = {0,0,0};

	for( int a = 0; a < nat; a++ )
	{
		cm1[0] += xyz1[3*a+0];
		cm1[1] += xyz1[3*a+1];
		cm1[2] += xyz1[3*a+2];
	}

	cm1[0] /= nat;
	cm1[1] /= nat;
	cm1[2] /= nat;

	for( int a = 0; a < nat; a++ )
	{
		xyz1[3*a+0] -= cm1[0];
		xyz1[3*a+1] -= cm1[1];
		xyz1[3*a+2] -= cm1[2];
	}

	double cm2[3] = {0,0,0};

	for( int a = 0; a < nat; a++ )
	{
		cm2[0] += xyz2[3*a+0];
		cm2[1] += xyz2[3*a+1];
		cm2[2] += xyz2[3*a+2];
	}

	cm2[0] /= nat;
	cm2[1] /= nat;
	cm2[2] /= nat;

	for( int a = 0; a < nat; a++ )
	{
		xyz2[3*a+0] -= cm2[0];
		xyz2[3*a+1] -= cm2[1];
		xyz2[3*a+2] -= cm2[2];
	}

/*
	for( int a = 0; a < nat; a++ )
		printf("C %lf %lf %lf\n", xyz1[3*a+0], xyz1[3*a+1], xyz1[3*a+2] );

	for( int a = 0; a < nat; a++ )
		printf("C %lf %lf %lf\n", xyz2[3*a+0], xyz2[3*a+1], xyz2[3*a+2] );

	exit(1);
*/
	double dcm[3] = { cm2[0]-cm1[0], cm2[1]-cm1[1], cm2[2]-cm1[2] };

	double S_xx = 0, S_xy = 0, S_xz = 0;	
	double S_yx = 0, S_yy = 0, S_yz = 0;	
	double S_zx = 0, S_zy = 0, S_zz = 0;	

	double G_A=0, G_B=0;

	for( int a = 0; a < nat; a++ )
	{
		G_A += xyz1[3*a+0] * xyz1[3*a+0] + xyz1[3*a+1] * xyz1[3*a+1];// + xyz1[3*a+2] * xyz1[3*a+2];
		G_B += xyz2[3*a+0] * xyz2[3*a+0] + xyz2[3*a+1] * xyz2[3*a+1];// + xyz2[3*a+2] * xyz2[3*a+2];

		S_xx += xyz1[3*a+0] * xyz2[3*a+0];
		S_xy += xyz1[3*a+0] * xyz2[3*a+1];
//		S_xz += xyz1[3*a+0] * xyz2[3*a+2];

		S_yx += xyz1[3*a+1] * xyz2[3*a+0];
		S_yy += xyz1[3*a+1] * xyz2[3*a+1];
//		S_yz += xyz1[3*a+1] * xyz2[3*a+2];

//		S_zx += xyz1[3*a+2] * xyz2[3*a+0];
//		S_zy += xyz1[3*a+2] * xyz2[3*a+1];
//		S_zz += xyz1[3*a+2] * xyz2[3*a+2];
	}

	double D = pow(S_xy * S_xy + S_xz * S_xz - S_yx * S_yx - S_zx * S_zx, 2.0);
	
	double E = ( - S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + S_yz * S_yz + S_zy * S_zy - 2 * ( S_yy * S_zz - S_yz * S_zy)) *
			( - S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + S_yz * S_yz + S_zy * S_zy + 2 * ( S_yy * S_zz - S_yz * S_zy));
	
	double F = ( -(S_xz + S_zx)*(S_yz - S_zy) + (S_xy - S_yx)*(S_xx - S_yy - S_zz)) *
			( -(S_xz - S_zx)*(S_yz + S_zy) + (S_xy - S_yx)*(S_xx - S_yy + S_zz));	

	double G = ( -(S_xz + S_zx)*(S_yz + S_zy) - (S_xy + S_yx)*(S_xx + S_yy - S_zz)) *
			( -(S_xz - S_zx)*(S_yz - S_zy) - (S_xy + S_yx)*(S_xx + S_yy + S_zz));

	double H = ((S_xy + S_yx) * (S_yz + S_zy) + (S_xz + S_zx)*(S_xx - S_yy + S_zz)) *
			(-(S_xy - S_yx) * (S_yz - S_zy) + (S_xz + S_zx)*(S_xx + S_yy + S_zz));

	double I = ((S_xy + S_yx)*(S_yz - S_zy) + (S_xz - S_zx) * (S_xx - S_yy - S_zz)) *
			(-(S_xy - S_yx)*(S_yz + S_zy) + (S_xz - S_zx) * (S_xx + S_yy - S_zz));

	double C2 = -2 *(S_xx * S_xx + S_xy * S_xy + S_xz * S_xz + S_yx * S_yx + S_yy * S_yy + S_yz * S_yz + S_zx * S_zx + S_zy * S_zy + S_zz * S_zz );
	
	double C1 = 8 * ( S_xx * S_yz * S_zy + S_yy * S_zx * S_xz + S_zz * S_xy * S_yx)
			- 8 *( S_xx * S_yy * S_zz + S_yz * S_zx * S_xy + S_zy * S_yx * S_xz);
	double C0 = D + E + F +G + H + I;

	double lam = 0.5 * (G_A + G_B);

	double lam_old;

//	printf("D: %lf E: %lf F: %lf G: %lf H: %lf I: %lf C0: %.12le C1: %.12le C2: %.12le\n",
//		D,E,F,G,H,I,C0,C1,C2 );

	int niter = 0;

#ifndef DISABLE_ALIGNMENT
	do {
		lam_old = lam;

		double P_lam = C0 + lam * (C1 + lam * ( C2 + lam * lam));
		double P_lam_d = C1 + lam * ( 2 * C2 + 4 * lam * lam );

		lam = lam - P_lam / P_lam_d; 	
		niter++;
		if( niter > 10000 )
		{
			printf("niter: %d %.12le %.12le.\n", niter, lam, lam_old );	
			exit(1);
		}
	} while( fabs(lam -lam_old) > 5e-7 );
#endif	

//	printf("G_A + G_B: %.12le\n", G_A + G_B );

//	printf("RMSD: %lf\n", (G_A + G_B - 2 * lam)/nat );


	double *K = (double *)malloc( sizeof(double) * 4 * 4 );

	K[0*4+0] = (S_xx + S_yy + S_zz);	K[0*4+1] = (S_yz - S_zy);	K[0*4+2] = (S_zx - S_xz);	K[0*4+3] = (S_xy - S_yx);
	K[1*4+0] = (S_yz - S_zy);		K[1*4+1] = (S_xx-S_yy-S_zz);	K[1*4+2] = (S_xy + S_yx);	K[1*4+3] = (S_xz + S_xz);
	K[2*4+0] = (S_zx - S_xz);		K[2*4+1] = (S_xy + S_yx);	K[2*4+2] = (-S_xx+S_yy-S_zz);	K[2*4+3] = (S_yz + S_zy);
	K[3*4+0] = (S_xy - S_yx);		K[3*4+1] = (S_zx+S_xz);		K[3*4+2] = (S_yz + S_zy);	K[3*4+3] = (-S_xx - S_yy + S_zz);	

	char jobz = 'V';
	char uplo = 'U';
	int order =4;
	double *ev = (double *)malloc( sizeof(double) * 4 );	
	double owork = 0;
	int lwork = -1;	
	int info = 0;

	dsyev(&jobz,&uplo,&order,K,&order,ev,&owork,&lwork,&info);
	double *work = (double *)malloc( sizeof(int)*(int)lround(owork));

	lwork = lround(owork);
	dsyev(&jobz,&uplo,&order,K,&order,ev,work,&lwork,&info);
	double qscale =  K[3*4+0];
	double q_vec[3] = { K[3*4+1], K[3*4+2], K[3*4+3] };
	
	quat rot_quat( qscale, q_vec);
	quat rot_quat_i;

	rot_quat_i = quat_invert(rot_quat);

	for( int i = 0; i < nat_tot; i++ )
	{
		double ovec[3] = {xyz2_in[3*i+0] - cm2[0], xyz2_in[3*i+1] - cm2[1], xyz2_in[3*i+2] - cm2[2]};

		double nvec[3] = {0,0,0};
		
		quat o_quat( 0, ovec );

		quat prod_quat = quat_imult( rot_quat, quat_mult( o_quat, rot_quat) );

#ifdef DISABLE_ALIGNMENT
		xyz2_in[3*i+0] = ovec[0] + cm1[0];
		xyz2_in[3*i+1] = ovec[1] + cm1[1];
		xyz2_in[3*i+2] = ovec[2] + cm1[2];
#else
		xyz2_in[3*i+0] = prod_quat.x + cm1[0];
		xyz2_in[3*i+1] = prod_quat.y + cm1[1];
		xyz2_in[3*i+2] = ovec[2] + cm1[2];//prod_quat.z + cm1[2];
#endif
	}

	free(xyz1);
	free(xyz2);
	free(K);
	free(ev);
	free(work);

	return (G_A + G_B - 2 * lam)/nat;
}
