#include "spline.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "lapack_we_use.h"

void AddPointToSpecificMatrix( double x, double f, double *ADD_A, double *ADD_B, int ind );


double A( double x, double x0, double x1 )
{
	return (x1-x)/(x1-x0);
}

double dA( double x, double x0, double x1 )
{
	return -1/(x1-x0);
}

double ddA( double x, double x0, double x1 )
{
	return 0;
}


double B( double x, double x0, double x1 )
{
	return 1.0 - A(x,x0,x1);
}

double dB( double x, double x0, double x1 )
{
	return - dA(x,x0,x1);
}

double ddB( double x, double x0, double x1 )
{
	return 0;
}
	

double C( double x, double x0, double x1 )
{
	double tA = A(x,x0,x1);

	return (1.0/6.0)*(tA*tA*tA-tA)*(x1-x0)*(x1-x0);
}

double dC( double x, double x0, double x1 )
{
	double tA = A(x,x0,x1);
	double tdA = dA(x,x0,x1);

	return (1.0/6.0) * (3 * tA * tA * tdA - tdA) * (x1-x0)*(x1-x0);
}

double ddC( double x, double x0, double x1 )
{
	double tA = A(x,x0,x1);
	double tdA = dA(x,x0,x1);
	double tddA = ddA(x,x0,x1);

	return (1.0/6.0) * (6 * tA * tdA * tdA + 3 * tA * tA * tddA - tddA ) * (x1-x0) * (x1-x0); 
}


double D( double x, double x0, double x1 )
{
	double tB = B(x,x0,x1);

	return (1.0/6.0)*(tB*tB*tB-tB)*(x1-x0)*(x1-x0);
}

double dD( double x, double x0, double x1 )
{
	double tB = B(x,x0,x1);
	double tdB = dB(x,x0,x1); 

	return (1.0/6.0)*(3*tB*tB*tdB-tdB)*(x1-x0)*(x1-x0);
}

double ddD( double x, double x0, double x1 )
{
	double tB = B(x,x0,x1);
	double tdB = dB(x,x0,x1); 
	double tddB = ddB(x,x0,x1); 

	return (1.0/6.0)*(6*tB*tdB*tdB + 3 * tB * tB * tddB -tddB)*(x1-x0)*(x1-x0);
}

#define MAX_INDEX	100

static double CUR_XMIN[MAX_INDEX];
static double CUR_XMAX[MAX_INDEX];
static int    CUR_NBINS[MAX_INDEX];
static double *CUR_A[MAX_INDEX];
static double *CUR_B[MAX_INDEX];
static double *CUR_SOLUTION[MAX_INDEX];
static int *CUR_PTS_IN_BIN[MAX_INDEX];
static int nparams[MAX_INDEX];
static int nparamstot[MAX_INDEX];
static int nboundaries[MAX_INDEX];
static int nconstraints[MAX_INDEX];
static int is_setup[MAX_INDEX] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static int is_periodic[MAX_INDEX] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

void setupSpline( double XMIN, double XMAX, int NBINS, int ind, int periodic)
{
	is_periodic[ind] = periodic;
	is_setup[ind] = 1;
	// two parameters at each bin.  value, derivative.

	if( is_periodic[ind] )
		nboundaries[ind]  = NBINS;
	else
		nboundaries[ind]  = NBINS+1;

	nconstraints[ind] = (nboundaries[ind]-2) + 2;

	nparams[ind] = 2 * nboundaries[ind];
	nparamstot[ind] = nparams[ind] + nconstraints[ind];	

	CUR_XMIN[ind] = XMIN;
	CUR_XMAX[ind] = XMAX;
	CUR_NBINS[ind] = NBINS;

	CUR_A[ind] = (double *)malloc( sizeof(double) * nparamstot[ind] * nparamstot[ind] );
	memset( CUR_A[ind], 0, sizeof(double) * nparamstot[ind] * nparamstot[ind] );
	CUR_B[ind] = (double *)malloc( sizeof(double) * nparamstot[ind] );
	memset( CUR_B[ind], 0, sizeof(double) * nparamstot[ind] );

	// set up the constraints on A.
	CUR_SOLUTION[ind] = (double *)malloc( sizeof(double) * nparamstot[ind] );
	memset(CUR_SOLUTION[ind], 0, sizeof(double) * nparamstot[ind] );
	CUR_PTS_IN_BIN[ind] = (int *)malloc( sizeof(int) * NBINS );
	memset( CUR_PTS_IN_BIN[ind], 0, sizeof(int) * NBINS );
	double dx = (XMAX-XMIN) / (double)NBINS;

	for( int c = 0; c < nconstraints[ind]; c++ )
	{
		int cp1 = c+1;
		int cp2 = c+2;

		if( !is_periodic[ind] && c >= nconstraints[ind]-2 )
			continue;
		if( is_periodic[ind] )
		{
			if( cp1 >= nconstraints[ind] )
				cp1 -= nconstraints[ind];
			if( cp2 >= nconstraints[ind] )
				cp2 -= nconstraints[ind];
		}	
		int p_y_m1   = c*2+0;
		int p_ydd_m1 = c*2+1;
		int p_y      = cp1*2+0;
		int p_ydd    = cp1*2+1;
		int p_y_p1   = cp2*2+0;
		int p_ydd_p1 = cp2*2+1;

		CUR_A[ind][(nparams[ind]+c)*nparamstot[ind]+p_y_m1] = -1.0 / (dx);
		CUR_A[ind][(nparams[ind]+c)*nparamstot[ind]+p_y]    =  2.0 / (dx);
		CUR_A[ind][(nparams[ind]+c)*nparamstot[ind]+p_y_p1] = -1.0 / (dx);
		
		CUR_A[ind][(nparams[ind]+c)+nparamstot[ind]*p_y_m1] = -1.0 / (dx);
		CUR_A[ind][(nparams[ind]+c)+nparamstot[ind]*p_y]    =  2.0 / (dx);
		CUR_A[ind][(nparams[ind]+c)+nparamstot[ind]*p_y_p1] = -1.0 / (dx);

		CUR_A[ind][(nparams[ind]+c)*nparamstot[ind]+p_ydd_m1] = (dx) / 6.0;	
		CUR_A[ind][(nparams[ind]+c)*nparamstot[ind]+p_ydd]    = 2*(dx) / 3.0;	
		CUR_A[ind][(nparams[ind]+c)*nparamstot[ind]+p_ydd_p1] = (dx) / 6.0;	
		
		CUR_A[ind][(nparams[ind]+c)+nparamstot[ind]*p_ydd_m1] = (dx) / 6.0;	
		CUR_A[ind][(nparams[ind]+c)+nparamstot[ind]*p_ydd]    = 2*(dx) / 3.0;	
		CUR_A[ind][(nparams[ind]+c)+nparamstot[ind]*p_ydd_p1] = (dx) / 6.0;	
	}

	if( !is_periodic[ind] )
	{
		int c_nat_left = nconstraints[ind]-2;
		int c_nat_right = nconstraints[ind]-1;

		CUR_A[ind][(nparams[ind]+c_nat_right)*nparamstot[ind]+(nparams[ind]-1)] = 1;
		CUR_A[ind][(nparams[ind]+c_nat_left)*nparamstot[ind]+1]            = 1;
		CUR_A[ind][(nparams[ind]+c_nat_right)+nparamstot[ind]*(nparams[ind]-1)] = 1;
		CUR_A[ind][(nparams[ind]+c_nat_left)+nparamstot[ind]*1]            = 1;
	}}
 

void AddPointToSpline( double x, double f, int ind)
{
	if( is_periodic[ind] )
	{
		while( x < CUR_XMIN[ind] ) x += (CUR_XMAX[ind]-CUR_XMIN[ind]);
		while( x >= CUR_XMAX[ind] ) x -= (CUR_XMAX[ind]-CUR_XMIN[ind]);
	}

	AddPointToSpecificMatrix( x, f, CUR_A[ind], CUR_B[ind], ind);
}

void AddPointToSpecificMatrix( double x, double f, double *ADD_A, double *ADD_B, int ind)
{
	int border_left  = CUR_NBINS[ind] * (x - CUR_XMIN[ind]) / (CUR_XMAX[ind]-CUR_XMIN[ind]);
	int border_right = border_left + 1;

	double dx = (CUR_XMAX[ind]-CUR_XMIN[ind])/CUR_NBINS[ind];

	double x0 = CUR_XMIN[ind] + border_left * dx;
	double x1 = CUR_XMIN[ind] + border_right * dx;
	
	if( is_periodic[ind] && border_right >= CUR_NBINS[ind] )
	{
		border_right -= CUR_NBINS[ind];
	}

	CUR_PTS_IN_BIN[ind][border_left] += 1;
//	printf("Adding point %lf mag %lf bin %d\n", x, f, border_left );

	double val_A = A(x,x0,x1);
	double val_B = B(x,x0,x1);
	double val_C = C(x,x0,x1);
	double val_D = D(x,x0,x1);
 
	// d f_s / d (y_m1) == A
	// d f_s / d (y_p1) == B
	// d f_s / d (ydd_m1) == C
	// d f_s / d (ydd_p1) == D

	int p_y_m1 = 2*border_left+0;
	int p_y_p1 = 2*border_right+0;
	int p_ydd_m1 = 2*border_left+1;
	int p_ydd_p1 = 2*border_right+1;

	ADD_B[p_y_m1] += f * val_A;
	ADD_B[p_y_p1] += f * val_B;
	ADD_B[p_ydd_m1] += f * val_C;
	ADD_B[p_ydd_p1] += f * val_D;

	ADD_A[p_y_m1*nparamstot[ind]+p_y_m1] += val_A * val_A;
	ADD_A[p_y_m1*nparamstot[ind]+p_y_p1] += val_B * val_A;
	ADD_A[p_y_m1*nparamstot[ind]+p_ydd_m1] += val_C * val_A;
	ADD_A[p_y_m1*nparamstot[ind]+p_ydd_p1] += val_D * val_A;
	
	ADD_A[p_y_p1*nparamstot[ind]+p_y_m1] += val_A * val_B;
	ADD_A[p_y_p1*nparamstot[ind]+p_y_p1] += val_B * val_B;
	ADD_A[p_y_p1*nparamstot[ind]+p_ydd_m1] += val_C * val_B;
	ADD_A[p_y_p1*nparamstot[ind]+p_ydd_p1] += val_D * val_B;
	
	ADD_A[p_ydd_m1*nparamstot[ind]+p_y_m1] += val_A * val_C;
	ADD_A[p_ydd_m1*nparamstot[ind]+p_y_p1] += val_B * val_C;
	ADD_A[p_ydd_m1*nparamstot[ind]+p_ydd_m1] += val_C * val_C;
	ADD_A[p_ydd_m1*nparamstot[ind]+p_ydd_p1] += val_D * val_C;
	
	ADD_A[p_ydd_p1*nparamstot[ind]+p_y_m1] += val_A * val_D;
	ADD_A[p_ydd_p1*nparamstot[ind]+p_y_p1] += val_B * val_D;
	ADD_A[p_ydd_p1*nparamstot[ind]+p_ydd_m1] += val_C * val_D;
	ADD_A[p_ydd_p1*nparamstot[ind]+p_ydd_p1] += val_D * val_D;
}

void SolveSpline( int ind)
{
	double *tA = (double *)malloc( sizeof(double) * nparamstot[ind] * nparamstot[ind] );
	memcpy( tA, CUR_A[ind], sizeof(double) * nparamstot[ind] * nparamstot[ind] );
	memcpy( CUR_SOLUTION[ind], CUR_B[ind], sizeof(double) * nparamstot[ind] );

	int nrhs = 1;
	int ipiv[nparamstot[ind]];
	int info = 0;

	// look for all-zero rows of A.

	double dx = (CUR_XMAX[ind]-CUR_XMIN[ind])/CUR_NBINS[ind];
	/*for( int b = 0; b < CUR_NBINS; b++ )
	{
		double bin_cen = CUR_XMIN + dx * (b+0.5);
		if( CUR_PTS_IN_BIN[b] == 0 )
			AddPointToSpecificMatrix( bin_cen, 0.0, tA, CUR_SOLUTION );
	}*/
	
	for( int t= 0; t < nparams[ind]; t++ )
		tA[t*nparamstot[ind]+t] += 1e-10;

	//dgesv( &nparams, &nrhs, tA, &nparamstot, ipiv, CUR_SOLUTION, &nparamstot, &info );
	dgesv( &nparamstot[ind], &nrhs, tA, &nparamstot[ind], ipiv, CUR_SOLUTION[ind], &nparamstot[ind], &info );

	if( info != 0 )
	{
		printf("Failed to solve for spline solution.\n");
		exit(1);
	}
	free(tA);
	
}

double evaluateSpline( double x, int ind)
{
	if( is_periodic[ind] )
	{
		if( is_periodic[ind] )
		{
			while( x < CUR_XMIN[ind] ) x += (CUR_XMAX[ind]-CUR_XMIN[ind]);
			while( x >= CUR_XMAX[ind] ) x -= (CUR_XMAX[ind]-CUR_XMIN[ind]);
		}
	}
	else
	{
		if( !( x >= CUR_XMIN[ind] && x <= CUR_XMAX[ind]) ) 
			return 0;
	}

	int border_left  = CUR_NBINS[ind] * (x - CUR_XMIN[ind]) / (CUR_XMAX[ind]-CUR_XMIN[ind]);
	int border_right = border_left + 1;
	

	int bin = border_left;

	if( CUR_PTS_IN_BIN[ind][bin] == 0 ) return 0;

	double dx = (CUR_XMAX[ind]-CUR_XMIN[ind])/CUR_NBINS[ind];

	double x0 = CUR_XMIN[ind] + border_left * dx;
	double x1 = CUR_XMIN[ind] + border_right * dx;
	
	if( is_periodic[ind] && border_right >= CUR_NBINS[ind] )
		border_right -= CUR_NBINS[ind];
	

	double val_A = A(x,x0,x1);
	double val_B = B(x,x0,x1);
	double val_C = C(x,x0,x1);
	double val_D = D(x,x0,x1);
	
	int p_y_m1   = 2*border_left+0;
	int p_y_p1   = 2*border_right+0;
	int p_ydd_m1 = 2*border_left+1;
	int p_ydd_p1 = 2*border_right+1;

	return val_A * CUR_SOLUTION[ind][p_y_m1] +
	       val_B * CUR_SOLUTION[ind][p_y_p1] +
	       val_C * CUR_SOLUTION[ind][p_ydd_m1] +
	       val_D * CUR_SOLUTION[ind][p_ydd_p1];
}


double splinedFdX( double x, int ind )
{
	if( is_periodic[ind] )
	{
		while( x < CUR_XMIN[ind] ) x += (CUR_XMAX[ind]-CUR_XMIN[ind]);
		while( x >= CUR_XMAX[ind] ) x -= (CUR_XMAX[ind]-CUR_XMIN[ind]);
	}
	else
	{
		if( !( x >= CUR_XMIN[ind] && x <= CUR_XMAX[ind]) ) 
			return 0;
	}

	int border_left  = CUR_NBINS[ind] * (x - CUR_XMIN[ind]) / (CUR_XMAX[ind]-CUR_XMIN[ind]);
	int border_right = border_left + 1;
	
	

	int bin = border_left;

	if( CUR_PTS_IN_BIN[ind][bin] == 0 ) return 0;

	double dx = (CUR_XMAX[ind]-CUR_XMIN[ind])/CUR_NBINS[ind];

	double x0 = CUR_XMIN[ind] + border_left * dx;
	double x1 = CUR_XMIN[ind] + border_right * dx;
	
	if( is_periodic[ind] && border_right >= CUR_NBINS[ind] )
		border_right -= CUR_NBINS[ind];

	double val_dA = dA(x,x0,x1);
	double val_dB = dB(x,x0,x1);
	double val_dC = dC(x,x0,x1);
	double val_dD = dD(x,x0,x1);
	
	int p_y_m1   = 2*border_left+0;
	int p_y_p1   = 2*border_right+0;
	int p_ydd_m1 = 2*border_left+1;
	int p_ydd_p1 = 2*border_right+1;

	return val_dA * CUR_SOLUTION[ind][p_y_m1] +
	       val_dB * CUR_SOLUTION[ind][p_y_p1] +
	       val_dC * CUR_SOLUTION[ind][p_ydd_m1] +
	       val_dD * CUR_SOLUTION[ind][p_ydd_p1];
}

double splined2FdX2( double x, int ind )
{
	if( is_periodic[ind] )
	{
		if( is_periodic[ind] )
		{
			while( x < CUR_XMIN[ind] ) x += (CUR_XMAX[ind]-CUR_XMIN[ind]);
			while( x >= CUR_XMAX[ind] ) x -= (CUR_XMAX[ind]-CUR_XMIN[ind]);
		}
	}
	else
	{
		if( !( x >= CUR_XMIN[ind] && x <= CUR_XMAX[ind]) ) 
			return 0;
	}
	int border_left  = CUR_NBINS[ind] * (x - CUR_XMIN[ind]) / (CUR_XMAX[ind]-CUR_XMIN[ind]);
	int border_right = border_left + 1;
	int bin = border_left;

	if( CUR_PTS_IN_BIN[ind][bin] == 0 ) return 0;

	double dx = (CUR_XMAX[ind]-CUR_XMIN[ind])/CUR_NBINS[ind];

	double x0 = CUR_XMIN[ind] + border_left * dx;
	double x1 = CUR_XMIN[ind] + border_right * dx;
	
	if( is_periodic[ind] && border_right >= CUR_NBINS[ind] )
		border_right -= CUR_NBINS[ind];

	double val_ddA = ddA(x,x0,x1);
	double val_ddB = ddB(x,x0,x1);
	double val_ddC = ddC(x,x0,x1);
	double val_ddD = ddD(x,x0,x1);
	
	int p_y_m1   = 2*border_left+0;
	int p_y_p1   = 2*border_right+0;
	int p_ydd_m1 = 2*border_left+1;
	int p_ydd_p1 = 2*border_right+1;

	return val_ddA * CUR_SOLUTION[ind][p_y_m1] +
	       val_ddB * CUR_SOLUTION[ind][p_y_p1] +
	       val_ddC * CUR_SOLUTION[ind][p_ydd_m1] +
	       val_ddD * CUR_SOLUTION[ind][p_ydd_p1];
}

