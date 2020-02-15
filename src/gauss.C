#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define INT_TOL (1e-9)
#define SOME_NUMBER 100
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "gsl/gsl_integration.h"
#include "lapack_we_use.h"

static gsl_integration_workspace *gsl_i_workspace = NULL; 

struct fblock
{
	double k;
	double v;
	int n;
	int type;
	double *a;
	double *b;
};

double eval_fn( double u, void *params )
{
	fblock *theBlock = (fblock *)params;

	double k = theBlock->k;
	double v = theBlock->v;

	double val = (1.0/(u+v)) * pow( k, -log2(u+v));
	int n = theBlock->n;

	double *a = theBlock->a;
	double *b = theBlock->b;

	double pnm1 = 0;
	double pn   = 1;

	for( int i = 0; i < theBlock->n; i++ )
	{
		double p = (u - a[i]) * pn - b[i] * pnm1;

		pnm1 = pn;
		pn = p;
	} 

	if( theBlock->type == 1 )	
	{
//		printf("type: %d n: %d u: %lf val: %le val: %le k: %le\n", 1, theBlock->n, u,u*pn * pn * val * val, val, k );
		return u * pn * pn /( val * val * val * val);
	}
	else
	{
//		printf("type: %d n: %d u: %lf val: %le val: %le\n", 1, theBlock->n, u, pn * pn * val * val, val );
		return pn * pn  /( val * val * val * val);
	}
}

double eval_fn_v( double v, void *params )
{
	fblock *theBlock = (fblock *)params;

	double k = theBlock->k;

	double val = - (-1+pow(k, 4.0 * log2(v))*v*v*v*v*v) * log(2.0) / (log(32.0)+4.0*log(k));
	int n = theBlock->n;

	double *a = theBlock->a;
	double *b = theBlock->b;

	double pnm1 = 0;
	double pn   = 1;

	for( int i = 0; i < theBlock->n; i++ )
	{
		double p = (v - a[i]) * pn - b[i] * pnm1;

		pnm1 = pn;
		pn = p;
	} 

	if( theBlock->type == 1 )	
	{
		return v * pn * pn * val;
	}
	else
	{
		return pn * pn * val;
	}
}

void get_u_pts( double *u_pts, double *w_pts, int npts, double eval, double v )
{
	if( gsl_i_workspace == NULL )
		gsl_i_workspace = gsl_integration_workspace_alloc( SOME_NUMBER );

	fblock theBlock;

	theBlock.a = (double *)malloc( sizeof(double) * (npts+1) );
	theBlock.b = (double *)malloc( sizeof(double) * (npts+1) );
	
	memset( theBlock.a, 0, sizeof(double) * (npts+1) );
	memset( theBlock.b, 0, sizeof(double) * (npts+1) );

	gsl_function thef;
	thef.params = &theBlock;

	thef.function = eval_fn;

	theBlock.k = eval;
	theBlock.n = 0;
	theBlock.v = v;

	double prev_a_den = 1;

	double mu0 = 0;

	for( int i = 0; i <= npts; i++ )
	{
		double f, aerr;

		theBlock.n = i;

		double a_num = 0, a_den = 0, b_num = 0, b_den = 0;
		
		theBlock.type = 1;
		// could use qagp since the singularity is at zero.
		gsl_integration_qags( &thef, 0, 1-v, INT_TOL, 0., SOME_NUMBER, gsl_i_workspace, &a_num, &aerr );

		theBlock.type = 0;
		gsl_integration_qags( &thef, 0, 1-v, INT_TOL, 0., SOME_NUMBER, gsl_i_workspace, &a_den, &aerr );

		theBlock.a[i] = a_num / a_den;

		if( i > 0 )
			theBlock.b[i] = a_den / prev_a_den;
		else
			theBlock.b[i] = 0;

		if( i == 0 )
			mu0 = a_den;

		prev_a_den = a_den;		

//		printf("i: %d num: %le den: %le\n", i, a_num, a_den );
//
//		printf("a[%d]: %le b[%d]: %le\n", i,theBlock.a[i],i, theBlock.b[i] ); 
	}
	
	double *GMat = (double *)malloc( sizeof(double) * npts*npts);

	memset( GMat, 0, sizeof(double) * npts * npts );
	int n = npts;
	for( int i = 0; i < n; i++ )
	{
		GMat[i*n+i] = theBlock.a[i];

		if( i > 0 )
		{
			GMat[(i-1)*n+i] = sqrt(theBlock.b[i]);
			GMat[i*n+i-1]   = sqrt(theBlock.b[i]);
		}
	}
	
	char jobz = 'V';
	char uplo = 'U';
	int N = npts;
	int LDA = N;
	double WR[npts];
	int LDVA = npts;
	int LWORK = 8*npts + npts*npts;
	double *work = (double *)malloc( sizeof(double) * LWORK );
	int info;
	dsyev( &jobz, &uplo, &N, GMat, &LDA, WR, work, &LWORK, &info ); 

	for( int i = 0; i < npts; i++ )
	{
		double nrm = 0;
	
		for( int j = 0; j< npts; j++ )
			nrm += GMat[j+npts*i] * GMat[j+npts*i];

		u_pts[i] = WR[i];
		w_pts[i] = mu0 * GMat[i*npts] * GMat[i*npts] / nrm;	
//		printf("i: %d nrm: %lf u: %le w: %le\n", i, nrm, u_pts[i], w_pts[i] );
		
	}

	free(work);
	free(GMat);
	free(theBlock.a);
	free(theBlock.b);	
}

void get_v_pts( double *v_pts, double *w_pts, int npts, double eval )
{
	if( gsl_i_workspace == NULL )
		gsl_i_workspace = gsl_integration_workspace_alloc( SOME_NUMBER );

	fblock theBlock;

	theBlock.a = (double *)malloc( sizeof(double) * (npts+1) );
	theBlock.b = (double *)malloc( sizeof(double) * (npts+1) );
	
	memset( theBlock.a, 0, sizeof(double) * (npts+1) );
	memset( theBlock.b, 0, sizeof(double) * (npts+1) );

	gsl_function thef;
	thef.params = &theBlock;

	thef.function = eval_fn_v;

	theBlock.k = eval;
	theBlock.n = 0;

	double prev_a_den = 1;

	double mu0 = 0;

	for( int i = 0; i <= npts; i++ )
	{
		double f, aerr;

		theBlock.n = i;

		double a_num = 0, a_den = 0, b_num = 0, b_den = 0;
		
		theBlock.type = 1;
		// could use qagp since the singularity is at zero.
		gsl_integration_qags( &thef, 0, 1, INT_TOL, 0., SOME_NUMBER, gsl_i_workspace, &a_num, &aerr );

		theBlock.type = 0;
		gsl_integration_qags( &thef, 0, 1, INT_TOL, 0., SOME_NUMBER, gsl_i_workspace, &a_den, &aerr );

		theBlock.a[i] = a_num / a_den;

		if( i > 0 )
			theBlock.b[i] = a_den / prev_a_den;
		else
			theBlock.b[i] = 0;

		if( i == 0 )
			mu0 = a_den;

		prev_a_den = a_den;		

//		printf("i: %d num: %le den: %le\n", i, a_num, a_den );
//
//		printf("a[%d]: %le b[%d]: %le\n", i,theBlock.a[i],i, theBlock.b[i] ); 
	}
	
	double *GMat = (double *)malloc( sizeof(double) * npts*npts);

	memset( GMat, 0, sizeof(double) * npts * npts );
	int n = npts;
	for( int i = 0; i < n; i++ )
	{
		GMat[i*n+i] = theBlock.a[i];

		if( i > 0 )
		{
			GMat[(i-1)*n+i] = sqrt(theBlock.b[i]);
			GMat[i*n+i-1]   = sqrt(theBlock.b[i]);
		}
	}
	
	char jobz = 'V';
	char uplo = 'U';
	int N = npts;
	int LDA = N;
	double WR[npts];
	int LDVA = npts;
	int LWORK = 8*npts + npts*npts;
	double *work = (double *)malloc( sizeof(double) * LWORK );
	int info;
	dsyev( &jobz, &uplo, &N, GMat, &LDA, WR, work, &LWORK, &info ); 

	for( int i = 0; i < npts; i++ )
	{
		double nrm = 0;
	
		for( int j = 0; j< npts; j++ )
			nrm += GMat[j+npts*i] * GMat[j+npts*i];

		v_pts[i] = WR[i];
		w_pts[i] = mu0 * GMat[i*npts] * GMat[i*npts] / nrm;	
//		printf("i: %d nrm: %lf u: %le w: %le\n", i, nrm, u_pts[i], w_pts[i] );
		
	}

	free(work);
	free(GMat);
	free(theBlock.a);
	free(theBlock.b);	
}
