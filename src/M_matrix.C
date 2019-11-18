#define __mmatrixh__
#include "interp.h"
#include "M_matrix.h"
#include <string.h>
static double **M = NULL;
static int mlow = 5;
static int mhigh = 7;
void generateSubdivisionMatrices( double **M, int mlow, int mhigh );

void getM( double ***M_out, int *mlow_out, int *mhigh_out )
{
	if( !M )
	{
		M = (double **)malloc( sizeof(double*) * (mhigh-mlow+1) );

		for( int mx = 0; mx <= mhigh-mlow; mx++ )
		{
			int val1 = 12;
			int val2 = 6 + mlow+mx;
			if( val2 > val1 )
				val1 = val2;
 
			M[mx] = (double *)malloc( sizeof(double) * 4 * val1 * val2 );	
		}

		generateSubdivisionMatrices( M, mlow, mhigh );

	}

	*M_out = M;
	*mlow_out = mlow;
	*mhigh_out = mhigh;
}

#define MAX_INV_VALENCE 12

void generateSubdivisionMatrices( double **M, int mlow, int mhigh )
{
	int xm = 0;
	for( int val = mlow; val <= mhigh; val++, xm++ )
	{
		double *put = M[xm];

		double w = 3.0/(8*val);
		double w6 = 3.0/(8*6);
		int ncoords_base  = 1 + val + 5; 

		// assemble the matrix:

		int ncoords_extra = 6;
	
		double *A = (double *)malloc( sizeof(double) * (ncoords_base+ncoords_extra) * ncoords_base ); 
		memset( A, 0, sizeof(double) * (ncoords_base+ncoords_extra)*ncoords_base );

		// point 0 (i), a carried vertex. 
		A[0*ncoords_base+0] = 1-val*w;
		for( int p = 1; p < 1 + val; p++ )
			A[0*ncoords_base+p] = w;
		
		A[1*ncoords_base+0]   = (3.0/8.0);	
		A[1*ncoords_base+1]   = (3.0/8.0);	
		A[1*ncoords_base+2]   = (1.0/8.0);	
		A[1*ncoords_base+val] = (1.0/8.0);	

		for( int p = 2; p < val; p++ )
		{
			A[p*ncoords_base+0]   = (3.0/8.0);
			A[p*ncoords_base+p]   = (3.0/8.0);
			A[p*ncoords_base+p-1] = (1.0/8.0);
			A[p*ncoords_base+p+1] = (1.0/8.0);
		}
		
		A[val*ncoords_base+0]     = (3.0/8.0);
		A[val*ncoords_base+val]   = (3.0/8.0);
		A[val*ncoords_base+val-1] = (1.0/8.0);
		A[val*ncoords_base+1]     = (1.0/8.0);

		A[(val+1)*ncoords_base+1] = (3.0/8.0);	
		A[(val+1)*ncoords_base+val] = (3.0/8.0);	
		A[(val+1)*ncoords_base+0] = (1.0/8.0);	
		A[(val+1)*ncoords_base+val+1] = (1.0/8.0);	

		// point 1 (k)			
		A[(val+2)*ncoords_base+1] = 1-6*w6;
			A[(val+2)*ncoords_base+(0)] = w6;
			A[(val+2)*ncoords_base+(val)] = w6;
			A[(val+2)*ncoords_base+(val+1)] = w6;
			A[(val+2)*ncoords_base+(val+2)] = w6;
			A[(val+2)*ncoords_base+(val+3)] = w6;
			A[(val+2)*ncoords_base+(2)] = w6;
		
		A[(val+3)*ncoords_base+1] = (3.0/8.0);	
		A[(val+3)*ncoords_base+2] = (3.0/8.0);	
		A[(val+3)*ncoords_base+0] = (1.0/8.0);	
		A[(val+3)*ncoords_base+(val+3)] = (1.0/8.0);	
		
		// point 2 (j)
		A[(val+4)*ncoords_base+2] = 1-6*w6;
			A[(val+4)*ncoords_base+(0)] = w6;
			A[(val+4)*ncoords_base+(1)] = w6;
			A[(val+4)*ncoords_base+(val+3)] = w6;
			A[(val+4)*ncoords_base+(val+4)] = w6;
			A[(val+4)*ncoords_base+(val+5)] = w6;
			A[(val+4)*ncoords_base+(3)] = w6;
		
		A[(val+5)*ncoords_base+2]       = (3.0/8.0);	
		A[(val+5)*ncoords_base+3]       = (3.0/8.0);	
		A[(val+5)*ncoords_base+0]       = (1.0/8.0);	
		A[(val+5)*ncoords_base+(val+5)] = (1.0/8.0);	


		// now, the extra coordinates.

		A[(val+6)*ncoords_base+1] = (3.0/8.0);
		A[(val+6)*ncoords_base+(val+1)] = (3.0/8.0); 
		A[(val+6)*ncoords_base+(val+2)] = (1.0/8.0); 
		A[(val+6)*ncoords_base+val] = (1.0/8.0); 
		
		A[(val+7)*ncoords_base+1] = (3.0/8.0);
		A[(val+7)*ncoords_base+(val+2)] = (3.0/8.0); 
		A[(val+7)*ncoords_base+(val+3)] = (1.0/8.0); 
		A[(val+7)*ncoords_base+(val+1)] = (1.0/8.0); 
		
		A[(val+8)*ncoords_base+1] = (3.0/8.0);
		A[(val+8)*ncoords_base+(val+3)] = (3.0/8.0); 
		A[(val+8)*ncoords_base+2] = (1.0/8.0); 
		A[(val+8)*ncoords_base+(val+2)] = (1.0/8.0); 
		
		A[(val+9)*ncoords_base+2] = (3.0/8.0);
		A[(val+9)*ncoords_base+(val+3)] = (3.0/8.0); 
		A[(val+9)*ncoords_base+1] = (1.0/8.0); 
		A[(val+9)*ncoords_base+(val+4)] = (1.0/8.0); 
		
		A[(val+10)*ncoords_base+2] = (3.0/8.0);
		A[(val+10)*ncoords_base+(val+4)] = (3.0/8.0); 
		A[(val+10)*ncoords_base+(val+5)] = (1.0/8.0); 
		A[(val+10)*ncoords_base+(val+3)] = (1.0/8.0); 
		
		A[(val+11)*ncoords_base+2] = (3.0/8.0);
		A[(val+11)*ncoords_base+(val+5)] = (3.0/8.0); 
		A[(val+11)*ncoords_base+3] = (1.0/8.0); 
		A[(val+11)*ncoords_base+(val+4)] = (1.0/8.0); 
			
		int nvmax = 12;
		if( ncoords_base > nvmax ) nvmax = ncoords_base;
		int pcycle[nvmax];

		for( int x = 0; x < nvmax; x++ )
			pcycle[x] = x;
	
		int max_result = 12;

		if( ncoords_base > max_result )
			max_result = ncoords_base;

		// the four subdivision triangles.
		for( int subM = 0; subM < 4; subM++ )
		{
			int nc = ncoords_base;

			if( subM == 1)
			{
				nc = 12;
				pcycle[0] = 2;
				pcycle[1] = val+3;
				pcycle[2] = val+4;
				pcycle[3] = val+5;
				pcycle[4] = 3;
				pcycle[5] = 0;
				pcycle[6] = 1;
				pcycle[7] = val+2;
				pcycle[8] = val+8;
				pcycle[9] = val+9;
				pcycle[10] =val+10;
				pcycle[11] = val+11;
			}
			else if( subM == 2)
			{
				nc = 12;
				pcycle[0] = 1;
				pcycle[1] = val+2;
				pcycle[2] = val+3;
				pcycle[3] = 2;
				pcycle[4] = 0;
				pcycle[5] = val;
				pcycle[6] = val+1;
				pcycle[7] = val+6;
				pcycle[8] = val+7;
				pcycle[9] = val+8;
				pcycle[10] =  val+9;
				pcycle[11] =  val+4;
			}
			else if( subM == 3 )
			{
				nc = 12;
				pcycle[0] = val+3; 
				pcycle[1] = 2;
				pcycle[2] = 1;
				pcycle[3] = val+2;
				pcycle[4] = val+8;
				pcycle[5] = val+9;
				pcycle[6] = val+4;
				pcycle[7] = val+5;
				pcycle[8] = 3;
				pcycle[9] = 0;
				pcycle[10] = val;
				pcycle[11] = val+1;
			}

			// get the subdivision matrix for this triangle.	
			for( int x = 0; x < nc; x++ )
			for( int y = 0; y < ncoords_base; y++ )
				put[subM*ncoords_base*max_result+x*ncoords_base+y] = A[pcycle[x]*ncoords_base+y];	
		}

		free(A);
	}
}
