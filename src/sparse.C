#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "util.h"
#include <sys/time.h>
#include "sparse.h"
#include "lapack_we_use.h"
double sp_fac( int n )
{
	double x = 1;

	for( int t = 1; t <= n; t++ )
		x *= t;
	return x;
}

void SparseMatrix::mvec( double *x, double *b )
{

	for( int i = 0; i < n; i++ )
	{
		double u = 0;
		
		for( int xj = 0; xj < nnz[i]; xj++ )
			u += nzv[i][xj] * x[nzl[i][xj]];

		b[i] = u;
	}
}

double SparseMatrix::melem( double *b )
{
	double t = 0;
	double u = 0;

	for( int i = 0; i < n; i++ )
	{
		u = 0;
		
		for( int xj = 0; xj < nnz[i]; xj++ )
			u += nzv[i][xj] * b[nzl[i][xj]];

		t += u * b[i];
	}
	
	return t;		
}

void SparseMatrix::debug( double *A )
{
	for( int i = 0; i < n; i++ )
	{
		int k = 0;
		
		for( int xj = 0; xj < nnz[i]; xj++ )
		{
			for( ; k < nzl[i][xj]; k++ )
			{
				if( fabs(A[i*n+k]) > 1e-10 )
				{
					printf("Non zero value not present in sparse matrix!!\n");
					exit(1);
				}
			}

			double del = fabs( A[i*n+nzl[i][xj]] - nzv[i][xj] );
			if( del > 1e-10 )
			{
				printf("Sparse matrix and A differ by %le!\n", del );
				exit(1);
			}

			k = nzl[i][xj] + 1;
		}
	}
}

void SparseMatrix::cgsolve( double *b, double *soln )
{
	double *b_use = (double *)malloc( sizeof(double) * n );
	memcpy( b_use, b, sizeof(double) * n );
	double *preconditioner = (double *)malloc( sizeof(double) * n );
	double b_norm = 1;

	
	double *r0 = (double *)malloc( sizeof(double) * n );
	double *r1 = (double *)malloc( sizeof(double) * n );
	double *p0 = (double *)malloc( sizeof(double) * n );
	double *p1 = (double *)malloc( sizeof(double) * n );
	double *x0 = (double *)malloc( sizeof(double) * n );
	double *x1 = (double *)malloc( sizeof(double) * n );
	double *z0 = (double *)malloc( sizeof(double) * n );
	double *z1 = (double *)malloc( sizeof(double) * n );
	
	for( int x = 0; x < n; x++ )
	{
		if( diage[x] == -1 )
		{
			printf("Trouble, matrix has zero on the diagonal at %d.\n", x );
			exit(1);
		}


		x0[x] = b_use[x];// / nzv[x][diage[x]];
//		x0[x] = (double)rand() / (double)RAND_MAX;

		preconditioner[x] = nzv[x][diage[x]];
	}


	int done = 0;

	mvec( x0, r0 );

	for( int x = 0; x < n; x++ )
	{
		r0[x] = b_use[x] - r0[x];
		z0[x] = (1.0 / (preconditioner[x]+1e-30)) * r0[x];
	}

	memcpy( p0, r0, sizeof(double) * n );

	int iter = 0;
	while( !done )
	{
		double zr = 0;
		for( int x = 0; x < n; x++ )
			zr += z0[x]*r0[x]; 

		double inner = melem( p0 );
		if( fabs(inner) < 1e-30 ) inner = 1e-30;
		double alpha = zr / inner;
		
		double del = 0;		
		for( int x = 0; x < n; x++ )
		{
			del += alpha * p0[x] * alpha * p0[x];
			x1[x] = x0[x] + alpha * p0[x];
		}	
	// temporary space	
		mvec( p0, p1 );	

		for( int x = 0; x < n; x++ )
			r1[x] = r0[x] - alpha * p1[x];	
		for( int x = 0; x < n; x++ )
			z1[x] = (1.0 / preconditioner[x]) * r1[x];

		double z1r1 = 0;
		for( int x = 0; x < n; x++ )
			z1r1 += z1[x] * r1[x];
		double beta = z1r1 / zr;

		for( int x = 0; x < n; x++ )
			p1[x] = z1[x] + beta * p0[x];
		
		memcpy( p0, p1, sizeof(double) * n );
		memcpy( x0, x1, sizeof(double) * n );
		memcpy( r0, r1, sizeof(double) * n );
		memcpy( z0, z1, sizeof(double) * n );
		
		if( iter % 50 == 0 )
		printf("iter %d del %.12le\n", iter, sqrt(del) );

		if( sqrt(del) < 1e-15 || iter > 2000 )
			break;

		iter++;
	}

	memcpy( soln, x1, sizeof(double) * n );

	mvec( soln, x1 );

//	for( int x = 0; x < n; x++ )
//		printf("%d %.12le %.12le %.12le\n", x, b[x], x1[x], b[x]-x1[x] );

	free(x0);
	free(x1);
	free(r0);
	free(r1);
	free(p0);
	free(p1);
	free(z0);
	free(z1);
	free(b_use);

	for( int x = 0; x < n; x++ )
		soln[x] *= b_norm; 
}

void SparseMatrix::solve2( double *b, double *soln, double omega )
{
	double *b_use = (double *)malloc( sizeof(double) * n );
	memcpy( b_use, b, sizeof(double) * n );

	double b_norm = 1;
	double *op = (double *)malloc( sizeof(double) * n );

	for( int x = 0; x < n; x++ )
	{
		if( diage[x] == -1 )
		{
			printf("Trouble, matrix has zero on the diagonal at %d.\n", x );
			exit(1);
		}

		op[x] = b_use[x] / nzv[x][diage[x]];
	}
	
	memset( op, 0, sizeof(double) * n );

	int done = 0;

	double *curs = (double *)malloc( sizeof(double) * n );

	int iter = 1;

	while( !done )
	{
		memset( curs, 0, sizeof(double) * n );

		mvec( op, curs );

		double del = 0;
		for( int i = 0; i < n; i++ )
		{
			del += (b_use[i] - curs[i]) * ( b_use[i] - curs[i] );
			op[i] += omega * b_use[i] - omega * curs[i]; 
		}

		iter++;

		printf("Iter %d del: %le\n", iter, del );

		if( del < 1e-10)
			break;
	}

	for( int i = 0; i < n; i++ )
	{	
		double tv = 0;

		for( int xj = 0; xj < nnz[i]; xj++ )
			tv += nzv[i][xj] * op[nzl[i][xj]];			
//		printf("%lf %lf\n", b[i], b_norm * tv );
	}

	free(b_use);
	free(curs);

	memcpy( soln, op, sizeof(double) * n );

	for( int x = 0; x < n; x++ )
		soln[x] *= b_norm; 
/*
	printf("soln:\n");
	for( int x = 0; x < n; x++  )
		printf("%lf ", soln[x] );
	printf("\n");*/
	free(op);
}
void SparseMatrix::solve( double *b, double *soln, double omega )
{
	double *b_use = (double *)malloc( sizeof(double) * n );
	memcpy( b_use, b, sizeof(double) * n );

	double b_norm = 1;
/*	for( int x = 0; x < n; x++ )
		b_norm += b_use[x] * b_use[x];
	for( int x = 0; x < n; x++ )
		b_use[x] /= b_norm;
*/
	double *op = (double *)malloc( sizeof(double) * n );

	for( int x = 0; x < n; x++ )
	{
		if( diage[x] == -1 )
		{
			printf("Trouble, matrix has zero on the diagonal at %d.\n", x );
			exit(1);
		}

		op[x] = b_use[x] / nzv[x][diage[x]];
	}
	
//	memset( op, 0, sizeof(double) * n );

	int done = 0;

	double *curs = (double *)malloc( sizeof(double) * n );

	int iter = 1;

	while( !done )
	{
		memset( curs, 0, sizeof(double) * n );

		for( int i = 0; i < n; i++ )
		{
			double sigma = 0;

			for( int xj = 0; xj < nnz[i]; xj++ )
			{
				if( nzl[i][xj] >= i )
					break;
				
				sigma += nzv[i][xj] * curs[nzl[i][xj]];
			}
			
			for( int xj = 0; xj < nnz[i]; xj++ )
			{
				if( nzl[i][xj] > i )
				{
				
					sigma += nzv[i][xj] * op[nzl[i][xj]];
				}
			}
			
			curs[i] = (1-omega) * op[i] + (omega / nzv[i][diage[i]]) * (b_use[i] - sigma);
		}

		double del = 0;
		
		for( int j = 0; j < n; j++ )
			del += (curs[j] - op[j]) * (curs[j] - op[j]);


		memcpy( op, curs, sizeof(double) * n );

		iter++;
		printf("Iter %d del: %le\n", iter, del );
		if( del/(omega*omega) < 1e-10)
			break;
	}

	for( int i = 0; i < n; i++ )
	{	
		double tv = 0;

		for( int xj = 0; xj < nnz[i]; xj++ )
			tv += nzv[i][xj] * op[nzl[i][xj]];			
//		printf("%lf %lf\n", b[i], b_norm * tv );
	}

	free(b_use);
	free(curs);

	memcpy( soln, op, sizeof(double) * n );

	for( int x = 0; x < n; x++ )
		soln[x] *= b_norm; 
/*
	printf("soln:\n");
	for( int x = 0; x < n; x++  )
		printf("%lf ", soln[x] );
	printf("\n");*/
	free(op);
}

void SparseMatrix::sort( void )
{
	for( int x = 0; x < n; x++ )
	{
		int done = 0;

		while( !done )
		{
			done = 1;

			for( int p = 0; p < nnz[x]-1; p++ )
			{
				if( nzl[x][p] > nzl[x][p+1] )
				{
					int it = nzl[x][p];
					double iv = nzv[x][p];

					nzl[x][p] = nzl[x][p+1];
					nzv[x][p] = nzv[x][p+1];
					nzl[x][p+1] = it;
					nzv[x][p+1] = iv;
					done = 0;
				}
			}
		}

		for( int p = 0; p < nnz[x]; p++ )
		{
			if( nzl[x][p] == x )
				diage[x] = p;
		}
	}	
}

void SparseMatrix::setNeedSource( void )
{
	for( int a = 0; a < n; a++ )
	{
		if( nnz[a] > 0 )
			need(a);
	
		for( int bx = 0; bx < nnz[a]; bx++ )
			source( nzl[a][bx] );
	}
	
}

void SparseMatrix::coupleParameters( int p1, int p2, double val, double tol )
{
	int got_it = 0;

	for( int x = 0; x < nnz[p1]; x++ )
	{
		if( nzl[p1][x] == p2 )
		{
			if( p1 == p2 )
				diagonal_element[p1] += val;

			nzv[p1][x] += val;
			got_it = 1;
			break;
		}
	}	

	if( !got_it && fabs(val) > tol )
	{
		if( nzl_space[p1] == nnz[p1] )
		{
			nzl_space[p1] *= 2;
			nzl[p1] = ( int *)realloc( nzl[p1], sizeof(int) * nzl_space[p1] );
			nzv[p1] = ( double *)realloc( nzv[p1], sizeof(double) * nzl_space[p1] );
		}

		if( p2 == p1 )
		{
			diage[p1] = nnz[p1];
			diagonal_element[p1] = val;
		}
		nzl[p1][nnz[p1]] = p2;
		nzv[p1][nnz[p1]] = val;

		nnz[p1] += 1;
	}

#if 0
	got_it = 0;

	for( int x = 0; x < nnz[p2]; x++ )
	{
		if( nzl[p2][x] == p1 )
		{
			nzv[p2][x] += val;
			got_it = 1;
			break;
		}
	}	
	
	if( !got_it )
	{
		if( nzl_space[p2] == nnz[p2] )
		{
			nzl_space[p2] *= 2;
			nzl[p2] = ( int *)realloc( nzl[p2], sizeof(int) * nzl_space[p2] );
			nzv[p2] = ( double *)realloc( nzv[p2], sizeof(double) * nzl_space[p2] );
		}
		
		if( p2 == p1 )
			diage[p1] = nnz[p1];

		nzl[p2][nnz[p2]] = p1;
		nzv[p2][nnz[p2]] = val;

		nnz[p2] += 1;
	}
#endif
}

void SparseMatrix::loadUnit( void )
{
	for( int t = 0; t < n; t++ )
	{
		nnz[t] = 1;
		nzl[t][0] = t;
		nzv[t][0] = 1;

		diagonal_element[t] = 1;
	}
}

void SparseMatrix::scale (double fac)
{
	for( int t = 0; t < n; t++ )
	{
		diagonal_element[t] *= fac;
		
		for( int x = 0; x < nnz[t]; x++ )
			nzv[t][x] *= fac;
	}
}

void SparseMatrix::clear( void )
{
	for( int t = 0; t < n; t++ )
	{
		nnz[t] = 0;
		diagonal_element[t] = 0;
	}
}

void SparseMatrix::init( int new_n )
{
	n = new_n;

	need_list = (int *)malloc( sizeof(int) * new_n );
	source_list = (int *)malloc( sizeof(int) * new_n );

	nneed = 0;
	nsource = 0;
	
	nnz = (int *)malloc( sizeof(int) * new_n );
	nzl = (int **)malloc( sizeof(int*) * new_n );
	nzl_space = (int *)malloc( sizeof(int*) * new_n );
	diage = (int *)malloc( sizeof(int) * new_n );
	diagonal_element = (double *)malloc( sizeof(double) * new_n );
	memset( diagonal_element, 0, sizeof(double) * new_n );
	nzv = (double **)malloc( sizeof(double*) * new_n );

	scratch = (double *)malloc( sizeof(double) * new_n );
	memset( nnz, 0, sizeof(int) * new_n );

	for( int x = 0; x < new_n; x++ )
	{
		diage[x] = -1;
		nzl_space[x] = 10;

		nzl[x] = (int *)malloc( sizeof(int) * nzl_space[x] );
		nzv[x] = (double *)malloc( sizeof(double) * nzl_space[x] );

	}
	
}

void SparseMatrix::freeMem( void )
{
	for( int x = 0; x < n; x++ )
	{
		free(nzl[x]);
		free(nzv[x]);
	}
	
	free(nnz);
	free(nzl);
	free(nzl_space);
	free(diage);
	free(nzv);
	free(scratch);
	
	free(diagonal_element);
}

void GramSchmidt( double *vecs, int start, int stop, int n )
{
	for( int v = start; v < stop; v++ )
	{
		for( int v2 = 0; v2 < v; v2++ )
		{
			double dp = 0;

			for( int k = 0; k < n; k++ )
				dp += vecs[v*n+k] * vecs[v2*n+k];

			for( int k = 0; k < n; k++ )
				vecs[v*n+k] -= dp * vecs[v2*n+k];
		} 
		
		double nrm = 0;
		for( int k = 0; k < n; k++ )
			nrm += vecs[v*n+k] * vecs[v*n+k];
		for( int k = 0; k < n; k++ )
			vecs[v*n+k] /= sqrt(nrm);
	}
}

void SparseMatrix::test( void )
{
	double *tmat = ( double *)malloc(sizeof(double) * n * n );
	for( int p = 0; p < n; p++ )
	{
		for( int q = p; q < n; q++ )
		{
			double t = (double)rand()/(double)RAND_MAX - 0.5;
	
			if( p == q )
				coupleParameters( p, q, t*0.5 );
			else
				coupleParameters( p, q, t );
			tmat[p*n+q] = t;
			tmat[q*n+p] = t;
		}
	}
	
	sort();

	for( int p = 0; p < n; p++ )
	{
		printf("{");
		for( int q = 0; q < n; q++ )
		{
			printf("%lf", tmat[p*n+q] );
			if( q != n-1) printf(",");
		}
		printf("}");
		if( p != n-1 )
			printf(",");
		printf("\n");
	}
		
	{
		char jobz = 'V';
		char uplo = 'U';
		int lwork = -1;
		double owork = 0;
		int info = 0;
		int order = n;
		double ev[n];
	
		dsyev(&jobz,&uplo,&order,tmat,&order,ev,&owork,&lwork,&info);
		double *work = (double *)malloc( sizeof(int)*(int)lround(owork));
		
		lwork = lround(owork);
	
		dsyev(&jobz,&uplo,&order,tmat,&order,ev,work,&lwork,&info);		
		free(work);

		printf("ev: %lf %lf %lf\n", ev[0], ev[1], ev[2] );
	}
	int nreq = 2;
	double out_vector[n*nreq];
	double out_ev[nreq];

	davidson( out_vector, out_ev, nreq );

	for( int t = 0; t < nreq; t++ )
	{
		printf("eigenvector: ");
		for( int p = 0; p < n; p++ )
			printf("%lf ", out_vector[t*n+p] );
		printf("\n");
		printf("evalue: %lf\n", out_ev[t] );
	}
}

void SparseMatrix::davidson( double *vectors, double *ev_out, int nvec )
{	
	FILE *theLog = fopen("davidson.log","w");

	double omega = 0.1;
	int max_subspace_add = 20*nvec;

	if( max_subspace_add + nvec > n )
		max_subspace_add = n - nvec;

	double *vspace = (double *)malloc( sizeof(double) * (nvec+max_subspace_add) * n );
	double *vspace_temp = (double *)malloc( sizeof(double) * (nvec+max_subspace_add) * n );
	double *approxe = (double*) malloc( sizeof(double) * (nvec+max_subspace_add) );
	double *residuals = (double *)malloc( sizeof(double) * (nvec+max_subspace_add) * n );
	double *correction = (double *)malloc( sizeof(double) * (nvec+max_subspace_add) * n );

	/* these are the guess eigenvectors */

	for( int v = 0; v < nvec; v++ )
	{
		double nrm = 0;
		for( int k = 0; k < n; k++ )
		{
			vspace[v*n+k] = (double)rand()/(double)RAND_MAX - 0.5;
			nrm += vspace[v*n+k] * vspace[v*n+k];
		}

		for( int k = 0; k < n; k++ )
			vspace[v*n+k] /= sqrt(nrm);
	}

	GramSchmidt( vspace, 0, nvec, n );

	// form the subspace matrix.

	double *W = (double *)malloc( sizeof(double) * n * (nvec+max_subspace_add) );
	double *H = (double *)malloc( sizeof(double) * (nvec+max_subspace_add) * (nvec+max_subspace_add) );
	double *ev = (double *)malloc( sizeof(double) * (nvec+max_subspace_add ) );
	double *temp = (double *)malloc( sizeof(double) *  n );
	
	int done = 0;

	int M = nvec;
	int cur_extra = 0;

	int iter = 0;
	while(!done) 
	{
		// H is the Rayleigh matrix.

		for( int v1 = 0; v1 < M; v1++ )
		{
			memset( temp, 0, sizeof(double) * n );
			mvec( vspace+v1*n, W + v1 * n );
	
			for( int v2 = 0; v2 < M; v2++ )
			{
				H[v1*M+v2] = 0;
				
				for( int k = 0; k < n; k++ )
					H[v1*M+v2] += vspace[v2*n+k] * W[v1*n+k];						
			}
		}
	
		// diagonalize the subspace matrix.
	
		char jobz = 'V';
		char uplo = 'U';
		int lwork = -1;
		double owork = 0;
		int info = 0;
		int order = M;
	
		dsyev(&jobz,&uplo,&order,H,&order,ev,&owork,&lwork,&info);
		double *work = (double *)malloc( sizeof(double)*(int)lround(owork));
		
		lwork = lround(owork);
	
		dsyev(&jobz,&uplo,&order,H,&order,ev,work,&lwork,&info);		
		free(work);

		fprintf(theLog, "cur_ev: ");
		for( int v1 = 0; v1 < nvec; v1++ )
		{
			fprintf(theLog, "%lf ", ev[v1] );
			ev_out[v1] = ev[v1];
		}
		fprintf(theLog,"\n");
		fflush(theLog);

		// H now contains the eigenvectors of the Rayleigh  matrix.

		for( int v1 = 0; v1 < nvec; v1++ )
		{
			memset( vspace_temp+v1*n, 0, sizeof(double) * n );
			approxe[v1] = ev[v1];

			for( int v2 = 0; v2 < M; v2++ )
			{
				for( int k = 0; k < n; k++ )
				{
					vspace_temp[v1*n+k] += H[v1*M+v2] * vspace[v2*n+k];
				}
			}
		}

		// vspace_temp now has the ritz vectors.

		for( int v1 = 0; v1 < nvec; v1++ )
		{
//			mvec( vspace_temp + v1 * n, residuals+v1 * n );
			memset( residuals+v1*n, 0, sizeof(double) * n );
			for( int m = 0; m < M; m++ )
			for( int k = 0; k < n; k++ )
				residuals[v1*n+k] -= W[m*n+k] * H[v1*M+m];

			// now has the ritz vector.

			for( int k = 0; k < n; k++ )
				residuals[v1*n+k] += ev[v1] * vspace_temp[v1*n+k]; 
		}

		double sum_del_2 = 0;
	
		if( M + nvec > nvec + max_subspace_add )
		{	
			fprintf(theLog,"restart\n");
			fflush(theLog);
			// restart

			for( int v1 = 0; v1 < nvec; v1++ )
			{
				//solve2( residuals+v1 * n, correction+v1*n, omega ); 
			
				for( int k = 0; k < n; k++ )
				{
					correction[v1*n+k] = residuals[v1*n+k] / (diage[k] - ev[v1]);  
					sum_del_2 += correction[v1*n+k] * correction[v1*n+k];
				}
				double nrm = 0;	
				for( int k = 0; k < n; k++ )
				{
					vspace[v1*n+k] = vspace_temp[v1*n+k];
					vspace[(v1+nvec)*n+k] = correction[v1*n+k];
				}
			}

			M = 2 * nvec;
		}
		else
		{ 
			// add the corrections to our subspace.

			for( int v1 = 0; v1 < nvec; v1++ )
			{
				//solve2( residuals+v1 * n, correction+v1*n, omega); 
				for( int k = 0; k < n; k++ )
				{
					correction[v1*n+k] = residuals[v1*n+k] / (nzv[k][diage[k]] - ev[v1]);  
					sum_del_2 += correction[v1*n+k] * correction[v1*n+k];
				}

				for( int k = 0; k < n; k++ )
					vspace[(M+v1)*n+k] = correction[v1*n+k];
			}

			M += nvec;
		}
			
		GramSchmidt( vspace, 0, M, n );
		fprintf(theLog,"residual sum: %le\n", sum_del_2 );
		if( sum_del_2 < (1e-7) )
			break;
	}
	
	memcpy( vectors, vspace, sizeof(double) * nvec * n );

	for( int v = 0; v < nvec; v++ )
	{
		printf("Eigenvalue: %le\n", ev_out[v] );
		printf("Eigenvector:");
		for( int k = 0; k < n; k++ )
			printf(" %lf", vectors[v*n+k] );
		printf("\n");
	}
	fclose(theLog);
}

void SparseMatrix::need( int m )
{
	// make sure it's not already there.

	int gotit = 0;

	for( int tm = 0; tm < nneed; tm++ )
	{
		if( need_list[tm] == m ) 
			gotit=1;
	}

	if(!gotit)
	{
		need_list[nneed] = m;
		nneed++;
	}	
}

void SparseMatrix::source( int m )
{
	// make sure it's not already there.

	int gotit = 0;

	for( int tm = 0; tm < nsource; tm++ )
	{
		if( source_list[tm] == m ) 
			gotit=1;
	}

	if(!gotit)
	{
		source_list[nsource] = m;
		nsource++;
	}

	if( nsource > n )
	{
		for( int x = 0; x < nsource; x++ )
		for( int y = x+1; y < nsource; y++ )
		{
			if( source_list[x] == source_list[y] )
			{
				printf("huh");
			}
		}
		printf("what.\n");
	}

}

void SparseMatrix::compress( void )
{
	int ax = 0;

	int *rev = (int *)malloc( sizeof(int) * n );
	
	for( int a = 0; a < n; a++ )
		rev[a] = -1;

	for( int a = 0; a < nsource; a++ )
		rev[source_list[a]] = a;

	for( int a = 0; a < n; a++ )
	{
		if( nnz[a] > 0 )
		{
			nnz[ax] = nnz[a];

			for( int x = 0; x < nnz[a]; x++ )
				nzl[a][x] = rev[nzl[a][x]];
			nzl[ax] = nzl[a];		
			nzv[ax] = nzv[a];		

			ax++;
		}
	}	
	

	n = nneed;
}

void SparseMatrix::compress_source_vector( double *full_B, double *compr_B, int nfast )
{
	for( int s = 0; s < nsource; s++ )
	{
		compr_B[s*nfast] = full_B[source_list[s]*nfast];
	}
}

void SparseMatrix::expand_target_vector( double *full_B, double *compr_B, int nfast )
{
	for( int s = 0; s < nneed; s++ )
	{
		full_B[need_list[s]*nfast] += compr_B[s*nfast];
	}
}

void SparseMatrix::mult(double *compr_A, double *compr_B, int nfast )
{	
	for( int s = 0; s < nneed; s++ )
	{
		compr_A[s*nfast] = 0;

		for( int xx = 0; xx < nnz[s]; xx++ )
		{
			compr_A[s*nfast] += nzv[s][xx] * compr_B[nzl[s][xx]*nfast];
		}
	}
}


/*
	Amat->coupleParameters( p_base_1+f1, p_base_2+f2, eps_boundary * fac * ( seg1->dR[d1]*seg2->dR[d2] ) );

	Amat->sort();
	double *solv = (double *)malloc( sizeof(double) * nparams );
	Amat->cgsolve( b, solv );

*/

void SparseMult( double *vec_out, double *vec_in, SparseMatrix *Mat )
{
	double *mcopy = (double *)malloc( sizeof(double) * Mat->nsource );
	double *vout  = (double *)malloc( sizeof(double) * Mat->nneed );
	memset( vout, 0, sizeof(double)*Mat->nneed);

	Mat->compress_source_vector( vec_in,   mcopy, 1 );
	Mat->mult( vout, mcopy, 1 );
	Mat->expand_target_vector( vec_out, vout, 1 );

	free(mcopy);
	free(vout);
}

void MMUL( SparseMatrix *in1, SparseMatrix *in2, SparseMatrix *out, double tol )
{
	for( int tA = 0; tA < in1->n; tA++ )
	{
		for( int ppA = 0; ppA < in1->nnz[tA]; ppA++ )
		{
			int tA2 = in1->nzl[tA][ppA];

			for( int ppB = 0; ppB < in2->nnz[tA2]; ppB++ )
			{
				int tB = in2->nzl[tA2][ppB];
			
				double val =  in1->nzv[tA][ppA] * in2->nzv[tA2][ppB];

				out->coupleParameters( tA, tB, val, tol );
			}
		}
	}
}

#if 0
void SparseMatrix::SquareRoot( SparseMatrix **root_output )
{
	double *preconditioner = (double *)malloc( sizeof(double) * n );

	double RMS = 0;

	for( int x = 0; x < n; x++ )
		RMS += diagonal_element[x]*diagonal_element[x];

	RMS /= n;
	RMS = sqrt(RMS);

	for( int x = 0; x < n; x++ )
		preconditioner[x] = 1.0 / sqrt( RMS );	

	// precondition the matrix.
	for( int t = 0; t < n; t++ )
		diagonal_element[t] *= preconditioner[t] * preconditioner[t];

	for( int t = 0; t < n; t++ )
	{
		for( int pp = 0; pp < nnz[t]; pp++ )
		{
			int b = nzl[t][pp];

			nzv[t][pp] *= preconditioner[t] * preconditioner[b];

			if( t == b )
				nzv[t][pp] -= 1;
		}
	}

	SparseMatrix *temp = (SparseMatrix *)malloc( sizeof(SparseMatrix ) );
	temp->init(n);
	
	SparseMatrix *temp2 = (SparseMatrix *)malloc( sizeof(SparseMatrix ) );
	temp2->init(n);

	SparseMatrix *result = (SparseMatrix *)malloc( sizeof(SparseMatrix) );
	result->init( n );

	/* ZEROth order: unit matrix */

	for( int t = 0; t < n; t++ )
		result->coupleParameters( t, t, 1.0 / preconditioner[t] );

	int order_limit = 2;

	SparseMatrix *output = NULL;
	SparseMatrix *cur = this;
	int o_switch = 0;

	for( int o = 1; o <= order_limit; o++ )
	{
		if( o > 1 )
		{
			cur->clear();
			MMUL( output, this, cur );
		}
		double factor = sp_fac(2*o) / ( (1 - 2. * o) * pow(sp_fac(o),2)*pow(4,o));

		if( o % 2 == 1 )
			factor *= -1;


		for( int t = 0; t < cur->n; t++ )
		{
			for( int pp = 0; pp < cur->nnz[t]; pp++ )
			{
				int t2 = cur->nzl[t][pp];
	
				// un-precondition.
				double pc = 1.0 / sqrt(preconditioner[t]) / sqrt( preconditioner[t2] );
				
				result->coupleParameters( t, t2, cur->nzv[t][pp] * pc * factor );
			}
		}

		o_switch = !o_switch;

		if( o == 1 )
		{
			output = this;
			cur = temp;
		}
		else if( o_switch )
		{ // o == 3
			output = temp2;
			cur = temp;
		}
		else
		{ // o == 2, 4
			output = temp;
			cur = temp2;
		}
	}
/*
	// SECOND order 
	
	MMUL( this, this, temp );

*/
	// replace matrix.
	
	for( int t = 0; t < n; t++ )
		diagonal_element[t] /= preconditioner[t] * preconditioner[t];

	for( int t = 0; t < n; t++ )
	{
		for( int pp = 0; pp < nnz[t]; pp++ )
		{
			int b = nzl[t][pp];

			if( t == b )
				nzv[t][pp] += 1;

			nzv[t][pp] /= preconditioner[t] * preconditioner[b];
		}
	}


	*root_output = result;

	temp->freeMem();
	temp2->freeMem();

	delete temp;
	delete temp2;

	free(preconditioner);
}
#endif

void SparseMatrix::SquareRoot( SparseMatrix **root_output )
{
	double *preconditioner = (double *)malloc( sizeof(double) * n );

	double maxv = 0;

	for( int x = 0; x < n; x++ )
	{
		if( diagonal_element[x] > maxv )
			maxv = diagonal_element[x];
	}

	for( int x = 0; x < n; x++ )
		preconditioner[x] = 1.0 / maxv;	

	// precondition the matrix.
	for( int t = 0; t < n; t++ )
		diagonal_element[t] *= preconditioner[t];

	for( int t = 0; t < n; t++ )
	{
		for( int pp = 0; pp < nnz[t]; pp++ )
		{
			int b = nzl[t][pp];

			nzv[t][pp] *= preconditioner[t];
		}
	}
	


	SparseMatrix *Yk = (SparseMatrix *)malloc( sizeof(SparseMatrix ) );
	Yk->init(n);
	SparseMatrix *Ykp1 = (SparseMatrix *)malloc( sizeof(SparseMatrix ) );
	Ykp1->init(n);
	SparseMatrix *Zk = (SparseMatrix *)malloc( sizeof(SparseMatrix ) );
	Zk->init(n);
	SparseMatrix *Zkp1 = (SparseMatrix *)malloc( sizeof(SparseMatrix ) );
	Zkp1->init(n);
	
	SparseMatrix *temp = (SparseMatrix *)malloc( sizeof(SparseMatrix ) );
	temp->init(n);

	Zk->loadUnit();
	MMUL( this, Zk, Yk );
	
	int done = 0;
	int iters = 0;

	double tol = 1e-5;

	while( !done )
	{
		Ykp1->clear();
		Zkp1->clear();
		temp->clear();

		MMUL( Zk, Yk, temp, tol );

		// compute 3I - Zk Yk
		
		for( int t = 0; t < temp->n; t++ )
		{
			temp->diagonal_element[t] = 3 - temp->diagonal_element[t];
			for( int x = 0; x < temp->nnz[t]; x++ )
			{
				int t2 = temp->nzl[t][x];

				if( t2 == t )
					temp->nzv[t][x] = 3 - temp->nzv[t][x];
				else
					temp->nzv[t][x] *= -1;
			}
		}

		MMUL( Yk, temp, Ykp1, tol );
		MMUL( temp, Zk, Zkp1, tol );

		Ykp1->scale(0.5);		
		Zkp1->scale(0.5);
		// swap

		SparseMatrix *ty = Yk;
		SparseMatrix *tz = Zk;

		Zk = Zkp1;
		Yk = Ykp1;

		Ykp1 = ty;
		Zkp1 = tz;

		iters++;

		if( iters > 5 )
			done = 1;
	}

	*root_output = Yk;

	temp->freeMem();
	Ykp1->freeMem();
	Zk->freeMem();
	Zkp1->freeMem();

	free(temp);
	free(Ykp1);
	free(Zk);
	free(Zkp1);


	// unprecondition the matrix.
	for( int t = 0; t < n; t++ )
	{
		Yk->diagonal_element[t] /= sqrt(preconditioner[t]);
		diagonal_element[t] /= preconditioner[t];
	}
	for( int t = 0; t < n; t++ )
	{
		for( int pp = 0; pp < nnz[t]; pp++ )
		{
			int b = nzl[t][pp];

			nzv[t][pp] /= preconditioner[t];
		}
	}
	
	for( int t = 0; t < Yk->n; t++ )
	{
		for( int pp = 0; pp < Yk->nnz[t]; pp++ )
		{
			int b = Yk->nzl[t][pp];

			Yk->nzv[t][pp] /= sqrt(preconditioner[t]);
		}
	}

	free(preconditioner);
}
