#ifndef __sparseh__
#define __sparseh__

struct SparseMatrix
{
	int n;

	int *need_list;
	int nneed;

	int *source_list;
	int nsource;

	int *diage;
	int *nnz;
	int **nzl;
	int *nzl_space;
	double *scratch;
	double **nzv;

	void sort( void );
	void init( int new_n );
	void solve( double *b, double *soln, double omega );
	void solve2( double *b, double *soln, double omega );
	void davidson( double *vectors, double *ev, int n );
	void cgsolve( double *b, double *soln );
	void coupleParameters( int p1, int p2, double val );
	void debug( double *A );
	void test( void );
	double melem( double *b );
	void mvec( double *x, double *b );
	void freeMem( void );

	void need( int m );
	void source( int m );
	void setNeedSource( void );
	void compress( void );

	void compress_source_vector( double *full_B, double *compr_B, int nfast );
	void expand_target_vector( double *full_A, double *compr_A, int nfast );
	void mult( double *compr_A, double *compr_B, int nfast  );
};

void SparseMult( double *vec_out, double *vec_in, SparseMatrix *Mat );

#endif
