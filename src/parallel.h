#ifndef __parallelh__
#define __parallelh__

#include "sparse.h"
struct surface;
struct pcomplex;

#define BASE_TASK 0

/*
	global parallel information

*/

typedef struct parallel_info
{
	int my_id;
	int nprocs;
	
	// list of generalized coordinates we are processing.
	int *genQ;
	int NQ;

	// list of faces we are computing
	int *faces;
	int nf;

	int *proc_for_vert;

	// list of vertices we are responsible for.
	int *verts;
	int nv;

	// list of complexes we are computing
	int *complexes;
	int nc;

	int *nsend_to;
	int **send_to;

	int *nget_from;
	int **get_from;
	
	int *psum_nsend_to;
	int **psum_send_to;

	int *psum_nget_from;
	int **psum_get_from;

} parallel_info;

#ifndef __parallelc__
extern parallel_info par_info;
extern int setup_for_parallel;
#endif

void quietParallel( void );
void setupParallel( surface *theSurface, pcomplex **allComplexes, int ncomplex, int NQ );
void setupSparseVertexPassing( SparseMatrix *EFFM, int nverts, int NQ  );
void ParallelSum( double *vec, int len );
void ParallelSyncComplexes( pcomplex **allComplexes, int ncomplexes );
void ParallelBroadcast( double *vec, int len );
void ParallelBroadcast( int *vec, int len );
void ParallelGather( double *vec, int len_per_proc );
void ParallelGather( int *vec, int len_per_proc );
void FullSyncVertices( double *total_vec );
void PartialSyncVertices( double *total_vec );
void PartialSumVertices( double *total_vec );
void PartialGenVertices( double *total_vec, int do_sum );

void FullSyncGenQ( double *total_vec );

void SparseCartMatVecIncrScale( double *vec_out, double *vec_in, double *Mat, double scale, int nv, int *sparse_use, int nv_use, double *alphas );
void AltSparseCartMatVecIncrScale( double *vec_out, double *vec_in, SparseMatrix *Mat, double scale, double *alphas );
void GenQMatVecIncrScale( double *vec_out, double *vec_in, SparseMatrix *Mat, double scale );

#endif
