#ifdef LAPACK_UNDERSCORE
	#define dgeev dgeev_
	#define dgetrf dgetrf_
	#define dgetri dgetri_
	#define dposv dposv_
	#define dgesv dgesv_
	#define dgemv dgemv_
	#define dgemm dgemm_
	#define dsyev dsyev_
	#define dsytrf dsytrf_
	#define dsytri dsytri_
	#define dpotri dpotri_
	#define dpotrf dpotrf_
#endif

extern "C" void dgeev( char *jobvl, char *jobvr, int *N, double *A, int *LDA, double *WR, double *WI, double *VL, int *LDVL, double *VR, int *LDVR, double *work, int *lwork, int *info );
extern "C" void dgetrf( int *, int *, double *, int *, int *, int * ); 
extern "C" void dgetri( int *, double *, int *, int *, double *, int *, int * ); 
extern "C" void dposv( char *uplo, int *n, int *nrhs, double *A, int *lda, double *b, int *ldb, int *info );
extern "C" void dgesv( int *n, int *nrhs, double *A, int *lda, int *ipiv, double *b, int *ldb, int *info );
extern "C" void dgemv( char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incrx, double *beta, double *y, int *incy );
extern "C" void dgemm( char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc );
extern "C" void dsyev( char *jobz, char *uplo, int *N, double *A, int *LDA, double *WR, double *work, int *lwork, int *info );
extern "C" void dsytrf( char *, int *, double *, int *, int *,  double *, int *, int *);
extern "C" void dsytri( char *, int *, double *, int *, int *,  double *, int *);
extern "C" void dpotri_( char *, int *, double *, int *, int *);
extern "C" void dpotrf_( char *, int *, double *, int *, int *);

extern "C" void dgeev( char *jobvl, char *jobvr, int *N, double *A, int *LDA, double *WR, double *WI, double *VL, int *LDVL, double *VR, int *LDVR, double *work, int *lwork, int *info );
