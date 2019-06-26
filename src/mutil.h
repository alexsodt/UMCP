#ifndef __mutilh__
#define __mutilh__
double dot( double *a, double *b );
double cross( double *dr1, double *dr2, double *cp );
double normalize( double *dr );
double length3( double *dr );
void d_triangle_area( double *der, double *pt1, double *pt2, double *pt3 );
double triangle_area( double *pt1, double *pt2, double *pt3 );
void dnrm( double *g, double *pt1, double *pt2, double *pt3 );
void d_tvol( double *dept, double *pt1, double *pt2, double *pt3, double *pt4 );
double tvol( double *pt1, double *pt2, double *pt3, double *pt4 );
int pointInTriangle( double *pt, double *pt1, double *pt2, double *pt3 , double tol );
void MinImage3D( double *dr1, double PBC_vec[3][3], double *put_vec, double *alphas=NULL );
void fprintPSFAtomExt( FILE *theFile, int atn, int res, const char *atname, const char *resname);
void writePSF( FILE *thePSF, int nv, char *atomNames, int *bond_list, int nbonds );
int sort3(int*ind);
int getcode( int *verts, int *sorted );
int iabs(int);
double Power( double, double );
double Power( double, int );
double Mul22( double *A, double *B, double *out );
void CartMatVecIncrScale( double *vec_out, double *vec_in, double *Mat, double scale, int nv, double *alphas  );
void MatVec( double *a, double *b, double *c, int m, int n);
#endif
