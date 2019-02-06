#ifndef __clusterh__
#define __clusterh__

double loopArea( double *pts, int npts );
double triangle_area( double *pt1, double *pt2, double *pt3 );
void get2DVoronoiConnectivity( double *pts, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds, double BoxL[2] );
void get2DHexVoronoiConnectivity( double *pts, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds, double BoxL[2] );
void getSphericalVoronoiConnectivity( double *pts, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds );
void get3DVoronoi( double *pts, int npts, const char *uniq, int *bonds, int *nbonds, int max_bonds, double PBC_vec[3][3], double max_cut );

#endif
