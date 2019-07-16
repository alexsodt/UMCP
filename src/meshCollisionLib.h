#ifndef __meshcollisionlibh__
#define __meshcollisionlibh__

void surfaceSurfaceCollisionForces( surface *surface1, surface *surface2, double *grad1, double *grad2, double alpha, double v0, double **M, int mlow, int mhigh );

#endif
