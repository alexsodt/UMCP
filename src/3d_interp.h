#ifndef __3D_INTERP__
#define __3D_INTERP__

void setup_rho( double *rho, int nx, int ny, int nz  );
double eval_rho(double fx, double fy, double fz );
void eval_drho(double fx, double fy, double fz, double *g );

#endif
