//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "l-bfgs.h"
#include "m_triangles.h"
//#define SKIP_OPT

#define FIX_VOLUME

surface *minSurface = NULL;
double thick = 13.8;
extern int do_print_area_ps;
extern double av_vol_err;
extern double KV_scale;
void fdiff_check( double *coords );
int print_energy = 0;
extern double tilt_scale;
surface *upperSurface;
surface *lowerSurface;
volumes *theVolumes; 
double surface_grad( double *r, double *g );
double surface_e( double *r );
	void curvatureGrad( double *g, surface *surface1, surface *surface2, double *r, double thick );
	double curvatureEnergy( surface *surface1, surface *surface2, double *r, double thick);
	void tiltGrad( double *g, surface *surface1, surface *surface2, double *r );
	void tiltGradExt( double *g, surface *surface1, surface *surface2, double *r );
	double tiltEnergy( surface *surface1, surface *surface2, double *r);
	double tiltEnergyExt( surface *surface1, surface *surface2, double *r, double *cons_energy);
double f(double *r);
double fdf( double *r, double *g);
void writeDualSurface( const char *fileNameXYZ, const char *fileNamePSF, surface *upperSurface, surface *lowerSurface, double *coords, double alpha );

int fix_alpha = 0;

int main( int argc, char **argv )
{
#ifdef FIXED_SEED
	srand(1);
#else
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);
#endif

	char buffer[4096];

	if( argc < 2 )
	{
		printf("Syntax: scale lattice [scale=1]\n");
		return -1;
	}
	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( argv[1], 0. );
	double scale = atof(argv[2]);

	theSurface->PBC_vec[0][0] *= scale;
	theSurface->PBC_vec[0][1] *= scale;
	theSurface->PBC_vec[0][2] *= scale;
	theSurface->PBC_vec[1][0] *= scale;
	theSurface->PBC_vec[1][1] *= scale;
	theSurface->PBC_vec[1][2] *= scale;
	theSurface->PBC_vec[2][0] *= scale;
	theSurface->PBC_vec[2][1] *= scale;
	theSurface->PBC_vec[2][2] *= scale;

	for( int i = 0; i < theSurface->nv; i++ )
	{
		theSurface->theVertices[i].r[0] *= scale;
		theSurface->theVertices[i].r[1] *= scale;
		theSurface->theVertices[i].r[2] *= scale;
	}

	theSurface->saveSurface("scale.mesh");
}
