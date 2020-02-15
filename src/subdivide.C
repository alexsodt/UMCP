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
		printf("Syntax: subdivide lattice [refine_at] [edge]\n");
		return -1;
	}
	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( argv[1], 0. );

	surface *sub_surface = (surface *)malloc( sizeof(surface) );
	surface *sub_surface2 = (surface *)malloc( sizeof(surface) );

	int refine_at = -1;
	int refine_edge = 0;

	if( argc > 2 )
	{
		refine_at = atoi(argv[2]);
		if( refine_at >= theSurface->nv )
		{
			printf("refine_at must be less than the number of vertices (zero-offset index).\n");
			exit(1);
		}

		if( argc > 3 )
			refine_edge = atoi(argv[3]);
	}
	if( refine_at >= 0 )
	{
		int new_center, new_edge;
		sub_surface->subdivideTriangle3( theSurface, refine_at, refine_edge, &new_center, &new_edge );
		printf("new center: %d new_edge: %d\n", new_center, new_edge );
	}
	else
		sub_surface->subdivideSurface( theSurface );

//	sub_surface->writeXYZSurface( "subdiv.xyz", "subdiv.psf", sub_surface);
	sub_surface->writeXYZandPSFPeriodic( "subdiv" );
	sub_surface->saveSurface("subdiv.mesh");
//	sub_surface->writeLimitingSurface("limit.xyz");
	sub_surface->generatePlan();
//	FILE* limit = fopen("limit.xyz","w");
//	sub_surface->writeStructure(limit);
//	fclose(limit);
	int nv = sub_surface->nv;
	int nc = 3 * nv+3;
	double *r = (double *)malloc( sizeof(double) * nc );
	double *g = (double *)malloc( sizeof(double) * nc );
	sub_surface->get(r);
	r[3*nv+0] = 1.;
	r[3*nv+1] = 1.;
	r[3*nv+2] = 1.;
	sub_surface->setg0(r);
	double eps = 1e-6;

	
}





