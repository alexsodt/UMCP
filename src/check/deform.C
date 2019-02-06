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

//#define FIXED_SEED


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
		printf("Syntax: deform lattice\n");
		return -1;
	}
	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( argv[1], 0. );

	double com1[3] = {0,0,0};

	theSurface->generatePlan();
	
	int nv1 = theSurface->nv;

	int nc1 = 3 * nv1+3;

	double *r1 = (double *)malloc( sizeof(double) * nc1 );

	theSurface->get(r1);
	r1[3*nv1+0] = 1.;
	r1[3*nv1+1] = 1.;
	r1[3*nv1+2] = 1.;


	for( int v = 0; v < nv1; v++ )
	{
		com1[0] += r1[3*v+0];
		com1[1] += r1[3*v+1];
		com1[2] += r1[3*v+2];
	}
	
	com1[0] /= nv1;
	com1[1] /= nv1;
	com1[2] /= nv1;
	
	for( int v = 0; v < nv1; v++ )
	{
		r1[3*v+0] -= com1[0];
		r1[3*v+1] -= com1[1];
		r1[3*v+2] -= com1[2];
	}
	
	double area0;
	double cur_area;
	theSurface->area(r1, -1, &cur_area, &area0 );

	double area_per_vertex = cur_area / theSurface->nv;

	double pert = sqrt(area_per_vertex) / 3.0;


	printf("PERT: %le\n", pert );

#define MAKE_SPIKY

#ifdef MAKE_SPIKY
	for( int v = 0; v < nv1; v++ )
	{
		if( rand()/(double)RAND_MAX < 0.5 )
		{
			r1[3*v+0] *= 1.5; 
			r1[3*v+1] *= 1.5; 
			r1[3*v+2] *= 1.5; 
		}
	}

#else
	for( int v = 0; v < nv1; v++ )
	{
		r1[3*v+0] += pert *2 *( rand()/(double)RAND_MAX - 0.5);
		r1[3*v+1] += pert *2 *( rand()/(double)RAND_MAX - 0.5);
		r1[3*v+2] += pert *2 *( rand()/(double)RAND_MAX - 0.5);
	}
#endif
	theSurface->put(r1);
	
	theSurface->saveSurface("deformed.mesh");

	
}





