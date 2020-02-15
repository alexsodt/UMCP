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
#include "M_matrix.h"
#define FIXED_SEED

void surfaceSurfaceCollisionForces( surface *surface1, surface *surface2, double *grad1, double *grad2, double alpha, double v0, double **M, int mlow, int mhigh );

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
		printf("Syntax: frc_sandbox mesh1 mesh2\n");
		return -1;
	}

	surface *theSurface1 =(surface *)malloc( sizeof(surface) );
	theSurface1->loadLattice( argv[1], 0. );
	
	surface *theSurface2 =(surface *)malloc( sizeof(surface) );
	theSurface2->loadLattice( argv[2], 0. );

	double com1[3] = {0,0,0};
	double com2[3] = {0,0,0};

	theSurface1->generatePlan();
	theSurface2->generatePlan();
	
	
	int nv1 = theSurface1->nv;
	int nv2 = theSurface2->nv;

	int nc1 = 3 * nv1+3;
	int nc2 = 3 * nv2+3;

	double *r1 = (double *)malloc( sizeof(double) * nc1 );
	double *r2 = (double *)malloc( sizeof(double) * nc2 );

	theSurface1->get(r1);
	r1[3*nv1+0] = 1.;
	r1[3*nv1+1] = 1.;
	r1[3*nv1+2] = 1.;

	theSurface2->get(r2);
	r2[3*nv2+0] = 1.;
	r2[3*nv2+1] = 1.;
	r2[3*nv2+2] = 1.;

	for( int v = 0; v < nv1; v++ )
	{
		com1[0] += r1[3*v+0];
		com1[1] += r1[3*v+1];
		com1[2] += r1[3*v+2];
	}
	
	for( int v = 0; v < nv2; v++ )
	{
		com2[0] += r2[3*v+0];
		com2[1] += r2[3*v+1];
		com2[2] += r2[3*v+2];
	}

	com1[0] /= nv1;
	com1[1] /= nv1;
	com1[2] /= nv1;
	
	com2[0] /= nv2;
	com2[1] /= nv2;
	com2[2] /= nv2;

	double sep = 700;

	for( int v = 0; v < nv1; v++ )
	{
		r1[3*v+0] -= com1[0];
		r1[3*v+1] -= com1[1];
		r1[3*v+2] -= com1[2];

		r1[3*v+0] += sep/2;
	}

	for( int v = 0; v < nv2; v++ )
	{
		r2[3*v+0] -= com2[0];
		r2[3*v+1] -= com2[1];
		r2[3*v+2] -= com2[2];

		r2[3*v+0] -= sep/2;
	}

	int nsteps = 100;
	int start_at_step = 15;
		
	for( int v = 0; v < nv1; v++ )
		r1[3*v+0] -= start_at_step * (double)(sep/2) / (double)nsteps;
	for( int v = 0; v < nv2; v++ )
		r2[3*v+0] += start_at_step * (double)(sep/2) / (double)nsteps;


	int mlow=5,mhigh=7;

	double **M;
	getM(&M,&mlow,&mhigh);

	int was_collision = 0;

	double *frc1 = (double *)malloc( sizeof(double) * 3 * nv1 );
	double *frc2 = (double *)malloc( sizeof(double) * 3 * nv2 );
	

	for( int step = start_at_step; step < nsteps && !was_collision; step++ )
	{	
		memset( frc1, 0, sizeof(double) * 3 * nv1 );
		memset( frc2, 0, sizeof(double) * 3 * nv2 );

		printf("%d \n", step);
		for( int v = 0; v < nv1; v++ )
			r1[3*v+0] -= (double)(sep/2) / (double)nsteps;
		for( int v = 0; v < nv2; v++ )
			r2[3*v+0] += (double)(sep/2) / (double)nsteps;

		theSurface1->put(r1);
		theSurface2->put(r2);


		surfaceSurfaceCollisionForces( theSurface1, theSurface2, frc1, frc2, 10.0, 1.0, M, mlow, mhigh ); 		

		printf("GRADIENTS\n");
		for( int v = 0; v < nv1; v++ )
			printf("S1 %d %le %le %le\n", v, frc1[3*v+0], frc1[3*v+1], frc1[3*v+2] );
		for( int v = 0; v < nv2; v++ )
			printf("S2 %d %le %le %le\n", v, frc2[3*v+0], frc2[3*v+1], frc2[3*v+2] );
	}

}





