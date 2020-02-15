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

#define FIXED_SEED

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
		printf("Syntax: coll_sandbox mesh1 mesh2\n");
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
	
	FILE* limit1 = fopen("surface1.xyz","w");
	FILE* limit2 = fopen("surface2.xyz","w");
//	FILE* psf1 = fopen("surface1.psf","w");
//	FILE* psf2 = fopen("surface2.psf","w");
		
//	theSurface1->writeLimitPSF(psf1);
//	theSurface2->writeLimitPSF(psf2);
	
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

//#define COLLISION_DEBUG
 
#ifdef COLLISION_DEBUG

	for( int v = 1; v < nv1; v++ )
	{	
		r1[3*v+0] -= r1[0];
		r1[3*v+1] -= r1[1];
		r1[3*v+2] -= r1[2];
	}
	
	r1[0] = 0;
	r1[1] = 1;
	r1[2] = 2;

	for( int v = 1; v < nv1; v++ )
	{	
		r2[3*v+0] -= r2[0];
		r2[3*v+1] -= r2[1];
		r2[3*v+2] -= r2[2];
	}

	r2[0] = 0;	
	r2[1] = 0;	
	r2[2] = 0;	

	for( int v = 0; v < nv1; v++ )
	{
		r1[3*v+0] += 0.02;
		r1[3*v+1] += 0.04;
		r1[3*v+2] -= 0.01;

		r1[3*v+0] = - r1[3*v+0];
		r1[3*v+1] =   r1[3*v+1];
		r1[3*v+2] = - r1[3*v+2];
	}

#endif

	int nsteps = 100;
	int start_at_step = 15;
		
	for( int v = 0; v < nv1; v++ )
		r1[3*v+0] -= start_at_step * (double)(sep/2) / (double)nsteps;
	for( int v = 0; v < nv2; v++ )
		r2[3*v+0] += start_at_step * (double)(sep/2) / (double)nsteps;


	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 ); 
	int mlow=5,mhigh=7;
	double *M[3]={M5,M6,M7};

	theSurface1->generateSubdivisionMatrices( M, 5, 7 );

	int was_collision = 0;
	for( int step = start_at_step; step < nsteps && !was_collision; step++ )
	{	
		printf("%d \n", step);
		for( int v = 0; v < nv1; v++ )
			r1[3*v+0] -= (double)(sep/2) / (double)nsteps;
		for( int v = 0; v < nv2; v++ )
			r2[3*v+0] += (double)(sep/2) / (double)nsteps;

		theSurface1->put(r1);
		theSurface2->put(r2);
		
		theSurface1->writeLimitingSurface(limit1);
		theSurface2->writeLimitingSurface(limit2);
		fflush(limit1);
		fflush(limit2);

		if(  collision( theSurface1, theSurface2, M, mlow, mhigh ) )
		{
			printf("COLLISION at frame %d! stopping.\n", step );
			was_collision = 1;
		} else {
			printf("nothing here");
		}
			
	}

	fclose(limit1);
	fclose(limit2);	
}





