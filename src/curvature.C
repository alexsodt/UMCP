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
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "Point.h"
#include "bd_functions.h"
#include <random>
#include "Bin.h"

#define FIXED_SEED


int main( int argc, char **argv )
{
          struct timeval tp;
#ifdef FIXED_SEED
	srand(1);
#else
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);
#endif

	char buffer[4096];


	if( argc < 2 )
	{
		printf("Syntax: bd_sandbox mesh\n");
		return -1;
	}

//	const char* input = argv[2];

	// setup surface
	surface *theSurface1 =(surface *)malloc( sizeof(surface) );
	theSurface1->loadLattice( argv[1], 0. );
	
	double com1[3] = {0,0,0};

	theSurface1->generatePlan();
	
	FILE* limit1 = fopen("surface1.xyz","w");
	FILE* psf1 = fopen("surface1.psf","w");
        theSurface1->writeLimitingSurface(limit1);
        theSurface1->writeLimitingSurfacePSF(psf1);
        fflush(limit1);
	fclose(psf1);

	int nv1 = theSurface1->nv;

	int nc1 = 3 * nv1+3;

	double *r1 = (double *)malloc( sizeof(double) * nc1 );

	theSurface1->get(r1);
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

	double *M5 = (double *)malloc( sizeof(double) * 4 * 11 * 12 ); 
	double *M6 = (double *)malloc( sizeof(double) * 4 * 12 * 12 ); 
	double *M7 = (double *)malloc( sizeof(double) * 4 * 13 * 13 ); 
	double *M[3] = { M5, M6, M7};
	int mlow = 5, mhigh = 7;
	theSurface1->generateSubdivisionMatrices( M, mlow, mhigh );

	double LA = theSurface1->PBC_vec[0][0];
	double LB = theSurface1->PBC_vec[1][1];
	double LC = theSurface1->PBC_vec[2][2];


	FILE *file;
	file = fopen("curvature.txt", "w+");
	theSurface1->checkCurvature(file);
	fclose(file);



}





