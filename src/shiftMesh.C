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
#include "parallel.h"
#include "simulation.h"
//#define SKIP_OPT



int main( int argc, char **argv )
{
	char buffer[4096];

	if( argc < 5 )
	{
		printf("Syntax: shiftMesh lattice dx dy dz\n");
		return -1;
	}

	double dx = atof(argv[2]);
	double dy = atof(argv[3]);
	double dz = atof(argv[4]);
	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( argv[1], 0. );

	theSurface->generatePlan();

	for( int i = 0; i < theSurface->nv; i++ )
	{
		theSurface->theVertices[i].r[0] += dx;
		theSurface->theVertices[i].r[1] += dy;
		theSurface->theVertices[i].r[2] += dz;
	}
	
	theSurface->saveSurface("shift.mesh");

}
