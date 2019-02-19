#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "min_distance.h"
#include "mutil.h"

using namespace std;

// full GJK
bool gjk_algorithm(double* a, int numA, double* b, int numB);
// gudrun's version
bool gjk_algorithm_gudrun(double* a, int numA, double* b, int numB);

//double* a;
//int numA;
//double* b;

bool quick_check( double *r1, int nv1, double *r2, double radius )
{
	// quick check.

	double com[3] = { 
		(r1[0]+r1[3]+r1[6])/3,
		(r1[1]+r1[4]+r1[7])/3,
		(r1[2]+r1[5]+r1[8])/3 };
	double dr[3] = { com[0]-r2[0],
			 com[1]-r2[1],
			 com[2]-r2[2] };

	normalize(dr);
	double mins=1e10;
	double maxs=-1e10;
	double dp2 = r2[0]*dr[0]+r2[1]*dr[1]+r2[2]*dr[2];

	for( int x = 0; x < nv1; x++ )
	{
		double ldp = r1[3*x+0] * dr[0] + r1[3*x+1] * dr[1] + r1[3*x+2] * dr[2];

		if( ldp < mins) mins = ldp;
		if( ldp > maxs ) maxs = ldp;
	}

//	printf("mins: %le maxs: %le ldp: %le radius: %le\n", mins, maxs, dp2, radius );

	if( dp2 - radius > maxs || dp2 + radius < mins )
		return false;

	return true;
}


bool minimum_distance(double* r1, int nv1, double* r2, double radius, int do_quick_abort ) {

	if( ! quick_check( r1, nv1, r2, radius ) )
		return false;

	Min_Distance test(r1, r2, nv1, radius, do_quick_abort);
//	test.write_Minkowski();
	bool trgjk =  test.run_GJK();

	return trgjk;
}

bool linear_collision_worker(double* r1, int nv1, double* pt1, double *pt2, double *pt,  double radius) {

	// we do a quick check from the midpoint between the two particles using the sphere.
	if( ! quick_check( r1, nv1, pt, radius ) )
		return false;

	double pts[6] = { pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2] };

	return gjk_algorithm( r1, nv1, pts, 2 );
}


