#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Point.h"
#include "Check.h"
//#include "Check2.h"

class Functions {
	private:
		// Array to hold Minkowski difference 
		Point* md; 

		// Number of points in Minkowski difference
		int num;

		// Number of points in simplex
		int elements;

		Point* simplex;

	public:
		Functions(double*& objA, double*& objB, int numA, int numB);

		~Functions();

		//inline double dot(Point& a, Point& b);
		inline double dot(Point& a, Point& b) {
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z); }

		Point cross(Point& a, Point& b);

		Point support(Point*& data, int num, Point& direction);

		Point findDirectionLine();

		Check findDirectionTriangle();

		Check findDirectionTetrahedral();

		Check NearestSimplex();
		
		bool run_GJK();

double near_point_on_triangle( double v1x, double v1y, double v1z,
                                double v2x, double v2y, double v2z,
                                double v3x, double v3y, double v3z, int *type, int *which );

};

#endif
