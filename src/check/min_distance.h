#ifndef MIN_DISTANCE_H
#define MIN_DISTANCE_H

#include <math.h>
#include "Point.h"
#include "Check.h"

class Min_Distance {
	private:
		Point* mem;
		Point part;
		int numM;

		// Number of points in simplex
		int elements;

		Point* simplex;

		double r;

	public:
		Min_Distance(double*& objA, double*& objB, int numA, double radius);

		~Min_Distance();

		double dot(Point& a, Point& b);

		Point cross(Point& a, Point& b);

		Point support(Point direction);

	inline Point normalize_direction(Point direction) {
		double magnitude = sqrt(pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2));
		return direction/magnitude;
		}

		Point findDirectionLine();

		Check findDirectionTriangle();

		Check findDirectionTetrahedral();

		Check NearestSimplex();
		
		bool run_GJK();


};
double near_point_on_triangle( double v1x, double v1y, double v1z,
                                double v2x, double v2y, double v2z,
                                double v3x, double v3y, double v3z, int *type, int *which );

#endif
