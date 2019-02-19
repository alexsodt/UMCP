#define min_distance_cpp
#include "min_distance.h"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "NearestSimplex.cpp"


using namespace std;

Min_Distance::Min_Distance(double*& objA, double*& objB, int numA, double radius) {
	elements = 0;
//	simplex = new Point[4];
	numM = numA;
	r = radius;
	
	if( numA > MAX_MD_POINTS )
	{
		printf("ERROR: MD points\n");
		exit(1);
	}

//	mem = new Point[numM];
	for (int i = 0; i < numM; i++) {
		mem[i].x = objA[3*i + 0];
		mem[i].y = objA[3*i + 1];
		mem[i].z = objA[3*i + 2];
	}
	
	part.x = objB[0];
	part.y = objB[1];
	part.z = objB[2];
}   

Min_Distance::~Min_Distance() {
//	delete[] mem;
//	delete[] simplex;
}

double Min_Distance::dot(Point& a, Point& b) {
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

Point Min_Distance::cross(Point& a, Point& b) {
	Point vec;
	vec.x = (a.y * b.z - a.z * b.y);
	vec.y = - (a.x * b.z - a.z * b.x);
	vec.z = (a.x * b.y - a.y * b.x);
	return vec;
}

Point Min_Distance::support(Point direction) {
//direction already normalized
	double mem_support = dot(mem[0], direction);
	int memIndex = 0;
	for (int i = 1; i < numM; i++) {
		double temp = dot(mem[i], direction);
		if (temp > mem_support) {
			mem_support = temp;
			memIndex = i;
		} 
	}	

	Point temp = (negative(direction) * r);
	Point sphere_support = temp + part;

	return mem[memIndex] - sphere_support;
}



Point Min_Distance::findDirectionLine() {
	// find vectors AB from new point to old and AO from new point to origin
	Point AB = simplex[0] - simplex[1];
	Point AO = negative(simplex[1]);

	if (dot(AB, AO) > 0) {
		Point cross1 = cross(AB, AO);

		return (cross(cross1, AB));
	}
	return AO;
}

Check Min_Distance::findDirectionTriangle() {
	Point AB = simplex[1] - simplex[2];
	Point AC = simplex[0] - simplex[2];
	Point AO = negative(simplex[2]);
	Point ABC = cross(AB, AC);
	
	Point temp = cross(ABC, AC);
	Point temp2 = cross(AB, ABC);

	if (dot(temp, AO) > 0) {
		if (dot(AC, AO) > 0) {
			//region 1
			simplex[1] = simplex[2];
			elements--;
			Point cross1 = cross(AC, AO);
			Check new_check = {false,  cross(cross1, AC)};
			return new_check;
		} else {
			if (dot(AB, AO) > 0) {
				//region 4
				simplex[0] = simplex[1];
				simplex[1] = simplex[2];
				elements--;
				Point cross2 = cross(AB, AO);
				Check new_check = {false, cross(cross2, AB)};
				return new_check;
		//	} else {
		//		//region 5
		//		cout << "weird" << endl;
		//		return AO;
			}
	
		}
	} else {
		if (dot(temp2, AO) > 0) {
            if (dot(AB, AO) > 0) {
                //region 4
                simplex[0] = simplex[1];
                simplex[1] = simplex[2];
		elements--;
                Point cross2 = cross(AB, AO);
                Check new_check = {false, cross(cross2, AB)};
                return new_check;
        //    } else {
         //       //region 5
		//		cout << "weird" << endl;
          //      return AO; 
            }   				
		} else {
			if (dot(ABC, AO) > 0) {
				Check new_check = {true, ABC};
				return new_check;
			} else {
				Check new_check = {true, negative(ABC)};
				return new_check;
			}
		}
	}
	Check new_check = {true, AO};
	return new_check;
}


Check Min_Distance::findDirectionTetrahedral() {
    Point AB = simplex[2] - simplex[3];
    Point AC = simplex[1] - simplex[3];
    Point AD = simplex[0] - simplex[3];
    Point AO = negative(simplex[3]);

    Point ABC = cross(AB, AC);
    Point ACD = cross(AC, AD);
    Point ADB = cross(AD, AB);

    double kABC = dot(ABC, simplex[3]);
    double kACD = dot(ACD, simplex[3]);
    double kADB = dot(ADB, simplex[3]);

    if (dot(ABC, simplex[0]) * kABC > 0) {
        //right face
        simplex[0] = simplex[1];
        simplex[1] = simplex[2];
        simplex[2] = simplex[3];

	if (kABC >= 0) {
	 Check new_check = {false, negative(ABC)};
	 return new_check;
	} else if (kABC < 0) {
	 Check new_check = {false, ABC};
	 return new_check;
	}

        Check new_check = {false, ABC};
        return new_check;
    } else if (dot(ACD, simplex[2]) * kACD > 0) {
        //left face
        simplex[2] = simplex[3];

        if (kACD >= 0) {
         Check new_check = {false, negative(ACD)};
         return new_check;
        } else if (kACD < 0) {
         Check new_check = {false, ACD};
         return new_check;
        }

        Check new_check = {false, ACD};
        return new_check;
    } else if (dot(ADB, simplex[1]) * kADB > 0) {
        //front face
        simplex[1] = simplex[2];
        simplex[2] = simplex[3];

        if (kADB >= 0) {
         Check new_check = {false, negative(ADB)};
         return new_check;
        } else if (kADB < 0) {
         Check new_check = {false, ADB};
         return new_check;
        }

        Check new_check = {false, ADB};
        return new_check;
    }
    Check new_check = {true, AO};
    return new_check;

}


bool Min_Distance::run_GJK() {


	// initial search direction (x-axis)
	Point direction = {mem[0].x - part.x, mem[0].y - part.y, mem[0].z - part.z};
	Point new_support = support(normalize_direction(direction));
	simplex[elements] = new_support;
	elements++;

    // find point in opposite direction
    direction = normalize_direction(negative(direction));

    bool intersect = true;
    int repeat = 0;
    while (intersect == true) {
        Point new_point = support(direction);
    
        if (dot(new_point, direction) < 0) {
           // cout << "No intersection." << endl;
            intersect = false;
			return false;
        } else {
            simplex[elements] = new_point;
            elements++;

            if (elements == 2) {
                direction = normalize_direction(findDirectionLine());
//		cout << "line" << endl;
            } else if (elements == 3) {
				Check searchTri = findDirectionTriangle();
//				if (searchTri.inside == true) {
					direction = normalize_direction(searchTri.direction);
//				} else {
				    direction = normalize_direction(searchTri.direction);
//					elements--;	
//				}

//		cout << "triangle" << endl;
            } else if (elements == 4) {
		//		cout << "HERE!!!" << endl;
				Check checkTetra = NearestSimplex();
				if (checkTetra.inside == true) {
					return true;
				} else {
			    	direction = normalize_direction(checkTetra.direction);
				}
			
            }   
        }   
    }   
	return false;
}

