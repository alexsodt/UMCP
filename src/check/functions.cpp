#define functions_cpp
#include "functions.h"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "NearestSimplex.cpp"


using namespace std;

Functions::Functions(double*& objA, double*& objB, int numA, int numB) {
	elements = 0;
	simplex = new Point[4];
	
	int temp_size = numA * numB;
    md = new Point[temp_size];
	num = 0;
    for (int i = 0; i < numA; i++) {
        for (int j = 0; j < numB; j++) {
            Point diff;
			diff.x = objA[3*i + 0] - objB[3*j + 0];
			diff.y = objA[3*i + 1] - objB[3*j + 1];
			diff.z = objA[3*i + 2] - objB[3*j + 2];  
    
                md[num] = diff;
                num++;
            }    
        }   
    }   

Functions::~Functions() {
	delete[] md;
	delete[] simplex;
}


Point Functions::cross(Point& a, Point& b) {
	Point vec;
	vec.x = (a.y * b.z - a.z * b.y);
	vec.y = - (a.x * b.z - a.z * b.x);
	vec.z = (a.x * b.y - a.y * b.x);
	return vec;
}

Point Functions::support(Point*& data, int num, Point& direction) {
	double max_value = dot(data[0], direction);
	int index = 0;
	for (int i = 1; i < num; i++) {
		double temp = dot(data[i], direction);
		if (temp > max_value) {
			max_value = temp;
			index = i;
		} 
	}	
	
	return data[index];
}

Point Functions::findDirectionLine() {
	// find vectors AB from new point to old and AO from new point to origin
	Point AB = simplex[0] - simplex[1];
	Point AO = negative(simplex[1]);

	if (dot(AB, AO) > 0) {
		Point cross1 = cross(AB, AO);

// RETURN NEGATIVE???
		return (cross(cross1, AB));
//	} else {
//		cout << "weird" << endl;
//		return AO;
	}
	return AO;
}

Check Functions::findDirectionTriangle() {
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
			Point cross1 = cross(AC, AO);
			Check new_check = {false,  cross(cross1, AC)};
			return new_check;
		} else {
			if (dot(AB, AO) > 0) {
				//region 4
				simplex[0] = simplex[1];
				simplex[1] = simplex[2];
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
                Point cross2 = cross(AB, AO);
                Check new_check = {false, cross(cross2, AB)};
                return new_check;
            } else {
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


Check Functions::findDirectionTetrahedral() {
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


bool Functions::run_GJK() {


	// initial search direction (x-axis)
	Point direction = {1, 0, 0};
	Point new_support = support(md, num, direction);
	simplex[elements] = new_support;
	elements++;

    // find point in opposite direction
    direction = negative(direction);

    bool intersect = true;
    int repeat = 0;
    while (intersect == true) {
        Point new_point = support(md, num, direction);
    
        if (dot(new_point, direction) < 0) {
           // cout << "No intersection." << endl;
            intersect = false;
			return false;
        } else {
            simplex[elements] = new_point;
            elements++;

            if (elements == 2) {
                direction = findDirectionLine();
//		cout << "line" << endl;
            } else if (elements == 3) {
				Check searchTri = findDirectionTriangle();
				if (searchTri.inside == true) {
					direction = searchTri.direction;
				} else {
				    direction = searchTri.direction;
					elements--;	
				}

//		cout << "triangle" << endl;
            } else if (elements == 4) {
			Check checkTetra = NearestSimplex();
			if (checkTetra.inside == true) {
				return true;
			} else {
			    	direction = checkTetra.direction;
			}
			
            }   
        }   
    }   
	return false;
}


