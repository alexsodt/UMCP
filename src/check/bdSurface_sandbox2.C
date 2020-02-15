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

// LJ sigma
double sigma = 1;

// LJ cutoff
double rc = pow(2, 1/6) * sigma;

//LJ epsilon
double eps = 1.0;

// diameter of particle
double d = 1.0;

// Length of the box
double L = 650.0 * d;
int l = 650 * d;

// K Spring constant
double k = 1;

// R0
double r0 = 1;

// Energy
double kbT = 2.0;

// Diffusion constant
double D = 7.7;

// Number of particles
int num_part = 1;

// Data holder
Point* data;
Point* temp_data;
Point* orig_data;

// time step
double del_t = .01;

// membrane variables
double* mem;
int numM;

//double*** boxes;
int xBox, yBox, zBox;
double defaultR = 25;

void read_data(const char*& file1) {
    ifstream ifs;
    ifs.open(file1);
    string line;
    getline(ifs, line);
    istringstream hold(line);
    hold >> num_part;

    getline(ifs, line);
    istringstream hold2(line);
    data = new Point[num_part];

    for (int i = 0; i < num_part; i++) {
        getline(ifs, line);
        istringstream hold3(line);
        string name;
        hold3 >> name;
        hold3 >> data[i].x;
        hold3 >> data[i].y;
        hold3 >> data[i].z;
    }  

}

void input_data() {
	data = new Point[num_part];
	BD_Functions bdFunc;
/**	Point part1 = {-10, -20, -180};
	Point part2 = {-70, 30, 230};
	data[0] = part1;
	data[1] = part2;
*/

    bool close = true;

    while( close )
	{
		close = false;
	    for (int i = 0; i < num_part; i++) {
	        data[i].x = rand() % 200 - 100;
	        data[i].y = rand() % 200 - 100;
	        data[i].z = rand() % 200 - 100;
	
	        for (int j = 0; j < i; j++) {
	            if (bdFunc.calc_distance(data[i], data[j]) < rc) {
			close = true;
		    }
		}
	    }
	}

	temp_data = new Point[num_part];
	orig_data = new Point[num_part];
	temp_data = data;
	orig_data = data;
}

void write_data(std::ofstream& file) {
    file << num_part << "\n";
    file << "//comment here" << "\n";

    for (int i = 0; i < num_part; i++) {
	file << "A" << " " << data[i].x << " " << data[i].y << " " << data[i].z << "\n";
    }
}

double * conversion(Point a) {
	double * newPoint = new double[3];
	newPoint[0] = a.x;
	newPoint[1] = a.y;
	newPoint[2] = a.z;
	return newPoint;
}


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

	theSurface1->generateSubdivisionMatrices( M5, M6 );

	double LA = theSurface1->PBC_vec[0][0];
	double LB = theSurface1->PBC_vec[1][1];
	double LC = theSurface1->PBC_vec[2][2];

	// start Brownian on particles
	input_data();
//	read_data(input);
	BD_Functions func;
	
	func.pbc(data, num_part, L);
	double potential = func.calc_potential(data, num_part, sigma, eps, L);
	
	double t = 0; 
	double tmax = 1;

	ofstream file, file2, file3, file4, file5, file6;
	file.open("coords.xyz");
	file3.open("coll.xyz");
	write_data(file);
	file2.open("energy.txt");
	file2 << t << " " << potential << "\n";

        gettimeofday( &tp, NULL );

	double time_1 = tp.tv_usec * (1e-6) + tp.tv_sec;

	theSurface1->box_system(L);

	theSurface1->returnSizeofBox(&xBox, &yBox, &zBox);
/**
	boxes = new double**[xBox];
        for (int i = 0; i < xBox; i++) {
                boxes[i] = new double*[yBox];
                for (int j = 0; j < yBox; j++) {
                        boxes[i][j] = new double[zBox];
                        for (int k = 0; k < zBox; k++) {
                                boxes[i][j][k] = 0;
                        }
                }
        }
*/

	int boxes[xBox][yBox][zBox];
        for (int i = 0; i < xBox; i++) {
                for (int j = 0; j < yBox; j++) {
                        for (int k = 0; k < zBox; k++) {
				boxes[i][j][k] = 0;
                                printf("Box %d %d %d: %d\n", i, j, k, boxes[i][j][k]);
                        }
                }
        }


	//radius to check
	double distance[num_part];
	for (int i = 0; i < num_part; i++) {
		double radius;

		//assume all points are at least 25 away from surface for now!! BAD CHOICE FOR LATER
		radius = defaultR;
		printf("%d: %lf\n", i,  radius);

		distance[i] = radius;
	}
       
	gettimeofday( &tp, NULL );
	
	double time_2 = tp.tv_usec * (1e-6) + tp.tv_sec;

                for (int a = 0; a < num_part; a++) {
                        double half = L/2;
                        int xB = (int) ((data[a].x + half) * xBox)/L;
                        int yB = (int) ((data[a].y + half) * yBox)/L;
                        int zB = (int) ((data[a].z + half) * zBox)/L;
                        boxes[xB][yB][zB] = boxes[xB][yB][zB] + 1;
                }

/**
        file6.open("begin.txt");
        for (int i = 0; i < xBox; i++) {
                for (int j = 0; j < yBox; j++) {
                        double average2D = 0;
                        for (int k = 0; k < zBox; k++) {
                                double average = (double)boxes[i][j][k];
                                average2D+=average;
                        }
                        file6 << i << " " << j << " " << average2D << "\n";

                }
        }
        file6.close();
*/


	printf("Inside/out initialization time: %le seconds.\n", time_2-time_1 );

	for (int i = 0; i < 10000; i++) {
		printf("%d\n", i);
        	t = t + del_t;
		temp_data = data;
//		printf("before particle move");
        	data = func.move_particles(data, num_part, sigma, eps, D, del_t, kbT, L, file3);
//		printf("after particle move");
   		func.pbc(data, num_part, L);


	//for loop to check whether our points collide with the surface
	// all points inside box
	bool allIn = true;
		for( int j = 0; j < num_part; j++ )
		{
		    double changeInDistance = func.calc_distance(data[j], orig_data[j]);
		    if (changeInDistance >= distance[j]) {

                        int f;
                        double u, v, radius;
			radius = theSurface1->returnRadius(conversion(data[j]), &f, &u, &v,  M5,  M6,  changeInDistance, distance[j],  L, defaultR);
			printf("before anything %lf %lf \n", u, v);
		//	printf("%lf   ", radius);
			
			if (radius > 0) {
				distance[j] = radius;
                                orig_data[j].x = data[j].x;
                                orig_data[j].y = data[j].y;
                                orig_data[j].z = data[j].z;
			} else {
	//			if ((u == 0 && v == 0)) {
	//				data[j].x = temp_data[j].x;
	//				data[j].y = temp_data[j].y;
	//				data[j].z = temp_data[j].z;
	//			} else {
//				int prob = rand() % 10; 
//				if (prob >=7 ) {
					while (prob >= 4) {	
						printf("interesting \n");
						printf("tricky part %lf %lf\n", u, v);
        	                                allIn = false;
                	                        double newPoint[3];
                        	                theSurface1->moveParticleonSurface(&f, &u, &v, newPoint);

                                	        printf("%lf %lf %lf %lf %lf %lf\n", data[j].x, data[j].y, data[j].z, newPoint[0], newPoint[1], newPoint[2]);
                                        	printf("after %lf %lf\n", u, v);

	                                        data[j].x = newPoint[0];
        	                                data[j].y = newPoint[1];
                	                        data[j].z = newPoint[2];
						write_data(file);
						prob = rand() % 10;
					}
				} else {
			        	data[j].x = temp_data[j].x;
                                	data[j].y = temp_data[j].y;
                                	data[j].z = temp_data[j].z;	
//				}
	//			}
			}
		    } 
			
		} 
	

		if ((i % 1 == 0)) {
			write_data(file);
           		//file2 << t << " " << func.calc_potential(data, num_part, sigma, eps, L) << "\n" ;
			
		}

		//printf("\n");
/**
		for (int a = 0; a < num_part; a++) {
		        double half = L/2;
        		int xB = (int) ((data[a].x + half) * xBox)/L;
		        int yB = (int) ((data[a].y + half) * yBox)/L;
        		int zB = (int) ((data[a].z + half) * zBox)/L;
        		boxes[xB][yB][zB] = boxes[xB][yB][zB] + 1;
		}
*/

	}
	
	
        gettimeofday( &tp, NULL );
	double time_3 = tp.tv_usec * (1e-6) + tp.tv_sec;

	printf("Brownian dynamics time: %le\n", time_3-time_2 );
/**
        file4.open("average.txt");
        file5.open("average3D.txt");
	for (int i = 0; i < xBox; i++) {
                for (int j = 0; j < yBox; j++) {
			double average2D = 0;
                        for (int k = 0; k < zBox; k++) {
                                double average = (double)boxes[i][j][k]/1000.0;
                                average2D+=average;
				file5 << i << " " << j << " " << k << " " << average << "\n";
                        }
			file4 << i << " " << j << " " << average2D << "\n";

                }
		file5 << "\n";
        }
*/


	file.close();
	file2.close();
	file3.close();
//	file4.close();
//	file5.close();
	fclose(limit1);
}





