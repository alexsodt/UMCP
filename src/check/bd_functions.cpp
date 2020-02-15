#include "bd_functions.h"
#include "Point.h"
#include <cmath>
#include <iostream>
#include <random>
#include <time.h>


using namespace std;

BD_Functions::BD_Functions() {

}

void BD_Functions::pbc(Point*& data, double num_part, double box_length) {
	double half = box_length/2;
	for (int i = 0; i < num_part; i++) {

		if (data[i].x > half) {
			data[i].x = data[i].x - box_length;
		} else if (data[i].x < -half) {
			data[i].x = data[i].x + box_length;
		}
		
		if (data[i].y > half) {
			data[i].y = data[i].y - box_length;
		} else if (data[i].y < -half) {
			data[i].y = data[i].y + box_length;
		}

		if (data[i].z > half) {
			data[i].z = data[i].z - box_length;
		} else if (data[i].z < -half) {
			data[i].z = data[i].z + box_length;
		}


	}	
}

Point BD_Functions::min_image(Point a, Point b, double box_length) {
	double half = box_length/2;
	double dx = b.x - a.x;
	double dy = b.y - a.y;
	double dz = b.z - a.z;
	if (dx > half) {
        b.x = b.x - box_length;
	} else if (dx < -half) {	
		b.x = b.x + box_length;
	}
	
	if (dy > half) {
	    b.y = b.y - box_length;
	} else if (dy < -half) {
		b.y = b.y + box_length;
	}

	if (dz > half) {
	    b.z = b.z - box_length;
	} else if (dz < -half) {
		b.z = b.z + box_length;
	}

	return b;
}

double BD_Functions::calc_potential(Point*& data, int num_part, double sigma, double epsilon, double box_length) {
	double potential = 0.0;
	for (int i = 0; i < num_part; i++) {
		for (int j = i + 1; j < num_part; j++) {
			Point b = min_image(data[i], data[j], box_length);
//			double r = calc_distance(data[i], data[j]);
			double r = calc_distance(data[i], b);
			potential = potential + 4*epsilon*(pow(sigma/r, 12) - pow(sigma/r, 6));
		}
	}
	return potential;
}  

double BD_Functions::calc_harm_potential(Point*& data, int num_part, double k, double r0, double box_length) {
	double potential = 0.0;
		for (int i = 0; i < num_part; i++) {
        for (int j = i + 1; j < num_part; j++) {
            Point b = min_image(data[i], data[j], box_length);
            double r = calc_distance(data[i], b);
            potential = potential + .5 * k * pow(r - r0, 2);
        }
    }
    return potential;	
}

double BD_Functions::calc_distance(Point a, Point b) {
	double r = sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2) + pow((a.z - b.z), 2));
	return r;
}



Point* BD_Functions::calc_force(Point*& data, int num_part, double sigma, double epsilon, double box_length) {
	Point * forces = new Point[num_part];

        for( int i = 0; i < num_part; i++ )
	{
		forces[i].x = 0;
		forces[i].y = 0;
		forces[i].z = 0;
	}

	for (int i = 0; i < num_part; i++) {
		for (int j = i + 1; j < num_part; j++) {
//			double r = calc_distance(data[i], data[j]);
//			double dx = (data[i].x - data[j].x);
//			double dy = (data[i].y - data[j].y);
	
			Point b = min_image(data[i], data[j], box_length);
			double dx = (data[i].x - b.x);
			double dy = (data[i].y - b.y);
			double dz = (data[i].z - b.z);
			double r = calc_distance(data[i], b);
			double dRdx = dx/r;
			double dRdy = dy/r;
			double dRdz = dz/r;
			double dV = (-4*epsilon*(12*pow(sigma/r, 12) - 6*pow(sigma/r, 6)))/r;
			double Fx = -dV * dRdx;
			double Fy = -dV * dRdy;
			double Fz = -dV * dRdz;
			forces[i].x += Fx;
			forces[i].y += Fy;
			forces[i].z += Fz;
			forces[j].x += -Fx;
			forces[j].y += -Fy;
			forces[j].z += -Fz;

		}	
//	std::cout << forces[i].x << " " << forces[i].y << std::endl;
	}
	return forces;
}


 
Point* BD_Functions::calc_harm_force(Point*& data, int num_part, double k, double r0, double box_length) {
    Point * forces = new Point[num_part];
        for( int i = 0; i < num_part; i++ )
	{
		forces[i].x = 0;
		forces[i].y = 0;
		forces[i].z = 0;
	}
    for (int i = 0; i < num_part; i++) {
        for (int j = i + 1; j < num_part; j++) {
            Point b = min_image(data[i], data[j], box_length);
            double dx = (b.x - data[i].x);
            double dy = (b.y - data[i].y);
            double dz = (b.z - data[i].z);
            double r = calc_distance(data[i], b); 
            double dRdx = dx/r;
            double dRdy = dy/r;
            double dRdz = dy/r;
            forces[i].x += k * (r - r0) * dRdx;
            forces[i].y += k * (r - r0) * dRdy;
            forces[i].z += k * (r - r0) * dRdz;
            forces[j].x += -k * (r - r0) * dRdx;
            forces[j].y += -k * (r - r0) * dRdy;
            forces[j].z += -k * (r - r0) * dRdz;

        }   
//  std::cout << forces[i].x << " " << forces[i].y << std::endl;
    }   
    return forces;
}




Point* BD_Functions::move_particles(Point*& data, int num_part, double sigma, double epsilon, double diffusion, double del_t, double kbT, double box_length, ofstream& file) {
	//std::random_device rdev;
	//std::default_random_engine generator {rdev()};
	std::random_device rdev;
	std::mt19937 generator(rdev());
	//std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 10);
	Point* new_data = new Point[num_part];
	//Point* forces = calc_force(data, num_part, sigma, epsilon, box_length); 	
//	Point* forces = calc_harm_force(data, num_part, sigma, epsilon, box_length); 	

	for (int i = 0; i < num_part; i++) {
		double rand1 = distribution(generator);
		double rand2 = distribution(generator);
		double rand3 = distribution(generator);
		//double dx = del_t * (kbT/diffusion) * forces[i].x + sqrt(2*diffusion*del_t) * rand1; 
		//double dy = del_t * (kbT/diffusion) * forces[i].y + sqrt(2*diffusion*del_t) * rand2;
		//double dz = del_t * (kbT/diffusion) * forces[i].z + sqrt(2*diffusion*del_t) * rand3;
		double dx = sqrt(2*diffusion*del_t) * rand1; 
		double dy = sqrt(2*diffusion*del_t) * rand2;
		double dz = sqrt(2*diffusion*del_t) * rand3;
//		dx = del_t * (kbT/diffusion) * forces[i].x; 
//		dy = del_t * (kbT/diffusion) * forces[i].y;
		new_data[i].x = data[i].x + dx;
		new_data[i].y = data[i].y + dy;
		new_data[i].z = data[i].z + dz;
//		file << forces[i].x << " " << forces[i].y << " " << forces[i].z << " " << sqrt(2*diffusion*del_t) * rand1 << " " << sqrt(2*diffusion*del_t) * rand2 << 
//			sqrt(2*diffusion*del_t) * rand3; 
	}

//	double mean = 0; 
//	double* rs = new double[num_part];
//	for (int i = 0; i < num_part; i++) {
//		for (int j = i+1; j < num_part; j++) {
//			file << calc_distance(data[i], data[j]);
//			mean += calc_distance(data[i], data[j]);
//		}
//	}
//	file << "\n";

//	mean = mean/num_part;
//	double stdev = 0;
//	for (int i = 0; i < num_part; i++) {
//		stdev += pow(rs[i] - mean, 2); 
//	}
//	stdev = sqrt(stdev/(num_part-1));
//	std::cout << mean << " " << stdev << std::endl;

/**
 //just to calc mean and stdev of random nums
	double mean;
	double stdev;
	double* num = new double[100];
	for (int i = 0; i < 100; i++) {
		num[i] = distribution(generator);
		mean += num[i];
	}
	mean = mean/100;
	for (int i = 0; i < 100; i++) {
		stdev += pow(num[i]-mean,2);
	}
	stdev = sqrt(stdev/99);
	std::cout << "here" << mean << " " << stdev << std::endl;
*/

	return new_data;

}
