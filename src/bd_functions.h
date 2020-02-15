#ifndef BD_FUNCTIONS_H
#define BD_FUNCTIONS_H

#include "Point.h"
#include <fstream>

using namespace std;
class BD_Functions {

	public:
		BD_Functions();

//		double ** input_data(const char*& filename);

		void pbc(Point*& data, double num_part, double box_length);

		Point min_image(Point a, Point b, double box_length);
	
		double calc_potential(Point*& data, int num_part, double sigma, double epsilon, double box_length);
		
		double calc_harm_potential(Point*& data, int num_part, double k, double r0, double box_length);
		
		double calc_distance(Point a, Point b);
		
		Point* calc_force(Point*& data, int num_part, double sigma, double epsilon, double box_length);
		Point* calc_harm_force(Point*& data, int num_part, double k, double r0, double box_length);

		Point* move_particles(Point*& data, int num_part, double sigma, double epsilon, double diffusion, double del_t, double kbT, double box_length, ofstream& file);

};

#endif
