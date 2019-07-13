#include "simulation.h"
#include "mutil.h"

void Simulation::wrapPBC( double *dr1, double *alphas)
{
	double put[3];
	MinImage3D( dr1, PBC_vec, put, alphas );
	
}
