#include "simulation.h"
#include "pcomplex.h"
#include "mutil.h"
#include "parallel.h"

void Simulation::wrapPBC( double *dr1, double *alphas)
{
	double put[3];
	MinImage3D( dr1, PBC_vec, put, alphas );
	
}

int Simulation::AddComplex( pcomplex *addMe )
{
	int free_spot = -1;
	int added = -1;

	for( int p = 0; p < ncomplex; p++ )
	{
		int nwatchers = allComplexes[p]->nwatchers;
#ifdef PARALLEL
		ParallelSum(&nwatchers,1);
#endif
		if( allComplexes[p]->disabled && nwatchers == 0 )
			free_spot = p;
	}

	if( free_spot )
	{
		added = free_spot;

		allComplexes[free_spot]->destroy();
		free(allComplexes[free_spot]);
		allComplexes[free_spot] = addMe;
	}
	else
	{
		if( ncomplexSpace == ncomplex )
		{
			ncomplexSpace += 10;
	
			allComplexes = (pcomplex **)realloc( allComplexes, sizeof(pcomplex *) * ncomplexSpace );
	
		}
		
		added = ncomplex;
			
		allComplexes[ncomplex] = addMe;
		ncomplex++;
	
		// assign this to a parallel process.
		int task=0;
	
#ifdef PARALLEL		
		// root process decides
		if( par_info.my_id == BASE_TASK )
			task = ncomplex % par_info.nprocs;
		ParallelBroadcast(&task,1);
#endif
	
		if( task == par_info.my_id )
		{
			par_info.nc += 1;
			par_info.complexes = (int *)realloc( par_info.complexes, par_info.nc );
			par_info.complexes[par_info.nc] = added;
		}
	}
	return added;
}

void Simulation::RemoveComplexDelayed( int id )
{
	allComplexes[id]->disable();
}

void Simulation::GarbageCollection( void )
{
}

