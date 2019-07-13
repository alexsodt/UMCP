#ifndef __p_p_h__
#define __p_p_h__

#include "interp.h"
#include "pcomplex.h"

/* potential */
double PP_V( Simulation *theSimulation );
/* gradient */
double PP_G( Simulation *theSimulation );
double Boxed_PP_G( Simulation *theSimulation );
double Boxed_PP_V( Simulation *theSimulation );
double timePrecedingElasticCollision( Simulation *theSimulation, pcomplex **allComplexes, int ncomplex, double time_step, int *, int *, int *, int* );
double handleElasticCollisions( Simulation *theSimulation, pcomplex **allComplexes, int ncomplex, double time_step );
double nElasticCollisions( Simulation *theSimulation, pcomplex **allComplexes, int ncomplex );

#endif
