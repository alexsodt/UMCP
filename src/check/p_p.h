#ifndef __p_p_h__
#define __p_p_h__

#include "interp.h"
#include "pcomplex.h"

/* potential */
double PP_V( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex );
/* gradient */
double PP_G( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, double *mesh_g);
double Boxed_PP_G( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, double *mesh_grad );
double Boxed_PP_V( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex );
void handleElasticCollisions( surface *theSurface, pcomplex **allComplexes, int ncomplex, double time_step );
double nElasticCollisions( surface *theSurface, pcomplex **allComplexes, int ncomplex );

#endif
