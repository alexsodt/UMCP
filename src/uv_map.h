#ifndef __uvmaph__
#define __uvmaph__
void disable_random_uv_step( void );
void enable_random_uv_step( void );
void randomDirection( surface *theSurface, double *r, int f, double u, double v, double len, double duv[2] );
#endif
