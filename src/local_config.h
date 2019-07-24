void printCompilationTime(void);

#include "../config.h"

#if defined(HAVE_FFTW3_H) && defined(HAVE_LIBFFTW3)
#define FFTW
#endif

#ifdef PARALLEL
#include "mpi.h"
#endif


//#define GPU
//#define GLOBAL_AREA
//#define MICRO_KA
//#define USE_KA4
//#define M_MIN -1 
//#define M_MAX 1
#define FFTW
#define DISABLE_ALIGNMENT
