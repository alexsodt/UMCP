#ifndef __sansh__
#define __sansh__

#define SANS_SAMPLE_NRM		0
#define SANS_SAMPLE_LAT		1

void loadBetaZ( const char *fileName, int load_to_spline );
void initSANS( parameterBlock *block, double **qvals, int *nq );
void writeSq( char *fileName, double *B_hist, double A2dz2_sampled, double sans_max_r, int nsans_bins, double qmin, double qmax, int nq, double *qvals=NULL );

#endif
