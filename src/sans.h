#ifndef __sansh__
#define __sansh__

#define SANS_SAMPLE_NRM		0
#define SANS_SAMPLE_LAT		1

void writeSq( char *fileName, double *B_hist, double A2dz2_sampled, double sans_max_r, int nsans_bins, double qmin, double qmax, int nq );

#endif
