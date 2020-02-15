#ifndef __lagrangeh__
#define __lagrangeh__

struct force_set
{	
	int npts;
	int *face_set;
	double *uv_set;
	double *frc_inverse;
	double *frc_coef;
	double *mass;
	int    *frc_coef_list;
	int    *frc_ncoef;
};

void clearForceSet( force_set *theForceSet );

#endif
