#ifndef __lipidcompositionh__
#define __lipidcompositionh__

struct localLipidComposition
{
	double A0;
	double A_inst;
	double *innerLeaflet;
	double *outerLeaflet;
	int    tracer; // identity of a lipid tracer particle to track diffusion.
};

struct bilayerComposition
{
	int nlipidTypes;
	double *sum_C;
	double *num_C;
	double *c0;
	double *APL;
	char **names;
	double *input_innerLeaflet;
	double *input_outerLeaflet;
};

#endif
