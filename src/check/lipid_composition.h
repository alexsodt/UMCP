#ifndef __lipidcompositionh__
#define __lipidcompositionh__

struct localLipidComposition
{
	double A0;
	double A_inst;
	double *innerLeaflet;
	double *outerLeaflet;
};

struct bilayerComposition
{
	int nlipidTypes;
	double *c0;
	double *APL;
	char **names;
};

#endif
