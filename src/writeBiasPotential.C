#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

double RMIN = 150.0;
double RMAX = 800.0;
int    NBINS = 100;

int main( int argc, char **argv )
{
	if( argc < 3 )
	{
		printf("Syntax: writeBiasPotential R0 meshFile [fraction prune]\n");
		return 0;
	}

	double R0 = atof(argv[1]);
	FILE *theFile = fopen(argv[2],"r");
	double fraction_prune = 0;
	if( argc > 3 )
		fraction_prune = atof(argv[3]);

	int nvals = 0;
	char buffer[4096];

	while( !feof(theFile) )
	{
		getLine(theFile, buffer );
		if( !feof(theFile) ) nvals++;
	}
	rewind(theFile);


	double hist[NBINS];
	memset( hist, 0, sizeof(double) * NBINS );

	double kT = 0.592;

	int nread = 0;

	while( !feof(theFile) )
	{	
		getLine( theFile, buffer );
		if( feof(theFile) ) break;
		double alpha;
		int nr = sscanf( buffer, "%lf", &alpha );

		if( nr !=  1 )
		{
			printf("Failed to read value from line '%s'.\n", buffer );
			exit(1);
		}

		double R = R0 * alpha;

		if( nread > fraction_prune * nvals )
		{
			int rbin = NBINS*(R-RMIN)/(RMAX-RMIN);
			if( rbin >=0 && rbin < NBINS )
				hist[rbin] += 1;
		}
		nread++;
	}

	for( int b = 0; b < NBINS; b++ )
	{
		if( hist[b] > 0 )
			printf("%lf %lf\n", RMIN + (b+0.5) * (RMAX-RMIN)/NBINS, -kT * log( hist[b] ) ); 
		else
			printf("%lf XX\n", RMIN + (b+0.5) * (RMAX-RMIN)/NBINS ); 
	}
}







