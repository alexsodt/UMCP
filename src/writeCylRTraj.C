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
	if( argc <3  )
	{
		printf("Syntax: writeBiasPotential R0 meshFile [fraction prune=0] [getstats=0]\n");
		return 0;
	}

	double R0 = atof(argv[1]);
	FILE *theFile = fopen(argv[2],"r");
	double fraction_prune =0;
	int get_stats = 0;
	if( argc > 3 )
		fraction_prune = atof(argv[3]);
	if( argc > 4 )
		get_stats = atoi(argv[4]);
	char buffer[4096];
	int nvals = 0;
	while( !feof(theFile) )
	{
		getLine(theFile, buffer );
		if( !feof(theFile) ) nvals++;
	}
	rewind(theFile);

	double hist[NBINS];
	memset( hist, 0, sizeof(double) * NBINS );

	double kT = 0.592;

	double av = 0;
	double av2 = 0;
	double nav = 0;
	int t = 0;
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
	
		if( t > fraction_prune * nvals )
		{
			av += R;
			av2 += R*R;
			nav += 1;
			if( !get_stats )
				printf("%d %lf\n", t, R );
		}
		t++;
	}
	av /= nav;
	av2 /= nav; 
	if( get_stats )
		printf("AV: %lf WIDTH: %lf SAMPLES: %lf\n", av, sqrt(av2-av*av), nav );
}







