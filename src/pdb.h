#ifndef __pdbh__
#define __pdbh__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct atom_rec
{
	int bead;
	int res;
	char *atname;
	char *resname;
	char *segid;
	char altloc;
	char chain;
	double x;
	double y;
	double z;

	double aux;
	double vdw;
	double charge;
	int segRes;
	void zap( void ) { free( atname ); free(resname); if( segid ) free( segid ); }
};

void readATOM( char *str, struct atom_rec *atrec );

void printATOM( FILE *toFile, int bead, int res, struct atom_rec *atrec, double aux=0, int write_hex = 0 );
void printATOM( char *toFile, int bead, int res, struct atom_rec *atrec, double aux=0, int write_hex = 0 );
void printCRYST( FILE *toFile, double LX, double LY, double LZ, double alpha, double beta, double gamma );

#endif
