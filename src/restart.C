#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "interp.h"
#include "pcomplex.h"
void surface::saveRestart( FILE *theFile, double *rsurf, double *pp, pcomplex **allComplexes, int ncomplex, int NQ, int seed )
{
	char *buf;

	saveRestart( &buf, rsurf, pp, allComplexes, ncomplex, NQ, seed );

	fprintf(theFile, "%s", buf );

	free(buf);
}

void surface::saveRestart( char **buf, double *rsurf, double *pp, pcomplex **allComplexes, int ncomplex, int NQ, int seed)
{
	int buffer_size = 4096;
	int put = 0;

	char *theBuffer = (char *)malloc( sizeof(char) * buffer_size );
	char *tbuf = (char *)malloc( sizeof(char) * 4096 );
	for( int v = 0; v < nv+1; v++ )
	{
		if( NQ == 0 )
			sprintf( tbuf, "%lf %lf %lf %lf %lf %lf\n", rsurf[3*v+0], rsurf[3*v+1], rsurf[3*v+2], pp[3*v+0], pp[3*v+1], pp[3*v+2] );
		else
			sprintf( tbuf, "%lf %lf %lf %lf %lf %lf\n", rsurf[3*v+0], rsurf[3*v+1], rsurf[3*v+2] );

		if( put + strlen(tbuf) >= buffer_size )
		{
			buffer_size *= 2;
			buffer_size += strlen(tbuf);

			theBuffer = (char *)realloc( theBuffer, buffer_size ); 
		}

		strcpy( theBuffer+put, tbuf );

		put += strlen(tbuf);
	}

	for( int c = 0; c < ncomplex; c++ )
	{
		// save at least 4096 characters for each complex.

		int len = 0;
		int tries = 0;

		while( tries < 3 && allComplexes[c]->saveComplex( theBuffer + put, &len, buffer_size - put ))
 		{
			buffer_size *= 2;

			theBuffer = (char *)realloc( theBuffer, buffer_size ); 

			tries++;
		}

		put += len;
	}

	for( int Q = 0; Q < NQ; Q++ )
	{
		sprintf( tbuf, "GQ %.14le\n", pp[Q] );

		if( put + strlen(tbuf) >= buffer_size )
		{
			buffer_size *= 2;
			buffer_size += strlen(tbuf);

			theBuffer = (char *)realloc( theBuffer, buffer_size ); 
		}

		strcpy( theBuffer+put, tbuf );
		put += strlen(tbuf);
	}
		
	sprintf( tbuf, "seed %d\n", seed );

	if( put + strlen(tbuf) >= buffer_size )
	{
		buffer_size *= 2;
		buffer_size += strlen(tbuf);

		theBuffer = (char *)realloc( theBuffer, buffer_size ); 
	}

	strcpy( theBuffer+put, tbuf );
	put += strlen(tbuf);

	free(tbuf);

	*buf = theBuffer;
}
