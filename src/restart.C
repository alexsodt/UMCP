#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "interp.h"
#include "pcomplex.h"
#include "util.h"
void Simulation::saveRestart( FILE *theFile, int seed ) 
{
	char *buf;

	saveRestart( &buf, seed );

	fprintf(theFile, "%s", buf );

	free(buf);
}

void Simulation::saveRestart( char **buf, int seed)
{
	int buffer_size = 4096;
	int put = 0;


	char *theBuffer = (char *)malloc( sizeof(char) * buffer_size );
	char *tbuf = (char *)malloc( sizeof(char) * 4096 );
	
	sprintf(tbuf, "%lf %lf %lf\n", alpha[0], alpha[1], alpha[2] );
		
	while( put + strlen(tbuf) >= buffer_size )
	{
		buffer_size *= 2;
		buffer_size += strlen(tbuf);
	
		theBuffer = (char *)realloc( theBuffer, buffer_size ); 
	}
	
	strcpy( theBuffer+put, tbuf );
	put += strlen(tbuf);
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		surface *theSurface = sRec->theSurface;
		double *rsurf = sRec->r;
		double *pp = sRec->pp;
		int nv = theSurface->nv;

		for( int v = 0; v < nv+1; v++ )
		{
			if( sRec->NQ == 0 )
				sprintf( tbuf, "%lf %lf %lf %lf %lf %lf\n", rsurf[3*v+0], rsurf[3*v+1], rsurf[3*v+2], pp[3*v+0], pp[3*v+1], pp[3*v+2] );
			else
				sprintf( tbuf, "%lf %lf %lf\n", rsurf[3*v+0], rsurf[3*v+1], rsurf[3*v+2] );
	
			if( put + strlen(tbuf) >= buffer_size )
			{
				buffer_size *= 2;
				buffer_size += strlen(tbuf);
	
				theBuffer = (char *)realloc( theBuffer, buffer_size ); 
			}
	
			strcpy( theBuffer+put, tbuf );
	
			put += strlen(tbuf);
		}
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

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		if( sRec->do_gen_q )
		{
			for( int Q = 0; Q < sRec->NQ; Q++ )
			{
				sprintf( tbuf, "GQ %.14le\n", sRec->pp[Q] );
		
				if( put + strlen(tbuf) >= buffer_size )
				{
					buffer_size *= 2;
					buffer_size += strlen(tbuf);
		
					theBuffer = (char *)realloc( theBuffer, buffer_size ); 
				}
		
				strcpy( theBuffer+put, tbuf );
				put += strlen(tbuf);
			}
		}
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

void Simulation::loadRestart( FILE *loadFile, int *seed )
{
	int buffer_size = 100000;
	int put = 0;

	char *buffer = (char *)malloc( sizeof(char) * buffer_size );

	getLine( loadFile, buffer );

	int nr = sscanf( buffer, "%lf %lf %lf", alpha+0, alpha+1, alpha+2 );

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		surface *theSurface = sRec->theSurface;
		double *rsurf = sRec->r;
		double *pp = sRec->pp;

		for( int v = 0; v < theSurface->nv+1; v++ )
		{
			getLine( loadFile, buffer );
			int nr = 0;
			if( sRec->NQ == 0 && pp )
				nr = sscanf( buffer, "%lf %lf %lf %lf %lf %lf\n", rsurf+3*v+0, rsurf+3*v+1, rsurf+3*v+2, pp+3*v+0, pp+3*v+1, pp+3*v+2 );
			else
				nr = sscanf( buffer, "%lf %lf %lf\n", rsurf+3*v+0, rsurf+3*v+1, rsurf+3*v+2 );
		}
	}
		
	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->loadComplex(loadFile,this,0);

	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		for( int Q = 0; Q < sRec->NQ; Q++ )
		{
			getLine(loadFile, buffer );
			sprintf( buffer, "GQ %.14le\n", sRec->pp[Q] );
		}
	}	
	nr = sscanf( buffer, "seed %d",seed );
	
	free(buffer);
}
