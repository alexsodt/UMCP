// by alex sodt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"
#include <stdarg.h>
int readNDoubles( char *buffer, double *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
	int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%lf", vals+x );

		if( nr != 1 )
			return nread;

		nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
                        s++;
        }

	return nread;
}



void getLine( FILE *theFile, char *theBuffer )
{
        int i = 0;

        while( !feof(theFile) )
        {
                char tc = fgetc(theFile);

                if( tc != '\n' && i < 50000 )
                {
                        theBuffer[i++] = tc;
                }
                else if( tc != '\n' && i >= 50000 )
		{
		}
		else
                        break;
        }

        theBuffer[i] = '\0';
}

int goToField( const char * buffer, int f )
{
	int ret = -1;
	int cf = 0;
	
	int cp = 0;
	while( buffer[cp] && (buffer[cp] == ' ' || buffer[cp] == '\t') )
		cp++;

	if( !buffer[cp] ) return -1;


	while( cf < f )
	{
		while( buffer[cp] && (buffer[cp] != ' ' && buffer[cp] != '\t') )
			cp++;
		
		if( !buffer[cp] ) return -1;

		while( buffer[cp] && (buffer[cp] == ' ' || buffer[cp] == '\t') )
			cp++;

		if( !buffer[cp] ) return -1;

		cf++;
	}
		

	return cp;
}

int readNInts( char *buffer, int *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
        int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {   
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%d", vals+x );

                if( nr != 1 ) 
                        return nread;

                nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
                        s++;
        }   

        return nread;
}

int my_isnan( double a )
{
	if( !(a > 0 || a < 1) )
		return 1;
	return 0;
}


static FILE *default_fp = 0;

void set_default_stream(FILE *fp)
{
    default_fp = fp;
}

int mprintf(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    if (default_fp == 0)
        default_fp = stdout;

    int rv = vfprintf(default_fp, fmt, args);

    va_end(args);
    return(rv);
 }

int decodeString( char *buf, char **out, int nmax )
{
	int done = 0;

	int ns = 0;
	int t = 0;
	int lt = 0;
	int lim = 256;
	char cur[lim+1];

	while( *(buf+t) )
	{
		while( *(buf+t) == '+' ) t++;

		while( *(buf+t) && *(buf+t) != '+' )
		{
			if( lt < lim )
			{
				cur[lt] = *(buf+t);	
				lt++;
			}

			t++;
		} 		

		cur[lt] = '\0';
 
		if( lt > 0 && ns < nmax )
		{
			out[ns] = (char *)malloc( sizeof(char) * (1 + lt ) );
			strcpy( out[ns], cur );
			ns++;
		}
			
		lt = 0;

	}

	return ns;
}

void print5( int val, char *str )
{
        if( val < 10 )
                sprintf(str, "0000%d", val );
        else if( val < 100 )
                sprintf(str, "000%d", val );
        else if( val < 1000 )
                sprintf(str, "00%d", val );
        else if( val < 10000 )
                sprintf(str, "0%d", val );
        else
                sprintf(str, "%d", val );
}
