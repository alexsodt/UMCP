
#ifndef __utilh__
#define __utilh__

#include <stdio.h>

void getLine( FILE *theFile, char *theBuffer );
int readNDoubles( char *buffer, double *vals, int nvalues );
int readNInts( char *buffer, int *vals, int nvalues );
int goToField( const char * buffer, int f );
int my_isnan( double val );

extern void set_default_stream(FILE *fp);
extern int  mprintf(const char *fmt, ...);
int decodeString( char *buf, char **out, int nmax );
void print5( int val, char *str );

#endif
