#include <stdio.h>
#include <stdlib.h>

void printCompilationTime( void )
{
	printf("Built on %s at %s.\n", __DATE__, __TIME__ );
}
