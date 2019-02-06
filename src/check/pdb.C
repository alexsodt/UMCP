// by alex sodt
#include <stdio.h>
#include "pdb.h"

void nospacecopy( char *str1, char *str2 )
{
	int tx = 0;
	for( int x = 0; x < strlen(str2); x++ )
	{
		if( str2[x] != ' ' )
		{
			str1[tx] = str2[x];
			tx++;
		}
	}
	str1[tx] = '\0';
}

void readATOM( char *str, struct atom_rec *atrec )
{
	char buffer[256];

	int i, j;
	
	j = 0;
	i = 0;

	for( i = 6; i < 6 + 5; i++,j++ )
		buffer[j] = str[i];
	buffer[11-6] = '\0';
	
	int atN = atoi(buffer);
	atrec->bead = atN;

	j = 0;
	for( i = 12; i < 16; i++,j++)		
		buffer[j] = str[i];
	buffer[16-12]='\0';

	atrec->altloc = str[16];

	atrec->atname = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );

	nospacecopy(atrec->atname, buffer );

	j = 0;

	for( i = 17; i < 21; i++,j++ )
		buffer[j] = str[i];
	buffer[21-17] = '\0';	

	atrec->resname = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );
	nospacecopy(atrec->resname, buffer );

	atrec->chain = str[21];

	j = 0;
	for( i = 22; i < 27; i++, j++ )
		buffer[j] = str[i];
	buffer[27-22] = '\0';

	sscanf(buffer, "%d", & (atrec->res) );

	j = 0;
	
	for( i = 30; i < 38; i++,j++ )
		buffer[j] = str[i];
	buffer[38-30] = '\0';
	
	atrec->x = atof(buffer);

	j = 0;
	
	for( i = 38; i < 46; i++,j++ )
		buffer[j] = str[i];
	buffer[46-38] = '\0';
	
	atrec->y = atof(buffer);
	j = 0;
	
	for( i = 46; i < 54; i++,j++ )
		buffer[j] = str[i];
	buffer[54-46] = '\0';
	
	atrec->z = atof(buffer);

	if( strlen(str) > 60 )
	{	
		j = 0;

		int tlen = strlen(str);
		for( i = 60; i < tlen; i++,j++ )
			buffer[j] = str[i];
		buffer[tlen-60] = '\0';
	
		atrec->vdw = atof(buffer);
	}
	else
		atrec->vdw = 0.0;

	int seg_len = strlen(str) - 72;

	if( seg_len > 0 )
	{
		atrec->segid = (char *)malloc( sizeof(char) * (seg_len+1) );
		nospacecopy( atrec->segid, str+72 );
	}
	else if( atrec->chain == 'A' )
	{
		atrec->segid = (char *)malloc( sizeof(char) * (5) );
		strcpy( atrec->segid, "PROA" );
	}
	else if( atrec->chain == 'B' )
	{
		atrec->segid = (char *)malloc( sizeof(char) * (5) );
		strcpy( atrec->segid, "PROB" );
	}
	else
	{
		atrec->segid = (char *)malloc( sizeof(char) * (1+strlen("UNSET") ) );
		sprintf(atrec->segid, "%s", "UNSET" );
//		atrec->segid = NULL;
	}
/*
	if( strlen(str) >= 73 && atrec->chain == ' ' )
		atrec->chain = str[72];
	if( strlen(str) >= 76 )
		atrec->chain = str[75];
*/
}
/*
void readATOM( char *str, struct atom_rec *atrec )
{
	char buffer[256];

	int i, j;
	
	j = 0;
	i = 0;

	for( i = 6; i < 6 + 5; i++,j++ )
		buffer[j] = str[i];
	buffer[11-6] = '\0';
	
	int atN = atoi(buffer);
	atrec->bead = atN;

	j = 0;
	for( i = 12; i < 16; i++,j++)		
		buffer[j] = str[i];
	buffer[16-12]='\0';

	atrec->altloc = str[16];

	atrec->atname = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );
	strcpy(atrec->atname, buffer );

	j = 0;

	for( i = 17; i < 20; i++,j++ )
		buffer[j] = str[i];
	buffer[20-17] = '\0';	

	atrec->resname = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );
	strcpy(atrec->resname, buffer );

	atrec->chain = str[21];

	j = 0;
	for( i = 22; i < 27; i++, j++ )
		buffer[j] = str[i];
	buffer[27-22] = '\0';

	sscanf(buffer, "%d", & (atrec->res) );

	j = 0;
	
	for( i = 30; i < 38; i++,j++ )
		buffer[j] = str[i];
	buffer[38-30] = '\0';
	
	atrec->x = atof(buffer);

	j = 0;
	
	for( i = 38; i < 46; i++,j++ )
		buffer[j] = str[i];
	buffer[46-38] = '\0';
	
	atrec->y = atof(buffer);
	j = 0;
	
	for( i = 46; i < 54; i++,j++ )
		buffer[j] = str[i];
	buffer[54-46] = '\0';
	
	atrec->z = atof(buffer);

	double d1, d2;
	int nr = sscanf( str + 54, " %lf %lf", &d1, &d2 );
	
	if( nr == 2 )
		atrec->aux = d2;
	else
		atrec->aux = 0;
}
*/
void printATOM( FILE *toFile, int bead, int res, struct atom_rec *atrec, double aux, int write_hex)
{
	fprintf(toFile,"ATOM  ");
	char tbuf[256];

	if( !write_hex )
		sprintf(tbuf, "%d", bead );
	else
		sprintf(tbuf, "%x", bead );

	if( strlen(tbuf) < 5 )
	for( int x = 0; x < 5 - strlen(tbuf); x++ )
		fprintf(toFile," ");
	fprintf(toFile,"%s",tbuf);
	fprintf(toFile," " );
	if( strlen(atrec->atname) < 4 )
		sprintf(tbuf, " %s", atrec->atname );
	else
		sprintf(tbuf, "%s", atrec->atname );
	fprintf(toFile,"%s",tbuf);
	if( strlen(tbuf) < 4 )
	for( int x = 0; x < 4 - strlen(tbuf); x++ )
		fprintf(toFile," ");
	fprintf(toFile, "%c", atrec->altloc );
	sprintf(tbuf, "%s", atrec->resname );
	fprintf(toFile,"%s",tbuf);
	if( strlen(tbuf) < 4 )
	for( int x = 0; x < 4 - strlen(tbuf); x++ )
		fprintf(toFile," " );
//	fprintf(toFile, " " );
	fprintf(toFile,"%c", atrec->chain );
	if( write_hex )
		sprintf(tbuf, "%x", res );
	else
		sprintf(tbuf, "%d", res );
	tbuf[4] = '\0';
	if( strlen(tbuf) < 4 )
	for( int x = 0; x < 4 - strlen(tbuf); x++ )
		fprintf(toFile," ");
	fprintf(toFile,"%s",tbuf);
	fprintf(toFile,"    ");
	sprintf(tbuf, "%.3lf", atrec->x );
	tbuf[8] = '\0';
	if( strlen(tbuf) < 8 )
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
		fprintf(toFile," ");
	fprintf(toFile,"%s",tbuf);
	sprintf(tbuf, "%.3lf", atrec->y );
	tbuf[8] = '\0';
	if( strlen(tbuf) < 8 )
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
		fprintf(toFile," ");
	fprintf(toFile,"%s",tbuf);
	sprintf(tbuf, "%.3lf", atrec->z );
	tbuf[8] = '\0';
	if( strlen(tbuf) < 8 )
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
		fprintf(toFile," ");
	fprintf(toFile,"%s",tbuf);
	//fprintf(toFile,"  1.00  %lf\n", aux);
	if( aux < 100 && aux > -100 )
		fprintf(toFile,"  1.00  %1.2lf", aux);
	else
		fprintf(toFile,"  1.00  0.00");
	if( atrec->segid )
	{
		fprintf(toFile, "      %s", atrec->segid ); 
	}	
	fprintf(toFile,"\n");
}

void printATOM( char *toFile_in, int bead, int res, struct atom_rec *atrec, double aux, int print_hex)
{
	char *toFile = toFile_in;
	sprintf(toFile,"ATOM  ");
	toFile += strlen(toFile);
	char tbuf[256];
	sprintf(tbuf, "%d", bead );

	if( strlen(tbuf) < 5 )
	for( int x = 0; x < 5 - strlen(tbuf); x++ )
	{
		sprintf(toFile," ");
		toFile += strlen(toFile);
	}
	sprintf(toFile,"%s",tbuf);
	toFile += strlen(toFile);
	
	sprintf(toFile," " );
	toFile += strlen(toFile);
	sprintf(tbuf, "%s", atrec->atname );
	sprintf(toFile,"%s",tbuf);

	if( strlen(tbuf) < 4 )
	for( int x = 0; x < 4 - strlen(tbuf); x++ )
	{
		sprintf(toFile," ");
		toFile += strlen(toFile);
	}

	toFile += strlen(toFile);
	sprintf(toFile, "%c", atrec->altloc );
	toFile += strlen(toFile);
	sprintf(tbuf, "%s", atrec->resname );
	if( strlen(tbuf) < 3 )
	for( int x = 0; x < 3 - strlen(tbuf); x++ )
	{
		sprintf(toFile," " );
		toFile += strlen(toFile);
	}
	sprintf(toFile,"%s",tbuf);
	toFile += strlen(toFile);
//	fprintf(toFile, " " );
	sprintf(toFile,"%c", atrec->chain );
	toFile += strlen(toFile);
	sprintf(tbuf, "%d", res );
	tbuf[4] = '\0';
	if( strlen(tbuf) < 4 )
	for( int x = 0; x < 4 - strlen(tbuf); x++ )
	{
		sprintf(toFile," ");
		toFile += strlen(toFile);
	}
	sprintf(toFile,"%s",tbuf);
	toFile += strlen(toFile);
	sprintf(toFile,"    ");
	toFile += strlen(toFile);
	sprintf(tbuf, "%.3lf", atrec->x );
	tbuf[8] = '\0';
	if( strlen(tbuf) < 8 )
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
	{
		sprintf(toFile," ");
		toFile += strlen(toFile);
	}
	sprintf(toFile,"%s",tbuf);
	toFile += strlen(toFile);
	sprintf(tbuf, "%.3lf", atrec->y );
	tbuf[8] = '\0';
	if( strlen(tbuf) < 8 )
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
	{
		sprintf(toFile," ");
		toFile += strlen(toFile);
	}
	sprintf(toFile,"%s",tbuf);
	toFile += strlen(toFile);
	sprintf(tbuf, "%.3lf", atrec->z );
	tbuf[8] = '\0';
	if( strlen(tbuf) < 8 )
	for( int x = 0; x < 8 - strlen(tbuf); x++ )
	{
		sprintf(toFile," ");
		toFile += strlen(toFile);
	}
	sprintf(toFile,"%s",tbuf);
	toFile += strlen(toFile);
	sprintf(toFile,"  1.00  0.00");
	if( atrec->segid )
	{
		sprintf(toFile, "      %s", atrec->segid ); 
	}	
	sprintf(toFile, "\n");
	toFile += strlen(toFile);
}



void printCRYST( FILE *toFile, double LX, double LY, double LZ, double alpha, double beta, double gamma )
{
        fprintf(toFile,"CRYST1");
        char tbuf[256];
        sprintf(tbuf, "%.3lf", LX);
        if( strlen(tbuf) < 9 ) 
        for( int x = 0; x < 9 - strlen(tbuf); x++ )
                fprintf(toFile," ");
        fprintf(toFile,"%s",tbuf);
    
        sprintf(tbuf, "%.3lf", LY);
        if( strlen(tbuf) < 9 ) 
        for( int x = 0; x < 9 - strlen(tbuf); x++ )
                fprintf(toFile," ");
        fprintf(toFile,"%s",tbuf);
    
        sprintf(tbuf, "%.3lf", LZ);
        if( strlen(tbuf) < 9 ) 
        for( int x = 0; x < 9 - strlen(tbuf); x++ )
                fprintf(toFile," ");
        fprintf(toFile,"%s",tbuf);


        sprintf(tbuf, "%.2lf",alpha);
        if( strlen(tbuf) < 7 ) 
        for( int x = 0; x < 7 - strlen(tbuf); x++ )
                fprintf(toFile," ");
        fprintf(toFile,"%s",tbuf);
        sprintf(tbuf, "%.2lf",beta);
        if( strlen(tbuf) < 7 ) 
        for( int x = 0; x < 7 - strlen(tbuf); x++ )
                fprintf(toFile," ");
        fprintf(toFile,"%s",tbuf);
        sprintf(tbuf, "%.2lf",gamma);
        if( strlen(tbuf) < 7 ) 
        for( int x = 0; x < 7 - strlen(tbuf); x++ )
                fprintf(toFile," ");
        fprintf(toFile,"%s",tbuf);
	fprintf(toFile,"\n");

}

