// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "lapack_we_use.h"
	
struct DCDheader1
{
	char length;
	char hdr_string[7];
	int nfile;
	int npriv;
	int nsavc;
	int nstep; // icntrl 4
	int unk1[4];
	int fixed; // icntrl 9
	int unk2[1];
	int qcrys;
	int dim4;
	int qcg;
	int unk3[7]; // icntrl12-20
	int itemp;
	int unk4[1];
	int ntitl;
};

struct DCDheader
{
	DCDheader1 header1;
	char *title;
	int natom;	
};


static char **res_name = NULL;
static char **at_name = NULL;
static char **seg_name = NULL;
static int  *at_num;
static double *at_q;
static int  *res_num;
static int success = 0;

static int *psf_dihedrals;
static int npsf_dihedrals;

static int *psf_carbons;
static int npsf_carbons;

static int psf_natoms = 0;
static int crd_type = 0;
static DCDheader curHeader;

static int doubleBondSpace = 0;
static int nDoubleBonds = 0;
static int *doubleBonds = NULL;

static int *bondOrder = NULL;

static double FracINV[9];
static double SMUSE[9];
static double SM_aligned[9];

double getXTLABC( double SMOUT[9] )
{
	memcpy( SMOUT, SMUSE, sizeof(double) * 9 );
	return 0;
}

double loadTransform( double FINV[9] )
{
	memcpy( FracINV, FINV, sizeof(double) * 9 );
	return 0;
}
double saveTransform( double FINV[9] )
{
	memcpy( FINV, FracINV, sizeof(double) * 9 );
	return 0;
}

int getNPSFDihedrals( void )
{
	return npsf_dihedrals;
}

int getNPSFCarbons( void )
{
	return npsf_carbons;
}

void putPSFDihedrals( int *dihe_buffer )
{
	memcpy( dihe_buffer, psf_dihedrals, sizeof(int) * 4 * npsf_dihedrals );
}

void putPSFCarbons( int *carbon_buffer )
{
	memcpy( carbon_buffer, psf_carbons, sizeof(int) * 5 * npsf_carbons );
}

int DCDsuccess(void)
{
	return success;
}
void printSingleCRD( FILE *theFile, struct atom_rec *at )
{
}

void printSingleCRD( char *theString, struct atom_rec *at )
{
	int x = 0;

	char printBuffer[1024];

	sprintf(printBuffer, "%d", at[x].bead );

	int field1 = 10;
	int field2 = 10;

	if( strlen(printBuffer) < field1 ) 
	{
                	for( int p = 0; p < field1- strlen(printBuffer); p++ )
				strcat( theString, " " );
//				fprintf(theFile, " ");
	}

//	if( strlen(printBuffer) >= field1 )
//		printf(" ");

	strcat(theString, printBuffer );
//	fprintf(theFile, "%s", printBuffer );
	
	sprintf(printBuffer, "%d", at[x].res );

	if( strlen(printBuffer) < field2 ) 
	{
                	for( int p = 0; p < field2- strlen(printBuffer); p++ )
				strcat(theString, " " );
			//fprintf(theFile, " ");
	}
	if( strlen(printBuffer) >= field2 )
		printf(" ");

//	fprintf(theFile, "%s ", printBuffer );
	strcat(theString, printBuffer);

//	fprintf(theFile, " ");	
	strcat(theString, "  ");
	sprintf(printBuffer, "%s", at[x].resname );

	//fprintf(theFile, "%s", printBuffer );
	strcat(theString, printBuffer);

	if( strlen(printBuffer) < 10 ) 
	{
                	for( int p = 0; p < 10 - strlen(printBuffer); p++ )
				strcat(theString, " ");
			//fprintf(theFile, " ");
	}
	
	sprintf(printBuffer, "%s", at[x].atname );
//		fprintf(theFile,"      ");
	//fprintf(theFile, "%s", printBuffer );
	strcat(theString, printBuffer);

	if( strlen(printBuffer) < 8 ) 
	{
                	for( int p = 0; p < 8 - strlen(printBuffer); p++ )
				strcat(theString, " ");
//			fprintf(theFile, " ");
	}

	{
		sprintf(printBuffer, "%.10lf", at[x].x );

		char tbuf[12];
		int t = 0;
		while( printBuffer[t] != '.' && printBuffer[t] )
		{
			tbuf[t] = printBuffer[t];
			t++;
			tbuf[t] = '\0';
		}
		t++;
		int s = 0;
		char pbuf[24];
		while( printBuffer[t] )
		{
			pbuf[s] = printBuffer[t];
			s++;
			t++;
			pbuf[s] = '\0';
		}

		if( strlen(tbuf) < 9 )
		{
			for( int s = 0; s < 9-strlen(tbuf);s++ )
				strcat(theString, " ");
//				fprintf(theFile, " ");
		} 

		strcat(theString, tbuf);
		strcat(theString, ".");
		strcat(theString, pbuf);
//		fprintf(theFile, "%s.", tbuf );
//		fprintf(theFile,"%s", pbuf );

		if( strlen(pbuf) < 10 )
		{
			for( int s = 0; s < 10-strlen(pbuf);s++ )
				strcat(theString, " ");
//				fprintf(theFile, " ");
		} 
	}
	
	{
		sprintf(printBuffer, "%.10lf", at[x].y );

		char tbuf[12];
		int t = 0;
		while( printBuffer[t] != '.' && printBuffer[t] )
		{
			tbuf[t] = printBuffer[t];
			t++;
			tbuf[t] = '\0';
		}
		t++;
		int s = 0;
		char pbuf[24];
		while( printBuffer[t] )
		{
			pbuf[s] = printBuffer[t];
			s++;
			t++;
			pbuf[s] = '\0';
		}

		if( strlen(tbuf) < 9 )
		{
			for( int s = 0; s < 9-strlen(tbuf);s++ )
				strcat(theString, " ");
//				fprintf(theFile, " ");
		} 

		strcat(theString, tbuf);
		strcat(theString, ".");
		strcat(theString, pbuf);
		//fprintf(theFile, "%s.", tbuf );
		//fprintf(theFile,"%s", pbuf );

		if( strlen(pbuf) < 10 )
		{
			for( int s = 0; s < 10-strlen(pbuf);s++ )
				strcat(theString, " ");
				//fprintf(theFile, " ");
		} 
	}
	
	{
		sprintf(printBuffer, "%.10lf", at[x].z );

		char tbuf[12];
		int t = 0;
		while( printBuffer[t] != '.' && printBuffer[t] )
		{
			tbuf[t] = printBuffer[t];
			t++;
			tbuf[t] = '\0';
		}
		t++;
		int s = 0;
		char pbuf[24];
		while( printBuffer[t] )
		{
			pbuf[s] = printBuffer[t];
			s++;
			t++;
			pbuf[s] = '\0';
		}

		if( strlen(tbuf) < 9 )
		{
			for( int s = 0; s < 9-strlen(tbuf);s++ )
				strcat( theString, " ");
//				fprintf(theFile, " ");
		} 

//		fprintf(theFile, "%s.", tbuf );
//		fprintf(theFile,"%s", pbuf );
		strcat(theString, tbuf);
		strcat(theString, ".");
		strcat(theString, pbuf);

		if( strlen(pbuf) < 10 )
		{
			for( int s = 0; s < 10-strlen(pbuf);s++ )
				strcat( theString, " ");
//				fprintf(theFile, " ");
		} 
	}

	strcat( theString, " ");
//	fprintf(theFile, " ");	
	
	sprintf(printBuffer, " %s", at[x].segid );

//	fprintf(theFile, "%s", printBuffer );
	strcat( theString, printBuffer );

	if( strlen(printBuffer) < 11 ) 
	{
                	for( int p = 0; p < 11- strlen(printBuffer); p++ )
				strcat( theString, " ");
		//	fprintf(theFile, " ");
	}
	
	sprintf(printBuffer, "%d", at[x].segRes );

	//fprintf(theFile, "%s", printBuffer );
	strcat( theString, printBuffer );

	if( strlen(printBuffer) < 16) 
	{
                	for( int p = 0; p < 16- strlen(printBuffer); p++ )
				strcat( theString, " ");
//			fprintf(theFile, " ");
	}

//	fprintf(theFile, "0.0000000000\n");
	strcat(theString, "0.0000000000\n");
}

void printCRDHeader( FILE *theFile, int nat)
{
	fprintf(theFile,
"* GENERATED BY ResetDihedral\n"
"* INPUT FILES FOR EQUILIBRATION\n"
"*  DATE:     5/25/12     16:55:15      CREATED BY USER: alex\n"
"*\n" );

	fprintf(theFile, "%d  EXT\n", nat );
}

void printCRD( FILE *theFile, struct atom_rec *at, int nat)
{
	fprintf(theFile,
"* GENERATED BY ResetDihedral\n"
"* INPUT FILES FOR EQUILIBRATION\n"
"*  DATE:     5/25/12     16:55:15      CREATED BY USER: alex\n"
"*\n" );

	fprintf(theFile, "%d  EXT\n", nat );


	int pres = -1;
	int cres = 0;
	for( int x = 0; x < nat ; x++ )
	{
		if( at[x].res != pres )
			cres++;
		pres = at[x].res;

		char printBuffer[1024];

		sprintf(printBuffer, "%d", x+1 );

		int field1 = 10;
		int field2 = 10;

		if( strlen(printBuffer) < field1 ) 
		{
                	for( int p = 0; p < field1- strlen(printBuffer); p++ )
				fprintf(theFile, " ");
		}

		if( strlen(printBuffer) >= field1 )
			printf(" ");

		fprintf(theFile, "%s", printBuffer );
		
		sprintf(printBuffer, "%d", cres );

		if( strlen(printBuffer) < field2 ) 
		{
                	for( int p = 0; p < field2- strlen(printBuffer); p++ )
				fprintf(theFile, " ");
		}
		if( strlen(printBuffer) >= field2 )
			printf(" ");

		fprintf(theFile, "%s ", printBuffer );

		fprintf(theFile, " ");	
		sprintf(printBuffer, "%s", at[x].resname );

		fprintf(theFile, "%s", printBuffer );

		if( strlen(printBuffer) < 4 ) 
		{
                	for( int p = 0; p < 4 - strlen(printBuffer); p++ )
				fprintf(theFile, " ");
		}
		
		sprintf(printBuffer, "%s", at[x].atname );
		fprintf(theFile,"      ");
		fprintf(theFile, "%s", printBuffer );

		if( strlen(printBuffer) < 8 ) 
		{
                	for( int p = 0; p < 8 - strlen(printBuffer); p++ )
				fprintf(theFile, " ");
		}

		{
			sprintf(printBuffer, "%.10lf", at[x].x );

			char tbuf[12];
			int t = 0;
			while( printBuffer[t] != '.' && printBuffer[t] )
			{
				tbuf[t] = printBuffer[t];
				t++;
				tbuf[t] = '\0';
			}
			t++;
			int s = 0;
			char pbuf[24];
			while( printBuffer[t] )
			{
				pbuf[s] = printBuffer[t];
				s++;
				t++;
				pbuf[s] = '\0';
			}

			if( strlen(tbuf) < 9 )
			{
				for( int s = 0; s < 9-strlen(tbuf);s++ )
					fprintf(theFile, " ");
			} 

			fprintf(theFile, "%s.", tbuf );
			fprintf(theFile,"%s", pbuf );

			if( strlen(pbuf) < 10 )
			{
				for( int s = 0; s < 10-strlen(pbuf);s++ )
					fprintf(theFile, " ");
			} 
		}
		
		{
			sprintf(printBuffer, "%.10lf", at[x].y );

			char tbuf[12];
			int t = 0;
			while( printBuffer[t] != '.' && printBuffer[t] )
			{
				tbuf[t] = printBuffer[t];
				t++;
				tbuf[t] = '\0';
			}
			t++;
			int s = 0;
			char pbuf[24];
			while( printBuffer[t] )
			{
				pbuf[s] = printBuffer[t];
				s++;
				t++;
				pbuf[s] = '\0';
			}

			if( strlen(tbuf) < 9 )
			{
				for( int s = 0; s < 9-strlen(tbuf);s++ )
					fprintf(theFile, " ");
			} 

			fprintf(theFile, "%s.", tbuf );
			fprintf(theFile,"%s", pbuf );

			if( strlen(pbuf) < 10 )
			{
				for( int s = 0; s < 10-strlen(pbuf);s++ )
					fprintf(theFile, " ");
			} 
		}
		
		{
			sprintf(printBuffer, "%.10lf", at[x].z );

			char tbuf[12];
			int t = 0;
			while( printBuffer[t] != '.' && printBuffer[t] )
			{
				tbuf[t] = printBuffer[t];
				t++;
				tbuf[t] = '\0';
			}
			t++;
			int s = 0;
			char pbuf[24];
			while( printBuffer[t] )
			{
				pbuf[s] = printBuffer[t];
				s++;
				t++;
				pbuf[s] = '\0';
			}

			if( strlen(tbuf) < 9 )
			{
				for( int s = 0; s < 9-strlen(tbuf);s++ )
					fprintf(theFile, " ");
			} 

			fprintf(theFile, "%s.", tbuf );
			fprintf(theFile,"%s", pbuf );

			if( strlen(pbuf) < 10 )
			{
				for( int s = 0; s < 10-strlen(pbuf);s++ )
					fprintf(theFile, " ");
			} 
		}
	
		fprintf(theFile, " ");	
		
		sprintf(printBuffer, " %s", at[x].segid );

		fprintf(theFile, "%s", printBuffer );

		if( strlen(printBuffer) < 11 ) 
		{
                	for( int p = 0; p < 11- strlen(printBuffer); p++ )
				fprintf(theFile, " ");
		}
		
		sprintf(printBuffer, "%d", at[x].segRes );

		fprintf(theFile, "%s", printBuffer );

		if( strlen(printBuffer) < 16) 
		{
                	for( int p = 0; p < 16- strlen(printBuffer); p++ )
				fprintf(theFile, " ");
		}

		fprintf(theFile, "0.0000000000\n");
	}
}

void writeNAMDBinary( FILE *theFile, struct atom_rec *at )
{
	fwrite( &psf_natoms, sizeof(int), 1, theFile );

	for( int x = 0; x < psf_natoms; x++ )
	{
		double xyz[3];

		xyz[0] = at[x].x;
		xyz[1] = at[x].y;
		xyz[2] = at[x].z;

		fwrite( xyz, sizeof(double), 3, theFile );		
	}
}

void copyNAMDBinary( FILE *theFile, struct atom_rec *at)
{

	int nat;

	fread( &nat, sizeof(int), 1, theFile );

	if( nat != psf_natoms )
	{
		printf("natoms psf mismatch cor: %d psf: %d.\n", nat, psf_natoms);
		exit(1);
	}

	for( int x = 0; x < nat; x++ )
	{
		double xyz[3];

		fread( xyz, sizeof(double), 3, theFile );		

		at[x].segRes = 0;
		
		at[x].bead = at_num[x];
		at[x].res = res_num[x];
		at[x].resname = (char *)malloc( sizeof(char) * (strlen(res_name[x])+1) );
		at[x].atname = (char *)malloc( sizeof(char) * (strlen(at_name[x])+1) );
		at[x].segid = (char *)malloc( sizeof(char) * (strlen(seg_name[x])+1) );
		at[x].aux = 0;
		at[x].segRes = 0;
		strcpy( at[x].resname, res_name[x] );
		strcpy( at[x].atname, at_name[x] );
		strcpy( at[x].segid, seg_name[x] );
		at[x].altloc = ' ';
		at[x].chain = ' ';

		at[x].x = xyz[0];
		at[x].y = xyz[1];
		at[x].z = xyz[2];
	}
}

void loadCRD( FILE *theFile, struct atom_rec *at)
{
	char buffer[4096];
	getLine( theFile, buffer );
	while( buffer[0] == '*' )
	{
		getLine(theFile, buffer );
	}

	int read_nat = 0;
	sscanf( buffer, "%d", &read_nat );

	if( read_nat != psf_natoms )
	{
		printf("natoms psf mismatch.\n");
		exit(1);
	}

	for( int x = 0; x < psf_natoms; x++ )
	{
		getLine( theFile, buffer );

		int atnum, resnum;
		char resname[1024];
		char atname[1024];
		double lx,ly,lz;
		int segRes;	
		char segName[1024];
		sscanf( buffer," %d %d %s %s %lf %lf %lf %s %d ",
			&atnum, &resnum, resname, atname, &lx, &ly, &lz, (char *)&segName, &segRes );

		at[x].bead = at_num[x];
		at[x].res = res_num[x];
		at[x].resname = (char *)malloc( sizeof(char) * (strlen(res_name[x])+1) );
		at[x].atname = (char *)malloc( sizeof(char) * (strlen(at_name[x])+1) );
		at[x].segid = (char *)malloc( sizeof(char) * (strlen(seg_name[x])+1) );
		at[x].aux = 0;
		at[x].segRes = segRes;
		strcpy( at[x].resname, res_name[x] );
		strcpy( at[x].atname, at_name[x] );
		strcpy( at[x].segid, seg_name[x] );
		at[x].altloc = ' ';
		at[x].chain = ' ';

		at[x].x = lx;
		at[x].y = ly;
		at[x].z = lz;

/*
		if( crd_type == 0 )
		{ // fractional
			at[x].x = SM[0] * lx + SM[1] * ly + SM[2] * lz;
			at[x].y = SM[3] * lx + SM[4] * ly + SM[5] * lz;
			at[x].z = SM[6] * lx + SM[7] * ly + SM[8] * lz;
		}
		else if( crd_type == 1 )
		{ // symmetric
			at[x].x = lx;
			at[x].y = ly;
			at[x].z = lz;
		}	
		else if( crd_type == 2 )
		{
			at[x].x = (SM_aligned[0] * lx + SM_aligned[1] * ly + SM_aligned[2] * lz);
			at[x].y = (SM_aligned[3] * lx + SM_aligned[4] * ly + SM_aligned[5] * lz);
			at[x].z = (SM_aligned[6] * lx + SM_aligned[7] * ly + SM_aligned[8] * lz);
		}
*/
	}
}

void loadPSFfromPDB( FILE *theFile )
{
	char buffer[1024];
	memset( buffer, 0, sizeof(char) * 1024 );

	int local_nat = 0;

	rewind(theFile);

	while( !feof(theFile) ) 
	{
		getLine( theFile, buffer );
		if( feof(theFile) ) break;

		if( !strncasecmp( buffer, "ATOM",4) ) 
			local_nat++;
	}

	rewind(theFile);

	res_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_q = (double *)malloc( sizeof( double ) * local_nat );	
	seg_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_num = (int *)malloc( sizeof(int) * local_nat );
	res_num = (int *)malloc( sizeof(int) * local_nat );

	psf_natoms = local_nat;

	int ncarbonSpace = 0;

	int x=0;
	while( !feof(theFile) ) 
	{
		getLine( theFile, buffer );
		
		if( feof(theFile) ) break;	
	
		if( !strncasecmp( buffer, "ATOM",4) )
		{ 
			struct atom_rec lat;

			readATOM( buffer, &lat );

		
	
			res_name[x] = (char *)malloc( sizeof(char) * ( strlen(lat.resname) +1 ) );
			at_name[x] = (char *)malloc( sizeof(char) * ( strlen(lat.atname) +1 ) );
			if( lat.segid )
				seg_name[x] = (char *)malloc( sizeof(char) * ( strlen(lat.segid) +1 ) );
			else
				seg_name[x] = (char *)malloc( sizeof(char) * ( strlen("NONE") + 1 ) );
	
			strcpy( res_name[x], lat.resname );	
			strcpy( at_name[x], lat.atname );	
			if( lat.segid ) 
				strcpy( seg_name[x], lat.segid );	
			else
				sprintf(seg_name[x], "NONE");
			at_num[x] = lat.bead;
			res_num[x] = lat.res;
			at_q[x] = lat.charge;

			lat.zap();
			x++;
		}
	}	
}


void loadPSF( FILE *theFile )
{
	char buffer[1024];
	memset( buffer, 0, sizeof(char) * 1024 );

	getLine( theFile, buffer );	
	getLine( theFile, buffer );	
	getLine( theFile, buffer );	
	int nlines = 0;
	sscanf( buffer, " %d ", &nlines );

	for( int x = 0; x < nlines; x++ )
		getLine( theFile, buffer );

	getLine( theFile, buffer );	
	getLine( theFile, buffer );	
	
	int local_nat = 0;
	sscanf( buffer, "%d ", &local_nat );
//	printf("Reading %d atoms from the PSF.\n", local_nat );

	res_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_q = (double *)malloc( sizeof( double ) * local_nat );	
	seg_name = (char **)malloc( sizeof( char *) * local_nat );	
	at_num = (int *)malloc( sizeof(int) * local_nat );
	res_num = (int *)malloc( sizeof(int) * local_nat );

	psf_natoms = local_nat;

	int ncarbonSpace = 0;

	for( int x = 0; x < local_nat; x++ )
	{
		getLine( theFile, buffer );
		
		char l_segname[16];
		char l_resname[16];
		int l_resnum;
		int l_atnum;
		int l_atid;
		double l_atq;
		char l_atname[16]; 
		
		sscanf( buffer, " %d %s %d %s %s %d %lf ",
			&l_atnum, l_segname, &l_resnum, l_resname, l_atname, &l_atid, &l_atq );
	
		res_name[x] = (char *)malloc( sizeof(char) * ( strlen(l_resname) +1 ) );
		at_name[x] = (char *)malloc( sizeof(char) * ( strlen(l_atname) +1 ) );
		seg_name[x] = (char *)malloc( sizeof(char) * ( strlen(l_segname) +1 ) );

		strcpy( res_name[x], l_resname );	
		strcpy( at_name[x], l_atname );	
		strcpy( seg_name[x], l_segname );	

		if( at_name[x][0] == 'C' ) ncarbonSpace += 1;

		at_num[x] = l_atnum;
		res_num[x] = l_resnum;
		at_q[x] = l_atq;
	}	
	
	// count the psf carbons.

	int npsf_c_space = 100;
	psf_carbons = (int *)malloc( sizeof(int) * npsf_c_space * 5 );
	npsf_carbons = 0;
	int *link = (int *)malloc( sizeof(int) * local_nat );

	
	for( int x = 0; x < local_nat; x++ )
	{
		link[x] = -1;

		if( at_name[x][0] == 'C' && strlen(at_name[x]) >= 3 )
		{
			char tstr[3];

			tstr[0] = at_name[x][1];
			tstr[1] = '\0';

			int num = atoi(tstr);

			if( num == 2 || num == 3 ) //chain!!
			{
				tstr[0] = at_name[x][2];
				tstr[1] = '\0';
				tstr[2] = '\0';
				if( at_name[x][3] != '\0' )
					tstr[1] = at_name[x][3];

				int seq = atoi(tstr);
				
				if( seq > 2 )
				{
					if( npsf_carbons == npsf_c_space )	
					{
						npsf_c_space *= 2;
						psf_carbons = (int *)realloc( psf_carbons, sizeof(int) * npsf_c_space * 5);	
					} 

					psf_carbons[5*npsf_carbons+0] = x; // the carbon.
					psf_carbons[5*npsf_carbons+1] = -1; // space for H
					psf_carbons[5*npsf_carbons+2] = -1; // space for H
					psf_carbons[5*npsf_carbons+3] = -1; // space for H
					psf_carbons[5*npsf_carbons+4] = -1; // space for H
					link[x] = npsf_carbons;
					npsf_carbons++;
				}
			} 	
		} 
	}

	getLine(theFile, buffer ); // blank line

	while( !feof(theFile) ) 
	{
		getLine(theFile, buffer ); // the read line

		int nunits;
		char unitName[256];
		int nr = sscanf(buffer, " %d !%s\n", &nunits, unitName );

		if( nr != 2 ) break;
		
		if( !strncasecmp( unitName, "NBOND", 5 ) )
		{
			// dihedrals.
			while( !feof(theFile) )
			{
				getLine( theFile, buffer );

				int units[8];
				int nr = sscanf(buffer, " %d %d %d %d %d %d %d %d",
					units+0,
					units+1,
					units+2,
					units+3,
					units+4,
					units+5,
					units+6,
					units+7 );

				for( int p = 0; p < nr; p++ )
					units[p] -= 1;

				for( int p = 0; p < nr/2; p++ )
				{
					int at_1 = units[2*p+0];
					int at_2 = units[2*p+1];

					if( link[at_1] >= 0 ) // && at_name[at_2][0] == 'H' )
					{
						for( int x = 1; x < 5; x++ )
						{
							if( psf_carbons[5*link[at_1]+x] == -1 )
							{
								psf_carbons[5*link[at_1]+x] = at_2;
								break;
							}
						}	
					}		
					if( link[at_2] >= 0 ) //&& at_name[at_1][0] == 'H' )
					{
						for( int x = 1; x < 5; x++ )
						{
							if( psf_carbons[5*link[at_2]+x] == -1 )
							{
								psf_carbons[5*link[at_2]+x] = at_1;
								break;
							}
						}	
					}		

				}

				if( feof(theFile) || strlen(buffer) < 3 ) break;
			}
		}
		else if( !strncasecmp( unitName, "NPHI", 4 ) )
		{
			psf_dihedrals = (int *)malloc( sizeof(int) * 4 * nunits );
			npsf_dihedrals = 0;
	
			// dihedrals.
			while( !feof(theFile) )
			{
				getLine( theFile, buffer );

				int units[8];
				int nr = sscanf(buffer, " %d %d %d %d %d %d %d %d",
					units+0,
					units+1,
					units+2,
					units+3,
					units+4,
					units+5,
					units+6,
					units+7 );

				if( nr == 4 )
				{
					psf_dihedrals[npsf_dihedrals*4+0] = units[0];
					psf_dihedrals[npsf_dihedrals*4+1] = units[1];
					psf_dihedrals[npsf_dihedrals*4+2] = units[2];
					psf_dihedrals[npsf_dihedrals*4+3] = units[3];
					npsf_dihedrals += 1;
				}	
				else if( nr == 8 )
				{
					psf_dihedrals[npsf_dihedrals*4+0] = units[0];
					psf_dihedrals[npsf_dihedrals*4+1] = units[1];
					psf_dihedrals[npsf_dihedrals*4+2] = units[2];
					psf_dihedrals[npsf_dihedrals*4+3] = units[3];
					npsf_dihedrals += 1;
					psf_dihedrals[npsf_dihedrals*4+0] = units[4];
					psf_dihedrals[npsf_dihedrals*4+1] = units[5];
					psf_dihedrals[npsf_dihedrals*4+2] = units[6];
					psf_dihedrals[npsf_dihedrals*4+3] = units[7];
					npsf_dihedrals += 1;
				}	

				if( feof(theFile) || strlen(buffer) < 3 ) break;
			}
		}
		else
		{
			while( !feof(theFile) )
			{
				getLine( theFile, buffer );

				if( feof(theFile) || strlen(buffer) < 3 ) break;
			}
		}		
	}

	bondOrder = (int*) malloc( sizeof(int) * local_nat );

	for( int x = 0; x < local_nat; x++ )
		bondOrder[x] = 0;

	for( int c = 0; c < npsf_carbons; c++ )
	{
		int na = 0;

		for( int p = 1; p < 5; p++ )
			if( psf_carbons[5*c+p] != -1 ) na++;	

		bondOrder[psf_carbons[5*c+0]] = na;
	}

	
	nDoubleBonds = 0;	
	doubleBondSpace = 100;
	doubleBonds = (int *)malloc( sizeof(int) * 4 * doubleBondSpace );

	for( int x = 0; x < npsf_dihedrals; x++ )
	{
		int at1 = psf_dihedrals[4*x+0]-1;
		int at2 = psf_dihedrals[4*x+1]-1;
		int at3 = psf_dihedrals[4*x+2]-1;
		int at4 = psf_dihedrals[4*x+3]-1;

		if( at_name[at1][0] == 'C' && at_name[at2][0] == 'C' && at_name[at3][0] == 'C' && at_name[at4][0] == 'C' &&
			bondOrder[at2] == 3 && bondOrder[at3] == 3 ) // a double bond
		{
			if( nDoubleBonds == doubleBondSpace )
			{
				doubleBondSpace *= 2;
				doubleBonds = (int *)realloc( doubleBonds, sizeof(int) * 4 * doubleBondSpace );
			}

			doubleBonds[4*nDoubleBonds+0] = at1;
			doubleBonds[4*nDoubleBonds+1] = at2;
			doubleBonds[4*nDoubleBonds+2] = at3;
			doubleBonds[4*nDoubleBonds+3] = at4;

			nDoubleBonds++;
		} 
	}
}

void getDoubleBonds( int **dbonds, int *nbonds )
{
	*dbonds = (int *)malloc( sizeof(int) * nDoubleBonds * 4 );
	*nbonds = nDoubleBonds;
	memcpy( *dbonds, doubleBonds, sizeof(int) * 4 * nDoubleBonds );
}

static double celv = 0;

double CellVolume( void )
{
	return celv;	
}

int pbc_set = 0;
double cstats[6] = { 0,0,0,0,0,0};

int PBCD( double *Lx, double *Ly, double *Lz, double *alpha, double *beta, double *gamma )
{

	*Lx = cstats[0];
	*Ly = cstats[1];
	*Lz = cstats[2];
	*alpha = cstats[3];
	*beta = cstats[4];
	*gamma = cstats[5];

	if( pbc_set )
		return 0;
	return 1;
}

void TransformFractional( double *dr )
{
	double tdr[3] = {0,0,0};

	tdr[0] = FracINV[0] * dr[0] + FracINV[1] * dr[1] + FracINV[2] * dr[2];
	tdr[1] = FracINV[3] * dr[0] + FracINV[4] * dr[1] + FracINV[5] * dr[2];
	tdr[2] = FracINV[6] * dr[0] + FracINV[7] * dr[1] + FracINV[8] * dr[2];
	
	dr[0] = tdr[0];
	dr[1] = tdr[1];
	dr[2] = tdr[2];
}

void TransformFractionalAligned( double *dr )
{
	double tdr[3] = {0,0,0};

	tdr[0] = FracINV[0] * dr[0] + FracINV[1] * dr[1] + FracINV[2] * dr[2];
	tdr[1] = FracINV[3] * dr[0] + FracINV[4] * dr[1] + FracINV[5] * dr[2];
	tdr[2] = FracINV[6] * dr[0] + FracINV[7] * dr[1] + FracINV[8] * dr[2];
	
	dr[0] = SM_aligned[0] * tdr[0] + SM_aligned[1] * tdr[1] + SM_aligned[2] * tdr[2];
	dr[1] = SM_aligned[3] * tdr[0] + SM_aligned[4] * tdr[1] + SM_aligned[5] * tdr[2];
	dr[2] = SM_aligned[6] * tdr[0] + SM_aligned[7] * tdr[1] + SM_aligned[8] * tdr[2];
}

int curNFrames( void )
{
//	printf("curHeader.header1.nstep : %d curHeader.header1.nsavc: %d\n", curHeader.header1.nstep, curHeader.header1.nsavc );
//	if( curHeader.header1.nstep < curHeader.header1.nsavc )
		return curHeader.header1.nstep;
	return curHeader.header1.nstep / curHeader.header1.nsavc;
}

int curNAtoms( void )
{
	return psf_natoms;
}

void setFractional( void )
{
	crd_type = 0;
}

void setSymmetric( void )
{
	crd_type = 1;
}

void setAligned( void )
{
	crd_type = 2;
}


void readDCDHeader( FILE *theFile )
{
	fread( &(curHeader.header1), sizeof(DCDheader1), 1, theFile ); 
//	printf("QCRYS: %d DIM4: %d QCG: %d\n", curHeader.header1.qcrys, curHeader.header1.dim4, curHeader.header1.qcg );	
//	printf("title len: %d\n", curHeader.header1.ntitl );

	curHeader.header1.ntitl *= 80;

	curHeader.title = (char *)malloc( sizeof(char) * (curHeader.header1.ntitl+1) );

	fread( (curHeader.title), sizeof(char), curHeader.header1.ntitl, theFile );
	curHeader.title[curHeader.header1.ntitl] = '\0';
//	printf("title: %s\n", curHeader.title );	
	int junk;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	fread( &(curHeader.natom), sizeof(int), 1, theFile );
	
	fread( &junk, sizeof(int), 1, theFile );

//	printf("natoms: %d\n", curHeader.natom );
}

void loadFrame( FILE *theFile, struct atom_rec *at )
{
	double a[3];
	double b[3];
	double c[3];
	
	int junk;

	if( curHeader.header1.qcrys )
	{
		fread( &junk, sizeof(int), 1, theFile );
		
		fread( &(a[0]), sizeof(double), 1, theFile );
		fread( &(a[1]), sizeof(double), 1, theFile );
		fread( &(b[1]), sizeof(double), 1, theFile );
		fread( &(a[2]), sizeof(double), 1, theFile );
		fread( &(b[2]), sizeof(double), 1, theFile );
		fread( &(c[2]), sizeof(double), 1, theFile );
	}
	else
	{
		a[0] = 1.0;
		a[1] = 0.0;
		a[2] = 0.0;
		b[1] = 1.0;
		b[2] = 0.0;
		c[2] = 1.0;
	}
	
	b[0] = a[1];
	c[0] = a[2];
	c[1] = b[2];
	
	double La,Lb,Lc;
	
	if( fabs(a[1]-90) < 1e-10 && fabs(b[2]-90) < 1e-10 )
	{
		La = a[0];
		Lb = b[1];
		Lc = c[2];

		a[1] = 0;
		a[2] = 0;
		b[0] = 0;
		b[2] = 0;
		c[0] = 0;
		c[1] = 0;
	}
	else
	{
		La = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
		Lb = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
		Lc = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
	}

	double gamma = (180.0/M_PI)*acos( (a[0]*b[0] + a[1]*b[1]+a[2]*b[2]) / (La*Lb) );
	if( fabs(a[1]) < 0.52 && fabs(a[1]) > 0.48 )
	{
		gamma = (180.0/M_PI)*acos( a[1] );
		double lorig = a[0];
		a[0] = lorig * cos( -15.0 * ( M_PI / 180.0) );
		a[1] = lorig * sin( -15.0 * ( M_PI / 180.0) );
		b[1] = lorig * cos( -15.0 * ( M_PI / 180.0) );
		b[0] = lorig * sin( -15.0 * ( M_PI / 180.0) );
	}
	double alpha = (180.0/M_PI)*acos( (c[0]*b[0] + c[1]*b[1]+c[2]*b[2]) / (Lc*Lb) );
	double beta  = (180.0/M_PI)*acos( (a[0]*c[0] + a[1]*c[1]+a[2]*c[2]) / (La*Lc) );

	cstats[0] = La;
	cstats[1] = Lb;
	cstats[2] = Lc;
	cstats[3] = alpha;
	cstats[4] = beta;
	cstats[5] = gamma;

	pbc_set = 1;	

	celv = La*Lb*Lc*sin((M_PI/180.0)*gamma); // currently for hexagonal or orthorhombic.


	double SM[9];


	FracINV[0] = SM[0] = a[0];
	FracINV[1] = SM[1] = a[1];
	FracINV[2] = SM[2] = a[2];
	FracINV[3] = SM[3] = b[0];
	FracINV[4] = SM[4] = b[1];
	FracINV[5] = SM[5] = b[2];
	FracINV[6] = SM[6] = c[0];
	FracINV[7] = SM[7] = c[1];
	FracINV[8] = SM[8] = c[2];


	memcpy( SMUSE, SM, sizeof(double) * 9 );

	char uplo = 'U';
	int N = 3;
	int IPIV[N];
	int lwork = N*N;
	double work[lwork];
	int info;
	dsytrf_( &uplo, &N, SM, &N, IPIV, work, &lwork, &info ); 
	dsytri_( &uplo, &N, SM, &N, IPIV, work, &info );

	SM[1] = SM[3];
	SM[2] = SM[6];
	SM[5] = SM[7];
/*
	for( int x = 0; x < 3; x++ )
	{
		for( int y = 0; y < 3; y++ )
		{
			printf("%lf ", SM[x*3+y] );
		}
		printf("\n");
	}

	printf("%lf %lf %lf :: %lf %lf %lf\n", La, Lb, Lc, alpha, beta, gamma );
*/
	// to align shape matrix, rotate around Z to diagonalize.

	double theta = atan2( -SM[1], SM[0] );

	SM_aligned[0] = cos(theta);
	SM_aligned[1] = sin(theta);
	SM_aligned[2] = 0;
	
	SM_aligned[3] = -sin(theta);
	SM_aligned[4] =  cos(theta);
	SM_aligned[5] = 0;
	
	SM_aligned[6] = 0;
	SM_aligned[7] = 0;
	SM_aligned[8] = 1;


#define TRAJ_TYPE float
	int natom = curHeader.natom;

	TRAJ_TYPE *loc_x = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_y = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_z = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	success = 1;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	size_t lread = fread( loc_x, sizeof(TRAJ_TYPE), natom, theFile );
        if( lread != natom )
		success = 0;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	lread = fread( loc_y, sizeof(TRAJ_TYPE), natom, theFile );
        if( lread != natom )
		success = 0;
	fread( &junk, sizeof(int), 1, theFile );
	fread( &junk, sizeof(int), 1, theFile );
	lread = fread( loc_z, sizeof(TRAJ_TYPE), natom, theFile );
        if( lread != natom )
		success = 0;



	for( int x = 0; x < natom; x++ )
	{
		double lx = loc_x[x];
		double ly = loc_y[x];
		double lz = loc_z[x];

		at[x].bead = at_num[x];
		at[x].charge = at_q[x];
		at[x].res = res_num[x];
		at[x].resname = (char *)malloc( sizeof(char) * (strlen(res_name[x])+1) );
		at[x].atname = (char *)malloc( sizeof(char) * (strlen(at_name[x])+1) );
		at[x].segid = (char *)malloc( sizeof(char) * (strlen(seg_name[x])+1) );
		at[x].aux = 0;

		strcpy( at[x].resname, res_name[x] );
		strcpy( at[x].atname, at_name[x] );
		strcpy( at[x].segid, seg_name[x] );
		at[x].altloc = ' ';
		at[x].chain = ' ';
		if( crd_type == 0 )
		{ // fractional
			at[x].x = SM[0] * lx + SM[1] * ly + SM[2] * lz;
			at[x].y = SM[3] * lx + SM[4] * ly + SM[5] * lz;
			at[x].z = SM[6] * lx + SM[7] * ly + SM[8] * lz;
		}
		else if( crd_type == 1 )
		{ // symmetric
			at[x].x = lx;
			at[x].y = ly;
			at[x].z = lz;
		}	
		else if( crd_type == 2 )
		{
			at[x].x = (SM_aligned[0] * lx + SM_aligned[1] * ly + SM_aligned[2] * lz);
			at[x].y = (SM_aligned[3] * lx + SM_aligned[4] * ly + SM_aligned[5] * lz);
			at[x].z = (SM_aligned[6] * lx + SM_aligned[7] * ly + SM_aligned[8] * lz);
		}
	}
		
	fread( &junk, sizeof(int), 1, theFile );

	free(loc_x);
	free(loc_y);
	free(loc_z);
}




/*
int main( int argc, char **argv )
{
	FILE *theFile = fopen(argv[1], "r");
	FILE *psfFile = fopen(argv[2], "r");

	loadPSF( psfFile );
	
	readDCDHeader( theFile );
		
	int nframes = curHeader.header1.nstep / curHeader.header1.nsavc;
	
	printf("natoms: %d\n", curHeader.natom );
	printf("nsteps: %d\n", curHeader.header1.nstep );
	printf("nsavc:  %d\n", curHeader.header1.nsavc );
	printf("nframes: %d\n", nframes );

	int natom = curHeader.natom;

	for( int x = 0; x < nframes; x++ )
	{
		struct atom_rec *at = (atom_rec *)malloc( sizeof(atom_rec) * natom );

		loadFrame( theFile, at );

		free(at);
	}

}DD*/

void loadPDB( FILE *theFile, struct atom_rec *at) 
{
        char buffer[4096];

        int a = 0; 

        while( !feof(theFile) )
        {    
                getLine( theFile, buffer );
                if( feof(theFile) ) break;

		if( !strncasecmp( buffer, "CRYST", 5 ) )
		{
			double LA,LB,LC;
			sscanf( buffer, "CRYST1 %lf %lf %lf", &LA, &LB, &LC );

			cstats[0] = LA;
			cstats[1] = LB;
			cstats[2] = LC;
			cstats[3] = 90;
			cstats[4] = 90;
			cstats[5] = 90;
	
			pbc_set = 1;	
		}

                if( !strncasecmp( buffer, "ATOM", 4 ) )
                {    
                        readATOM( buffer, at+a );

			if( res_name && res_name[a])
			{
				free(at[a].resname);
				at[a].resname = (char *)malloc( sizeof(char) * (strlen(res_name[a])+1) );
				strcpy( at[a].resname, res_name[a] );
			}
			
			if( at_name && at_name[a])
			{
				free(at[a].atname);
				at[a].atname = (char *)malloc( sizeof(char) * (strlen(at_name[a])+1) );
				strcpy( at[a].atname, at_name[a] );
			}
			
			if( seg_name && seg_name[a])
			{
				free(at[a].segid);
				at[a].segid = (char *)malloc( sizeof(char) * (strlen(seg_name[a])+1) );
				strcpy( at[a].segid, seg_name[a] );
			}

                        a++; 
                }     
                if( a >= psf_natoms) break;     
        }    
}

int loadPDB( FILE *theFile, struct atom_rec *at, int max_space) 
{
        char buffer[4096];

        int a = 0; 

        while( !feof(theFile) )
        {    
                getLine( theFile, buffer );
                if( feof(theFile) ) break;

		if( !strncasecmp( buffer, "CRYST", 5 ) )
		{
			double LA,LB,LC;
			sscanf( buffer, "CRYST1 %lf %lf %lf", &LA, &LB, &LC );

			cstats[0] = LA;
			cstats[1] = LB;
			cstats[2] = LC;
			cstats[3] = 90;
			cstats[4] = 90;
			cstats[5] = 90;
			
			pbc_set = 1;	
		}

                if( !strncasecmp( buffer, "ATOM", 4 ) )
                {    
                        readATOM( buffer, at+a );
                        a++; 
                }     
                if( a >= max_space ) break;     
        }

	return a;    
}

void copyDCDHeader( FILE *theFile1, FILE *theFile2, int activateQCRYS )
{
	int nr = fread( &(curHeader.header1), sizeof(DCDheader1), 1, theFile1 ); 

	if( activateQCRYS) 
		curHeader.header1.qcrys = 1;

	fwrite( &(curHeader.header1), sizeof(DCDheader1), 1, theFile2 ); 

	curHeader.title = (char *)malloc( sizeof(char) * (curHeader.header1.ntitl*80+1) );
	nr = fread( (curHeader.title), sizeof(char), curHeader.header1.ntitl*80, theFile1 );
	fwrite( (curHeader.title), sizeof(char), curHeader.header1.ntitl*80, theFile2 );

	int junk;
	 nr = fread( &junk, sizeof(int), 1, theFile1 );
	printf("header junk1: %d\n", junk );
	fwrite( &junk, sizeof(int), 1, theFile2);
	nr = fread( &junk, sizeof(int), 1, theFile1);
	printf("header junk2: %d\n", junk );
	fwrite( &junk, sizeof(int), 1, theFile2 );

	nr= fread( &(curHeader.natom), sizeof(int), 1, theFile1 );
	fwrite( &(curHeader.natom), sizeof(int), 1, theFile2 );
	
	nr = fread( &junk, sizeof(int), 1, theFile1 ); // natom footer
	printf("header junk3: %d\n", junk );
	fwrite( &junk, sizeof(int), 1, theFile2 );
	
}


void copyFrameNewCoords( FILE *theFile1, FILE *theFile2, double *coords )
{
	static double use_abc[6];
	static int abc_init = 0;

	double a[3];
	double b[3];
	double c[3];
	
	int junk;

	int nr ;
	double La=1.0;
	double Lb=1.0;	
	double Lc=1.0;

	double gamma = 0;
	double beta  = 0;
	double alpha = 0;
/*
	// the length of vec_A
	fread( &(a[0]), sizeof(double), 1, theFile1 );
	// the angle between vec_A and vec_B
	fread( &(a[1]), sizeof(double), 1, theFile1 );
	// the length of vec_B
	fread( &(b[1]), sizeof(double), 1, theFile1 );
	// the angle between ? .. let's call it vec_A and vec_C.
	fread( &(a[2]), sizeof(double), 1, theFile1 );
	// the angle between ? .. let's call it vec_B and vec_C.
	fread( &(b[2]), sizeof(double), 1, theFile1 );
	// the length of vec_C
	fread( &(c[2]), sizeof(double), 1, theFile1 );

	double xtlabc[6] = { a[0], a[1], b[1], a[2], b[2], c[2] };
	
	if( abc_init == 0 )
		memcpy( use_abc, xtlabc, sizeof(double) * 6 );
	abc_init = 1;
	memcpy( xtlabc, use_abc, sizeof(double) * 6 );
	fwrite( xtlabc, sizeof(double), 6, theFile2 );
	*/
	if( curHeader.header1.qcrys )
	{
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
		
		fread( &(a[0]), sizeof(double), 1, theFile1 );
		fwrite( &(a[0]), sizeof(double), 1, theFile2 );
		fread( &(a[1]), sizeof(double), 1, theFile1 );
		fwrite( &(a[1]), sizeof(double), 1, theFile2 );
		fread( &(b[1]), sizeof(double), 1, theFile1 );
		fwrite( &(b[1]), sizeof(double), 1, theFile2 );
		fread( &(a[2]), sizeof(double), 1, theFile1 );
		fwrite( &(a[2]), sizeof(double), 1, theFile2 );
		fread( &(b[2]), sizeof(double), 1, theFile1 );
		fwrite( &(b[2]), sizeof(double), 1, theFile2 );
		fread( &(c[2]), sizeof(double), 1, theFile1 );
		fwrite( &(c[2]), sizeof(double), 1, theFile2 );
	}
	
#define TRAJ_TYPE float

	int natom = curHeader.natom;

	TRAJ_TYPE *loc_x = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_y = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	TRAJ_TYPE *loc_z = (TRAJ_TYPE *)malloc( sizeof(TRAJ_TYPE) * natom );
	{	
	int junk[3][2];
	
	fread( &(junk[0][0]), sizeof(int), 1, theFile1 );
	fread( &(junk[0][1]), sizeof(int), 1, theFile1 );
	size_t lread = fread( loc_x, sizeof(TRAJ_TYPE), natom, theFile1 );
	fread( &(junk[1][0]), sizeof(int), 1, theFile1 );
	fread( &(junk[1][1]), sizeof(int), 1, theFile1 );
	lread = fread( loc_y, sizeof(TRAJ_TYPE), natom, theFile1 );
	fread( &(junk[2][0]), sizeof(int), 1, theFile1 );
	fread( &(junk[2][1]), sizeof(int), 1, theFile1 );
	lread = fread( loc_z, sizeof(TRAJ_TYPE), natom, theFile1 );

	for( int a = 0; a < natom; a++ )
	{
		loc_x[a] = (float)coords[3*a+0]; 
		loc_y[a] = (float)coords[3*a+1]; 
		loc_z[a] = (float)coords[3*a+2]; 
	}
		
	fwrite( &(junk[0][0]), sizeof(int), 1, theFile2 );
	fwrite( &(junk[0][1]), sizeof(int), 1, theFile2 );
	fwrite( loc_x, sizeof(TRAJ_TYPE), natom, theFile2 );
	fwrite( &(junk[1][0]), sizeof(int), 1, theFile2 );
	fwrite( &(junk[1][1]), sizeof(int), 1, theFile2 );
	fwrite( loc_y, sizeof(TRAJ_TYPE), natom, theFile2 );
	fwrite( &(junk[2][0]), sizeof(int), 1, theFile2 );
	fwrite( &(junk[2][1]), sizeof(int), 1, theFile2 );
	fwrite( loc_z, sizeof(TRAJ_TYPE), natom, theFile2 );
	}
	free(loc_x);
	free(loc_y);
	free(loc_z);

/*	
	for( int x = 0; x < sizeof(int); x++ )
	{
		char cjunk;
		nr = fread( &cjunk, sizeof(char), 1, theFile1 );
		fwrite( &cjunk, sizeof(char), 1, theFile2 );
	}
*/
		fread( &junk, sizeof(int), 1, theFile1 );
		fwrite( &junk, sizeof(int), 1, theFile2 );
}




void dcd_MinImage( double *f_target, double *f_pt, double cell[3][3] )
{
	int done = 0;

	double cur_dr[3] = { f_target[0] - f_pt[0], f_target[1] - f_pt[1], f_target[2] - f_pt[2] };
	
	double cur_r2 = cur_dr[0]*cur_dr[0]+cur_dr[1]*cur_dr[1]+cur_dr[2]*cur_dr[2];

	while( !done )
	{
		done = 1;

		double displ[4][3] =
		{
			{cell[0][0], cell[0][1], cell[0][2]},
			{cell[1][0], cell[1][1], cell[1][2]},
			{cell[1][0]+cell[0][0], cell[1][1]+cell[0][1], cell[1][2]+cell[0][2]},
			{cell[2][0], cell[2][1], cell[2][2]},
		};

		for( int sgn = -1; sgn <= 1; sgn++ )
		for( int attr = 0; attr < 4; attr++ )
		{
			if( sgn == 0 ) continue;

			double trial_dr[3] = { cur_dr[0] + sgn * displ[attr][0],
						 cur_dr[1] + sgn * displ[attr][1],
						 cur_dr[2] + sgn * displ[attr][2] };

			double trial_r2 = trial_dr[0]*trial_dr[0]+trial_dr[1]*trial_dr[1]+trial_dr[2]*trial_dr[2];

			if( trial_r2 < cur_r2 - (1e-8) )
			{
				done = 0;

				cur_dr[0] += sgn * displ[attr][0];
				cur_dr[1] += sgn * displ[attr][1];
				cur_dr[2] += sgn * displ[attr][2];
				cur_r2 = trial_r2;
			}
		}
	}

	f_target[0] = cur_dr[0] + f_pt[0];
	f_target[1] = cur_dr[1] + f_pt[1];
	f_target[2] = cur_dr[2] + f_pt[2];
}

static char *cached_header = NULL;
static int cached_pt = 0;
void cacheDCDHeader( FILE *theFile1 )
{
	if( cached_header )
	{
		free(cached_header);
		cached_pt = 0;

	}
	int nr = fread( &(curHeader.header1), sizeof(DCDheader1), 1, theFile1 ); 


	cached_header = (char *)malloc( sizeof(char) * ( sizeof(DCDheader1) + curHeader.header1.ntitl*80 + sizeof(int) * 4 + 1000 ) );
	memcpy( cached_header+cached_pt, &(curHeader.header1), sizeof(DCDheader1) );
 //	fwrite( &(curHeader.header1), sizeof(DCDheader1), 1, theFile2 ); 
	cached_pt += sizeof(DCDheader1);

	curHeader.title = (char *)malloc( sizeof(char) * (curHeader.header1.ntitl*80+1) );
	nr = fread( (curHeader.title), sizeof(char), curHeader.header1.ntitl*80, theFile1 );

	memcpy( cached_header+cached_pt,  (curHeader.title), curHeader.header1.ntitl*80 );
	cached_pt += curHeader.header1.ntitl*80;

//	fwrite( (curHeader.title), sizeof(char), curHeader.header1.ntitl*80, theFile2 );

	int junk;

	nr = fread( &junk, sizeof(int), 1, theFile1 );
	memcpy( cached_header+cached_pt, &junk, sizeof(int) );
	cached_pt += sizeof(int);

	//fwrite( &junk, sizeof(int), 1, theFile2);
	
	nr = fread( &junk, sizeof(int), 1, theFile1);
	memcpy( cached_header+cached_pt, &junk, sizeof(int) );
	cached_pt += sizeof(int);
	//fwrite( &junk, sizeof(int), 1, theFile2 );

	nr= fread( &(curHeader.natom), sizeof(int), 1, theFile1 );
	memcpy( cached_header+cached_pt, &(curHeader.natom), sizeof(int) );
	cached_pt += sizeof(int);
	//fwrite( &(curHeader.natom), sizeof(int), 1, theFile2 );
	
	nr = fread( &junk, sizeof(int), 1, theFile1 ); // natom footer
	memcpy( cached_header+cached_pt, &junk, sizeof(int) );
	cached_pt += sizeof(int);
	//fwrite( &junk, sizeof(int), 1, theFile2 );
	
}

void uncacheDCDHeader( FILE *theFile1 )
{
	fwrite( cached_header, 1, cached_pt, theFile1 );
}

/*
int fetchPSFCycles( cycle **outCycles )
{

} */
