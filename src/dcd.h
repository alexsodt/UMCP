#ifndef __dcd_h__
#define __dcd_h__

#include "pdb.h"

void getBonds( int **bonds, int *nbonds );
void getDoubleBonds( int **dbonds, int *nbonds );
int getNPSFCarbons( void );
void putPSFCarbons( int *carbon_buffer );
void loadPSF( FILE *theFile );
void loadPSFfromPDB( FILE *theFile );
int PBCD( double *Lx, double *Ly, double *Lz, double *alpha, double *beta, double *gamma );
int curNAtoms( void );
int curNFrames( void );
void setFractional( void );
void setSymmetric( void );
void setAligned( void );
void readDCDHeader( FILE *theFile );
void loadFrame( FILE *theFile, struct atom_rec *at );
void TransformFractional( double *dr );
void TransformFractionalAligned( double *dr );
double CellVolume( void );
int DCDsuccess(void);
double loadTransform( double *);
double saveTransform( double *);
int getNPSFDihedrals( void );
void putPSFDihedrals( int *dihe_buffer );
void loadCRD( FILE *theFile, struct atom_rec *at);
void loadPDB( FILE *theFile, struct atom_rec *at);
void printCRD( FILE *theFile, struct atom_rec *at, int nat);
void printSingleCRD( FILE *theFile, struct atom_rec *at );
void printSingleCRD( char *theString, struct atom_rec *at );
void printCRDHeader( FILE *theFile, int nat);
double getXTLABC( double SMOUT[9] );
void copyNAMDBinary( FILE *theFile, struct atom_rec *at);
void writeNAMDBinary( FILE *theFile, struct atom_rec *at );
int loadPDB( FILE *theFile, struct atom_rec *at, int max_space);


void copyDCDHeader( FILE *theFile1, FILE *theFile2, int activateQCRYS );
void copyFrameImage( FILE *theFile1, FILE *theFile2 );
void dcd_MinImage( double *f_target, double *f_pt, double cell[3][3] );
void copyFrameNewCoords( FILE *theFile1, FILE *theFile2, double *coords );
void cacheDCDHeader( FILE *theFile1 );
void uncacheDCDHeader( FILE *theFile1 );

#endif
