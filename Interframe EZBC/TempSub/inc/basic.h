#ifndef __BASICH__
#define __BASICH__

#include <stdio.h>

#define U8	    unsigned char       /* unsigned 8 bit integer        (FF hex) */
#define U16     unsigned short int      /* unsigned 16 bit integer (FFFF hex) */
#define U32	    unsigned int        /* unsigned 32 bit integer       (FFFFFFFF hex) */
#define S32	    signed   int        /* signed 32 bit integer (80000000 hex) */
#define R32	    float       /*32 bit real number     (7F800000 hex (+ infinity)) */
#define ASCII	char            /*(NULL  00 hex) */

typedef unsigned char CharImg;
typedef int IntImg;

#define MY_SQUARE(A) (A)*(A)
#define MY_SQUARE_NORM(A,B) (MY_SQUARE(A)) + (MY_SQUARE(B))
#define MY_ABS(A) (((A) < 0) ? ((-1)*(A)) : (A))
#define MY_ABS_DIFF(A,B) (MY_ABS((A) - (B)))
#define MY_MAX(A,B) ((A) > (B) ? (A) : (B))
#define MY_MIN(A,B) ((A) > (B) ? (B) : (A))
#define MY_SIGN(A)  ((A) < 0 ? -1 : 1)
#define MY_DELTA(A,B) (((A) == (B)) ? (1) : (0))
#define nint(a) (((a)<0) ? ((int)(a - 0.5)) : ((int)(a + 0.5))) // Ô¶ÁãÈ¡Õû

int scanargs( int argc, char *argv[], char *format, ... );
struct rasterfile *read_charimg( char *name, unsigned char **img );
struct rasterfile *read_floatimg( char *name, float **img );
char *getarray( int num, int siz, char *name );
void write_charimg( char *name, struct rasterfile *rhead,
                    unsigned char *data );
void write_floatimg1( char *name, struct rasterfile *rhead, float *data );
void write_floatimg4( char *name, struct rasterfile *rhead, float *data );
void cnv_data_1_4( unsigned char *y, float *x, int n );
void cnv_data_4_1( float *y, unsigned char *x, int n );
void cnv_data_1_2( unsigned char *y, int *x, int n );
void cnv_data_2_1( int *y, unsigned char *x, int n );
float **getmatrix( int nrl, int nrh, int ncl, int nch );
void freematrix( float **m, int nrl, int nrh, int ncl, int nch );
void write_intimg1( char *name, struct rasterfile *rhead, int *data );
struct rasterfile *read_intimg( char *name, int **img );
void writeFloatImg1( char *name, struct rasterfile *rhead, float *data );
/*  Std ANSI prototypes (missing in SunOS include files)  */

#endif
