// this file contains functions for down/up-sampling chrominance compant by MPEG half band filters.
/* f444_422 f444_420 f420_444 f422_444 up/down-sample chrominance components. 
The size of chromiance component will change correspondingly, but is not returned to its parent function*/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "structN.h"
#include "basic.h"
#include "util_filtering.h"


void
YUV2RGB( float Y, float U, float V, float *R, float *G, float *B )
{
  *R = ( float )( 1.000 * Y - 0.0009 * U + 1.4017 * V );
  *G = ( float )( 1.000 * Y - 0.3437 * U - 0.7142 * V );
  *B = ( float )( 1.000 * Y + 1.7722 * U + 0.0010 * V );
/*  *R = 1.000*Y - 0.0009*V + 1.4017*U; 
  *G = 1.000*Y - 0.3437*V - 0.7142*U; 
  *B = 1.000*Y + 1.7722*V + 0.0010*U;*/

}

void
RGB2YUV( float R, float G, float B, float *Y, float *U, float *V )
{
  *Y = ( float )( 0.299 * R + 0.587 * G + 0.114 * B );
  *U = ( float )( -0.169 * R - 0.331 * G + 0.500 * B );
  *V = ( float )( 0.500 * R - 0.419 * G - 0.081 * B );
}


// 240M
void
YUV2RGB_( float Y, float U, float V, float *R, float *G, float *B )
{
  *R = ( float )( 1.000 * Y - 0.0007 * U + 1.5758 * V );
  *G = ( float )( 1.000 * Y - 0.2264 * U - 0.4765 * V );
  *B = ( float )( 1.000 * Y + 1.8260 * U - 0.0004 * V );

}

//240M
void
RGB2YUV_( float R, float G, float B, float *Y, float *U, float *V )
{
  *Y = ( float )( 0.212 * R + 0.701 * G + 0.087 * B );
  *U = ( float )( -0.116 * R - 0.384 * G + 0.500 * B );
  *V = ( float )( 0.500 * R - 0.445 * G - 0.055 * B );
}





void
f444_422( videoinfo info, YUVimage * frame )
{
  int i, j, chor, cver, chor_half;
  int filter_length, half_length, new_length;
  float *h = NULL, *line = NULL, *extension = NULL, *outline = NULL;
  float *U = NULL, *V = NULL;

  filter_length = 7;
  h = ( float * )getarray( filter_length, sizeof( float ), "h" );
  h = &h[filter_length / 2];
  h[-3] = h[3] = -29 / 256.;
  h[-2] = h[2] = 0.;
  h[-1] = h[1] = 88 / 256.;
  h[0] = 138 / 256.;

  if( info.cwidth % 2 ) {
    printf( "can not handle this case!\n" );
    exit( 1 );
  }
  chor = info.cwidth;
  chor_half = chor / 2;
  cver = info.cheight;

  U = ( float * )getarray( chor_half * cver, sizeof( float ), "U" );
  V = ( float * )getarray( chor_half * cver, sizeof( float ), "V" );

  outline = ( float * )getarray( chor, sizeof( float ), "outline" );


  half_length = ( filter_length - 1 ) / 2;
  new_length = chor + 2 * half_length;
  extension = ( float * )getarray( new_length, sizeof( float ), "extension" );
  extension = &extension[half_length];



  for( i = 0; i < cver; i++ ) {
    line = &( frame->U[i * chor] );
    line_convolve( line, extension, chor, h, filter_length, ANAL, outline );
    for( j = 0; j < chor_half; j++ )
      U[i * chor_half + j] = outline[j * 2];

    line = &( frame->V[i * chor] );
    line_convolve( line, extension, chor, h, filter_length, ANAL, outline );
    for( j = 0; j < chor_half; j++ )
      V[i * chor_half + j] = outline[j * 2];

  }

  free( frame->U );
  free( frame->V );
  free( ( char * )&extension[-half_length] );

  frame->U = U;
  frame->V = V;

  free( outline );

  free( ( char * )&h[-filter_length / 2] );
}


void
f444_420( videoinfo info, YUVimage * frame )
{
  int i, j, chor, cver, chor_half, cver_half;
  int filter_length, half_length, new_length;
  float *h = NULL, *line = NULL, *outline = NULL;
  float *U = NULL, *V = NULL;
  float *extension;

  filter_length = 7;
  h = ( float * )getarray( filter_length, sizeof( float ), "h" );
  h = &h[filter_length / 2];
  h[-3] = h[3] = -29 / 256.;
  h[-2] = h[2] = 0.;
  h[-1] = h[1] = 88 / 256.;
  h[0] = 138 / 256.;

  if( info.cwidth % 2 || info.cheight % 2 ) {
    printf( "can not handle this case!\n" );
    exit( 1 );
  }
  chor = info.cwidth;
  chor_half = chor / 2;
  cver = info.cheight;
  cver_half = cver / 2;


  U = ( float * )getarray( chor_half * cver, sizeof( float ), "U" );
  V = ( float * )getarray( chor_half * cver, sizeof( float ), "V" );

  outline =
    ( float * )getarray( MY_MAX( chor, cver ), sizeof( float ), "outline" );
  half_length = ( filter_length - 1 ) / 2;
  new_length = MY_MAX( chor, cver ) + 2 * half_length;
  extension = ( float * )getarray( new_length, sizeof( float ), "extension" );
  extension = &extension[half_length];



  for( i = 0; i < cver; i++ ) {
    line = &( frame->U[i * chor] );
    line_convolve( line, extension, chor, h, filter_length, ANAL, outline );
    for( j = 0; j < chor_half; j++ )
      U[i * chor_half + j] = outline[j * 2];
  }

  line = ( float * )getarray( cver, sizeof( float ), "line" );
  for( j = 0; j < chor_half; j++ ) {
    for( i = 0; i < cver; i++ )
      line[i] = U[i * chor_half + j];
    line_convolve( line, extension, cver, h, filter_length, ANAL, outline );
    for( i = 0; i < cver_half; i++ )
      U[i * chor_half + j] = outline[i * 2];
  }
  free( line );



  for( i = 0; i < cver; i++ ) {
    line = &( frame->V[i * chor] );
    line_convolve( line, extension, chor, h, filter_length, ANAL, outline );
    for( j = 0; j < chor_half; j++ )
      V[i * chor_half + j] = outline[j * 2];

  }

  line = ( float * )getarray( cver, sizeof( float ), "line" );
  for( j = 0; j < chor_half; j++ ) {
    for( i = 0; i < cver; i++ )
      line[i] = V[i * chor_half + j];
    line_convolve( line, extension, cver, h, filter_length, ANAL, outline );
    for( i = 0; i < cver_half; i++ )
      V[i * chor_half + j] = outline[i * 2];
  }
  free( line );

  free( frame->U );
  free( frame->V );
  free( outline );
  free( ( char * )&extension[-half_length] );

  frame->U = U;
  frame->V = V;


  free( ( char * )&h[-filter_length / 2] );
}

void
f422_444( videoinfo info, YUVimage * frame )
{
  int i, j, chor, cver, chor2;
  int filter_length, half_length, new_length;
  float *h = NULL, *line = NULL, *extension = NULL, *outline = NULL;
  float *U = NULL, *V = NULL;

  filter_length = 7;
  h = ( float * )getarray( filter_length, sizeof( float ), "h" );
  h = &h[filter_length / 2];
  h[-3] = h[3] = -12 / 256.;
  h[-2] = h[2] = 0.;
  h[-1] = h[1] = 140 / 256.;
  h[0] = 256 / 256.;

  chor = info.cwidth;
  chor2 = chor * 2;
  cver = info.cheight;

  line = ( float * )getarray( chor2, sizeof( float ), "line" );
  outline = ( float * )getarray( chor2, sizeof( float ), "outline" );

  half_length = ( filter_length - 1 ) / 2;
  new_length = chor2 + 2 * half_length;
  extension = ( float * )getarray( new_length, sizeof( float ), "extension" );
  extension = &extension[half_length];

  U = ( float * )getarray( chor2 * cver, sizeof( float ), "U" );
  V = ( float * )getarray( chor2 * cver, sizeof( float ), "V" );


  for( i = 0; i < cver; i++ ) {
    for( j = 0; j < chor2; j++ ) {
      if( j % 2 )
        line[j] = 0.;
      else
        line[j] = frame->U[i * chor + j / 2];
    }
    line_convolve( line, extension, chor2, h, filter_length, ANAL, outline );

    for( j = 0; j < chor2; j++ )
      U[i * chor2 + j] = outline[j];

  }

  for( i = 0; i < cver; i++ ) {
    for( j = 0; j < chor2; j++ ) {
      if( j % 2 )
        line[j] = 0;
      else
        line[j] = frame->V[i * chor + j / 2];
    }
    line_convolve( line, extension, chor2, h, filter_length, ANAL, outline );

    for( j = 0; j < chor2; j++ )
      V[i * chor2 + j] = outline[j];

  }

  free( frame->U );
  free( frame->V );

  frame->U = U;
  frame->V = V;
  free( outline );


  free( line );
  free( ( char * )&extension[-half_length] );
  free( ( char * )&h[-filter_length / 2] );
}



void
f420_422( videoinfo info, YUVimage * frame )
{
  int i, j, chor, cver, cver2;
  float *U = NULL, *V = NULL;

  chor = info.cwidth;
  cver = info.cheight;
  cver2 = cver * 2;

  U = ( float * )getarray( chor * cver2, sizeof( float ), "U" );
  V = ( float * )getarray( chor * cver2, sizeof( float ), "V" );

  // transpose U and V
  for( i = 0; i < cver; i++ ) {
    for( j = 0; j < chor; j++ ) {
      U[j * cver + i] = frame->U[i * chor + j];
      V[j * cver + i] = frame->V[i * chor + j];
    }
  }
  info.cwidth = cver;
  info.cheight = chor;

  free( frame->U );
  free( frame->V );

  frame->U = U;
  frame->V = V;

  f422_444( info, frame );


  U = ( float * )getarray( chor * cver2, sizeof( float ), "U" );
  V = ( float * )getarray( chor * cver2, sizeof( float ), "V" );

  // transpose U and V
  for( i = 0; i < info.cheight; i++ ) {
    for( j = 0; j < info.cwidth * 2; j++ ) {
      U[j * info.cheight + i] = frame->U[i * info.cwidth * 2 + j];
      V[j * info.cheight + i] = frame->V[i * info.cwidth * 2 + j];
    }
  }

  free( frame->U );
  free( frame->V );

  frame->U = U;
  frame->V = V;
}




void
f420_444( videoinfo info, YUVimage * frame )
{
  //upsample vertically
  f420_422( info, frame );

  info.cheight *= 2;
  //upsample horizontally
  f422_444( info, frame );
}
