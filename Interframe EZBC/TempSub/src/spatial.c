#include "stdio.h"
#include "stdlib.h"
#include "basic.h"
#include "structN.h"
#include "Choisubband.h"
#include <assert.h>


  /*************************************************/
  /*            in              out                */
  /*        -----------     -----------            */
  /*        |         |     |    |    |            */
  /*        |         |     | LL | HL |            */
  /*        |         |     -----------            */
  /*        |         |     |    |    |            */
  /*        |         |     | LH | HH |            */
  /*        -----------     -----------            */
  /*************************************************/
/*
 * spatial_anal()
 * subband decomposition 
 * full: input image for analysis
 * (hor, ver): size of the full
 */

void
spatial_anal( float *full, int hor, int ver, float *ll, float *lh, float *hl,
              float *hh, int filterType )
{
  int i, k, l, pos, ppos, sx, sy;
  int half_hor, half_ver;
  float *image, *subband[4];    // image[comp], subband[comp][sub]

  if( ( hor % 2 != 0 ) || ( ver % 2 != 0 ) ) {
    printf( "can not handle this case(encoder.c)\n" );
    exit( 1 );
  }

  half_hor = hor / 2;
  half_ver = ver / 2;

  image = ( float * )getarray( hor * ver, sizeof( float ), "image" );
  for( i = 0; i < hor * ver; i++ )
    image[i] = full[i];

  subband[0] = ll;
  subband[1] = lh;
  subband[2] = hl;
  subband[3] = hh;

  analysis( image, 0, 0, hor, ver, hor, ver, filterType );      // subband.c

  for( i = 0; i < 4; i++ ) {
    sx = i / 2;
    sy = i % 2;
    //printf("sx = %d, sy = %d (encoder.c)\n", sx, sy);
    for( k = 0; k < half_ver; k++ ) {
      pos = k * half_hor;
      ppos = ( sy * half_ver + k ) * hor + sx * half_hor;

      for( l = 0; l < half_hor; l++ ) {
        subband[i][pos + l] = image[ppos + l];
      }
    }
  }
  free( image );
}

void
spatial_anal_frame( YUVimage_ptr fr, videoinfo info, YUVimage_ptr low,
                    YUVimage_ptr high0, YUVimage_ptr high1,
                    YUVimage_ptr high2, int filterType )
{

  spatial_anal( fr->Y, info.ywidth, info.yheight, low->Y,
                high0->Y, high1->Y, high2->Y, filterType );


  if( info.cwidth && info.cheight ) {
    spatial_anal( fr->U, info.cwidth, info.cheight, low->U,
                  high0->U, high1->U, high2->U, filterType );

    spatial_anal( fr->V, info.cwidth, info.cheight, low->V,
                  high0->V, high1->V, high2->V, filterType );
  }

}

/*
 * spatial_syn()
 * subband synthesis
 * full: result of synthesis
 * (hor, ver): size of the full
 */
void
spatial_syn( float *full, int hor, int ver, float *ll, float *lh, float *hl,
             float *hh, int filterType )
{
  int i, k, l, pos, ppos, sx, sy;
  int half_hor, half_ver;
  float *subband[4];            // image[comp], subband[comp][sub]

  if( ( hor % 2 != 0 ) || ( ver % 2 != 0 ) ) {
    printf( "can not handle this case(encoder.c)\n" );
    exit( 1 );
  }

  half_hor = hor / 2;
  half_ver = ver / 2;


  subband[0] = ll;
  subband[1] = lh;
  subband[2] = hl;
  subband[3] = hh;

  for( i = 0; i < 4; i++ ) {
    sx = i / 2;
    sy = i % 2;
    //printf("sx = %d, sy = %d (encoder.c)\n", sx, sy);
    for( k = 0; k < half_ver; k++ ) {
      pos = k * half_hor;
      ppos = ( sy * half_ver + k ) * hor + sx * half_hor;

      for( l = 0; l < half_hor; l++ ) {
        full[ppos + l] = subband[i][pos + l];
      }
    }
  }

  synthesis( full, 0, 0, hor, ver, hor, ver, filterType );      // subband.c
}


void
spatial_syn_frame( YUVimage_ptr fr, videoinfo info, YUVimage_ptr low,
                   YUVimage_ptr high0, YUVimage_ptr high1, YUVimage_ptr high2,
                   int filterType )
{
  spatial_syn( fr->Y, info.ywidth, info.yheight, low->Y,
               high0->Y, high1->Y, high2->Y, filterType );

  if( info.cwidth && info.cheight ) {
    spatial_syn( fr->U, info.cwidth, info.cheight, low->U,
                 high0->U, high1->U, high2->U, filterType );

    spatial_syn( fr->V, info.cwidth, info.cheight, low->V,
                 high0->V, high1->V, high2->V, filterType );
  }
}
