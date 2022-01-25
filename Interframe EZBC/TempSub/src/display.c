#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#define EXTERN extern
#include "rasterfile.h"
#include "basic.h"
#include "structN.h"
#include "coderN.h"
#include "bmeN.h"
#include "memoryN.h"
#include "miscN.h"
#include "dpx.h"
#include "unix_pc.h"
#include "ioN.h"
#include "analsyn.h"
#include "util_filtering.h"
#include <assert.h>

#define OUTPUT_RAS

void copyframe( YUVimage * source, YUVimage * dest, videoinfo info );
void snr_frame( float *ysnr, float *usnr, float *vsnr, YUVimage_ptr codeframe,
                YUVimage_ptr inframe, videoinfo info );
void write_frame( YUVimage frame, videoinfo info, char *inname, int index,
                  enum FORMAT format );
void yuv2RGB( unsigned char *RGBframe, int yhor, int yver, float *Yframe,
              int chor, int cver, float *Uframe, float *Vframe );
struct rasterfile *exchange_rasheader( struct rasterfile *rhead );
struct rasterfile make_header( int hor, int ver, int depth );
void write_ras( char *filename, struct rasterfile header,
                unsigned char *frame );


void
writefloatimg4( char *name, struct rasterfile *rhead, float *data )
{
  FILE *fp;
  int i, nn;
  U32 *p;
  struct rasterfile temphead;

  rhead->ras_depth = 32;
  rhead->ras_length = rhead->ras_height * rhead->ras_width * sizeof( float );
  rhead->ras_maptype = RMT_NONE;
  rhead->ras_maplength = 0;

  if( !name ) {
    fp = stdout;
  } else if( ( fp = fopen( name, "wb" ) ) == NULL ) {
    ( void )fprintf( stderr, "write_floatimg4: can't open %s for write\n",
                     name );
    exit( -1 );
  }
  nn = rhead->ras_width * rhead->ras_height;


  /* write the file header */
  temphead = *rhead;
  p = ( U32 * ) & temphead;
  for( i = 0; i < 8; i++ ) {
    p[i] = exchange_4byte_order( p[i] );
  }
  if( fwrite( p, sizeof( U32 ), 8, fp ) != 8 ) {
    printf( "write_ras: can't write header to %s\n", name );
    exit( 1 );
  }

  p = ( U32 * ) data;
  for( i = 0; i < nn; i++ ) {
    p[i] = exchange_4byte_order( p[i] );
  }

  if( fwrite( p, sizeof( U32 ), nn, fp ) != ( unsigned )nn ) {
    printf( "write_ras: can't write header to %s\n", name );
    exit( 1 );
  }

  ( void )fclose( fp );
  return;
}


void
block2pixel( float *mvx, float *mvy, vector_ptr fmv, int cx, int cy, int xblk,
             int yblk, int hor, int ver )
{
  int i, j, xblock, yblock, pos;

  /* change the structure of motion vectors */
  /* from the block-based to the pixel-based */
  /* write the motion vector of the block recursively */

  if( fmv->child ) {
    block2pixel( mvx, mvy, fmv->child0, cx, cy, xblk / 2, yblk / 2, hor,
                 ver );
    block2pixel( mvx, mvy, fmv->child1, cx + xblk / 2, cy, xblk / 2, yblk / 2,
                 hor, ver );
    block2pixel( mvx, mvy, fmv->child2, cx, cy + yblk / 2, xblk / 2, yblk / 2,
                 hor, ver );
    block2pixel( mvx, mvy, fmv->child3, cx + xblk / 2, cy + yblk / 2,
                 xblk / 2, yblk / 2, hor, ver );
  } else {
    /* consider the small block around the boundaries */
    xblock = ( cx + xblk <= hor ) ? xblk : hor - cx;
    yblock = ( cy + yblk <= ver ) ? yblk : ver - cy;

    if( xblock <= 0 || yblock <= 0 ) {
      /*      printf("xblock<=0 || yblock<=0 in block2pixel2() !\n");*/
      return;
    }

    for( i = cy; i < cy + yblock; i++ ) {
      for( j = cx; j < cx + xblock; j++ ) {
        pos = i * hor + j;
        mvx[pos] = fmv->mvx;
        mvy[pos] = fmv->mvy;
      }
    }
  }
}



void write_highband(YUVimage frame, int GOPIndex, int level, int index, 
                    videoinfo info)
{
  int ysize, pos;
  float weight, *data; 
  char name[512];
  unsigned char *char_data = NULL;
  struct rasterfile rhead;

  ysize = info.ywidth * info.yheight;

  temporal_filter();
  weight = (float) pow( HPW1[1], (float) level);
  
  data = 
    (float *) getarray ( ysize, sizeof(float), "data_write_highband" );
  char_data =
    (unsigned char *) getarray ( ysize, sizeof(unsigned char), "char_data" );
  
  for( pos = 0; pos < ysize; pos++ ) {
    data[pos] = frame.Y[pos] / weight;
  }

  rhead = make_header(info.ywidth, info.yheight, GRAYDATA);
  cnv_data_4_1(data, char_data, ysize);

  sprintf( name, "%s_statistics/GOP%02d_%s%01d_fr%03d.ras",
           info.bitname, GOPIndex, "H", level, index );  
  write_ras(name, rhead, char_data);
  
  free(char_data); 
  free(data);
}


void
write_lowband( YUVimage frame, int GOPIndex, int level, int index,
               videoinfo info )
{
  int ysize, csize, pos;
  float weight; 
  char name[512];
  YUVimage temp;

  ysize = info.ywidth * info.yheight;
  csize = info.cwidth * info.cheight;

  temporal_filter();
  weight = (float) pow( LPW4[1], (float) level);

  frame_alloc( &temp, info );
  
  for( pos = 0; pos < ysize; pos++ ) {
    temp.Y[pos] = frame.Y[pos] / weight;
  }
  
  for( pos = 0; pos < csize; pos++ ) {
    temp.U[pos] = frame.U[pos] / weight;
    temp.V[pos] = frame.V[pos] / weight;
  }
 
#ifdef OUTPUT_RAS
  sprintf( name, "%s_statistics/GOP%02d_%s%01d_fr%03d.ras", 
           info.bitname, GOPIndex, "L", level, index );
  write_frame(temp, info, name, index, RAS);
#else 
  sprintf( name, "%s_statistics/GOP%02d_%s%01d_fr%03d.yuv", 
           info.bitname, GOPIndex, "L", level, index );
  write_frame(temp, info, name, index, YUV);
#endif
  
  free_frame( temp );
}

