#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "basic.h"
#include "rasterfile.h"
#include "structN.h"
#include "dpx.h"
#include "unix_pc.h"
#include "util_filtering.h"
#define EXTERN extern
#include "coderN.h"

struct rasterfile make_header( int hor, int ver, int depth );

/*****************************************************************************/
/*                               make_header()                               */
/*****************************************************************************/
struct rasterfile
make_header( int hor, int ver, int depth )
{
  struct rasterfile header;

  header.ras_magic = RAS_MAGIC;
  header.ras_width = hor;
  header.ras_height = ver;
  header.ras_depth = depth;
  header.ras_type = RT_STANDARD;

  if( depth == GRAYDATA ) {
    header.ras_length = hor * ver;
    header.ras_maptype = RMT_EQUAL_RGB;
    header.ras_maplength = 256 * 3;
  } else if( depth == RGBDATA ) {
    header.ras_length = 3 * hor * ver;
    header.ras_maptype = 0;
    header.ras_maplength = 0;
    /*header.ras_type = RT_FORMAT_RGB; */
  }
  return header;
}

/*****************************************************************************/
/*                              write_ras()                                  */
/*****************************************************************************/
void
write_ras( char *filename, struct rasterfile header, unsigned char *frame )
{
  FILE *fp;
  int i, ver, hor;
  int write;
  unsigned char map[256];
  int length;
  U32 *p;
  struct rasterfile temphead;

  /* open the file */
  if( ( fp = fopen( filename, "wb" ) ) == NULL ) {
    printf( "write_ras: can't open %s for writing\n", filename );
    exit( 1 );
  }

  /* write the file header */
  temphead = header;
  p = ( U32 * ) & temphead;
  for( i = 0; i < 8; i++ ) {
    p[i] = exchange_4byte_order( p[i] );
  }
  /*if((write = fwrite(header, sizeof(struct rasterfile), 1, fp)) != 1) {
     printf("write_ras: can't write header to %s\n", filename);
     exit(1);
     } */

  if( ( write = fwrite( p, sizeof( U32 ), 8, fp ) ) != 8 ) {
    printf( "write_ras: can't write header to %s\n", filename );
    exit( 1 );
  }

  /* write the color map */
  if( header.ras_maplength ) {
    for( i = 0; i < 256; i++ )
      map[i] = ( unsigned char )i;
    for( i = 0; i < 3; i++ ) {
      if( fwrite( map, sizeof( unsigned char ), 256, fp ) != 256 ) {
        printf( "write_ras: can't write colormap to %s(gray)\n", filename );
        exit( 1 );
      }
    }
  }

  /* write the image data */
  ver = header.ras_height;
  hor = header.ras_width;
  if( ( header.ras_type == RT_STANDARD )
      || ( header.ras_type == RT_FORMAT_RGB ) ) {
    if( header.ras_depth == GRAYDATA )
      length = hor * ver;
    else if( header.ras_depth == RGBDATA )
      length = 3 * hor * ver;

    if( ( write =
          fwrite( frame, sizeof( unsigned char ), length, fp ) ) != length ) {
      printf( "write_ras: can't write data of %s write %d length %d\n",
              filename, write, length );
      exit( 1 );
    }
  } else {
    printf( "write_ras: input file is not standard file\n" );
    exit( 1 );
  }

  /* file close */
  fclose( fp );
}

/*****************************************************************************/
/*                              read_ras()                                   */
/*****************************************************************************/
struct rasterfile *
read_ras( char *filename, unsigned char **frame )
{
  FILE *fp;
  int i, ver, hor;
  unsigned char map[256], index;
  int length;
  U32 *p;
  struct rasterfile *header;

  /* open the file */
  if( ( fp = fopen( filename, "rb" ) ) == NULL ) {
    printf( "read_ras: can't open %s for reading\n", filename );
    exit( 1 );
  }

  p = ( U32 * ) getarray( 8, sizeof( U32 ), "p" );
  /* read the file header */
  if( fread( p, sizeof( U32 ), 8, fp ) != 8 ) {
    printf( "read_ras: can't read header to %s (read_ras)\n", filename );
    exit( 1 );
  }
  for( i = 0; i < 8; i++ ) {
    p[i] = exchange_4byte_order( p[i] );
  }
  header = ( struct rasterfile * )p;

  //printf("maplength %d\n", header->ras_maplength);
  /* read the color map */
  if( header->ras_maplength ) {
    for( i = 0; i < header->ras_maplength / 3; i++ )
      map[i] = ( unsigned char )i;
    for( i = 0; i < 3; i++ ) {
      if( fread( map, sizeof( unsigned char ), header->ras_maplength / 3, fp )
          != ( unsigned )header->ras_maplength / 3 ) {
        printf( "read_ras: can't read colormap to %s(gray)\n", filename );
        exit( 1 );
      }
    }
  }

  /* read the image data */
  ver = header->ras_height;
  hor = header->ras_width;
  //if((header->ras_type == RT_STANDARD) || (header->ras_type == RT_FORMAT_RGB)){
  if( 1 ) {
    if( header->ras_depth == GRAYDATA )
      length = hor * ver;
    else if( header->ras_depth == RGBDATA )
      length = 3 * hor * ver;

    *frame =
      ( unsigned char * )getarray( length, sizeof( unsigned char ),
                                   "*frame" );
    if( fread( *frame, sizeof( unsigned char ), length, fp ) !=
        ( unsigned )length ) {
      printf( "read_ras: can't read data of %s(read_ras)\n", filename );
      exit( 1 );
    }
    if( header->ras_depth == GRAYDATA ) {       //printf("gray(ras_util.c)");
      for( i = 0; i < length; i++ ) {
        index = ( *frame )[i];
        ( *frame )[i] = map[index];
      }
    }
  } else {
    printf( "read_ras: input file is not standard file\n" );
    exit( 1 );
  }

  /* file close */
  fclose( fp );

  return ( header );
}


struct rasterfile *
read_charimg24( char *name, unsigned char **img )
{
  struct rasterfile *rhead = NULL;
  FILE *fp;
  int nn_3, i, tt;
  unsigned char *buffer = NULL;
  U32 *p;



  /* open file stream to read */
  if( !name ) {
    fp = stdin;
  } else if( ( fp = fopen( name, "rb" ) ) == NULL ) {
    fprintf( stderr, "read_charimg24: can't open %s for read\n", name );
    exit( 1 );
  }


  /* read in ras header */
  rhead =
    ( struct rasterfile * )getarray( 1, sizeof( struct rasterfile ),
                                     "rhead" );
  if( fread( ( char * )rhead, sizeof( struct rasterfile ), 1, fp ) != 1 ) {
    fprintf( stderr, "read_charimg: can't read ras header of %s\n", name );
    exit( 1 );
  }
  p = ( U32 * ) rhead;
  for( i = 0; i < 8; i++ ) {
    p[i] = exchange_4byte_order( p[i] );
  }
  rhead = ( struct rasterfile * )p;

  if( rhead->ras_type != RT_STANDARD || rhead->ras_depth != 24 ) {
    fprintf( stderr, "read_charimg24: Mismatch image read types\n" );
    exit( -1 );
  }
  if( rhead->ras_maplength != 0 ) {
    fseek( fp, ( long )rhead->ras_maplength, 1 );
  }
  nn_3 = 3 * rhead->ras_width * rhead->ras_height;


  /* allocate memory for image */
  buffer =
    ( unsigned char * )getarray( nn_3, sizeof( unsigned char ), "buffer" );


  /* read in image into buffer */
  if( ( tt = fread( buffer, sizeof( unsigned char ), nn_3, fp ) ) != nn_3 ) {
    fprintf( stderr, "Error reading %s %d %d \n", name, tt, nn_3 );
    exit( -1 );
  }
  fclose( fp );


  *img = buffer;
  return ( rhead );
}
