#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "structN.h"
#include "rasterfile.h"
#include "basic.h"
#include "dpx.h"
#include "unix_pc.h"
#include "chrom.h"
#include "miscN.h"
#define EXTERN extern
#include "coderN.h"
#include "util_filtering.h"
extern enum FLAG **scene_change;
extern enum FLAG **dec_scene_change;

#define RGBYUV //do RGB to YUV conversion
//#define OUTPUT_RAS


/*****************************************************************************/
/*                                 clipping()                                */
/*****************************************************************************/
int
clipping( float ipix )
{
  int opix;

  if( ipix > 255.0 )
    opix = 255;
  else if( ipix < 0.0 )
    opix = 0;
  else
    opix = nint( ipix );

  return ( opix );
}

/*****************************************************************************/
/*                                 clipping2()                               */
/*****************************************************************************/
float
clipping2( float ipix )
{
  float opix;

  if( ipix > 255.0 )
    opix = 255.;
  else if( ipix < 0.0 )
    opix = 0.;
  else
    opix = ipix;

  return ( opix );
}


/*****************************************************************************/
/*                                 read_child()                              */
/*****************************************************************************/
void
read_child( vector_ptr fmv, FILE * fpio )
{
  /* if child is NULL, then there is nothing to do */
  /* else, allocate the memory */
  /*       and read the 4 child motion vectors */
  /*       and call again */

  if( fmv->child ) {
    fmv->child0 = ( vector_ptr ) calloc( 1, sizeof( vector ) );
    fmv->child1 = ( vector_ptr ) calloc( 1, sizeof( vector ) );
    fmv->child2 = ( vector_ptr ) calloc( 1, sizeof( vector ) );
    fmv->child3 = ( vector_ptr ) calloc( 1, sizeof( vector ) );

    fread( fmv->child0, sizeof( vector ), 1, fpio );
    fread( fmv->child1, sizeof( vector ), 1, fpio );
    fread( fmv->child2, sizeof( vector ), 1, fpio );
    fread( fmv->child3, sizeof( vector ), 1, fpio );

    fmv->child0->parent = fmv;
    fmv->child1->parent = fmv;
    fmv->child2->parent = fmv;
    fmv->child3->parent = fmv;

    read_child( fmv->child0, fpio );
    read_child( fmv->child1, fpio );
    read_child( fmv->child2, fpio );
    read_child( fmv->child3, fpio );
  } else {
  }
}


/*****************************************************************************/
/*                               write_child()                               */
/*****************************************************************************/
void
write_child( vector_ptr fmv, FILE * fpio )
{
  /* if child is NULL, then there is nothing to do */
  /* else, write the 4 child motion vectors and call again */

  if( fmv->child ) {
    fwrite( fmv->child0, sizeof( vector ), 1, fpio );
    fwrite( fmv->child1, sizeof( vector ), 1, fpio );
    fwrite( fmv->child2, sizeof( vector ), 1, fpio );
    fwrite( fmv->child3, sizeof( vector ), 1, fpio );

    write_child( fmv->child0, fpio );
    write_child( fmv->child1, fpio );
    write_child( fmv->child2, fpio );
    write_child( fmv->child3, fpio );
  } else {
  }
}


/*****************************************************************************/
/*                                   print_mv()                              */
/*****************************************************************************/
void
print_mv( vector_ptr fmv, int cx, int cy, int xblk, int yblk, int hor,
          int ver )
{
  float mvx, mvy, mad;

  /* this routine draws the variable block motion vectors */

  if( fmv->child ) {

    if( xblk >= 32 ) {
      mvx = fmv->mvx;
      mvy = fmv->mvy;
      mad = fmv->mad;
      printf( "x %d y %d xblk %d yblk %d mvx %.1f mvy %.1f mad %.1f\n", cx,
              cy, xblk, yblk, mvx, mvy, mad );
    }
    print_mv( fmv->child0, cx, cy, xblk / 2, yblk / 2, hor, ver );
    print_mv( fmv->child1, cx + xblk / 2, cy, xblk / 2, yblk / 2, hor, ver );
    print_mv( fmv->child2, cx, cy + yblk / 2, xblk / 2, yblk / 2, hor, ver );
    print_mv( fmv->child3, cx + xblk / 2, cy + yblk / 2, xblk / 2, yblk / 2,
              hor, ver );
  } else {
    if( cx < hor && cy < ver && xblk >= 32 ) {
      mvx = fmv->mvx;
      mvy = fmv->mvy;
      mad = fmv->mad;
      printf( "x %d y %d xblk %d yblk %d mvx %.1f mvy %.1f mad %.1f\n", cx, cy,
              xblk, yblk, mvx, mvy, mad );
    }

    return;
  }
}




/*****************************************************************************
 **                                                                         ** 
 **  ras2yuv                                                                **
 **                                                                         ** 
 **  -. convert the sun rasterfile to YUV file (YUV with 4:1:1 or 4:2:2)    **
 **                                                                         ** 
 **  author : Seung-Jong Choi                                               **
 **                                                                         ** 
 *****************************************************************************/

/* sub422=0; *//* 4:1:1 subsampling */
/*sub422=1;  *//* 4:2:2 subsampling */
void
ras2yuv( char *rasname, int yhor, int yver, float *Yfr, int chor, int cver,
         float *Ufr, float *Vfr )
{
  unsigned char *RGBfr, *temp, map[768];
  FILE *fpin;
  int i, j, hor, ver, mapsize, sub422;
  float Ydata, Udata, Vdata, Rdata, Gdata, Bdata;
  struct rasterfile *header;
  U32 *p;

  if( !( fpin = fopen( rasname, "rb" ) ) ) {
    printf( "ras2yuv: can't open rasterfile %s.\n", rasname );
    exit( 1 );
  }

  header =
    ( struct rasterfile * )getarray( 1, sizeof( struct rasterfile ),
                                     "header" );
  /* read the header of rasterfile */
  if( fread( header, sizeof( struct rasterfile ), 1, fpin ) != 1 ) {
    printf( "\n" );
    printf( "ras2yuv: can't read header from %s.\n", rasname );
    exit( 1 );
  }

  p = ( U32 * ) header;
  for( i = 0; i < 8; i++ ) {
    p[i] = exchange_4byte_order( p[i] );
  }

  header = ( struct rasterfile * )p;



  hor = header->ras_width;
  ver = header->ras_height;
  if( ( yhor != hor ) || ( yver != ver ) ) {
    printf
      ( "size of image error info.width %d, info.height %d, header->ras_width %d, header->ras_height %d (ioN.c)!\n",
        yhor, yver, hor, ver );
    header->ras_depth = 24;
    hor = yhor;
    ver = yver;
    /*exit(0); */
  }


  /************ allocate memory and read the data of rasterfile **********/
  /*
    if(sub422){ 
    RGBfr = (unsigned char*) calloc(ver*hor*3,   sizeof(unsigned char));
    Yfr   = (float*) calloc(hor*ver,     sizeof(float));
    Ufr   = (float*) calloc((hor/2)*ver, sizeof(float));
    Vfr   = (float*) calloc((hor/2)*ver, sizeof(float));
    }
    else {   
    RGBfr = (unsigned char*) calloc(ver*hor*3,  sizeof(unsigned char));
    Yfr   = (float*) calloc(hor*ver,    sizeof(float));
    Ufr   = (float*) calloc(hor*ver/4,  sizeof(float));
    Vfr   = (float*) calloc(hor*ver/4,  sizeof(float));
    }
  */

  /* read the map */
  mapsize = header->ras_maplength;
  if( ( mapsize != 768 ) && ( mapsize != 0 ) ) {
    printf( "the size of maplength %d is not 768(ioN.c)\n", mapsize );
    exit( 0 );
  }

  if( mapsize ) {
    if( fread( map, sizeof( unsigned char ), mapsize, fpin ) !=
        ( unsigned )mapsize ) {
      printf( "ras2yuv: error in reading of map from %s.\n", rasname );
      exit( 1 );
    }
  }

  /********************* read the data *****************/
  if( header->ras_depth == GRAYDATA ) {
    if( ( chor != 0 ) || ( cver != 0 ) ) {
      printf( " the size of chrom component error(ras2yuv.c)!\n" );
      exit( 0 );
    }
    temp = ( unsigned char * )calloc( hor * ver, sizeof( unsigned char ) );
    if( fread( temp, sizeof( unsigned char ), hor * ver, fpin ) !=
        ( unsigned )( hor * ver ) ) {
      printf( "ras2yuv: can't read the Y data from %s.\n", rasname );
      exit( 1 );
    }
    for( i = 0; i < hor * ver; i++ )
      Yfr[i] = ( float )temp[i];
    free( temp );
  } else if( header->ras_depth == RGBDATA ) {

    if( ( yhor == 2 * chor ) && ( yver == 2 * cver ) )
      sub422 = 0;
    else if( ( yhor == 2 * chor ) && ( yver == cver ) )
      sub422 = 1;
    else if( ( yhor == chor ) && ( yver == cver ) )
      sub422 = 2;
    else {
      printf( "size of image error(ras2yuv.c)\n" );
      exit( 0 );
    }

    RGBfr =
      ( unsigned char * )calloc( ver * hor * 3, sizeof( unsigned char ) );
    if( fread( RGBfr, sizeof( unsigned char ), 3 * hor * ver, fpin ) !=
        ( unsigned )( 3 * hor * ver ) ) {
      printf( "ras2yuv: can't read the RGB data from %s.\n", rasname );
      exit( 1 );
    }
    for( i = 0; i < ver; i++ ) {
      for( j = 0; j < hor; j++ ) {
        /* read the data and change to the float */
        Bdata = ( float )RGBfr[i * ( 3 * hor ) + 3 * j];
        Gdata = ( float )RGBfr[i * ( 3 * hor ) + 3 * j + 1];
        Rdata = ( float )RGBfr[i * ( 3 * hor ) + 3 * j + 2];

        RGB2YUV( Rdata, Gdata, Bdata, &Ydata, &Udata, &Vdata );
        Udata += 128.;
        Vdata += 128.;

        if( sub422 == 2 ) {
          Yfr[i * hor + j] = Ydata;
          Ufr[i * hor + j] = Udata;
          Vfr[i * hor + j] = Vdata;
        } else if( sub422 == 1 ) {
          Yfr[i * hor + j] = Ydata;
          if( !( j % 2 ) ) {
            Ufr[i * ( hor / 2 ) + ( j / 2 )] = Udata;
            Vfr[i * ( hor / 2 ) + ( j / 2 )] = Vdata;
          }
        } else {
          Yfr[i * hor + j] = Ydata;
          if( !( j % 2 ) && !( i % 2 ) ) {
            Ufr[( i / 2 ) * ( hor / 2 ) + j / 2] = Udata;
            Vfr[( i / 2 ) * ( hor / 2 ) + j / 2] = Vdata;
          }
        }
      }                         /* for j */
    }                           /* for i */
    free( RGBfr );

  } else {
    printf( "ras2yuv: %s is not supported.\n", rasname );
    exit( 1 );
  }


  fclose( fpin );

  free( header );

}



/*****************************************************************************
 **                                                                         ** 
 **  yuv2ras                                                                **
 **                                                                         ** 
 **  -. convert the YUV file to sun rasterfile (YUV with 4:1:1 or 4:2:2)    **
 **                                                                         ** 
 **  author : Seung-Jong Choi                                               **
 **                                                                         ** 
 *****************************************************************************/

void
yuv2ras( char *rasname, int yhor, int yver, float *Yframe, int chor, int cver,
         float *Uframe, float *Vframe )
{
  int i, j, ci, cj;
  float Ydata, Udata, Vdata, Rdata, Gdata, Bdata;
  struct rasterfile header;
  unsigned char *RGBframe, *temp;

  ci = ( cver ) ? yver / cver : 0;
  cj = ( chor ) ? yhor / chor : 0;

  /* allocate memory */

  /*
    Yframe    = (unsigned char*) calloc(yhor*yver,   sizeof(unsigned char));
    Uframe    = (unsigned char*) calloc(chor*cver,   sizeof(unsigned char));
    Vframe    = (unsigned char*) calloc(chor*cver,   sizeof(unsigned char));
  */

  /* read the data from SIF or CCIR */
  /*
    if(!(fpin  = fopen(yuvname, "rb"))){
    printf("yuv2ras: can't open %s for reading.\n", yuvname);
    exit(1);
    }

    if(fread(Yframe, sizeof(unsigned char), yhor*yver, fpin) != yhor*yver){
    printf("yuv2ras: can't read Y from %s\n", yuvname);
    exit(1);
    }
    if(fread(Uframe, sizeof(unsigned char), chor*cver, fpin) != chor*cver){
    printf("yuv2ras: can't read U from %s\n", yuvname);
    exit(1);
    }
    if(fread(Vframe, sizeof(unsigned char), chor*cver, fpin) != chor*cver){
    printf("yuv2ras: can't read V from %s\n", yuvname);
    exit(1);
    }
    fclose(fpin);
  */

  /* convert YUV to BGR */
  if( chor * cver ) {
    RGBframe =
      ( unsigned char * )calloc( yhor * yver * 3, sizeof( unsigned char ) );
    for( i = 0; i < yver; i++ ) {
      for( j = 0; j < yhor; j++ ) {
        Ydata = clipping2( Yframe[i * yhor + j] );
        Udata = clipping2( Uframe[( i / ci ) * chor + ( j / cj )] );
        Vdata = clipping2( Vframe[( i / ci ) * chor + ( j / cj )] );

        /* YUV -> RGB */
        Udata -= 128.;
        Vdata -= 128.;
        YUV2RGB( Ydata, Udata, Vdata, &Rdata, &Gdata, &Bdata );

        /* clipping */
        if( Rdata >= 255. )
          Rdata = 255.;
        if( Rdata <= 0. )
          Rdata = 0.;
        if( Gdata >= 255. )
          Gdata = 255.;
        if( Gdata <= 0. )
          Gdata = 0.;
        if( Bdata >= 255. )
          Bdata = 255.;
        if( Bdata <= 0. )
          Bdata = 0.;

        /* formatting */
        RGBframe[i * 3 * yhor + 3 * j] = ( unsigned char )nint( Bdata );
        RGBframe[i * 3 * yhor + 3 * j + 1] = ( unsigned char )nint( Gdata );
        RGBframe[i * 3 * yhor + 3 * j + 2] = ( unsigned char )nint( Rdata );
      }
    }
    header = make_header( yhor, yver, RGBDATA );

    write_ras( rasname, header, RGBframe );
    free( RGBframe );
  } else {                      /* gray file */
    temp = ( unsigned char * )calloc( yhor * yver, sizeof( unsigned char ) );
    header = make_header( yhor, yver, GRAYDATA );
    for( i = 0; i < yhor * yver; i++ ) {
      if( Yframe[i] >= 255. )
        temp[i] = 255;
      else if( Yframe[i] <= 0. )
        temp[i] = 0;
      else
        temp[i] = ( unsigned char )Yframe[i];
    }
    write_ras( rasname, header, temp );
    free( temp );
  }



}




/*****************************************************************************
 **                                                                         ** 
 **  yuv2RGB                                                                **
 **                                                                         ** 
 **  -. convert the YUV file to RGB  (YUV with 4:1:1 or 4:2:2)              **
 **                                                                         ** 
 **                                                                         **
 **                                                                         ** 
 *****************************************************************************/

void
yuv2RGB( unsigned char *RGBframe, int yhor, int yver, float *Yframe, int chor,
         int cver, float *Uframe, float *Vframe )
{
  int i, j, ci, cj;
  float Ydata, Udata, Vdata, Rdata, Gdata, Bdata;

  ci = ( cver ) ? yver / cver : 0;
  cj = ( chor ) ? yhor / chor : 0;

  /* allocate memory */

  /* convert YUV to BGR */
  if( chor * cver ) {
    /*      RGBframe  = (unsigned char*) calloc(yhor*yver*3, sizeof(unsigned char));*/
    for( i = 0; i < yver; i++ ) {
      for( j = 0; j < yhor; j++ ) {
        Ydata = clipping2( Yframe[i * yhor + j] );      /* QUEST */
        Udata = clipping2( Uframe[( i / ci ) * chor + ( j / cj )] );
        Vdata = clipping2( Vframe[( i / ci ) * chor + ( j / cj )] );

        /* YUV -> RGB */
        Udata -= 128.;
        Vdata -= 128.;
        YUV2RGB( Ydata, Udata, Vdata, &Rdata, &Gdata, &Bdata );

        /* clipping */
        if( Rdata >= 255. )
          Rdata = 255.;
        if( Rdata <= 0. )
          Rdata = 0.;
        if( Gdata >= 255. )
          Gdata = 255.;
        if( Gdata <= 0. )
          Gdata = 0.;
        if( Bdata >= 255. )
          Bdata = 255.;
        if( Bdata <= 0. )
          Bdata = 0.;

        /* formatting */
        RGBframe[i * 3 * yhor + 3 * j] = ( unsigned char )nint( Bdata );
        RGBframe[i * 3 * yhor + 3 * j + 1] = ( unsigned char )nint( Gdata );
        RGBframe[i * 3 * yhor + 3 * j + 2] = ( unsigned char )nint( Rdata );
      }
    }

  } else {                      /* gray file */

    for( i = 0; i < yhor * yver; i++ ) {
      if( Yframe[i] >= 255. )
        RGBframe[i] = 255;
      else if( Yframe[i] <= 0. )
        RGBframe[i] = 0;
      else
        RGBframe[i] = ( unsigned char )Yframe[i];
    }
  }

}


/*****************************************************************************/
/*                              read_frame()                                 */
/*****************************************************************************/
void
read_frame( YUVimage_ptr frame, videoinfo info, char *inname, int index,
            enum FORMAT format )
{
  int i, size, ret;
#ifdef SUPPORT_DPX
  int j, row, col, pos;
  float Rdata, Gdata, Bdata;
  U16 *RGBframe;
  struct dpx *dpxheader;
#endif
  int frameskip = 0; 
  _int64 offset; 
  unsigned char *temp;
  char frame_name[250];
  //float *LRGB;
  FILE *fpio;

  offset = ( frameskip + 1 ) * index * ( info.ywidth*info.yheight + 2 * info.cwidth*info.cheight ); 
//  printf("\noffset = %ld\n",offset);

  switch ( format ) {
  case RAS:
    sprintf( frame_name, inname, index );
    ras2yuv( frame_name, info.ywidth, info.yheight, frame->Y, info.cwidth,
             info.cheight, frame->U, frame->V );
    break;
  case YUV:
    /* assume that the order of data file is Y[total], U[total], V[total] */
  
    if ((strchr (inname, '%'))!=0) //check the filename for frame or video mode 
      {
        sprintf(frame_name, inname, index);   // for frame-mode only
        if (!(fpio = fopen(frame_name, "rb")))
          {
            printf("read_frame: %s\n", frame_name);
            exit(1);
          }
      }
    else
      { // for container-mode (video)
        sprintf(frame_name, inname); 
        if (!(fpio = fopen(inname, "rb")))
          {
            printf("read_Video: %s\n", inname);
            exit(1);
          }
//        fseek (fpio, offset, SEEK_SET); //seek to the picture in the file
		_fseeki64(fpio, offset, SEEK_SET);
        /* if ( fseek (fpio, offset, SEEK_SET)); //seek to the picture in the file with error check
           {
           printf("Error while seeking to picture");
           exit (1);  
           } */
      }

    size = ( info.ywidth ) * ( info.yheight );
    temp =
      ( unsigned char * )getarray( size, sizeof( unsigned char ), "temp" );
    if( ( ret = fread( temp, sizeof( unsigned char ), size, fpio ) ) != size ) {

      printf
        ( "error for reading file %s.Y in read_frame(): return size %d, target size %d\n",
          frame_name, ret, size );
      exit( 1 );
    }
    for( i = 0; i < size; i++ )
      frame->Y[i] = ( float )temp[i];

    size = ( info.cwidth ) * ( info.cheight );
    if( fread( temp, sizeof( unsigned char ), size, fpio ) !=
        ( unsigned )( size ) ) {
      printf( "error for reading file %s.U in read_frame()\n", frame_name );
      exit( 1 );
    }
    for( i = 0; i < size; i++ )
      frame->U[i] = ( float )temp[i];

    if( fread( temp, sizeof( unsigned char ), size, fpio ) !=
        ( unsigned )( size ) ) {
      printf( "error for reading file %s.V in read_frame()\n", frame_name );
      exit( 1 );
    }
    for( i = 0; i < size; i++ )
      frame->V[i] = ( float )temp[i];

    free( temp );
    fclose( fpio );
    break;
#ifdef SUPPORT_DPX
  case DPX:
    //printf("index %d %d\n", index, 0);
    sprintf( frame_name, inname, index );
    RGBframe = read_dpx( frame_name, &dpxheader, &row, &col );

    if( ( row != info.yheight ) || ( col != info.ywidth ) ) {
      printf( " the size of image does not match with the info\n" );
      exit( 1 );
    }


    for( i = 0; i < info.yheight; i++ ) {
      for( j = 0; j < info.ywidth; j++ ) {
        Bdata = RGBframe[i * ( 3 * info.ywidth ) + 3 * j];
        Gdata = RGBframe[i * ( 3 * info.ywidth ) + 3 * j + 1];
        Rdata = RGBframe[i * ( 3 * info.ywidth ) + 3 * j + 2];

        if( info.coding_domain != LOG ) {
          Bdata = convertLOG( PEAK, Bdata, info.coding_domain );
          Gdata = convertLOG( PEAK, Gdata, info.coding_domain );
          Rdata = convertLOG( PEAK, Rdata, info.coding_domain );
        }
        /* convert RGB to YUV */
        pos = i * info.ywidth + j;

#ifdef RGBYUV
        RGB2YUV( Rdata, Gdata, Bdata, &( frame->Y[pos] ), &( frame->U[pos] ),
                 &( frame->V[pos] ) );

        if( info.coding_domain != LOG ) {
          frame->U[pos] += ( PEAK + 1 ) / 2;
          frame->V[pos] += ( PEAK + 1 ) / 2;
        } else {
          frame->U[pos] += 1024 / 2;
          frame->V[pos] += 1024 / 2;
        }
#else
        frame->Y[pos] = Rdata;
        frame->U[pos] = Gdata;
        frame->V[pos] = Bdata;
#endif
      }
    }


    //      sprintf(dpxname, "%s%04d%s", "temp", index, ".dpx");
    //      write_dpx(dpxname, *frame, info);

    free( RGBframe );
    free( dpxheader );

    break;
#endif
  default:
    printf( "image format error    format = %d(ioN.c)\n", format );
    exit( 0 );

  }


}

/*****************************************************************************/
/*                         read_jp2k_frame()                                 */
/*****************************************************************************/
void
read_jp2k_frame( YUVimage_ptr frame, videoinfo info, char *inname, int index,
            enum FORMAT format )
{
  int i, size, ret;
#ifdef SUPPORT_DPX
  int j, row, col, pos;
  float Rdata, Gdata, Bdata;
  U16 *RGBframe;
  struct dpx *dpxheader;
#endif
  int frameskip = 0; 
  long offset; 
  float *temp;
  char frame_name[250];
  //float *LRGB;
  FILE *fpio;

  offset = ( sizeof( float ) / sizeof( unsigned char ) ) * ( frameskip + 1 ) * index * ( info.ywidth*info.yheight + 2 * info.cwidth*info.cheight ); 

  switch ( format ) {
  case RAS:
    sprintf( frame_name, inname, index );
    ras2yuv( frame_name, info.ywidth, info.yheight, frame->Y, info.cwidth,
             info.cheight, frame->U, frame->V );
    break;
  case YUV:
    /* assume that the order of data file is Y[total], U[total], V[total] */
  
    if ((strchr (inname, '%'))!=0) //check the filename for frame or video mode 
      {
        sprintf(frame_name, inname, index);   // for frame-mode only
        if (!(fpio = fopen(frame_name, "rb")))
          {
            printf("read_frame: %s\n", frame_name);
            exit(1);
          }
      }
    else
      { // for container-mode (video)
        sprintf(frame_name, inname); 
        if (!(fpio = fopen(inname, "rb")))
          {
            printf("read_Video: %s\n", inname);
            exit(1);
          }
        fseek (fpio, offset, SEEK_SET); //seek to the picture in the file
        /* if ( fseek (fpio, offset, SEEK_SET)); //seek to the picture in the file with error check
           {
           printf("Error while seeking to picture");
           exit (1);  
           } */
      }

    size = ( info.ywidth ) * ( info.yheight );
    temp =
      ( float * )getarray( size, sizeof( float ), "temp" );
    if( ( ret = fread( temp, sizeof( float ), size, fpio ) ) != size ) {

      printf
        ( "error for reading file %s.Y in read_frame(): return size %d, target size %d\n",
          frame_name, ret, size );
      exit( 1 );
    }
    for( i = 0; i < size; i++ )
      frame->Y[i] = ( float )temp[i];

    size = ( info.cwidth ) * ( info.cheight );
    if( fread( temp, sizeof( float ), size, fpio ) !=
        ( size ) ) {
      printf( "error for reading file %s.U in read_frame()\n", frame_name );
      exit( 1 );
    }
    for( i = 0; i < size; i++ )
      frame->U[i] = ( float )temp[i];

    if( fread( temp, sizeof( float ), size, fpio ) !=
        ( size ) ) {
      printf( "error for reading file %s.V in read_frame()\n", frame_name );
      exit( 1 );
    }
    for( i = 0; i < size; i++ )
      frame->V[i] = ( float )temp[i];

    free( temp );
    fclose( fpio );
    break;
#ifdef SUPPORT_DPX
  case DPX:
    //printf("index %d %d\n", index, 0);
    sprintf( frame_name, inname, index );
    RGBframe = read_dpx( frame_name, &dpxheader, &row, &col );

    if( ( row != info.yheight ) || ( col != info.ywidth ) ) {
      printf( " the size of image does not match with the info\n" );
      exit( 1 );
    }


    for( i = 0; i < info.yheight; i++ ) {
      for( j = 0; j < info.ywidth; j++ ) {
        Bdata = RGBframe[i * ( 3 * info.ywidth ) + 3 * j];
        Gdata = RGBframe[i * ( 3 * info.ywidth ) + 3 * j + 1];
        Rdata = RGBframe[i * ( 3 * info.ywidth ) + 3 * j + 2];

        if( info.coding_domain != LOG ) {
          Bdata = convertLOG( PEAK, Bdata, info.coding_domain );
          Gdata = convertLOG( PEAK, Gdata, info.coding_domain );
          Rdata = convertLOG( PEAK, Rdata, info.coding_domain );
        }
        /* convert RGB to YUV */
        pos = i * info.ywidth + j;

#ifdef RGBYUV
        RGB2YUV( Rdata, Gdata, Bdata, &( frame->Y[pos] ), &( frame->U[pos] ),
                 &( frame->V[pos] ) );

        if( info.coding_domain != LOG ) {
          frame->U[pos] += ( PEAK + 1 ) / 2;
          frame->V[pos] += ( PEAK + 1 ) / 2;
        } else {
          frame->U[pos] += 1024 / 2;
          frame->V[pos] += 1024 / 2;
        }
#else
        frame->Y[pos] = Rdata;
        frame->U[pos] = Gdata;
        frame->V[pos] = Bdata;
#endif
      }
    }


    //      sprintf(dpxname, "%s%04d%s", "temp", index, ".dpx");
    //      write_dpx(dpxname, *frame, info);

    free( RGBframe );
    free( dpxheader );

    break;
#endif
  default:
    printf( "image format error    format = %d(ioN.c)\n", format );
    exit( 0 );

  }


}

/*****************************************************************************/
/*                               write_frame()                               */
/*****************************************************************************/
void
write_frame( YUVimage frame, videoinfo info, char *inname, int index,
             enum FORMAT format )
{
  FILE *fpio;
  int i, size;
  char frame_name[250];
  unsigned char *temp;
  long int offset;
#ifdef SUPPORT_DPX
  int j, pos;
  float Ydata, Udata, Vdata;
#endif
  float *LRGB = 0;
  int frameskip = 0; 
  
//  printf("Enter write frame!\n");

  offset = ( frameskip + 1 ) * index * ( info.ywidth*info.yheight + 2 * info.cwidth*info.cheight ); 

  switch ( format ) {
  case RAS:
    sprintf( frame_name, inname, index );
    yuv2ras( frame_name, info.ywidth, info.yheight, frame.Y, info.cwidth,
             info.cheight, frame.U, frame.V );
    break;

  case YUV:
    /* assume that the order of data file is Y[total], U[total], V[total] */
    if ((strchr (inname, '%'))!=0)
      {
		assert(0);
        sprintf(frame_name, inname, index);
      
        if (!(fpio = fopen(frame_name, "wb"))) // open a new frame file
          {
            printf("write_frame: %s\n", frame_name);
            exit(1);
          }
      }
    else
      {
        if ( index == info.start ) // open new video file
          {  
            if (!(fpio = fopen(inname, "wb")))
              {
                printf("write_Video: %s\n", inname);
                exit(1);
              }
          }
        else // append new frame to existing video file
          {
            if (!(fpio = fopen(inname, "r+b")))
              {
                printf("write_Video: %s\n", inname);
                exit(1);
              }
          }
		  fseek (fpio, offset, SEEK_SET);
      }

    size = ( info.ywidth ) * ( info.yheight );
    temp = ( unsigned char * )calloc( size, sizeof( unsigned char ) );
    for( i = 0; i < size; i++ )
      temp[i] = ( unsigned char )clipping( frame.Y[i] );
    if( fwrite( temp, sizeof( unsigned char ), size, fpio ) !=
        ( unsigned )( size ) ) {
      printf( "error in write_frame()\n" );
      exit( 1 );
    }

    size = ( info.cwidth ) * ( info.cheight );  /* QUEST */
    for( i = 0; i < size; i++ )
      temp[i] = ( unsigned char )clipping( frame.U[i] );
    if( fwrite( temp, sizeof( unsigned char ), size, fpio ) !=
        ( unsigned )( size ) ) {
      printf( "error in write_frame()\n" );
      exit( 1 );
    }

    for( i = 0; i < size; i++ )
      temp[i] = ( unsigned char )clipping( frame.V[i] );
    if( fwrite( temp, sizeof( unsigned char ), size, fpio ) !=
        ( unsigned )( size ) ) {
      printf( "error in write_frame()\n" );
      exit( 1 );
    }
    free( temp );
    fclose( fpio );

    break;

  case DPX:

#ifdef DPX_SUPPORT
    LRGB =
      ( float * )getarray( info.ywidth * info.yheight * 3, sizeof( float ),
                           "LRGB" );

    for( i = 0; i < info.yheight; i++ ) {
      for( j = 0; j < info.ywidth; j++ ) {
        pos = i * info.ywidth + j;
        Ydata = frame.Y[pos];
        Udata = frame.U[pos];
        Vdata = frame.V[pos];
        /* YUV -> RGB */
        if( info.coding_domain != LOG ) {
          Udata -= ( PEAK + 1 ) / 2.;
          Vdata -= ( PEAK + 1 ) / 2.;
        } else {
          Udata -= 1024 / 2.;
          Vdata -= 1024 / 2.;
        }
        YUV2RGB( Ydata, Udata, Vdata, &( LRGB[pos * 3 + 2] ),
                 &( LRGB[pos * 3 + 1] ), &( LRGB[pos * 3] ) );
      }
    }
    //printf("Check the DPX header (ioN.c) \n"); getchar();

    sprintf( frame_name, inname, index );

    save_dpx_header( info, index );     // save the header from the original dpx file for the reconstructed frame

    write_dpx( frame_name, LRGB, info );
#else
    assert(0);
#endif

#ifdef OUTPUT_RAS
    unsigned char *RGBframe;
    int data;
    struct rasterfile header;
    char rasname[250];

    RGBframe =
      ( unsigned char * )getarray( info.ywidth * info.yheight * 3,
                                   sizeof( unsigned char ), "LRGB" );

    strcpy( rasname, frame_name );
    strcat( rasname, ".ras" );

    header = make_header( info.ywidth, info.yheight, RGBDATA );

    switch ( info.coding_domain ) {
    case VIDEO:
      for( i = 0; i < info.ywidth * info.yheight * 3; i++ ) {

        data = nint( LRGB[i] * 255. / PEAK );
        if( data < 0 )
          RGBframe[i] = 0;
        else if( data > 255 )
          RGBframe[i] = 255;
        else
          RGBframe[i] = data;
      }
      break;

    case LINEAR:               // so we need to do Gamma correction
      for( i = 0; i < info.ywidth * info.yheight * 3; i++ ) {

        LRGB[i] = LRGB[i] / PEAK;
        if( LRGB[i] > 1 )
          LRGB[i] = 1;
        else if( LRGB[i] < 0 )
          LRGB[i] = 0;

        data = nint( 255 * pow( LRGB[i], 1. / Gamma ) );
        if( data < 0 )
          RGBframe[i] = 0;
        else if( data > 255 )
          RGBframe[i] = 255;
        else
          RGBframe[i] = data;
      }
      break;

    case LOG:
      for( i = 0; i < info.ywidth * info.yheight * 3; i++ ) {

        //LRGB[i] = LRGB[i]*1023/PEAK;

        if( LRGB[i] < 0 )
          LRGB[i] = 0;
        else if( LRGB[i] > 1023 )
          LRGB[i] = 1023;

        RGBframe[i] = nint( convertLOG( 255, LRGB[i], VIDEO ) );
      }
      break;

    default:
      printf( "coding_domain error (ioN.c)\n" );
      exit( 1 );
    }

    write_ras( rasname, header, RGBframe );

    free( RGBframe );
#endif

    free( LRGB );
    break;

  default:
    printf( "image format error format = %d(ioN.c)\n", format );
    exit( 0 );

  }

}

/*****************************************************************************/
/*                          write_jp2k_frame()                               */
/*****************************************************************************/
void
write_jp2k_frame( YUVimage frame, videoinfo info, char *inname, int index,
             enum FORMAT format )
{
  FILE *fpio;
  int i, size;
  char frame_name[250];
  float *temp;
#ifdef SUPPORT_DPX
  int j, pos;
  float Ydata, Udata, Vdata;
#endif
  float *LRGB = 0;

  switch ( format ) {
  case RAS:
    sprintf( frame_name, inname, index );
    yuv2ras( frame_name, info.ywidth, info.yheight, frame.Y, info.cwidth,
             info.cheight, frame.U, frame.V );
    break;

  case YUV:
    /* assume that the order of data file is Y[total], U[total], V[total] */
    if ((strchr (inname, '%'))!=0)
      {
        sprintf(frame_name, inname, index);
      
        if (!(fpio = fopen(frame_name, "wb"))) // open a new frame file
          {
            printf("write_frame: %s\n", frame_name);
            exit(1);
          }
      }
    else
      {
        if ( index == info.start ) // open new video file
          {  
            if (!(fpio = fopen(inname, "wb")))
              {
                printf("write_Video: %s\n", inname);
                exit(1);
              }
          }
        else // append new frame to existing video file
          {
            if (!(fpio = fopen(inname, "ab")))
              {
                printf("write_Video: %s\n", inname);
                exit(1);
              }
          }
      }

    size = ( info.ywidth ) * ( info.yheight );
    temp = ( float * )calloc( size, sizeof( float ) );
    for( i = 0; i < size; i++ )
      temp[i] = ( frame.Y[i] );
    if( fwrite( temp, sizeof( float ), size, fpio ) !=
        ( size ) ) {
      printf( "error in write_frame()\n" );
      exit( 1 );
    }

    size = ( info.cwidth ) * ( info.cheight );  /* QUEST */
    for( i = 0; i < size; i++ )
      temp[i] = ( frame.U[i] );
    if( fwrite( temp, sizeof( float ), size, fpio ) !=
       ( size ) ) {
      printf( "error in write_frame()\n" );
      exit( 1 );
    }

    for( i = 0; i < size; i++ )
      temp[i] = ( frame.V[i] );
    if( fwrite( temp, sizeof( float ), size, fpio ) !=
        ( size ) ) {
      printf( "error in write_frame()\n" );
      exit( 1 );
    }
    free( temp );
    fclose( fpio );

    break;

  case DPX:

#ifdef DPX_SUPPORT
    LRGB =
      ( float * )getarray( info.ywidth * info.yheight * 3, sizeof( float ),
                           "LRGB" );

    for( i = 0; i < info.yheight; i++ ) {
      for( j = 0; j < info.ywidth; j++ ) {
        pos = i * info.ywidth + j;
        Ydata = frame.Y[pos];
        Udata = frame.U[pos];
        Vdata = frame.V[pos];
        /* YUV -> RGB */
        if( info.coding_domain != LOG ) {
          Udata -= ( PEAK + 1 ) / 2.;
          Vdata -= ( PEAK + 1 ) / 2.;
        } else {
          Udata -= 1024 / 2.;
          Vdata -= 1024 / 2.;
        }
        YUV2RGB( Ydata, Udata, Vdata, &( LRGB[pos * 3 + 2] ),
                 &( LRGB[pos * 3 + 1] ), &( LRGB[pos * 3] ) );
      }
    }
    //printf("Check the DPX header (ioN.c) \n"); getchar();

    sprintf( frame_name, inname, index );

    save_dpx_header( info, index );     // save the header from the original dpx file for the reconstructed frame

    write_dpx( frame_name, LRGB, info );
#else
    assert(0);
#endif

#ifdef OUTPUT_RAS
    unsigned char *RGBframe;
    int data;
    struct rasterfile header;
    char rasname[250];

    RGBframe =
      ( unsigned char * )getarray( info.ywidth * info.yheight * 3,
                                   sizeof( unsigned char ), "LRGB" );

    strcpy( rasname, frame_name );
    strcat( rasname, ".ras" );

    header = make_header( info.ywidth, info.yheight, RGBDATA );

    switch ( info.coding_domain ) {
    case VIDEO:
      for( i = 0; i < info.ywidth * info.yheight * 3; i++ ) {

        data = nint( LRGB[i] * 255. / PEAK );
        if( data < 0 )
          RGBframe[i] = 0;
        else if( data > 255 )
          RGBframe[i] = 255;
        else
          RGBframe[i] = data;
      }
      break;

    case LINEAR:               // so we need to do Gamma correction
      for( i = 0; i < info.ywidth * info.yheight * 3; i++ ) {

        LRGB[i] = LRGB[i] / PEAK;
        if( LRGB[i] > 1 )
          LRGB[i] = 1;
        else if( LRGB[i] < 0 )
          LRGB[i] = 0;

        data = nint( 255 * pow( LRGB[i], 1. / Gamma ) );
        if( data < 0 )
          RGBframe[i] = 0;
        else if( data > 255 )
          RGBframe[i] = 255;
        else
          RGBframe[i] = data;
      }
      break;

    case LOG:
      for( i = 0; i < info.ywidth * info.yheight * 3; i++ ) {

        //LRGB[i] = LRGB[i]*1023/PEAK;

        if( LRGB[i] < 0 )
          LRGB[i] = 0;
        else if( LRGB[i] > 1023 )
          LRGB[i] = 1023;

        RGBframe[i] = nint( convertLOG( 255, LRGB[i], VIDEO ) );
      }
      break;

    default:
      printf( "coding_domain error (ioN.c)\n" );
      exit( 1 );
    }

    write_ras( rasname, header, RGBframe );

    free( RGBframe );
#endif

    free( LRGB );
    break;

  default:
    printf( "image format error format = %d(ioN.c)\n", format );
    exit( 0 );

  }

}

void
scaling( int curr, YUVimage ** pyrTemp, YUVimage * pyrFrs, videoinfo info,
         enum FLAG first_GOP, int remaining_frs )
{
  int i, nfrs, t_level;
  float weight;

  temporal_filter ();
  t_level = MY_MIN (info.t_level, info.tPyrLev);
   
  if ( remaining_frs != 0 ){ 
    nfrs = info.GOPsz;
  } else {
    nfrs = info.eff_GOPsz;
  }
  nfrs = nfrs >> t_level;
  if ( first_GOP == YES ) nfrs = 1;

  weight = ( float )pow( 1 / LPW4[1], t_level );
  
  if( t_level < info.tPyrLev ) {
    for( i = 0; i < nfrs; i++ ) {
      wcopyframe( &pyrTemp[t_level][i], &pyrTemp[t_level][i], weight, info );
    }
  } else { // t_level = info.tPyrLev
    wcopyframe( &pyrFrs[0], &pyrFrs[0], weight, info );
  }
}

void
write_GOP( int curr, YUVimage ** pyrTemp, YUVimage * pyrFrs, videoinfo info,
           enum FLAG first_GOP, int remaining_frs )
{
  int i, nfrs, dist, t_level; 

  // assert (info.t_level <= info.tPyrLev);

  t_level = MY_MIN (info.t_level, info.tPyrLev);
  dist    = 1 << ( t_level );

  if( info.tPyrLev > 0 && info.t_level > 0 ) 
    scaling( curr, pyrTemp, pyrFrs, info, first_GOP, remaining_frs );

  if( t_level < info.tPyrLev ) { // default case
  
    if( first_GOP == YES ){
      write_frame( pyrTemp[t_level][0], info, info.decname, curr, info.format );
    } else {
      if ( remaining_frs != 0 ){
        nfrs = info.GOPsz;
      } else {
        nfrs = info.eff_GOPsz;
      }

      for( i = 0; i < nfrs >> t_level; i++ ) {
        if ( (curr - info.GOPsz + (i+1) * dist) <= info.last ){
          write_frame( pyrTemp[t_level][i], info, info.decname, curr - info.GOPsz + (i+1) * dist, info.format );
          // printf(" write_frame( pyrTemp[%d][%d], curr - info.GOPsz + (i+1) * dist %d );\n", 
          //        t_level, i, curr - info.GOPsz + (i+1) * dist );
        }
      }
    }
    
  } else {
    if ( curr <= info.last ){
      write_frame( pyrFrs[0], info, info.decname, curr, info.format );
    }
  }
  
}


void  rec_write_block_mode_motion_vector( vector_ptr fmv, int cx, int cy, 
										   int xblk, int yblk, int hor, int ver, 
										   int t_level, videoinfo info, FILE *fpmvio, int  frame_type )
{
  int xblock, yblock;

  if( fmv->child ) {
    rec_write_block_mode_motion_vector( fmv->child0, cx, cy, 
										xblk/2, yblk/2, hor,  ver, t_level,  info, fpmvio, frame_type );

    rec_write_block_mode_motion_vector( fmv->child1, cx+xblk/2, cy, 
										xblk/2, yblk/2, hor,  ver, t_level,  info, fpmvio, frame_type );

    rec_write_block_mode_motion_vector( fmv->child2, cx, cy+yblk/2, 
										xblk/2, yblk/2, hor,  ver, t_level,  info, fpmvio, frame_type );

    rec_write_block_mode_motion_vector( fmv->child3, cx+xblk/2, cy+yblk/2, 
										xblk/2, yblk/2, hor,  ver, t_level,  info, fpmvio, frame_type );

  } else {
    /* consider the small block around the boundaries */
    xblock = ( cx + xblk <= hor ) ? xblk : hor - cx;
    yblock = ( cy + yblk <= ver ) ? yblk : ver - cy;

    if( xblock <= 0 || yblock <= 0 )		return;
 
    // the format of block_mode_motion_info file 
	fprintf(fpmvio, "%d\t %d\t %d\t %d\t %d\t", cx, cy, xblock, yblock, fmv->bi_mode);
	if (frame_type==1 || frame_type==2) // uni-left or uni-right
	{
		if ( fmv->bi_mode != DIRECTIONAL_IBLOCK )
			//            indicator of mv, mvx, mvy
			fprintf(fpmvio, "%d\t %f\t %f\n", 1, fmv->mvx, fmv->mvy);
		else
			fprintf(fpmvio, "%d\t %f\t %f\n", 0, fmv->iblock_spatial_mode*1.0, 0);
	}

	if (frame_type==0) // bi-directional
	{
		if ( fmv->lifting_mode==CONNECTED || fmv->lifting_mode==PREDICTED )
			fprintf(fpmvio, "%d\t %f\t %f\n", 1, fmv->mvx, fmv->mvy);
		else if  ( fmv->lifting_mode==SPATIAL_PREDICTED)
			fprintf(fpmvio, "%d\t %f\t %f\n", 0, fmv->iblock_spatial_mode*1.0, 0);
		else
			fprintf(fpmvio, "%d\t %f\t %f\n", 0, 0, 0);
	}

  }
}


void write_block_mode_motion_vector(char *direction, int GOP_counter, int count, int  frame_type, 
									vector_ptr fmv,  videoinfo info, int t_level)
{

  int x, y, X, Y, yhor, yver;
  int xnum, ynum, xblk, yblk;
  char file_name[80];
  FILE *fpmvio; 

  sprintf(file_name, "%s_block_mode_mv_GOP%03d_count%03d_frametype%1d.txt", 
					direction, GOP_counter, count, frame_type); 
  fpmvio = fopen(file_name, "wt"); 
  yhor = info.ywidth;
  yver = info.yheight;
  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];
  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
    for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
      rec_write_block_mode_motion_vector( &fmv[Y * xnum + X], x, y, xblk, yblk, yhor, yver, 
										  t_level, info, fpmvio, frame_type );
    }
  }
  fclose(fpmvio); 

}



// write high temporal frames into files
void write_frame_into_file(YUVimage frame, videoinfo info,  int index, int GOP_counter)
{
	char framenameY[80], framenameU[80], framenameV[80];
	int x, y, yhor, yver;
	FILE *fpimgY, *fpimgU, *fpimgV;
    
	yhor = info.ywidth;  yver = info.yheight;

	// write high temporal frame into text file to have a look
	sprintf(framenameY, "highframeY_GOP%03d_index%d.txt", GOP_counter, index);
	sprintf(framenameU, "highframeU_GOP%03d_index%d.txt", GOP_counter, index);
	sprintf(framenameV, "highframeV_GOP%03d_index%d.txt", GOP_counter, index);
	fpimgY = fopen( framenameY, "wt");
	fpimgU = fopen( framenameU, "wt");
	fpimgV = fopen( framenameV, "wt");
	for (y=0; y<yver; y++)
	{
		for (x=0; x<yhor;  x++)	
		{
			fprintf(fpimgY, "%f\t", frame.Y[y*yhor+x]);
			if (x%2==0 && y%2==0)
			{
				fprintf(fpimgU, "%f\t", frame.U[y/2*yhor/2+x/2]);
				fprintf(fpimgV, "%f\t", frame.V[y/2*yhor/2+x/2]);

			}
		}
		fprintf(fpimgY, "\n");
		fprintf(fpimgU, "\n");
		fprintf(fpimgV, "\n");
	}
    fclose(fpimgY);
	fclose(fpimgU);
	fclose(fpimgV);

}


