#include "basic.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "rasterfile.h"
#include "structN.h"
#include "dpx.h"
#include "unix_pc.h"
#include "chrom.h"
#include "memoryN.h"

#define Detector 0xffc00000
float clipping2( float ipix );

void dpx2ras( char *dpxname, char *outname );
U16 colorLUT[1024];


U16
Release( U32 * data, U8 bitdepth )
{
  U32 temp;
  U16 outdata;

  outdata = 0;

  temp = ( *data ) & Detector;
  outdata = temp >> ( 32 - bitdepth );
  ( *data ) = ( *data ) << bitdepth;

  return outdata;
}


void
Pack( U32 * data, int comp, U8 bitdepth )
{
  ( *data ) <<= bitdepth;
  comp = comp & 0x3ff;
  ( *data ) = ( *data ) | comp;
}



U16 *
read_dpx( char *dpxname, struct dpx **dpxheader, int *nr, int *nc )
{
  int nn;
  int i, j, pos;
  FILE *fp = NULL;
  U32 *data = NULL, rgbU32, offset;
  U16 *RGBframe;
  U8 bitdepth;

  if( ( fp = fopen( dpxname, "rb" ) ) == NULL ) {
    fprintf( stderr, "Can't open %s for read\n", dpxname );
    exit( 1 );
  }

  *dpxheader =
    ( struct dpx * )getarray( 1, sizeof( struct dpx ), "*dpxheader" );

  if( fread( *dpxheader, sizeof( struct dpx ), 1, fp ) != 1 ) {
    printf( "read dpx file header error \n" );
    exit( 1 );
  }

  *nc = exchange_4byte_order( ( *dpxheader )->imageinfo.pixels_per_line );
  *nr = exchange_4byte_order( ( *dpxheader )->imageinfo.lines_per_image_ele );

  offset = exchange_4byte_order( ( *dpxheader )->fileinfo.offset );
  printf( "offset = %d\n", offset );
  //total = exchange_4byte_order((*dpxheader)->fileinfo.file_size); printf("total file size= %d\n", total);
/* 
  if( (nr!=1080) || (nc!=1920) ){
	  printf("the height of image %d and the width of image are error\n", nr, nc);
	  exit(0);
  }
*/
  //*nr = 540; *nc = 960;
  nn = ( *nr ) * ( *nc );
  //printf("row %d col %d pixels %d\n", *nr, *nc, nn);


  fseek( fp, offset, SEEK_SET );


  data = ( U32 * ) getarray( nn, sizeof( U32 ), "data" );
  if( fread( data, sizeof( U32 ), nn, fp ) != ( unsigned )nn ) {
    printf( "read image data error\n" );
    exit( 1 );
  }
  fclose( fp );

  for( i = 0; i < nn; i++ ) {
    data[i] = exchange_4byte_order( data[i] );
  }

  /*  seperate 4 byte data into R, G, B */
  bitdepth = ( *dpxheader )->imageinfo.bit_size;
  RGBframe = ( U16 * ) getarray( nn * 3, sizeof( U16 ), "RGBframe" );

  for( i = 0; i < ( *nr ); i++ ) {
    for( j = 0; j < ( *nc ); j++ ) {
      pos = i * ( *nc ) + j;

      rgbU32 = data[pos];

      RGBframe[i * 3 * ( *nc ) + 3 * j + 2] = Release( &rgbU32, bitdepth );

      RGBframe[i * 3 * ( *nc ) + 3 * j + 1] = Release( &rgbU32, bitdepth );

      RGBframe[i * 3 * ( *nc ) + 3 * j] = Release( &rgbU32, bitdepth );
    }
  }
  free( data );
  return RGBframe;
}




struct dpx *
make_dpxheader( videoinfo info )
{
  struct dpx *dpxheader;

  dpxheader = read_dpx_header( "dpxheader" );

  dpxheader->imageinfo.pixels_per_line = exchange_4byte_order( info.ywidth );
  dpxheader->imageinfo.lines_per_image_ele =
    exchange_4byte_order( info.yheight );
  //dpxheader->imageinfo.bit_size = 10;
  return ( dpxheader );
}

#ifdef DPX_SUPPORT
void
write_dpx( char *dpxname, float *LRGB, videoinfo info )
{
  int i, j, pos;
  int temp;
  U8 bitdepth;
  int hor, ver, nn;
  FILE *fp = NULL;
  U32 *data = NULL;
  struct dpx *dpxheader;
  //*RGBframe,


  ver = info.yheight;
  hor = info.ywidth;
  nn = ver * hor;

  dpxheader = make_dpxheader( info );

  data = ( U32 * ) getarray( nn, sizeof( U32 ), "data" );

  if( !( fp = fopen( dpxname, "wb" ) ) ) {
    printf( "write_frame: %s\n", dpxname );
    exit( 1 );
  }
//        RGBframe = Nonlinearize(LRGB, ver, hor, info);

  for( i = 0; i < ver; i++ ) {
    for( j = 0; j < hor; j++ ) {
      pos = i * hor + j;

      bitdepth = dpxheader->imageinfo.bit_size; // no need to exchange order. it is U8 data 
      data[pos] = 0;

      temp = nint( LRGB[pos * 3 + 2] ); // R 
      if( temp > info.max )
        temp = info.max;
      if( temp < info.min )
        temp = info.min;
      Pack( &data[pos], temp, bitdepth );

      temp = nint( LRGB[pos * 3 + 1] ); // G
      if( temp > info.max )
        temp = info.max;
      if( temp < info.min )
        temp = info.min;
      Pack( &data[pos], temp, bitdepth );

      temp = nint( LRGB[pos * 3] );     // B
      if( temp > info.max )
        temp = info.max;
      if( temp < info.min )
        temp = info.min;
      Pack( &data[pos], temp, bitdepth );

      data[pos] <<= 2;

      data[pos] = exchange_4byte_order( data[pos] );
    }
  }



  if( fwrite( dpxheader, sizeof( struct dpx ), 1, fp ) != 1 ) {
    printf( "error in write_frame()\n" );
    exit( 1 );
  }

  if( fwrite( data, sizeof( U32 ), nn, fp ) != ( unsigned )nn ) {
    printf( "error in write_frame()\n" );
    exit( 1 );
  }

  fclose( fp );
  free( data );                 //free(RGBframe);
  free( dpxheader );
}
#endif


void
dpx2ras( char *dpxname, char *outname, enum CODING_DOMAIN coding_domain )
{
  int nr, nc, pos, downnr, downnc, downpos, t;
  int i, j, k, r;
  U8 *temp;
  U16 *RGBframe;
  U16 Bmin, Bmax, Gmin, Gmax, Rmin, Rmax;
  struct dpx *dpxheader;
  struct rasterfile rhead;

  RGBframe = read_dpx( dpxname, &dpxheader, &nr, &nc );


  //printf("%d %d\n", nr, nc);


  dpxrange( RGBframe, nr, nc, &Bmin, &Bmax, &Gmin, &Gmax, &Rmin, &Rmax );
  printf( "image min(R %d G %d B %d), image max(R %d G %d B %d) \n", Rmin,
          Gmin, Bmin, Rmax, Gmax, Bmax );


  r = 1;
  downnr = nr / r;
  downnc = nc / r;
  temp = ( U8 * ) getarray( downnr * downnc * 3, sizeof( U8 ), "temp" );
  for( i = 0; i < downnr; i++ ) {
    for( j = 0; j < downnc; j++ ) {
      pos = ( i * r ) * ( nc * 3 ) + ( j * r ) * 3;     /* 3 comes from R, G, B */
      downpos = i * ( downnc * 3 ) + j * 3;
      for( k = 0; k < 3; k++ ) {

        if( coding_domain != LOG )
          t = colorLUT[RGBframe[pos + k]];
        else
          t = ( RGBframe[pos + k] + 2 ) >> 2;

        if( t > 255 )
          temp[downpos + k] = 255;
        else if( t < 0 )
          temp[downpos + k] = 0;
        else
          temp[downpos + k] = ( U8 ) t;
      }
    }
  }

  rhead = make_header( downnc, downnr, RGBDATA );
  write_ras( outname, rhead, temp );

  free( temp );
  free( dpxheader );
  free( RGBframe );
  return;
}

#ifdef DPX_SUPPORT
void
dpx2dpx( char *dpxname, char *outname )
{
  int nr, nc, pos, pos2;
  int i, j;
  U16 *RGBframe;
  float Ydata, Udata, Vdata;
  struct dpx *dpxheader;
  YUVimage YUV;
  float *LRGB;                  //*RGBframe,

  videoinfo info;

  RGBframe = read_dpx( dpxname, &dpxheader, &nr, &nc );


  info.ywidth = nc;
  info.yheight = nr;
  info.cwidth = nc;
  info.cheight = nr;

  frame_alloc( &YUV, info );

  for( i = 0; i < nr; i++ ) {
    for( j = 0; j < nc; j++ ) {
      pos = i * ( nc * 3 ) + j * 3;     /* 3 comes from R, G, B */
      pos2 = i * nc + j;
      RGB2YUV( ( float )RGBframe[pos + 2], ( float )RGBframe[pos + 1],
               ( float )RGBframe[pos], &( YUV.Y[pos2] ), &( YUV.V[pos2] ),
               &( YUV.U[pos2] ) );
      YUV.U[pos2] += 512.;
      YUV.V[pos2] += 512.;
    }
  }


  info.coding_domain = LOG;
  info.max = 1023;
  info.min = 0;



  LRGB =
    ( float * )getarray( info.ywidth * info.yheight * 3, sizeof( float ),
                         "LRGB" );

  for( i = 0; i < info.yheight; i++ ) {
    for( j = 0; j < info.ywidth; j++ ) {
      pos = i * info.ywidth + j;
      Ydata = YUV.Y[pos];
      Udata = YUV.U[pos];
      Vdata = YUV.V[pos];
      /* YUV -> RGB */
      Udata -= 1024 / 2.;
      Vdata -= 1024 / 2.;
      YUV2RGB( Ydata, Udata, Vdata, &( LRGB[pos * 3 + 2] ),
               &( LRGB[pos * 3 + 1] ), &( LRGB[pos * 3] ) );
    }
  }



  write_dpx( outname, LRGB, info );



  free_frame( YUV );
  free( LRGB );
  free( dpxheader );
  free( RGBframe );
  return;
}
#endif

void
dpxrange( U16 * RGBframe, int row, int col, U16 * Bmin, U16 * Bmax,
          U16 * Gmin, U16 * Gmax, U16 * Rmin, U16 * Rmax )
{
  int i;

  *Bmax = 0;
  *Bmin = 1024;

  *Gmax = 0;
  *Gmin = 1024;

  *Rmax = 0;
  *Rmin = 1024;

  for( i = 0; i < row * col; i++ ) {

    if( RGBframe[i * 3] < ( *Bmin ) ) {
      *Bmin = RGBframe[i * 3];
    }
    if( RGBframe[i * 3] > ( *Bmax ) ) {
      *Bmax = RGBframe[i * 3];
    }

    if( RGBframe[i * 3 + 1] < ( *Gmin ) ) {
      *Gmin = RGBframe[i * 3 + 1];
    }
    if( RGBframe[i * 3 + 1] > ( *Gmax ) ) {
      *Gmax = RGBframe[i * 3 + 1];
    }
    if( RGBframe[i * 3 + 2] < ( *Rmin ) ) {
      *Rmin = RGBframe[i * 3 + 2];
    }
    if( RGBframe[i * 3 + 2] > ( *Rmax ) ) {
      *Rmax = RGBframe[i * 3 + 2];
    }
  }
//  printf("min = %d, max = %d (dpxio.c)\n", *min, *max);
}

/*
 * convert 10 bits LOG data to "coding_domain"
 * #define Gamma 2.2
 * #define Refblack  95
 * #define Refwhite  685
 * #define KODAK
 */
float
convertLOG( int peak, float data, enum CODING_DOMAIN coding_domain )
{
  float newdata;
  //int i;
  double linear;                //, shift, gain

#ifdef KODAK
  k = 0.002 / 0.6;

  if( data <= Refblack )
    newdata = 0;
  else if( data >= Refwhite )
    newdata = peak;
  else {

    range = pow( 10, Refwhite * k ) - pow( 10, Refblack * k );

    /* direct from The Cineon Digital Film System
       shift = data-Refwhite;
       shift *= k;
       shift  = pow(10, shift);
       gain = 1/(1-pow(10, (Refblack-Refwhite)*k));

       linear = shift*gain-(gain-1);
     */
    linear = ( pow( 10, data * k ) - pow( 10, Refblack * k ) ) / range; // equivalent to the upper one, but easy to understand

    if( coding_domain == LINEAR )
      newdata = peak * linear;
    else if( coding_domain == VIDEO )
      newdata = peak * pow( linear, 1. / Gamma );
    else {
      printf( "coding_domain error (dpxio.c)\n" );
      exit( 1 );
    }

  }

  return newdata;

#else //cintel
  double black, white;

  black = pow( 2., Refblack / 102.4 ) / 1024;
  white = pow( 2., Refwhite / 102.4 ) / 1024;

  linear = ( pow( 2., data / 102.4 ) / 1024 - black ) / ( white - black );
  if( linear > 1 )
    linear = 1.;
  if( linear < 0 )
    linear = 0;

  if( coding_domain == LINEAR )
    newdata = ( float )( peak * linear );
  else if( coding_domain == VIDEO )
    newdata = ( float )( peak * pow( linear, 1. / Gamma ) );
  else {
    printf( "coding_domain error (dpxio.c)\n" );
    exit( 1 );
  }

  return newdata;

#endif


}

void
compute_colorLUT( int bits, enum CODING_DOMAIN coding_domain )
{
  int i;
  int peak;

  peak = ( 1 << bits ) - 1;     //printf("peak = %d \n", peak);

  if( coding_domain == LOG )
    return;

  for( i = 0; i < 1024; i++ ) {
    colorLUT[i] = nint( convertLOG( peak, ( float )i, coding_domain ) );
    //printf("i %d, LUT %d\n", i, colorLUT[i]); getchar();
  }
}

/*float *linearize(U16 *RGBframe, int row, int col, videoinfo info)
{
	int i;
    U16 Bmin, Bmax, Gmin, Gmax, Rmin, Rmax;
	float *LRGB;

	LRGB = (float *)getarray(row*col*3, sizeof(float), "LRGB");

	if( (info.log == YES) && (info.linear ==YES) ){
    	dpxrange(RGBframe, row, col, &Bmin, &Bmax, &Gmin, &Gmax, &Rmin, &Rmax);
        if( (Rmin < info.min) || (Rmax > info.max) ||(Gmin < info.min) || (Gmax > info.max) 
		     	||(Bmin < info.min) || (Bmax > info.max)){
		   printf("image min(R%d G%d B%d) != info.min(%d), image max(R%d G%d B%d) != info.max(%d)\n", Rmin, Gmin, Bmin, info.min, Rmax, Gmin, Bmin, info.max);
           exit(1);
		}
	}
	else{
        for(i=0; i<row*col*3; i++){
		    LRGB[i] = (float)RGBframe[i];
		}
	}
	return LRGB;
}
*/


/*float *Nonlinearize(float *LRGB, int row, int col, videoinfo info)
{
	int i;
	double shift, k, range, tmp;
    float *RGBframe;

    	
	k = 0.002/0.6/Gamma;

		range = info.max-info.min;
        range *= k;
        range = pow(10, range)-1;

	RGBframe = (float *)getarray(row*col*3, sizeof(float), "RGBframe");
			   
    if( (info.log == YES) && (info.linear ==YES) ){
		for(i=0; i<row*col*3; i++){

          tmp = LRGB[i]*range/1024;
		  if(tmp > 0){
 		     shift = log10(tmp+1);
		  }
		  else{
			shift = 0;
		  }
		  RGBframe[i] = (float)(shift/k + info.min);
		}
	}
	else{
		for(i=0; i<row*col*3; i++){
			RGBframe[i] = LRGB[i]; 
		}
	}


   return RGBframe;
}
*/

/* save the dpx header copied from the original dpx file */
void
save_dpx_header( videoinfo info, int index )
{
  char frame_name[250];
  struct dpx origheader;
  FILE *fp;
  sprintf( frame_name, info.inname, index );
  //sprintf(frame_name, info.inname);
  if( !( fp = fopen( frame_name, "rb" ) ) ) {
    printf( "can not open: %s (dpxio.c)\n", frame_name );
    exit( 1 );
  }

  if( fread( &origheader, sizeof( struct dpx ), 1, fp ) != 1 ) {
    printf( "read dpx file header error \n" );
    exit( 1 );
  }
  fclose( fp );

  if( !( fp = fopen( "dpxheader", "wb" ) ) ) {
    printf( "can not open file (in dpxio.c)\n" );
    exit( 1 );
  }

  if( fwrite( &origheader, sizeof( struct dpx ), 1, fp ) != 1 ) {
    printf( "write dpx file header error \n" );
    exit( 1 );
  }
  fclose( fp );
}


struct dpx *
read_dpx_header( char *filename )
{
  FILE *fp;
  struct dpx *origheader;

  origheader =
    ( struct dpx * )getarray( 1, sizeof( struct dpx ), "origheader" );

  fp = fopen( filename, "rb" );
  if( fp == NULL ) {
    fp = fopen( "z:\\3d_ezbc\\exp\\dpxheader", "rb" );
    if( fp == NULL ) {
      printf( "can not find the dpx header file ( in dpxio.c )\n" );
      exit( 1 );
    } else {
      printf( "use the default dpx header\n" );
    }
  }

  if( fread( origheader, sizeof( struct dpx ), 1, fp ) != 1 ) {
    printf( "read dpx file header error \n" );
    exit( 1 );
  }
  fclose( fp );

  return ( origheader );
}
