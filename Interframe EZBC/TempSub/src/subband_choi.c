#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "basic.h"
#include "structN.h"
#include "Choisubband.h"
#define   LL     1
#define   HL     2
#define   LH     3
#define   HH     4


int filter_coeff( float **lpf, float **hpf, int filterType, int analsyn );
void filter_2d( float *obuf, float *ibuf, float *buf1, float *buf2, int lenx,
                int leny, float *lpf, float *hpf, int tap, int mode );
void interpol_2d( float *obuf, float *in, float *ibuf, float *buf1,
                  float *buf2, float *buf3, float *lpf, float *hpf, int tap,
                  int lenx, int leny, int sx, int sy, int hor, int mode );
void filter_1d( float *out, float *in, int length, float *filter, int tap,
                int start );
void filter_1d_normal( float *out, float *in, int length, float *filter,
                       int tap, int start );


/*****************************************************************************/
/*                              analysis()                                   */
/*****************************************************************************/
void
analysis( float *in, int sx, int sy, int lenx, int leny, int hor, int ver,
          int filterType )
{
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
  int x, y, pos, ipos, tap;
  float *lpf, *hpf;
  float *ibuf, *obuf, *buf1, *buf2;

  tap = filter_coeff( &lpf, &hpf, filterType, 0 );      /* filter */

  ibuf = ( float * )calloc( lenx * leny, sizeof( float ) );
  obuf = ( float * )calloc( lenx * leny, sizeof( float ) );
  buf1 = ( float * )calloc( leny, sizeof( float ) );
  buf2 = ( float * )calloc( leny, sizeof( float ) );

  /* in--> ibuf */
  /* since inframe is filtered 4 times */
  /* output is saved in inframe, so we need buffer */
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      pos = y * lenx + x;
      ipos = ( y + sy ) * hor + ( x + sx );
      ibuf[pos] = in[ipos];     /* int -> float */
    }
  }

  /* LL */
  filter_2d( obuf, ibuf, buf1, buf2, lenx, leny, lpf, hpf, tap, LL );
  for( y = 0; y < leny; y += 2 ) {      /* subsampling */
    for( x = 0; x < lenx; x += 2 ) {
      pos = y * lenx + x;
      ipos = ( y / 2 + sy ) * hor + ( x / 2 + sx );
      in[ipos] = obuf[pos];
    }
  }

  // goto stop;
  
  /* HL */
  filter_2d( obuf, ibuf, buf1, buf2, lenx, leny, lpf, hpf, tap, HL );
  for( y = 0; y < leny; y += 2 ) {      /* subsampling */
    for( x = 1; x < lenx; x += 2 ) {
      pos = y * lenx + x;
      ipos = ( y / 2 + sy ) * hor + ( x / 2 + sx + lenx / 2 );
      in[ipos] = obuf[pos];
    }
  }

  /* LH */
  filter_2d( obuf, ibuf, buf1, buf2, lenx, leny, lpf, hpf, tap, LH );
  for( y = 1; y < leny; y += 2 ) {      /* subsampling */
    for( x = 0; x < lenx; x += 2 ) {
      pos = y * lenx + x;
      ipos = ( y / 2 + sy + leny / 2 ) * hor + ( x / 2 + sx );
      in[ipos] = obuf[pos];
    }
  }

  /* HH */
  filter_2d( obuf, ibuf, buf1, buf2, lenx, leny, lpf, hpf, tap, HH );
  for( y = 1; y < leny; y += 2 ) {      /* subsampling */
    for( x = 1; x < lenx; x += 2 ) {
      pos = y * lenx + x;
      ipos = ( y / 2 + sy + leny / 2 ) * hor + ( x / 2 + sx + lenx / 2 );
      in[ipos] = obuf[pos];
    }
  }

  // stop: printf("only do LL filtering (subband.c)\n");

  free( ibuf );
  free( obuf );
  free( buf1 );
  free( buf2 );
  free( lpf );
  free( hpf );

  return;
}


/*****************************************************************************/
/*                              synthesis()                                  */
/*****************************************************************************/
void
synthesis( float *in, int sx, int sy, int lenx, int leny, int hor, int ver,
           int filterType )
{
  /*************************************************/
  /*            in              out                */
  /*        -----------     -----------            */
  /*        |    |    |     |         |            */
  /*        | LL | HL |     |         |            */
  /*        -----------     |         |            */
  /*        |    |    |     |         |            */
  /*        | LH | HH |     |         |            */
  /*        -----------     -----------            */
  /*************************************************/

  int x, y, pos, ipos, tap;
  float *lpf, *hpf;
  float *ibuf, *obuf, *tbuf, *buf1, *buf2, *buf3;
  float *fptr1, *fptr2;

  tap = filter_coeff( &lpf, &hpf, filterType, 1 );      /* filter coefficient */

  ibuf = ( float * )calloc( lenx * leny, sizeof( float ) );
  obuf = ( float * )calloc( lenx * leny, sizeof( float ) );
  tbuf = ( float * )calloc( lenx * leny, sizeof( float ) );
  buf1 = ( float * )calloc( lenx + tap - 1, sizeof( float ) );
  buf2 = ( float * )calloc( leny + tap - 1, sizeof( float ) );
  buf3 = ( float * )calloc( leny, sizeof( float ) );

  /* LL */
  interpol_2d( obuf, in, ibuf, buf1, buf2, buf3, lpf, hpf, tap, lenx, leny,
               sx, sy, hor, LL );
  fptr1 = tbuf;
  fptr2 = obuf;
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      *fptr1++ = 4 * ( *fptr2++ );
    }
  }

  /* HL */
  interpol_2d( obuf, in, ibuf, buf1, buf2, buf3, lpf, hpf, tap, lenx, leny,
               sx, sy, hor, HL );
  fptr1 = tbuf;
  fptr2 = obuf;
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      *fptr1++ += 4 * ( *fptr2++ );
    }
  }

  /* LH */
  interpol_2d( obuf, in, ibuf, buf1, buf2, buf3, lpf, hpf, tap, lenx, leny,
               sx, sy, hor, LH );
  fptr1 = tbuf;
  fptr2 = obuf;
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      *fptr1++ += 4 * ( *fptr2++ );
    }
  }

  /* HH */
  interpol_2d( obuf, in, ibuf, buf1, buf2, buf3, lpf, hpf, tap, lenx, leny,
               sx, sy, hor, HH );
  fptr1 = tbuf;
  fptr2 = obuf;
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      *fptr1++ += 4 * ( *fptr2++ );
    }
  }

  /* tbuf --> in */
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      pos = y * lenx + x;
      ipos = ( y + sy ) * hor + ( x + sx );
      in[ipos] = tbuf[pos];
    }
  }

  free( ibuf );
  free( obuf );
  free( tbuf );
  free( buf1 );
  free( buf2 );
  free( buf3 );
  free( lpf );
  free( hpf );

  return;
}

/****************************************************************************/
/*                            analysis_child()                              */
/****************************************************************************/
void
analysis_child( float *frame, int sx, int sy, int lenx, int leny, int hor,
                int ver, int filterType )
{

  analysis( frame, sx, sy, lenx / 2, leny / 2, hor, ver, filterType );
  analysis( frame, sx + lenx / 2, sy, lenx / 2, leny / 2, hor, ver,
            filterType );
  analysis( frame, sx, sy + leny / 2, lenx / 2, leny / 2, hor, ver,
            filterType );
  analysis( frame, sx + lenx / 2, sy + leny / 2, lenx / 2, leny / 2, hor, ver,
            filterType );
}

/****************************************************************************/
/*                            synthesis_child()                             */
/****************************************************************************/
void
synthesis_child( float *frame, int sx, int sy, int lenx, int leny, int hor,
                 int ver, int filterType )
{
  synthesis( frame, sx, sy, lenx / 2, leny / 2, hor, ver, filterType );
  synthesis( frame, sx + lenx / 2, sy, lenx / 2, leny / 2, hor, ver,
             filterType );
  synthesis( frame, sx, sy + leny / 2, lenx / 2, leny / 2, hor, ver,
             filterType );
  synthesis( frame, sx + lenx / 2, sy + leny / 2, lenx / 2, leny / 2, hor,
             ver, filterType );
}


/*****************************************************************************/
/*                            filter_2d()                                    */
/*****************************************************************************/
void
filter_2d( float *obuf, float *ibuf, float *buf1, float *buf2, int lenx,
           int leny, float *lpf, float *hpf, int tap, int mode )
{
  int x, y;

  /* horizontal direction */
  for( y = 0; y < leny; y++ ) {
    if( mode == LL || mode == LH )
      filter_1d( &obuf[y * lenx], &ibuf[y * lenx], lenx, lpf, tap, 0 );
    else
      filter_1d( &obuf[y * lenx], &ibuf[y * lenx], lenx, hpf, tap, 1 );
  }

  /* vertical direction */
  for( x = 0; x < lenx; x++ ) {
    for( y = 0; y < leny; y++ )
      buf1[y] = obuf[y * lenx + x];
    if( mode == LL || mode == HL )
      filter_1d( buf2, buf1, leny, lpf, tap, 0 );
    else
      filter_1d( buf2, buf1, leny, hpf, tap, 1 );
    for( y = 0; y < leny; y++ )
      obuf[y * lenx + x] = buf2[y];
  }
}



/*****************************************************************************/
/*                           interpol_2d()                                   */
/*****************************************************************************/
void
interpol_2d( float *obuf, float *in, float *ibuf, float *buf1, float *buf2,
             float *buf3, float *lpf, float *hpf, int tap, int lenx, int leny,
             int sx, int sy, int hor, int mode )
{
  int x, y;

  /* horizontal interpolation and filtering */
  for( y = 0; y < leny / 2; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      buf1[x + tap / 2] = 0;
      if( mode == LL && !( x % 2 ) )
        buf1[x + tap / 2] = in[( y + sy ) * hor + ( x / 2 + sx )];
      else if( mode == HL && ( x % 2 ) )
        buf1[x + tap / 2] = in[( y + sy ) * hor + ( x / 2 + sx + lenx / 2 )];
      else if( mode == LH && !( x % 2 ) )
        buf1[x + tap / 2] = in[( y + sy + leny / 2 ) * hor + ( x / 2 + sx )];
      else if( mode == HH && ( x % 2 ) )
        buf1[x + tap / 2] =
          in[( y + sy + leny / 2 ) * hor + ( x / 2 + sx + lenx / 2 )];
    }
    if( mode == LL || mode == LH )
      filter_1d_normal( &ibuf[y * lenx], buf1 + tap / 2, lenx, lpf, tap, 0 );
    else
      filter_1d_normal( &ibuf[y * lenx], buf1 + tap / 2, lenx, hpf, tap, 1 );
  }

  /* vertical interpolation and filtering */
  for( x = 0; x < lenx; x++ ) {
    for( y = 0; y < leny; y++ ) {
      buf2[y + tap / 2] = 0;
      if( mode == LL || mode == HL ) {
        if( !( y % 2 ) )
          buf2[y + tap / 2] = ibuf[( y / 2 ) * lenx + x];
      } else {
        if( ( y % 2 ) )
          buf2[y + tap / 2] = ibuf[( y / 2 ) * lenx + x];
      }
    }
    if( mode == LL || mode == HL )
      filter_1d_normal( buf3, buf2 + tap / 2, leny, lpf, tap, 0 );
    else
      filter_1d_normal( buf3, buf2 + tap / 2, leny, hpf, tap, 1 );

    for( y = 0; y < leny; y++ )
      obuf[y * lenx + x] = buf3[y];
  }
}

/*****************************************************************************/
/*                            filter_1d()                                    */
/*****************************************************************************/
void
filter_1d( float *out, float *in, int length, float *filter, int tap,
           int start )
{
  int i, j, k;
  float tempf;

  /* 1 dimensional filtering by convolution */
  /* the symmetric extension for the boundary data */
  /* Example, for filter tap = 7 and data length = 10 */
  /*    -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 */
  /*     3  2  1 0 1 2 3 4 5 6 7 8 9 8  7  6 */
  /*                                 18 18 18 */

  for( i = start; i < length; i += 2 ) {        /* considering the subsampling process */
    tempf = 0.;
    for( j = 0; j < tap; j++ ) {
      k = i - j + tap / 2;
      if( k < 0 )
        k = -k;
      else if( k >= length )
        k = ( length - 1 ) * 2 - k;

      if( k < 0 || k >= length ) {
        printf( "output range k %d\n", k );
        exit( 0 );
      }
      tempf += in[k] * filter[j];
    }
    out[i] = tempf;
  }
}

/*****************************************************************************/
/*                         filter_1d_normal()                                */
/*****************************************************************************/
void
filter_1d_normal( float *out, float *in, int length, float *filter, int tap,
                  int start )
{
  int i, j, jstart;
  float tempf, *fptr1, *fptr2;

  /* 1 dimensional filtering by convolution */
  /* the symmetric extension for the boundary data */
  /* Example, for filter tap = 7 and data length = 10 */
  /*    -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 */
  /*     3  2  1 0 1 2 3 4 5 6 7 8 9 8  7  6 */
  /*                                 18 18 18 */

  /* reflection for boundary */
  for( i = 1; i <= tap / 2; i++ ) {
    in[-i] = in[i];
    in[length - 1 + i] = in[length - 1 - i];
  }

  /* we don't have to multiply zeros which is filled for upsampling */
  for( i = 0; i < length; i++ ) {
    tempf = 0.;

    fptr1 = &in[i - tap / 2];
    fptr2 = filter;
    jstart = 0;
    if( ( ( i % 2 ) && !start ) || ( !( i % 2 ) && start ) ) {
      fptr1++;
      fptr2++;
      jstart = 1;
    }

    for( j = jstart; j < tap; j += 2 ) {
      tempf += ( *fptr1 ) * ( *fptr2 );
      fptr1 += 2;
      fptr2 += 2;
    }
    out[i] = tempf;
  }
}

/*****************************************************************************/
/*                              filter_coeff()                               */
/*****************************************************************************/
int
filter_coeff( float **lpf, float **hpf, int FILTER, int analsyn )
{
  int i, tap;
  float *h0, *h1, *g0, *g1;
  FILE *fp;

  if( FILTER == 0 ) {
    tap = 9;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );

    h0[4] = ( float )0.5645751; /* Adelson and Simoncelli */
    h0[3] = h0[5] = ( float )0.2927051; /* 9 tap filter */
    h0[2] = h0[6] = ( float )-0.05224239;       /* Senoo and Girod '92 paper */
    h0[1] = h0[7] = ( float )-0.04270508;
    h0[0] = h0[8] = ( float )0.01995484;

    tap = 9;
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    h1[4] = ( float )0.5645751;
    h1[3] = h1[5] = ( float )-0.2927051;
    h1[2] = h1[6] = ( float )-0.05224239;
    h1[1] = h1[7] = ( float )0.04270508;
    h1[0] = h1[8] = ( float )0.01995484;
  } else if( FILTER == 1 ) {
    tap = 9;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );

    h0[4] = ( float )0.635577;  /* Nguyen and Vaidyanathan */
    h0[3] = h0[5] = ( float )0.307180;
    h0[2] = h0[6] = ( float )-0.164740;
    h0[1] = h0[7] = ( float )0.002902;
    h0[0] = h0[8] = ( float )0.036870;

    tap = 9;
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    h1[4] = ( float )0.490663;
    h1[3] = h1[5] = ( float )-0.309311;
    h1[2] = h1[6] = ( float )0.004668;
    h1[1] = h1[7] = ( float )0.059311;
    h1[0] = h1[8] = ( float )0.0;

  } else if( FILTER == 2 ) {    /* Daubechies */
    tap = 9;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );

    h0[4] = ( float )0.602949;
    h0[3] = h0[5] = ( float )0.266864;
    h0[2] = h0[6] = ( float )-0.078223;
    h0[1] = h0[7] = ( float )-0.016864;
    h0[0] = h0[8] = ( float )0.026749;

    tap = 9;
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    h1[4] = ( float )0.557543;
    h1[3] = h1[5] = ( float )-0.295636;
    h1[2] = h1[6] = ( float )-0.028772;
    h1[1] = h1[7] = ( float )0.045636;
    h1[0] = h1[8] = ( float )0.0;
  } else if( FILTER == 3 ) {    /* MPEG half band filter */
    tap = 9;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );

    h0[4] = ( float )138. / 256;
    h0[3] = h0[5] = ( float )88. / 256;
    h0[2] = h0[6] = ( float )0;
    h0[1] = h0[7] = ( float )-29. / 256;
    h0[0] = h0[8] = ( float )0.;

    tap = 9;
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    h1[4] = ( float )0.;
    h1[3] = h1[5] = ( float )0.;
    h1[2] = h1[6] = ( float )0.;
    h1[1] = h1[7] = ( float )0.;
    h1[0] = h1[8] = ( float )0.0;
  } else if( FILTER == 4 ) {    /* new MPEG half band filter */
    tap = 13;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );

    h0[6] = ( float )26. / 64;
    h0[5] = h0[7] = ( float )19. / 64;
    h0[4] = h0[8] = ( float )5. / 64;
    h0[3] = h0[9] = ( float )-3. / 64;
    h0[2] = h0[10] = ( float )-4. / 64;
    h0[1] = h0[11] = ( float )0.;
    h0[0] = h0[12] = ( float )2. / 64;

    tap = 13;
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    h1[6] = ( float )0.;
    h1[5] = h1[7] = ( float )0.;
    h1[4] = h1[8] = ( float )0.;
    h1[3] = h1[9] = ( float )0.;
    h1[2] = h1[10] = ( float )0.;
    h1[1] = h1[11] = ( float )0.;
    h1[0] = h1[12] = ( float )0.;
  } else if( FILTER == 5 ) {    /* Daubechies  exchange the analytic filter and synthetic filter */
    tap = 9;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );

    h0[4] = ( float )0.557543;
    h0[3] = h0[5] = ( float )0.295636;
    h0[2] = h0[6] = ( float )-0.028772;
    h0[1] = h0[7] = ( float )-0.045636;
    h0[0] = h0[8] = ( float )0.0;

    tap = 9;
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    h1[4] = ( float )0.602949;
    h1[3] = h1[5] = ( float )-0.266864;
    h1[2] = h1[6] = ( float )-0.078223;
    h1[1] = h1[7] = ( float )0.016864;
    h1[0] = h1[8] = ( float )0.026749;

  } else if( FILTER == 6 ) {
    tap = 33;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    if( ( fp = fopen( "z:\\3d_ezbc\\Filter\\halfband.txt", "rb" ) ) == NULL ) {
      printf( "can not open z:\\3d_ezbc\\Filter\\halfband.txt\n" );
      exit( 1 );
    }
    for( i = 0; i < tap; i++ ) {
      fscanf( fp, "%f\n", &h0[i] );
    }
    fclose( fp );

  } else if( FILTER == 7 ) {    /* Daubechies - scaled to sqrt(2) */
    // printf ("use Daubechies spatial filters ( %d )\n", FILTER);
    tap = 9;
    h0 = ( float * )getarray( tap, sizeof( float ), "h0" );

    h0[4] = ( float )0.852699;
    h0[3] = h0[5] = ( float ) 0.377403;
    h0[2] = h0[6] = ( float )-0.110624;
    h0[1] = h0[7] = ( float )-0.023849;
    h0[0] = h0[8] = ( float ) 0.037829;

    tap = 9;
    h1 = ( float * )getarray( tap, sizeof( float ), "h1" );

    h1[4] = ( float )0.788485;
    h1[3] = h1[5] = ( float )-0.418092;
    h1[2] = h1[6] = ( float )-0.040690;
    h1[1] = h1[7] = ( float ) 0.064539;
    h1[0] = h1[8] = ( float ) 0.0;
  } else {
    printf( "error in filter_coeff()\n" );
    exit( 1 );
  }

  if( !analsyn ) {
    *lpf = h0;
    *hpf = h1;
  } else {
    if( FILTER == 3 ) {
      tap = 9;
      g0 = ( float * )getarray( tap, sizeof( float ), "g0" );

      g0[4] = ( float )256. / 256 / 2;
      g0[3] = g0[5] = ( float )140. / 256 / 2;
      g0[2] = g0[6] = ( float )0.;
      g0[1] = g0[7] = ( float )-12. / 256 / 2;
      g0[0] = g0[8] = ( float )0;

      tap = 9;
      g1 = ( float * )getarray( tap, sizeof( float ), "g1" );

      g1[4] = h1[4];
      g1[3] = g1[5] = -h1[3];
      g1[2] = g1[6] = h1[2];
      g1[1] = g1[7] = -h1[1];
      g1[0] = g1[8] = h1[0];
    } else if( FILTER == 4 ) {
      printf( "not available\n" );
      exit( 1 );
    }  else if( FILTER == 7 ) {    /* Daubechies - scaled to sqrt(2) */
      // printf ("use Daubechies spatial filters ( %d )\n", FILTER);
      /*
      tap = 9;
      g0 = ( float * )getarray( tap, sizeof( float ), "g0" );
      
      g0[4] = ( float )0.852699;
      g0[3] = g0[5] = ( float )0.418092;
      g0[2] = g0[6] = ( float )-0.110624;
      g0[1] = g0[7] = ( float )-0.064539;
      g0[0] = g0[8] = ( float )0.037829;
      
      tap = 9;
      g1 = ( float * )getarray( tap, sizeof( float ), "g1" );
      
      g1[4] = ( float )0.788485;
      g1[3] = g1[5] = ( float )-0.377403;
      g1[2] = g1[6] = ( float )-0.040690;
      g1[1] = g1[7] = ( float )0.023849;
      g1[0] = g1[8] = ( float )0.0;
      */
      float alpha = 2.0;
      tap = 9;
      g0 = ( float * )getarray( tap, sizeof( float ), "g0" );
      
      g0[4] = h1[4] / alpha;
      g0[3] = g0[5] = -h1[3] / alpha;
      g0[2] = g0[6] =  h1[2] / alpha;
      g0[1] = g0[7] = -h1[1] / alpha;
      g0[0] = g0[8] =  h1[0] / alpha;
      
      tap = 9;
      g1 = ( float * )getarray( tap, sizeof( float ), "g1" );
      
      g1[4] = h0[4] / alpha;
      g1[3] = g1[5] = -h0[3] / alpha;
      g1[2] = g1[6] =  h0[2] / alpha;
      g1[1] = g1[7] = -h0[1] / alpha;
      g1[0] = g1[8] =  h0[0] / alpha;
    } else {
      // printf ("use inverted reconstruction filters ( %d )\n", FILTER);
      tap = 9;
      g0 = ( float * )getarray( tap, sizeof( float ), "g0" );
      
      g0[4] = h1[4];
      g0[3] = g0[5] = -h1[3];
      g0[2] = g0[6] = h1[2];
      g0[1] = g0[7] = -h1[1];
      g0[0] = g0[8] = h1[0];

      tap = 9;
      g1 = ( float * )getarray( tap, sizeof( float ), "g1" );

      g1[4] = h0[4];
      g1[3] = g1[5] = -h0[3];
      g1[2] = g1[6] = h0[2];
      g1[1] = g1[7] = -h0[1];
      g1[0] = g1[8] = h0[0];
    }
    *lpf = g0;
    *hpf = g1;

    free( h0 );
    free( h1 );
  }

  return tap;
}

/****************************************************************************/
/*                          analysis_1frame()                               */
/****************************************************************************/
void
analysis_1frame( float *frame, int band, int hor, int ver, int filterType )
{

  /******************************/
  /*      -----------------     */
  /*      | 0 | 1 | 4 | 5 |     */
  /*      -----------------     */
  /*      | 2 | 3 | 6 | 7 |     */
  /*      -----------------     */
  /*      | 8 | 9 | 12| 13|     */
  /*      -----------------     */
  /*      | 10| 11| 14| 15|     */
  /*      -----------------     */
  /******************************/

  /**********************************/
  /*      ---------------------     */
  /*      |0 1 |    |         |     */
  /*      |2 3 |  4 |         |     */
  /*      -----------    7    |     */
  /*      |    |    |         |     */
  /*      | 5  |  6 |         |     */
  /*      ---------------------     */
  /*      |         |         |     */
  /*      |         |         |     */
  /*      |    8    |    9    |     */
  /*      |         |         |     */
  /*      |         |         |     */
  /*      ---------------------     */
  /**********************************/


  switch ( band ) {
  default:
    printf( "error in analysis_1frame()\n" );
    exit( 1 );
    break;
  case 4:
    analysis( frame, 0, 0, hor, ver, hor, ver, filterType );
    break;
  case 7:
    analysis( frame, 0, 0, hor, ver, hor, ver, filterType );
    analysis( frame, 0, 0, hor / 2, ver / 2, hor, ver, filterType );
    break;
  case 10:
    analysis( frame, 0, 0, hor, ver, hor, ver, filterType );
    analysis( frame, 0, 0, hor / 2, ver / 2, hor, ver, filterType );
    analysis( frame, 0, 0, hor / 4, ver / 4, hor, ver, filterType );
    break;
  case 16:
    analysis( frame, 0, 0, hor, ver, hor, ver, filterType );
    analysis_child( frame, 0, 0, hor, ver, hor, ver, filterType );
    break;
  }
}


/****************************************************************************/
/*                          synthesis_1frame()                               */
/****************************************************************************/
void
synthesis_1frame( float *frame, int band, int hor, int ver, int filterType )
{

  switch ( band ) {
  default:
    printf( "error in synthesis_1frame()\n" );
    exit( 1 );
    break;
  case 4:
    synthesis( frame, 0, 0, hor, ver, hor, ver, filterType );
    break;
  case 7:
    synthesis( frame, 0, 0, hor / 2, ver / 2, hor, ver, filterType );
    synthesis( frame, 0, 0, hor, ver, hor, ver, filterType );
    break;
  case 10:
    synthesis( frame, 0, 0, hor / 4, ver / 4, hor, ver, filterType );
    synthesis( frame, 0, 0, hor / 2, ver / 2, hor, ver, filterType );
    synthesis( frame, 0, 0, hor, ver, hor, ver, filterType );
    break;
  case 16:
    synthesis_child( frame, 0, 0, hor, ver, hor, ver, filterType );
    synthesis( frame, 0, 0, hor, ver, hor, ver, filterType );
    break;
  }
}
