
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//                     I M A G E   C L A S S

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//           > > > >    C++ version 11.06 -  02/01/96   < < < <

// Amir Said - amir@densis.fee.unicamp.br
// University of Campinas (UNICAMP)
// Campinas, SP 13081, Brazil

// William A. Pearlman - pearlman@ecse.rpi.edu
// Rensselaer Polytechnic Institute
// Troy, NY 12180, USA

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// Copyright (c) 1995, 1996 Amir Said & William A. Pearlman

// This program is Copyright (c) by Amir Said & William A. Pearlman.
// It may not be redistributed without the consent of the copyright
// holders. In no circumstances may the copyright notice be removed.
// The program may not be sold for profit nor may they be incorporated
// in commercial programs without the written permission of the copyright
// holders. This program is provided as is, without any express or
// implied warranty, without even the warranty of fitness for a
// particular purpose.

// - - Inclusion - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include <string.h>
#include "general.h"
#include "image_bw.h"
#include "structN.h"
#include <math.h>

// - - Constants - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static char *M_MSG = "< Image_SP >";

static char *R_MSG = "< Image_BW > cannot read from file";

static char *W_MSG = "< Image_BW > cannot write to file";

static char *L_MSG = "< Image_BW > larger than specified dimension";


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//  Auxiliary functions

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//ostream & operator << ( ostream & out, Image_Coord & in_ord ) {
//  cout << "( " << in_ord.x << ", " << in_ord.y << " )";
// return out;
//}

#ifdef LOSSLESS

static void
SP_Transform( int m, int in[], int l[], int h[] )
{
  int i, k, d1, d2, mm = m - 1;

  for( i = k = 0; i < m; i++, k += 2 ) {
    l[i] = ( in[k] + in[k + 1] ) >> 1;
    h[i] = in[k] - in[k + 1];
  }

  h[0] -= ( d2 = l[0] - l[1] ) >> 2;
  for( i = 1; i < mm; i++ ) {
    d1 = d2;
    d2 = l[i] - l[i + 1];
    h[i] -= ( ( ( d1 + d2 - h[i + 1] ) << 1 ) + d2 + 3 ) >> 3;
  }
  h[i] -= d2 >> 2;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void
SP_Recover( int m, int l[], int h[], int out[] )
{
  int i, k, d1, d2, t;

  t = ( h[m - 1] += ( d1 = l[m - 2] - l[m - 1] ) >> 2 );
  for( i = m - 2; i > 0; i-- ) {
    d2 = d1;
    d1 = l[i - 1] - l[i];
    t = ( h[i] += ( ( ( d1 + d2 - t ) << 1 ) + d2 + 3 ) >> 3 );
  }
  h[0] += d1 >> 2;

  for( i = k = 0; i < m; i++, k += 2 ) {
    out[k] = l[i] + ( ( h[i] + 1 ) >> 1 );
    out[k + 1] = out[k] - h[i];
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#else
     // 下面是进行小波变换的系数
static const float SmoothingFactor = ( float )0.8;

static const int NumbTap = 4;

static const float T_LowPass[5] =
  { ( float )0.852699, ( float )0.377403, ( float )-0.110624,
( float )-0.023849, ( float )0.037829 };

static const float T_HighPass[5] =
  { ( float )0.788485, ( float )-0.418092, ( float )-0.040690,
( float )0.064539, ( float )0.0 };

static const float R_LowPass[5] =
  { ( float )0.852699, ( float )0.418092, ( float )-0.110624,
( float )-0.064539, ( float )0.037829 };

static const float R_HighPass[5] =
  { ( float )0.788485, ( float )-0.377403, ( float )-0.040690,
( float )0.023849, ( float )0.0 };

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline float
Filter_L( const float *f, float *v )
{
  return f[0] * v[0] +
    f[1] * ( v[1] + v[-1] ) + f[2] * ( v[2] + v[-2] ) +
    f[3] * ( v[3] + v[-3] ) + f[4] * ( v[4] + v[-4] );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline float
Filter_H( const float *f, float *v )
{
  return f[0] * v[0] +
    f[1] * ( v[1] + v[-1] ) + f[2] * ( v[2] + v[-2] ) +
    f[3] * ( v[3] + v[-3] );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void
Reflection( float *h, float *t )
{
  for( int i = 1; i <= NumbTap; i++ ) {
    h[-i] = h[i];
    t[i] = t[-i];
  }
}

#endif

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//  Member functions of the class  < Image_BW >

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// - - Private functions - - - - - - - - - - - - - - - - - - - - - - - -

int
Image_BW::max_levels( int n ) // 能进行二维小波分解的最大的次数
{
  int l1, l2;
  for( l1 = 0; !( n & 1 ); l1++ ) // 最小的n的偶数，表示能进行的最大的二维小波分解的次数
    n >>= 1;
#ifdef LOSSLESS
  for( l2 = l1 - 3; n; l2++ )
    n >>= 1;
#else
  for( l2 = l1 - 4; n; l2++ )
    n >>= 1;                    // this makes sure the size of the lowest subband will be larger or equal to 8
#endif
  return ( l1 < l2 ? l1 : l2 );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::assign_mem( Image_Coord d, int b )
{
  if( ( b < 1 ) || ( b > 2 ) )
    Error( "Invalid number of < Image_BW > bytes" );
  if( ( levels >= 0 ) && ( dim.x == d.x ) && ( dim.y == d.y ) )
    return;
  free_mem(  );
  if( ( d.x < 64 ) || ( d.y < 64 ) ) {
    printf( "d.x = %d, d.y = %d\n", d.x, d.y );
    Error( "< Image_BW > dimension is too small or negative" );
  }
  dim = d;

#ifdef EXPAND
  // expand image to garantee the number of spatial decomposition levels (refering to max_levels())
  pdim.x = ( (d.x < 256) ? ( d.x + 7 ) & 0x3FF8 : ( d.x + 15 ) & 0x3FF0 );
  pdim.y = ( (d.y < 256) ? ( d.y + 7 ) & 0x3FF8 : ( d.y + 15 ) & 0x3FF0 );
#else
  pdim = dim;
#endif




  NEW_VECTOR( coeff, pdim.x, Pel_Type *, M_MSG );
  for( int i = 0; i < pdim.x; i++ ) {
    NEW_VECTOR( coeff[i], pdim.y, Pel_Type, M_MSG );
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::free_mem( void )
{
  if( levels >= 0 ) {
    for( int i = pdim.x - 1; i >= 0; i-- )
      delete[]coeff[i];
    delete[]coeff;
  }
  bytes = dim.x = dim.y = 0;
  levels = -1;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::extend( void )
{
  int i, j;
  for( j = dim.y - 1; j < pdim.y - 1; j++ ) {
    coeff[0][j + 1] = ( coeff[0][j] + coeff[1][j] ) / 2;
    coeff[dim.x - 1][j + 1] =
      ( coeff[dim.x - 1][j] + coeff[dim.x - 2][j] ) / 2;
    for( i = dim.x - 2; i > 0; i-- )
      coeff[i][j + 1] =
        ( coeff[i - 1][j] + coeff[i][j] + coeff[i + 1][j] ) / 3;
  }
  for( i = dim.x - 1; i < pdim.x - 1; i++ ) {
    coeff[i + 1][0] = ( coeff[i][0] + coeff[i][1] ) / 2;
    coeff[i + 1][pdim.y - 1] =
      ( coeff[i][pdim.y - 1] + coeff[i][pdim.y - 2] ) / 2;
    for( j = pdim.y - 2; j > 0; j-- )
      coeff[i + 1][j] =
        ( coeff[i][j - 1] + coeff[i][j] + coeff[i][j + 1] ) / 3;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - Public functions  - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::read_pic( Image_Coord d, char *file_name, int b )
{
  assign_mem( d, b );
  mean = levels = 0;
  bytes = b;

  FILE *in_file = fopen( file_name, "rb" );
  if( in_file == NULL )
    Error( R_MSG );

  int i, j, k, p, c;
  for( i = 0; i < dim.x; i++ )
    for( j = 0; j < dim.y; j++ ) {
      for( p = k = 0; k < bytes; k++ ) {
        if( ( c = getc( in_file ) ) == EOF )
          Error( R_MSG );
        p = ( p << 8 ) | c;
      }
      coeff[i][j] = ( float )p;
    }
  if( getc( in_file ) != EOF )
    Error( L_MSG );
  fclose( in_file );

  extend(  );
}

void
Image_BW::read_float_image( Image_Coord d, float *image )
{
  assign_mem( d, 1 );
  mean = levels = 0;
  bytes = 1;

  int i, j, k = 0;
  for( i = 0; i < dim.x; i++ )
    for( j = 0; j < dim.y; j++ )
      //coeff[i][j] = *image++;
      coeff[i][j] = image[k++];
  
  extend(  );

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

float
Image_BW::compare( char *file_name )
{
  if( levels )
    Error( "cannot compare < Image_BW >" );

  FILE *in_file = fopen( file_name, "rb" );
  if( in_file == NULL )
    Error( R_MSG );

  double mse = 0.0;
  int i, j, k, p, c, t;
  for( i = 0; i < dim.x; i++ )
    for( j = 0; j < dim.y; j++ ) {
#ifdef LOSSLESS
      t = coeff[i][j];
#else
      t = int ( floor( 0.499 + coeff[i][j] ) );
#endif
      if( t < 0 )
        t = 0;
      if( ( bytes == 1 ) && ( t > 255 ) )
        t = 255;
      for( p = k = 0; k < bytes; k++ ) {
        if( ( c = getc( in_file ) ) == EOF )
          Error( R_MSG );
        p = ( p << 8 ) | c;
      }
      mse += Sqr( ( float )( p - t ) );
    }
  if( getc( in_file ) != EOF )
    Error( L_MSG );
  fclose( in_file );

  return ( float )( ( mse / dim.x ) / dim.y );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::write_pic( char *file_name )
{
  if( levels )
    Error( "cannot write < Image_BW >" );

  FILE *out_file = fopen( file_name, "wb" );
  if( out_file == NULL )
    Error( W_MSG );

  int i, j, k;
  for( i = 0; i < dim.x; i++ )
    for( j = 0; j < dim.y; j++ ) {
#ifdef LOSSLESS
      k = coeff[i][j];
#else
      k = int ( floor( 0.499 + coeff[i][j] ) );
#endif
      if( k < 0 )
        k = 0;
      if( bytes == 2 ) {
        if( putc( k >> 8, out_file ) == EOF )
          Error( W_MSG );
      } else {
        if( k > 255 )
          k = 255;
      }
      if( putc( k & 0xFF, out_file ) == EOF )
        Error( W_MSG );
    }

  fclose( out_file );
}

void
Image_BW::write_float_image( float *image )
{
  int i, j;

  for( i = 0; i < dim.x; i++ )
    for( j = 0; j < dim.y; j++ )
      *image++ = coeff[i][j];


}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::reset( Image_Coord d )
{
  assign_mem( d, 1 );
  bytes = 1;
  mean = shift = smoothing = 0;



  levels = Min( max_levels( pdim.x ), max_levels( pdim.y ) );

  int i, j;
  for( i = 0; i < pdim.x; i++ )
    for( j = 0; j < pdim.y; j++ )
      coeff[i][j] = 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::reset( int m, int shf, int smt )
{
  mean = m;
  shift = shf;
  smoothing = smt;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::reset( Image_Coord d, int b, int m, int shf, int smt )
{
  reset( d );
  bytes = b;
  mean = m;
  shift = shf;
  smoothing = smt;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::transform( int sm_numb ) // sm_numb 平滑系数
{
  if( levels )
    Error( "cannot transform < Image_BW >" );
  if( ( dim.x < 32 ) || ( dim.y < 32 ) )
    Error( "< Image_BW > is too small to transform" );

  float sum;

  // Chronometer cpu_time;
  // cpu_time.start("\n  Starting image transformation...");
  if(high_band == NO){
	levels = Min( max_levels( pdim.x ), max_levels( pdim.y ) ); //levels = 5 for luma, levels = 4 for chroma bands
//	printf("LOW band!\n");
  }
  else{
	assert(high_band == YES);
	levels = 2;
//	printf("HIGH band!\n");
  }

  if( levels <= 0 )
    Error( "invalid < Image_BW > dimension" );

  smoothing = sm_numb;
  if( ( sm_numb < 0 ) || ( sm_numb > 7 ) )
    Error( "invalid < Image_BW > smoothing factor" );

#ifdef LOSSLESS
  int i = 0, j = Max( pdim.x, pdim.y ), k = j << 1;
#else
  int i = NumbTap, j = 0, k = Max( pdim.x, pdim.y ) + ( i << 1 );
  //float sm_mt;
#endif
  CREATE_VECTOR( temp_line, k, Pel_Type, M_MSG );
  Pel_Type *t, *in_line = temp_line + i; // *out_line = in_line + j;

// hierarchical wavelet or S+P transformation

  int lv, nx, ny, mx = pdim.x, my = pdim.y;
//  printf("pdimx = %d, pdimy = %d, dimx = %d, dimy = %d\n",pdim.x,pdim.y,dim.x,dim.y);

  //printf("subband decomposition levels = %d(image_bw.cpp)\n", levels);

  for( lv = 0; lv < levels; lv++ ) { // 子带进行几次分解

    // shifts are halved, multiplier is updated

    nx = mx;
    mx >>= 1;
    ny = my;
    my >>= 1;

#ifndef LOSSLESS
    float sm_mt = 1 + smoothing * SmoothingFactor / ( 2 + lv * lv );
#endif
    // transformation of columns

    for( j = 0; j < ny; j++ ) 
	{
      for( i = 0; i < nx; i++ )
        in_line[i] = coeff[i][j];  // 读取一列系数
#ifdef LOSSLESS
      SP_Transform( mx, in_line, out_line, out_line + mx );
      for( i = 0; i < nx; i++ )
        coeff[i][j] = out_line[i];
    }
#else

#ifdef Haar
      for( i = 0; i < mx; i++ ) {
        coeff[i][j] =
          sm_mt * ( in_line[i * 2] + in_line[i * 2 + 1] ) / sqrt( 2. );
        coeff[i + mx][j] =
          ( in_line[i * 2] - in_line[i * 2 + 1] ) / sqrt( 2. );
      }
    }
#else
      Reflection( in_line, in_line + nx - 1 );
      for( i = 0, t = in_line; i < mx; i++ ) {
        coeff[i][j] = sm_mt * Filter_L( T_LowPass, t++ );
        coeff[i + mx][j] = Filter_H( T_HighPass, t++ );
      }
    }
#endif

#endif

    // transformation of rows

    for( i = 0; i < nx; i++ ) {
      memcpy( in_line, coeff[i], ny * sizeof( Pel_Type ) );
#ifdef LOSSLESS
      SP_Transform( my, in_line, coeff[i], coeff[i] + my );
    }
  }
#else

#ifdef Haar
      for( j = 0; j < my; j++ ) {
        coeff[i][j] =
          sm_mt * ( in_line[j * 2] + in_line[j * 2 + 1] ) / sqrt( 2. );
        coeff[i][j + my] =
          ( in_line[j * 2] - in_line[j * 2 + 1] ) / sqrt( 2. );
      }
    }
  }
#else
      Reflection( in_line, in_line + ny - 1 );
      for( j = 0, t = in_line; j < my; j++ ) {
        coeff[i][j] = sm_mt * Filter_L( T_LowPass, t++ );
        coeff[i][j + my] = Filter_H( T_HighPass, t++ );
      }
    }
  }
#endif

#endif

#ifdef   ROLL_STRUCTURE_TWO  
  // for frequency roll off structure two, we want small subbands in this level
  int dx, dy, lenx, leny; 
  // decompose LH subband further into 4 small bands 
  dx   = 0;         dy   = pdim.y/4; 
  lenx = pdim.x/4;  leny = pdim.y/4;  
  one_step_decompose(dx, dy, lenx, leny);
  // decompose HL subband further into 4 small bands
  dx   = pdim.x/4;  dy   = 0; 
  one_step_decompose(dx, dy, lenx, leny);
  // decompose HH subband further into 4 small bands 
  dx   = pdim.x/4;  dy   = pdim.y/4; 
  one_step_decompose(dx, dy, lenx, leny);
#endif 


 // calculate and subtract mean
  float s = 0;
  for( i = 0; i < mx; i++ )
    for( j = 0; j < my; j++ )
      s += coeff[i][j];
  s /= float ( mx ) * float ( my );

//  for (shift = 0; s > 1e3; shift++) s *= 0.25;
  for( shift = 0; abs( ( int )s ) > 1e3; shift++ )
    s *= 0.25;
  mean = int ( 0.5 + s );

#ifdef LOSSLESS
  int tm = mean << ( shift + shift );
#else
  float tm = ( float )( mean * pow( 4, shift ) );
#endif
  for( i = 0; i < mx; i++ )
    for( j = 0; j < my; j++ )
      coeff[i][j] -= tm;

#ifdef MY_SCALE
  for (i = 0; i < pdim.x; i++)
    for (j = 0; j < pdim.y; j++)
      coeff[i][j] *= MY_SCALE;
#endif

  delete[]temp_line;
//  cpu_time.display(" Image transformed in");

/*
  sum = 0;
  printf("Transformed band:\n");
  for( i = 0; i < dim.x; i++ ){
    for( j = 0; j < dim.y; j++ ){
//	  printf("%f\t",coeff[i][j]*coeff[i][j]);
      sum += coeff[i][j]*coeff[i][j];
	}
//	printf("\n");
  }
//  printf("\n\n");

//  printf("sum = %f\n",sum);
  sum = sum / (dim.x * dim.y);
  printf("avg = %f\n",sum);

  mse = sum;
*/
//  printf("pixel num = %d\n",(dim.x * dim.y));
}


void Image_BW::one_step_decompose(int dx, int dy, int lenx, int leny)
{
  int i = NumbTap, j = 0, k = Max(lenx, leny) + (i << 1);
  CREATE_VECTOR(temp_line, k, Pel_Type, M_MSG);
  Pel_Type * t, * in_line = temp_line + i, * out_line = in_line + j;

  int nx, ny, mx = lenx, my = leny;
 
  nx = mx;  mx >>= 1;  ny = my;  my >>= 1;
  // transformation of columns
  for (j = 0; j < ny; j++) {
     for (i = 0; i < nx; i++) in_line[i] = coeff[dx+i][dy+j];
     Reflection(in_line, in_line + nx - 1);
     for (i = 0, t = in_line; i < mx; i++) {
       coeff[dx+i][dy+j]    = Filter_L(T_LowPass, t++);
       coeff[dx+i+mx][dy+j] = Filter_H(T_HighPass, t++); 
	 } 
  }
  // transformation of rows
  for (i = 0; i < nx; i++) {
	 for (j=0; j < ny; j++ ) in_line[j] = coeff[dx+i][dy+j]; 
     Reflection(in_line, in_line + ny - 1);
     for (j = 0, t = in_line; j < my; j++) {
       coeff[dx+i][dy+j]    = Filter_L(T_LowPass,  t++);
       coeff[dx+i][dy+j+my] = Filter_H(T_HighPass, t++); 
	 } 
   }

  delete [] temp_line;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::scale_down( void )
{
  int i, j;
  int mx = pdim.x, my = pdim.y;

  mx >>= 1;
  my >>= 1;

  for( i = mx; i < pdim.x; i++ )
    for( j = 0; j < my; j++ )
      coeff[i][j] /= 2;

  for( i = mx; i < pdim.x; i++ )
    for( j = my; j < pdim.y; j++ )
      coeff[i][j] /= 2;

  for( i = 0; i < mx; i++ )
    for( j = my; j < pdim.y; j++ )
      coeff[i][j] /= 2;


}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::scale_up( void )
{

  int i, j;
  int mx = pdim.x, my = pdim.y;

  mx >>= 1;
  my >>= 1;

  for( i = mx; i < pdim.x; i++ )
    for( j = 0; j < my; j++ )
      coeff[i][j] *= 2;

  for( i = mx; i < pdim.x; i++ )
    for( j = my; j < pdim.y; j++ )
      coeff[i][j] *= 2;

  for( i = 0; i < mx; i++ )
    for( j = my; j < pdim.y; j++ )
      coeff[i][j] *= 2;


}




// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Image_BW::recover( void )
{
  if( levels <= 0 )
    Error( "cannot recover < Image_BW >" );

  // Chronometer cpu_time;
  // cpu_time.start("\n  Starting inverse transformation...");
  // printf("subband decomposition levels = %d(image_bw.cpp)\n", levels);

#ifdef LOSSLESS
  int i = 0, j = Max( pdim.x, pdim.y ), k = j << 1;
#else
  int i = NumbTap, j = 0, k = Max( pdim.x, pdim.y ) + ( i << 1 );
  //float sm_mt;
#endif
  CREATE_VECTOR( temp_line, k, Pel_Type, M_MSG );
  Pel_Type *t, *in_line = temp_line + i; // *out_line = in_line + j;

  if(high_band == NO){
	levels = Min( max_levels( pdim.x ), max_levels( pdim.y ) ); //levels = 5 for luma, levels = 4 for chroma bands
//	printf("LOW band!\n");
  }else{
	assert(high_band == YES);
	levels = 2;
//	printf("HIGH band!\n");
  }

  int lv, mx, my, nx = pdim.x >> levels, ny = pdim.y >> levels;

#ifdef MY_SCALE
  for (i = 0; i < pdim.x; i++)
    for (j = 0; j < pdim.y; j++)
      coeff[i][j] /= MY_SCALE;
#endif

// add mean

#ifdef FREQUENCY_ROLL_OFF
	FILE *froll;
	float roll_factors[15]; 
    int   factor_num; 
    
	if ( (s_level==1 || s_level==2) && frame_res==HD)
    {
#ifdef ROLL_STRUCTURE_ONE
	  froll=fopen("roll_factors_structure_one.txt", "rt"); 
	  for (factor_num=0; factor_num<6; factor_num++)
		  fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
      if (s_level==1)
		  for (i=0; i<pdim.x/2; i++)
			  for (j=0; j<pdim.y/2; j++)
			  {
				  coeff[i         ][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[0]));
				  coeff[i+pdim.x/2][j         ] *= (float)(1/pow(2, roll_factors[1]));
				  coeff[i+pdim.x/2][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[2]));
			  }

	   if (s_level==2)
		  for (i=0; i<pdim.x/2; i++)
			  for (j=0; j<pdim.y/2; j++)
			  {
				  coeff[i         ][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[3])) ;
				  coeff[i+pdim.x/2][j         ] *= (float)(1/pow(2, roll_factors[4])); 
				  coeff[i+pdim.x/2][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[5])); 
			  }
#endif 

#ifdef ROLL_STRUCTURE_TWO
	  if (s_level==1 || s_level==2)
	  {
		  froll=fopen("roll_factors_structure_two.txt", "rt"); 
		  for (factor_num=0; factor_num<15; factor_num++)
			  fscanf(froll, "%f", &(roll_factors[factor_num]));
		  fclose(froll); 
	  }

	 if (s_level==1)
	   for (i=0; i<pdim.x/4; i++)
		  for (j=0; j<pdim.y/4; j++)
		  {
			  coeff[i                  ][j+pdim.y/2]               *= (float)(1/pow(2, roll_factors[0])); 
			  coeff[i                  ][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[1])); 
			  coeff[i+pdim.x/4         ][j+pdim.y/2          ]     *= (float)(1/pow(2, roll_factors[2])); 
			  coeff[i+pdim.x/4         ][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[3])); 

			  coeff[i+pdim.x/2         ][j                   ]     *= (float)(1/pow(2, roll_factors[4])); 
			  coeff[i+pdim.x/2         ][j+pdim.y/4          ]     *= (float)(1/pow(2, roll_factors[5])); 
			  coeff[i+pdim.x/2+pdim.x/4][j                   ]     *= (float)(1/pow(2, roll_factors[6])); 
			  coeff[i+pdim.x/2+pdim.x/4][j+pdim.y/4          ]     *= (float)(1/pow(2, roll_factors[7])); 

			  coeff[i+pdim.x/2         ][j+pdim.y/2          ]     *= (float)(1/pow(2, roll_factors[8])); 
			  coeff[i+pdim.x/2         ][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[9])); 
			  coeff[i+pdim.x/2+pdim.x/4][j+pdim.y/2          ]     *= (float)(1/pow(2, roll_factors[10])); 
			  coeff[i+pdim.x/2+pdim.x/4][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[11])); 
		  }

	   if (s_level==2)
		  for (i=0; i<pdim.x/2; i++)
			  for (j=0; j<pdim.y/2; j++)
			  {
				  coeff[i         ][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[12])) ;
				  coeff[i+pdim.x/2][j         ] *= (float)(1/pow(2, roll_factors[13])); 
				  coeff[i+pdim.x/2][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[14])); 
			  }

	   int dx, dy, lenx, leny; 
	   if (s_level<2)
	   {
		   dx   = pdim.x/(1<<(2-s_level));  dy = 0; 
		   lenx = pdim.x/(1<<(2-s_level));  leny = pdim.y/(1<<(2-s_level)); 
	       one_step_recover(dx, dy, lenx, leny);
		   dx   = 0;         dy = pdim.y/(1<<(2-s_level)); 
	       one_step_recover(dx, dy, lenx, leny);
	       dx   = pdim.x/(1<<(2-s_level));  dy = pdim.y/(1<<(2-s_level)); 
	       one_step_recover(dx, dy, lenx, leny);
	   }


#endif 

   }

   if (  s_level==1 && (frame_res == SD || frame_res == SD2) )
   {
#ifdef ROLL_STRUCTURE_ONE
	  froll=fopen("roll_factors_structure_one.txt", "rt"); 
	  for (factor_num=0; factor_num<6; factor_num++)
		  fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
      if (s_level==1)
		  for (i=0; i<pdim.x/2; i++)
			  for (j=0; j<pdim.y/2; j++)
			  {
				  coeff[i         ][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[0]));
				  coeff[i+pdim.x/2][j         ] *= (float)(1/pow(2, roll_factors[1]));
				  coeff[i+pdim.x/2][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[2]));
			  }

	   if (s_level==2)
		  for (i=0; i<pdim.x/2; i++)
			  for (j=0; j<pdim.y/2; j++)
			  {
				  coeff[i         ][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[3])) ;
				  coeff[i+pdim.x/2][j         ] *= (float)(1/pow(2, roll_factors[4])); 
				  coeff[i+pdim.x/2][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[5])); 
			  }
#endif 

#ifdef ROLL_STRUCTURE_TWO
	  if (s_level==1 || s_level==2)
	  {
		  froll=fopen("roll_factors_structure_two.txt", "rt"); 
		  for (factor_num=0; factor_num<15; factor_num++)
			  fscanf(froll, "%f", &(roll_factors[factor_num]));
		  fclose(froll); 
	  }

	 if (s_level==1)
	   for (i=0; i<pdim.x/4; i++)
		  for (j=0; j<pdim.y/4; j++)
		  {
			  coeff[i                  ][j+pdim.y/2]               *= (float)(1/pow(2, roll_factors[0])); 
			  coeff[i                  ][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[1])); 
			  coeff[i+pdim.x/4         ][j+pdim.y/2          ]     *= (float)(1/pow(2, roll_factors[2])); 
			  coeff[i+pdim.x/4         ][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[3])); 

			  coeff[i+pdim.x/2         ][j                   ]     *= (float)(1/pow(2, roll_factors[4])); 
			  coeff[i+pdim.x/2         ][j+pdim.y/4          ]     *= (float)(1/pow(2, roll_factors[5])); 
			  coeff[i+pdim.x/2+pdim.x/4][j                   ]     *= (float)(1/pow(2, roll_factors[6])); 
			  coeff[i+pdim.x/2+pdim.x/4][j+pdim.y/4          ]     *= (float)(1/pow(2, roll_factors[7])); 

			  coeff[i+pdim.x/2         ][j+pdim.y/2          ]     *= (float)(1/pow(2, roll_factors[8])); 
			  coeff[i+pdim.x/2         ][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[9])); 
			  coeff[i+pdim.x/2+pdim.x/4][j+pdim.y/2          ]     *= (float)(1/pow(2, roll_factors[10])); 
			  coeff[i+pdim.x/2+pdim.x/4][j+pdim.y/4+pdim.y/2 ]     *= (float)(1/pow(2, roll_factors[11])); 
		  }

	   if (s_level==2)
		  for (i=0; i<pdim.x/2; i++)
			  for (j=0; j<pdim.y/2; j++)
			  {
				  coeff[i         ][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[12])) ;
				  coeff[i+pdim.x/2][j         ] *= (float)(1/pow(2, roll_factors[13])); 
				  coeff[i+pdim.x/2][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[14])); 
			  }

	   int dx, dy, lenx, leny; 
	   if (s_level<2)
	   {
		   dx   = pdim.x/(1<<(2-s_level));  dy = 0; 
		   lenx = pdim.x/(1<<(2-s_level));  leny = pdim.y/(1<<(2-s_level)); 
	       one_step_recover(dx, dy, lenx, leny);
		   dx   = 0;         dy = pdim.y/(1<<(2-s_level)); 
	       one_step_recover(dx, dy, lenx, leny);
	       dx   = pdim.x/(1<<(2-s_level));  dy = pdim.y/(1<<(2-s_level)); 
	       one_step_recover(dx, dy, lenx, leny);
	   }


#endif 

   }

   if ( s_level==1 && frame_res == CIF ) 
   {
	  froll=fopen("roll_factors_structure_one.txt", "rt"); 
	  for (factor_num=0; factor_num<3; factor_num++)
		  fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
		for (i=0; i<pdim.x/2; i++)
		  for (j=0; j<pdim.y/2; j++)
		  {
			  coeff[i         ][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[0]));
			  coeff[i+pdim.x/2][j         ] *= (float)(1/pow(2, roll_factors[1]));
			  coeff[i+pdim.x/2][j+pdim.y/2] *= (float)(1/pow(2, roll_factors[2]));
		  }
	}

#endif


#ifdef LOSSLESS
  int tm = mean << ( shift + shift );
#else
  float tm = ( float )( mean * pow( 4, shift ) );
#endif
  for( i = 0; i < nx; i++ )
    for( j = 0; j < ny; j++ )
      coeff[i][j] += tm;

// inverse hierarchical wavelet or S+P transformation

  for( lv = levels - 1; lv >= 0; lv-- ) {

    // shifts are doubled, multiplier is updated

    mx = nx;
    nx <<= 1;
    my = ny;
    ny <<= 1;

#ifndef LOSSLESS
    float sm_mt = 1 / ( 1 + smoothing * SmoothingFactor / ( 2 + lv * lv ) );
#endif

    // inverse transformation of rows

    for( i = 0; i < nx; i++ ) {
#ifdef LOSSLESS
      memcpy( in_line, coeff[i], ny * sizeof( Pel_Type ) );
      SP_Recover( my, in_line, in_line + my, coeff[i] );
    }
#else
      for( j = 0, t = in_line; j < my; j++ ) {
        *( t++ ) = sm_mt * coeff[i][j];
        *( t++ ) = coeff[i][j + my];
      }

#ifdef Haar
      for( j = 0; j < ny; ) {
        coeff[i][j] = ( in_line[j] + in_line[j + 1] ) / sqrt( 2. );
        coeff[i][j + 1] = ( in_line[j] - in_line[j + 1] ) / sqrt( 2. );
        j = j + 2;
      }
    }
#else
      Reflection( in_line, in_line + ny - 1 );
      for( j = 0, t = in_line; j < ny; ) {
        coeff[i][j++] = Filter_H( R_HighPass, t++ );
        coeff[i][j++] = Filter_L( R_LowPass, t++ );
      }
    }
#endif

#endif

    // inverse transformation of columns

    for( j = 0; j < ny; j++ ) {
#ifdef LOSSLESS
      for( i = 0; i < nx; i++ )
        in_line[i] = coeff[i][j];
      SP_Recover( mx, in_line, in_line + mx, out_line );
      for( i = 0; i < nx; i++ )
        coeff[i][j] = out_line[i];
    }
  }
#else
      for( i = 0, t = in_line; i < mx; i++ ) {
        *( t++ ) = sm_mt * coeff[i][j];
        *( t++ ) = coeff[i + mx][j];
      }

#ifdef Haar
      for( i = 0; i < nx; ) {
        coeff[i][j] = ( in_line[i] + in_line[i + 1] ) / sqrt( 2. );
        coeff[i + 1][j] = ( in_line[i] - in_line[i + 1] ) / sqrt( 2. );
        i = i + 2;
      }
    }
  }
#else
      Reflection( in_line, in_line + nx - 1 );
      for( i = 0, t = in_line; i < nx; ) {
        coeff[i++][j] = Filter_H( R_HighPass, t++ );
        coeff[i++][j] = Filter_L( R_LowPass, t++ );
      }
    }
  }
#endif

#endif

  levels = 0;
  delete[]temp_line;
//  cpu_time.display(" Image transformed in");
}

void Image_BW::one_step_recover(int dx, int dy, int lenx, int leny)
{

  int i = NumbTap, j = 0, k = Max(lenx, leny) + (i << 1);

  CREATE_VECTOR(temp_line, k, Pel_Type, M_MSG);
  Pel_Type * t, * in_line = temp_line + i, * out_line = in_line + j;

  int mx, my, nx = lenx/2, ny = leny/2;

  mx = nx;  nx <<= 1;  my = ny;  ny <<= 1;

  // inverse transformation of rows
  for (i = 0; i < nx; i++) {
	  for (j = 0, t = in_line; j < my; j++) {
		   *(t++) = coeff[dx+i][dy+j];  
		   *(t++) = coeff[dx+i][dy+j+my]; 
	   }
       Reflection(in_line, in_line + ny - 1);
       for (j = 0, t = in_line; j < ny;){
         coeff[dx+i][dy+j] = Filter_H(R_HighPass, t++);
		 j++;
         coeff[dx+i][dy+j] = Filter_L(R_LowPass, t++); 
		 j++;
	   } 
  }

   // inverse transformation of columns
   for (j = 0; j < ny; j++) {
      for (i = 0, t = in_line; i < mx; i++) {
        *(t++) = coeff[dx+i][dy+j];  
		*(t++) = coeff[dx+i+mx][dy+j]; 
	  }
      Reflection(in_line, in_line + nx - 1);
      for (i = 0, t = in_line; i < nx;) {
        coeff[dx+i][dy+j] = Filter_H(R_HighPass, t++);
		i++;
        coeff[dx+i][dy+j] = Filter_L(R_LowPass, t++); 
		i++;
	  } 
   } 

  delete [] temp_line;

}


void Image_BW::write_raw( char *file_name )
{
  int i;

  FILE *out_file = fopen( file_name, "wb" );
  if( out_file == NULL ) {
    char msg[512];
    sprintf(msg, "%s: %s\n", W_MSG, file_name);
    Error(msg);
  }

  for( i = 0; i < dim.x; i++ ) {
    if(fwrite( coeff[i], sizeof( Pel_Type ), dim.y, out_file ) - dim.y != 0) {
      char msg[512];
      sprintf(msg, "%s: %s\n", W_MSG, file_name);
      Error(msg);
    }
  }

  fclose( out_file );
}

void Image_BW::read_raw( Image_Coord d, char *file_name )
{
  int i;

  assign_mem( d, 1 );
  mean = levels = 0;
  bytes = 1;

  FILE *in_file = fopen( file_name, "rb" );
  if( in_file == NULL ) {
    char msg[512];
    sprintf(msg, "%s: %s\n", R_MSG, file_name);
    Error(msg);
  }

  for( i = 0; i < dim.x; i++ ) {
    if(fread( coeff[i], sizeof( Pel_Type ), dim.y, in_file ) - dim.y != 0) {
      char msg[512];
      sprintf(msg, "%s: %s\n", R_MSG, file_name);
      Error(msg);
    }
  }

  fclose( in_file );
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// end of file  < Image_BW.C >


