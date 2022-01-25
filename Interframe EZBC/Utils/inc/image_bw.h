
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//                     I M A G E   C L A S S
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//           > > > >    C++ version 11.04 -  02/01/96   < < < <
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
// - - External definitions  - - - - - - - - - - - - - - - - - - - - - -
#ifndef _IMAGE_BW_
#define _IMAGE_BW_

#include <iostream>
#ifdef LOSSLESS

typedef int Pel_Type;
#define ABS abs

#else

typedef float Pel_Type;
#define ABS fabs

#endif


struct Image_Coord
{
  int x, y;
    Image_Coord( void )
  {
    x = y = 0;
  }
  Image_Coord( int i, int j )
  {
    x = i;
    y = j;
  }
  Image_Coord( int i )
  {
    x = y = i;
  }
  Image_Coord & operator -= ( Image_Coord in_ord ) {
    x -= in_ord.x;
    y -= in_ord.y;
    return *this;
  }
  Image_Coord & operator += ( Image_Coord in_ord ) {
    x += in_ord.x;
    y += in_ord.y;
    return *this;
  }
  Image_Coord operator + ( Image_Coord in_ord )
  {
    in_ord.x += x;
    in_ord.y += y;
    return in_ord;
  }
  Image_Coord operator - ( Image_Coord in_ord )
  {
    in_ord.x = x - in_ord.x;
    in_ord.y = y - in_ord.y;
    return in_ord;
  }
  Image_Coord t(  )
  {
    return Image_Coord( y, x );
  }
};

//ostream & operator << ( ostream & out, Image_Coord & in_ord );

// - - Class definition  - - - - - - - - - - - - - - - - - - - - - - - -

enum FRAME_RES
{ SD, SD2, CIF, HD, INVALID_RES };


class Image_BW
{
  // . private data .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  Image_Coord dim, pdim;  // 表示图像的大小 coordinates坐标

  float mse; //Added on 11.26.2017

  int levels, bytes, mean, shift, smoothing; // level：空域进行小波变换的次数

  Pel_Type **coeff;

  // . private functions  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  int max_levels( int );

  void assign_mem( Image_Coord, int );

  void free_mem( void );

  void extend( void );

  // . constructor and destructor  .  .  .  .  .  .  .  .  .  .  .  .  .

public:

  int s_level;
  enum FRAME_RES frame_res; 
  int high_band; //Added on 08.14.2018
  Image_BW( void )
  {
    levels = -1;
  }
//    Image_BW(const Image_BW&);
   ~Image_BW( void )
  {
    free_mem(  );
  }

  // . public functions   .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  Pel_Type & operator[]( const Image_Coord & c ) {
    return coeff[c.x][c.y];
  }

  Pel_Type & operator(  )( int i, int j ) {
    return coeff[i][j];
  }

  Pel_Type *address( const Image_Coord & c )
  {
    return coeff[c.x] + c.y;
  }

  Image_Coord dimension( void )
  {
    return dim;
  }

  Image_Coord pyramid_dim( void )
  {
    return pdim;
  }

  int transform_mean( void )
  {
    return mean;
  }

  int mean_shift( void )
  {
    return shift;
  }

  int smoothing_factor( void )
  {
    return smoothing;
  }

  int pyramid_levels( void )
  {
    return levels;
  }                             // accurately this subband decomposition level

  int pixel_bytes( void )
  {
    return bytes;
  }

  float compare( char *file_name );

  void read_pic( Image_Coord, char *file_name, int nbytes = 1 );

  void read_float_image( Image_Coord, float *image );

  void write_pic( char *file_name );

  void write_float_image( float *image );

  void reset( Image_Coord );

  void reset( int m, int mshift, int smoothing_factor = 0 );

  void reset( Image_Coord, int nbytes, int m, int mshift,
              int smoothing_factor = 0 );

  void transform( int smoothing_factor = 0 );

  void scale_down( void );

  void scale_up( void );

  void subtract_mean(  );

  void intensity_density(  );

  void density_intensity(  );

  void add_mean(  );

  void recover( void );

  void one_step_decompose(int dx, int dy, int lenx, int leny);
  void one_step_recover(int dx, int dy, int lenx, int leny);

  void dispose( void )
  {
    free_mem(  );
  }

  float get_mse( ) //Added on 11.26.2017
  {
     return mse;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  void write_raw( char *filename );

  void read_raw( Image_Coord d, char *filename );

};                              // end definition of class  < Image_BW >

#endif
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// end of file  < Image_BW.H >
