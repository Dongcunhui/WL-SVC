/* ========================================================================= */
/* Description: wavelet decomposotion structure                              */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#ifndef _SUBBAND_
#define _SUBBAND_

#include "ifc.h"
#include "general.h"
#include "image_bw.h"

typedef ifc_int SUB_COEFF_TYPE;

class SubbandLayer;             //forward def.
class SubSet;
class SubSetLayer;

typedef SubbandLayer SUBBAND_TYPE;

const SUB_COEFF_TYPE MAX_SUB_COEFF = MAX_IFC_INT;
const SUB_COEFF_TYPE MIN_SUB_COEFF = MIN_IFC_INT;
const SUB_COEFF_TYPE SIGN_BIT = MIN_IFC_INT;
const SUB_COEFF_TYPE MAG_MASK = ~SIGN_BIT;
const std_short PARITY_MASK_SHT = ~( ( std_short ) 1 );

// old definition
const SUB_COEFF_TYPE sign = 1 << ( 8 * sizeof( SUB_COEFF_TYPE ) - 1 );
const SUB_COEFF_TYPE SIGN = 1 << ( 8 * sizeof( SUB_COEFF_TYPE ) - 1 );

#define EXTRA_BIT 1
#define NO_CHILD -1

enum ORIENTATION
{ LL, HL, LH, HH };

class Subband
{
protected:
  int band_idx, transpose_flag;
  Image_Coord dim;
  Image_Coord mem_org;
  Image_Coord mem_dim;
  SUB_COEFF_TYPE **coeff;
  SubSet *wp;

  void alloc_coeff( void );
  void free_coeff( void );
  /*static const SUB_COEFF_TYPE SIGN = 1 << (8 * sizeof(SUB_COEFF_TYPE) - 1); */

public:
    Subband( void )
  {
    coeff = NULL;
  }
  Subband( Image_Coord dimen );
  Subband( Image_Coord dimen, Image_Coord memo_org, Image_Coord memo_dim );
  virtual ~ Subband( void )
  {
    free_coeff(  );
  }
  int initialization( Image_Coord dimen, Image_Coord memo_org,
                      Image_Coord memo_dim );
  int initialization( Image_Coord dimen )
  {
    return initialization( dimen, 0, dimen );
  }
  void transpose( void );
  void remove( void )
  {
    free_coeff(  );
  }
  void reset_coeff(  );
  SUB_COEFF_TYPE *get_coeff_mem_head_ptr(  )
  {
    return &( coeff[mem_org.x][mem_org.x] );
  }
  SUB_COEFF_TYPE *get_coeff_org_ptr(  )
  {
    return &( coeff[0][0] );
  }
  SUB_COEFF_TYPE **get_coeff(  )
  {
    return coeff;
  }
  SUB_COEFF_TYPE & operator[]( const Image_Coord & c ) {
    return coeff[c.x][c.y];
  }
  SUB_COEFF_TYPE & operator(  )( int i, int j ) {
    return coeff[i][j];
  }
  SUB_COEFF_TYPE *address( const Image_Coord & c )
  {
    return coeff[c.x] + c.y;
  }
  Image_Coord get_dim( void )
  {
    return dim;
  }
  Image_Coord get_mem_dim( void )
  {
    return mem_dim;
  }
  Image_Coord get_mem_org( void )
  {
    return mem_org;
  }
  int get_transpose_flag( void )
  {
    return transpose_flag;
  }
  void set_wp( SubSet * subset )
  {
    wp = subset;
  }
  SubSet *get_wp( void )
  {
    return wp;
  }
  void set_band_idx( int idx )
  {
    band_idx = idx;
  }
  int get_band_idx(  )
  {
    return band_idx;
  }
  int get_band_level(  );
  int get_child_band_idx(  );
  int get_grand_child_band_idx(  );
  int get_band_size( void )
  {
    return dim.x * dim.y;
  }
};

class SubbandLayer:public Subband  // 子带层
{
protected:
  Image_Coord org;
  int par_ind, max_msb;

  //static variable
  static int lsb;
public:

  SubbandLayer( void ):Subband(  )
  {
  }
  SubbandLayer( Image_BW & image_pyr, Image_Coord org, Image_Coord dim,
                int par );
  SubbandLayer( Image_Coord org, Image_Coord dim, int par, int msb );
  ~SubbandLayer( void )
  {
  }
  void remove( void )
  {
    free_coeff(  );
  }
  Image_Coord get_org( void )
  {
    return org;
  }
  int get_par( void )
  {
    return par_ind;
  }
  int get_max_msb( void )
  {
    return max_msb;
  }
  void set_max_msb( int msb_digit )
  {
    max_msb = msb_digit;
  }
  double subband_energy(  );
  int get_lsb( void )
  {
    return lsb;
  }
  int get_sign( int i, int j )
  {
    return ( coeff[i][j] & SIGN ) != 0;
  }
  void set_sign( int bit, int i, int j )
  {
    if( bit )
      coeff[i][j] |= SIGN;
  }
  SUB_COEFF_TYPE get_mag( int i, int j )
  {
    return coeff[i][j] & MAG_MASK;
  }
  int read_bit( SUB_COEFF_TYPE mask, int i, int j )
  {
    return coeff[i][j] & mask;
  }
  void write_bit( int bit, SUB_COEFF_TYPE mask, int i, int j )
  {
    if( bit )
      coeff[i][j] |= mask;
  }
  enum ORIENTATION get_orientation(  );
};


class SubSet
{
protected:
  int pyr_levels; // 空域分解的次数
  int nband; // 空域分解得到的子带的数量 nband = 1 + 3 * pyr_levels
  Image_Coord pdim;
  SUBBAND_TYPE **bands;

public:
    SubSet( void )
  {
    bands = NULL;
    nband = 0;
  }
   ~SubSet( void )
  {
    if( bands ) {
      for( int i = 0; i < nband; i++ )
        delete bands[i];
      delete[]bands;
    }
  }
  SUBBAND_TYPE *operator[] ( int k )
  {
    return bands[k];
  }
  SUBBAND_TYPE *get_band_ptr( int k )
  {
    return bands[k];
  }
  int get_pyr_levels( void )
  {
    return pyr_levels;
  }
  int get_nband( void )
  {
    return nband;
  }
  Image_Coord get_pdim( void )
  {
    return pdim;
  }
};


class SubSetLayer:public SubSet
{
protected:
  int max_msb;

public:
    SubSetLayer( void ):SubSet(  )
  {
    max_msb = -1;
  }
  //msb_layers = NULL: initialize for encoding
  //OW  nitialize for decoding
  SubSetLayer( Image_BW & image_pyr, int s_level, int *msb_layers = NULL );
  void SubSetLayer2Image_BW( Image_BW & image_pyr );
  int get_max_msb( void )
  {
    return max_msb;
  }
  void reset_max_msb( void );
  int cal_band_level( int band_idx );
};
#endif
