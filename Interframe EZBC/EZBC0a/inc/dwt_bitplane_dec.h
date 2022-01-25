/* ========================================================================= */
/* Description: definition for class DecSubband and DecSubbandTree           */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */
#ifndef _DWT_DEC_H
#define _DWT_DEC_H

#include "dwt_bitplane_codec.h"

typedef Decoder DECODER_TYPE;

class DecSubband;               //forward def.
class DecSubbandTree;
class EzbcDec3d;

typedef DecSubband DEC_SUBBAND_TYPE;
typedef DecSubbandTree DEC_SUBBAND_TREE_TYPE;


class DecSubband:public SubbandCodec
{
  friend class DecSubbandTree;
  friend class EzbcDec3d;
protected:
    DecSubbandTree * dec_band_tree;
  DECODER_TYPE *sub_decoder;
  std_int *LSP_break_pt;

  void create_max_and_cxt_qtrees( void );
  void setup_cxt_qtrees( void );
  void set_functions( void );
  void reset_band_dec( SUBBAND_TYPE * subband );

  void decode_LIP_cxt_AC( void );       //dwt_bitplane_dec_cxt_AC.C
  void decode_sig_leaf_cxt_AC( std_int cur_coord );
  void decode_sig_node_cxt_AC( std_int cur_coord, int lev );
  void decode_LIS_leaves_cxt_AC( void );
  void decode_cur_qtree_level_cxt_AC( void );
  void decode_LSP_cxt_AC( void );
  void decode_LSP_cxt_AC_and_bit_idx( void );

  //dwt_bitplane_enc_pos_dep_cxt_AC.C
  void decode_sig_leaf_pos_dep_cxt_AC( std_int cur_coord );
  void decode_sig_node_pos_dep_cxt_AC( std_int cur_coord, int lev );

  void ( DecSubband::*decode_LIP ) ( void );
  void ( DecSubband::*decode_LSP ) ( void );
  void ( DecSubband::*decode_sig_node ) ( std_int cur_coord, int lev );
  void ( DecSubband::*decode_sig_leaf ) ( std_int cur_coord );
  void ( DecSubband::*decode_LIS_leaves ) ( void );
  void ( DecSubband::*decode_cur_qtree_level ) ( void );
  void decode_LIS_nodes( void );
  void decode_LIS_stack(  );
public:
    DecSubband( void )
  {
    sub_decoder = NULL;
    dec_band_tree = NULL;
  }
  // ~DecSubband(void){if (sub_decoder) delete [] sub_decoder; }
  DecSubband( SUBBAND_TYPE * band, DecSubbandTree * dec_tree,
              DECODER_TYPE * dec )
  {
    initialize( band, dec_tree, dec );
  }
  void initialize( SUBBAND_TYPE * subband, DecSubbandTree * dec_tree,
                   DECODER_TYPE * dec = NULL );

  void initialize2( SUBBAND_TYPE * subband, DecSubbandTree * dec_tree );


  int start_dec_subband( void );
  void rec_subband(  );
};

class DecSubbandTree:public SubbandTreeCodec
{
  friend class DecSubband;
  friend class EzbcDec3d;
protected:
    DECODER_TYPE * decoder;
  DEC_SUBBAND_TYPE *dec_subs;
public:
    DecSubbandTree( void )
  {
    dec_subs = NULL;
  }
  DecSubbandTree( SUBBAND_TREE_TYPE * subs, DECODER_TYPE * dec );
  ~DecSubbandTree( void )
  {
    if( dec_subs )
      delete[]dec_subs;
  }
  void set_decoder( DECODER_TYPE * dec )
  // int set_decoder( DECODER_TYPE * dec )
  {
    decoder = dec;
  }
  void rec_subbands( void );
  void initialize( SUBBAND_TREE_TYPE * subs, DECODER_TYPE * dec );
  void initialize( SUBBAND_TREE_TYPE * subs, int *subband_ACcoder, DECODER_TYPE * enc );        //062602
  void reset_tree_dec( SUBBAND_TREE_TYPE * subs );
};

#endif
