/* ========================================================================= */
/* Description: definition for classes EncSubband and EncSubbandTree         */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#ifndef _DWT_ENC_H
#define _DWT_ENC_H

#include "dwt_bitplane_codec.h"

typedef Encoder ENCODER_TYPE;

class EncSubband;               //forward def. 
class EncSubbandTree;
class EzbcEnc3d;

typedef EncSubband ENC_SUBBAND_TYPE;
typedef EncSubbandTree ENC_SUBBAND_TREE_TYPE;

class EncSubband:public SubbandCodec    // class SubbandCodec defined in dwt_bitplane_codec.h
{
  friend class EncSubbandTree;
  friend class EzbcEnc3d;
protected:
    SUB_COEFF_TYPE mag_mask;    //typedef ifc_int SUB_COEFF_TYPE;Utils/subband.h
  EncSubbandTree *enc_band_tree;
  ENCODER_TYPE *sub_encoder;

  void create_max_and_cxt_qtrees( void );
  void setup_max_and_cxt_qtrees( void );
  void set_functions( void );

  void encode_LIP_cxt_AC( void );       //dwt_bitplane_enc_cxt_AC.C
  void encode_sig_leaf_cxt_AC( std_int cur_coord );
  void encode_sig_node_cxt_AC( std_int cur_coord, int lev );
  void encode_LIS_leaves_cxt_AC( void );
  void encode_cur_qtree_level_cxt_AC( void );
  void encode_LSP_cxt_AC( void );
  void encode_LSP_cxt_AC_and_bit_idx( void );

  //dwt_bitplane_enc_pos_dep_cxt_AC.C
  void encode_sig_leaf_pos_dep_cxt_AC( std_int cur_coord );
  void encode_sig_node_pos_dep_cxt_AC( std_int cur_coord, int lev );

  void ( EncSubband::*encode_LIP ) ( void );    // see EncSubband::set_functions
  void ( EncSubband::*encode_sig_node ) ( std_int cur_coord, int lev );
  void ( EncSubband::*encode_sig_leaf ) ( std_int cur_coord );
  void ( EncSubband::*encode_LIS_leaves ) ( void );
  void ( EncSubband::*encode_cur_qtree_level ) ( void );
  void ( EncSubband::*encode_LSP ) ( void );
  void encode_LIS_nodes( void );
  void encode_LIS_stack( void );
  void set_mag_mask( void )
  {
    mag_mask = ( SUB_COEFF_TYPE ) ( MAX_SUB_COEFF & ( ( -1 ) << bit_idx ) );
  }

public:
    EncSubband( void ):SubbandCodec(  )
  {
    sub_encoder = NULL;
    enc_band_tree = NULL;
  }
//  EncSubband(SUBBAND_TYPE* band, EncSubbandTree *enc_tree , ENCODER_TYPE *enc) //comment out 062602
//    {initialize(band, enc_tree, enc);}
//  ~EncSubband(void){ if (sub_encoder) delete [] sub_encoder;}

  void initialize( SUBBAND_TYPE * subband, EncSubbandTree * enc_tree,
                   ENCODER_TYPE * enc = NULL );

  //void initialize2(SUBBAND_TYPE* subband, EncSubbandTree *enc_tree); //062602

  int start_enc_subband( void );
  void reset_band_enc( SUBBAND_TYPE * band );
};


class EncSubbandTree:public SubbandTreeCodec
{

  friend class EncSubband;
  friend class EzbcEnc3d;
protected:
    ENCODER_TYPE * encoder;
  ENC_SUBBAND_TYPE *enc_subs;
public:
    EncSubbandTree( void ):SubbandTreeCodec(  )
  {
    enc_subs = NULL;
    encoder = NULL;
  }
  EncSubbandTree( SUBBAND_TREE_TYPE * subs, ENCODER_TYPE * enc );
  ~EncSubbandTree( void )
  {
    if( enc_subs )
      delete[]enc_subs;
  }
  void set_encoder( ENCODER_TYPE * enc )
  {
    encoder = enc;
  }
  void initialize( SUBBAND_TREE_TYPE * subs, int *subband_ACcoder, ENCODER_TYPE * enc );        //062602
  void initialize( SUBBAND_TREE_TYPE * subs, ENCODER_TYPE * enc );
  void reset_tree_enc( SUBBAND_TREE_TYPE * subs );
};

#endif
