/* ========================================================================= */
/* Description: member functions for class EncSubbandTree                    */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */


#include <math.h>
#include <assert.h>
#include <string.h>
#include "dwt_bitplane_enc.h"

//=============================================================================
//-----------------------------------------------------------------------------

// EncSubbandTree
// class EncSubbandTree:public SubbandTreeCodec
//-----------------------------------------------------------------------------
//comment out 062602
/*EncSubbandTree::
EncSubbandTree(SUBBAND_TREE_TYPE *subs, ENCODER_TYPE *enc)
  : SubbandTreeCodec(subs), encoder(enc)
{
  tree_msb = subband_tree->get_max_msb();
  tree_lsb = subband_tree->get_band_ptr(0)->get_lsb();

  NEW_VECTOR(enc_subs, subband_tree->get_nband(), ENC_SUBBAND_TYPE,
            "SubbandEnc");
  int nbands = subband_tree->get_nband();
  for(int k = 0; k < nbands; k++){
    enc_subs[k].initialize(subband_tree->get_band_ptr(k), this, enc);
  }
}
*/
void
EncSubbandTree::initialize( SUBBAND_TREE_TYPE * subs, ENCODER_TYPE * enc )
{
  SubbandTreeCodec::initialize( subs );
  encoder = enc;

  tree_msb = subband_tree->get_max_msb(  );     // SUBBAND_TREE_TYPE *subband_tree; 
  tree_lsb = subband_tree->get_band_ptr( 0 )->get_lsb(  );

  NEW_VECTOR( enc_subs, subband_tree->get_nband(  ), ENC_SUBBAND_TYPE,
              "SubbandEnc" );
  int nbands = subband_tree->get_nband(  );

  for( int k = 0; k < nbands; k++ ) {
    enc_subs[k].initialize( subband_tree->get_band_ptr( k ), this, encoder );
  }

}
void
EncSubbandTree::

initialize( SUBBAND_TREE_TYPE * subs, int *subband_ACcoder,
            ENCODER_TYPE * enc )
{

//  SubbandTreeCodec::initialize(subs);

  tree_msb = subband_tree->get_max_msb(  );     // SUBBAND_TREE_TYPE *subband_tree; 
  tree_lsb = subband_tree->get_band_ptr( 0 )->get_lsb(  );

  int nbands = subband_tree->get_nband(  );

  NEW_VECTOR( enc_subs, nbands, ENC_SUBBAND_TYPE, "SubbandEnc" );

  for( int k = 0; k < nbands; k++ ) {
    enc_subs[k].initialize( subband_tree->get_band_ptr( k ), this,
                            &enc[subband_ACcoder[k]] );
  }

}


void
EncSubbandTree::reset_tree_enc( SUBBAND_TREE_TYPE * subs )
{
  reset_tree_codec( subs );
  tree_msb = subband_tree->get_max_msb(  );
  int nbands = subband_tree->get_nband(  );
  for( int k = 0; k < nbands; k++ ) {
    enc_subs[k].reset_band_enc( subband_tree->get_band_ptr( k ) );
  }
}
