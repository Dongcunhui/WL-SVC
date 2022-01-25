/* ========================================================================= */
/* Description: menber functions for class DecSubbandTree                    */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <math.h>
#include <assert.h>
#include "dwt_bitplane_dec.h"

DecSubbandTree::DecSubbandTree( SUBBAND_TREE_TYPE * subs, DECODER_TYPE * dec )
:SubbandTreeCodec( subs ), decoder( dec )
{

  int nbands = subband_tree->get_nband(  );
  NEW_VECTOR( dec_subs, nbands, DEC_SUBBAND_TYPE, "SubbandDec" );
  for( int k = 0; k < nbands; k++ )
    dec_subs[k].initialize( subband_tree->get_band_ptr( k ), this, dec );
}

void
DecSubbandTree::initialize( SUBBAND_TREE_TYPE * subs, DECODER_TYPE * dec )
{
  SubbandTreeCodec::initialize( subs );
  decoder = dec;

  int nbands = subband_tree->get_nband(  );
  NEW_VECTOR( dec_subs, nbands, DEC_SUBBAND_TYPE, "SubbandDec" );
  for( int k = 0; k < nbands; k++ )
    dec_subs[k].initialize( subband_tree->get_band_ptr( k ), this, dec );

}

void
DecSubbandTree::

initialize( SUBBAND_TREE_TYPE * subs, int *subband_ACcoder,
            DECODER_TYPE * dec )
{
//  SubbandTreeCodec::initialize(subs);

  int nbands = subband_tree->get_nband(  );
  NEW_VECTOR( dec_subs, nbands, DEC_SUBBAND_TYPE, "SubbandDec" );


  for( int k = 0; k < nbands; k++ )
    dec_subs[k].initialize( subband_tree->get_band_ptr( k ), this,
                            &dec[subband_ACcoder[k]] );

}

void
DecSubbandTree::reset_tree_dec( SUBBAND_TREE_TYPE * subs )
{
  reset_tree_codec( subs );
  int nbands = subband_tree->get_nband(  );
  for( int k = 0; k < nbands; k++ ) {
    dec_subs[k].reset_band_dec( subband_tree->get_band_ptr( k ) );
  }
}



//-----------------------------------------------------------------------------

//      rec_subbands(void)

//-----------------------------------------------------------------------------
void
DecSubbandTree::rec_subbands( void )
{
  int pyr_lev, nband;

  pyr_lev = subband_tree->get_pyr_levels(  );
  nband = subband_tree->get_nband(  );

  for( int k = 0; k < nband; dec_subs[k++].rec_subband(  ) );
}
