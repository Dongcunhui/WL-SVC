/* ========================================================================= */
/* Description: menber functions for class EncSubband*/
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <math.h>
#include <assert.h>
#include <string.h>
#include "dwt_bitplane_enc.h"

//-----------------------------------------------------------------------------

//   encode_sig_leaf_pos_dep_cxt_AC()

//-----------------------------------------------------------------------------

void
EncSubband::encode_sig_leaf_pos_dep_cxt_AC( std_int cur_coord )
{
  // std_short dim_r = qtree.dims[0].r;
  // std_short dim_c = qtree.dims[0].c;
  std_short row_gap = base_coeff[1] - base_coeff[0];
  std_short cxt_row_gap = cxt_qtree.base_cxt[1] - cxt_qtree.base_cxt[0];

  std_short r = ( cur_coord >> 16 ) << 1;
  std_short c = ( cur_coord & 0xFFFF ) << 1;

  SUB_COEFF_TYPE *sp = base_coeff[r] + c;
  PEL_CXT_TYPE *cp, *cxt_sp = cxt_qtree.base_cxt[r] + c;


  MODEL_TYPE *sig_models, *cxt_models = cxt_qtree.cxt_models;
  int *jsig_offsets = cxt_qtree.jsig_offsets[0];
  std_byte *cxt_tab, **jsig_tabs = cxt_qtree.jsig_tabs[0];

  //sign related data........
  MODEL_TYPE *sign_models = cxt_qtree.cxt_models + cxt_qtree.sign_offset;
  SUB_COEFF_TYPE *sign_cxt_tab = cxt_qtree.sign_tab;
  SIGN_CXT_TYPE **base_sign_cxt = cxt_qtree.base_sign_cxt;
  SIGN_CXT_TYPE *sign_cp, *sign_cxt_sp = base_sign_cxt[r] + c;
  std_short sign_cxt_row_gap = base_sign_cxt[1] - base_sign_cxt[0];
  SUB_COEFF_TYPE sign_predict, sign_bit;

  //interband related data
#ifdef INTERBANDS
  int cd_base_cxt_row_gap, cd_cxt_node_row_gap;
  PEL_CXT_TYPE *cd_base_cp, *cd_base_cxt_sp, **cd_base_cxt;
  NODE_CXT_TYPE *cd_node_cp, *cd_cxt_node_sp, **cd_cxt_node;

  if( child_cxt_qtree ) {
    cd_base_cxt = child_cxt_qtree->base_cxt;
    cd_base_cxt_row_gap = cd_base_cxt[1] - cd_base_cxt[0];
    cd_base_cxt_sp = cd_base_cxt[r << 1] + ( c << 1 );
    cd_cxt_node = child_cxt_qtree->cxt_nodes[1];
    cd_cxt_node_row_gap = cd_cxt_node[1] - cd_cxt_node[0];
    cd_cxt_node_sp = cd_cxt_node[r] + c;
  }
#endif

  int nbits = 0;
  cur_coord <<= 1;

  //#ifdef DEBUG
  //fprintf( stderr, "LEAF, coord = %8x, (%d,%d)\n", cur_coord, r, c );
  //#endif


  sig_models = cxt_models + jsig_offsets[0];
  cxt_tab = jsig_tabs[0];

  if( *sp & mag_mask ) {

    sub_encoder->code_symbol( 1, sig_models[cxt_tab[*cxt_sp & ZC_MASK]] );
    UPDATE_CXT( cxt_sp, cxt_row_gap );

#ifdef INTERBANDS
    if( child_cxt_qtree ) {
      UPDATE_4_CD_NODES( cd_base_cxt_sp, cd_base_cxt_row_gap );
      UPDATE_1_CD_NODE( cd_cxt_node_sp );
    }
#endif
    sign_predict = sign_cxt_tab[*sign_cxt_sp & SIGN_CXT_MASK];
    sign_bit = *sp & SIGN_BIT;

    sub_encoder->code_symbol( ( sign_predict & SIGN_BIT ) ^ sign_bit,
                              sign_models[sign_predict & SIGN_CXT_MASK] );
    UPDATE_SIGN_CXT( sign_cxt_sp, sign_cxt_row_gap, sign_bit );


    *++node_list.LSP_end = cur_coord;
    nbits++;

  } else {
    sub_encoder->code_symbol( 0, sig_models[cxt_tab[*cxt_sp & ZC_MASK]] );
    *--node_list.LIP_end = cur_coord;
  }
  cp = cxt_sp + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {

#ifdef SEPARATE_JSIG_MODELS
    sig_models = cxt_models + jsig_offsets[1];
    cxt_tab = jsig_tabs[1];
#else
    if( nbits ) {
      sig_models = cxt_models + jsig_offsets[3];
      cxt_tab = jsig_tabs[3];
    } else {
      sig_models = cxt_models + jsig_offsets[1];
      cxt_tab = jsig_tabs[1];
    }
#endif

    if( sp[1] & mag_mask ) {

      sub_encoder->code_symbol( 1, sig_models[cxt_tab[*cp & ZC_MASK]] );
      UPDATE_CXT( cp, cxt_row_gap );

#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_base_cp = cd_base_cxt_sp + 2;
        UPDATE_4_CD_NODES( cd_base_cp, cd_base_cxt_row_gap );
        cd_node_cp = cd_cxt_node_sp + 1;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif
      sign_cp = sign_cxt_sp + 1;
      sign_predict = sign_cxt_tab[*sign_cp & SIGN_CXT_MASK];
      sign_bit = sp[1] & SIGN_BIT;

      sub_encoder->code_symbol( ( sign_predict & SIGN_BIT ) ^ sign_bit,
                                sign_models[sign_predict & SIGN_CXT_MASK] );
      UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

      *++node_list.LSP_end = cur_coord | 0x1;
      nbits++;

    } else {

      sub_encoder->code_symbol( 0, sig_models[cxt_tab[*cp & ZC_MASK]] );
      *--node_list.LIP_end = cur_coord | 0x1;

    }
  }

  cp = cxt_sp + cxt_row_gap;
  if( !( *cp & OUT_OF_BOUNDS ) ) {

#ifdef SEPARATE_JSIG_MODELS
    sig_models = cxt_models + jsig_offsets[2];
    cxt_tab = jsig_tabs[2];
#else
    if( nbits ) {
      sig_models = cxt_models + jsig_offsets[3];
      cxt_tab = jsig_tabs[3];
    } else {
      sig_models = cxt_models + jsig_offsets[2];
      cxt_tab = jsig_tabs[2];
    }
#endif

    if( sp[row_gap] & mag_mask ) {

      sub_encoder->code_symbol( 1, sig_models[cxt_tab[*cp & ZC_MASK]] );
      UPDATE_CXT( cp, cxt_row_gap );

#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_base_cp = cd_base_cxt_sp + ( cd_base_cxt_row_gap << 1 );
        UPDATE_4_CD_NODES( cd_base_cp, cd_base_cxt_row_gap );
        cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif
      sign_cp = sign_cxt_sp + sign_cxt_row_gap;
      sign_predict = sign_cxt_tab[*sign_cp & SIGN_CXT_MASK];
      sign_bit = sp[row_gap] & SIGN_BIT;


      sub_encoder->code_symbol( ( sign_predict & SIGN_BIT ) ^ sign_bit,
                                sign_models[sign_predict & SIGN_CXT_MASK] );
      UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

      *++node_list.LSP_end = cur_coord | 0x10000;
      nbits++;

    } else {

      sub_encoder->code_symbol( 0, sig_models[cxt_tab[*cp & ZC_MASK]] );
      *--node_list.LIP_end = cur_coord | 0x10000;

    }
  }

  cp = cxt_sp + cxt_row_gap + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    sig_models = cxt_models + jsig_offsets[3];
    cxt_tab = jsig_tabs[3];

    if( nbits ) {
      if( sp[row_gap + 1] & mag_mask ) {

        sub_encoder->code_symbol( 1, sig_models[cxt_tab[*cp & ZC_MASK]] );
        UPDATE_CXT( cp, cxt_row_gap );

#ifdef INTERBANDS
        if( child_cxt_qtree ) {
          cd_base_cp = cd_base_cxt_sp + ( cd_base_cxt_row_gap << 1 ) + 2;
          UPDATE_4_CD_NODES( cd_base_cp, cd_base_cxt_row_gap );
          cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap + 1;
          UPDATE_1_CD_NODE( cd_node_cp );
        }
#endif

        sign_cp = sign_cxt_sp + sign_cxt_row_gap + 1;
        sign_predict = sign_cxt_tab[*sign_cp & SIGN_CXT_MASK];
        sign_bit = sp[row_gap + 1] & SIGN_BIT;

        sub_encoder->code_symbol( ( sign_predict & SIGN_BIT ) ^ sign_bit,
                                  sign_models[sign_predict & SIGN_CXT_MASK] );
        UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

        *++node_list.LSP_end = cur_coord | 0x10001;

      } else {

        sub_encoder->code_symbol( 0, sig_models[cxt_tab[*cp & ZC_MASK]] );
        *--node_list.LIP_end = cur_coord | 0x10001;
      }
    } else {
      UPDATE_CXT( cp, cxt_row_gap );

#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_base_cp = cd_base_cxt_sp + ( cd_base_cxt_row_gap << 1 ) + 2;
        UPDATE_4_CD_NODES( cd_base_cp, cd_base_cxt_row_gap );
        cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap + 1;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif
      sign_cp = sign_cxt_sp + sign_cxt_row_gap + 1;
      sign_predict = sign_cxt_tab[*sign_cp & SIGN_CXT_MASK];
      sign_bit = sp[row_gap + 1] & SIGN_BIT;

      sub_encoder->code_symbol( ( sign_predict & SIGN_BIT ) ^ sign_bit,
                                sign_models[sign_predict & SIGN_CXT_MASK] );
      UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

      *++node_list.LSP_end = cur_coord | 0x10001;
    }
  }
}

//-----------------------------------------------------------------------------

//  encode_sig_node_pos_dep_cxt_AC()

//    sp[0] is always assumed to be within the image boundaries.
//-----------------------------------------------------------------------------

void
EncSubband::encode_sig_node_pos_dep_cxt_AC( std_int cur_coord, int lev )
{
  std_short row_gap = qtree.nodes[lev][1] - qtree.nodes[lev][0];
  std_short cxt_row_gap =
    cxt_qtree.cxt_nodes[lev][1] - cxt_qtree.cxt_nodes[lev][0];
  std_short r = ( cur_coord >> 16 ) << 1;
  std_short c = ( cur_coord & 0xFFFF ) << 1;
  SUB_COEFF_TYPE *sp = qtree.nodes[lev][r] + c;
  NODE_CXT_TYPE *cp, *cxt_sp = cxt_qtree.cxt_nodes[lev][r] + c;

  MODEL_TYPE *sig_models, *cxt_models = cxt_qtree.cxt_models;
  int *jsig_offsets = cxt_qtree.jsig_offsets[lev];
  std_byte *sig_cxt_tab, **jsig_tabs = cxt_qtree.jsig_tabs[lev];

  std_int coord_buf[4];

  //......................................

#ifdef INTERBANDS
  int cd_cxt_node_row_gap;
  NODE_CXT_TYPE *cd_node_cp, *cd_cxt_node_sp, **cd_cxt_node;

  if( child_cxt_qtree ) {
    cd_cxt_node = child_cxt_qtree->cxt_nodes[lev + 1];
    cd_cxt_node_row_gap = cd_cxt_node[1] - cd_cxt_node[0];
    cd_cxt_node_sp = cd_cxt_node[r] + c;
  }
#endif

  int nbits = 0;
  cur_coord <<= 1;

  sig_models = cxt_models + jsig_offsets[0];
  sig_cxt_tab = jsig_tabs[0];

  if( *sp & mag_mask ) {

    sub_encoder->code_symbol( 1, sig_models[sig_cxt_tab[*cxt_sp & ZC_MASK]] );
    coord_buf[nbits++] = cur_coord;
    UPDATE_CXT( cxt_sp, cxt_row_gap );
#ifdef INTERBANDS
    if( child_cxt_qtree ) {
      UPDATE_1_CD_NODE( cd_cxt_node_sp );
    }
#endif

  } else {

    sub_encoder->code_symbol( 0, sig_models[sig_cxt_tab[*cxt_sp & ZC_MASK]] );
    *( ++node_list.LIS_end[lev] ) = cur_coord;

  }

  cp = cxt_sp + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {

#ifdef SEPARATE_JSIG_MODELS
    sig_models = cxt_models + jsig_offsets[1];
    sig_cxt_tab = jsig_tabs[1];
#else
    if( nbits ) {
      sig_models = cxt_models + jsig_offsets[3];
      sig_cxt_tab = jsig_tabs[3];
    } else {
      sig_models = cxt_models + jsig_offsets[1];
      sig_cxt_tab = jsig_tabs[1];
    }
#endif

    if( sp[1] & mag_mask ) {

      sub_encoder->code_symbol( 1, sig_models[sig_cxt_tab[*cp & ZC_MASK]] );
      coord_buf[nbits++] = cur_coord | 0x1;
      UPDATE_CXT( cp, cxt_row_gap );
#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_node_cp = cd_cxt_node_sp + 1;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif

    } else {

      sub_encoder->code_symbol( 0, sig_models[sig_cxt_tab[*cp & ZC_MASK]] );
      *( ++node_list.LIS_end[lev] ) = cur_coord | 0x1;

    }
  }

  cp = cxt_sp + cxt_row_gap;
  if( !( *cp & OUT_OF_BOUNDS ) ) {

#ifdef SEPARATE_JSIG_MODELS
    sig_models = cxt_models + jsig_offsets[2];
    sig_cxt_tab = jsig_tabs[2];
#else
    if( nbits ) {
      sig_models = cxt_models + jsig_offsets[3];
      sig_cxt_tab = jsig_tabs[3];
    } else {
      sig_models = cxt_models + jsig_offsets[2];
      sig_cxt_tab = jsig_tabs[2];
    }
#endif

    if( sp[row_gap] & mag_mask ) {

      sub_encoder->code_symbol( 1, sig_models[sig_cxt_tab[*cp & ZC_MASK]] );
      coord_buf[nbits++] = cur_coord | 0x10000;
      UPDATE_CXT( cp, cxt_row_gap );
#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif

    } else {

      sub_encoder->code_symbol( 0, sig_models[sig_cxt_tab[*cp & ZC_MASK]] );
      *( ++node_list.LIS_end[lev] ) = cur_coord | 0x10000;

    }
  }

  cp = cxt_sp + cxt_row_gap + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    sig_models = cxt_models + jsig_offsets[3];
    sig_cxt_tab = jsig_tabs[3];

    if( nbits ) {
      if( sp[row_gap + 1] & mag_mask ) {
        sub_encoder->code_symbol( 1, sig_models[sig_cxt_tab[*cp & ZC_MASK]] );
        coord_buf[nbits++] = cur_coord | 0x10001;

        UPDATE_CXT( cp, cxt_row_gap );
#ifdef INTERBANDS
        if( child_cxt_qtree ) {
          cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap + 1;
          UPDATE_1_CD_NODE( cd_node_cp );
        }
#endif

      } else {
        sub_encoder->code_symbol( 0, sig_models[sig_cxt_tab[*cp & ZC_MASK]] );
        *( ++node_list.LIS_end[lev] ) = cur_coord | 0x10001;

      }
    } else {
      coord_buf[nbits++] = cur_coord | 0x10001;

      UPDATE_CXT( cp, cxt_row_gap );
#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap + 1;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif
    }
  }

  while( nbits ) {
    ( ++node_list.LIS_stack_top )->node = coord_buf[--nbits];
    node_list.LIS_stack_top->level = lev;
  }
}
