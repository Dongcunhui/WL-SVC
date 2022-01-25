/* ========================================================================= */
/* Description: menber functions for class DecSubband with contex modeling   */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <math.h>
#include <assert.h>
#include <string.h>
#include "dwt_bitplane_dec.h"

//-----------------------------------------------------------------------------

//  decode_LIP_cxt_AC()

//-----------------------------------------------------------------------------
void
DecSubband::decode_LIP_cxt_AC(  )
{
  std_short r, c;
  std_int cur_coord;
  PEL_CXT_TYPE *cxt_sp;
  std_int *LIP_cur, *LIP_end, *LIP_old_end;


  std_int *LSP_end = node_list.LSP_end;
  // std_int *LSP = node_list.LSP;
  std_int *LIP = node_list.LIP;
  PEL_CXT_TYPE **base_cxt = cxt_qtree.base_cxt;
  std_short row_gap = base_cxt[1] - base_cxt[0];
  MODEL_TYPE *LIP_models = cxt_qtree.cxt_models + cxt_qtree.LIP_offset;
  std_byte *cxt_tab = cxt_qtree.LIP_tab;

  //sign related data........
  MODEL_TYPE *sign_models = cxt_qtree.cxt_models + cxt_qtree.sign_offset;
  SUB_COEFF_TYPE *sign_cxt_tab = cxt_qtree.sign_tab;
  SIGN_CXT_TYPE *sign_cxt_sp, **base_sign_cxt = cxt_qtree.base_sign_cxt;
  std_short sign_cxt_row_gap = base_sign_cxt[1] - base_sign_cxt[0];
  SUB_COEFF_TYPE sign_predict, sign_bit;

  //interband related data
#ifdef INTERBANDS
  int cd_base_cxt_row_gap;
  PEL_CXT_TYPE *cd_base_cxt_sp, **cd_base_cxt;
  NODE_CXT_TYPE *cd_cxt_node_sp, **cd_cxt_node;

  if( child_cxt_qtree ) {
    cd_base_cxt = child_cxt_qtree->base_cxt;
    cd_base_cxt_row_gap = cd_base_cxt[1] - cd_base_cxt[0];
    cd_cxt_node = child_cxt_qtree->cxt_nodes[1];
  }
#endif

#ifdef INITIALIZE_LIP_MODELS_FROM_PAR
  MODEL_TYPE *par_LIP_models;

  if( ( bit_idx == max_bit_idx - 1 ) && ( par_cxt_qtree ) ) {
    assert( cxt_qtree.LIP_cxts == par_cxt_qtree->LIP_cxts );
    par_LIP_models = par_cxt_qtree->cxt_models + par_cxt_qtree->LIP_offset;
    for( int i = cxt_qtree.LIP_cxts - 1; i >= 0;
         LIP_models[i].reset( par_LIP_models[i] ),
         LIP_models[i--].taub_scale(  ) );
  }
#endif

  LIP_cur = LIP_end = LIP + node_list.LIP_sz;
  LIP_old_end = node_list.LIP_end;

  while( LIP_cur != LIP_old_end ) {
    cur_coord = *( --LIP_cur );
    r = cur_coord >> 16;
    c = cur_coord & 0xFFFF;
    cxt_sp = base_cxt[r] + c;

    if( sub_decoder->decode_symbol( LIP_models[cxt_tab[*cxt_sp & ZC_MASK]] ) ) {

#ifdef TEST_CODING_RATE         // June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif


      if( ++LSP_end == LIP_old_end )    //LSP_pos runs beyond the end of the LSP
        *( LIP_cur++ ) = *( LIP_old_end++ );
      *LSP_end = cur_coord;

      base_coeff[r][c] |= bit_idx_mask;
      UPDATE_CXT( cxt_sp, row_gap );

#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_base_cxt_sp = cd_base_cxt[r << 1] + ( c << 1 );
        UPDATE_4_CD_NODES( cd_base_cxt_sp, cd_base_cxt_row_gap );
        cd_cxt_node_sp = cd_cxt_node[r] + c;
        UPDATE_1_CD_NODE( cd_cxt_node_sp );
      }
#endif

      sign_cxt_sp = base_sign_cxt[r] + c;
      sign_predict = sign_cxt_tab[*sign_cxt_sp & SIGN_CXT_MASK];
      sign_bit =
        sub_decoder->
        decode_symbol( sign_models[sign_predict & SIGN_CXT_MASK] )
        ? ( sign_predict & SIGN_BIT ) ^ SIGN_BIT : sign_predict & SIGN_BIT;

      if( sign_bit )
        base_coeff[r][c] |= SIGN_BIT;
      UPDATE_SIGN_CXT( sign_cxt_sp, sign_cxt_row_gap, sign_bit );

    } else
      *( --LIP_end ) = cur_coord;

#ifdef TEST_CODING_RATE
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
      node_list.LSP_end = LSP_end;
      return;
    }
#endif

  }
  node_list.LSP_end = LSP_end;
  node_list.LIP_end = LIP_end;

#ifdef LSP_BIT_IDX
  node_list.LSP_ids[bit_idx][0] = LSP_end;
#endif
}

//-----------------------------------------------------------------------------

//  decode_sig_leaf_cxt_AC()

//-----------------------------------------------------------------------------
void
  DecSubband::
//bool DecSubband::   
decode_sig_leaf_cxt_AC( std_int cur_coord )
{
  // std_short dim_r = qtree.dims[0].r;
  // std_short dim_c = qtree.dims[0].c;
  std_short
    row_gap =
    base_coeff[1] -
    base_coeff[0];
  std_short
    cxt_row_gap =
    cxt_qtree.
    base_cxt[1] -
    cxt_qtree.
    base_cxt[0];
  std_short
    r = ( cur_coord >> 16 ) <<
    1;
  std_short
    c = ( cur_coord & 0xFFFF ) <<
    1;

  SUB_COEFF_TYPE *
    sp =
    base_coeff[r] +
    c;
  PEL_CXT_TYPE *
  cp, *
    cxt_sp =
    cxt_qtree.
    base_cxt[r] +
    c;
  MODEL_TYPE *
    sig_models =
    cxt_qtree.
    cxt_models +
    cxt_qtree.
    sig_offsets[0];
  std_byte *
    cxt_tab =
    cxt_qtree.
    sig_tabs[0];

  //sign related data........
  MODEL_TYPE *
    sign_models =
    cxt_qtree.
    cxt_models +
    cxt_qtree.
    sign_offset;
  SUB_COEFF_TYPE *
    sign_cxt_tab =
    cxt_qtree.
    sign_tab;
  SIGN_CXT_TYPE **
    base_sign_cxt =
    cxt_qtree.
    base_sign_cxt;
  SIGN_CXT_TYPE *
  sign_cp, *
    sign_cxt_sp =
    base_sign_cxt[r] +
    c;
  std_short
    sign_cxt_row_gap =
    base_sign_cxt[1] -
    base_sign_cxt[0];
  SUB_COEFF_TYPE
    sign_predict,
    sign_bit;

  //interband related data
#ifdef INTERBANDS
  int
    cd_base_cxt_row_gap,
    cd_cxt_node_row_gap;
  PEL_CXT_TYPE *
  cd_base_cp, *
  cd_base_cxt_sp, **
    cd_base_cxt;
  NODE_CXT_TYPE *
  cd_node_cp, *
  cd_cxt_node_sp, **
    cd_cxt_node;

  if( child_cxt_qtree ) {
    cd_base_cxt = child_cxt_qtree->base_cxt;
    cd_base_cxt_row_gap = cd_base_cxt[1] - cd_base_cxt[0];
    cd_base_cxt_sp = cd_base_cxt[r << 1] + ( c << 1 );
    cd_cxt_node = child_cxt_qtree->cxt_nodes[1];
    cd_cxt_node_row_gap = cd_cxt_node[1] - cd_cxt_node[0];
    cd_cxt_node_sp = cd_cxt_node[r] + c;
  }
#endif

  int
    nbits =
    0;
  cur_coord <<= 1;

  if( sub_decoder->decode_symbol( sig_models[cxt_tab[*cxt_sp & ZC_MASK]] ) ) {


#ifdef TEST_CODING_RATE         // June15
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif

    *sp |= bit_idx_mask;
    UPDATE_CXT( cxt_sp, cxt_row_gap );

#ifdef INTERBANDS
    if( child_cxt_qtree ) {
      UPDATE_4_CD_NODES( cd_base_cxt_sp, cd_base_cxt_row_gap );
      UPDATE_1_CD_NODE( cd_cxt_node_sp );
    }
#endif

    sign_predict = sign_cxt_tab[*sign_cxt_sp & SIGN_CXT_MASK];
    sign_bit =
      sub_decoder->
      decode_symbol( sign_models[sign_predict & SIGN_CXT_MASK] )
      ? ( sign_predict & SIGN_BIT ) ^ SIGN_BIT : sign_predict & SIGN_BIT;
    if( sign_bit )
      *sp |= SIGN_BIT;
    UPDATE_SIGN_CXT( sign_cxt_sp, sign_cxt_row_gap, sign_bit );

    *++node_list.LSP_end = cur_coord;
    nbits++;
  } else {
    *--node_list.LIP_end = cur_coord;
  }



#ifdef TEST_CODING_RATE         // June15
  if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
    return;
#endif




  cp = cxt_sp + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    if( sub_decoder->decode_symbol( sig_models[cxt_tab[*cp & ZC_MASK]] ) ) {

#ifdef TEST_CODING_RATE         // June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif

      sp[1] |= bit_idx_mask;
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
      sign_bit =
        sub_decoder->
        decode_symbol( sign_models[sign_predict & SIGN_CXT_MASK] )
        ? ( sign_predict & SIGN_BIT ) ^ SIGN_BIT : sign_predict & SIGN_BIT;
      if( sign_bit )
        sp[1] |= SIGN_BIT;
      UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

      *++node_list.LSP_end = cur_coord | 0x1;
      nbits++;
    } else {
      *--node_list.LIP_end = cur_coord | 0x1;
    }

#ifdef TEST_CODING_RATE         // June15
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif
  }




  cp = cxt_sp + cxt_row_gap;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    if( sub_decoder->decode_symbol( sig_models[cxt_tab[*cp & ZC_MASK]] ) ) {


#ifdef TEST_CODING_RATE         // June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif


      sp[row_gap] |= bit_idx_mask;
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
      sign_bit =
        sub_decoder->
        decode_symbol( sign_models[sign_predict & SIGN_CXT_MASK] )
        ? ( sign_predict & SIGN_BIT ) ^ SIGN_BIT : sign_predict & SIGN_BIT;
      if( sign_bit )
        sp[row_gap] |= SIGN_BIT;
      UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

      *++node_list.LSP_end = cur_coord | 0x10000;
      nbits++;
    } else {
      *--node_list.LIP_end = cur_coord | 0x10000;
    }
#ifdef TEST_CODING_RATE         // June15
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif
  }







  cp = cxt_sp + cxt_row_gap + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    if( nbits ) {
      if( sub_decoder->decode_symbol( sig_models[cxt_tab[*cp & ZC_MASK]] ) ) {


#ifdef TEST_CODING_RATE         // June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
          return;
#endif


        sp[row_gap + 1] |= bit_idx_mask;
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
        sign_bit =
          sub_decoder->
          decode_symbol( sign_models[sign_predict & SIGN_CXT_MASK] )
          ? ( sign_predict & SIGN_BIT ) ^ SIGN_BIT : sign_predict & SIGN_BIT;
        if( sign_bit )
          sp[row_gap + 1] |= SIGN_BIT;
        UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

        *++node_list.LSP_end = cur_coord | 0x10001;
      } else {
        *--node_list.LIP_end = cur_coord | 0x10001;

      }
    } else {
      sp[row_gap + 1] |= bit_idx_mask;
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
      sign_bit =
        sub_decoder->
        decode_symbol( sign_models[sign_predict & SIGN_CXT_MASK] )
        ? ( sign_predict & SIGN_BIT ) ^ SIGN_BIT : sign_predict & SIGN_BIT;

      if( sign_bit )
        sp[row_gap + 1] |= SIGN_BIT;
      UPDATE_SIGN_CXT( sign_cp, sign_cxt_row_gap, sign_bit );

      *++node_list.LSP_end = cur_coord | 0x10001;
    }
  }


  return;
}

//-----------------------------------------------------------------------------

//  decode_sig_node_cxt_AC()

//-----------------------------------------------------------------------------

void
DecSubband::decode_sig_node_cxt_AC( std_int cur_coord, int lev )
{
  std_short cxt_row_gap =
    cxt_qtree.cxt_nodes[lev][1] - cxt_qtree.cxt_nodes[lev][0];
  std_short r = ( cur_coord >> 16 ) << 1;
  std_short c = ( cur_coord & 0xFFFF ) << 1;
  NODE_CXT_TYPE *cp, *cxt_sp = cxt_qtree.cxt_nodes[lev][r] + c;
  MODEL_TYPE *sig_models = cxt_qtree.cxt_models + cxt_qtree.sig_offsets[lev];
  std_byte *sig_cxt_tab = cxt_qtree.sig_tabs[lev];

  //new stuffs
  std_int coord_buf[4];

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

  if( sub_decoder->
      decode_symbol( sig_models[sig_cxt_tab[*cxt_sp & ZC_MASK]] ) ) {

#ifdef TEST_CODING_RATE         //June15
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif



    coord_buf[nbits++] = cur_coord;

    UPDATE_CXT( cxt_sp, cxt_row_gap );
#ifdef INTERBANDS
    if( child_cxt_qtree ) {
      UPDATE_1_CD_NODE( cd_cxt_node_sp );
    }
#endif

    //nbits++;
  } else
    *( ++node_list.LIS_end[lev] ) = cur_coord;


#ifdef TEST_CODING_RATE         //June15
  if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
    return;
#endif


  cp = cxt_sp + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    if( sub_decoder->decode_symbol( sig_models[sig_cxt_tab[*cp & ZC_MASK]] ) ) {

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif


      coord_buf[nbits++] = cur_coord | 0x1;

      UPDATE_CXT( cp, cxt_row_gap );
#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_node_cp = cd_cxt_node_sp + 1;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif

    } else
      *( ++node_list.LIS_end[lev] ) = cur_coord | 0x1;

#ifdef TEST_CODING_RATE         //June15
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif
  }


  cp = cxt_sp + cxt_row_gap;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    if( sub_decoder->decode_symbol( sig_models[sig_cxt_tab[*cp & ZC_MASK]] ) ) {

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif


      coord_buf[nbits++] = cur_coord | 0x10000;

      UPDATE_CXT( cp, cxt_row_gap );
#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap;
        UPDATE_1_CD_NODE( cd_node_cp );
      }
#endif

    } else
      *( ++node_list.LIS_end[lev] ) = cur_coord | 0x10000;

#ifdef TEST_CODING_RATE         //June15
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif
  }


  cp = cxt_sp + cxt_row_gap + 1;
  if( !( *cp & OUT_OF_BOUNDS ) ) {
    if( nbits ) {
      if( sub_decoder->
          decode_symbol( sig_models[sig_cxt_tab[*cp & ZC_MASK]] ) ) {

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
          return;
#endif


        coord_buf[nbits++] = cur_coord | 0x10001;

        UPDATE_CXT( cp, cxt_row_gap );
#ifdef INTERBANDS
        if( child_cxt_qtree ) {
          cd_node_cp = cd_cxt_node_sp + cd_cxt_node_row_gap + 1;
          UPDATE_1_CD_NODE( cd_node_cp );
        }
#endif
      } else
        *( ++node_list.LIS_end[lev] ) = cur_coord | 0x10001;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif

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

#ifdef TEST_CODING_RATE         //June15
  if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
    return;
#endif


  while( nbits ) {
    ( ++node_list.LIS_stack_top )->node = coord_buf[--nbits];
    node_list.LIS_stack_top->level = lev;
  }

}

//-----------------------------------------------------------------------------

//  decode_LIS_leaves_cxt_AC(void)

//-----------------------------------------------------------------------------
void
DecSubband::decode_LIS_leaves_cxt_AC( void )
{
  int i;
  std_short r, c;
  std_int cur_coord, *LIS_cur, *LIS_end, *LIS_end_old;
  NODE_CXT_TYPE *cxt_sp;
  NODE_CXT_TYPE **cxt_nodes = cxt_qtree.cxt_nodes[1];
  std_short cxt_row_gap = cxt_nodes[1] - cxt_nodes[0];
  MODEL_TYPE *leaf_models = cxt_qtree.cxt_models + cxt_qtree.node_offsets[1];
  std_byte *leaf_cxt_tab = cxt_qtree.node_tabs[1];

#ifdef INTERBANDS
  NODE_CXT_TYPE *cd_cxt_node_sp, **cd_cxt_node;
  if( child_cxt_qtree ) {
    cd_cxt_node = child_cxt_qtree->cxt_nodes[2];
  }
#endif


#ifdef SCALE_SIGN_MODELS_AT_LEVEL_ENDS
  MODEL_TYPE *sign_models = cxt_qtree.cxt_models + cxt_qtree.sign_offset;
  for( i = cxt_qtree.sign_cxts - 1; i >= 0; sign_models[i--].taub_scale(  ) );
#endif


#ifdef INITIALIZE_NODE_MODELS_FROM_PAR
  if( ( bit_idx == max_bit_idx - 1 ) && ( par_cxt_qtree ) ) {
    MODEL_TYPE *par_leaf_models;

    assert( cxt_qtree.node_cxts[1] == par_cxt_qtree->node_cxts[1] );
    par_leaf_models = par_cxt_qtree->cxt_models +
      par_cxt_qtree->node_offsets[1];
    for( i = cxt_qtree.node_cxts[1] - 1; i >= 0;
         leaf_models[i].reset( par_leaf_models[i] ),
         leaf_models[i--].taub_scale(  ) );
  }
#endif

  LIS_cur = LIS_end = node_list.LIS[1] - 1;
  LIS_end_old = node_list.LIS_end[1];

  for( ; LIS_cur < LIS_end_old; ) {
    cur_coord = *( ++LIS_cur );
    r = cur_coord >> 16;
    c = cur_coord & 0xFFFF;
    cxt_sp = cxt_nodes[r] + c;

    if( sub_decoder->
        decode_symbol( leaf_models[leaf_cxt_tab[*cxt_sp & ZC_MASK]] ) ) {


#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif


      ( this->*decode_sig_leaf ) ( cur_coord );

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif




      UPDATE_CXT( cxt_sp, cxt_row_gap );
#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_cxt_node_sp = cd_cxt_node[r] + c;
        UPDATE_1_CD_NODE( cd_cxt_node_sp );
      }
#endif
    } else
      *++LIS_end = cur_coord;

#ifdef TEST_CODING_RATE
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif

  }
  node_list.LIS_end[1] = LIS_end;

#ifdef LSP_BIT_IDX
  node_list.LSP_ids[bit_idx][1] = node_list.LSP_end;
#endif
}

//-----------------------------------------------------------------------------

//  decode_cur_qtree_level_cxt_AC(void)

//-----------------------------------------------------------------------------

void
DecSubband::decode_cur_qtree_level_cxt_AC( void )
{
  assert( cur_lev < qtree.depth );

  int i;
  std_short r, c;
  std_int cur_coord, *pLIS_cur, *pLIS_end, *pLIS_end_old;       //, *cLIS_end
  NODE_CXT_TYPE *cxt_sp, **cxt_nodes = cxt_qtree.cxt_nodes[cur_lev];
  std_short cxt_row_gap = cxt_nodes[1] - cxt_nodes[0];
  MODEL_TYPE *node_models =
    cxt_qtree.cxt_models + cxt_qtree.node_offsets[cur_lev];
  std_byte *node_cxt_tab = cxt_qtree.node_tabs[cur_lev];

#ifdef INTERBANDS
  NODE_CXT_TYPE *cd_cxt_node_sp, **cd_cxt_node;
  if( child_cxt_qtree ) {
    cd_cxt_node = child_cxt_qtree->cxt_nodes[cur_lev + 1];
  }
#endif

#ifdef SCALE_SIGN_MODELS_AT_LEVEL_ENDS
  MODEL_TYPE *sign_models = cxt_qtree.cxt_models + cxt_qtree.sign_offset;
  for( i = cxt_qtree.sign_cxts - 1; i >= 0; sign_models[i--].taub_scale(  ) );
#endif


  pLIS_cur = pLIS_end = node_list.LIS[cur_lev] - 1;
  pLIS_end_old = node_list.LIS_end[cur_lev];
  while( pLIS_cur < pLIS_end_old ) {
    cur_coord = *( ++pLIS_cur );
    r = cur_coord >> 16;
    c = cur_coord & 0xFFFF;
    cxt_sp = cxt_nodes[r] + c;
    if( sub_decoder->
        decode_symbol( node_models[node_cxt_tab[*cxt_sp & ZC_MASK]] ) ) {

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif


      UPDATE_CXT( cxt_sp, cxt_row_gap );
#ifdef INTERBANDS
      if( child_cxt_qtree ) {
        cd_cxt_node_sp = cd_cxt_node[r] + c;
        UPDATE_1_CD_NODE( cd_cxt_node_sp );
      }
#endif

      ( this->*decode_sig_node ) ( cur_coord, cur_lev - 1 );

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif


      decode_LIS_stack(  );

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif
    } else
      *++pLIS_end = cur_coord;

#ifdef TEST_CODING_RATE
    if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
      return;
#endif

  }
  node_list.LIS_end[cur_lev] = pLIS_end;

#ifdef LSP_BIT_IDX
  node_list.LSP_ids[bit_idx][cur_lev] = node_list.LSP_end;
#endif
  cur_lev++;
}

//-----------------------------------------------------------------------------

//  decode_LSP_cxt_AC()

//-----------------------------------------------------------------------------
void
DecSubband::decode_LSP_cxt_AC(  )
{

  std_short r, c;
  std_int *last, *sp, LSP_offset;
  SUB_COEFF_TYPE mag;
  std_int LSP_plane = node_list.LSP_plane;
  std_int *LSP_mark = node_list.LSP_mark;
  MODEL_TYPE *LSP_models = cxt_qtree.cxt_models + cxt_qtree.LSP_offset;
#ifdef LSP_LUT
  std_byte *LSP_cxt_tab = cxt_qtree.LSP_tab;
#endif
  PEL_CXT_TYPE **base_cxt = cxt_qtree.base_cxt;

  LSP_mark[LSP_plane] = node_list.LSP_end - node_list.LSP;

  if( LSP_plane ) {
    last = node_list.LSP + LSP_mark[LSP_plane - 1];
//    for (sp = node_list.LSP; sp <= last; sp++){
    for( sp = last; sp >= node_list.LSP; sp-- ) {
      r = *sp >> 16;
      c = *sp & 0xFFFF;
      mag = ( base_coeff[r][c] & MAG_MASK ) >> bit_idx;

      LSP_offset = 0;

#ifdef LSP_LUT
      if( mag < 4 ) {
        LSP_offset = LSP_cxt_tab[base_cxt[r][c]] + 1;
      }
#else
      if( mag < 4 ) {
        if( base_cxt[r][c] & 0x0FFF )
          LSP_offset++;
      } else
        LSP_offset += 2;
#endif

      //  if(sub_decoder->decode_symbol(*LSP_models))
      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif

    }
  }
  ++( node_list.LSP_plane );
  LSP_break_pt = node_list.LSP;
}

//-----------------------------------------------------------------------------

//  decode_LSP_cxt_AC_and_bit_idx()

//-----------------------------------------------------------------------------
void
DecSubband::decode_LSP_cxt_AC_and_bit_idx(  )
{


  std_short r, c;
  std_int *last, *sp, LSP_offset;
  std_int LSP_plane = node_list.LSP_plane;
  std_int *LSP_mark = node_list.LSP_mark;

  std_int **LSP_bit_idx_marks = node_list.LSP_bit_idx_marks;
  std_int ***LSP_ids = node_list.LSP_ids;
  int i, depth = qtree.depth;

  MODEL_TYPE *LSP_models = cxt_qtree.cxt_models + cxt_qtree.LSP_offset;
#ifdef LSP_LUT
  std_byte *LSP_cxt_tab = cxt_qtree.LSP_tab;
#endif
  PEL_CXT_TYPE **base_cxt = cxt_qtree.base_cxt;

  LSP_mark[LSP_plane] = node_list.LSP_end - node_list.LSP;

  LSP_bit_idx_marks[bit_idx] = node_list.LSP_end;

//new stuffs
#ifdef INITIALIZE_LSP_MODELS_FROM_PAR
  MODEL_TYPE *par_LSP_models;

  if( ( bit_idx == max_bit_idx - 1 ) && ( par_cxt_qtree ) ) {
    assert( cxt_qtree.LSP_cxts == par_cxt_qtree->LSP_cxts );
    par_LSP_models = par_cxt_qtree->cxt_models + par_cxt_qtree->LSP_offset;
    for( i = cxt_qtree.LSP_cxts - 1; i >= 0;
         LSP_models[i].reset( par_LSP_models[i] ),
         LSP_models[i--].taub_scale(  ) );
  }
#endif

#ifdef LSP_LUT

  int LSP_lev, LSP_set_offset;
  if( bit_idx < max_bit_idx ) {
#ifdef MERGE_HIGH_IDX_SETS_OF_TOP_LSP_PLANE

    LSP_lev = depth - 1;
#ifdef MERGE_ALL_LSP_PLANES
    LSP_set_offset = 0;
#else
    LSP_set_offset = 20;
#endif

#else
    LSP_lev = ( depth < 6 ) ? ( depth - 1 ) : 5;
#endif

    for( sp = LSP_bit_idx_marks[bit_idx + 1]; LSP_lev > 1; LSP_lev-- ) {

#ifndef MERGE_HIGH_IDX_SETS_OF_TOP_LSP_PLANE
      LSP_set_offset = LSP_lev * 10;
#endif
      for( last = LSP_ids[bit_idx + 1][LSP_lev - 1]; sp > last; sp-- ) {
        //decoding unit.................................
        r = *sp >> 16;
        c = *sp & 0xFFFF;
        LSP_offset = LSP_set_offset + LSP_cxt_tab[base_cxt[r][c]];

        if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
          base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
          LSP_break_pt = sp + 1;
          return;
        }
#endif

        //end of decoding unit.................................
      }

/*      June15
#ifdef TEST_CODING_RATE
   if(sub_decoder->bytes_used() > sub_decoder->byte_budget){
     LSP_break_pt = sp + 1;
     return;
   }
#endif
*/

#ifdef SCALE_AT_LSP_BREAK_PTS
      for( i = cxt_qtree.LSP_cxts - 1; i >= 0;
           LSP_models[i--].taub_scale(  ) );
#endif
    }


    if( bit_idx < max_bit_idx - 1 ) {
      //from leaves
      last = LSP_ids[bit_idx + 1][0];

#ifndef MERGE_SETS_OF_TOP_LSP_PLANE
      LSP_set_offset = 10;
#else
      //for(i = cxt_qtree.LSP_cxts - 1; i >= 0; LSP_models[i--].taub_scale());
#endif

      for( ; sp > last; sp-- ) {
        //decoding unit.................................
        r = *sp >> 16;
        c = *sp & 0xFFFF;
        LSP_offset = LSP_set_offset + LSP_cxt_tab[base_cxt[r][c]];

        if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
          base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
          LSP_break_pt = sp + 1;
          return;
        }
#endif


        //end of decoding unit.................................
      }
/*  June15
#ifdef TEST_CODING_RATE
   if(sub_decoder->bytes_used() > sub_decoder->byte_budget){
     LSP_break_pt = sp + 1;
     return;
   }
#endif
*/

#ifdef MERGE_SETS_0_1_OF_TOP_PLANE

#ifdef SCALE_AT_LSP_BREAK_PTS
      for( i = cxt_qtree.LSP_cxts - 1; i >= 0;
           LSP_models[i--].taub_scale(  ) );
#endif

#else
      LSP_set_offset = 0;
#endif

      last = LSP_bit_idx_marks[bit_idx + 2];
      for( ; sp > last; sp-- ) {
        //decoding unit.................................
        r = *sp >> 16;
        c = *sp & 0xFFFF;
        LSP_offset = LSP_set_offset + LSP_cxt_tab[base_cxt[r][c]];

        if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
          base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
          LSP_break_pt = sp + 1;
          return;
        }
#endif


        //end of decoding unit.................................
      }

/*  June15
#ifdef TEST_CODING_RATE
   if(sub_decoder->bytes_used() > sub_decoder->byte_budget){
     LSP_break_pt = sp + 1;
     return;
   }
#endif
*/

    }
  }

  if( bit_idx < max_bit_idx - 1 ) {


#ifdef SEPARATE_SETS_OF_PLANE_1

#ifndef MERGE_LSP_PLANE_1
    LSP_set_offset = 60;
#else

#ifdef SCALE_AT_LSP_BREAK_PTS
    for( i = cxt_qtree.LSP_cxts - 1; i >= 0; LSP_models[i--].taub_scale(  ) );
#endif

#endif //MERGE_LSP_PLANE_1

    LSP_lev = depth - 1;
    for( ; LSP_lev >= 0; LSP_lev-- ) {
      if( LSP_lev )
        last = LSP_ids[bit_idx + 2][LSP_lev - 1];
      else
        last = ( bit_idx == max_bit_idx - 2 ) ?
          node_list.LSP - 1 : LSP_bit_idx_marks[bit_idx + 3];

      for( ; sp > last; sp-- ) {
        //decoding unit.................................
        r = *sp >> 16;
        c = *sp & 0xFFFF;
        LSP_offset = LSP_set_offset + LSP_cxt_tab[base_cxt[r][c]];

        if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
          base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
          LSP_break_pt = sp + 1;
          return;
        }
#endif


        //end of decoding unit.................................
      }

/*  June15
#ifdef TEST_CODING_RATE
   if(sub_decoder->bytes_used() > sub_decoder->byte_budget){
     LSP_break_pt = sp + 1;
     return;
   }
#endif
*/
      for( i = cxt_qtree.LSP_cxts - 1; i >= 0;
           LSP_models[i--].taub_scale(  ) );

    }

#else //SEPARATE_SETS_OF_PLANE_1

    last = ( bit_idx == max_bit_idx - 2 ) ?
      node_list.LSP - 1 : LSP_bit_idx_marks[bit_idx + 3];
    LSP_set_offset = 60;
    for( ; sp > last; sp-- ) {
      //decoding unit.................................
      r = *sp >> 16;
      c = *sp & 0xFFFF;
      LSP_offset = LSP_set_offset + LSP_cxt_tab[base_cxt[r][c]];

      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif
      //end of decoding unit.................................
    }

/*  June15
#ifdef TEST_CODING_RATE
   if(sub_decoder->bytes_used() > sub_decoder->byte_budget){
     LSP_break_pt = sp + 1;
     return;
   }
#endif
*/

#endif
  }

  if( bit_idx < max_bit_idx - 2 ) {

#ifndef MERGE_LOWER_LSP_PLANE
    LSP_set_offset = 70;

#else

#ifdef SCALE_AT_LSP_BREAK_PTS
    for( i = cxt_qtree.LSP_cxts - 1; i >= 0; LSP_models[i--].taub_scale(  ) );
#endif

#endif

    last = node_list.LSP - 1;
    for( ; sp > last; sp-- ) {
      //decoding unit.................................
      r = *sp >> 16;
      c = *sp & 0xFFFF;
      LSP_offset = LSP_set_offset + LSP_cxt_tab[base_cxt[r][c]];

      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif
      //end of decoding unit.................................
    }

/*  June15
#ifdef TEST_CODING_RATE
   if(sub_decoder->bytes_used() > sub_decoder->byte_budget){
     LSP_break_pt = sp + 1;
     return;
   }
#endif
*/
  }
#else
  int LSP_lev;
  if( bit_idx < max_bit_idx ) {
    LSP_lev = ( depth < 6 ) ? ( depth - 1 ) : 5;
    for( sp = LSP_bit_idx_marks[bit_idx + 1]; LSP_lev > 1; LSP_lev-- ) {
      LSP_offset = LSP_lev;
      for( last = LSP_ids[bit_idx + 1][LSP_lev - 1]; sp > last; sp-- ) {
        //decoding unit.................................
        r = *sp >> 16;
        c = *sp & 0xFFFF;

        if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
          base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
          LSP_break_pt = sp + 1;
          return;
        }
#endif


        //end of decoding unit.................................
      }
    }


    if( bit_idx < max_bit_idx - 1 ) {
      //from leaves
      last = LSP_ids[bit_idx + 1][0];
      LSP_offset = 1;
      for( ; sp > last; sp-- ) {
        //decoding unit.................................
        r = *sp >> 16;
        c = *sp & 0xFFFF;

        if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
          base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
          LSP_break_pt = sp + 1;
          return;
        }
#endif


        //end of decoding unit.................................
      }

      last = LSP_bit_idx_marks[bit_idx + 2];
      LSP_offset = 0;
      for( ; sp > last; sp-- ) {
        //decoding unit.................................
        r = *sp >> 16;
        c = *sp & 0xFFFF;

        if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
          base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
        if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
          LSP_break_pt = sp + 1;
          return;
        }
#endif
        //end of decoding unit.................................
      }
    }
  }


  if( bit_idx < max_bit_idx - 1 ) {
    last = ( bit_idx == max_bit_idx - 2 ) ?
      node_list.LSP - 1 : LSP_bit_idx_marks[bit_idx + 3];
    LSP_offset = 10;
    for( ; sp > last; sp-- ) {
      //decoding unit.................................
      r = *sp >> 16;
      c = *sp & 0xFFFF;

      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif

      //end of decoding unit.................................
    }
  }

  if( bit_idx < max_bit_idx - 2 ) {
    last = node_list.LSP - 1;
    LSP_offset = 11;
    for( ; sp > last; sp-- ) {
      //decoding unit.................................
      r = *sp >> 16;
      c = *sp & 0xFFFF;

      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif


      //end of decoding unit.................................
    }
  }
  //3 seperate levs.........................................

#ifdef 3LEVS_AC0
  if( bit_idx < max_bit_idx ) {
    last = ( bit_idx == max_bit_idx - 1 ) ?
      node_list.LSP - 1 : LSP_bit_idx_marks[bit_idx + 2];
    LSP_offset = 0;
    for( sp = LSP_bit_idx_marks[bit_idx + 1]; sp > last; sp-- ) {
      //decoding unit.................................
      r = *sp >> 16;
      c = *sp & 0xFFFF;

      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif


      //end of decoding unit.................................
    }
  }

  if( bit_idx < max_bit_idx - 1 ) {
    last = ( bit_idx == max_bit_idx - 2 ) ?
      node_list.LSP - 1 : LSP_bit_idx_marks[bit_idx + 3];
    LSP_offset = 1;
    for( ; sp > last; sp-- ) {
      //decoding unit.................................
      r = *sp >> 16;
      c = *sp & 0xFFFF;

      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif

      //end of decoding unit.................................
    }
  }

  if( bit_idx < max_bit_idx - 2 ) {
    last = node_list.LSP - 1;
    LSP_offset = 2;
    for( ; sp > last; sp-- ) {
      //decoding unit.................................
      r = *sp >> 16;
      c = *sp & 0xFFFF;

      if( sub_decoder->decode_symbol( LSP_models[LSP_offset] ) )
        base_coeff[r][c] |= bit_idx_mask;

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget ) {
        LSP_break_pt = sp + 1;
        return;
      }
#endif

      //end of decoding unit.................................
    }
  }
#endif //3LEVS_AC0

#endif //LSP_LUT

  ++( node_list.LSP_plane );
  LSP_break_pt = node_list.LSP;
}
