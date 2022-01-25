/* ========================================================================= */
/* Description: menber functions for class DecSubband                        */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <math.h>
#include <assert.h>
#include "dwt_bitplane_dec.h"


void
DecSubband::

initialize( SUBBAND_TYPE * subband,
            DecSubbandTree * dec_tree, DECODER_TYPE * dec )
{
  SubbandCodec::initialize( subband );
  dec_band_tree = dec_tree;
  sub_decoder = dec;
  create_max_and_cxt_qtrees(  );
  initialize_cxt_models(  );
  set_functions(  );
}

/*void DecSubband::
initialize2(SUBBAND_TYPE* subband,
                           DecSubbandTree *dec_tree)
{
  SubbandCodec::initialize(subband);
  dec_band_tree = dec_tree;


  NEW_VECTOR(sub_decoder, 1, DECODER_TYPE, "sub_decoder (dwt_bitplane_enc.cpp)") ; 

  create_max_and_cxt_qtrees();
  initialize_cxt_models();
  set_functions();
}
*/

void
DecSubband::reset_band_dec( SUBBAND_TYPE * subband )
{
  reset_band_codec( subband );
}

//-----------------------------------------------------------------------------

//  create_max_and_cxt_qtrees()

//-----------------------------------------------------------------------------
void
DecSubband::create_max_and_cxt_qtrees(  )
{
  NODE_CXT_TYPE ***cxt_nodes;
  Image_Coord_Sht *dims;
  int depth;

  Image_Coord band_dim = band->get_dim(  );
  int rows = band_dim.x;
  int cols = band_dim.y;

  depth = qtree.depth = Max( log2( rows - 1 ), log2( cols - 1 ) ) + 1;
  
  NEW_VECTOR( qtree.dims, depth + 1, Image_Coord_Sht, "qtree.nodes" );
  dims = qtree.dims;
  dims[0].r = rows;
  dims[0].c = cols;


  NEW_VECTOR( cxt_qtree.cxt_nodes, depth + 1, NODE_CXT_TYPE **,
              "qtree.nodes" );
  cxt_nodes = cxt_qtree.cxt_nodes;

  //compute dims
  int lev, r, c, cxt_r_sz, cxt_node_sz, node_sz;        //r_sz,
  Image_Coord_Sht *mem_dims;

  NEW_VECTOR( mem_dims, depth + 1, Image_Coord_Sht, "qtree.nodes" );

  node_sz = cxt_r_sz = cxt_node_sz = 0;
  for( lev = 1, r = rows, c = cols; lev <= depth; lev++ ) {
    dims[lev].r = r = 1 + ( ( r - 1 ) >> 1 );
    dims[lev].c = c = 1 + ( ( c - 1 ) >> 1 );

    mem_dims[lev].r = r = ( ( r + 1 ) >> 1 ) << 1;      // make r, c always even
    mem_dims[lev].c = c = ( ( c + 1 ) >> 1 ) << 1;

    cxt_r_sz += r + WIDTH_OF_CXT_BDY_X_2;
    node_sz += r * c;
    cxt_node_sz +=
      ( r + WIDTH_OF_CXT_BDY_X_2 ) * ( c + WIDTH_OF_CXT_BDY_X_2 );
  }
  qtree.qtree_sz = node_sz;
  qtree.nodes = NULL;
  cxt_qtree.cxt_qtree_sz = cxt_node_sz;

  NEW_VECTOR( cxt_nodes[1], cxt_r_sz, NODE_CXT_TYPE *, "NODE_CXT_TYPE*" );
  cxt_nodes[1] += WIDTH_OF_CXT_BDY;
  NEW_VECTOR( cxt_nodes[1][-WIDTH_OF_CXT_BDY], cxt_node_sz, NODE_CXT_TYPE,
              "NODE_CXT_TYPE" );
  cxt_nodes[1][-WIDTH_OF_CXT_BDY] += WIDTH_OF_CXT_BDY;

  for( lev = 2; lev <= depth; lev++ ) {
    r = mem_dims[lev - 1].r;
    c = mem_dims[lev - 1].c;
    cxt_nodes[lev] = cxt_nodes[lev - 1] + r + WIDTH_OF_CXT_BDY_X_2;
    cxt_nodes[lev][-WIDTH_OF_CXT_BDY] =
      cxt_nodes[lev - 1][-WIDTH_OF_CXT_BDY] + ( r +
                                                WIDTH_OF_CXT_BDY_X_2 ) * ( c +
                                                                           WIDTH_OF_CXT_BDY_X_2 );
  }

  for( lev = 1; lev <= depth; lev++ ) {
    int c_cxt = mem_dims[lev].c + WIDTH_OF_CXT_BDY_X_2;
    NODE_CXT_TYPE **node_cxt_dptr = cxt_nodes[lev] - WIDTH_OF_CXT_BDY;
    NODE_CXT_TYPE **node_cxt_end_dptr = node_cxt_dptr +
      mem_dims[lev].r + WIDTH_OF_CXT_BDY_X_2 - 1;
    for( ; node_cxt_dptr != node_cxt_end_dptr; node_cxt_dptr++ )
      node_cxt_dptr[1] = *node_cxt_dptr + c_cxt;
  }
  // even dim. all the time
  r = ( ( ( rows + 1 ) >> 1 ) << 1 ) + WIDTH_OF_CXT_BDY_X_2;
  c = ( ( ( cols + 1 ) >> 1 ) << 1 ) + WIDTH_OF_CXT_BDY_X_2;
  int rc = r * c;

  NEW_VECTOR( cxt_qtree.base_cxt, r, PEL_CXT_TYPE *, "base_cxt*" );
  cxt_qtree.base_cxt += WIDTH_OF_CXT_BDY;
  NEW_VECTOR( cxt_qtree.base_cxt[-WIDTH_OF_CXT_BDY], rc, PEL_CXT_TYPE,
              "base_cxt" );
  cxt_qtree.base_cxt[-WIDTH_OF_CXT_BDY] += WIDTH_OF_CXT_BDY;

  PEL_CXT_TYPE **pel_cxt_dptr = cxt_qtree.base_cxt - WIDTH_OF_CXT_BDY;
  PEL_CXT_TYPE **pel_cxt_end_dptr = pel_cxt_dptr + r - 1;
  for( ; pel_cxt_dptr != pel_cxt_end_dptr; pel_cxt_dptr++ )
    pel_cxt_dptr[1] = *pel_cxt_dptr + c;

  NEW_VECTOR( cxt_qtree.base_sign_cxt, r, SIGN_CXT_TYPE *, "base_sign_cxt*" );
  cxt_qtree.base_sign_cxt += WIDTH_OF_CXT_BDY;
  NEW_VECTOR( cxt_qtree.base_sign_cxt[-WIDTH_OF_CXT_BDY], rc, SIGN_CXT_TYPE,
              "base_sign_cxt" );
  cxt_qtree.base_sign_cxt[-WIDTH_OF_CXT_BDY] += WIDTH_OF_CXT_BDY;

  SIGN_CXT_TYPE **base_sign_cxt_dptr = cxt_qtree.base_sign_cxt -
    WIDTH_OF_CXT_BDY;
  SIGN_CXT_TYPE **base_sign_cxt_end_dptr = base_sign_cxt_dptr + r - 1;
  for( ; base_sign_cxt_dptr != base_sign_cxt_end_dptr; base_sign_cxt_dptr++ )
    base_sign_cxt_dptr[1] = *base_sign_cxt_dptr + c;
  //for inter subbands correlations

#ifdef INTERBANDS

  if( !( band->get_band_idx(  ) ) || ( band->get_band_level(  ) == 1 ) )
    child_cxt_qtree = NULL;     //bottom bands have no kids and the top level
  //only use info within subbands
  else
    child_cxt_qtree =
      &( dec_band_tree->dec_subs[band->get_child_band_idx(  )].cxt_qtree );

#endif


#ifdef GET_PARENT_MODELS
  int par_idx = band->get_par(  );
  if( ( par_idx < 0 )
      || ( ( par_idx == 0 ) && ( band->get_orientation(  ) == HH ) ) )
    par_cxt_qtree = NULL;
  else
    par_cxt_qtree = &( dec_band_tree->dec_subs[par_idx].cxt_qtree );

#endif

  delete[]mem_dims;

}

//-----------------------------------------------------------------------------

//  setup_cxt_qtrees()

//-----------------------------------------------------------------------------
void
DecSubband::setup_cxt_qtrees(  )
{
  int r, c, rc, lev;
  PEL_CXT_TYPE *pel_cxt_ptr;
  NODE_CXT_TYPE *node_cxt_ptr;
  NODE_CXT_TYPE ***cxt_nodes = cxt_qtree.cxt_nodes;
  PEL_CXT_TYPE **base_cxt = cxt_qtree.base_cxt;
  Image_Coord_Sht *dims = qtree.dims;
  int depth = qtree.depth;
  Image_Coord band_dim = band->get_dim(  );
  int rows = band_dim.x;
  int cols = band_dim.y;

  rc = ( ( ( ( rows + 1 ) >> 1 ) << 1 ) + WIDTH_OF_CXT_BDY_X_2 ) *
    ( ( ( ( cols + 1 ) >> 1 ) << 1 ) + WIDTH_OF_CXT_BDY_X_2 );

  pel_cxt_ptr = &( base_cxt[-WIDTH_OF_CXT_BDY][-WIDTH_OF_CXT_BDY] );
  PEL_CXT_TYPE *pel_cxt_end_ptr = pel_cxt_ptr + rc;
  while( pel_cxt_end_ptr != pel_cxt_ptr )
    *( pel_cxt_ptr++ ) = OUT_OF_BOUNDS;

  for( r = rows - 1; r >= 0; r-- ) {    //initialize cxt = 0 for pels within the bounds
    pel_cxt_ptr = base_cxt[r];
    memset( pel_cxt_ptr, 0, sizeof( PEL_CXT_TYPE ) * cols );
  }

  memset( &( cxt_qtree.base_sign_cxt[-WIDTH_OF_CXT_BDY][-WIDTH_OF_CXT_BDY] ),
          0, sizeof( SIGN_CXT_TYPE ) * rc );

  node_cxt_ptr = &( cxt_nodes[1][-WIDTH_OF_CXT_BDY][-WIDTH_OF_CXT_BDY] );
  NODE_CXT_TYPE *node_cxt_end_ptr = node_cxt_ptr + cxt_qtree.cxt_qtree_sz;
  while( node_cxt_ptr != node_cxt_end_ptr )
    *( node_cxt_ptr++ ) = OUT_OF_BOUNDS;

  //setup cxt
  for( lev = 1; lev <= depth; lev++ ) {
    for( c = dims[lev].c, r = dims[lev].r - 1; r >= 0; r-- ) {
      node_cxt_ptr = cxt_nodes[lev][r];
      memset( node_cxt_ptr, 0, sizeof( NODE_CXT_TYPE ) * c );
    }
  }
}

//-----------------------------------------------------------------------------

//  set_functions(void)

//-----------------------------------------------------------------------------
void
DecSubband::set_functions( void )
{


  decode_LIP = &DecSubband::decode_LIP_cxt_AC;
  decode_LIS_leaves = &DecSubband::decode_LIS_leaves_cxt_AC;
  decode_cur_qtree_level = &DecSubband::decode_cur_qtree_level_cxt_AC;

#ifdef LSP_BIT_IDX
  decode_LSP = &DecSubband::decode_LSP_cxt_AC_and_bit_idx;
#else
  decode_LSP = &DecSubband::decode_LSP_cxt_AC;
#endif

  decode_sig_node = &DecSubband::decode_sig_node_pos_dep_cxt_AC;
  decode_sig_leaf = &DecSubband::decode_sig_leaf_pos_dep_cxt_AC;
}

//-----------------------------------------------------------------------------

// decode_LIS_stack()

//-----------------------------------------------------------------------------
void
DecSubband::decode_LIS_stack(  )
{
  do {

//added by Peisong
    if( node_list.LIS_stack > node_list.LIS_stack_top ) {
      //cout<<"=======>LIS_stack="<<node_list.LIS_stack<<", LIS_stack_top="<<node_list.LIS_stack_top<<endl;               
      break;
    }
//added by Peisong


    if( node_list.LIS_stack_top->level > 1 ) {
      std_int node = node_list.LIS_stack_top->node;
      int level = ( node_list.LIS_stack_top-- )->level - 1;
      ( this->*decode_sig_node ) ( node, level );
      //      (this->*decode_sig_node)(node_list.LIS_stack_top->node,
      //                               (node_list.LIS_stack_top--)->level-1);
      // Rusert, 13.09.02
      
#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif
    } else {
      ( this->*decode_sig_leaf ) ( ( node_list.LIS_stack_top-- )->node );

#ifdef TEST_CODING_RATE         //June15
      if( sub_decoder->bytes_used(  ) > sub_decoder->byte_budget )
        return;
#endif

    }

  } while( node_list.LIS_stack <= node_list.LIS_stack_top );
}


//-----------------------------------------------------------------------------

//  decode_LIS_nodes()

//-----------------------------------------------------------------------------

void
DecSubband::decode_LIS_nodes(  )
{
#ifdef INITIALIZE_NODE_MODELS_FROM_PAR
  if( ( bit_idx == max_bit_idx - 1 ) && ( par_cxt_qtree ) ) {
    int lev, cxts, par_cxts;
    MODEL_TYPE *par_node_models, *node_models;

    cxts = par_cxts = 0;
    for( lev = qtree.depth - 2; lev > 1; lev-- ) {
      //parent band with depth = qtree.depth - 1
      cxts += cxt_qtree.node_cxts[lev];
      par_cxts += par_cxt_qtree->node_cxts[lev];
      assert( cxts == par_cxts );
    }
    node_models = cxt_qtree.cxt_models + cxt_qtree.node_offsets[2];
    par_node_models =
      par_cxt_qtree->cxt_models + par_cxt_qtree->node_offsets[2];
    for( int i = cxts - 1; i >= 0; ) {
      node_models[i].reset( par_node_models[i] );
      node_models[i--].taub_scale(  );
    }
  }
#endif

#ifdef SCALE_JSIG_MODELS_AT_LEVEL_ENDS
  int jsig_cxts = 0;
  MODEL_TYPE *jsig_models = cxt_qtree.cxt_models + cxt_qtree.sig_offsets[0];
  for( int lev = qtree.depth - 1; lev >= 0;
       jsig_cxts += cxt_qtree.sig_cxts[lev--] );
#endif //SCALE_JSIG_MODELS_AT_LEVEL_ENDS




#ifdef SCALE_NODE_MODELS_AT_LEVEL_ENDS
  int node_cxts;
  MODEL_TYPE *node_models;

  node_models = cxt_qtree.cxt_models + cxt_qtree.node_offsets[2];
  node_cxts = cxt_qtree.node_cxts[2];
#endif //SCALE_JSIG_MODELS_AT_LEVEL_ENDS

  cur_lev = 2;
  while( cur_lev < qtree.depth ) {

#ifdef SCALE_JSIG_MODELS_AT_LEVEL_ENDS
    for( int i = jsig_cxts - 1; i >= 0; jsig_models[i--].taub_scale(  ) );
#endif


#ifdef SCALE_NODE_MODELS_AT_LEVEL_ENDS
    for( int i = node_cxts - 1; i >= 0; node_models[i--].taub_scale(  ) );
#endif

    ( this->*decode_cur_qtree_level ) (  );
  }
}

//-----------------------------------------------------------------------------

//  start_dec_subband()

//-----------------------------------------------------------------------------
int
DecSubband::start_dec_subband(  )
{
  assert( band );
  delete_coding_stats(  );
  //for now
  //max_bit_idx = bit_idx = sub_decoder->decode_bits(4) - band->get_lsb();
  max_bit_idx = sub_decoder->decode_bits( 5 ) - band->get_lsb(  ) - EXTRA_BIT;
  //printf("max_bit_idx %d (dwt_bitplane_dec.cpp)\n", max_bit_idx);
  bit_idx = max_bit_idx + 1;
  no_bitplanes = max_bit_idx - EXTRA_BIT + 1;
  if( no_bitplanes <= 0 ) {
    no_bitplanes = 1;
  }                             //printf("error\n");
  min_bit_idx = EXTRA_BIT;

  band->set_max_msb( max_bit_idx + band->get_lsb(  ) );
  setup_cxt_qtrees(  );
  initialize_node_list(  );
  create_coding_stats(  );
  LSP_break_pt = node_list.LSP;

  return 0;
}


//-----------------------------------------------------------------------------

//            rec_subband()

//-----------------------------------------------------------------------------

void
DecSubband::rec_subband(  )
{
  std_int *last, *sp;
  SUB_COEFF_TYPE mask = bit_idx_mask >> 1;

  for( sp = node_list.LSP_end; sp >= LSP_break_pt; sp-- )
    base_coeff[*sp >> 16][*sp & 0xFFFF] |= mask;

  mask = bit_idx_mask;

  for( last = node_list.LSP; sp >= last; sp-- )
    base_coeff[*sp >> 16][*sp & 0xFFFF] |= mask;

}
