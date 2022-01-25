/* ========================================================================= */
/* Description: member functions for class EncSubband                        */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <math.h>
#include <assert.h>
#include <string.h>
#include "dwt_bitplane_enc.h"

void
EncSubband::

initialize( SUBBAND_TYPE * subband, EncSubbandTree * enc_tree,
            ENCODER_TYPE * enc )
{
  SubbandCodec::initialize( subband );
  enc_band_tree = enc_tree;
  sub_encoder = enc;

  create_max_and_cxt_qtrees(  );
  initialize_cxt_models(  );
  set_functions(  );

}

/*void EncSubband::
initialize2(SUBBAND_TYPE* subband, EncSubbandTree *enc_tree)
{
  SubbandCodec::initialize(subband);
  enc_band_tree = enc_tree;

  NEW_VECTOR(sub_encoder, 1, ENCODER_TYPE, "sub_encoder (dwt_bitplane_enc.cpp)") ; 

  sub_encoder->new_open_file(); 

  create_max_and_cxt_qtrees();
  initialize_cxt_models();
  set_functions();
}*/


void
EncSubband::reset_band_enc( SUBBAND_TYPE * subband )
{
  reset_band_codec( subband );
}

//-----------------------------------------------------------------------------

//  create_max_and_cxt_qtrees()

//    The dims of qtree.nodes are always even.
//    The dims of cxt_qtree.nodes = qtree.nodes + WIDTH_OF_CXT_BDY_X_2.
//    The results of this def. is that the dim of the root = 2 x 2
//
//-----------------------------------------------------------------------------
void
EncSubband::create_max_and_cxt_qtrees(  )
{
  SUB_COEFF_TYPE ***nodes;
  NODE_CXT_TYPE ***cxt_nodes;
  Image_Coord_Sht *dims;
  int depth;

  Image_Coord band_dim = band->get_dim(  );
  int rows = band_dim.x;
  int cols = band_dim.y;

  depth = qtree.depth = Max( log2( rows - 1 ), log2( cols - 1 ) ) + 1;

  assert( depth < 15 );         //dim is only short int

  // 0 - depth
  NEW_VECTOR( qtree.dims, depth + 1, Image_Coord_Sht, "qtree.nodes" );
  NEW_VECTOR( qtree.nodes, depth + 1, SUB_COEFF_TYPE **, "qtree.nodes" );

  dims = qtree.dims;
  nodes = qtree.nodes;

  nodes[0] = band->get_coeff(  );
  dims[0].r = rows;
  dims[0].c = cols;

  NEW_VECTOR( cxt_qtree.cxt_nodes, depth + 1, NODE_CXT_TYPE **,
              "qtree.nodes" );
  cxt_nodes = cxt_qtree.cxt_nodes;

  //compute dims
  int lev, r, c, cxt_r_sz, r_sz, cxt_node_sz, node_sz;
  Image_Coord_Sht *mem_dims;

  NEW_VECTOR( mem_dims, depth + 1, Image_Coord_Sht, "qtree.nodes" );

  r_sz = node_sz = cxt_r_sz = cxt_node_sz = 0;
  for( lev = 1, r = rows, c = cols; lev <= depth; lev++ ) {
    dims[lev].r = r = 1 + ( ( r - 1 ) >> 1 );
    dims[lev].c = c = 1 + ( ( c - 1 ) >> 1 );

    mem_dims[lev].r = r = ( ( r + 1 ) >> 1 ) << 1;      // make r, c always even
    mem_dims[lev].c = c = ( ( c + 1 ) >> 1 ) << 1;
    r_sz += r;
    cxt_r_sz += r + WIDTH_OF_CXT_BDY_X_2;
    node_sz += r * c;
    cxt_node_sz +=
      ( r + WIDTH_OF_CXT_BDY_X_2 ) * ( c + WIDTH_OF_CXT_BDY_X_2 );
  }

  qtree.qtree_sz = node_sz;
  cxt_qtree.cxt_qtree_sz = cxt_node_sz;

  //memory allocation

  NEW_VECTOR( nodes[1], r_sz, SUB_COEFF_TYPE *, "qtree.nodes[1]" );
  NEW_VECTOR( nodes[1][0], node_sz, SUB_COEFF_TYPE, "qtree.nodes[1][0]" );

  NEW_VECTOR( cxt_nodes[1], cxt_r_sz, NODE_CXT_TYPE *, "NODE_CXT_TYPE*" );
  cxt_nodes[1] += WIDTH_OF_CXT_BDY;
  NEW_VECTOR( cxt_nodes[1][-WIDTH_OF_CXT_BDY], cxt_node_sz, NODE_CXT_TYPE,
              "NODE_CXT_TYPE" );
  cxt_nodes[1][-WIDTH_OF_CXT_BDY] += WIDTH_OF_CXT_BDY;

  for( lev = 2; lev <= depth; lev++ ) {
    r = mem_dims[lev - 1].r;
    c = mem_dims[lev - 1].c;
    nodes[lev] = nodes[lev - 1] + r;
    nodes[lev][0] = nodes[lev - 1][0] + r * c;

    cxt_nodes[lev] = cxt_nodes[lev - 1] + r + WIDTH_OF_CXT_BDY_X_2;
    cxt_nodes[lev][-WIDTH_OF_CXT_BDY] =
      cxt_nodes[lev - 1][-WIDTH_OF_CXT_BDY] + ( r +
                                                WIDTH_OF_CXT_BDY_X_2 ) * ( c +
                                                                           WIDTH_OF_CXT_BDY_X_2 );
  }


  for( lev = 1; lev <= depth; lev++ ) {
    c = mem_dims[lev].c;
    r = mem_dims[lev].r;
    SUB_COEFF_TYPE **node_dptr = nodes[lev];
    SUB_COEFF_TYPE **node_end_dptr = node_dptr + r - 1;
    for( ; node_dptr != node_end_dptr; node_dptr++ )
      node_dptr[1] = *node_dptr + c;

    int c_cxt = c + WIDTH_OF_CXT_BDY_X_2;
    NODE_CXT_TYPE **node_cxt_dptr = cxt_nodes[lev] - WIDTH_OF_CXT_BDY;
    NODE_CXT_TYPE **node_cxt_end_dptr = node_cxt_dptr
      + r + WIDTH_OF_CXT_BDY_X_2 - 1;

    for( ; node_cxt_dptr != node_cxt_end_dptr; node_cxt_dptr++ )
      node_cxt_dptr[1] = *node_cxt_dptr + c_cxt;
  }

  //memory allocation for pel context

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
      &( enc_band_tree->enc_subs[band->get_child_band_idx(  )].cxt_qtree );

#endif


#ifdef GET_PARENT_MODELS
  int par_idx = band->get_par(  );
  if( ( par_idx < 0 )
      || ( ( par_idx == 0 ) && ( band->get_orientation(  ) == HH ) ) )
    par_cxt_qtree = NULL;
  else
    par_cxt_qtree = &( enc_band_tree->enc_subs[par_idx].cxt_qtree );

#endif

  delete[]mem_dims;
}

//-----------------------------------------------------------------------------

//void EncSubband::setup_max_and_cxt_qtrees()

//-----------------------------------------------------------------------------
void
EncSubband::setup_max_and_cxt_qtrees(  )
{
  int r, c, rc, lev;
  PEL_CXT_TYPE *pel_cxt_ptr;
  NODE_CXT_TYPE *node_cxt_ptr;
  SUB_COEFF_TYPE *node_ptr, *cnode_ptr;

  SUB_COEFF_TYPE ***nodes = qtree.nodes;
  NODE_CXT_TYPE ***cxt_nodes = cxt_qtree.cxt_nodes;
  PEL_CXT_TYPE **base_cxt = cxt_qtree.base_cxt;
  Image_Coord_Sht *dims = qtree.dims;
  int depth = qtree.depth;

  Image_Coord band_dim = band->get_dim(  );
  int rows = band_dim.x;
  int cols = band_dim.y;

  rc = ( ( ( ( rows + 1 ) >> 1 ) << 1 ) + WIDTH_OF_CXT_BDY_X_2 ) *
    ( ( ( ( cols + 1 ) >> 1 ) << 1 ) + WIDTH_OF_CXT_BDY_X_2 );

  //pel cxt state initialization
  pel_cxt_ptr = &( base_cxt[-WIDTH_OF_CXT_BDY][-WIDTH_OF_CXT_BDY] );

  PEL_CXT_TYPE *pel_cxt_end_ptr = pel_cxt_ptr + rc;

  while( pel_cxt_end_ptr != pel_cxt_ptr )
    *( pel_cxt_ptr++ ) = OUT_OF_BOUNDS;

  for( r = rows - 1; r >= 0; r-- ) {    //initialize cxt = 0 for pels within the bounds
    memset( base_cxt[r], 0, sizeof( PEL_CXT_TYPE ) * cols );
  }

  // OUT_OF_BOUNDS for all node_cxt
  node_cxt_ptr = &( cxt_nodes[1][-WIDTH_OF_CXT_BDY][-WIDTH_OF_CXT_BDY] );
  NODE_CXT_TYPE *node_cxt_end_ptr = node_cxt_ptr + cxt_qtree.cxt_qtree_sz;
  while( node_cxt_ptr != node_cxt_end_ptr )
    *( node_cxt_ptr++ ) = OUT_OF_BOUNDS;

  //reset node states
  memset( nodes[1][0], 0, qtree.qtree_sz * sizeof( SUB_COEFF_TYPE ) );
  nodes[0] = base_coeff;

  for( lev = 1; lev <= depth; lev++ ) {
    c = dims[lev].c;
    int crow_gap = nodes[lev - 1][1] - nodes[lev - 1][0];
    for( r = dims[lev].r - 1; r >= 0; r-- ) {
      node_ptr = nodes[lev][r];
      cnode_ptr = nodes[lev - 1][r << 1];
      node_cxt_ptr = cxt_nodes[lev][r];
      for( SUB_COEFF_TYPE * node_end_ptr = node_ptr + c;
           node_ptr != node_end_ptr; cnode_ptr += 2 ) {
        *( node_ptr++ ) = ( cnode_ptr[0] | cnode_ptr[1] | cnode_ptr[crow_gap]
                            | cnode_ptr[crow_gap + 1] ) & MAG_MASK;
        *( node_cxt_ptr++ ) = 0;
      }
    }
  }
  // sign cxt state initialization
  memset( &( cxt_qtree.base_sign_cxt[-WIDTH_OF_CXT_BDY][-WIDTH_OF_CXT_BDY] ),
          0, sizeof( SIGN_CXT_TYPE ) * rc );
}

//-----------------------------------------------------------------------------

//  set_functions(void)

//-----------------------------------------------------------------------------
void
EncSubband::set_functions( void )
{

  encode_LIP = &EncSubband::encode_LIP_cxt_AC;
  encode_LIS_leaves = &EncSubband::encode_LIS_leaves_cxt_AC;
  encode_cur_qtree_level = &EncSubband::encode_cur_qtree_level_cxt_AC;

#ifdef LSP_BIT_IDX
  encode_LSP = &EncSubband::encode_LSP_cxt_AC_and_bit_idx;
#else
  encode_LSP = &EncSubband::encode_LSP_cxt_AC;

#endif

  encode_sig_leaf = &EncSubband::encode_sig_leaf_pos_dep_cxt_AC;
  encode_sig_node = &EncSubband::encode_sig_node_pos_dep_cxt_AC;
}

//-----------------------------------------------------------------------------

//     start_enc_subband()

//-----------------------------------------------------------------------------
int
EncSubband::start_enc_subband(  )
{

  assert( band );
  delete_coding_stats(  );

  setup_max_and_cxt_qtrees(  );
  max_bit_idx = band->get_max_msb(  ) - band->get_lsb(  );
  bit_idx = max_bit_idx + 1;
  no_bitplanes = max_bit_idx - EXTRA_BIT + 1;
  if( no_bitplanes <= 0 ) {
    no_bitplanes = 1;
  }                             //printf("error\n");
  min_bit_idx = EXTRA_BIT;

  initialize_node_list(  );
  create_coding_stats(  );

  //for now
  sub_encoder->code_bits( 5, band->get_max_msb(  ) + EXTRA_BIT );
  //printf("band->get_max_msb %d \n", band->get_max_msb());
  return 0;
}


//-----------------------------------------------------------------------------

//   enc_LIS_stack()

//-----------------------------------------------------------------------------

void
EncSubband::encode_LIS_stack(  )
{

  do {
    if( node_list.LIS_stack_top->level > 1 ) {
       std_int node = node_list.LIS_stack_top->node;
       int level = ( node_list.LIS_stack_top-- )->level - 1;
       ( this->*encode_sig_node ) ( node, level );
 
       //(this->*encode_sig_node)(node_list.LIS_stack_top->node,
       //(node_list.LIS_stack_top--)->level-1);
       // Rusert, 13.09.02
    }
    else
      ( this->*encode_sig_leaf ) ( ( node_list.LIS_stack_top-- )->node );
  } while( node_list.LIS_stack <= node_list.LIS_stack_top );
}

//-----------------------------------------------------------------------------

// encode_LIS_nodes()

//-----------------------------------------------------------------------------
void
EncSubband::encode_LIS_nodes(  )
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
#endif //SCALE_NODE_MODELS_AT_LEVEL_ENDS

  cur_lev = 2;

  while( cur_lev < qtree.depth ) {

#ifdef SCALE_JSIG_MODELS_AT_LEVEL_ENDS
    for( int i = jsig_cxts - 1; i >= 0;
         // jsig_models[i--].scale());
         jsig_models[i--].taub_scale(  ) );
#endif


#ifdef SCALE_NODE_MODELS_AT_LEVEL_ENDS
    for( int i = node_cxts - 1; i >= 0; node_models[i--].taub_scale(  ) );
#endif

    ( this->*encode_cur_qtree_level ) (  );
  }
}
