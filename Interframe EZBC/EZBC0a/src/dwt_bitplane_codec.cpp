/* ========================================================================= */
/* Description: menber functions for class SubbandCodec                      */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */
#include <math.h>
#include <assert.h>
#include "dwt_bitplane_codec.h"

#ifdef NDEBUG
#define myassert(seq) seq
#else
#define myassert(seq) assert(seq)
#endif

//long  SubbandCodec::BYTE_BUDGET = MAX_STD_INT;

SubbandCodec::~SubbandCodec(  )
{

  delete_coding_stats(  );

  if( qtree.nodes ) {
    delete[]qtree.nodes[1][0];
    delete[]qtree.nodes[1];
    DELETE_VECTOR( qtree.nodes );
  }

  if( cxt_qtree.cxt_nodes ) {
    cxt_qtree.cxt_nodes[1][-WIDTH_OF_CXT_BDY] -= WIDTH_OF_CXT_BDY;
    delete[]cxt_qtree.cxt_nodes[1][-WIDTH_OF_CXT_BDY];
    cxt_qtree.cxt_nodes[1] -= WIDTH_OF_CXT_BDY;
    delete[]cxt_qtree.cxt_nodes[1];
    DELETE_VECTOR( cxt_qtree.cxt_nodes );
  }
  if( cxt_qtree.base_cxt ) {
    delete[]( cxt_qtree.base_cxt[-WIDTH_OF_CXT_BDY] - WIDTH_OF_CXT_BDY );
    cxt_qtree.base_cxt -= WIDTH_OF_CXT_BDY;
    DELETE_VECTOR( cxt_qtree.base_cxt );
  }
  if( cxt_qtree.base_sign_cxt ) {
    delete[]( cxt_qtree.base_sign_cxt[-WIDTH_OF_CXT_BDY] - WIDTH_OF_CXT_BDY );
    cxt_qtree.base_sign_cxt -= WIDTH_OF_CXT_BDY;
    DELETE_VECTOR( cxt_qtree.base_sign_cxt );
  }

  clear_node_list(  );

  DELETE_VECTOR( qtree.dims );
  DELETE_VECTOR( cxt_qtree.sig_cxts );
  DELETE_VECTOR( cxt_qtree.sig_offsets );
  DELETE_VECTOR( cxt_qtree.sig_tabs );

  if( cxt_qtree.node_cxts ) {
    cxt_qtree.node_cxts++;
    DELETE_VECTOR( cxt_qtree.node_cxts );
  }
  if( cxt_qtree.node_offsets ) {
    cxt_qtree.node_offsets++;
    DELETE_VECTOR( cxt_qtree.node_offsets );
  }
  if( cxt_qtree.node_tabs ) {
    cxt_qtree.node_tabs++;
    DELETE_VECTOR( cxt_qtree.node_tabs );
  }
  DELETE_VECTOR( cxt_qtree.cxt_models );

  if( cxt_qtree.jsig_cxts ) {
    delete[]cxt_qtree.jsig_cxts[0];
    DELETE_VECTOR( cxt_qtree.jsig_cxts );
  }
  if( cxt_qtree.jsig_offsets ) {
    delete[]cxt_qtree.jsig_offsets[0];
    DELETE_VECTOR( cxt_qtree.jsig_offsets );
  }
  if( cxt_qtree.jsig_tabs ) {
    delete[]cxt_qtree.jsig_tabs[0];
    DELETE_VECTOR( cxt_qtree.jsig_tabs );
  }
}

void
SubbandCodec::delete_coding_stats( void )
{
  if( coding_stats.bitplanes ) {
    for( int n = 0; n < no_bitplanes; n++ )
      delete[]coding_stats.bitplanes[EXTRA_BIT + n].passes;
    coding_stats.bitplanes += EXTRA_BIT;
    DELETE_VECTOR( coding_stats.bitplanes );
  }
}

//-----------------------------------------------------------------------------

//     initialize()

//-----------------------------------------------------------------------------
int
SubbandCodec::initialize( SUBBAND_TYPE * subband )
{
  memset( this, 0, sizeof( SubbandCodec ) );

  //assert(band = subband);
  myassert( band = subband );
  if( subband->get_orientation(  ) == HL )
    subband->transpose(  );
  base_coeff = band->get_coeff(  );

  fpout = stdout;

  return 0;
}

void
SubbandCodec::reset_band_codec( SUBBAND_TYPE * subband )
{
  //assert(band = subband);
  myassert( band = subband );
  if( subband->get_orientation(  ) == HL )
    subband->transpose(  );
  base_coeff = band->get_coeff(  );
  cur_pass = cur_lev = 0;
  for( int i = cxt_qtree.total_cxts - 1; i >= 0;
       cxt_qtree.cxt_models[i--].reset(  ) );
}

//-----------------------------------------------------------------------------

//             initialize_node_list()

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_node_list( void )
{
  clear_node_list(  );

  int i, depth = qtree.depth;
  Image_Coord_Sht *dims = qtree.dims;
  std_int band_size = dims[0].r * dims[0].c;
  int tree_sz = qtree.qtree_sz;


  node_list.LIP_sz = ( std_int ) band_size;
  NEW_VECTOR( node_list.LIP, band_size, std_int, "node_list.LIP" );
  node_list.LSP = node_list.LIP;
  NEW_VECTOR( node_list.LIS, depth, std_int *, "node_list.LIS" );
  ( node_list.LIS )--;

  std_int **LIS = node_list.LIS;

  NEW_VECTOR( LIS[1], tree_sz, std_int, "node_list.LIS[1]" );
  int lev;
  for( lev = 1; lev < depth; lev++ ) {
    LIS[lev + 1] = LIS[lev] + dims[lev].r * dims[lev].c;
  }
  NEW_VECTOR( node_list.LIS_end, depth, std_int *, "node_list.LIS_end" );
  ( node_list.LIS_end )--;
  for( lev = 1; lev <= depth; lev++ ) {
    node_list.LIS_end[lev] = node_list.LIS[lev] - 1;
  }
  //*(++node_list.LIS_end[depth]) = 0;  //root  node

  NEW_VECTOR( node_list.LIS_prev_mark, depth, std_int *,
              "node_list.LIS_prev_mark" );
  ( node_list.LIS_prev_mark )--;

  NEW_VECTOR( node_list.LIS_old_end, depth, std_int *,
              "node_list.LIS_old_end" );
  ( node_list.LIS_old_end )--;

  NEW_VECTOR( node_list.LIS_stack, depth << 2, Node_Entry,
              "node_list.LIS_stack" );

  node_list.LIS_stack_top = node_list.LIS_stack - 1;

  node_list.LSP_plane = 0;
  node_list.LSP_end = node_list.LSP - 1;
  node_list.LIP_end = node_list.LIP + band_size;

  NEW_VECTOR( node_list.LSP_bit_idx_marks, no_bitplanes, std_int *,
              "LSP_bit_idx_marks" );
  memset( node_list.LSP_bit_idx_marks, 0,
          sizeof( std_int * ) * no_bitplanes );
  node_list.LSP_bit_idx_marks -= EXTRA_BIT;

  NEW_VECTOR( node_list.LSP_ids, no_bitplanes, std_int **, "LSP_ids" );
  node_list.LSP_ids -= EXTRA_BIT;

  NEW_VECTOR( node_list.LSP_ids[EXTRA_BIT], no_bitplanes * depth,
              std_int *, "LSP_ids[bit_idx]" );
  memset( node_list.LSP_ids[EXTRA_BIT], 0,
          sizeof( std_int * ) * no_bitplanes * depth );
  for( i = EXTRA_BIT; i < max_bit_idx; i++ )
    node_list.LSP_ids[i + 1] = node_list.LSP_ids[i] + depth;

}

//-----------------------------------------------------------------------------

//        clear_node_list()

//-----------------------------------------------------------------------------
void
SubbandCodec::clear_node_list(  )
{

  DELETE_VECTOR( node_list.LIP );
  if( node_list.LIS ) {
    DELETE_VECTOR( node_list.LIS[1] );
    ( node_list.LIS )++;
    DELETE_VECTOR( node_list.LIS );
  }
  if( node_list.LIS_end ) {
    ( node_list.LIS_end )++;
    DELETE_VECTOR( node_list.LIS_end );
  }
  if( node_list.LIS_prev_mark ) {
    ( node_list.LIS_prev_mark )++;
    DELETE_VECTOR( node_list.LIS_prev_mark );
  }
  if( node_list.LIS_old_end ) {
    ( node_list.LIS_old_end )++;
    DELETE_VECTOR( node_list.LIS_old_end );
  }
  DELETE_VECTOR( node_list.LIS_stack );

  if( node_list.LSP_ids ) {
    //delete [] node_list.LSP_ids[EXTRA_BIT];
    DELETE_VECTOR( node_list.LSP_ids[EXTRA_BIT] );      //Peisong


    node_list.LSP_ids += EXTRA_BIT;
    DELETE_VECTOR( node_list.LSP_ids );
  }
  if( node_list.LSP_bit_idx_marks ) {
    node_list.LSP_bit_idx_marks += EXTRA_BIT;
    DELETE_VECTOR( node_list.LSP_bit_idx_marks );
  }
}

//-----------------------------------------------------------------------------

//           create_coding_stats()

//-----------------------------------------------------------------------------
void
SubbandCodec::create_coding_stats(  )
{

  // delete_coding_stats();

#ifdef LSP_BIT_IDX
  no_passes = qtree.depth + 5;  //LIP, 1 ~ depth-1, LSP 0-4
#else
  no_passes = qtree.depth + 1;  //LIP, 1 ~ depth-1, LSP
#endif
  //for now
  no_passes = 32;

  NEW_VECTOR( coding_stats.bitplanes, no_bitplanes, BitPlaneStats,
              "plane_stats" );

  coding_stats.bitplanes -= EXTRA_BIT;

  for( int n = 0; n < no_bitplanes; n++ ) {
    NEW_VECTOR( coding_stats.bitplanes[EXTRA_BIT + n].passes,
                no_passes, PassStats, "PassStat" );
    memset( coding_stats.bitplanes[EXTRA_BIT + n].passes, 0,
            sizeof( PassStats ) * no_passes );
  }
}

//-----------------------------------------------------------------------------

//  initialize_cxt_models()

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_cxt_models(  )
{
  int i, depth = qtree.depth;

  //0 ~ depth-1
  NEW_VECTOR( cxt_qtree.sig_cxts, depth, int, "sig_cxts" );
  NEW_VECTOR( cxt_qtree.sig_offsets, depth, int, "sig_offsets" );
  NEW_VECTOR( cxt_qtree.sig_tabs, depth, std_byte *, "sig_tabs" );

  //0 ~ depth-1
  NEW_VECTOR( cxt_qtree.jsig_cxts, depth, int *, "sig_cxts" );
  NEW_VECTOR( cxt_qtree.jsig_offsets, depth, int *, "sig_offsets" );
  NEW_VECTOR( cxt_qtree.jsig_tabs, depth, std_byte **, "sig_tabs" );

  NEW_VECTOR( cxt_qtree.jsig_cxts[0], depth << 2, int, "sig_cxts" );
  NEW_VECTOR( cxt_qtree.jsig_offsets[0], depth << 2, int, "sig_offsets" );
  NEW_VECTOR( cxt_qtree.jsig_tabs[0], depth << 2, std_byte *, "sig_tabs" );

  for( i = 1; i < depth; i++ ) {
    cxt_qtree.jsig_cxts[i] = cxt_qtree.jsig_cxts[i - 1] + 4;
    cxt_qtree.jsig_offsets[i] = cxt_qtree.jsig_offsets[i - 1] + 4;
    cxt_qtree.jsig_tabs[i] = cxt_qtree.jsig_tabs[i - 1] + 4;
  }
  //1 ~ depth-1, lev0 is represented by LIP
  NEW_VECTOR( cxt_qtree.node_cxts, depth - 1, int, "node_cxts" );
  cxt_qtree.node_cxts--;
  NEW_VECTOR( cxt_qtree.node_offsets, depth - 1, int, "node_offsets" );
  cxt_qtree.node_offsets--;
  NEW_VECTOR( cxt_qtree.node_tabs, depth - 1, std_byte *, "node_tabs" );
  cxt_qtree.node_tabs--;

  create_cxt_models(  );

  int total = cxt_qtree.total_cxts;

  NEW_VECTOR( cxt_qtree.cxt_models, total, MODEL_TYPE, "MODEL_TYPE" );

  for( i = total - 1; i >= 0; cxt_qtree.cxt_models[i--].set_symbols( 2 ) );
}


//-----------------------------------------------------------------------------

//  create_cxt_models()

//-----------------------------------------------------------------------------
void
SubbandCodec::create_cxt_models( void )
{
  int i, lev, total;
  int node_zc_cxts, LIP_zc_cxts, LSP_cxts, sign_cxts;
  std_byte *node_zc_lut, *LIP_zc_lut, *LSP_lut;
  SUB_COEFF_TYPE *sign_lut;
  // Image_Coord_Sht *dims = qtree.dims;

  int sig_zc_cxts[4], sig0_zc_cxts[4];
  std_byte *sig_zc_lut[4], *sig0_zc_lut[4];

  int depth = qtree.depth;

  if( band->get_orientation(  ) == HH ) {

    //assert(LIP_zc_lut = sm0_sp_diag_node_lut);
    myassert( LIP_zc_lut = sm0_sp_diag_node_lut );
    LIP_zc_cxts = sm0_sp_diag_node_cxts;

    //assert(node_zc_lut = mid_inter_and_sp_diag_node_lut);
    myassert( node_zc_lut = mid_inter_and_sp_diag_node_lut );
    node_zc_cxts = mid_inter_and_sp_diag_node_cxts;

    //assert(sig0_zc_lut[0] = mid_diag_jsig0_00_lut);
    myassert( sig0_zc_lut[0] = mid_diag_jsig0_00_lut );
    sig0_zc_cxts[0] = mid_diag_jsig0_00_cxts;

    //assert(sig_zc_lut[0] = mid_diag_jsig_00_lut);
    myassert( sig_zc_lut[0] = mid_diag_jsig_00_lut );
    sig_zc_cxts[0] = mid_diag_jsig_00_cxts;


    //assert(sig0_zc_lut[1] = mid_diag_jsig0_01_lut);
    myassert( sig0_zc_lut[1] = mid_diag_jsig0_01_lut );
    sig0_zc_cxts[1] = mid_diag_jsig0_01_cxts;
    //assert(sig_zc_lut[1] = mid_diag_jsig_01_lut);
    myassert( sig_zc_lut[1] = mid_diag_jsig_01_lut );
    sig_zc_cxts[1] = mid_diag_jsig_01_cxts;


    //assert(sig0_zc_lut[2] = mid_diag_jsig0_10_lut);
    myassert( sig0_zc_lut[2] = mid_diag_jsig0_10_lut );
    sig0_zc_cxts[2] = mid_diag_jsig0_10_cxts;
    //assert(sig_zc_lut[2] = mid_diag_jsig_10_lut);
    myassert( sig_zc_lut[2] = mid_diag_jsig_10_lut );
    sig_zc_cxts[2] = mid_diag_jsig_10_cxts;

    //assert(sig0_zc_lut[3] = mid_diag_jsig0_11_lut);
    myassert( sig0_zc_lut[3] = mid_diag_jsig0_11_lut );
    sig0_zc_cxts[3] = mid_diag_jsig0_11_cxts;
    //assert(sig_zc_lut[3] = mid_diag_jsig_11_lut);
    myassert( sig_zc_lut[3] = mid_diag_jsig_11_lut );
    sig_zc_cxts[3] = mid_diag_jsig_11_cxts;

  } else {

    //assert(LIP_zc_lut = sm0_sp_main_node_lut);
    myassert( LIP_zc_lut = sm0_sp_main_node_lut );
    LIP_zc_cxts = sm0_sp_main_node_cxts;

    //assert(node_zc_lut = mid_inter_and_sp_main_node_lut);
    myassert( node_zc_lut = mid_inter_and_sp_main_node_lut );
    node_zc_cxts = mid_inter_and_sp_main_node_cxts;

    //assert(sig0_zc_lut[0] = mid_main_jsig0_00_lut);
    myassert( sig0_zc_lut[0] = mid_main_jsig0_00_lut );
    sig0_zc_cxts[0] = mid_main_jsig0_00_cxts;
    //assert(sig_zc_lut[0] = mid_main_jsig_00_lut);
    myassert( sig_zc_lut[0] = mid_main_jsig_00_lut );
    sig_zc_cxts[0] = mid_main_jsig_00_cxts;


    //assert(sig0_zc_lut[1] = mid_main_jsig0_01_lut);
    myassert( sig0_zc_lut[1] = mid_main_jsig0_01_lut );
    sig0_zc_cxts[1] = mid_main_jsig0_01_cxts;
    //assert(sig_zc_lut[1] = mid_main_jsig_01_lut);
    myassert( sig_zc_lut[1] = mid_main_jsig_01_lut );
    sig_zc_cxts[1] = mid_main_jsig_01_cxts;

    //assert(sig0_zc_lut[2] = mid_main_jsig0_10_lut);
    myassert( sig0_zc_lut[2] = mid_main_jsig0_10_lut );
    sig0_zc_cxts[2] = mid_main_jsig0_10_cxts;
    //assert(sig_zc_lut[2] = mid_main_jsig_10_lut);
    myassert( sig_zc_lut[2] = mid_main_jsig_10_lut );
    sig_zc_cxts[2] = mid_main_jsig_10_cxts;

    //assert(sig0_zc_lut[3] = mid_main_jsig0_11_lut);
    myassert( sig0_zc_lut[3] = mid_main_jsig0_11_lut );
    sig0_zc_cxts[3] = mid_main_jsig0_11_cxts;
    //assert(sig_zc_lut[3] = mid_main_jsig_11_lut);
    myassert( sig_zc_lut[3] = mid_main_jsig_11_lut );
    sig_zc_cxts[3] = mid_main_jsig_11_cxts;

  }

  //assignments of LUTS for LSP coding

#ifdef LSP_LUT

  if( band->get_orientation(  ) == HH ) {

#ifdef TWO_SIGNIF_BITS
    //assert(LSP_lut = sm3_2_signif_bits_diag_LSP_lut);
    myassert( LSP_lut = sm3_2_signif_bits_diag_LSP_lut );
    LSP_cxts = sm3_2_signif_bits_diag_LSP_cxts;
#else
    //assert(LSP_lut = mid_diag_LSP_lut);
    myassert( LSP_lut = mid_diag_LSP_lut );
    LSP_cxts = mid_diag_LSP_cxts;
#endif
  } else {

#ifdef TWO_SIGNIF_BITS
    //assert(LSP_lut = sm3_2_signif_bits_main_LSP_lut);
    myassert( LSP_lut = sm3_2_signif_bits_main_LSP_lut );
    LSP_cxts = sm3_2_signif_bits_main_LSP_cxts;
#else
    //assert(LSP_lut = mid_main_LSP_lut);
    myassert( LSP_lut = mid_main_LSP_lut );
    LSP_cxts = mid_main_LSP_cxts;
#endif
  }
#endif //LSP_LUT

  //sign luts

  if( band->get_orientation(  ) == HH ) {
    //assert(sign_lut = mid_diag_sign_lut);
    myassert( sign_lut = mid_diag_sign_lut );
    sign_cxts = mid_diag_sign_cxts;

  } else {
    //assert(sign_lut = mid_main_sign_lut);
    myassert( sign_lut = mid_main_sign_lut );
    sign_cxts = mid_main_sign_cxts;
  }

  cxt_qtree.sig_offsets[0] = total = 0;
  for( i = 0; i < 4; i++ ) {
    cxt_qtree.jsig_offsets[0][i] = total;
    cxt_qtree.jsig_cxts[0][i] = sig0_zc_cxts[i];
    cxt_qtree.jsig_tabs[0][i] = sig0_zc_lut[i];
    total += sig0_zc_cxts[i];
  }
  cxt_qtree.sig_cxts[0] = total - cxt_qtree.sig_offsets[0];

  for( lev = 1; lev < 3; lev++ ) {
    cxt_qtree.sig_offsets[lev] = total;
    for( i = 0; i < 4; i++ ) {
      cxt_qtree.jsig_offsets[lev][i] = total;
      cxt_qtree.jsig_cxts[lev][i] = sig_zc_cxts[i];
      cxt_qtree.jsig_tabs[lev][i] = sig_zc_lut[i];
      total += sig_zc_cxts[i];
    }
    cxt_qtree.sig_cxts[lev] = total - cxt_qtree.sig_offsets[lev];
  }



  for( lev = 3; lev < depth; lev++ ) {  // 3 ~ depth-1
    cxt_qtree.sig_offsets[lev] = cxt_qtree.sig_offsets[2];
    for( i = 0; i < 4; i++ ) {
      cxt_qtree.jsig_offsets[lev][i] = cxt_qtree.jsig_offsets[2][i];
      cxt_qtree.jsig_cxts[lev][i] = 0;
      cxt_qtree.jsig_tabs[lev][i] = cxt_qtree.jsig_tabs[2][i];
    }
    cxt_qtree.sig_cxts[lev] = 0;
  }
  cxt_qtree.LIP_offset = total;
  cxt_qtree.LIP_cxts = LIP_zc_cxts;
  cxt_qtree.LIP_tab = LIP_zc_lut;
  total += cxt_qtree.LIP_cxts;

  for( lev = 1; lev < 3; lev++ ) {
    cxt_qtree.node_offsets[lev] = total;
    cxt_qtree.node_cxts[lev] = node_zc_cxts;
    cxt_qtree.node_tabs[lev] = node_zc_lut;
    total += cxt_qtree.node_cxts[lev];
  }

  for( lev = 3; lev < depth; lev++ ) {  // 3 ~ depth-1
    cxt_qtree.node_offsets[lev] = cxt_qtree.node_offsets[2];
    cxt_qtree.node_cxts[lev] = 0;
    cxt_qtree.node_tabs[lev] = cxt_qtree.node_tabs[2];
  }

  cxt_qtree.LSP_offset = total;

#ifdef LSP_BIT_IDX


#ifdef LSP_LUT

#ifndef MERGE_ALL_LSP_PLANES
  cxt_qtree.LSP_cxts = 79;
#else
  cxt_qtree.LSP_cxts = LSP_cxts;
#endif //MERGE_ALL_LSP_PLANES

  cxt_qtree.LSP_tab = LSP_lut;
#else //LSP_LUT
  cxt_qtree.LSP_cxts = 15;
#endif //LSP_LUT


#else //LSP_BIT_IDX

#ifdef LSP_LUT
  cxt_qtree.LSP_cxts = LSP_cxts + 1;
  cxt_qtree.LSP_tab = LSP_lut;
#else
  cxt_qtree.LSP_cxts = MAG_CONTEXTS;
#endif

#endif //LSP_BIT_IDX

  total += cxt_qtree.LSP_cxts;

  cxt_qtree.sign_offset = total;

  cxt_qtree.sign_cxts = sign_cxts;
  cxt_qtree.sign_tab = sign_lut;

  total += cxt_qtree.sign_cxts;

  cxt_qtree.total_cxts = total;

}

//-----------------------------------------------------------------------------

//   reset_cxt_models()

//-----------------------------------------------------------------------------
void
SubbandCodec::reset_cxt_models( void )
{

  for( int i = cxt_qtree.total_cxts - 1; i >= 0;
       cxt_qtree.cxt_models[i--].taub_scale(  ) );
}

//-----------------------------------------------------------------------------

//   reset_node_cxts(void)

//-----------------------------------------------------------------------------
void
SubbandCodec::update_node_cxts( void )
{
  std_int *LIP_cur, *LIP_end, *LIS_cur, *LIS_end, coord_cur;
  PEL_CXT_TYPE *base_cxt_sp, **base_cxt;
  NODE_CXT_TYPE *cxt_node_sp, **cxt_node;

  LIP_cur = node_list.LSP + node_list.LIP_sz;
  LIP_end = node_list.LIP_end;
  base_cxt = cxt_qtree.base_cxt;


  while( LIP_cur > LIP_end ) {
    coord_cur = *( --LIP_cur );
    base_cxt_sp = base_cxt[coord_cur >> 16] + ( coord_cur & 0xFFFF );

    if( *base_cxt_sp & TC_SIG )
      *base_cxt_sp |= TC2_SIG;

    if( *base_cxt_sp & CL_SIG )
      *base_cxt_sp |= CL2_SIG;

    if( *base_cxt_sp & CR_SIG )
      *base_cxt_sp |= CR2_SIG;

    if( *base_cxt_sp & BC_SIG )
      *base_cxt_sp |= BC2_SIG;
  }

  int depth = qtree.depth;
  for( int lev = 1; lev < depth; lev++ ) {
    LIS_cur = node_list.LIS[lev] - 1;
    LIS_end = node_list.LIS_end[lev];
    cxt_node = cxt_qtree.cxt_nodes[lev];
    while( LIS_cur < LIS_end ) {
      coord_cur = *( ++LIS_cur );
      cxt_node_sp = cxt_node[coord_cur >> 16] + ( coord_cur & 0xFFFF );

      if( *cxt_node_sp & TC_SIG )
        *cxt_node_sp |= TC2_SIG;

      if( *cxt_node_sp & CL_SIG )
        *cxt_node_sp |= CL2_SIG;

      if( *cxt_node_sp & CR_SIG )
        *cxt_node_sp |= CR2_SIG;

      if( *cxt_node_sp & BC_SIG )
        *cxt_node_sp |= BC2_SIG;
    }
  }

  std_int *LSP_cur = node_list.LSP - 1;
  std_int *LSP_end = node_list.LSP_end;

  while( LSP_cur < LSP_end ) {
    coord_cur = *( ++LSP_cur );
    base_cxt_sp = base_cxt[coord_cur >> 16] + ( coord_cur & 0xFFFF );

    if( *base_cxt_sp & TC_SIG )
      *base_cxt_sp |= TC2_SIG;

    if( *base_cxt_sp & CL_SIG )
      *base_cxt_sp |= CL2_SIG;

    if( *base_cxt_sp & CR_SIG )
      *base_cxt_sp |= CR2_SIG;

    if( *base_cxt_sp & BC_SIG )
      *base_cxt_sp |= BC2_SIG;
  }

}

//-----------------------------------------------------------------------------

//       setup_luts(void)

//-----------------------------------------------------------------------------

void
SubbandCodec::setup_luts( void )
{
  initialize_zc_luts(  );       /* in setup_cxt_tables.C */
  initialize_sc_lut(  );
}

//-----------------------------------------------------------------------------

//  Member function of class SubbandTreeCodec()

//-----------------------------------------------------------------------------

SubbandTreeCodec::SubbandTreeCodec( SUBBAND_TREE_TYPE * subs )
:fpout( stdout )
{
  //assert(subband_tree = subs);
  myassert( subband_tree = subs );
  total_byte_budget = MAX_STD_INT;
}

//-----------------------------------------------------------------------------

//  ~SubbandTreeCodec()

//-----------------------------------------------------------------------------
SubbandTreeCodec::~SubbandTreeCodec( void )
{
}

void
SubbandTreeCodec::initialize( SUBBAND_TREE_TYPE * subs )
{
  //assert(subband_tree = subs);
  myassert( subband_tree = subs );
  total_byte_budget = MAX_STD_INT;

  fpout = stdout;
}

void
SubbandTreeCodec::reset_tree_codec( SUBBAND_TREE_TYPE * subs )
{
  //assert(subband_tree = subs);
  myassert( subband_tree = subs );
}
