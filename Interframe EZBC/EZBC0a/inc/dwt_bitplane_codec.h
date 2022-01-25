/* ========================================================================= */
/* Description: definition of  classes SubbandCodec and SubbandTreeCodec     */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */
#ifndef _DWT_CODEC_H
#define _DWT_CODEC_H

#include <string.h>
#include "ar_code.h"
#include "subband.h"
#include "context_status_consts.h"

//-----------------------------------------------------------------------------
typedef HistoBiModel MODEL_TYPE;        //type of context prob models
typedef HistoBi STAT_MODEL_TYPE;        //type of context prob models

typedef SubSetLayer SUBBAND_TREE_TYPE;  //class SubSetLayer: public SubSet in Utils/subband.h


typedef std_short PEL_CXT_TYPE;
typedef std_short NODE_CXT_TYPE;
typedef std_byte SIGN_CXT_TYPE;


struct Image_Coord_Sht
{
  std_short r, c;
};

struct Quad_Tree
{
  int depth;
  int qtree_sz;
  Image_Coord_Sht *dims;
  SUB_COEFF_TYPE ***nodes;
};

struct Node_Entry
{
  int level;
  std_int node;
};

struct List
{
  std_int *LIP;
  std_int *LIP_end;
  std_int *LIP_prev_mark;
  std_int *LIP_old_end;
  std_int LIP_sz;
  std_int *LSP;
  std_int *LSP_end;
  std_int LSP_mark[sizeof( SUB_COEFF_TYPE ) * 8];
  int LSP_plane;
  std_int **LIS;
  std_int **LIS_end;
  std_int **LIS_prev_mark;
  std_int **LIS_old_end;
  Node_Entry *LIS_stack;
  Node_Entry *LIS_stack_top;

  std_int **LSP_bit_idx_marks;
  std_int ***LSP_ids;
};

struct PassStats
{
  std_int cumulative_bytes;
  double delta_mse;
  float rd_slope;
  std_int pass_bytes;
};

struct BitPlaneStats
{
  PassStats *passes;
};

struct CodingStats
{
  BitPlaneStats *bitplanes;
};

struct Cxt_States
{
  int cxt_qtree_sz;
  NODE_CXT_TYPE ***cxt_nodes;
  PEL_CXT_TYPE **base_cxt;
  SIGN_CXT_TYPE **base_sign_cxt;
  int *sig_cxts;
  int *sig_offsets;
  std_byte **sig_tabs;
  int *node_cxts;
  int *node_offsets;
  std_byte **node_tabs;
  int LIP_cxts;
  int LIP_offset;
  std_byte *LIP_tab;
  int LSP_cxts;
  int LSP_offset;
  std_byte *LSP_tab;
  int sign_cxts;
  int sign_offset;
  SUB_COEFF_TYPE *sign_tab;     //SUB_COEFF_TYPE used,
  //see setup_cxt_tables() for details
  int total_cxts;
  MODEL_TYPE *cxt_models;

  int **jsig_cxts;
  int **jsig_offsets;
  std_byte ***jsig_tabs;
};

class SubbandCodec
{
protected:
  int bit_idx, max_bit_idx, min_bit_idx, no_bitplanes, no_passes, cur_pass,
    cur_lev;
  SUB_COEFF_TYPE **base_coeff, bit_idx_mask;
  SUBBAND_TYPE *band;
  List node_list;
  Quad_Tree qtree;
  Cxt_States cxt_qtree, *child_cxt_qtree, *grand_child_cxt_qtree;
#ifdef GET_PARENT_MODELS
  Cxt_States *par_cxt_qtree;
#endif
  CodingStats coding_stats;
  FILE *fpout;

  void initialize_node_list( void );
  void clear_node_list( void );
  void create_coding_stats( void );
  void delete_coding_stats( void );
  void set_bit_idx_mask( void )
  {
    bit_idx_mask = ( SUB_COEFF_TYPE ) ( 0x1 << bit_idx );
  }
  void initialize_cxt_models( void );   //cxt related
  void create_cxt_models_AC0( void );
  void create_cxt_models_EBCOT( void );
  void create_cxt_models( void );
  void reset_cxt_models( void );
  void update_node_cxts( void );


//=============================================================================

//static variables, cxt lookup tables

//===========================================================================
//  static long BYTE_BUDGET;
  static std_byte *sm0_sp_main_node_lut, *sm0_sp_diag_node_lut;

  static std_byte *mid_inter_and_sp_main_node_lut;
  static std_byte *mid_inter_and_sp_diag_node_lut;
  static int sm0_sp_main_node_cxts, sm0_sp_diag_node_cxts;

  static int mid_inter_and_sp_main_node_cxts;
  static int mid_inter_and_sp_diag_node_cxts;

  static std_byte *mid_main_jsig_00_lut, *mid_diag_jsig_00_lut;
  static std_byte *mid_main_jsig_01_lut, *mid_diag_jsig_01_lut;
  static std_byte *mid_main_jsig_10_lut, *mid_diag_jsig_10_lut;
  static std_byte *mid_main_jsig_11_lut, *mid_diag_jsig_11_lut;

  static int mid_main_jsig_00_cxts, mid_diag_jsig_00_cxts;
  static int mid_main_jsig_01_cxts, mid_diag_jsig_01_cxts;
  static int mid_main_jsig_10_cxts, mid_diag_jsig_10_cxts;
  static int mid_main_jsig_11_cxts, mid_diag_jsig_11_cxts;

  static std_byte *mid_main_jsig0_00_lut, *mid_diag_jsig0_00_lut;
  static std_byte *mid_main_jsig0_01_lut, *mid_diag_jsig0_01_lut;
  static std_byte *mid_main_jsig0_10_lut, *mid_diag_jsig0_10_lut;
  static std_byte *mid_main_jsig0_11_lut, *mid_diag_jsig0_11_lut;

  static int mid_main_jsig0_00_cxts, mid_diag_jsig0_00_cxts;
  static int mid_main_jsig0_01_cxts, mid_diag_jsig0_01_cxts;
  static int mid_main_jsig0_10_cxts, mid_diag_jsig0_10_cxts;
  static int mid_main_jsig0_11_cxts, mid_diag_jsig0_11_cxts;

  static std_byte *mid_main_LSP_lut, *mid_diag_LSP_lut;
  static int mid_main_LSP_cxts, mid_diag_LSP_cxts;

#ifdef TWO_SIGNIF_BITS

  static std_byte *sm3_2_signif_bits_main_LSP_lut;
  static std_byte *sm3_2_signif_bits_diag_LSP_lut;
  static int sm3_2_signif_bits_main_LSP_cxts, sm3_2_signif_bits_diag_LSP_cxts;
#endif


  static SUB_COEFF_TYPE *mid_main_sign_lut, *mid_diag_sign_lut;
  static int mid_main_sign_cxts, mid_diag_sign_cxts;

public:
  SubbandCodec( void )
  {
    fpout = stdout;
  }
  SubbandCodec( SUBBAND_TYPE * subband ) {
    initialize( subband );
  }
  ~SubbandCodec( void );
  int initialize( SUBBAND_TYPE * subband );
  void reset_band( SUBBAND_TYPE * subband )
  {
    band = subband;
  }
  //int get_no_bitplanes(void){return no_bitplanes;}
  int get_max_bit_idx(  )
  {
    return max_bit_idx;
  }
  void set_fpout( FILE * fp )
  {
    fpout = fp;
  }
  void reset_band_codec( SUBBAND_TYPE * subband );
  int get_qtree_depth( void )
  {
    return qtree.depth;
  }


  static void setup_luts( void );
  static void delete_cxt_tables( void );
  /*static void set_byte_budget(long assigned_bytes){
     BYTE_BUDGET = assigned_bytes;} */
  static void initialize_zc_luts( void );
  static void initialize_sc_lut( void );

  static void initialize_sm0_sp_node_luts( void );
  static void initialize_inter_and_sp_node_luts( void );

  static void initialize_mid_LSP_luts( void );

  static void initialize_mid_sign_luts( void );

#ifdef TWO_SIGNIF_BITS
  static void initialize_sm3_2_signif_bits_LSP_luts( void );
#endif

  static void initialize_mid_jsig_luts( void );

  static void initialize_mid_jsig0_00_luts( void );
  static void initialize_mid_jsig0_01_luts( void );
  static void initialize_mid_jsig0_10_luts( void );
  static void initialize_mid_jsig0_11_luts( void );
};

class SubbandTreeCodec
{
protected:
  SUBBAND_TREE_TYPE * subband_tree;
  int tree_msb, tree_lsb;
  long total_byte_budget;
  int subblock_dim;
  FILE *fpout;
public:
    SubbandTreeCodec( void )
  {
    subband_tree = NULL;
    fpout = stdout;
  }
  SubbandTreeCodec( SUBBAND_TREE_TYPE * subs );
  ~SubbandTreeCodec( void );
  int get_tree_msb(  )
  {
    return tree_msb;
  }
  int get_tree_lsb(  )
  {
    return tree_lsb;
  }
  void set_total_byte_budget( long total_byte )
  {
    total_byte_budget = total_byte;
  }
  void set_fpout( FILE * fp )
  {
    fpout = fp;
  }
  void initialize( SUBBAND_TREE_TYPE * subs );
  void reset_tree_codec( SUBBAND_TREE_TYPE * subs );
};

#define UPDATE_CXT(cp, row_gap){ (cp)[-1] |= CR_SIG;\
  (cp)[1] |= CL_SIG;   (cp)[-row_gap] |= BC_SIG;     (cp)[row_gap] |= TC_SIG;\
  (cp)[-row_gap-1] |= BR_SIG;  (cp)[-row_gap+1] |= BL_SIG;\
  (cp)[row_gap-1]  |= TR_SIG;  (cp)[row_gap+1]  |= TL_SIG;}



#ifdef UPDATE_SIGN_DIAG_CXT

#define UPDATE_SIGN_CXT(cp, row_gap, sign_bit){ \
  if(sign_bit){ (cp)[-1] |= H_NVE_MASK; (cp)[1] |= H_NVE_MASK;\
      (cp)[-row_gap] |= V_NVE_MASK; (cp)[row_gap] |= V_NVE_MASK;\
      (cp)[-row_gap-1] |= NW_NVE_MASK; (cp)[-row_gap+1] |= NE_NVE_MASK;\
      (cp)[row_gap+1] |= NW_NVE_MASK; (cp)[row_gap-1] |= NE_NVE_MASK;}\
  else{ (cp)[-1] |= H_PVE_MASK; (cp)[1] |= H_PVE_MASK;\
      (cp)[-row_gap] |= V_PVE_MASK; (cp)[row_gap] |= V_PVE_MASK;\
      (cp)[row_gap+1] |= NW_PVE_MASK; (cp)[row_gap-1] |= NE_PVE_MASK;}}

#else


#define UPDATE_SIGN_CXT(cp, row_gap, sign_bit){ \
  if(sign_bit){ (cp)[-1] |= H_NVE_MASK; (cp)[1] |= H_NVE_MASK;\
      (cp)[-row_gap] |= V_NVE_MASK; (cp)[row_gap] |= V_NVE_MASK;}\
  else{ (cp)[-1] |= H_PVE_MASK; (cp)[1] |= H_PVE_MASK;\
      (cp)[-row_gap] |= V_PVE_MASK; (cp)[row_gap] |= V_PVE_MASK;}}

#endif


#define UPDATE_1_CD_NODE(cd_node){*cd_node |= PA_SIG;}
#define UPDATE_4_CD_NODES(cd_node, cd_node_row_gap){\
  cd_node[0] |= PA_SIG;                 cd_node[1] |= PA_SIG; \
  cd_node[cd_node_row_gap] |= PA_SIG; cd_node[cd_node_row_gap+1] |= PA_SIG;}

#endif // file dwt_bitplane_codec.h
