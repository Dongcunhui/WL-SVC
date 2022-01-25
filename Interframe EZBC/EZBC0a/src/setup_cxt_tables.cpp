/* ========================================================================= */
/* Description: luts for context-based AR. coding*/
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */


#include "dwt_bitplane_codec.h"


std_byte *
  SubbandCodec::sm0_sp_main_node_lut =
  NULL;
std_byte *
  SubbandCodec::sm0_sp_diag_node_lut =
  NULL;

std_byte *
  SubbandCodec::mid_inter_and_sp_main_node_lut =
  NULL;
std_byte *
  SubbandCodec::mid_inter_and_sp_diag_node_lut =
  NULL;

std_byte *
  SubbandCodec::mid_main_LSP_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_LSP_lut =
  NULL;

SUB_COEFF_TYPE *
  SubbandCodec::mid_main_sign_lut =
  NULL;
SUB_COEFF_TYPE *
  SubbandCodec::mid_diag_sign_lut =
  NULL;

int
  SubbandCodec::sm0_sp_main_node_cxts =
  4;
int
  SubbandCodec::sm0_sp_diag_node_cxts =
  3;

int
  SubbandCodec::mid_inter_and_sp_main_node_cxts =
  10;
int
  SubbandCodec::mid_inter_and_sp_diag_node_cxts =
  9;

int
  SubbandCodec::mid_main_LSP_cxts =
  5;
int
  SubbandCodec::mid_diag_LSP_cxts =
  4;

int
  SubbandCodec::mid_main_sign_cxts =
  9;
int
  SubbandCodec::mid_diag_sign_cxts =
  9;

#ifdef TWO_SIGNIF_BITS

std_byte *
  SubbandCodec::sm3_2_signif_bits_main_LSP_lut =
  NULL;
std_byte *
  SubbandCodec::sm3_2_signif_bits_diag_LSP_lut =
  NULL;

int
  SubbandCodec::sm3_2_signif_bits_main_LSP_cxts =
  3;
int
  SubbandCodec::sm3_2_signif_bits_diag_LSP_cxts =
  3;
#endif

std_byte *
  SubbandCodec::mid_main_jsig_00_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig_00_lut =
  NULL;
std_byte *
  SubbandCodec::mid_main_jsig0_00_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig0_00_lut =
  NULL;

std_byte *
  SubbandCodec::mid_main_jsig_01_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig_01_lut =
  NULL;
std_byte *
  SubbandCodec::mid_main_jsig0_01_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig0_01_lut =
  NULL;

std_byte *
  SubbandCodec::mid_main_jsig_10_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig_10_lut =
  NULL;
std_byte *
  SubbandCodec::mid_main_jsig0_10_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig0_10_lut =
  NULL;

std_byte *
  SubbandCodec::mid_main_jsig_11_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig_11_lut =
  NULL;
std_byte *
  SubbandCodec::mid_main_jsig0_11_lut =
  NULL;
std_byte *
  SubbandCodec::mid_diag_jsig0_11_lut =
  NULL;

int
  SubbandCodec::mid_main_jsig_00_cxts =
  10;
int
  SubbandCodec::mid_diag_jsig_00_cxts =
  10;
int
  SubbandCodec::mid_main_jsig0_00_cxts =
  10;
int
  SubbandCodec::mid_diag_jsig0_00_cxts =
  10;

int
  SubbandCodec::mid_main_jsig_01_cxts =
  10;
int
  SubbandCodec::mid_diag_jsig_01_cxts =
  10;
int
  SubbandCodec::mid_main_jsig0_01_cxts =
  10;
int
  SubbandCodec::mid_diag_jsig0_01_cxts =
  10;

int
  SubbandCodec::mid_main_jsig_10_cxts =
  30;
int
  SubbandCodec::mid_diag_jsig_10_cxts =
  30;
int
  SubbandCodec::mid_main_jsig0_10_cxts =
  30;
int
  SubbandCodec::mid_diag_jsig0_10_cxts =
  30;

int
  SubbandCodec::mid_main_jsig_11_cxts =
  11;
int
  SubbandCodec::mid_diag_jsig_11_cxts =
  11;
int
  SubbandCodec::mid_main_jsig0_11_cxts =
  11;
int
  SubbandCodec::mid_diag_jsig0_11_cxts =
  11;

//===========================================================================

void
SubbandCodec::initialize_zc_luts( void )
{
  initialize_sm0_sp_node_luts(  );

  initialize_inter_and_sp_node_luts(  );

  initialize_mid_LSP_luts(  );

  initialize_sm3_2_signif_bits_LSP_luts(  );

  initialize_mid_jsig_luts(  );

  initialize_mid_jsig0_00_luts(  );
  initialize_mid_jsig0_01_luts(  );
  initialize_mid_jsig0_10_luts(  );
  initialize_mid_jsig0_11_luts(  );

}

//-----------------------------------------------------------------------------

//  initialize_sc_lut()

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_sc_lut( void )
{
  initialize_mid_sign_luts(  );
}

//-----------------------------------------------------------------------------

//STATIC delete_cxt_tables()

//-----------------------------------------------------------------------------
void
SubbandCodec::delete_cxt_tables(  )
{
  DELETE_VECTOR( sm0_sp_main_node_lut );
  DELETE_VECTOR( sm0_sp_diag_node_lut );

  DELETE_VECTOR( mid_inter_and_sp_main_node_lut );
  DELETE_VECTOR( mid_inter_and_sp_diag_node_lut );

  DELETE_VECTOR( mid_main_sign_lut );
  DELETE_VECTOR( mid_diag_sign_lut );

  DELETE_VECTOR( mid_main_LSP_lut );
  DELETE_VECTOR( mid_diag_LSP_lut );

#ifdef TWO_SIGNIF_BITS
  DELETE_VECTOR( sm3_2_signif_bits_main_LSP_lut );
  DELETE_VECTOR( sm3_2_signif_bits_diag_LSP_lut );
#endif


  DELETE_VECTOR( mid_main_jsig_00_lut );
  DELETE_VECTOR( mid_diag_jsig_00_lut );
  DELETE_VECTOR( mid_main_jsig0_00_lut );
  DELETE_VECTOR( mid_diag_jsig0_00_lut );

  DELETE_VECTOR( mid_main_jsig_01_lut );
  DELETE_VECTOR( mid_diag_jsig_01_lut );
  DELETE_VECTOR( mid_main_jsig0_01_lut );
  DELETE_VECTOR( mid_diag_jsig0_01_lut );

  DELETE_VECTOR( mid_main_jsig_10_lut );
  DELETE_VECTOR( mid_diag_jsig_10_lut );
  DELETE_VECTOR( mid_main_jsig0_10_lut );
  DELETE_VECTOR( mid_diag_jsig0_10_lut );

  DELETE_VECTOR( mid_main_jsig_11_lut );
  DELETE_VECTOR( mid_diag_jsig_11_lut );
  DELETE_VECTOR( mid_main_jsig0_11_lut );
  DELETE_VECTOR( mid_diag_jsig0_11_lut );
}



//-----------------------------------------------------------------------------

//  initialize_sm0_sp_node_luts()

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_sm0_sp_node_luts( void )
{
  std_short idx, ctxt, v1, v2, v3, diag_ctxt;

  NEW_VECTOR( sm0_sp_main_node_lut, ZC_MASK + 1, std_byte,
              "sm0_sp_main_node_lut" );

  NEW_VECTOR( sm0_sp_diag_node_lut, ZC_MASK + 1, std_byte,
              "sm0_sp_diag_node_lut" );

  for( idx = 0; idx <= ZC_MASK; idx++ ) {
    /* First, form the context map for the horizontal and vertical bands.
       These both use the same context map, because the horizontally
       high-pass band is physically transposed before encoding. */

    v1 = ( ( idx >> CL_POS ) & 1 ) + ( ( idx >> CR_POS ) & 1 );
    v2 = ( ( idx >> TC_POS ) & 1 ) + ( ( idx >> BC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 ) + ( ( idx >> BR_POS ) & 1 );


    if( ( v1 == 0 ) && ( v2 == 0 ) && ( v3 == 1 ) ) {
      ctxt = 0;
      diag_ctxt = 0;
    } else if( ( v1 == 0 ) && ( v2 == 0 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else if( ( v1 == 0 ) && ( v2 == 1 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else if( ( v1 == 1 ) && ( v2 == 0 ) ) {
      ctxt = 2;
      diag_ctxt = 1;
    } else if( ( ( v1 == 0 ) && ( v2 == 2 ) )
               || ( ( v1 == 2 ) && ( v2 == 0 ) ) ) {
      ctxt = 3;
      diag_ctxt = 2;
    } else if( ( v1 + v2 ) > 2 ) {
      ctxt = 3;
      diag_ctxt = 2;
    } else {                    //(v1 == 1) && (v2 == 1)
      ctxt = 3;
      diag_ctxt = 2;
    }

    assert( ctxt < sm0_sp_main_node_cxts );
    sm0_sp_main_node_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < sm0_sp_diag_node_cxts );
    sm0_sp_diag_node_lut[idx] = ( std_byte ) diag_ctxt;
  }
}

//-----------------------------------------------------------------------------

// initialize_inter_and_sp_node_luts()

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_inter_and_sp_node_luts( void )
{

  std_short idx, ctxt, v1, v2, v3, diag_ctxt, v0;

  NEW_VECTOR( mid_inter_and_sp_main_node_lut, ZC_MASK + 1, std_byte,
              "mid_inter_and_sp_main_node_lut" );
  NEW_VECTOR( mid_inter_and_sp_diag_node_lut, ZC_MASK + 1, std_byte,
              "mid_inter_and_sp_diag_node_lut" );


  for( idx = 0; idx <= ZC_MASK; idx++ ) {

    v0 = ( ( idx >> PA_POS ) & 1 );
    v1 = ( ( idx >> CL_POS ) & 1 ) + ( ( idx >> CR_POS ) & 1 );
    v2 = ( ( idx >> TC_POS ) & 1 ) + ( ( idx >> BC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 ) + ( ( idx >> BR_POS ) & 1 );

    if( v0 ) {
      if( ( v1 == 0 ) && ( v2 == 0 ) && ( v3 == 1 ) ) {
        ctxt = 0;
        diag_ctxt = 0;
      } else if( ( v1 == 0 ) && ( v2 == 0 ) ) {
        ctxt = 1;
        diag_ctxt = 1;
      } else if( ( v1 == 0 ) && ( v2 == 1 ) ) {
        ctxt = 1;
        diag_ctxt = ( v3 < 3 ) ? 1 : 2;
      } else if( ( v1 == 1 ) && ( v2 == 0 ) ) {
        ctxt = 2;
        diag_ctxt = ( v3 < 3 ) ? 1 : 2;
      } else if( ( ( v1 == 0 ) && ( v2 == 2 ) )
                 || ( ( v1 == 2 ) && ( v2 == 0 ) ) ) {
        ctxt = 3;
        diag_ctxt = 2;
      } else if( ( v1 + v2 ) > 2 ) {
        ctxt = 4;
        diag_ctxt = 3;
      } else {                  //(v1 == 1) && (v2 == 1)
        ctxt = 3;
        diag_ctxt = ( v3 < 2 ) ? 1 : 2;
      }
    } else {
      if( ( v1 == 0 ) && ( v2 == 0 ) ) {
        ctxt = 5;
        diag_ctxt = 4;
      } else if( ( v1 + v2 == 1 ) ) {
        ctxt = ( v3 < 3 ) ? 5 : 6;
        diag_ctxt = ( v3 < 3 ) ? 4 : 5;
      } else if( ( ( v1 == 0 ) && ( v2 == 2 ) )
                 || ( ( v1 == 2 ) && ( v2 == 0 ) ) ) {
        ctxt = 6;
        diag_ctxt = 5;
      } else if( ( v1 + v2 ) > 2 ) {
        ctxt = 7;
        diag_ctxt = 6;
      } else {                  //(v1 == 1) && (v2 == 1)
        ctxt = 6;
        diag_ctxt = 5;
      }
    }

    assert( ctxt < mid_inter_and_sp_main_node_cxts );
    mid_inter_and_sp_main_node_lut[idx] = ( std_byte ) ctxt;
    assert( diag_ctxt < mid_inter_and_sp_diag_node_cxts );
    mid_inter_and_sp_diag_node_lut[idx] = ( std_byte ) diag_ctxt;


  }
}


//-----------------------------------------------------------------------------

//  initialize_mid_sign_luts(void)

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_mid_sign_luts( void )
{

  ifc_int idx, vpos, vneg, hpos, hneg, ctxt, predict, v1, v2;
  ifc_int diag_ctxt, diag_predict;
  ifc_int nwpos, nwneg, nepos, neneg, v3, v4;

  NEW_VECTOR( mid_main_sign_lut, SIGN_CXT_MASK + 1, SUB_COEFF_TYPE,
              "mid_main_sign_lut" );
  NEW_VECTOR( mid_diag_sign_lut, SIGN_CXT_MASK + 1, SUB_COEFF_TYPE,
              "mid_diag_sign_lut" );

  for( idx = 0; idx <= SIGN_CXT_MASK; idx++ ) {

    vpos = ( idx >> ( V_PVE_BIT_POS - SIGN_BIT_POS ) ) & 1;
    vneg = ( idx >> ( V_NVE_BIT_POS - SIGN_BIT_POS ) ) & 1;
    hpos = ( idx >> ( H_PVE_BIT_POS - SIGN_BIT_POS ) ) & 1;
    hneg = ( idx >> ( H_NVE_BIT_POS - SIGN_BIT_POS ) ) & 1;

    v1 = hpos - hneg;
    v2 = vpos - vneg;

    if( v2 > 0 ) {
      diag_predict = predict = MIN_IFC_INT;
    } else {
      diag_predict = predict = 0;
      if( v2 < 0 ) {
        v1 = -v1;
        v2 = -v2;
      }
    }

    if( v2 == 0 ) {
      if( v1 < 0 ) {
        predict = MIN_IFC_INT;
        ctxt = 1, diag_ctxt = 1;
      } else if( v1 > 0 ) {
        diag_predict = MIN_IFC_INT;
        ctxt = 1, diag_ctxt = 1;
      } else {
        ctxt = 0;
        diag_ctxt = 0;
      }
    } else {
      ctxt = 3 + v1;
      diag_ctxt = 3 + v1;
    }

    if( ctxt == 0 ) {
      nwpos = ( idx >> ( NW_PVE_BIT_POS - SIGN_BIT_POS ) ) & 1;
      nwneg = ( idx >> ( NW_NVE_BIT_POS - SIGN_BIT_POS ) ) & 1;

      nepos = ( idx >> ( NE_PVE_BIT_POS - SIGN_BIT_POS ) ) & 1;
      neneg = ( idx >> ( NE_NVE_BIT_POS - SIGN_BIT_POS ) ) & 1;

      v3 = nwpos - nwneg;
      v4 = nepos - neneg;

      if( v3 < 0 ) {
        v3 = -v3;
        v4 = -v4;
        diag_predict = predict = MIN_IFC_INT;
      }

      if( v3 == 0 ) {
        if( v4 < 0 ) {
          diag_predict = MIN_IFC_INT;
          ctxt = 5, diag_ctxt = 5;
          //predict = MIN_IFC_INT;
        } else if( v4 > 0 ) {
          predict = MIN_IFC_INT;
          ctxt = 5, diag_ctxt = 5;
        } else {
          ctxt = 0;
          diag_ctxt = 0;
        }
      } else {
        ctxt = 5;
        diag_ctxt = 5;
      }
    }
    assert( ( ctxt >= 0 ) && ( ctxt < mid_main_sign_cxts ) );
    mid_main_sign_lut[idx] = ctxt | predict;

    assert( ( diag_ctxt >= 0 ) && ( diag_ctxt < mid_diag_sign_cxts ) );
    mid_diag_sign_lut[idx] = diag_ctxt | diag_predict;
  }

}


//-----------------------------------------------------------------------------

//  initialize_mid_LSP_lut()

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_mid_LSP_luts( void )
{
  std_short idx, ctxt, v1, v2, v3, diag_ctxt;

  NEW_VECTOR( mid_main_LSP_lut, ZC_MASK + 1, std_byte, "mid_main_LSP_lut" );

  NEW_VECTOR( mid_diag_LSP_lut, ZC_MASK + 1, std_byte, "mid_main_LSP_lut" );

  for( idx = 0; idx <= ZC_MASK; idx++ ) {
    /* First, form the context map for the horizontal and vertical bands.
       These both use the same context map, because the horizontally
       high-pass band is physically transposed before encoding. */

    v1 = ( ( idx >> CL_POS ) & 1 ) + ( ( idx >> CR_POS ) & 1 );
    v2 = ( ( idx >> TC_POS ) & 1 ) + ( ( idx >> BC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 ) + ( ( idx >> BR_POS ) & 1 );

    if( ( v1 == 0 ) && ( v2 == 0 ) ) {
      ctxt = 0;
      diag_ctxt = 0;
    } else if( ( v1 == 0 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else if( ( v2 == 0 ) ) {
      ctxt = 2;
      diag_ctxt = 1;
    } else if( ( v1 == 2 ) && ( v2 == 2 ) ) {
      ctxt = 4;
      diag_ctxt = 3;
    } else {
      ctxt = 3;
      diag_ctxt = 2;
    }

    assert( ctxt < mid_main_LSP_cxts );
    mid_main_LSP_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_LSP_cxts );
    mid_diag_LSP_lut[idx] = ( std_byte ) diag_ctxt;

  }
}

#ifdef TWO_SIGNIF_BITS

//-----------------------------------------------------------------------------

//  initialize_sm3_2_signif_bits_LSP_luts(void)

//-----------------------------------------------------------------------------
void
SubbandCodec::initialize_sm3_2_signif_bits_LSP_luts( void )
{
  std_short idx, ctxt, v1, v2, v3, vv1, vv2, diag_ctxt;

  NEW_VECTOR( sm3_2_signif_bits_main_LSP_lut, ZC_MASK + 1, std_byte,
              "sm3_2_signif_bits_main_LSP_lut" );

  NEW_VECTOR( sm3_2_signif_bits_diag_LSP_lut, ZC_MASK + 1, std_byte,
              "sm3_2_signif_bits_diag_LSP_lut" );

  for( idx = 0; idx <= ZC_MASK; idx++ ) {
    /* First, form the context map for the horizontal and vertical bands.
       These both use the same context map, because the horizontally
       high-pass band is physically transposed before encoding. */

    vv1 = ( ( idx >> CL2_POS ) & 1 ) + ( ( idx >> CR2_POS ) & 1 );
    vv2 = ( ( idx >> TC2_POS ) & 1 ) + ( ( idx >> BC2_POS ) & 1 );

    v1 = ( ( idx >> CL_POS ) & 1 ) + ( ( idx >> CR_POS ) & 1 );
    v2 = ( ( idx >> TC_POS ) & 1 ) + ( ( idx >> BC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 ) + ( ( idx >> BR_POS ) & 1 );


    if( ( v1 == 0 ) && ( v2 == 0 ) ) {
      ctxt = 0;
      diag_ctxt = 0;
    } else if( ( vv1 == 0 ) && ( vv2 == 0 ) && ( v3 < 2 ) ) {
      ctxt = 0;
      diag_ctxt = 0;
    } else if( ( vv1 == 0 ) && ( vv2 == 0 ) ) {
      ctxt = 0;
      diag_ctxt = 0;
    } else if( ( vv1 == 0 ) ) {
      ctxt = 2;
      diag_ctxt = 1;
    } else if( ( vv2 == 0 ) ) {
      ctxt = 2;
      diag_ctxt = 1;
    } else if( ( vv1 == 2 ) && ( vv2 == 2 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else if( ( vv1 == 2 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else if( ( vv2 == 2 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else {
      ctxt = 2;
      diag_ctxt = 1;
    }


    assert( ctxt < sm3_2_signif_bits_main_LSP_cxts );
    sm3_2_signif_bits_main_LSP_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < sm3_2_signif_bits_diag_LSP_cxts );
    sm3_2_signif_bits_diag_LSP_lut[idx] = ( std_byte ) diag_ctxt;
  }
}
#endif

void
SubbandCodec::initialize_mid_jsig_luts( void )
{
  std_short idx, ctxt, diag_ctxt, v1, v2, v3, p, l, r, t, bc, tr;       //v0, 

  NEW_VECTOR( mid_main_jsig_00_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig_00_lut" );
  NEW_VECTOR( mid_diag_jsig_00_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig_00_lut" );

  NEW_VECTOR( mid_main_jsig_01_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig_01_lut" );
  NEW_VECTOR( mid_diag_jsig_01_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig_01_lut" );

  NEW_VECTOR( mid_main_jsig_10_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig_10_lut" );
  NEW_VECTOR( mid_diag_jsig_10_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig_10_lut" );

  NEW_VECTOR( mid_main_jsig_11_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig_11_lut" );
  NEW_VECTOR( mid_diag_jsig_11_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig_11_lut" );


  for( idx = 0; idx <= ZC_MASK; idx++ ) {

    p = ( ( idx >> PA_POS ) & 1 );
    l = ( ( idx >> CL_POS ) & 1 );
    t = ( ( idx >> TC_POS ) & 1 );
    v3 =
      ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 );

    if( p == 0 ) {
      if( ( l && t ) ) {
        ctxt = 3;
        diag_ctxt = 3;
      } else if( l || t ) {
        ctxt = 1;
        diag_ctxt = 1;
      } else {
        ctxt = 0;
        diag_ctxt = 0;
      }
    } else {
      if( ( l && t ) ) {
        ctxt = 3;
        diag_ctxt = 3;
      } else if( l || t ) {
        ctxt = 2;
        diag_ctxt = 2;
      } else if( v3 > 1 ) {
        ctxt = 1;
        diag_ctxt = 2;
      } else {
        ctxt = 1;
        diag_ctxt = 1;
      }
    }

    assert( ctxt < mid_main_jsig_00_cxts );
    mid_main_jsig_00_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig_00_cxts );
    mid_diag_jsig_00_lut[idx] = ( std_byte ) diag_ctxt;


    r = ( ( idx >> CR_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BR_POS ) & 1 );


    if( l == 0 ) {              //sig
      if( p == 0 ) {
        if( r && t ) {
          ctxt = 2;
          diag_ctxt = 2;
        } else if( r || t ) {
          ctxt = 1;
          diag_ctxt = 1;
        } else if( v3 > 1 ) {
          ctxt = 0;
          diag_ctxt = 1;
        } else {
          ctxt = 0;
          diag_ctxt = 0;
        }
      } else {
        if( r || t ) {
          ctxt = 2;
          diag_ctxt = 2;
        } else if( v3 > 1 ) {
          ctxt = 1;
          diag_ctxt = 2;
        } else {
          ctxt = 1;
          diag_ctxt = 1;
        }
      }
    } else {
      if( p == 0 ) {
        if( r && t ) {
          ctxt = 7;
          diag_ctxt = 7;
        } else if( r || t ) {
          ctxt = 4;
          diag_ctxt = 4;
        } else if( v3 > 1 ) {
          ctxt = 3;
          diag_ctxt = 4;
        } else {
          ctxt = 3;
          diag_ctxt = 3;
        }
      } else {
        if( r && t ) {
          ctxt = 7;
          diag_ctxt = 7;
        } else if( r || t ) {
          ctxt = 6;
          diag_ctxt = 6;
        } else if( v3 > 1 ) {
          ctxt = 5;
          diag_ctxt = 6;
        } else {
          ctxt = 5;
          diag_ctxt = 5;
        }
      }
    }

    assert( ctxt < mid_main_jsig_01_cxts );
    mid_main_jsig_01_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig_01_cxts );
    mid_diag_jsig_01_lut[idx] = ( std_byte ) diag_ctxt;


    bc = ( ( idx >> BC_POS ) & 1 );
    tr = ( ( idx >> TR_POS ) & 1 );
    v2 = ( ( idx >> TC_POS ) & 1 ) + ( ( idx >> BC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 ) + ( ( idx >> BR_POS ) & 1 );

    if( ( t == 0 ) && ( tr == 0 ) ) {
      if( p == 0 ) {
        if( bc && l ) {
          ctxt = 3;
          diag_ctxt = 3;
        } else if( bc || l ) {
          ctxt = 1;
          diag_ctxt = 2;
        } else if( v3 ) {
          ctxt = 0;
          diag_ctxt = ( v3 > 1 ) ? 1 : 0;
        } else {
          ctxt = 0;
          diag_ctxt = 0;
        }
      } else {
        if( l ) {
          ctxt = 3;
          diag_ctxt = 3;
        } else if( bc ) {
          ctxt = 2;
          diag_ctxt = 3;
        } else if( v3 > 1 ) {
          ctxt = 1;
          diag_ctxt = 3;
        } else {
          ctxt = 1;
          diag_ctxt = 2;
        }
      }
    } else {
      if( p == 0 ) {
        if( t ) {
          if( l && bc ) {
            ctxt = 8;
            diag_ctxt = 9;
          } else if( l || bc ) {
            ctxt = 7;
            diag_ctxt = 8;
          } else if( v3 ) {
            ctxt = 5;
            diag_ctxt = ( v3 > 1 ) ? 8 : 7;
          } else {
            ctxt = 5;
            diag_ctxt = 5;
          }
        } else {
          if( l || bc ) {
            ctxt = 6;
            diag_ctxt = 6;
          } else if( v3 > 1 ) {
            ctxt = 5;
            diag_ctxt = 6;
          } else {
            ctxt = 4;
            diag_ctxt = 4;
          }
        }
      } else {
        if( t ) {
          if( l && bc ) {
            ctxt = 12;
            diag_ctxt = 13;
          } else if( l || bc ) {
            ctxt = 11;
            diag_ctxt = 12;
          } else if( v3 ) {
            ctxt = 9;
            diag_ctxt = ( v3 > 1 ) ? 12 : 10;
          } else {
            ctxt = 9;
            diag_ctxt = 10;
          }
        } else {
          if( l && bc ) {
            ctxt = 12;
            diag_ctxt = 13;
          } else if( l || bc ) {
            ctxt = 10;
            diag_ctxt = 11;
          } else if( v3 > 1 ) {
            ctxt = 9;
            diag_ctxt = 11;
          } else {
            ctxt = 9;
            diag_ctxt = 10;
          }
        }
      }
    }

    assert( ctxt < mid_main_jsig_10_cxts );
    mid_main_jsig_10_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig_10_cxts );
    mid_diag_jsig_10_lut[idx] = ( std_byte ) diag_ctxt;


    v1 = ( ( idx >> CL_POS ) & 1 ) + ( ( idx >> CR_POS ) & 1 );

    if( p ) {
      if( ( v1 == 0 ) && ( v2 == 0 ) && ( v3 == 1 ) ) {
        ctxt = 0;
        diag_ctxt = 0;
      } else if( ( v1 == 0 ) && ( v2 == 0 ) ) {
        ctxt = 1;
        diag_ctxt = 1;
      } else if( ( v1 == 0 ) && ( v2 == 1 ) ) {
        ctxt = 1;
        diag_ctxt = 1;
      } else if( ( v1 == 1 ) && ( v2 == 0 ) ) {
        ctxt = 2;
        diag_ctxt = 1;
      } else if( ( v1 + v2 ) > 2 ) {
        ctxt = 4;
        diag_ctxt = 3;
      } else {                  //(v1 == 1) && (v2 == 1)
        ctxt = 3;
        diag_ctxt = 2;
      }
    } else {
      if( ( v1 == 0 ) && ( v2 == 0 ) && ( v3 == 1 ) ) {
        ctxt = 5;
        diag_ctxt = 4;
      } else if( ( v1 == 0 ) && ( v2 == 0 ) ) {
        ctxt = 6;
        diag_ctxt = 5;
      } else if( ( v1 == 0 ) && ( v2 == 1 ) ) {
        ctxt = 6;
        diag_ctxt = 5;
      } else if( ( v1 == 1 ) && ( v2 == 0 ) ) {
        ctxt = 6;
        diag_ctxt = 5;
      } else if( ( v1 + v2 ) > 2 ) {
        ctxt = 8;
        diag_ctxt = 7;
      } else {                  //(v1 == 1) && (v2 == 1)
        ctxt = 7;
        diag_ctxt = 6;
      }
    }

    assert( ctxt < mid_main_jsig_11_cxts );
    mid_main_jsig_11_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig_11_cxts );
    mid_diag_jsig_11_lut[idx] = ( std_byte ) diag_ctxt;

  }
}


//-----------------------------------------------------------------------------

// initialize_mid_jsig0_00_luts(void)

//----------------------------------------------------------------------------

void
SubbandCodec::initialize_mid_jsig0_00_luts( void )
{
  std_short idx, ctxt, v0, v1, v2, v3, diag_ctxt;

  NEW_VECTOR( mid_main_jsig0_00_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig0_00_lut" );

  NEW_VECTOR( mid_diag_jsig0_00_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig0_00_lut" );

  for( idx = 0; idx <= ZC_MASK; idx++ ) {

    v0 = ( ( idx >> PA_POS ) & 1 );
    v1 = ( ( idx >> CL_POS ) & 1 );
    v2 = ( ( idx >> TC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 );

    if( v1 && v2 ) {
      ctxt = 5;
      diag_ctxt = 3;
    } else if( v1 ) {
      ctxt = ( v0 ) ? 4 : 2;
      diag_ctxt = 2;
    } else if( v2 ) {
      ctxt = ( v0 ) ? 3 : 2;
      diag_ctxt = 2;
    } else if( v3 > 1 ) {
      ctxt = 1;
      diag_ctxt = 2;
    } else if( v3 ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else {
      ctxt = 0;
      diag_ctxt = 0;
    }

    assert( ctxt < mid_main_jsig0_00_cxts );
    mid_main_jsig0_00_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig0_00_cxts );
    mid_diag_jsig0_00_lut[idx] = ( std_byte ) diag_ctxt;
  }

}


void
SubbandCodec::initialize_mid_jsig0_01_luts( void )
{
  std_short idx, ctxt, p, l, r, t, v3, diag_ctxt;

  NEW_VECTOR( mid_main_jsig0_01_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig0_01_lut" );

  NEW_VECTOR( mid_diag_jsig0_01_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig0_01_lut" );

  for( idx = 0; idx <= ZC_MASK; idx++ ) {
    p = ( ( idx >> PA_POS ) & 1 );
    l = ( ( idx >> CL_POS ) & 1 );
    r = ( ( idx >> CR_POS ) & 1 );
    t = ( ( idx >> TC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BR_POS ) & 1 );

    if( l == 0 ) {
      if( r && t ) {
        ctxt = 4;
        diag_ctxt = 3;
      } else if( r ) {
        ctxt = ( p ) ? 3 : 2;
        diag_ctxt = 2;
      } else if( t ) {
        ctxt = ( p ) ? 3 : 2;
        diag_ctxt = 2;
      } else if( v3 > 1 ) {
        ctxt = 1;
        diag_ctxt = 2;
      } else if( v3 ) {
        ctxt = 1;
        diag_ctxt = 1;
      } else {
        ctxt = 0;
        diag_ctxt = 0;
      }
    } else {
      if( r && t ) {
        ctxt = 7;
        diag_ctxt = 6;
      } else if( r || t ) {
        ctxt = ( p ) ? 7 : 6;
        diag_ctxt = ( p ) ? 6 : 5;
      } else if( v3 > 1 ) {
        ctxt = 5;
        diag_ctxt = ( p ) ? 6 : 5;
      } else {
        ctxt = 5;
        diag_ctxt = 4;
      }
    }

    assert( ctxt < mid_main_jsig0_01_cxts );
    mid_main_jsig0_01_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig0_01_cxts );
    mid_diag_jsig0_01_lut[idx] = ( std_byte ) diag_ctxt;

  }
}



void
SubbandCodec::initialize_mid_jsig0_10_luts( void )
{
  std_short idx, ctxt, p, l, t, tr, bc, v2, v3, diag_ctxt;      // r,

  NEW_VECTOR( mid_main_jsig0_10_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig0_10_lut" );

  NEW_VECTOR( mid_diag_jsig0_10_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig0_10_lut" );

  for( idx = 0; idx <= ZC_MASK; idx++ ) {

    p = ( ( idx >> PA_POS ) & 1 );
    l = ( ( idx >> CL_POS ) & 1 );
    t = ( ( idx >> TC_POS ) & 1 );
    bc = ( ( idx >> BC_POS ) & 1 );
    tr = ( ( idx >> TR_POS ) & 1 );

    v2 = ( ( idx >> TC_POS ) & 1 ) + ( ( idx >> BC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 ) + ( ( idx >> BR_POS ) & 1 );

    if( ( t == 0 ) && ( tr == 0 ) ) {
      if( p == 0 ) {
        if( bc && l ) {
          ctxt = 5;
          diag_ctxt = 4;
        } else if( bc || l ) {
          ctxt = 1;
          diag_ctxt = 1;
        } else if( v3 ) {
          ctxt = 2;
          diag_ctxt = ( v3 > 1 ) ? 1 : 2;
        } else {
          ctxt = 0;
          diag_ctxt = 0;
        }
      } else {
        if( bc && l ) {
          ctxt = 5;
          diag_ctxt = 4;
        } else if( l ) {
          ctxt = 4;
          diag_ctxt = 3;
        } else if( bc ) {
          ctxt = 3;
          diag_ctxt = 3;
        } else if( v3 > 1 ) {
          ctxt = 2;
          diag_ctxt = 3;
        } else if( v3 ) {
          ctxt = 2;
          diag_ctxt = 2;
        } else {
          ctxt = 0;
          diag_ctxt = 0;
        }
      }
    } else {
      if( p == 0 ) {
        if( t ) {
          if( l && bc ) {
            ctxt = 19;
            diag_ctxt = 18;
          } else if( l ) {
            ctxt = 17;
            diag_ctxt = 17;
          } else if( bc ) {
            ctxt = 17;
            diag_ctxt = 17;
          } else if( v3 ) {
            ctxt = 15;
            diag_ctxt = ( v3 > 1 ) ? 17 : 15;
          } else {
            ctxt = 15;
            diag_ctxt = 14;
          }
        } else {
          if( l && bc ) {
            ctxt = 12;
            diag_ctxt = 12;
          } else if( l ) {
            ctxt = 12;
            diag_ctxt = 12;
          } else if( bc ) {
            ctxt = 12;
            diag_ctxt = 12;
          } else if( v3 > 1 ) {
            ctxt = 15;
            diag_ctxt = 12;
          } else {
            ctxt = 10;
            diag_ctxt = 10;
          }
        }
      } else {
        if( t ) {
          if( l && bc ) {
            ctxt = 24;
            diag_ctxt = 23;
          } else if( l ) {
            ctxt = 27;
            diag_ctxt = 27;
          } else if( bc ) {
            ctxt = 27;
            diag_ctxt = 27;
          } else if( v3 ) {
            ctxt = 25;
            diag_ctxt = ( v3 > 1 ) ? 27 : 24;
          } else {
            ctxt = 25;
            diag_ctxt = 24;
          }
        } else {
          if( l && bc ) {
            ctxt = 24;
            diag_ctxt = 23;
          } else if( l ) {
            ctxt = 22;
            diag_ctxt = 22;
          } else if( bc ) {
            ctxt = 22;
            diag_ctxt = 22;
          } else if( v3 > 1 ) {
            ctxt = 25;
            diag_ctxt = ( v3 == 2 ) ? 22 : 22;
          } else {
            ctxt = 25;
            diag_ctxt = 24;
          }
        }
      }
    }

    assert( ctxt < mid_main_jsig0_10_cxts );
    mid_main_jsig0_10_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig0_10_cxts );
    mid_diag_jsig0_10_lut[idx] = ( std_byte ) diag_ctxt;
  }
}



void
SubbandCodec::initialize_mid_jsig0_11_luts( void )
{
  std_short idx, ctxt, v1, v2, v3, diag_ctxt, v0;


  NEW_VECTOR( mid_main_jsig0_11_lut, ZC_MASK + 1, std_byte,
              "mid_main_jsig0_11_lut" );
  NEW_VECTOR( mid_diag_jsig0_11_lut, ZC_MASK + 1, std_byte,
              "mid_diag_jsig0_11_lut" );


  for( idx = 0; idx <= ZC_MASK; idx++ ) {

    v0 = ( ( idx >> PA_POS ) & 1 );
    v1 = ( ( idx >> CL_POS ) & 1 ) + ( ( idx >> CR_POS ) & 1 );
    v2 = ( ( idx >> TC_POS ) & 1 ) + ( ( idx >> BC_POS ) & 1 );
    v3 = ( ( idx >> TL_POS ) & 1 ) + ( ( idx >> TR_POS ) & 1 ) +
      ( ( idx >> BL_POS ) & 1 ) + ( ( idx >> BR_POS ) & 1 );

    if( ( v1 == 0 ) && ( v2 == 0 ) && ( v3 == 1 ) ) {
      ctxt = 0;
      diag_ctxt = 0;
    } else if( ( v1 == 0 ) && ( v2 == 0 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else if( ( v1 == 0 ) && ( v2 == 1 ) ) {
      ctxt = 1;
      diag_ctxt = 1;
    } else if( ( v1 == 1 ) && ( v2 == 0 ) ) {
      ctxt = 2;
      diag_ctxt = 1;
    } else if( ( v1 == 0 ) && ( v2 == 2 ) ) {
      ctxt = 3;
      diag_ctxt = 2;
    } else if( ( v1 == 2 ) && ( v2 == 0 ) ) {
      ctxt = 3;
      diag_ctxt = 2;
    } else if( ( v1 + v2 ) > 2 ) {
      ctxt = 4;
      diag_ctxt = 3;
    } else {                    //(v1 == 1) && (v2 == 1)
      ctxt = 3;
      diag_ctxt = 2;
    }

    assert( ctxt < mid_main_jsig0_11_cxts );
    mid_main_jsig0_11_lut[idx] = ( std_byte ) ctxt;

    assert( diag_ctxt < mid_diag_jsig0_11_cxts );
    mid_diag_jsig0_11_lut[idx] = ( std_byte ) diag_ctxt;

  }
}
