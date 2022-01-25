/* ========================================================================= */
/* Description: definition and configuration of EZBC                         */
/* Notes: modified from ebcot_constants.h                                    */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

/* ========================================================================= */
/*****************************************************************************/
/* Copyright 1998, Hewlett-Packard Company                                   */
/* All rights reserved                                                       */
/* File: "ebcot_constants.h"                                                 */
/* Description: Constants which configure the EBCOT encoder and decoder      */
/* Author: David Taubman                                                     */
/* Affiliation: Hewlett-Packard and                                          */
/*              The University of New South Wales, Australia                 */
/* Version: V3.1A                                                            */
/* Last Revised: 14 January, 1999                                            */
/*****************************************************************************/
#ifndef _CONTEXT_CONSTANTS_H
#define _CONTEXTCONSTANTS_H
#include <ifc.h>

/* ========================================================================= */
/* ---------------------------- Context Assignment ------------------------- */
/* ========================================================================= */

#define TL_POS    0             /* Bit-pos for significance of top-left neighbour        */
#define TC_POS    1             /* Bit-pos for significance of top-centre neighbour      */

#define TC2_POS   2             /*2nd Bit_pos for significance of top-centre neighbour   */

#define TR_POS    3             /* Bit-pos for significance of top-right neighbour       */
#define CL_POS    4             /* Bit-pos for significance of centre-left neighbour     */
#define CL2_POS   5             /*2nd Bit-pos for significance of centre-left neighbour  */
#define CR_POS    6             /* Bit-pos for significance of centre-right neighbour    */
#define CR2_POS   7             /*2nd Bit-pos for significance of centre-right neighbour */
#define BL_POS    8             /* Bit-pos for significance of bottom-left neighbour     */
#define BC_POS    9             /* Bit-pos for significance of bottom-centre neighbour   */
#define BC2_POS  10             /*2nd Bit-pos for significance of bottom-centre neighbour */
#define BR_POS   11             /* Bit-pos for significance of bottom-right neighbour    */

#define OUT_OF_BOUNDS_POS 15    /* May be used to identify context words which
                                   lie beyond the boundaries of the code block */

#define TL_SIG        ((std_short)(1<<TL_POS))
#define TC_SIG        ((std_short)(1<<TC_POS))
#define TR_SIG        ((std_short)(1<<TR_POS))
#define CL_SIG        ((std_short)(1<<CL_POS))
#define CR_SIG        ((std_short)(1<<CR_POS))
#define BL_SIG        ((std_short)(1<<BL_POS))
#define BC_SIG        ((std_short)(1<<BC_POS))
#define BR_SIG        ((std_short)(1<<BR_POS))

#define OUT_OF_BOUNDS ((std_short)(1<<OUT_OF_BOUNDS_POS))

/* ========================================================================= */
//-----------------------------------------------------------------------------

// Choose operation modes

//-----------------------------------------------------------------------------

//#define DEBUG

#define TEST_CODING_RATE        //test if target rate reached in subplane pass

#define INTERBANDS  // INTERBANDS

#define GET_PARENT_MODELS // INTERBANDS


#define INITIALIZE_JSIG_MODELS_FROM_PAR // INTERBANDS

#define INITIALIZE_LIP_MODELS_FROM_PAR // INTERBANDS
#define INITIALIZE_NODE_MODELS_FROM_PAR // INTERBANDS


//sign

#define SCALE_SIGN_MODELS_AT_LEVEL_ENDS
#define UPDATE_SIGN_DIAG_CXT
//#define MID_SIGN_LUTS

//LSP
#define INITIALIZE_LSP_MODELS_FROM_PAR // INTERBANDS
#define MERGE_ALL_LSP_PLANES
#define TWO_SIGNIF_BITS



#ifdef MERGE_ALL_LSP_PLANES

#define LSP_LUT
#define LSP_BIT_IDX
#define SCALE_AT_LSP_BREAK_PTS
#define MERGE_HIGH_IDX_SETS_OF_TOP_LSP_PLANE
#define MERGE_SETS_OF_TOP_LSP_PLANE
#define MERGE_SETS_0_1_OF_TOP_PLANE
#define SEPARATE_SETS_OF_PLANE_1
#define MERGE_LOWER_LSP_PLANE
#define MERGE_LSP_PLANE_1

#endif


#ifdef NOT_INITIALIZE_FROM_PAR

#undef INITIALIZE_JSIG_MODELS_FROM_PAR
#undef INITIALIZE_LIP_MODELS_FROM_PAR
#undef INITIALIZE_NODE_MODELS_FROM_PAR
#undef INITIALIZE_LSP_MODELS_FROM_PAR
//#undef INITIALIZE_SIGN_MODELS_FROM_PAR

#endif

/* ========================================================================= */

#define WIDTH_OF_CXT_BDY 2
#define WIDTH_OF_CXT_BDY_X_2 (WIDTH_OF_CXT_BDY << 1)


#define V_PVE_BIT_POS 0
#define V_NVE_BIT_POS 1
#define H_PVE_BIT_POS 2
#define H_NVE_BIT_POS 3

#define NW_PVE_BIT_POS 4
#define NW_NVE_BIT_POS 5
#define NE_PVE_BIT_POS 6
#define NE_NVE_BIT_POS 7

#define V_PVE_MASK ((std_byte)(1 << V_PVE_BIT_POS))
#define V_NVE_MASK ((std_byte)(1 << V_NVE_BIT_POS))
#define H_PVE_MASK ((std_byte)(1 << H_PVE_BIT_POS))
#define H_NVE_MASK ((std_byte)(1 << H_NVE_BIT_POS))

#define NW_PVE_MASK ((std_byte)(1 << NW_PVE_BIT_POS))
#define NW_NVE_MASK ((std_byte)(1 << NW_NVE_BIT_POS))
#define NE_PVE_MASK ((std_byte)(1 << NE_PVE_BIT_POS))
#define NE_NVE_MASK ((std_byte)(1 << NE_NVE_BIT_POS))

#define SIGN_BIT_POS V_PVE_BIT_POS

#define SC_MASK_EBCOT ((std_byte)(~((-1) << 4)))        /*Sign-coding LUT index mask. */

#define PA_POS 12               //Bit-pos for significance of parent pel


#define PA_SIG ((std_short)(1<<PA_POS))

#define TC2_SIG ((std_short)(1<<TC2_POS))
#define CL2_SIG ((std_short)(1<<CL2_POS))
#define CR2_SIG ((std_short)(1<< CR2_POS))
#define BC2_SIG ((std_short)(1<< BC2_POS))

#define ZC_CXT_BITS 12


/* Zero-coding LUT index mask. */
#define ZC_MASK ((std_short)(~((-1) << (ZC_CXT_BITS+1))))

#define SIGN_CXT_BITS 8
#define SIGN_CXT_MASK ((std_byte)(~((-1) << SIGN_CXT_BITS)))

#define DIGITS_11 ((std_short)0x3)

#endif /* EBCOT_CONSTANTS_H */
