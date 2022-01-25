#include "structN.h"

EXTERN videoinfo info;
EXTERN videoheader header;

EXTERN YUVimage efr;            /* prediction error frame */
EXTERN YUVimage pcode;          /* coded previous frame */

EXTERN FILE *fpbit;

EXTERN vector_ptr *yfmv;        /* motion vector */
EXTERN vector_ptr *yfmv_bigGOP;
EXTERN vector_ptr *mv_ref;      /* arrays of MV references */
EXTERN vector_ptr *mv_ref_bigGOP;

EXTERN YUVimage_ptr pyrFrs;     /*frames after temporal decomposition */
EXTERN YUVimage_ptr pyrFrs_first; /* first L-frames in each GOP */
EXTERN YUVimage_ptr pyrFrs_bigGOP;
EXTERN YUVimage_ptr *pyrTemp;   /*temporary memory for temporally filtering */
