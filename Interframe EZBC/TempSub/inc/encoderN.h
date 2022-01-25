#include "structN.h"

EXTERN videoinfo info;
EXTERN videoheader header;

EXTERN YUVimage efr;            /* prediction error frame */
EXTERN YUVimage pcode;          /* coded previous frame */

EXTERN FILE *fpbit;
EXTERN Rate FrsRate;

EXTERN vector_ptr tmp_yfmv;        /* motion vector */
EXTERN vector_ptr *yfmv;           /* motion vector */
EXTERN vector_ptr *yfmv_bigGOP;    /* motion vector for bigGOP */
EXTERN YUVimage_ptr pyrFrs;        /* frames after temporal decomposition */
EXTERN YUVimage_ptr pyrFrs_bigGOP; /* frames of bigGOP */
EXTERN YUVimage_ptr *pyrTemp;      /* temporary memory for temporally filtering */
EXTERN float **frY;                /* pointer for coding */
EXTERN float **frU;
EXTERN float **frV;
