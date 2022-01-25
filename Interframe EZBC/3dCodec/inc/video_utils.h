/* ========================================================================= */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */
#ifndef _VIDEO_UTILS_H
#define _VIDEO_UTILS_H


#include "struct_sht.h"

int alloc_pyr_frames( YUVimage * pyrFrs, videoinfo & info );
void free_frames( YUVimage * pyrFrs, int n );
int get_mvBytes( FILE * fp_mv );
 
#endif
