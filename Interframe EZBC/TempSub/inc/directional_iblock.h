
#ifndef DIRECTIONAL_IBLOCK_H
#define DIRECTIONAL_IBLOCK_H


#include "structN.h"

//Added by Yuan Liu
void  rec_directional_iblock_analysis(vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
									  int hor, int ver, 
									YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
					                YUVimage_ptr  fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
									vector_ptr    fmv1, vector_ptr fmv2, vector_ptr fmv3,
									vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
									int t_level, int remaining_frs, videoinfo info,
									vector_ptr fmv_root,
									ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray);
//Added by Yuan Liu

void directional_iblock_detection(videoinfo info, vector_ptr fmv1_array, float *fr_cur,  
								  int cx, int cy, int xblock, int yblock, int hor, int ver, int t_level, 
								  vector *tmv1,  float *fadditional_penalty );

void  rec_directional_iblock_analysis_with_OBMC(vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
									  int hor, int ver, 
									  YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
					                  YUVimage_ptr  fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
									  vector_ptr    fmv1, vector_ptr fmv2, vector_ptr fmv3,
									  vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
									  int t_level, int remaining_frs, videoinfo info,
									  vector_ptr fmv_root, 
									  ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray);

void rec_directional_iblock_synthesis_with_OBMC( vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
									   int hor, int ver, 
									   YUVimage_ptr fr0, YUVimage_ptr fr1, YUVimage H0,
									   YUVimage L1, YUVimage H1, YUVimage frp,
                                       vector_ptr fmv0, vector_ptr fmv1, vector_ptr fmv2,
                                       vector_ptr mv_ref0, vector_ptr mv_ref1, 
                                       vector_ptr mv_ref2, int level, videoinfo info, 
									   vector_ptr fmv_root,
									   ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray);

#endif