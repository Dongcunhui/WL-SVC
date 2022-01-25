#ifndef MCTF_OBMC_H
#define MCTF_OBMC_H


#include "structN.h"

void get_mv_side_information(videoinfo info, vector_ptr fmv2, vector_ptr fmv3, 
							 vector_ptr mv_ref2, vector_ptr mv_ref3, 
							 ImageMEinfo *imagemeinfo, Varblkarrayinfo *imageblkarray, 
							 int *total_blk, int i, int j, 
							 enum FLAG left_scene, enum FLAG right_scene, int encoder_sign);

void mv_weight_info(videoinfo info, ImageMEinfo *pixelmeinfo, Varblkarrayinfo *blockinfoarray,
					int total_varblk, int t_level, int tindex);

#endif