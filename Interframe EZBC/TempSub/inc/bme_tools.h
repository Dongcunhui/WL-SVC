#ifndef BME_TOOLS_H
#define BME_TOOLS_H

#include "structN.h"

void get_dmv(float *dmvx, float *dmvy, enum FLAG *is_predictor,
             vector_ptr fmv, int x_dst, int y_dst, videoinfo info, int t_level, int blk_thresh);

void get_predictor(float *pmvx, float *pmvy, enum FLAG *is_predictor,
                   vector_ptr fmv, int x_dest, int y_dest, videoinfo info, int t_level,
			       enum BiMode *block_mode, int *propagate_iblk, int blk_thresh);  
                  // by Yongjun Wu: get the block mode

void clear_predictors(vector_ptr fmv, videoinfo info, int t_level);

float get_bit_cost(float lambda, float mvx, float mvy, float pmvx, float pmvy,
                   int ctx_x, int ctx_y, int subpel);

void get_median_predictor(float *pmvx, float *pmvy, vector_ptr fmv, vector_ptr prev_fmv,
						  vector_ptr prev_fmv2, int x_pos, int y_pos, int xblock, int yblock,
                          videoinfo info, int t_level, int blk_thresh);

//Added by Yuan Liu
int enc_trans_info(int aff_idx);

#endif
