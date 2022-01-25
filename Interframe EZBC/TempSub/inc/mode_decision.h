#ifndef MODE_DECISION_H
#define MODE_DECISION_H

#include "structN.h"

/*void find_best_mode(vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2, 
                    float *fr_cur, float *fr_ref1, float *fr_ref2, 
                    float *upframe1, float *upframe2,
                    int x, int y, int xblk, int yblk, int maxx, int maxy, 
                    int hor, int ver, int t_level, videoinfo info,
                    int ctx1x, int ctx1y, int ctx2x, int ctx2y, 
                    float *pmv1x, float *pmv1y, float *pmv2x, float *pmv2y, int dec);*/

void find_best_mode(vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3_array, vector_ptr fmv4_array,
                    float *fr_cur, float *fr_ref1, float *fr_ref2, 
                    float *upframe1, float *upframe2,
                    int x, int y, int xblk, int yblk, int maxx, int maxy, 
                    int hor, int ver, int t_level, videoinfo info,
                    int ctx1x, int ctx1y, int ctx2x, int ctx2y, 
                    float *pmv1x, float *pmv1y, float *pmv2x, float *pmv2y, int dec);

vector_ptr find_block(int x, int y, vector_ptr root, videoinfo info, int t_level, float *mvx, 
					  float *mvy, int *get_xblk, int *get_yblk, int get_addx, int get_addy);

int compare_mv(vector_ptr left1, vector_ptr right1, vector_ptr left2, vector_ptr right2);

#endif // MODE_DECISION_H
