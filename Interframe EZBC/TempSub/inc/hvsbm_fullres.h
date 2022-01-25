#ifndef HVSBM_FULLRES_H
#define HVSBM_FULLRES_H

#include "structN.h"

float hvsbm_fullres(vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3, vector_ptr fmv4, 
                    float *fr_cur, float *fr_ref1, float *fr_ref2, 
                    float *upframe1, float *upframe2,
                    int t_level, int dist, videoinfo info, int dec);

void rdme(vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3, vector_ptr fmv4, 
          YUVimage *fr_cur, YUVimage *fr_ref1, YUVimage *fr_ref2, 
          enum FLAG *sc1, enum FLAG *sc2, int t_level, int dist, videoinfo info, 
          float *upframe1, float *upframe2, float *upsamp_x );

int get_searchrange(int dist, videoinfo info);

#endif // HVSBM_FULLRES_H

