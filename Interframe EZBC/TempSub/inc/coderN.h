#pragma once

#ifndef CODERN_H
#define CODERN_H

#define ADAPT_LEVEL 5

EXTERN YUVimage end_of_lastGOP;
//Added on 02.03.2018
EXTERN YUVimage buff_end_of_lastGOP;

EXTERN FILE *fpbit;
EXTERN Rate FrsRate;

EXTERN vector_ptr tmp_yfmv;
EXTERN vector_ptr *yfmv;           /* motion vector */
EXTERN vector_ptr *yfmv_bigGOP;    /* only help for different GOP */
EXTERN vector_ptr *mv_ref;         /* arrays of MV references */
EXTERN vector_ptr *mv_ref_bigGOP;

//Added on 02.03.2018
EXTERN vector_ptr buff_tmp_yfmv;
EXTERN vector_ptr *buff_yfmv;           /* motion vector */
EXTERN vector_ptr *buff_yfmv_bigGOP;    /* only help for different GOP */
EXTERN vector_ptr *buff_mv_ref;         /* arrays of MV references */
EXTERN vector_ptr *buff_mv_ref_bigGOP;

EXTERN vector_ptr *dec_yfmv;           /* motion vector */
EXTERN vector_ptr *dec_yfmv_bigGOP;    /* only help for different GOP */
EXTERN vector_ptr *dec_mv_ref_bigGOP;

//Added on 02.03.2018
EXTERN vector_ptr *buff_dec_yfmv;           /* motion vector */
EXTERN vector_ptr *buff_dec_yfmv_bigGOP;    /* only help for different GOP */
EXTERN vector_ptr *buff_dec_mv_ref_bigGOP;

EXTERN YUVimage_ptr pyrFrs;        /* frames after temporal decomposition */
EXTERN YUVimage_ptr pyrFrs_bigGOP; /* only help for different GOP */
EXTERN YUVimage_ptr pyrFrs_first;  /* first L-Frame by decoding */

//Added on 02.03.2018
EXTERN YUVimage_ptr buff_pyrFrs;        /* frames after temporal decomposition */
EXTERN YUVimage_ptr buff_pyrFrs_bigGOP; /* only help for different GOP */
EXTERN YUVimage_ptr buff_pyrFrs_first;  /* first L-Frame by decoding */

EXTERN YUVimage_ptr dec_pyrFrs;        /* frames after temporal decomposition 解码端的Frs*/
EXTERN YUVimage_ptr dec_pyrFrs_bigGOP; /* only help for different GOP */
EXTERN YUVimage_ptr dec_pyrFrs_first;  /* first L-Frame by decoding */

//Added on 02.03.2018
EXTERN YUVimage_ptr buff_dec_pyrFrs;        /* frames after temporal decomposition */
EXTERN YUVimage_ptr buff_dec_pyrFrs_bigGOP; /* only help for different GOP */
EXTERN YUVimage_ptr buff_dec_pyrFrs_first;  /* first L-Frame by decoding */

EXTERN enum FLAG **scene_change;  // 第一维表示时域分解的第几层， 第二维表示一层时域分解内部的每个帧
//Added on 02.03.2018
EXTERN enum FLAG **buff_scene_change;

EXTERN enum FLAG **simul_scene_change;

EXTERN enum FLAG scene_change_help;
//Added on 02.03.2018
EXTERN enum FLAG buff_scene_change_help;

EXTERN enum FLAG **dec_scene_change;
//Added on 02.03.2018
EXTERN enum FLAG **buff_dec_scene_change;

EXTERN YUVimage_ptr *pyrTemp;             /* temporary memory for temporally filtering */
//Added on 02.03.2018
EXTERN YUVimage_ptr *buff_pyrTemp;             /* temporary memory for temporally filtering */

EXTERN YUVimage_ptr Four_GOP, Four_bigGOP; /* only used if info.denoise_flag == YES */
EXTERN YUVimage_ptr spatial_high[3];      /* only used if info.denoise_flag == YES */

EXTERN YUVimage_ptr *dec_pyrTemp;   
//Added on 02.03.2018
EXTERN YUVimage_ptr *buff_dec_pyrTemp;   

EXTERN float HPW1[3], HPW2[3], HPW3[3], HPW4[3], /* filter for temporal decomposition */ // 依次为双向滤波系数、前向、后向、
  HPW1_pred[3], HPW2_pred[3], HPW3_pred[3], HPW4_pred[3],
  LPW1[3], LPW2[3], LPW3[3], LPW4[3]; // 后面四个是低频帧的系数，共四种模式，双向、前向、后向、帧内

EXTERN float FIR14[FIR_LEN], FIR12[FIR_LEN], FIR34[FIR_LEN], FIR18[FIR_LEN], FIR38[FIR_LEN], FIR58[FIR_LEN],
  FIR78[FIR_LEN], FIR1hex[FIR_LEN], FIR3hex[FIR_LEN], FIR5hex[FIR_LEN], FIR7hex[FIR_LEN], FIR9hex[FIR_LEN],
  FIR11hex[FIR_LEN], FIR13hex[FIR_LEN], FIR15hex[FIR_LEN];

/* Check */
EXTERN float **frY;             /*pointer for coding */

EXTERN float **frU;
EXTERN float **frV;

EXTERN int connected, unreferred, pre_connected, multi_connected,
  next_connected, case4;
EXTERN long int totalY[5], totalU[5], totalV[5], totalmap[5], totalMV[5];
EXTERN float avgpsnr[3][16];    /* avgpsnr[y,u,v][frames in several levels] averaged psnr over all GOPS */
EXTERN float avgpsnr_cod[3][16];        /* avgpsnr_cod[y,u,v][frames in several levels] averaged psnr over all GOPS */
EXTERN float avgvar[3][16];     /* avgvar[y,u,v][frames in several levels]  averaged variance over all GOPS */

EXTERN long int *unconnectedL;


// by Yongjun Wu 
// for directional spatial prediction from RPI
#define  IBLOCK_MAX_SIZE 8      // the maximum block size for directional IBLOCK
EXTERN float neighborA[IBLOCK_MAX_SIZE];
EXTERN float neighborB[IBLOCK_MAX_SIZE];
EXTERN float neighborC[IBLOCK_MAX_SIZE];
EXTERN float neighborD[IBLOCK_MAX_SIZE];
EXTERN float neighborD1[IBLOCK_MAX_SIZE];
EXTERN float neighborE[IBLOCK_MAX_SIZE];
EXTERN float neighborF[IBLOCK_MAX_SIZE];
EXTERN float neighborG[IBLOCK_MAX_SIZE];
EXTERN float neighborH;
EXTERN float predict_blkY[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE]; 
EXTERN float neighborUA[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUB[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUC[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUD[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUD1[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUE[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUF[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUG[IBLOCK_MAX_SIZE/2];
EXTERN float neighborUH;
EXTERN float predict_blkU[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE/4]; 
EXTERN float neighborVA[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVB[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVC[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVD[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVD1[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVE[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVF[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVG[IBLOCK_MAX_SIZE/2];
EXTERN float neighborVH;
EXTERN float predict_blkV[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE/4]; 

// for IBLOCK in OBMC framework 
EXTERN float neighbor_predictY[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE]; 
EXTERN float neighbor_predictU[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE/4]; 
EXTERN float neighbor_predictV[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE/4]; 
EXTERN float self_weight_matrix[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE];
EXTERN float self_weight_matrixUV[IBLOCK_MAX_SIZE*IBLOCK_MAX_SIZE/4]; 

EXTERN Varblkarrayinfo  *varblkarray;
EXTERN ImageMEinfo      *frameMEinfo; 

EXTERN ImageMEinfo      *buff_frameMEinfo;  //Added on 05.22.2018
EXTERN float gop_block; //Added on 06.28.2018
EXTERN int frame_cnt;  //Added on 06.28.2018
EXTERN int simul_skip; //Added on 06.30.2018

EXTERN int skip_frame; //Added on 06.30.2018
EXTERN float orig_lambda[ADAPT_LEVEL];

/* Check */

#endif // CODERN_H
