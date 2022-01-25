#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <fstream>
#define EXTERN  extern
#include "basic.h"
#include "rasterfile.h"
#include "structN.h"
#include "coderN.h"
#include "analsyn.h"
#include "ioN.h"
#include "miscN.h"
#include "memoryN.h"
#include "chrom.h"
#include "mvcodingN.h"
#include "ezbc_dec_3d.h"
#include "bmeN.h"
#include "pstatN.h"
#include "video_utils.h"

#include "directional_iblock.h"

#include "obmc_varblk.h"

#include "util_filtering.h"


void calsnr_frame( YUVimage_ptr fr0, YUVimage_ptr fr1, int curr,
                   videoinfo info );

void ezbc3d_dec_GOP( YUVimage * pyrFrs, videoinfo info, long total_bytes_past,
                     long int GOP_counter, int curr );


void
decode_MV( long int *total_bytes_past, videoinfo info, int GOP_counter, int simul_dec, int theo_dec )
{
  // this function should go to mvcodingN.c
  // KH, 2003-12-02

  long int starting_pos; 

  int mvBytes, GOPheader_bytes, s_level;
  s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
  
  if( s_level > 0 ) { // 调整分辨率的大小
    info.ywidth  <<= s_level;
    info.yheight <<= s_level;
    info.cwidth  <<= s_level;
    info.cheight <<= s_level;
  }

  if( info.tPyrLev >= 1 ) {  
    if( !( fpbit = fopen( info.bitname, "rb" ) ) ) { // 打开bit文件
      printf( "init_dec: %s\n", info.bitname );
      exit( 1 );
    }
	
    fseek( fpbit, *total_bytes_past, SEEK_SET );// 定位文件位置
    GOPheader_bytes = read_GOPheader( dec_scene_change, info );     // 标记每一帧是否有场景改变
    *total_bytes_past += GOPheader_bytes;
	starting_pos = *total_bytes_past; //后移GOPheader_bytes的定位
    mvBytes = mv_decoding( info, dec_yfmv, GOP_counter, starting_pos, simul_dec, theo_dec );// 解码mv
    *total_bytes_past += mvBytes;       
    fclose( fpbit );
  }
}

// 合成一层的重建，不是一帧
void
temporal_synthesis( YUVimage_ptr fr0, YUVimage_ptr fr1, YUVimage H0,
                    YUVimage L1, YUVimage H1, YUVimage frp,
                    vector_ptr fmv0, vector_ptr fmv1, vector_ptr fmv2,
                    vector_ptr mv_ref0, vector_ptr mv_ref1, 
                    vector_ptr mv_ref2, int level, videoinfo info )
     /* fr1--current, fr0--reference */
{
  int yhor, yver, chor, cver;
  float *ymvx0, *ymvy0, *ymvx1, *ymvy1, *ymvx2, *ymvy2;
  float *cmvx0, *cmvy0, *cmvx1, *cmvy1, *cmvx2, *cmvy2;
  float *ymvx0_int, *ymvy0_int, *ymvx1_int, *ymvy1_int, *ymvx2_int, *ymvy2_int;
  float *cmvx0_int, *cmvy0_int, *cmvx1_int, *cmvy1_int, *cmvx2_int, *cmvy2_int;

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;
   
  blockmv2pixelmv( fmv0, &ymvx0, &ymvy0, &cmvx0, &cmvy0, CONNECTED, info, level );
  blockmv2pixelmv( fmv1, &ymvx1, &ymvy1, &cmvx1, &cmvy1, CONNECTED, info, level );
  blockmv2pixelmv( fmv2, &ymvx2, &ymvy2, &cmvx2, &cmvy2, CONNECTED, info, level );
  blockmv2pixelmv( fmv0, &ymvx0_int, &ymvy0_int, &cmvx0_int, &cmvy0_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv1, &ymvx1_int, &ymvy1_int, &cmvx1_int, &cmvy1_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv2, &ymvx2_int, &ymvy2_int, &cmvx2_int, &cmvy2_int, PREDICTED, info, level );

  mc_synthesis( fr0->Y, fr1->Y, H0.Y, L1.Y, H1.Y, frp.Y, 
                ymvx0, ymvy0, ymvx1, ymvy1, ymvx2, ymvy2, 
                ymvx0_int, ymvy0_int, ymvx1_int, ymvy1_int,
                ymvx2_int, ymvy2_int, mv_ref0, mv_ref1, mv_ref2, 
                yhor, yver, level, info );
    
  if( info.cwidth && info.cheight ) 
    {
      mc_synthesis( fr0->U, fr1->U, H0.U, L1.U, H1.U, frp.U, 
                    cmvx0, cmvy0, cmvx1, cmvy1, cmvx2, cmvy2,
                    cmvx0_int, cmvy0_int, cmvx1_int, cmvy1_int, 
                    cmvx2_int, cmvy2_int, mv_ref0, mv_ref1, mv_ref2, 
                    chor, cver, level, info );
      mc_synthesis( fr0->V, fr1->V, H0.V, L1.V, H1.V, frp.V, 
                    cmvx0, cmvy0, cmvx1, cmvy1, cmvx2, cmvy2, 
                    cmvx0_int, cmvy0_int, cmvx1_int, cmvy1_int, 
                    cmvx2_int, cmvy2_int, mv_ref0, mv_ref1, mv_ref2, 
                    chor, cver, level, info );
    }

  free( ymvx0 );
  free( ymvy0 );
  free( cmvx0 );
  free( cmvy0 );

  free( ymvx1 );
  free( ymvy1 );
  free( cmvx1 );
  free( cmvy1 );

  free( ymvx2 );
  free( ymvy2 );
  free( cmvx2 );
  free( cmvy2 );

  free( ymvx0_int );
  free( ymvy0_int );
  free( cmvx0_int );
  free( cmvy0_int );

  free( ymvx1_int );
  free( ymvy1_int );
  free( cmvx1_int );
  free( cmvy1_int );

  free( ymvx2_int );
  free( ymvy2_int );
  free( cmvx2_int );
  free( cmvy2_int );

}


void
temporal_synthesis_with_OBMC( YUVimage_ptr fr0, YUVimage_ptr fr1, YUVimage H0,
                    YUVimage L1, YUVimage H1, YUVimage frp,
                    vector_ptr fmv0, vector_ptr fmv1, vector_ptr fmv2,
                    vector_ptr mv_ref0, vector_ptr mv_ref1, 
                    vector_ptr mv_ref2, int level, videoinfo info,
					ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
     /* fr1--current, fr0--reference */
{
  int yhor, yver, chor, cver;
  float *ymvx0, *ymvy0, *ymvx1, *ymvy1, *ymvx2, *ymvy2;
  float *cmvx0, *cmvy0, *cmvx1, *cmvy1, *cmvx2, *cmvy2;
  float *ymvx0_int, *ymvy0_int, *ymvx1_int, *ymvy1_int, *ymvx2_int, *ymvy2_int;
  float *cmvx0_int, *cmvy0_int, *cmvx1_int, *cmvy1_int, *cmvx2_int, *cmvy2_int;

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;
   
  blockmv2pixelmv( fmv0, &ymvx0, &ymvy0, &cmvx0, &cmvy0, CONNECTED, info, level );
  blockmv2pixelmv( fmv1, &ymvx1, &ymvy1, &cmvx1, &cmvy1, CONNECTED, info, level );
  blockmv2pixelmv( fmv2, &ymvx2, &ymvy2, &cmvx2, &cmvy2, CONNECTED, info, level );
  blockmv2pixelmv( fmv0, &ymvx0_int, &ymvy0_int, &cmvx0_int, &cmvy0_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv1, &ymvx1_int, &ymvy1_int, &cmvx1_int, &cmvy1_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv2, &ymvx2_int, &ymvy2_int, &cmvx2_int, &cmvy2_int, PREDICTED, info, level );

  mc_synthesis_with_OBMC( fr0->Y, fr1->Y, H0.Y, L1.Y, H1.Y, frp.Y, 
                ymvx0, ymvy0, ymvx1, ymvy1, ymvx2, ymvy2, 
                ymvx0_int, ymvy0_int, ymvx1_int, ymvy1_int,
                ymvx2_int, ymvy2_int, mv_ref0, mv_ref1, mv_ref2, 
                yhor, yver, level, info, frameMEinfo, varblkarray, 0, 0);
    
  if( info.cwidth && info.cheight ) 
    {
      mc_synthesis_with_OBMC( fr0->U, fr1->U, H0.U, L1.U, H1.U, frp.U, 
                    cmvx0, cmvy0, cmvx1, cmvy1, cmvx2, cmvy2,
                    cmvx0_int, cmvy0_int, cmvx1_int, cmvy1_int, 
                    cmvx2_int, cmvy2_int, mv_ref0, mv_ref1, mv_ref2, 
                    chor, cver, level, info, frameMEinfo, varblkarray, 1, 1);
      mc_synthesis_with_OBMC( fr0->V, fr1->V, H0.V, L1.V, H1.V, frp.V, 
                    cmvx0, cmvy0, cmvx1, cmvy1, cmvx2, cmvy2, 
                    cmvx0_int, cmvy0_int, cmvx1_int, cmvy1_int, 
                    cmvx2_int, cmvy2_int, mv_ref0, mv_ref1, mv_ref2, 
                    chor, cver, level, info, frameMEinfo, varblkarray, 1, 0);
    }

  free( ymvx0 );
  free( ymvy0 );
  free( cmvx0 );
  free( cmvy0 );

  free( ymvx1 );
  free( ymvy1 );
  free( cmvx1 );
  free( cmvy1 );

  free( ymvx2 );
  free( ymvy2 );
  free( cmvx2 );
  free( cmvy2 );

  free( ymvx0_int );
  free( ymvy0_int );
  free( cmvx0_int );
  free( cmvy0_int );

  free( ymvx1_int );
  free( ymvy1_int );
  free( cmvx1_int );
  free( cmvy1_int );

  free( ymvx2_int );
  free( ymvy2_int );
  free( cmvx2_int );
  free( cmvy2_int );

}



void
directional_iblock_synthesis_with_OBMC( YUVimage_ptr fr0, YUVimage_ptr fr1, YUVimage H0,
                    YUVimage L1, YUVimage H1, YUVimage frp,
                    vector_ptr fmv0, vector_ptr fmv1, vector_ptr fmv2,
                    vector_ptr mv_ref0, vector_ptr mv_ref1, 
                    vector_ptr mv_ref2, int t_level, videoinfo info,
					ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
     /* fr1--current, fr0--reference */
{
	int x, y, X, Y, yhor, yver;
	int xnum, ynum, xblk, yblk;
	vector_ptr fmv; 
	int s_level;

	// left scene change and right scene change, this is an intra-frame
	if( mv_ref0 == NULL && mv_ref1 == NULL ) return;  

	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));

	yhor = info.ywidth  << s_level;
	yver = info.yheight << s_level;
	xnum = info.xnum[t_level];
	ynum = info.ynum[t_level];
	xblk = info.xblk[t_level];
	yblk = info.yblk[t_level];

	if (mv_ref0 != NULL ) 
		fmv = fmv0;
	else if ( mv_ref1 != NULL )
		fmv = fmv1; 
	else
		assert(0); 

	for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
		for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
		// fmv is for current high temporal frame H0 ( non-NULL motion vector quad-tree ) 
		// fmv[Y * xnum + X] is current Macro-block 
		rec_directional_iblock_synthesis_with_OBMC( &fmv[Y * xnum + X], x, y, xblk, yblk, yhor, yver, 
										  fr0, fr1, H0, L1,  H1, frp,
										  fmv0, fmv1, fmv2,  mv_ref0, mv_ref1, 
										  mv_ref2, t_level, info, fmv,
										  frameMEinfo, varblkarray);
		}
	}

}


// by Yongjun Wu
// In decoder mv_ref0, mv_ref1 and mv_ref2 (whether it's NULL ) can exactly indicate the scene change 
// information and the existence of motion vectors in left and/or rigth directions.
// However, in encoder mv_ref1, mv_ref2 and mv_ref3 can not exactly indicate the information. 
// Sometimes scene is changed but mv_refi is still not NULL.
// hence we have to use the correpsonding scene_change[][] to indicate the information.

/*
 *                               synscheme 3
 */
void
synscheme3( int curr, videoinfo info, enum FLAG first_GOP,
            enum FLAG Level_change, int remaining_frs, int simul_dec, int theo_dec )
{
  int i, j, k, tPyrLev, t_level, dist, yhor, yver, chor, cver;
  int bigGOP, GOPsz, half_bigGOP, half_GOPsz, GOPIndex;  
  int ind_mv_GOP[20], ind_pyrFrs_bigGOP[20], ind_pyrFrs_GOP[20];
  int eff_GOPsz[20]; // eff_GOPsz in level i
  YUVimage dumb;

  YUVimage middle;

  int total_varblk;  // for MCTF with OBMC 

  frame_alloc(&dumb, info); 

  frame_alloc(&middle, info); 
  
  /*****************/
  /* Base Settings */
  /*****************/

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;
  tPyrLev     = info.tPyrLev; 
  t_level     = 0;
  bigGOP      = info.bigGOP;       // 31, 15, 7, 3
  half_bigGOP = info.bigGOP / 2;   // 15,  7, 3, 1
  GOPsz       = info.GOPsz;        // 16,  8, 4, 2
  half_GOPsz  = info.GOPsz  / 2;   //  8,  4, 2, 1
  GOPIndex    = curr / GOPsz;

  // determine effective GOP size in level i
  for( i = 0; i < info.tPyrLev; i++ ) {
    if ( Level_change == YES )
      eff_GOPsz[i] = (int) ceil ((double)info.eff_GOPsz / (pow (2, i)));  
    else
      eff_GOPsz[i] = GOPsz;  
    GOPsz /= 2;
  }
  GOPsz = info.GOPsz;

  dist = 1 << ( info.tPyrLev - 1 );
  
  // note: indices for bigGop are not the same as in analysis!
  // MV indices for transmission, e.g. for tPyrLev = 4:
  // Level0: 18-34; Level1: 9-17; Level2: 4-8; Level3: 1-3
  // the first set of motion vectors are always NULL
  // such as 1 in level3, 4 in level 2, ...
  ind_mv_GOP       [tPyrLev] = -1; // 1, 4, 9 , 18, ...

  // Indices of BigGOP-Array
  ind_pyrFrs_bigGOP[tPyrLev] =  0; // 0, 1, 3, 6, 11, ...
  ind_pyrFrs_GOP   [tPyrLev] =  1; // 0, 1, 2, 4, 8 , ...

  // calculation of base indices for MV-sets and frames
  for( i = tPyrLev - 1; i >= 0; i-- ){
    ind_mv_GOP[i]        = (int)(pow (2, tPyrLev-i-1)) + ind_mv_GOP[i+1] + 1;
    ind_pyrFrs_bigGOP[i] = (int)(pow (2, tPyrLev-i-2)) + ind_pyrFrs_bigGOP[i+1] + 1; 
    ind_pyrFrs_GOP[i]    = (int)(pow (2, tPyrLev-i-1));
  }

//  printf("Enter scene change\n");
  // free MV in case of scene changes
  for( i = 0; i < tPyrLev; i++ ){
    for( j = 0; j <= GOPsz; j++ ){ // effective GOPsz!
      if( dec_scene_change[i][j] == YES ){
        free_vector( dec_yfmv[j + ind_mv_GOP[i]], info );
        mv_ref[j + ind_mv_GOP[i]] = NULL;
      } else {
        mv_ref[j + ind_mv_GOP[i]] = dec_yfmv[j + ind_mv_GOP[i]];
      }
    }
    GOPsz /= 2;
  }
  GOPsz = info.GOPsz;

//  printf("Out of scene change\n");
  
  // not really necessary?
  free_vector( dec_yfmv[1], info );  
  mv_ref[1] = NULL; 

  /**********************/
  /* Temporal Filtering */
  /**********************/

  if( first_GOP == YES ){  // *** first_GOP = YES *** // 合成一个gop
    // decode first frame (..., LLL, LL, L, A)
    GOPsz = 1; 
    for( i = tPyrLev - 1; i >= t_level; i-- ){ // 逐级的合成
      // for synthesis all MV's must be given 为了synthesis，所有的MV必须要给
      // ==> initialize missing MV's and set mode = PREDICTED 初始化丢失的MV，并且设置mode = PREDICTED，mvx & mvy = (float)HUGE_VAL
      // ==> mvx & mvy = (float)HUGE_VAL in temporal_synthesis
      if( i == tPyrLev - 1 ){ // *** last level ***  最下面一层，先进行这里
#ifdef   MCTF_WITH_OBMC 
		// no need for MCTF with OBMC here, since H0 is dumb 这里不需要带OBMC的MCTF，因为H0是哑的
		// no need for directional IBLOCK reconstruction since H0 is dumb 这里不需要directional IBLOCK reconstruction，因为H0是哑的
		printf("i == tPyrLev - 1  level and frame: i = %d, j = %d\n",i,j);
        temporal_synthesis_with_OBMC( 
			                &dumb,          // fr0 (B1), dumb 
                            &dec_pyrTemp[i][0], // fr1 (A1)
                            dumb,           // H0, dumb 
                            dec_pyrFrs[0],      // L1
                            dec_pyrFrs[1],      // H1
                            dumb,           // frp (A0), dumb
                            dec_yfmv[1],        // mv0, free
                            dec_yfmv[1],        // mv1, free
                            dec_yfmv[2],        // mv2
                            mv_ref[1],      // NULL
                            mv_ref[1],      // NULL
                            mv_ref[2],      // ref_mv2
                            i, info, NULL, NULL );

#else
//		printf("i == tPyrLev - 1  level and frame: i = %d, j = %d\n",i,j);
        temporal_synthesis( &dumb,          // fr0 (B1), dumb 
                            &dec_pyrTemp[i][0], // fr1 (A1)
                            dumb,           // H0, dumb 
                            dec_pyrFrs[0],      // L1
                            dec_pyrFrs[1],      // H1
                            dumb,           // frp (A0), dumb
                            dec_yfmv[1],        // mv0, free
                            dec_yfmv[1],        // mv1, free
                            dec_yfmv[2],        // mv2
                            mv_ref[1],      // NULL
                            mv_ref[1],      // NULL
                            mv_ref[2],      // ref_mv2
                            i, info );
	   
#endif
      } 
	  else { // *** next levels *** 
#ifdef   MCTF_WITH_OBMC
		// no need for MCTF with OBMC here, since H0 is dumb
		// no need for directional IBLOCK reconstruction since H0 is dumb

		printf("next levels, level and frame: i = %d, j = %d\n",i,j);
        temporal_synthesis_with_OBMC( 
			                &dumb,             // fr0 (B1), dumb  
                            &dec_pyrTemp[i][0],    // fr1 (A1)
                            dumb,              // H0, dumb 
                            dec_pyrTemp[i + 1][0], // L1 
                            dec_pyrFrs[ind_pyrFrs_GOP[i]], // H1
                            dumb,              // frp (A0), dumb
                            dec_yfmv[1],           // mv0, free
                            dec_yfmv[1],           // mv1, free
                            dec_yfmv[ind_mv_GOP[i] + 1],   // mv2 
                            mv_ref[1],         // NULL
                            mv_ref[1],         // NULL
                            mv_ref[ind_mv_GOP[i] + 1], // ref_mv2
                            i, info, NULL, NULL );  

#else
//		printf("next levels, level and frame: i = %d, j = %d\n",i,j);
        temporal_synthesis( &dumb,             // fr0 (B1), dumb  
                            &dec_pyrTemp[i][0],    // fr1 (A1)
                            dumb,              // H0, dumb 
                            dec_pyrTemp[i + 1][0], // L1 
                            dec_pyrFrs[ind_pyrFrs_GOP[i]], // H1
                            dumb,              // frp (A0), dumb
                            dec_yfmv[1],           // mv0, free
                            dec_yfmv[1],           // mv1, free
                            dec_yfmv[ind_mv_GOP[i] + 1],   // mv2 
                            mv_ref[1],         // NULL
                            mv_ref[1],         // NULL
                            mv_ref[ind_mv_GOP[i] + 1], // ref_mv2
                            i, info );  
		
#endif
      }
      GOPsz *= 2;  
    }
  }
  else { // *** first_GOP = NO ***
     
    // copy first frame / MVF of new GOP to last position of preceeding (big)GOP
    // in order to be able to decode the preceeding GOP completely 为了完全解码前面一个gop，复制新gop的第一帧/MVF到先前gop的最后一个位置
    if ( remaining_frs != 0 ){  
      GOPsz = 2;
      for( i = tPyrLev - 1; i >= t_level; i-- ){
       
        // bigGOP-Array
        copyframe( &dec_pyrFrs[ind_pyrFrs_GOP[i]], &dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + GOPsz / 2 ], info );
        free_vector( dec_yfmv_bigGOP[ind_mv_GOP[i] + GOPsz], info );
        mv_copy( dec_yfmv[ind_mv_GOP[i] + 1], dec_yfmv_bigGOP[ind_mv_GOP[i] + GOPsz], info );
        dec_mv_ref_bigGOP[ind_mv_GOP[i] + GOPsz] = mv_ref[ind_mv_GOP[i] + 1] ? dec_yfmv_bigGOP[ind_mv_GOP[i] + GOPsz] : NULL;
        
        GOPsz *= 2;
      } 
      copyframe( &dec_pyrFrs[0], &dec_pyrFrs_bigGOP[0], info ); // lowpass copy  
    }

    // decode GOP

    // *** last level ***

    for( k = 1; k <= 3; k++ ){
      assert ((dec_mv_ref_bigGOP[k] == dec_yfmv_bigGOP[k]) || (dec_mv_ref_bigGOP[k] == NULL));
    }

#ifdef   MCTF_WITH_OBMC 
	if (tPyrLev - 1 >= info.t_level) {
		get_mv_side_information(info, 
								dec_yfmv_bigGOP[1],            // mv0
								dec_yfmv_bigGOP[2],            // mv1
								dec_mv_ref_bigGOP[1],          // ref_mv0
								dec_mv_ref_bigGOP[2],          // ref_mv1
  								frameMEinfo,  varblkarray, 
								&total_varblk, tPyrLev - 1, 0,
								NO, NO, 0);
		mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, tPyrLev - 1, 0);
	}

	printf("i == tPyrLev - 1 , next GOP, level and frame: i = %d, j = %d\n",i,j);
	temporal_synthesis_with_OBMC( 
		                &dec_pyrTemp[tPyrLev - 1][0],  // fr0 (B1)
                        &dec_pyrTemp[tPyrLev - 1][1],  // fr1 (A1)
                        dec_pyrFrs_bigGOP[1],          // H0
                        dec_pyrFrs_bigGOP[0],          // L1 
                        dec_pyrFrs_bigGOP[2],          // H1
                        dec_pyrFrs_first[tPyrLev - 1], // frp (A0)
                        dec_yfmv_bigGOP[1],            // mv0
                        dec_yfmv_bigGOP[2],            // mv1
                        dec_yfmv_bigGOP[3],            // mv2
                        dec_mv_ref_bigGOP[1],          // ref_mv0
                        dec_mv_ref_bigGOP[2],          // ref_mv1
                        dec_mv_ref_bigGOP[3],          // ref_mv2
                        tPyrLev - 1, info,
						frameMEinfo, varblkarray);    

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // synthesis IBLOCK 
    directional_iblock_synthesis_with_OBMC( 
		                &dec_pyrTemp[tPyrLev - 1][0],  // fr0 (B1)
                        &dec_pyrTemp[tPyrLev - 1][1],  // fr1 (A1)
                        dec_pyrFrs_bigGOP[1],          // H0
                        dec_pyrFrs_bigGOP[0],          // L1 
                        dec_pyrFrs_bigGOP[2],          // H1
                        dec_pyrFrs_first[tPyrLev - 1], // frp (A0)
                        dec_yfmv_bigGOP[1],            // mv0
                        dec_yfmv_bigGOP[2],            // mv1
                        dec_yfmv_bigGOP[3],            // mv2
                        dec_mv_ref_bigGOP[1],          // ref_mv0
                        dec_mv_ref_bigGOP[2],          // ref_mv1
                        dec_mv_ref_bigGOP[3],          // ref_mv2
                        tPyrLev - 1, info,
						frameMEinfo, varblkarray);  
#endif 


#else
	if (tPyrLev - 1 >= info.t_level) {
		get_mv_side_information(info, 
								dec_yfmv_bigGOP[1],            // mv0
								dec_yfmv_bigGOP[2],            // mv1
								dec_mv_ref_bigGOP[1],          // ref_mv0
								dec_mv_ref_bigGOP[2],          // ref_mv1
  								frameMEinfo,  varblkarray, 
								&total_varblk, tPyrLev - 1, 0,
								NO, NO, 0);
		mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, tPyrLev - 1, 0);
	}

//    printf("i == tPyrLev - 1 , next GOP, level and frame: i = %d, j = %d\n",i,j);
    temporal_synthesis( &dec_pyrTemp[tPyrLev - 1][0],  // fr0 (B1)
                        &dec_pyrTemp[tPyrLev - 1][1],  // fr1 (A1)
                        dec_pyrFrs_bigGOP[1],          // H0
                        dec_pyrFrs_bigGOP[0],          // L1 
                        dec_pyrFrs_bigGOP[2],          // H1
                        dec_pyrFrs_first[tPyrLev - 1], // frp (A0)
                        dec_yfmv_bigGOP[1],            // mv0
                        dec_yfmv_bigGOP[2],            // mv1
                        dec_yfmv_bigGOP[3],            // mv2
                        dec_mv_ref_bigGOP[1],          // ref_mv0
                        dec_mv_ref_bigGOP[2],          // ref_mv1
                        dec_mv_ref_bigGOP[3],          // ref_mv2
                        tPyrLev - 1, info );    

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // synthesis IBLOCK 

    directional_iblock_synthesis_with_OBMC( 
		                &dec_pyrTemp[tPyrLev - 1][0],  // fr0 (B1)
                        &dec_pyrTemp[tPyrLev - 1][1],  // fr1 (A1)
                        dec_pyrFrs_bigGOP[1],          // H0
                        dec_pyrFrs_bigGOP[0],          // L1 
                        dec_pyrFrs_bigGOP[2],          // H1
                        dec_pyrFrs_first[tPyrLev - 1], // frp (A0)
                        dec_yfmv_bigGOP[1],            // mv0
                        dec_yfmv_bigGOP[2],            // mv1
                        dec_yfmv_bigGOP[3],            // mv2
                        dec_mv_ref_bigGOP[1],          // ref_mv0
                        dec_mv_ref_bigGOP[2],          // ref_mv1
                        dec_mv_ref_bigGOP[3],          // ref_mv2
                        tPyrLev - 1, info,
						frameMEinfo, varblkarray);  
#endif 

#endif 

    // *** next levels ***
    GOPsz = 2; 
    for( i = tPyrLev - 2; i >= t_level; i-- ){
      for( j = 0; j < GOPsz; j++ ){
                
        for( k = (ind_mv_GOP[i] + 2 * j); k <= (ind_mv_GOP[i] + 2 * j + 2); k++ ){
          assert ((dec_mv_ref_bigGOP[k] == dec_yfmv_bigGOP[k]) || (dec_mv_ref_bigGOP[k] == NULL));
        }
        if( j == 0 ){
#ifdef   MCTF_WITH_OBMC 
		  if (i >= info.t_level) {
			  get_mv_side_information(info, 
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],         // mv0
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],     // mv1
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
  						   	  		  frameMEinfo,  varblkarray, 
									  &total_varblk, i, j,
									  NO, NO, 0);

			  mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);
		  }

		  printf("Next level, next GOP, level and frame: i = %d, j = %d\n",i,j);
          temporal_synthesis_with_OBMC( 
			                  &dec_pyrTemp[i][2 * j],                         // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                     // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],    // H0
                              dec_pyrTemp[i + 1][j],                          // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1],// H1
                              dec_pyrFrs_first[i],                            // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],         // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],     // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],     // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],   // ref_mv2
                              i, info,
							  frameMEinfo, varblkarray);   
#ifdef DIRECTIONAL_IBLOCK_EMPLOYED
          directional_iblock_synthesis_with_OBMC( 
			                  &dec_pyrTemp[i][2 * j],                         // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                     // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],    // H0
                              dec_pyrTemp[i + 1][j],                          // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1],// H1
                              dec_pyrFrs_first[i],                            // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],         // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],     // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],     // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],   // ref_mv2
                              i, info,
							  frameMEinfo, varblkarray);  
#endif 

		  if(i == EXTRACT_TMP && EXT){
			  if( dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j] != NULL || dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1] != NULL ){
				wcopyframe(&dec_pyrTemp[i][2 * j], &middle, 1/pow(LPW1[1],EXTRACT_TMP), info);
				write_frame(middle,info,info.jp2k_decname,curr/2 + 2*j, info.format);
			  }else{
//				wcopyframe(&dec_pyrTemp[i][2 * j], &middle, 1, info);
//			    assert(0);
			  }

			  if( dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1] != NULL || dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2] != NULL ){
				wcopyframe(&dec_pyrTemp[i][2 * j + 1], &middle, 1/pow(LPW1[1],EXTRACT_TMP), info);
				write_frame(middle,info,info.jp2k_decname,curr/2 + 2*j + 1, info.format);
			  }else{
//				wcopyframe(&dec_pyrTemp[i][2 * j + 1], &middle, 1, info);
//				assert(0);
			  }

		  }

#else
		  if (i >= info.t_level) {
			  get_mv_side_information(info, 
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],         // mv0
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],     // mv1
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
  						   	  		  frameMEinfo,  varblkarray, 
									  &total_varblk, i, j,
									  NO, NO, 0);

			  mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);
		  }

//		  printf("Next level, next GOP, level and frame: i = %d, j = %d\n",i,j);
          temporal_synthesis( &dec_pyrTemp[i][2 * j],                         // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                     // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],    // H0
                              dec_pyrTemp[i + 1][j],                          // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1],// H1
                              dec_pyrFrs_first[i],                            // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],         // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],     // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],     // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],   // ref_mv2
                              i, info );

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED

          directional_iblock_synthesis_with_OBMC( 
			                  &dec_pyrTemp[i][2 * j],                         // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                     // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],    // H0
                              dec_pyrTemp[i + 1][j],                          // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1],// H1
                              dec_pyrFrs_first[i],                            // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],         // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],     // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],     // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],   // ref_mv2
                              i, info,
							  frameMEinfo, varblkarray);  
#endif 

#endif 
        } else {

#ifdef MCTF_WITH_OBMC
		  if (i >= info.t_level) {
			  get_mv_side_information(info, 
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],          // mv0
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],      // mv1
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
  						   	  		  frameMEinfo,  varblkarray, 
									  &total_varblk, i, j,
									  NO, NO, 0);

			  mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);
		  }

		  printf("level and frame: i = %d, j = %d\n",i,j);
          temporal_synthesis_with_OBMC( 
			                  &dec_pyrTemp[i][2 * j],						   // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                      // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],     // H0
                              dec_pyrTemp[i + 1][j],                           // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1], // H1
                              dec_pyrTemp[i][2 * j - 1], // statt pyrFrs_first // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],          // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],      // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],      // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],        // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],    // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],    // ref_mv2
                              i, info,
							  frameMEinfo, varblkarray); 

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED
          directional_iblock_synthesis_with_OBMC( 
			                  &dec_pyrTemp[i][2 * j],                          // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                      // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],     // H0
                              dec_pyrTemp[i + 1][j],                           // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1], // H1
                              dec_pyrTemp[i][2 * j - 1], // statt pyrFrs_first // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],          // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],      // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],      // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],        // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],    // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],    // ref_mv2
                              i, info,
							  frameMEinfo, varblkarray); 
#endif 

#else
		  if (i >= info.t_level) {
			  get_mv_side_information(info, 
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],          // mv0
									  dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],      // mv1
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],       // ref_mv0
									  dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],   // ref_mv1
  						   	  		  frameMEinfo,  varblkarray, 
									  &total_varblk, i, j,
									  NO, NO, 0);

			  mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);
		  }

//		  printf("level and frame: i = %d, j = %d\n",i,j);
          temporal_synthesis( &dec_pyrTemp[i][2 * j],						   // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                      // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],     // H0
                              dec_pyrTemp[i + 1][j],                           // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1], // H1
                              dec_pyrTemp[i][2 * j - 1], // statt pyrFrs_first // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],          // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],      // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],      // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],        // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],    // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],    // ref_mv2
                              i, info );

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED

          directional_iblock_synthesis_with_OBMC( 
			                  &dec_pyrTemp[i][2 * j],                          // fr0 (B1)
                              &dec_pyrTemp[i][2 * j + 1],                      // fr1 (A1)
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j],     // H0
                              dec_pyrTemp[i + 1][j],                           // L1 
                              dec_pyrFrs_bigGOP[ind_pyrFrs_bigGOP[i] + j + 1], // H1
                              dec_pyrTemp[i][2 * j - 1], // statt pyrFrs_first // frp (A0)
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j],          // mv0
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 1],      // mv1
                              dec_yfmv_bigGOP[ind_mv_GOP[i] + 2 * j + 2],      // mv2
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j],        // ref_mv0
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 1],    // ref_mv1
                              dec_mv_ref_bigGOP[ind_mv_GOP[i] + 2 * j + 2],    // ref_mv2
                              i, info,
							  frameMEinfo, varblkarray); 
#endif 

#endif
        }
      }
      GOPsz *= 2;  
    }
    //  last GOP: remaining_frs == 0 !
    
  }
  
  /********************/
  /* Save information 为下一个gop保存信息 */
  /*   for next GOP   */
  /********************/
  if ( remaining_frs != 0 && theo_dec == NO ){ // 最后不满一个gop也会进来，来保存信息，因为下一个gop不读fream的

    GOPsz = info.GOPsz >> t_level; 

    for( i = t_level; i < tPyrLev; i++ ){
      // save dec_pyrTemp[i][0] for next GOP  保存第一帧
      if ( first_GOP == YES )
        copyframe( &dec_pyrTemp[i][0], &dec_pyrFrs_first[i], info );  
      else { 
        copyframe( &dec_pyrTemp[i][GOPsz - 1], &dec_pyrFrs_first[i], info );
      }	  

      for( j = 0; j < GOPsz / 2; j++ ){ //eff_GOPsz[i]
        // copy frames 保存帧
        copyframe( &dec_pyrFrs[j + ind_pyrFrs_GOP[i]], &dec_pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]], info );  
      }
	  
      for( j = 0; j <= GOPsz; j++ ){
        if (j + ind_mv_GOP[i] != 1){ 
          // copy MVs 保存mv
          free_vector( dec_yfmv_bigGOP[j + ind_mv_GOP[i] - 1], info );
          mv_copy( dec_yfmv[j + ind_mv_GOP[i]], dec_yfmv_bigGOP[j + ind_mv_GOP[i] - 1], info );
          dec_mv_ref_bigGOP[j + ind_mv_GOP[i] - 1] = mv_ref[j + ind_mv_GOP[i]] ?
            dec_yfmv_bigGOP[j + ind_mv_GOP[i] - 1] : NULL;
        }
      }
      GOPsz /= 2;
    }
  }// if simul_dec

}


void
denoise_mctf_syn_ezbc( int curr, int GOP_counter, long int *total_bytes_past,
                       videoinfo info, enum FLAG first_GOP, enum FLAG Level_change,
                       int remaining_frs )
{
  int i, j, t_level, nfrs, dist;
  int filterType = DENOISE_FILTER;
  YUVimage tempfr, tempfr2;

  t_level = MY_MIN (info.t_level, info.tPyrLev);
  dist = 1 << ( t_level );
  
  if ( remaining_frs != 0 ){
    nfrs = info.GOPsz;
  } else {
    nfrs = info.eff_GOPsz;
  }
 
  if ( first_GOP == YES ) nfrs = 1;
  
  /*************/
  /* MV + EZBC */
  /*************/
  if ( remaining_frs != 0 ){
    decode_MV( total_bytes_past, info, GOP_counter, 0, 0 );
    ezbc3d_dec_GOP( Four_GOP, info, *total_bytes_past, GOP_counter, curr );
      
    if ( info.s_level == 0 ){
      if ( first_GOP == YES ){
        for( i = 1; i < ((YUV420 == 1) ? 2 : 4); i++ ) {
          for( j = 0; j < info.GOPsz; j++ ) {
            spatial_high[(i-1)][j] = Four_GOP[i * info.GOPsz + j];
          }
        }
      } else { // first_GOP == NO
        for( i = 1; i < ((YUV420 == 1) ? 2 : 4); i++ ) {
          copyframe( &Four_GOP[i * info.GOPsz], &Four_bigGOP[(i+1) * info.GOPsz - (1 << t_level)], info );      
          // printf("copyframe( &Four_GOP[%d], &Four_bigGOP[%d] );\n", i * info.GOPsz, (i+1) * info.GOPsz - (1 << t_level));      
           
          for( j = 0; j < info.GOPsz; j++ ) {
            spatial_high[(i-1)][j] = Four_bigGOP[i * info.GOPsz + j];
          }
        }
      }
    }
    
    for( j = 0; j < info.GOPsz; j++ ) {
      dec_pyrFrs[j] = Four_GOP[j];
    }

  } 
 
  /********/
  /* MCTF */
  /********/ 
  mc_syn( curr, total_bytes_past, info, first_GOP,
          Level_change, remaining_frs, 0, 0 );

  /*************/
  /* spatial   */
  /* synthesis */
  /*************/
  if( info.s_level == 0 ){ // Quarter == 0
     
    info.ywidth  *= 2;
    info.yheight *= 2;
    info.cwidth  *= 2;
    info.cheight *= 2;

    // for( i = 0; i < info.eff_GOPsz; i++ ) {
    for( i = 0; i < nfrs; i += dist ) {  
      frame_alloc( &tempfr, info );
      frame_alloc( &tempfr2, info );
      
      if ( first_GOP == YES ) dist = 1;

      if( YUV420 == 1 ) {
        spatial_syn( tempfr.Y, info.ywidth, info.yheight, dec_pyrTemp[0][i + dist - 1].Y,
                     spatial_high[0][i].Y, spatial_high[0][i].U,
                     spatial_high[0][i].V, filterType );

        for( j = 0; j < info.cwidth * info.cheight / 4; j++ ) {
          tempfr.U[j] = dec_pyrTemp[0][i + dist - 1].U[j] * 4;
          tempfr.V[j] = dec_pyrTemp[0][i + dist - 1].V[j] * 4;
        }
        read_frame( &tempfr2, info, info.inname, curr + i, info.format );

        f444_420( info, &tempfr2 );
        info.cwidth  /= 2;  
        info.cheight /= 2;
        calsnr_frame( &tempfr, &tempfr2, curr + i, info );

        f420_444( info, &tempfr );
        info.cwidth  *= 2;
        info.cheight *= 2;

      } else { // YUV420 == 0
        spatial_syn_frame( &tempfr, info, &dec_pyrTemp[0][i + dist - 1], &spatial_high[0][i],
                           &spatial_high[1][i], &spatial_high[2][i], filterType );
        // printf("spatial_syn_frame( &dec_pyrTemp[0][%d], &spatial_high[0][%d], &spatial_high[1][%d], &spatial_high[2][%d]);\n",
        //       i >> t_level, i, i, i);
      }

      if ( first_GOP == YES )
        write_frame( tempfr, info, info.decname, curr, info.format );
      else{
        if ( (curr - info.GOPsz + i + dist) <= info.last ){
          write_frame( tempfr, info, info.decname, curr - info.GOPsz + i + dist, info.format );
          // printf(" write_frame( tempfr, curr - info.GOPsz + i + dist %d );\n",  curr - info.GOPsz + i + dist );
        }
      }
                   
      free_frame( tempfr );
      free_frame( tempfr2 );
    }

    info.ywidth  /= 2;
    info.yheight /= 2;
    info.cwidth  /= 2;
    info.cheight /= 2;

  } else { // Quarter == 1
    // printf( "denoise_mctf_syn_ezbc: only output spatial LL bands\n" );
    
    for( i = 0; i < nfrs; i += dist ) {  
      if ( filterType == 7 ){ 
        // scale down frame (correct scaling of the spatial transform) 
        wcopyframe( &dec_pyrTemp[0][i], &dec_pyrTemp[0][i], ( float )pow( 2, -1 ), info );
      }
      
      if ( first_GOP == YES )
        write_frame( dec_pyrTemp[0][0], info, info.decname, curr, info.format );
      else {
        if ( (curr - info.GOPsz + i + dist) <= info.last ){
          write_frame( dec_pyrTemp[0][i + dist - 1], info, info.decname, curr - info.GOPsz + i + dist, info.format );
          // printf(" write_frame( dec_pyrTemp[0][%d], curr - info.GOPsz + i + dist %d );\n", i, curr - info.GOPsz + i + dist );
        }
      }
    }
    
  }

  /********************/
  /* Save information */
  /*   for next GOP   */
  /********************/
  if ( remaining_frs != 0 ){
    for( i = 1; i < ((YUV420 == 1) ? 2 : 4); i++ ) {
      for( j = 1 << t_level; j < info.eff_GOPsz; j+=dist ) {
        copyframe( &Four_GOP[i * info.GOPsz + j], &Four_bigGOP[i * info.GOPsz + j - (1 << t_level)], info );      
        // printf("-> copyframe( &Four_GOP[%d], &Four_bigGOP[%d] );\n", i * info.GOPsz + j, i * info.GOPsz + j - (1 << t_level)); 
      }
    }
  }
  
}


/*
 *                              mc_syn  mc综合网络
 */
void
mc_syn( int curr, long int *total_bytes_past, videoinfo info,
        enum FLAG first_GOP, enum FLAG Level_change, int remaining_frs, int simul_dec, int theo_dec )
{
 
  // GOP header and motion fields
  if( info.tPyrLev >= 1 ) {
  
    synscheme3( curr, info, first_GOP, Level_change, remaining_frs, simul_dec, theo_dec );
    
    free_mvs( dec_yfmv, info );
  }
  else { // necessary??? KH, 2003-12-12
    copyframe( &dec_pyrFrs[0], &dec_pyrTemp[0][0], info );
  }
}


void
mctf_syn_ezbc( int curr, int GOP_counter, long int *total_bytes_past,
               videoinfo info, enum FLAG first_GOP, enum FLAG Level_change,
               int remaining_frs, int simul_dec, int theo_dec )
{
  int t_level, nfrs, dist;
  int i,j;
  float small_y=0, small_u=0, small_v=0, large_y=0, large_u=0, large_v=0;

  /*************/
  /* MV + EZBC */
  /*************/
  if (remaining_frs != 0) { // 如果是最后一帧，这次就不读了
	  decode_MV(total_bytes_past, info, GOP_counter, simul_dec, theo_dec); // 解码MV
#ifdef CNN_wavelet
  {  // 以文本形式读入
	//std::cout << "open file" << info.hpcbindata << std::endl;
	//  std::ifstream myfile(info.hpcbindata);
	//  //std::ifstream myfile("BasketballDrill_832x480_50FTFS_GOP1_000.data");
	//  if (myfile.is_open() == NULL)
	//  {
	//	  std::cout << "file open failed!";
	//	  exit(0);
	//  }
	//  float temp = 0;
	//  for (int pic = 0; pic < info.GOPsz; pic++)
	//  {
	//	  for (int i = 0; i < info.yheight; i++)
	//	  {
	//		  for (int j = 0; j < info.ywidth; j++)
	//		  {
	//			  myfile >> temp;
	//			  dec_pyrFrs[pic].Y[i*info.ywidth + j] = temp;
	//		  }
	//	  }
	//  }
	//  for (int pic = 0; pic < info.GOPsz; pic++)
	//  {
	//	  for (int i = 0; i < info.cheight; i++)
	//	  {
	//		  for (int j = 0; j < info.cwidth; j++)
	//		  {
	//			  myfile >> temp;
	//			  dec_pyrFrs[pic].U[i*info.cwidth + j] = temp;
	//		  }
	//	  }
	//  }
	//  for (int pic = 0; pic < info.GOPsz; pic++)
	//  {
	//	  for (int i = 0; i < info.cheight; i++)
	//	  {
	//		  for (int j = 0; j < info.cwidth; j++)
	//		  {
	//			  myfile >> temp;
	//			  dec_pyrFrs[pic].V[i*info.cwidth + j] = temp;
	//		  }
	//	  }
	//  }
	//  myfile.close();
  }
 {// 读二进制文件
	  //char name[200] = "BasketballDrill_832x480_50_FTFS_GOP1_000.bindata";
	  //std::cout << "open file" << name << std::endl;
	  //std::ifstream myfile(name, std::ios::binary);
	  //if (myfile.is_open() == NULL)
	  //{
		 // std::cout << "open file failed" << std::endl;
		 // exit(0);
	  //}
	  //short int temp = 0;
	  //myfile.seekg(sizeof(short int) * info.cwidth * info.cheight * 6 * GOP_counter * info.GOPsz);
	  //std::cout << "GOP_counter" << GOP_counter << std::endl;
	  //for (int pic = 0; pic < info.GOPsz; pic++)
	  //{
		 // for (int i = 0; i < info.yheight; i++)
		 // {
			//  for (int j = 0; j < info.ywidth; j++)
			//  {
			//	  myfile.read((char*)(&temp), sizeof(short int));
			//	  dec_pyrFrs[pic].Y[i*info.ywidth + j] = temp;
			//  }
		 // }
	  //}
	  //for (int pic = 0; pic < info.GOPsz; pic++)
	  //{
		 // for (int i = 0; i < info.cheight; i++)
		 // {
			//  for (int j = 0; j < info.cwidth; j++)
			//  {
			//	  myfile.read((char*)(&temp), sizeof(short int));
			//	  dec_pyrFrs[pic].U[i*info.cwidth + j] = temp;
			//  }
		 // }
	  //}
	  //for (int pic = 0; pic < info.GOPsz; pic++)
	  //{
		 // for (int i = 0; i < info.cheight; i++)
		 // {
			//  for (int j = 0; j < info.cwidth; j++)
			//  {
			//	  myfile.read((char*)(&temp), sizeof(short int));
			//	  dec_pyrFrs[pic].V[i*info.cwidth + j] = temp;
			//  }
		 // }
	  //}
	  //myfile.close();
 }
	  {// 读二进制文件  float形式读取 读Y
		 // //char name_y[200] = "hpcbindata/Tango2_3840x2160_60FTFS_GOP1_000.hpcbindata";
		 // //char name_u[200] = "RaceHorses_352x288_30FTFS_GOP1_000_post_u.bindata";
		 // //char name_v[200] = "RaceHorses_352x288_30FTFS_GOP1_000_post_v.bindata";
		 // std::cout << "open file" << info.hpcbindata << std::endl;
		 // std::ifstream myfile_y(info.hpcbindata, std::ios::binary);
		 // //std::ifstream myfile_u(name_u, std::ios::binary);
		 // //std::ifstream myfile_v(name_v, std::ios::binary);
		 // if (myfile_y.is_open() == NULL)
		 // {
			//  std::cout << "open file_y failed" << std::endl;
			//  exit(0);
		 // }
		 // //if (myfile_u.is_open() == NULL)
		 // //{
			// // std::cout << "open file_u failed" << std::endl;
			// // exit(0);
		 // //}
		 // //if (myfile_v.is_open() == NULL)
		 // //{
			// // std::cout << "open file_v failed" << std::endl;
			// // exit(0);
		 // //}
		 // myfile_y.seekg(sizeof(float) * info.ywidth * info.yheight * GOP_counter * info.GOPsz);
		 // //myfile_u.seekg(sizeof(float) * info.cwidth * info.cheight * GOP_counter * info.GOPsz);
		 // //myfile_v.seekg(sizeof(float) * info.cwidth * info.cheight * GOP_counter * info.GOPsz);
		 // std::cout << "GOP_counter" << GOP_counter << std::endl;
		 // float temp = 0;
		 // for (int pic = 0; pic < info.GOPsz; pic++)
		 // {
			//  for (int i = 0; i < info.yheight; i++)
			//  {
			//	  for (int j = 0; j < info.ywidth; j++)
			//	  {
			//		  myfile_y.read((char*)(&temp), sizeof(float));
			//		  dec_pyrFrs[pic].Y[i*info.ywidth + j] = temp;
			//	  }
			//  }
		 // }
		 // //for (int pic = 0; pic < info.GOPsz; pic++)
		 // //{
			// // for (int i = 0; i < info.cheight; i++)
			// // {
			//	//  for (int j = 0; j < info.cwidth; j++)
			//	//  {
			//	//	  myfile_u.read((char*)(&temp), sizeof(float));
			//	//	  dec_pyrFrs[pic].U[i*info.cwidth + j] = temp;
			//	//  }
			// // }
		 // //}
		 // //for (int pic = 0; pic < info.GOPsz; pic++) 
		 // //{
			// // for (int i = 0; i < info.cheight; i++)
			// // {
			//	//  for (int j = 0; j < info.cwidth; j++)
			//	//  {
			//	//	  myfile_v.read((char*)(&temp), sizeof(float));
			//	//	  dec_pyrFrs[pic].V[i*info.cwidth + j] = temp;
			//	//  }
			// // }
		 // //}
		 // myfile_y.close();
		 // //myfile_u.close();
		 // //myfile_v.close();
	  }
	  {// 读二进制文件  float形式读取 读Y
		  //std::cout << "open file" << info.hpcbindata << std::endl;
		  //std::ifstream myfile(info.hpcbindata, std::ios::binary);
		  //if (myfile.is_open() == NULL)
		  //{
			 // std::cout << "open file_y failed" << std::endl;
			 // exit(0);
		  //}
		  //
		  ////myfile.seekg(sizeof(short int) * info.cwidth * info.cheight * 6 * GOP_counter * info.GOPsz);
		  //myfile.seekg(sizeof(float) * info.cwidth * info.cheight * 4 * GOP_counter * info.GOPsz);
		  //std::cout << "GOP_counter" << GOP_counter << std::endl;
		  ////short int temp = 0;
		  //float temp = 0;
		  //for (int pic = 0; pic < info.GOPsz; pic++)
		  //{
			 // for (int i = 0; i < info.yheight; i++)
			 // {
				//  for (int j = 0; j < info.ywidth; j++)
				//  {
				//	  //myfile.read((char*)(&temp), sizeof(short int));
				//	  myfile.read((char*)(&temp), sizeof(float));
				//	  dec_pyrFrs[pic].Y[i*info.ywidth + j] = temp;
				//  }
			 // }
		  //}
		  //myfile.close();
	  }
	  {// 读二进制文件  float形式读取 读YUV
		  std::cout << "open file" << info.hpcbindata << std::endl;
		  std::ifstream myfile(info.hpcbindata, std::ios::binary);
		  if (myfile.is_open() == NULL)
		  {
			  std::cout << "open file_y failed" << std::endl;
			  exit(0);
		  }
		  
		  //myfile.seekg(sizeof(short int) * info.cwidth * info.cheight * 6 * GOP_counter * info.GOPsz);
		  myfile.seekg(sizeof(float) * info.cwidth * info.cheight * 6 * GOP_counter * info.GOPsz);
		  std::cout << "GOP_counter" << GOP_counter << std::endl;
		  //short int temp = 0;
		  float temp = 0;
		  for (int pic = 0; pic < info.GOPsz; pic++)
		  {
			  for (int i = 0; i < info.yheight; i++)
			  {
				  for (int j = 0; j < info.ywidth; j++)
				  {
					  //myfile.read((char*)(&temp), sizeof(short int));
					  myfile.read((char*)(&temp), sizeof(float));
					  //std::cout << temp << std::endl;
					  dec_pyrFrs[pic].Y[i*info.ywidth + j] = temp;
				  }
			  }
		  }
		  for (int pic = 0; pic < info.GOPsz; pic++)
		  {
			  for (int i = 0; i < info.cheight; i++)
			  {
				  for (int j = 0; j < info.cwidth; j++)
				  {
					  //myfile.read((char*)(&temp), sizeof(short int));
					  myfile.read((char*)(&temp), sizeof(float));
					  dec_pyrFrs[pic].U[i*info.cwidth + j] = temp;
				  }
			  }
		  }
		  for (int pic = 0; pic < info.GOPsz; pic++)
		  {
			  for (int i = 0; i < info.cheight; i++)
			  {
				  for (int j = 0; j < info.cwidth; j++)
				  {
					  //myfile.read((char*)(&temp), sizeof(short int));
					  myfile.read((char*)(&temp), sizeof(float));
					  dec_pyrFrs[pic].V[i*info.cwidth + j] = temp;
				  }
			  }
		  }
		  myfile.close();
	  }
 //   ezbc3d_dec_GOP( dec_pyrFrs, info, *total_bytes_past, GOP_counter, curr );// 获得重建帧
	//std::ofstream myfile("dec_Y_float_3.txt", std::ios::app);
	//for (int pic = 0; pic < info.GOPsz; pic++)
	//{
	//	for (int i = 0; i < info.yheight; i++)
	//	{
	//		for (int j = 0; j < info.ywidth; j++)
	//		{
	//			float temp = dec_pyrFrs[pic].Y[i*info.ywidth + j];
	//			myfile << temp << " ";
	//		}
	//		myfile << std::endl;
	//	}
	//}
	//myfile.close();
#else
    ezbc3d_dec_GOP( dec_pyrFrs, info, *total_bytes_past, GOP_counter, curr );// 获得重建帧

	//{// 写文本文件
	//	char name[256], data_file_name[256], start[10];
	//	strncpy(name, info.bitname, strlen(info.bitname) - 4);
	//	name[strlen(info.bitname) - 4] = '_';
	//	snprintf(start, 9, "%03d", info.start);
	//	for (int i = 0; i < 3; i++)
	//	{
	//		name[strlen(info.bitname) - 3 + i] = start[i];
	//	}
	//	name[strlen(info.bitname)] = '\0';               //
	//	sprintf(data_file_name, "%s.decdata", name);
	//	std::cout << data_file_name << std::endl;
	//	std::ofstream myfile(data_file_name, std::ios::app);
	//	for (int pic = 0; pic < info.GOPsz; pic++)
	//	{
	//		for (int i = 0; i < info.yheight; i++)
	//		{
	//			for (int j = 0; j < info.ywidth; j++)
	//			{
	//				//int temp = round(pyrFrs[pic].Y[i*info.ywidth + j]);
	//				float temp = (dec_pyrFrs[pic].Y[i*info.ywidth + j]);
	//				myfile << temp << " ";
	//			}
	//			myfile << std::endl;
	//		}
	//	}
	//	for (int pic = 0; pic < info.GOPsz; pic++)
	//	{
	//		for (int i = 0; i < info.cheight; i++)
	//		{
	//			for (int j = 0; j < info.cwidth; j++)
	//			{
	//				//int temp = round(pyrFrs[pic].U[i*info.cwidth + j]);
	//				float temp = (dec_pyrFrs[pic].U[i*info.cwidth + j]);
	//				myfile << temp << " ";
	//			}
	//			myfile << std::endl;
	//		}
	//	}
	//	for (int pic = 0; pic < info.GOPsz; pic++)
	//	{
	//		for (int i = 0; i < info.cheight; i++)
	//		{
	//			for (int j = 0; j < info.cwidth; j++)
	//			{
	//				//int temp = round(pyrFrs[pic].V[i*info.cwidth + j]);
	//				float temp = (dec_pyrFrs[pic].V[i*info.cwidth + j]);
	//				myfile << temp << " ";
	//			}
	//			myfile << std::endl;
	//		}
	//	}
	//	myfile.close();
	//}


#endif
  }

  /********/
  /* MCTF */
  /********/

  mc_syn( curr, total_bytes_past, info, first_GOP,
	  Level_change, remaining_frs, simul_dec, theo_dec );

  /****************/
  /* write frames */
  /****************/

  t_level = MY_MIN (info.t_level, info.tPyrLev);
  dist = 1 << ( t_level );
  
  if ( remaining_frs != 0 ){
    nfrs = info.GOPsz;
  } else {
    nfrs = info.eff_GOPsz;
  }
 
  if( first_GOP == YES ){
      write_frame( dec_pyrTemp[0][0], info, info.decname, curr, info.format );
      // printf(" write_frame( dec_pyrTemp[0][0], curr %d );\n", curr );
  } 
  else 
  {
    for( i = 0; i < nfrs; i += dist ) 
	{  
      if ( (curr - info.GOPsz + i + dist) <= info.last ){
        write_frame( dec_pyrTemp[0][i + dist - 1], info, info.decname, curr - info.GOPsz + i + dist, info.format );
        // printf(" write_frame( dec_pyrTemp[0][%d], curr - info.GOPsz + i + dist %d );\n", i, curr - info.GOPsz + i + dist );
      }
    }
  }
  
}

