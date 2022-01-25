#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#define EXTERN extern
#include "structN.h"
#include "coderN.h"
#include "basic.h"
#include "util_filtering.h"
#include "memoryN.h"
#include "general.h"
#include <iostream>

/*****************************************************************************/
/*                         init_dec()                                        */
/*****************************************************************************/
void
init_dec( videoinfo info )
{
  int i, j, /*xnum, ynum, xblk, yblk, hor, ver,*/ GOPsz;
  int num_mv_GOP, num_pyrFrs_bigGOP, rel_s_level;
  /*int   picture_rate, totalrate, illbit, pllbit, lhbit, hbit; */

//#ifdef   MCTF_WITH_OBMC
  int hor, ver, s_level;
  s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
  hor = info.ywidth<<s_level;
  ver = info.yheight<<s_level; 
  // allocate memory for frameMEinfo, varblkarray
  // pixel based motion vector side information, and side information for each block
  // they are always in full resolution ( original resolution )
  frameMEinfo  = (ImageMEinfo*) getarray(hor*ver, sizeof(ImageMEinfo), "frameMEinfo");
  varblkarray  = (Varblkarrayinfo*) getarray(hor/SMALLEST_BLKSIZE*ver/SMALLEST_BLKSIZE, 
		sizeof(Varblkarrayinfo), "varblkarray");
//#endif


  interpolate_filter(  );
  temporal_filter(  );
  rel_s_level = (info.s_level == 0 && info.denoise_flag == YES);

  printf( "---------------------------------------------------------------\n" );
  printf( " Bit file         : %s\n", info.bitname );
  printf( " Decoded sequence : %s\n", info.decname );
  printf( " Sequence length  : %03d frames ( %03d - %03d )\n", 
          info.last - info.start + 1, info.start, info.last );
  for (i = 0; i < info.tPyrLev; i++) {
    printf( " Motion Level %d   : HVSBM with Bx %d By %d ... Bx %d By %d and ", 
            i, info.xblk[i], info.yblk[i],
            info.xblk[i] >> (info.level[i] - 1), info.yblk[i] >> (info.level[i] - 1));
    
    if( info.subpel[i] == 4 )
      printf( "1/16-pixel accuracy\n" );
    else if( info.subpel[i] == 3 )
      printf( "1/8-pixel accuracy\n" );
    else if( info.subpel[i] == 2 )
      printf( "quarter-pixel accuracy\n" );
    else if( info.subpel[i] == 1 )
      printf( "half-pixel accuracy\n" );
    else
      printf( "full-pixel accuracy\n" );
  }
  printf( " Frame rate       : %0.2f frames/sec (t_level %d)\n",
          ((double)(info.framerate) / pow( (float)2, (int)info.t_level)), info.t_level );
  printf( " Frame size       : (Y) %d x %d (C) %d x %d (s_level %d)\n",
          info.ywidth << rel_s_level, info.yheight << rel_s_level,
          info.cwidth << rel_s_level, info.cheight << rel_s_level, 
          info.s_level );
  if ( info.denoise_flag == YES )
    printf( " Denoise          : YES\n");
#ifdef MY_SCALE
  printf( " MY_SCALE         : %d (%d extra bits)\n", MY_SCALE, (int)log2(MY_SCALE));
#endif
  printf( "---------------------------------------------------------------\n" );
  /****************** initialize the sequence information ****************/
  GOPsz = info.GOPsz;
 
  // same number of MVs for GOP and bigGOP !
  num_mv_GOP        = 1; // 4, 9, 18, 35, ...
  num_pyrFrs_bigGOP = 1; // 3, 6, 11, 20, ...
     
  // calculation of base indices for MV-sets and frames
  for( i = info.tPyrLev-1; i >= 0; i-- ){
    num_mv_GOP += info.GOPsz / (int)( pow( (float)2, (int)i)) + 1;
  }
  num_pyrFrs_bigGOP = info.GOPsz + info.tPyrLev + 1; // 为什么这样？？？？？？？

  printf("num_pyrFrs_bigGOP = %d, num_mv_GOP = %d\n", num_pyrFrs_bigGOP, num_mv_GOP);
  
  /***************** allocate the 1st layer motion vectors ***************/
  dec_yfmv_bigGOP =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "dec_yfmv_bigGOP_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ){
    dec_yfmv_bigGOP[j] = 
      ( vector_ptr ) getarray( info.maxMBnum, sizeof( vector ), "dec_yfmv_bigGOP" );
    alloc_vector( dec_yfmv_bigGOP[j], info );
  }
  
  dec_yfmv =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "dec_yfmv_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ) {
    dec_yfmv[j] =
      ( vector_ptr ) getarray( info.maxMBnum, sizeof( vector ), "dec_yfmv" );
    alloc_vector( dec_yfmv[j], info );
    }
  
  dec_mv_ref_bigGOP =
    ( vector_ptr * ) getarray( num_mv_GOP, sizeof( vector_ptr ),
                               "dec_mv_ref_bigGOP" );
  mv_ref =
    ( vector_ptr * ) getarray( num_mv_GOP, sizeof( vector_ptr ),
                               "mv_ref" );
  
  /************ allocate the frame *************/
  frame_alloc( &end_of_lastGOP, info );
  
  if( info.tPyrLev >= 1 ) {
    dec_pyrFrs_bigGOP =
      ( YUVimage_ptr ) getarray( num_pyrFrs_bigGOP, sizeof( YUVimage ),
                                 "dec_pyrFrs_bigGOP_ptr" );    
    for(i = 0; i < num_pyrFrs_bigGOP; i++){
      frame_alloc( &dec_pyrFrs_bigGOP[i], info );
    }    

    dec_pyrFrs =
      ( YUVimage_ptr ) getarray( info.GOPsz, sizeof( YUVimage ),
                                 "dec_pyrFrs_ptr" );

    if( info.denoise_flag == NO ) {
      for(i = 0; i < info.GOPsz; i++){
        frame_alloc( &dec_pyrFrs[i], info );
      }    
    }
    
    dec_pyrFrs_first =
      ( YUVimage_ptr ) getarray( info.tPyrLev, sizeof( YUVimage ),
                                 "dec_pyrFrs_first_ptr" );
    for( i = 0; i < info.tPyrLev; i++ ){
      frame_alloc( &dec_pyrFrs_first[i], info ); 
    }
  
    dec_pyrTemp =
      ( YUVimage_ptr * ) getarray( info.tPyrLev, sizeof( YUVimage_ptr ),
                                   "dec_pyrTemp" );
    
    GOPsz = info.GOPsz;          
    for( i = 0; i < info.tPyrLev; i++ ) {  
      dec_pyrTemp[i] = ( YUVimage_ptr ) getarray( GOPsz, sizeof( YUVimage ),
                                              "dec_pyrTemp[]" );
      for(j = 0; j < GOPsz; j++){
        frame_alloc( &dec_pyrTemp[i][j], info );
      }
      GOPsz /= 2;
    }
    
  } else {
    assert(0);
  }

  if( info.denoise_flag == YES ) {
    Four_bigGOP =
      ( YUVimage * ) getarray( ((YUV420 == 1) ? 2 : 4) * info.GOPsz,
                               sizeof( YUVimage ), "Four_bigGOP" );
    for( j = 0; j < ((YUV420 == 1) ? 2 : 4) * info.GOPsz; j++ ) {
      frame_alloc( &Four_bigGOP[j], info ); 
    }
 
    Four_GOP =
      ( YUVimage * ) getarray( ((YUV420 == 1) ? 2 : 4) * info.GOPsz,
                               sizeof( YUVimage ), "Four_GOP" );
    for( j = 0; j < ((YUV420 == 1) ? 2 : 4) * info.GOPsz; j++ ) {
      frame_alloc( &Four_GOP[j], info ); 
    }
    
    for( i = 0; i < ((YUV420 == 1) ? 1 : 3); i++ ) { 
      spatial_high[i] =
        ( YUVimage * ) getarray( info.GOPsz, sizeof( YUVimage ),
                                 "spatial_high[i]" );
    }

  }
 
  frY = ( float ** )getarray( info.GOPsz, sizeof( float * ), "*frY" );
  frU = ( float ** )getarray( info.GOPsz, sizeof( float * ), "*frU" );
  frV = ( float ** )getarray( info.GOPsz, sizeof( float * ), "*frV" );
  
	dec_scene_change =
    ( enum FLAG ** )getarray( info.tPyrLev, sizeof( enum FLAG * ),
                              "dec_scene_change_ptr" );
    
  GOPsz = info.GOPsz;
  for( i = 0; i < info.tPyrLev; i++ ) {
    GOPsz += 1;
    dec_scene_change[i] =
      ( enum FLAG * )getarray( GOPsz, sizeof( enum FLAG ), "dec_scene_change" );
    GOPsz /= 2;
  }



    
}
