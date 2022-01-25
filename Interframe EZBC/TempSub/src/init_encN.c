#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#define EXTERN extern
#include "basic.h"
#include "structN.h"
#include "coderN.h"
#include "init_encN.h"
#include "memoryN.h"
#include "dpx.h"
#include "util_filtering.h"
#include "structN.h"
#include "general.h"
#include <iostream>

extern U16 colorLUT[1024];

#define  Yes      1
#define  No       0

extern long int totalY[5], totalU[5], totalV[5], totalmap[5], totalMV[5];
int rateAlloc( Rate_ptr FrsRate, int nfrs, videoinfo * info );

/*****************************************************************************/
/*                             init()                                        */
/*****************************************************************************/
void
init_enc( videoinfo * info )
{
  int i, j, xblk, yblk, hor, ver, bigGOP, rel_s_level;
  int num_mv_bigGOP, num_mv_GOP, num_pyrFrs_bigGOP;//, num_scene_change;
  char name[512], strtmp[512];
  FILE *fpstat, *fpbit, *fp, *fp_mv;

  int GOPsz;

  interpolate_filter(  );
  temporal_filter(  );
  rel_s_level = (info->denoise_flag == YES);

  for( i = 0; i < 5; i++ ) {
    totalY[i] = 0;
    totalU[i] = 0;
    totalV[i] = 0;
    totalmap[i] = 0;
    totalMV[i] = 0;
  }

  /************ initialize the sequence information ***************/
  hor  = info->ywidth;
  ver  = info->yheight;
  info->maxMBnum = 0;
  for (i = 0; i < info->tPyrLev; i++) {
    xblk = info->xblk[i];
    yblk = info->yblk[i];
    info->xnum[i] = ( !( hor % xblk ) ) ? hor / xblk : hor / xblk + 1;
    info->ynum[i] = ( !( ver % yblk ) ) ? ver / yblk : ver / yblk + 1;

    if(info->xnum[i] * info->ynum[i] > info->maxMBnum) {
      info->maxMBnum = info->xnum[i] * info->ynum[i];
    }
  }

  num_mv_bigGOP     = 1; // 4, 11, 26, 57, ...
  num_mv_GOP        = 1; // 4,  9, 18, 35, ...
  num_pyrFrs_bigGOP = 1; // 2,  4, 12, 27, ...
    
  // calculation of base indices for MV-sets and frames
  for( i = info->tPyrLev-1; i >= 0; i-- ){
    num_pyrFrs_bigGOP += info->bigGOP / (int)(pow ( (float)2, (int)(i + 1))); 
    num_mv_bigGOP     += info->bigGOP / (int)(pow ( (float)2, (int)i)); 
    num_mv_GOP        += info->GOPsz  / (int)(pow ( (float)2, (int)i)) + 1;
  }
 
  printf("num_pyrFrs_bigGOP = %d, num_mv_GOP = %d, num_mv_bigGOP = %d\n", num_pyrFrs_bigGOP, num_mv_GOP, num_mv_bigGOP);
  /************ allocate the 1st layer motion vectors *************/
  yfmv_bigGOP =
    ( vector_ptr * ) getarray( num_mv_bigGOP + 1, sizeof( vector_ptr ),
                               "yfmv_bigGOP_ptr" );
  for( j = 1; j <= num_mv_bigGOP; j++ ){
    yfmv_bigGOP[j] = 
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "yfmv_bigGOP" );
    alloc_vector( yfmv_bigGOP[j], *info );
  }
//Shadow
  buff_yfmv_bigGOP =
    ( vector_ptr * ) getarray( num_mv_bigGOP + 1, sizeof( vector_ptr ),
                               "buff_yfmv_bigGOP_ptr" );
  for( j = 1; j <= num_mv_bigGOP; j++ ){
    buff_yfmv_bigGOP[j] = 
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "buff_yfmv_bigGOP" );
    alloc_vector( buff_yfmv_bigGOP[j], *info );
  }
  
  yfmv =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "yfmv_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ) {
    yfmv[j] =
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "yfmv" );
    alloc_vector( yfmv[j], *info );
  }
//Shadow
  buff_yfmv =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "buff_yfmv_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ) {
    buff_yfmv[j] =
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "buff_yfmv" );
    alloc_vector( buff_yfmv[j], *info );
  }
  
  tmp_yfmv =
    ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "tmp_yfmv" );
  alloc_vector( tmp_yfmv, *info );
//Shadow
  buff_tmp_yfmv =
    ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "buff_tmp_yfmv" );
  alloc_vector( buff_tmp_yfmv, *info );
  
  mv_ref_bigGOP =
    ( vector_ptr * ) getarray( num_mv_bigGOP, sizeof( vector_ptr ),
                               "mv_ref_bigGOP" );
//Shadow
  buff_mv_ref_bigGOP =
    ( vector_ptr * ) getarray( num_mv_bigGOP, sizeof( vector_ptr ),
                               "buff_mv_ref_bigGOP" );

// for simul decoding 01.09.2018
  mv_ref =
    ( vector_ptr * ) getarray( num_mv_GOP, sizeof( vector_ptr ),
                               "mv_ref" );
//Shadow
  buff_mv_ref =
    ( vector_ptr * ) getarray( num_mv_GOP, sizeof( vector_ptr ),
                               "buff_mv_ref" );

  dec_yfmv_bigGOP =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "dec_yfmv_bigGOP_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ){
    dec_yfmv_bigGOP[j] = 
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "dec_yfmv_bigGOP" );
    alloc_vector( dec_yfmv_bigGOP[j], *info );
  }
//Shadow
  buff_dec_yfmv_bigGOP =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "buff_dec_yfmv_bigGOP_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ){
    buff_dec_yfmv_bigGOP[j] = 
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "buff_dec_yfmv_bigGOP" );
    alloc_vector( buff_dec_yfmv_bigGOP[j], *info );
  }
  
  dec_yfmv =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "dec_yfmv_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ) {
    dec_yfmv[j] =
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "dec_yfmv" );
    alloc_vector( dec_yfmv[j], *info );
    }
//Shadow
  buff_dec_yfmv =
    ( vector_ptr * ) getarray( num_mv_GOP + 1, sizeof( vector_ptr ),
                               "buff_dec_yfmv_ptr" );
  for( j = 1; j <= num_mv_GOP; j++ ) {
    buff_dec_yfmv[j] =
      ( vector_ptr ) getarray( info->maxMBnum, sizeof( vector ), "buff_dec_yfmv" );
    alloc_vector( buff_dec_yfmv[j], *info );
    }
  
  dec_mv_ref_bigGOP =
    ( vector_ptr * ) getarray( num_mv_GOP, sizeof( vector_ptr ),
                               "dec_mv_ref_bigGOP" );
//Shadow
  buff_dec_mv_ref_bigGOP =
    ( vector_ptr * ) getarray( num_mv_GOP, sizeof( vector_ptr ),
                               "buff_dec_mv_ref_bigGOP" );

// for simul decoding 01.09.2018

//#ifdef   MCTF_WITH_OBMC
	// allocate memory for frameMEinfo, varblkarray
	// pixel based motion vector side information, and side information for each block
	frameMEinfo  = (ImageMEinfo*) getarray(hor*ver, sizeof(ImageMEinfo), "frameMEinfo");
	varblkarray  = (Varblkarrayinfo*) getarray(hor/SMALLEST_BLKSIZE*ver/SMALLEST_BLKSIZE, 
		sizeof(Varblkarrayinfo), "varblkarray");
//#endif




  /*********************** allocate the frames *********************/
  frame_alloc( &end_of_lastGOP, *info );
//Shadow
  frame_alloc( &buff_end_of_lastGOP, *info );
  
  if( info->tPyrLev >= 1 ){
  
    pyrFrs_bigGOP = 
      ( YUVimage_ptr ) getarray( num_pyrFrs_bigGOP, sizeof( YUVimage ),
                                 "pyrFrs_bigGOP__ptr" );                               
    for(i = 0; i < num_pyrFrs_bigGOP; i++){
      frame_alloc( &pyrFrs_bigGOP[i], *info );
    }    
//Shadow
	buff_pyrFrs_bigGOP = 
      ( YUVimage_ptr ) getarray( num_pyrFrs_bigGOP, sizeof( YUVimage ),
                                 "buff_pyrFrs_bigGOP_ptr" );                               
    for(i = 0; i < num_pyrFrs_bigGOP; i++){
      frame_alloc( &buff_pyrFrs_bigGOP[i], *info );
    } 

    pyrFrs =
      ( YUVimage_ptr ) getarray( info->GOPsz, sizeof( YUVimage ),
                                 "pyrFrs_ptr" );
    for(i = 0; i < info->GOPsz; i++){
      frame_alloc( &pyrFrs[i], *info );
    }    
//Shadow
	buff_pyrFrs =
      ( YUVimage_ptr ) getarray( info->GOPsz, sizeof( YUVimage ),
                                 "buff_pyrFrs_ptr" );
    for(i = 0; i < info->GOPsz; i++){
      frame_alloc( &buff_pyrFrs[i], *info );
    }   
    
    pyrFrs_first = 
      ( YUVimage_ptr ) getarray( info->tPyrLev, sizeof( YUVimage ),
                                 "pyrFrs_first_ptr" );
    for( i = 0; i < info->tPyrLev; i++ ){
      frame_alloc( &pyrFrs_first[i], *info ); 
    }
//Shadow
	buff_pyrFrs_first = 
      ( YUVimage_ptr ) getarray( info->tPyrLev, sizeof( YUVimage ),
                                 "buff_pyrFrs_first_ptr" );
    for( i = 0; i < info->tPyrLev; i++ ){
      frame_alloc( &buff_pyrFrs_first[i], *info ); 
    }
  
    pyrTemp =
      ( YUVimage_ptr * ) getarray( info->tPyrLev, sizeof( YUVimage_ptr ),
                                   "pyrTemp" );
    bigGOP = info->bigGOP;          
    for( i = 0; i < info->tPyrLev; i++ ) {  
      pyrTemp[i] = ( YUVimage_ptr ) getarray( bigGOP, sizeof( YUVimage ),
                                              "pyrTemp[]" );
      for(j = 0; j < bigGOP; j++){
        frame_alloc( &pyrTemp[i][j], *info );
      }
      bigGOP /= 2;
    }
//Shadow
	buff_pyrTemp =
      ( YUVimage_ptr * ) getarray( info->tPyrLev, sizeof( YUVimage_ptr ),
                                   "buff_pyrTemp" );
    bigGOP = info->bigGOP;          
    for( i = 0; i < info->tPyrLev; i++ ) {  
      buff_pyrTemp[i] = ( YUVimage_ptr ) getarray( bigGOP, sizeof( YUVimage ),
                                              "buff_pyrTemp[]" );
      for(j = 0; j < bigGOP; j++){
        frame_alloc( &buff_pyrTemp[i][j], *info );
      }
      bigGOP /= 2;
    }
//	For simul decoding *************************************************
	dec_pyrFrs_bigGOP =
      ( YUVimage_ptr ) getarray( num_pyrFrs_bigGOP, sizeof( YUVimage ),
                                 "dec_pyrFrs_bigGOP_ptr" );    
    for(i = 0; i < num_pyrFrs_bigGOP; i++){
      frame_alloc( &dec_pyrFrs_bigGOP[i], *info );
    }    
//Shadow
	buff_dec_pyrFrs_bigGOP =
      ( YUVimage_ptr ) getarray( num_pyrFrs_bigGOP, sizeof( YUVimage ),
                                 "buff_dec_pyrFrs_bigGOP_ptr" );    
    for(i = 0; i < num_pyrFrs_bigGOP; i++){
      frame_alloc( &buff_dec_pyrFrs_bigGOP[i], *info );
    }  

    dec_pyrFrs =
      ( YUVimage_ptr ) getarray( info->GOPsz, sizeof( YUVimage ),
                                 "dec_pyrFrs_ptr" );
    if( info->denoise_flag == NO ) {
      for(i = 0; i < info->GOPsz; i++){
        frame_alloc( &dec_pyrFrs[i], *info );
      }    
    }
//Shadow
	buff_dec_pyrFrs =
      ( YUVimage_ptr ) getarray( info->GOPsz, sizeof( YUVimage ),
                                 "buff_dec_pyrFrs_ptr" );
    if( info->denoise_flag == NO ) {
      for(i = 0; i < info->GOPsz; i++){
        frame_alloc( &buff_dec_pyrFrs[i], *info );
      }    
    }
    
    dec_pyrFrs_first =
      ( YUVimage_ptr ) getarray( info->tPyrLev, sizeof( YUVimage ),
                                 "dec_pyrFrs_first_ptr" );
    for( i = 0; i < info->tPyrLev; i++ ){
      frame_alloc( &dec_pyrFrs_first[i], *info ); 
    }
//Shadow
	buff_dec_pyrFrs_first =
      ( YUVimage_ptr ) getarray( info->tPyrLev, sizeof( YUVimage ),
                                 "buff_dec_pyrFrs_first_ptr" );
    for( i = 0; i < info->tPyrLev; i++ ){
      frame_alloc( &buff_dec_pyrFrs_first[i], *info ); 
    }
  
    dec_pyrTemp =
      ( YUVimage_ptr * ) getarray( info->tPyrLev, sizeof( YUVimage_ptr ),
                                   "dec_pyrTemp" );    
    GOPsz = info->GOPsz;          
    for( i = 0; i < info->tPyrLev; i++ ) {  
      dec_pyrTemp[i] = ( YUVimage_ptr ) getarray( GOPsz, sizeof( YUVimage ),
                                              "dec_pyrTemp[]" );
      for(j = 0; j < GOPsz; j++){
        frame_alloc( &dec_pyrTemp[i][j], *info );
      }
      GOPsz /= 2;
    }
//Shadow
	buff_dec_pyrTemp =
      ( YUVimage_ptr * ) getarray( info->tPyrLev, sizeof( YUVimage_ptr ),
                                   "buff_dec_pyrTemp" );
	GOPsz = info->GOPsz;          
    for( i = 0; i < info->tPyrLev; i++ ) {  
      buff_dec_pyrTemp[i] = ( YUVimage_ptr ) getarray( GOPsz, sizeof( YUVimage ),
                                              "buff_dec_pyrTemp[]" );
      for(j = 0; j < GOPsz; j++){
        frame_alloc( &buff_dec_pyrTemp[i][j], *info );
      }
      GOPsz /= 2;
    }
  
  } else { // info->tPyrLev == 0; relevant?
    assert(0);
  }

  if( info->denoise_flag == YES ) {
    Four_GOP =
      ( YUVimage * ) getarray( 4 * info->GOPsz, sizeof( YUVimage ),
                               "Four_GOP" );
    
    for( i = 0; i < ((YUV420 == 1) ? 1 : 3); i++ ) { 
      spatial_high[i] =
        ( YUVimage * ) getarray( info->bigGOP, sizeof( YUVimage ),
                                 "spatial_high[i]" );
      for( j = 0; j < info->bigGOP; j++ ) {
        frame_alloc( &spatial_high[i][j], *info ); 
      }
    }

  }
  
  /********************* rate allocation **********************/
  rateAlloc( &FrsRate, num_mv_GOP + 1, info );

  frY = ( float ** )getarray( info->GOPsz, sizeof( float * ), "*frY" );
  frU = ( float ** )getarray( info->GOPsz, sizeof( float * ), "*frU" );
  frV = ( float ** )getarray( info->GOPsz, sizeof( float * ), "*frV" );


  /**************** output encoder information **************/
  printf( "---------------------------------------------------------------\n" );
  printf( " Input sequence   : %s\n", info->inname );
  // printf( " Coded sequence : %s\n", info->decname);
  printf( " Motion vectors   : %s\n", info->mvname );
  printf( " Bit file         : %s\n", info->bitname );
  printf( " Status log       : %s\n", info->statname );
  printf( " Sequence length  : %03d frames ( %03d - %03d )\n", 
          info->last - info->start + 1, info->start, info->last );
  printf( " Frame rate       : %0.2f frames/sec\n",
          (double)(info->framerate));
  printf( " Frame size       : (Y) %d x %d (C) %d x %d\n",
          info->ywidth << rel_s_level, info->yheight << rel_s_level,
          info->cwidth << rel_s_level, info->cheight << rel_s_level);
  if ( info->denoise_flag == YES )
    printf( " Denoise          : YES\n");
  if( info->ME == 3 ){
    for (i = 0; i < info->tPyrLev; i++) {
      printf( " Motion Level %d   : HVSBM with Bx %d By %d ... Bx %d By %d / ", 
              i, info->xblk[i], info->yblk[i], 
              info->xblk[i] >> (info->level[i] - 1), info->yblk[i] >> (info->level[i] - 1));
      if( info->subpel[i] == 4 )
        printf( "1/16-pixel accuracy\n" );
      else if( info->subpel[i] == 3 )
        printf( "1/8-pixel accuracy\n" );
      else if( info->subpel[i] == 2 )
        printf( "quarter-pixel accuracy\n" );
      else if( info->subpel[i] == 1 )
        printf( "half-pixel accuracy\n" );
      else
        printf( "full-pixel accuracy\n" );
    }
  } else {
    printf( " Motion           : unknown /  " );
  }
  printf( " Searchrange      : %d\n", info->searchrange );
  printf( " Max Searchrange  : %d\n", info->maxsearchrange );
  printf( " Lambda           : ");
  for (i = 0; i < info->tPyrLev; i++) printf("%.1f ", info->lambda[i]);
  printf("\n");
#ifdef MY_SCALE
  printf( " MY_SCALE         : %d (%d extra bits)\n", MY_SCALE, (int)log2(MY_SCALE));
#endif
  // printf( "   initial bits for GOP = %d bits\n", info->GOPbit); 
  printf( "---------------------------------------------------------------\n" );

  /******************** initialize the statefile **************/
  if( !( fpstat = fopen( info->statname, "wb" ) ) ) {
    printf( "init: %s\n", info->statname );
    exit( 1 );
  }
  fprintf( fpstat, "---------------------------------------------------------------\n" );
  fprintf( fpstat, " Input sequence   : %s\n", info->inname );
  fprintf( fpstat, " Motion vectors   : %s\n", info->mvname );
  fprintf( fpstat, " Bit file         : %s\n", info->bitname );
  fprintf( fpstat, " Status log       : %s\n", info->statname );
  fprintf( fpstat, " Sequence length  : %03d frames ( %03d - %03d )\n", 
          info->last - info->start + 1, info->start, info->last );
  fprintf( fpstat, " Frame rate       : %0.2f frames/sec\n",
          (double)(info->framerate));
  fprintf( fpstat, " Frame size       : (Y) %d x %d (C) %d x %d\n",
          info->ywidth << rel_s_level, info->yheight << rel_s_level,
          info->cwidth << rel_s_level, info->cheight << rel_s_level);
  if ( info->denoise_flag == YES )
    fprintf( fpstat, " Denoise          : YES\n");
  fprintf( fpstat, " Searchrange      : %d\n", info->searchrange );
  fprintf( fpstat, " Max searchrange  : %d\n", info->maxsearchrange );
  fprintf( fpstat, " Lambda           : ");
  for (i = 0; i < info->tPyrLev; i++) 
    fprintf(fpstat, "%.1f ", info->lambda[i]);
  fprintf(fpstat, "\n");
#ifdef MY_SCALE
  fprintf( fpstat, " MY_SCALE         : %d (%d extra bits)\n", MY_SCALE, (int)log2(MY_SCALE));
#endif
  // fprintf( fpstat, "   initial bits for GOP = %d bits\n", info->GOPbit); 
  fprintf( fpstat, "---------------------------------------------------------------\n" );
  fclose( fpstat );

  /******************* newly open the bitfile ******************/
  if( !( fpbit = fopen( info->bitname, "wb" ) ) ) {
    printf( "init: %s\n", info->bitname );
    exit( 1 );
  }
  fclose( fpbit );

  if( !( fp_mv = fopen( info->mvstatname, "wb" ) ) ) {
    printf( "can not open %s\n", info->mvstatname );
    exit( 1 );
  }
  fclose( fp_mv );

  strncpy( strtmp, info->bitname, strlen( info->bitname ) - 4 );        // Hanke, 16.09.02 
  strtmp[strlen( info->bitname ) - 4] = '\0';   // 
  sprintf( name, "%s.rd_sample_dat", strtmp );  //
  if( !( fp = fopen( name, "wb" ) ) ) {
    printf( "init: can not open %s\n", name );
    exit( 1 );
  }
  fclose( fp );

  /******************* allocation of scene changes ******************/
  // scene change description:
  scene_change =
    ( enum FLAG ** )getarray( info->tPyrLev, sizeof( enum FLAG * ),
                              "scene_change" );
  bigGOP = info->bigGOP;
  for( i = 0; i < info->tPyrLev; i++ ) {
    scene_change[i] =
      ( enum FLAG * )getarray( bigGOP + 1, sizeof( enum FLAG ), "scene_change" );
    bigGOP /= 2;
  }
//Shadow
  buff_scene_change =
    ( enum FLAG ** )getarray( info->tPyrLev, sizeof( enum FLAG * ),
                              "buff_scene_change" );
  bigGOP = info->bigGOP;
  for( i = 0; i < info->tPyrLev; i++ ) {
    buff_scene_change[i] =
      ( enum FLAG * )getarray( bigGOP + 1, sizeof( enum FLAG ), "buff_scene_change" );
    bigGOP /= 2;
  }

//For simul decoding
  dec_scene_change =
    ( enum FLAG ** )getarray( info->tPyrLev, sizeof( enum FLAG * ),
                              "dec_scene_change_ptr" );
  GOPsz = info->GOPsz;
  for( i = 0; i < info->tPyrLev; i++ ) {
    GOPsz += 1;
    dec_scene_change[i] =
      ( enum FLAG * )getarray( GOPsz, sizeof( enum FLAG ), "dec_scene_change" );
    GOPsz /= 2;
  }
//Shadow
  buff_dec_scene_change =
    ( enum FLAG ** )getarray( info->tPyrLev, sizeof( enum FLAG * ),
                              "buff_dec_scene_change_ptr" );
  GOPsz = info->GOPsz;
  for( i = 0; i < info->tPyrLev; i++ ) {
    GOPsz += 1;
    buff_dec_scene_change[i] =
      ( enum FLAG * )getarray( GOPsz, sizeof( enum FLAG ), "buff_dec_scene_change" );
    GOPsz /= 2;
  }
  
}
