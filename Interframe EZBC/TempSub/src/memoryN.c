#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#define EXTERN extern
#include "basic.h"
#include "structN.h"
#include "coderN.h"
#include "memoryN.h"
#include <iostream>

/****************************************************************************/
/*                              file: getarray.c                            */
/****************************************************************************/
char *
getarray( int num, int siz, char *name )
{
  char *ptr;
  //printf("allocating %d bytes for %s (%d pieces of %d bytes).\n", num * siz, name, num, siz);

  ptr = ( char * )calloc( num, siz );

  //printf("get_array! num = %d, siz = %d, %s \n",num,siz,name);

  if( ptr == NULL ) {
	std::cout<< "num:"<<num<<" siz:"<<siz<<"name:"<<name<<std::endl;
    printf( "getarray: can't allocate <%s>\n", name );
    exit( 1 );
  }
  return ( ptr );
}

int
rateAlloc( Rate_ptr FrsRate, int nfrs, videoinfo * info )
{
  int layer_num; 

  FrsRate->map  = ( int * )getarray( nfrs, sizeof( int ), "FrsR.map" );

  FrsRate->mv = (int **)getarray(LAYER_NUM, sizeof(int **), "FrsRate->mv");
  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
	  FrsRate->mv[layer_num] = (int *)getarray(nfrs, sizeof(int *), "FrsRate->mv[]");

  FrsRate->ybit = ( int * )getarray( nfrs, sizeof( int ), "FrsR.ybit" );
  FrsRate->ubit = ( int * )getarray( nfrs, sizeof( int ), "FrsR.ubit" );
  FrsRate->vbit = ( int * )getarray( nfrs, sizeof( int ), "FrsR.vbit" );
  FrsRate->yvar = ( float * )getarray( nfrs, sizeof( float ), "FrsR.yvar" );
  FrsRate->uvar = ( float * )getarray( nfrs, sizeof( float ), "FrsR.uvar" );
  FrsRate->vvar = ( float * )getarray( nfrs, sizeof( float ), "FrsR.vvar" );


  // for scalable motion vector coding: AGP 
  FrsRate->submv = (int ***)getarray(LAYER_NUM, sizeof(int **), "FrsRate->submv");
  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
	FrsRate->submv[layer_num] = (int **)getarray(MAX_AGP_LEVEL, sizeof(int *), "FrsRate->submv[]");
  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
  for (int i=0; i<MAX_AGP_LEVEL; i++)
	  FrsRate->submv[layer_num][i] = (int *)getarray(nfrs, sizeof(int), "FrsRate->submv[][]");


  return ( 0 );
}

int
free_rate( Rate_ptr FrsRate )
{
  free( FrsRate->map );
  free( FrsRate->mv );
  free( FrsRate->ybit );
  free( FrsRate->ubit );
  free( FrsRate->vbit );
  free( FrsRate->yvar );
  free( FrsRate->uvar );
  free( FrsRate->vvar );
  return ( 0 );
}

/*****************************************************************************/
/*                             transition filtering                          */
/*****************************************************************************/
void
alloc_transmap( transmap ** map, int hor, int ver )
{
  *map = ( transmap * )getarray( 1 , sizeof( transmap ), "transmap" ); 
  (*map)->cmode =
      ( int * )getarray( hor * ver, sizeof( int ), "transmap_cmode" );
  (*map)->dist =
      ( int * )getarray( hor * ver, sizeof( int ), "transmap_dist" );
  (*map)->tmode =
      ( int * )getarray( hor * ver, sizeof( int ), "transmap_tmode" ); 
}

void
free_transmap( transmap * map )
{
  free ( map->cmode );
  free ( map->dist  );
  free ( map->tmode );
  free ( map );
}


/*****************************************************************************/
/*                               motion vectors                              */
/*****************************************************************************/
void
alloc_vector( vector_ptr fmv, videoinfo info )
{
  int i, j, size;

  //  size = ( info.xnum ) * ( info.ynum );
  size = info.maxMBnum;
  for( i = 0; i < size; i++ ) { 
    fmv[i].child  = 0;
    fmv[i].parent = NULL;
    fmv[i].bi_mode = UNDEFINED;
    fmv[i].lifting_mode = IGNORED;
    fmv[i].is_predictor = NO;
    fmv[i].merge  = YES;
	fmv[i].merge_sign = 0; 
    for (j = 0; j < NUMBER_OF_BI_MODES; j++) {
      fmv[i].mode_info[j] = invalid_mode_info;
    }
  }
}

void
free_child( vector_ptr fmv )
{
  // if there is no child, then cut itself      
  // otherwise cut the children first and parent

  if( fmv->child == 1 ) { /* having children */
    
    /* move to free_vector() */
    free_child( fmv->child0 );
    free_child( fmv->child1 );
    free_child( fmv->child2 );
    free_child( fmv->child3 );
    
    free( fmv->child0 );
    free( fmv->child1 );
    free( fmv->child2 );
    free( fmv->child3 );
  } 
}

void
free_vector( vector_ptr fmv, videoinfo info )
{
  int i, j, size;

  /* quad tree structured should be cleaned before starting the next */
  /* so, clean the structure except 1st layer */
  /* fmv->child and fmv->lifting_mode should be initialized to CONNECTED */

  //  size = ( info.xnum ) * ( info.ynum );
  size = info.maxMBnum;
  for( i = 0; i < size; i++ ) {
    free_child( &fmv[i] );

    fmv[i].child = 0;         // reset after free_child not before
    fmv[i].merge = YES;
	fmv[i].merge_sign = 0; 
    fmv[i].mvx   = (float)HUGE_VAL;  // free_vector also used in case of
    fmv[i].mvy   = (float)HUGE_VAL;  // scene_changes (cmp. anal.c)
    fmv[i].lifting_mode = IGNORED;
    fmv[i].is_predictor = NO;
	fmv[i].bi_mode = UNDEFINED;  // added by Yongjun Wu, June 8, 2005
    for (j = 0; j < NUMBER_OF_BI_MODES; j++) {
      fmv[i].mode_info[j] = invalid_mode_info;
    }

	fmv[i].med_idx = -1;	//Added on 01.17.2016
	///////////  Added by Yuan Liu   /////////////
	fmv[i].direct_idx = NONAFF;
	fmv[i].aff_idx = -1;
	fmv[i].merge_idx = INEFFECTIVE;
	fmv[i].merge_dir = MUDA;
	fmv[i].trans_pred_idx = MUDA;

	fmv[i].aff_mrg = NO;
	fmv[i].two_comp_src = 0;
	fmv[i].skip_sign = NO;

	fmv[i].aff1_mvx   = (float)HUGE_VAL; 
    fmv[i].aff1_mvy   = (float)HUGE_VAL;
	fmv[i].aff2_mvx   = (float)HUGE_VAL; 
    fmv[i].aff2_mvy   = (float)HUGE_VAL;
	fmv[i].aff3_mvx   = (float)HUGE_VAL; 
    fmv[i].aff3_mvy   = (float)HUGE_VAL;

	fmv[i].aff1_dmvx   = (float)HUGE_VAL; 
    fmv[i].aff1_dmvy   = (float)HUGE_VAL;
	fmv[i].aff2_dmvx   = (float)HUGE_VAL; 
    fmv[i].aff2_dmvy   = (float)HUGE_VAL;
	fmv[i].aff3_dmvx   = (float)HUGE_VAL; 
    fmv[i].aff3_dmvy   = (float)HUGE_VAL;

	//////////////////////////////////////////////
  }
}

void
free_mvs( vector_ptr * yfmv, videoinfo info )
{
  int j, num_mv_GOP;

  num_mv_GOP = 1;  // 4, 9, 19, 35, ...
  
  for( j = info.tPyrLev-1; j >= 0; j-- )
    num_mv_GOP += info.GOPsz / (int)(pow ( (float)2, (int)j)) + 1;

  for( j = 1; j <= num_mv_GOP; j++ )
    free_vector( yfmv[j], info );
}


/*****************************************************************************/
/*                                   frames                                  */
/*****************************************************************************/
void
frame_alloc( YUVimage_ptr frame, videoinfo info )
{
  int i, ysize, csize;
  
  ysize =  info.ywidth * info.yheight;
  csize =  info.cwidth * info.cheight;
 
  frame->Y = ( float * )getarray( ysize, sizeof( float ), "frame" );
  frame->U = ( float * )getarray( csize, sizeof( float ), "frame" );
  frame->V = ( float * )getarray( csize, sizeof( float ), "frame" );

  for( i = 0; i < ysize; i++ ) frame->Y[i] = 0.;
  for( i = 0; i < csize; i++ ) frame->U[i] = 0.;
  for( i = 0; i < csize; i++ ) frame->V[i] = 0.;
}

void
free_frame( YUVimage frame )
{
  free( frame.Y );
  free( frame.U );
  free( frame.V );
}


/*****************************************************************************/
/*                                 destructors                               */
/*****************************************************************************/
void
enc_destructor( videoinfo info )
{
  int i, j, num_mv_bigGOP, num_mv_GOP, num_pyrFrs_bigGOP;

  num_mv_bigGOP     = 1; // 4, 11, 26, 57, ...
  num_mv_GOP        = 1; // 4,  9, 19, 35, ...
  num_pyrFrs_bigGOP = 1; // 2,  4, 12, 27, ...
    
  // calculation of base indices for MV-sets and frames
  for( i = info.tPyrLev-1; i >= 0; i-- ){
    num_pyrFrs_bigGOP += info.bigGOP / (int)(pow ( (float)2, (int)(i + 1) )); 
    num_mv_bigGOP     += info.bigGOP / (int)(pow ( (float)2, (int)i)); 
    num_mv_GOP        += info.GOPsz  / (int)(pow ( (float)2, (int)i)) + 1;
  }
 
  // printf("num_mv_bigGOP %d, num_mv_GOP %d, num_pyrFrs_bigGOP %d, num_pyrFrs %d\n",
  //        num_mv_bigGOP, num_mv_GOP, num_pyrFrs_bigGOP, info.GOPsz);
 
  /******************** free the motion vectors ********************/
    
  for( j = 1; j <= num_mv_bigGOP; j++ ){
    free_vector( yfmv_bigGOP[j], info );
    free( yfmv_bigGOP[j] );
  }
//Shadow
  for( j = 1; j <= num_mv_bigGOP; j++ ){
    free_vector( buff_yfmv_bigGOP[j], info );
    free( buff_yfmv_bigGOP[j] );
  }
  
  for( j = 1; j <= num_mv_GOP; j++ ){
    free_vector( yfmv[j], info );
    free( yfmv[j] );
  }
//Shadow
  for( j = 1; j <= num_mv_GOP; j++ ){
    free_vector( buff_yfmv[j], info );
    free( buff_yfmv[j] );
  }
  
  free_vector( tmp_yfmv, info );
//Shadow
  free_vector( buff_tmp_yfmv, info );

  /************************* free the frames ***********************/
  free_frame( end_of_lastGOP );
//Shadow
  free_frame( buff_end_of_lastGOP );

  if( info.tPyrLev >= 1 ){
  
    for(i = 0; i < num_pyrFrs_bigGOP; i++)
      free_frame( pyrFrs_bigGOP[i] );
    free( pyrFrs_bigGOP );  
//Shadow
	for(i = 0; i < num_pyrFrs_bigGOP; i++)
      free_frame( buff_pyrFrs_bigGOP[i] );
    free( buff_pyrFrs_bigGOP );  
    
    for(i = 0; i < info.GOPsz; i++)
      free_frame( pyrFrs[i] );
    free( pyrFrs );
//Shadow
	for(i = 0; i < info.GOPsz; i++)
      free_frame( buff_pyrFrs[i] );
    free( buff_pyrFrs );
    
    for( i = 0; i < info.tPyrLev; i++ )
      free_frame( pyrFrs_first[i] );
    free( pyrFrs_first );
//Shadow
	for( i = 0; i < info.tPyrLev; i++ )
      free_frame( buff_pyrFrs_first[i] );
    free( buff_pyrFrs_first );
    
  } else {  } // info->tPyrLev == 0, nie verwendet (hoffentlich!)

  if( info.denoise_flag == YES ) {
    for( i = 0; i < ((YUV420 == 1) ? 1 : 3); i++ )
      free( spatial_high[i] );
  }
    
  /************************* free scene changes *********************/
  for( i = 0; i < info.tPyrLev; i++ ) {
    free ( scene_change[i] );
  }
  free ( scene_change );
//Shadow
  for( i = 0; i < info.tPyrLev; i++ ) {
    free ( buff_scene_change[i] );
  }
  free ( buff_scene_change );

}

void
dec_destructor( videoinfo info )
{
  int i, j, num_mv_GOP, num_pyrFrs_bigGOP;

  // same number of MVs for GOP and bigGOP !
  num_mv_GOP        = 1; // 4, 9, 18, 35, ...
  num_pyrFrs_bigGOP = 1; // 3, 6, 11, 20, ...
     
  // calculation of base indices for MV-sets and frames
  for( i = info.tPyrLev-1; i >= 0; i-- ){
    num_mv_GOP += info.GOPsz / (int)(pow ( (float)2, (int)i)) + 1;
  }
  num_pyrFrs_bigGOP = info.GOPsz + info.tPyrLev;
  
  /******************** free the motion vectors ********************/
    
  for( j = 1; j <= num_mv_GOP; j++ ){
    free_vector( dec_yfmv_bigGOP[j], info );
    free( dec_yfmv_bigGOP[j] );
    free_vector( dec_yfmv[j], info );
    free( dec_yfmv[j] );
  }

  /************************* free the frames ***********************/
  free_frame( end_of_lastGOP );

  if( info.tPyrLev >= 1 ){
  
    for(i = 0; i < num_pyrFrs_bigGOP; i++)
      free_frame( dec_pyrFrs_bigGOP[i] );
    free( dec_pyrFrs_bigGOP );  
    
    for(i = 0; i < info.GOPsz; i++)
      free_frame( dec_pyrFrs[i] );
    free( dec_pyrFrs );
    
    for( i = 0; i < info.tPyrLev; i++ )
      free_frame( dec_pyrFrs_first[i] );
    free( dec_pyrFrs_first );
    
    for( i = info.tPyrLev - 1; i >= 0; i-- ) {  
      for(j = 0; j < (info.GOPsz / (int)pow( (float)2, (int)i)); j++){
        free_frame( dec_pyrTemp[i][j] );
      }
      free( dec_pyrTemp[i] ); 
    }
    free( dec_pyrTemp );
    
  }
 
  if( info.denoise_flag == YES ) {
    for( i = 0; i < ((YUV420 == 1) ? 1 : 3); i++ )
      free( spatial_high[i] );
    free( Four_GOP );
  }

  /************************* free scene changes *********************/
  for( i = 0; i < info.tPyrLev; i++ ) {
    free ( dec_scene_change[i] );
  }
  free ( dec_scene_change );

}



