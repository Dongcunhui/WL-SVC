#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define EXTERN  extern
#include "rasterfile.h"
#include "basic.h"
#include "structN.h"
#include "coderN.h"
#include "analsyn.h"
#include "ioN.h"
#include "miscN.h"
#include "memoryN.h"
#include "chrom.h"
#include "dwt_bitplane_enc.h"
#include "ezbc_enc_3d.h"
#include "bmeN.h"
#include "mvcodingN.h"
#include "pstatN.h"
#include "util_filtering.h"
#include <fstream>

#include "directional_iblock.h"
#include "obmc_varblk.h"
#include "layer_mv.h"

#define LAMBDA_SWITCH_POINT 208
#define LAMBDA_SWITCH_POINT2 64

EXTERN unsigned long int *gop_mv;

int simul_alert;

float image_entropy(float *image, int hor, int ver){
	int ii,jj,k;
	int ss,tt;
	float *tem;
	int cnt;
	double p;
	float H;

	tem = (float *)getarray( hor*ver,sizeof(float),"tem" );

	for(tt = 1;tt < ver-1; tt ++)
		for(ss = 1;ss < hor-1; ss ++)
			tem[tt*hor + ss] = (image[(tt-1)*hor + ss] + image[(tt+1)*hor + ss] + image[tt*hor + (ss-1)] + image[tt*hor + (ss+1)]) /4;

	H = 0;

	for(ii = 0;ii < 255; ii ++){
		for(jj = 0;jj < 255; jj ++){
			cnt = 0;
			for(tt = 0;tt < ver; tt ++){
				for(ss = 0;ss < hor; ss ++){
					if( floor(image[tt*hor + ss]) == ii && floor(tem[tt*hor + ss]) == jj )
						cnt ++;
				}
			}

			if(cnt > 0){
				p = ((double)cnt)/((double)(ver*hor));
				H = H - p * log((double)p)/log((double)2);
			}

		}
	}

	free(tem);
	return H;
}

void save_enc_status(videoinfo *info){
  int i, j, xblk, yblk, hor, ver, bigGOP, rel_s_level;
  int num_mv_bigGOP, num_mv_GOP, num_pyrFrs_bigGOP, GOPsz;//, num_scene_change;

  num_mv_bigGOP     = 1; // 4, 11, 26, 57, ...
  num_mv_GOP        = 1; // 4,  9, 18, 35, ...
  num_pyrFrs_bigGOP = 1; // 2,  4, 12, 27, ...
    
  // calculation of base indices for MV-sets and frames
  for( i = info->tPyrLev-1; i >= 0; i-- ){
    num_pyrFrs_bigGOP += info->bigGOP / (int)(pow ( (float)2, (int)(i + 1))); 
    num_mv_bigGOP     += info->bigGOP / (int)(pow ( (float)2, (int)i)); 
    num_mv_GOP        += info->GOPsz  / (int)(pow ( (float)2, (int)i)) + 1;
  }

//enc MV sets
  for( j = 1; j <= num_mv_bigGOP; j++ )
	  if( yfmv_bigGOP[j] != NULL ){
		free_vector(buff_yfmv_bigGOP[j], *info);
		mv_copy(yfmv_bigGOP[j], buff_yfmv_bigGOP[j], *info);
	  }

  for( j = 1; j <= num_mv_GOP; j++ )
	  if( yfmv[j] != NULL ){
		free_vector(buff_yfmv[j], *info);
		mv_copy(yfmv[j], buff_yfmv[j], *info);
	  }


  if( tmp_yfmv != NULL ){
	free_vector(buff_tmp_yfmv, *info);
	mv_copy(tmp_yfmv, buff_tmp_yfmv, *info);
  }

  for( j = 0; j < num_mv_bigGOP; j++ )
	buff_mv_ref_bigGOP[j] = mv_ref_bigGOP[j];

//Simul dec MV sets
  for( j = 0; j < num_mv_GOP; j++ )
		buff_mv_ref[j] = mv_ref[j];

  for( j = 1; j <= num_mv_GOP; j++ )
	  if( dec_yfmv_bigGOP[j] != NULL ){
		free_vector(buff_dec_yfmv_bigGOP[j], *info);
		mv_copy(dec_yfmv_bigGOP[j], buff_dec_yfmv_bigGOP[j], *info);
	  }

  for( j = 1; j <= num_mv_GOP; j++ )
	  if( dec_yfmv[j] != NULL ){
		free_vector(buff_dec_yfmv[j], *info);
		mv_copy(dec_yfmv[j], buff_dec_yfmv[j], *info);
	  }

  for( j = 0; j < num_mv_GOP; j++ )
	  buff_dec_mv_ref_bigGOP[j] = dec_mv_ref_bigGOP[j];

//enc frame sets
  copyframe(&end_of_lastGOP, &buff_end_of_lastGOP, *info);

  if( info->tPyrLev >= 1 ){
	  for(i = 0; i < num_pyrFrs_bigGOP; i++)
		  copyframe(&pyrFrs_bigGOP[i], &buff_pyrFrs_bigGOP[i], *info);

	  for(i = 0; i < info->GOPsz; i++)
		  copyframe(&pyrFrs[i], &buff_pyrFrs[i], *info);

	  for( i = 0; i < info->tPyrLev; i++ )
		  copyframe(&pyrFrs_first[i], &buff_pyrFrs_first[i], *info);

	  bigGOP = info->bigGOP;          
      for( i = 0; i < info->tPyrLev; i++ ) {  
        for(j = 0; j < bigGOP; j++){
          copyframe( &pyrTemp[i][j], &buff_pyrTemp[i][j], *info );
        }
        bigGOP /= 2;
      }

//dec frame sets
	  for(i = 0; i < num_pyrFrs_bigGOP; i++)
		  copyframe( &dec_pyrFrs_bigGOP[i], &buff_dec_pyrFrs_bigGOP[i], *info );

	  for(i = 0; i < info->GOPsz; i++)
		  copyframe( &dec_pyrFrs[i], &buff_dec_pyrFrs[i], *info );

	  for( i = 0; i < info->tPyrLev; i++ )
		  copyframe( &dec_pyrFrs_first[i], &buff_dec_pyrFrs_first[i], *info );

	  GOPsz = info->GOPsz;          
      for( i = 0; i < info->tPyrLev; i++ ) {  
        for(j = 0; j < GOPsz; j++){
          copyframe( &dec_pyrTemp[i][j], &buff_dec_pyrTemp[i][j], *info );
        }
        GOPsz /= 2;
      }

//scene change info
	  bigGOP = info->bigGOP;
	  for( i = 0; i < info->tPyrLev; i++ ) {
		  for( j = 0; j <= bigGOP ; j++ ){
				buff_scene_change[i][j] = scene_change[i][j];
		  }
		  bigGOP /= 2;
	  }

	  bigGOP = info->bigGOP;
	  for( i = 0; i < info->tPyrLev; i++ ) {
		  GOPsz += 1;
		  for( j = 0; j < GOPsz ; j++ ){
			buff_dec_scene_change[i][j] = dec_scene_change[i][j];
		  }
		  bigGOP /= 2;
	  }

  }else{
	assert(0);
  }

}

void resume_enc_status(videoinfo *info){
  int i, j, xblk, yblk, hor, ver, bigGOP, rel_s_level;
  int num_mv_bigGOP, num_mv_GOP, num_pyrFrs_bigGOP, GOPsz;//, num_scene_change;

  num_mv_bigGOP     = 1; // 4, 11, 26, 57, ...
  num_mv_GOP        = 1; // 4,  9, 18, 35, ...
  num_pyrFrs_bigGOP = 1; // 2,  4, 12, 27, ...
    
  // calculation of base indices for MV-sets and frames
  for( i = info->tPyrLev-1; i >= 0; i-- ){
    num_pyrFrs_bigGOP += info->bigGOP / (int)(pow ( (float)2, (int)(i + 1))); 
    num_mv_bigGOP     += info->bigGOP / (int)(pow ( (float)2, (int)i)); 
    num_mv_GOP        += info->GOPsz  / (int)(pow ( (float)2, (int)i)) + 1;
  }

//enc MV sets
  for( j = 1; j <= num_mv_bigGOP; j++ )
	  if( yfmv_bigGOP[j] != NULL ){
		free_vector(yfmv_bigGOP[j], *info);
		mv_copy(buff_yfmv_bigGOP[j], yfmv_bigGOP[j], *info);
	  }

  for( j = 1; j <= num_mv_GOP; j++ )
	  if( yfmv[j] != NULL ){
		free_vector(yfmv[j], *info);
		mv_copy(buff_yfmv[j], yfmv[j], *info);
	  }


  if( tmp_yfmv != NULL ){
	free_vector(tmp_yfmv, *info);
	mv_copy(buff_tmp_yfmv, tmp_yfmv, *info);
  }

  for( j = 0; j < num_mv_bigGOP; j++ )
		mv_ref_bigGOP[j] = buff_mv_ref_bigGOP[j];

//  printf("\nhere!\n");

//Simul dec MV sets
  for( j = 0; j < num_mv_GOP; j++ )
		mv_ref[j] = buff_mv_ref[j];

  for( j = 1; j <= num_mv_GOP; j++ )
	  if( dec_yfmv_bigGOP[j] != NULL ){
		free_vector(dec_yfmv_bigGOP[j], *info);
		mv_copy(buff_dec_yfmv_bigGOP[j], dec_yfmv_bigGOP[j], *info);
	  }

  for( j = 1; j <= num_mv_GOP; j++ )
	  if( dec_yfmv[j] != NULL ){
		free_vector(dec_yfmv[j], *info);
		mv_copy(buff_dec_yfmv[j], dec_yfmv[j], *info);
	  }

  for( j = 0; j < num_mv_GOP; j++ )
	  dec_mv_ref_bigGOP[j] = buff_dec_mv_ref_bigGOP[j];

//enc frame sets
  copyframe(&buff_end_of_lastGOP, &end_of_lastGOP, *info);

  if( info->tPyrLev >= 1 ){
	  for(i = 0; i < num_pyrFrs_bigGOP; i++)
		  copyframe(&buff_pyrFrs_bigGOP[i], &pyrFrs_bigGOP[i], *info);

	  for(i = 0; i < info->GOPsz; i++)
		  copyframe(&buff_pyrFrs[i], &pyrFrs[i], *info);

	  for( i = 0; i < info->tPyrLev; i++ )
		  copyframe(&buff_pyrFrs_first[i], &pyrFrs_first[i], *info);

	  bigGOP = info->bigGOP;          
      for( i = 0; i < info->tPyrLev; i++ ) {  
        for(j = 0; j < bigGOP; j++){
          copyframe( &buff_pyrTemp[i][j], &pyrTemp[i][j], *info );
        }
        bigGOP /= 2;
      }

//dec frame sets
	  for(i = 0; i < num_pyrFrs_bigGOP; i++)
		  copyframe( &buff_dec_pyrFrs_bigGOP[i], &dec_pyrFrs_bigGOP[i], *info );

	  for(i = 0; i < info->GOPsz; i++)
		  copyframe( &buff_dec_pyrFrs[i], &dec_pyrFrs[i], *info );

	  for( i = 0; i < info->tPyrLev; i++ )
		  copyframe( &buff_dec_pyrFrs_first[i], &dec_pyrFrs_first[i], *info );

	  GOPsz = info->GOPsz;          
      for( i = 0; i < info->tPyrLev; i++ ) {  
        for(j = 0; j < GOPsz; j++){
          copyframe( &buff_dec_pyrTemp[i][j], &dec_pyrTemp[i][j], *info );
        }
        GOPsz /= 2;
      }

//scene change info
	  bigGOP = info->bigGOP;
	  for( i = 0; i < info->tPyrLev; i++ ) {
		  for( j = 0; j <= bigGOP ; j++ ){
			scene_change[i][j] = buff_scene_change[i][j];
		  }
		  bigGOP /= 2;
	  }

	  bigGOP = info->bigGOP;
	  for( i = 0; i < info->tPyrLev; i++ ) {
		  GOPsz += 1;
		  for( j = 0; j < GOPsz ; j++ ){
			dec_scene_change[i][j] = buff_dec_scene_change[i][j];
		  }
		  bigGOP /= 2;
	  }

  }else{
	assert(0);
  }

}

//Added on 09.10.2017
void lambda_revise2(float *lambda, int t_level, float add_num){
	*lambda += add_num;

	if( *lambda <= 1 )
		*lambda -= add_num;
}

void anal_file_copy(char *dest_name, char *source_name){

	FILE * file1,*file2;  
    //使用二进制模式打开文件   
    file1 = fopen(source_name,"rb"); // rb 表示读   
    file2 = fopen(dest_name,"wb"); // wb 表示写   
    if(!file1)  
    {  
        printf("文件%s打开失败！",source_name);  
        return;  
    }  
    char c;  
    int index = 0;  
    fseek(file1,0,SEEK_END);        //将源文件定位到文件尾   
    int length = ftell(file1);      //获取当前位置，即文件大小（按字节算）   
    //printf("%d\n",length);        //此处可输出字节数，以进行验证   
    if(!length)  
        return;  
    while(!fseek(file1,index,SEEK_SET)) //循环定位文件，向后移动一个字节   
    {  
        fread(&c,1,1,file1);            //从源文件读取一个字节的内容到 中间变量 c   
        fwrite(&c,1,1,file2);           //将这个字节的内容写入目标文件   
        if(index == length - 1)         //如果已经读到文件尾，则跳出循环   
        {  
            break;  
        }  
        index++;                        //往后推进一个字节   
    }  
    fclose(file1);                      //关闭源文件   
    fclose(file2);                      //关闭目标文件   
}


// use left_scene and right_scene for connection checking 使用左右场景进行连接检查
void
directional_iblock_analysis_with_OBMC( YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
                   YUVimage_ptr fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
                   vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3,
                   vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
                   int t_level, int remaining_frs, videoinfo info,
				   ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray,
				   enum FLAG left_scene, enum FLAG right_scene)
{
  int x, y, X, Y, yhor, yver;
  int xnum, ynum, xblk, yblk;
  vector_ptr fmv; 

  // left scene change and right scene change
  // this is a intra coded frame 
  if( left_scene == YES &&  right_scene == YES ) return; 

  yhor = info.ywidth;
  yver = info.yheight;
  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];

  if (left_scene == NO ) 
	  fmv = fmv2;
  else if ( right_scene == NO )
	  fmv = fmv3; 
  else
	  assert(0); // validity checking 

  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
    for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {  // 逐个块进行
	  // fmv is non-NULL quad-tree root for current high temporal frame H1 
	  // fmv[Y * xnum + X] is current Macro-block 
      rec_directional_iblock_analysis_with_OBMC( &fmv[Y * xnum + X], x, y, xblk, yblk, yhor, yver, 
												L1,   H1,   H0, fr1,  fr2,  fr3, 
												fmv1,  fmv2,  fmv3, mv_ref1,  mv_ref2,  mv_ref3,
												t_level,  remaining_frs, info, fmv,
												frameMEinfo, varblkarray);
    }
  }
}


/*
 *                       temporal_analysis                              
 */

void
temporal_analysis( YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
                   YUVimage_ptr fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
                   vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3,
                   vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
                   int level, int remaining_frs, videoinfo info )
     /* fr0--previous, fr1--current(reference), fr2--next */
{
  int yhor, yver, chor, cver;
  float *ymvx1, *ymvy1, *ymvx2, *ymvy2, *ymvx3, *ymvy3;
  float *cmvx1, *cmvy1, *cmvx2, *cmvy2, *cmvx3, *cmvy3; 
  float *ymvx1_int, *ymvy1_int, *ymvx2_int, *ymvy2_int, *ymvx3_int, *ymvy3_int;
  float *cmvx1_int, *cmvy1_int, *cmvx2_int, *cmvy2_int, *cmvx3_int, *cmvy3_int; 

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;
    
  blockmv2pixelmv( fmv1, &ymvx1, &ymvy1, &cmvx1, &cmvy1, CONNECTED, info, level );
  blockmv2pixelmv( fmv2, &ymvx2, &ymvy2, &cmvx2, &cmvy2, CONNECTED, info, level );
  blockmv2pixelmv( fmv3, &ymvx3, &ymvy3, &cmvx3, &cmvy3, CONNECTED, info, level );
  blockmv2pixelmv( fmv1, &ymvx1_int, &ymvy1_int, &cmvx1_int, &cmvy1_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv2, &ymvx2_int, &ymvy2_int, &cmvx2_int, &cmvy2_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv3, &ymvx3_int, &ymvy3_int, &cmvx3_int, &cmvy3_int, PREDICTED, info, level );

  /* Y, U, V */
  mc_analysis( L1->Y, H1->Y, H0->Y, fr1->Y, fr2->Y, fr3->Y, 
               ymvx1, ymvy1, ymvx2, ymvy2, ymvx3, ymvy3,
               ymvx1_int, ymvy1_int, ymvx2_int, ymvy2_int,
               ymvx3_int, ymvy3_int, mv_ref1, mv_ref2, mv_ref3, 
               yhor, yver, level, remaining_frs, info );
  
  if( info.cwidth && info.cheight ) 
    { 
      mc_analysis( L1->U, H1->U, H0->U, fr1->U, fr2->U, fr3->U, 
                   cmvx1, cmvy1, cmvx2, cmvy2, cmvx3, cmvy3, 
                   cmvx1_int, cmvy1_int, cmvx2_int, cmvy2_int, 
                   cmvx3_int, cmvy3_int, mv_ref1, mv_ref2, mv_ref3, 
                   chor, cver, level, remaining_frs, info );
      mc_analysis( L1->V, H1->V, H0->V, fr1->V, fr2->V, fr3->V, 
                   cmvx1, cmvy1, cmvx2, cmvy2, cmvx3, cmvy3, 
                   cmvx1_int, cmvy1_int, cmvx2_int, cmvy2_int, 
                   cmvx3_int, cmvy3_int, mv_ref1, mv_ref2, mv_ref3, 
                   chor, cver, level, remaining_frs, info );
    }

  free( ymvx1 );
  free( ymvy1 );
  free( cmvx1 );
  free( cmvy1 );
  
  free( ymvx2 );
  free( ymvy2 );
  free( cmvx2 );
  free( cmvy2 );
  
  free( ymvx3 );
  free( ymvy3 );
  free( cmvx3 );
  free( cmvy3 );

  free( ymvx1_int );
  free( ymvy1_int );
  free( cmvx1_int );
  free( cmvy1_int );
  
  free( ymvx2_int );
  free( ymvy2_int );
  free( cmvx2_int );
  free( cmvy2_int );
  
  free( ymvx3_int );
  free( ymvy3_int );
  free( cmvx3_int );
  free( cmvy3_int );
}


/*
 *                       temporal_analysis_with_OBMC         带obmc的时域分解                     
 */

void
temporal_analysis_with_OBMC( YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
                   YUVimage_ptr fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
                   vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3,
                   vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
                   int level, int remaining_frs, videoinfo info,
				   ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray,
				   enum FLAG left_scene, enum FLAG right_scene)
     /* fr0--previous, fr1--current(reference), fr2--next */
{
  int    yhor, yver, chor, cver;
  float *ymvx1, *ymvy1, *ymvx2, *ymvy2, *ymvx3, *ymvy3;
  float *cmvx1, *cmvy1, *cmvx2, *cmvy2, *cmvx3, *cmvy3; 
  float *ymvx1_int, *ymvy1_int, *ymvx2_int, *ymvy2_int, *ymvx3_int, *ymvy3_int;
  float *cmvx1_int, *cmvy1_int, *cmvx2_int, *cmvy2_int, *cmvx3_int, *cmvy3_int; 

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;
    
  // connected使用预测加更新    predicted只使用预测，不更新
  blockmv2pixelmv( fmv1, &ymvx1, &ymvy1, &cmvx1, &cmvy1, CONNECTED, info, level );// 块mv转为像素mv，其中色度直接复用的亮度mv
  blockmv2pixelmv( fmv2, &ymvx2, &ymvy2, &cmvx2, &cmvy2, CONNECTED, info, level );
  blockmv2pixelmv( fmv3, &ymvx3, &ymvy3, &cmvx3, &cmvy3, CONNECTED, info, level );
  blockmv2pixelmv( fmv1, &ymvx1_int, &ymvy1_int, &cmvx1_int, &cmvy1_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv2, &ymvx2_int, &ymvy2_int, &cmvx2_int, &cmvy2_int, PREDICTED, info, level );
  blockmv2pixelmv( fmv3, &ymvx3_int, &ymvy3_int, &cmvx3_int, &cmvy3_int, PREDICTED, info, level );

  // Y componet MCTF analysis with OBMC Y分量MCTF分解
  mc_analysis_with_OBMC( L1->Y, H1->Y, H0->Y, fr1->Y, fr2->Y, fr3->Y, 
               ymvx1, ymvy1, ymvx2, ymvy2, ymvx3, ymvy3,
               ymvx1_int, ymvy1_int, ymvx2_int, ymvy2_int,
               ymvx3_int, ymvy3_int, mv_ref1, mv_ref2, mv_ref3, 
               yhor, yver, level, remaining_frs, info, frameMEinfo, varblkarray, 0,
			   left_scene, right_scene, 0);
  
  if( info.cwidth && info.cheight ) 
  { 
      mc_analysis_with_OBMC( L1->U, H1->U, H0->U, fr1->U, fr2->U, fr3->U, 
                   cmvx1, cmvy1, cmvx2, cmvy2, cmvx3, cmvy3, 
                   cmvx1_int, cmvy1_int, cmvx2_int, cmvy2_int, 
                   cmvx3_int, cmvy3_int, mv_ref1, mv_ref2, mv_ref3, 
                   chor, cver, level, remaining_frs, info, frameMEinfo, varblkarray, 1,
				   left_scene, right_scene, 1);
      mc_analysis_with_OBMC( L1->V, H1->V, H0->V, fr1->V, fr2->V, fr3->V, 
                   cmvx1, cmvy1, cmvx2, cmvy2, cmvx3, cmvy3, 
                   cmvx1_int, cmvy1_int, cmvx2_int, cmvy2_int, 
                   cmvx3_int, cmvy3_int, mv_ref1, mv_ref2, mv_ref3, 
                   chor, cver, level, remaining_frs, info, frameMEinfo, varblkarray, 1,
				   left_scene, right_scene, 0);
   }

  free( ymvx1 );
  free( ymvy1 );
  free( cmvx1 );
  free( cmvy1 );
  
  free( ymvx2 );
  free( ymvy2 );
  free( cmvx2 );
  free( cmvy2 );
  
  free( ymvx3 );
  free( ymvy3 );
  free( cmvx3 );
  free( cmvy3 );

  free( ymvx1_int );
  free( ymvy1_int );
  free( cmvx1_int );
  free( cmvy1_int );
  
  free( ymvx2_int );
  free( ymvy2_int );
  free( cmvx2_int );
  free( cmvy2_int );
  
  free( ymvx3_int );
  free( ymvy3_int );
  free( cmvx3_int );
  free( cmvy3_int );
}


// by Yongjun Wu
// In decoder mv_ref0, mv_ref1 and mv_ref2 (whether it's NULL ) can exactly indicate the scene change 
// information and the existence of motion vectors in left and/or rigth directions.
// However, in encoder mv_ref1, mv_ref2 and mv_ref3 can not exactly indicate the information. 
// Sometimes scene is changed but mv_refi is still not NULL.
// hence we have to use the correpsonding scene_change[][] to indicate the information.


/*********************************************************************/
/*                                                                   */
/*                            analscheme3 进行一个gop的时域分解，完成所有次数的时域分解                           */
/*                                                                   */
/*********************************************************************/
void
analscheme3( int curr, videoinfo info, enum FLAG first_GOP,
             enum FLAG Level_change, int remaining_frs, int simul_enc, float *upframe1, float *upframe2, float *upsamp_x )
{ 
  int i, j, k, tPyrLev, dist, yhor, yver, chor, cver;
  int num_pyrFrs_bigGOP, bigGOP, GOPsz, half_bigGOP, half_GOPsz;
  int leaving_frs, eff_GOPsz; // leaving_frs表示是否有前面bigGOP剩余的已经处理的帧
  
  int ind_mv_bigGOP[20], ind_mv_GOP[20], ind_pyrFrs_bigGOP[20];

#ifdef  DEBUG_BLOCK_MODE_MV_INFO
  int curr_mv_index, frame_type;
#endif

  enum FLAG *sc_ref1, *sc_ref2; // 是否参考标志1、2
  vector_ptr mv_ref1, mv_ref2, mv_ref3, mv_ref4;
  YUVimage *fr_cur, *fr_ref1, *fr_ref2;// 当前帧、参考帧1、2
  int total_varblk;

  float buff_lambda;

  printf("remaining_frs = %d\n",remaining_frs);
 
  // scene_change[i][j] shows whether frame i has scene change compared 
  // with its previous frame i-1, so whether we can do mctf depends on
  // entries at odd-index places.
  temporal_filter(); // 初始化系数
  
  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;
  
  tPyrLev     = info.tPyrLev; 
  bigGOP      = info.bigGOP;       // 31, 15, 7, 3
  half_bigGOP = info.bigGOP / 2;   // 15,  7, 3, 1
  GOPsz       = info.GOPsz;        // 16,  8, 4, 2
  half_GOPsz  = info.GOPsz  / 2;   //  8,  4, 2, 1
  eff_GOPsz   = info.eff_GOPsz;    // effective GOP size (cmp. encoderN.c)
  dist = 1; // frame difference of a pair of frames 或 temporal distance between fr1 and fr0
  
  // by Yongjun Wu
  // the starting points for the motion vector set in each level of bigGOP (size 31)在每一层bigGOP中，为运动向量集的起始点
  // An example of 4 level MCTF: 
  // level 3 starts from 1  and ends at 3,  level 2 starts from 4  and ends at 10
  // level 1 starts from 11 and ends at 25, level 0 starts from 26 and ends at 56

  ind_mv_bigGOP    [tPyrLev] =  0; // 1, 4, 11, 26, ...  这几个设置与encodern.c 1334处公式 设置相同

  // by Yongjun Wu
  // MV indices for transmission, e.g. for total tPyrLev = 4:
  // Level 0: 18-34; Level 1: 9-17; Level 2: 4-8; Level 3: 1-3
  // here the first set of motion vectors in each level, 
  // i.e. motion vector sets 1, 4, 9, 18 come from previous GOP (not bigGOP)  运动向量1,、4/9/18来自前面的gop
  ind_mv_GOP       [tPyrLev] = -1; // 1, 4,  9, 18, ...
  
  // indices of highpass frames in each temporal level
  ind_pyrFrs_bigGOP[tPyrLev] =  1; // 1, 2,  5, 12, ...

  // lowpass of last level plus number of all highpass frames
  num_pyrFrs_bigGOP          =  1; // 2, 5, 12, 27, ...  
   
  // calculation of base indices for MV-sets and frames:  计算这几个mv的每一层时域滤波的索引
  for( i = tPyrLev - 1; i >= 0; i-- ){
    ind_mv_bigGOP[i] = ind_mv_bigGOP[i+1] + bigGOP / (int)(pow (2, i+1));  // 
    ind_mv_GOP   [i] = 1+ ind_mv_GOP[i+1] + GOPsz  / (int)(pow (2, i+1));
    ind_pyrFrs_bigGOP[i] = ind_pyrFrs_bigGOP[i+1] + half_bigGOP / (int)(pow (2, i+1));
    num_pyrFrs_bigGOP   += half_bigGOP / (int)(pow (2, i)); 
  }

  ///////////////////////////////////
  printf("curr = %d\n\n",curr);
/*
  if(curr >= 0 && curr <= LAMBDA_SWITCH_POINT2){
	  for(i = 0; i < tPyrLev; i ++)
		lambda_revise2(&info.lambda[i],i,-1);
  }
*/
/*
  if(curr >= LAMBDA_SWITCH_POINT){
	  for(i = 0; i < tPyrLev; i ++)
		lambda_revise2(&info.lambda[i],i,5);
  }
*/  ///////////////////////////////////

  /***************
    Base Settings
  ****************/
  // 如果不是first_GOP，就会将上一个bigGOP的后半部分复制到当前的bigGOP内
  if( first_GOP == YES ) // *** first GOP ***
  { 
      leaving_frs = 0;
	  simul_alert = 0;
    
      // initialize scene changes TR 20041302
      bigGOP = info.bigGOP;
      for( i = 0; i < tPyrLev; i++ ) // 设置默认的场景改变，每一层的第一帧默认场景改变，后续帧默认场景不改变
	  {
        scene_change[i][0] = YES;  // 第一帧默认场景改变，后续帧默认场景不改变
        for( j = 1; j <= bigGOP; j++ ) { 
          scene_change[i][j] = NO; // video_scene_change[curr + j];
        }
        bigGOP /= 2;
        
        free_vector( yfmv_bigGOP[ind_mv_bigGOP[i]], info );  // the first set of motion vector does not exist always in the first GOP
        mv_ref_bigGOP[ind_mv_bigGOP[i]] = NULL;              // the first set of motion vector does not exist always in the first GOP
      }
      bigGOP = info.bigGOP;
      
      // scene_change_help / tmp_yfmv initialization
      scene_change_help = YES;
      free_vector(tmp_yfmv, info);
      mv_ref_bigGOP[1] = NULL;      // the first set of motion vector does not exist always in the first GOP 
    
  }
  else // *** next GOPs ***
  {
      bigGOP = info.bigGOP;
      leaving_frs = half_bigGOP; // = info.bigGOP - info.GOPsz;  
      // no complete GOP available, thus adapt indices for next (smaller) GOP

      for( i = 0; i < tPyrLev - 1; i++ ){

        scene_change[i][0] = YES; // do not transmit MVs 场景改变则不传输mv。如果...大于一个阈值，就会判为场景改变。

        // get scene_change information from last GOP  从上一个gop获取场景改变信息。下面两个for使得gop之间的scene_change连续
        for( j = 1; j < leaving_frs; j++ )
		{	 
          scene_change[i][j] = scene_change[i][j + GOPsz];
          // printf( "scene_change[%d][%d] = scene_change[%d][%d + GOPsz %d = %d];\n",i,j,i,j,GOPsz,j+GOPsz ); 
        }
        for( j = leaving_frs; j <= bigGOP; j++ ) 
		{ // 后半个gop的场景改变设为不改变
          scene_change[i][j] = NO; // video_scene_change[curr + j]; 
        }
    
        // put last processed L frames to correct position  by Yongjun Wu 下面的if for把前面的gop的产生的高低频帧存起来
//		if(simul_alert == 0)
			copyframe( &pyrFrs_first[i], &pyrTemp[i + 1][leaving_frs / 2 - 1], info );
        // copy last calculated highpass frames to the beginning of pyrFrs_bigGOP 
        for( j = 0; j < leaving_frs / 2; j++ ){ // (30 - 16) Frames
//			if(simul_alert == 0)
				copyframe( &pyrFrs_bigGOP[j + half_GOPsz + ind_pyrFrs_bigGOP[i]], &pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]], info );
        }
        
        for( k = 1; k < leaving_frs; k++ ){
          // BTW: number of vectors (bigGOP) per level: 2, 6, 14, 30, ... = 2x half_bigGOP
		  // copy the motion vectors from previous bigGOP to current bigGOP  by Yongjun Wu  从前面的biggop的信息拷贝到当前biggop
//		  if(simul_alert == 0){
			  free_vector ( yfmv_bigGOP[k + ind_mv_bigGOP[i]], info );
			  mv_copy( yfmv_bigGOP[k + ind_mv_bigGOP[i] + GOPsz], yfmv_bigGOP[k + ind_mv_bigGOP[i]], info );
			  // mv_ref_bigGOP[k + ind_mv_bigGOP[i]] = mv_ref_bigGOP[k + ind_mv_bigGOP[i] + GOPsz];
			  mv_ref_bigGOP[k + ind_mv_bigGOP[i]] = mv_ref_bigGOP[k + ind_mv_bigGOP[i] + GOPsz] ?
				yfmv_bigGOP[k + ind_mv_bigGOP[i]] : NULL;
//		   }//if simul
        }                  
        half_bigGOP /= 2;
        half_GOPsz /= 2;
        GOPsz /= 2;
        leaving_frs/= 2;
        bigGOP /= 2;
      }
      half_bigGOP = info.bigGOP / 2; 
      half_GOPsz = info.GOPsz / 2;
      GOPsz = info.GOPsz;  
      leaving_frs = half_bigGOP; // = info.bigGOP - info.GOPsz;
      bigGOP = info.bigGOP;
    }

  printf("leaving_frs = %d, GOPsz = %d, info.bigGOP = %d, half_bigGOP = %d\n",leaving_frs,GOPsz,info.bigGOP,half_bigGOP);
 

  /********************
    Temporal Filtering
  *********************/
  
  for( i = 0; i < tPyrLev - 1; i++ )// 逐层次的进行1.时域补偿，2.时域分解
  { 
    printf("\n");
    printf("temporal pyramid level %d (mc_anal)\n", i+1 );
       
    // mark end of sequence with a scene_changes  level_change是true表示是最后一个gop了
    if ( Level_change == YES ) {
      for (j = eff_GOPsz; j <= bigGOP; j++ )
	    scene_change[i][j] = YES;
      printf("scene_change[%d][eff_GOPsz=%d] = YES\n", i, eff_GOPsz);
      // corresponding MV will be cleared during ME
	}

	printf("bigGOP = %d, level_change = %d, eff_GOPsz = %d\n",bigGOP, Level_change, eff_GOPsz );
     
    /*** motion estimation ***/  // 运动估计。感觉就是设置了一点东西，在最后面进行了块的匹配
    assert((leaving_frs + (first_GOP==YES)) % 2 != 0);// 剩余帧是奇数或者是第一个gop
    for( j = leaving_frs + (first_GOP==YES); j < bigGOP; j += 2 ) // 从1或leaving_frs开始，到bigGop结束，一帧一帧进行处理
	{
      if (j == bigGOP - 1 && scene_change[i][j] == NO) {
		assert(0);  // 这种情况是错误的，退出

        printf("ME (uni) for frames %.3d / %.3d\n", 
               curr + (j - 1) * dist, curr + j * dist);
        sc_ref1 = &scene_change[i][j];
        sc_ref2 = NULL;
        mv_ref1 = yfmv_bigGOP[j + ind_mv_bigGOP[i]];
        mv_ref2 = NULL;
        fr_cur  = &pyrTemp[i][j];
        fr_ref1 = &pyrTemp[i][j - 1];
        fr_ref2 = NULL;
        mv_ref_bigGOP[j + ind_mv_bigGOP[i]] = yfmv_bigGOP[j + ind_mv_bigGOP[i]];
        mv_ref_bigGOP[j + 1 + ind_mv_bigGOP[i]] = NULL;
      } 
	  else if (scene_change[i][j] == NO && scene_change[i][j + 1] == NO) // 左右场景都没有改变，可以双向参考
	  {
        printf("ME (bi) for frames %.3d / %.3d / %.3d\n", 
               curr + (j - 1) * dist, curr + j * dist, curr + (j + 1) * dist);
        sc_ref1 = &scene_change[i][j];
        sc_ref2 = &scene_change[i][j + 1];
        mv_ref1 = yfmv_bigGOP[j + ind_mv_bigGOP[i]]; // 
        mv_ref2 = yfmv_bigGOP[j + 1 + ind_mv_bigGOP[i]];//

		if( (j - 2 + ind_mv_bigGOP[i]) - ind_mv_bigGOP[i] >= 0 )
			mv_ref3 = yfmv_bigGOP[j - 2 + ind_mv_bigGOP[i]];
		else
			mv_ref3 = NULL;
/*
		if( ( (first_GOP==NO) && (j == leaving_frs + 2) ) || (j == (int)( pow(2, info.tPyrLev - i) + 1 ) && (first_GOP==YES) ) )
			mv_ref3 = NULL;
*/
		if( (j - 1 + ind_mv_bigGOP[i]) - ind_mv_bigGOP[i] >= 0 )
			mv_ref4 = yfmv_bigGOP[j - 1 + ind_mv_bigGOP[i]];
		else
			mv_ref4 = NULL;
/*		
		if( ( (first_GOP==NO) && (j == leaving_frs + 2) ) || (j == (int)( pow(2, info.tPyrLev - i) + 1 ) && (first_GOP==YES) ) )
			mv_ref4 = NULL;
*/
        fr_cur  = &pyrTemp[i][j];
        fr_ref1 = &pyrTemp[i][j - 1];
        fr_ref2 = &pyrTemp[i][j + 1]; 
        mv_ref_bigGOP[j + ind_mv_bigGOP[i]] = yfmv_bigGOP[j + ind_mv_bigGOP[i]];
        mv_ref_bigGOP[j + 1 + ind_mv_bigGOP[i]] = yfmv_bigGOP[j + 1 + ind_mv_bigGOP[i]];

		printf("(int)( pow(2, info.tPyrLev - i) ) = %d\n",(int)( pow(2, info.tPyrLev - i) ));
		printf("prev left = %d, prev right = %d, ind_mv_bigGOP[i] = %d, ind_mv_GOP[i] = %d\n"
			,j - 2 + ind_mv_bigGOP[i], j - 1 + ind_mv_bigGOP[i], ind_mv_bigGOP[i], ind_mv_GOP[i]);
		printf("i = %d, j = %d, leaving_frs = %d, first_GOP = %d\n",i,j,leaving_frs, first_GOP);
      } 
	  else if (scene_change[i][j] == NO) {
        printf("ME (left) for frames %.3d / %.3d / %.3d\n", 
               curr + (j - 1) * dist, curr + j * dist, curr + (j + 1) * dist);
        sc_ref1 = &scene_change[i][j];
        sc_ref2 = NULL;
        mv_ref1 = yfmv_bigGOP[j + ind_mv_bigGOP[i]];
        mv_ref2 = NULL;

		if( (j - 2 + ind_mv_bigGOP[i]) - (ind_mv_bigGOP[i]) >= 0 )
			mv_ref3 = yfmv_bigGOP[j - 2 + ind_mv_bigGOP[i]];
		else
			mv_ref3 = NULL;
/*
		if( ( (first_GOP==NO) && (j == leaving_frs + 2) ) || (j == (int)( pow(2, info.tPyrLev - i) + 1 ) && (first_GOP==YES) ) )
			mv_ref3 = NULL;
*/
		mv_ref4 = NULL;

        fr_cur  = &pyrTemp[i][j];
        fr_ref1 = &pyrTemp[i][j - 1];
        fr_ref2 = NULL; 
        mv_ref_bigGOP[j + ind_mv_bigGOP[i]] = yfmv_bigGOP[j + ind_mv_bigGOP[i]];
        mv_ref_bigGOP[j + 1 + ind_mv_bigGOP[i]] = NULL;
        assert(scene_change[i][j + 1] == YES);

		printf("prev left = %d, prev right = %d, ind_mv_bigGOP[i] = %d, ind_mv_GOP[i] = %d\n"
			,j - 2 + ind_mv_bigGOP[i], j - 1 + ind_mv_bigGOP[i], ind_mv_bigGOP[i], ind_mv_GOP[i]);
		printf("i = %d, j = %d, leaving_frs = %d, first_GOP = %d\n",i,j,leaving_frs, first_GOP);
      } 
	  else if (scene_change[i][j + 1] == NO) {
//		assert(0);

        printf("ME (right) for frames %.3d / %.3d / %.3d\n", 
               curr + (j - 1) * dist, curr + j * dist, curr + (j + 1) * dist);
        sc_ref1 = &scene_change[i][j + 1];
        sc_ref2 = NULL;
        mv_ref1 = yfmv_bigGOP[j + 1 + ind_mv_bigGOP[i]];
        mv_ref2 = NULL;
        fr_cur  = &pyrTemp[i][j];
        fr_ref1 = &pyrTemp[i][j + 1];
        fr_ref2 = NULL;
        mv_ref_bigGOP[j + ind_mv_bigGOP[i]] = NULL;
        mv_ref_bigGOP[j + 1 + ind_mv_bigGOP[i]] = yfmv_bigGOP[j + 1 + ind_mv_bigGOP[i]];
        assert(scene_change[i][j] == YES);
        // free_vector( yfmv_bigGOP[j + ind_mv_bigGOP[i]], info );  
      } 
	  else {
        sc_ref1 = sc_ref2 = NULL;
        mv_ref1 = mv_ref2 = NULL;
		mv_ref3 = NULL;
		mv_ref4 = NULL;
        fr_cur = fr_ref1 = fr_ref2 = NULL; 
        mv_ref_bigGOP[j + ind_mv_bigGOP[i]] = NULL;
        mv_ref_bigGOP[j + 1 + ind_mv_bigGOP[i]] = NULL;
      }

      // free_vector before motion estimation
	  free_vector( yfmv_bigGOP[j + ind_mv_bigGOP[i]], info );       
	  free_vector( yfmv_bigGOP[j + 1 + ind_mv_bigGOP[i]], info ); 

#ifndef VBR_DEC
	  if( j == (leaving_frs + (first_GOP==YES)) && simul_skip == YES){
		buff_lambda = info.lambda[i];

		if(i <= 2){
			info.lambda[i] = orig_lambda[i];
		}else{
			info.lambda[i] = orig_lambda[2];
		}

	  }
#endif

      if (sc_ref1 != NULL) { // 左向参考存在，进行划分，补偿
        block_matching(mv_ref1, mv_ref2, mv_ref3, mv_ref4, fr_cur, fr_ref1, fr_ref2, 
                       sc_ref1, sc_ref2, info, i, dist, info.subpel[i], upframe1, upframe2, upsamp_x);
      }
#ifndef VBR_DEC
	  if( j == (leaving_frs + (first_GOP==YES)) && simul_skip == YES){
		info.lambda[i] = buff_lambda;
	  }
#endif
    } 
    
    /*** temporal analysis 时域分解***/
	// 从1或leaving_frs开始，到half bigGop结束，一帧一帧进行处理完本层的时域小波的所有帧
    for( j = (half_GOPsz - 1)*(first_GOP == NO); j < half_bigGOP; j++ )
	{	
#ifdef   MCTF_WITH_OBMC
		// get the side information for each block: mode, mv and so on
		// in encoder use "scene_change[i][2 * j + 1], scene_change[i][2 * j + 2], 1"
		// for connection checking
		get_mv_side_information(info, // 获得mv信息，放在了imageneinfo中
	                           yfmv_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],    // fmv2 
                               yfmv_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],    // fmv3
		                       mv_ref_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],  // mv_ref2
                               mv_ref_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],  // mv_ref3
  							   frameMEinfo,    varblkarray, 
							   &total_varblk, i, j,
							   scene_change[i][2 * j + 1], scene_change[i][2 * j + 2], 1);
		// get the weighting coefficients for self mv and neighbor mvs 获得mv的权重系数
		mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
		// for the analysis of directional iblock 
		// for convenience we separate from the function of temporal_analysis()
#ifdef DEGUB_DIRECTIONAL_IBLOCK_OBMC
		FILE *fiblk = fopen("iblock_info_encoder.txt", "at");
		fprintf(fiblk, "i=%d\t j=%d\n", i, j); 
		fclose(fiblk); 
#endif 

		directional_iblock_analysis_with_OBMC( 
			             &pyrTemp[i + 1][j],               // L1
						 &pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]],     // H1
						 &pyrFrs_bigGOP[j - 1 + ind_pyrFrs_bigGOP[i]], // H0 &dump for j=0
						 &pyrTemp[i][2 * j],                           // fr1
						 &pyrTemp[i][2 * j + 1],                       // fr2
						 &pyrTemp[i][2 * j + 2],                       // fr3
						 yfmv_bigGOP[2 * j + ind_mv_bigGOP[i]],        // fmv1
						 yfmv_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],    // fmv2 
						 yfmv_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],    // fmv3
						 mv_ref_bigGOP[2 * j + ind_mv_bigGOP[i]],      // mv_ref1
						 mv_ref_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],  // mv_ref2
						 mv_ref_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],  // mv_ref3
						 i, remaining_frs, info, 
						 frameMEinfo, varblkarray,
						 scene_change[i][2 * j + 1], scene_change[i][2 * j + 2]);
#ifdef DEGUB_DIRECTIONAL_IBLOCK_OBMC
		fiblk = fopen("iblock_info_encoder.txt", "at");
		fprintf(fiblk, "\n\n"); 
		fclose(fiblk); 
#endif 

#endif 
		printf("level and frame: i = %d, j = %d\n",i,j);

        temporal_analysis_with_OBMC( 
		                 &pyrTemp[i + 1][j],                           // L1 // 下一层的L1放入这个里面
                         &pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]],     // H1  // 生成的H1
                         &pyrFrs_bigGOP[j - 1 + ind_pyrFrs_bigGOP[i]], // &dump for j=0
                         &pyrTemp[i][2 * j],                           // fr1
                         &pyrTemp[i][2 * j + 1],                       // fr2
                         &pyrTemp[i][2 * j + 2],                       // fr3
                         yfmv_bigGOP[2 * j + ind_mv_bigGOP[i]],        // fmv1
                         yfmv_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],    // fmv2
                         yfmv_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],    // fmv3
                         mv_ref_bigGOP[2 * j + ind_mv_bigGOP[i]],      // mv_ref1
                         mv_ref_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],  // mv_ref2
                         mv_ref_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],  // mv_ref3
                         i, remaining_frs, info,
						 frameMEinfo, varblkarray,
						 scene_change[i][2 * j + 1], scene_change[i][2 * j + 2]);
#else 
	  // NO OBMC for MCTF
	  // there are always the first layer quad-tree in yfmv_bigGOP
      // so we need mv_ref_bigGOP to indicate whether the motion vector set is NULL or not
	  printf("level and frame: i = %d, j = %d\n",i,j);

	  get_mv_side_information(info, 
	                           yfmv_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],    // fmv2 
                               yfmv_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],    // fmv3
		                       mv_ref_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],  // mv_ref2
                               mv_ref_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],  // mv_ref3
  							   frameMEinfo,    varblkarray, 
							   &total_varblk, i, j,
							   scene_change[i][2 * j + 1], scene_change[i][2 * j + 2], 1);
	  // get the weighting coefficients for self mv and neighbor mvs 
	  mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
	  // for the analysis of directional iblock 
	  // for convenience we separate from the function of temporal_analysis()

	  // get the side information for each block: mode, mv and so on
	  // in encoder use "scene_change[i][2 * j + 1], scene_change[i][2 * j + 2], 1"
	  // for connection checking

	  directional_iblock_analysis_with_OBMC( 
			             &pyrTemp[i + 1][j],               // L1
						 &pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]],     // H1
						 &pyrFrs_bigGOP[j - 1 + ind_pyrFrs_bigGOP[i]], // H0 &dump for j=0
						 &pyrTemp[i][2 * j],                           // fr1
						 &pyrTemp[i][2 * j + 1],                       // fr2
						 &pyrTemp[i][2 * j + 2],                       // fr3
						 yfmv_bigGOP[2 * j + ind_mv_bigGOP[i]],        // fmv1
						 yfmv_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],    // fmv2 
						 yfmv_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],    // fmv3
						 mv_ref_bigGOP[2 * j + ind_mv_bigGOP[i]],      // mv_ref1
						 mv_ref_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],  // mv_ref2
						 mv_ref_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],  // mv_ref3
						 i, remaining_frs, info, 
						 frameMEinfo, varblkarray,
						 scene_change[i][2 * j + 1], scene_change[i][2 * j + 2]);

#endif

      temporal_analysis( &pyrTemp[i + 1][j],                           // L1
                         &pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]],     // H1
                         &pyrFrs_bigGOP[j - 1 + ind_pyrFrs_bigGOP[i]], // H0 &dump for j=0
                         &pyrTemp[i][2 * j],                           // fr1
                         &pyrTemp[i][2 * j + 1],                       // fr2
                         &pyrTemp[i][2 * j + 2],                       // fr3
                         yfmv_bigGOP[2 * j + ind_mv_bigGOP[i]],        // fmv1
                         yfmv_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],    // fmv2
                         yfmv_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],    // fmv3
                         mv_ref_bigGOP[2 * j + ind_mv_bigGOP[i]],
                         mv_ref_bigGOP[2 * j + 1 + ind_mv_bigGOP[i]],
                         mv_ref_bigGOP[2 * j + 2 + ind_mv_bigGOP[i]],
                         i, remaining_frs, info );

#endif

      if (scene_change[i][2 * j] == YES && scene_change[i][2 * j + 1] == YES) 
	  {

		  printf("enter here! i = %d, j = %d\n",i,j);

#ifdef COPYCOMPENSATION_WEIGHTING 
        wcopyframe(&pyrTemp[i][2 * j], &pyrTemp[i + 1][j], 
                   LPW4[1] * copycomp_weight_low[i], info);
#else
        // propagate L frame with factor 1
        copyframe(&pyrTemp[i][2 * j], &pyrTemp[i + 1][j], info);
#endif

        // propagate double scene change
        scene_change[i + 1][j] = YES;
        scene_change[i + 1][j + 1] = YES;
      }

      if (scene_change[i][2 * j + 1] == YES && scene_change[i][2 * j + 2] == YES) 
	  {
#ifdef COPYCOMPENSATION_WEIGHTING 
        wcopyframe(&pyrTemp[i][2 * j + 1], 
                   &pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]], 
                   HPW4[1] * copycomp_weight_high[i], info);
#else
        // propagate H frame with factor 1
        copyframe(&pyrTemp[i][2 * j + 1], &pyrFrs_bigGOP[j + ind_pyrFrs_bigGOP[i]], info);
#endif
      }
    }//for j 一层的所有帧
    
    // save frames for spatial analysis  只保存高频帧.从后往前放，因此越低频帧越在前面
    for( k = 0; k < half_GOPsz; k++ ) {
      copyframe( &pyrFrs_bigGOP[k + ind_pyrFrs_bigGOP[i]],
                 &pyrFrs[k + half_GOPsz], info );  // these high temporal frames belong to current GOP, which will be coded for this GOP
	}
    
    // copy MV-sets, the motion vector sets in current GOP, which will be encoded for this GOP, by Yongjun复制当前GOPsz的mv
    for( k = 0; k <= GOPsz; k++ ) {
      free_vector ( yfmv[k + ind_mv_GOP[i]], info );
      mv_copy( yfmv_bigGOP[k + ind_mv_bigGOP[i]], yfmv[k + ind_mv_GOP[i]], info ); 
    }

    bigGOP /= 2;  
    GOPsz  /= 2;
    half_bigGOP /= 2;
    half_GOPsz  /= 2;
    leaving_frs /= 2; 
    dist *= 2; 

    if ( Level_change == YES ) {
      eff_GOPsz = (int) ceil ((double)remaining_frs / (pow (2, i + 1)));  
    } 
	else {
      eff_GOPsz = bigGOP;
    }
    
  } // i  
  
  /**************/ 
  /* last level */
  /**************/        
  printf("\n");
  printf("*** last level ***\n");
  printf("temporal pyramid level %d (mc_anal)\n", tPyrLev );

  printf("bigGOP = %d, level_change = %d, eff_GOPsz = %d\n",bigGOP, Level_change, eff_GOPsz );

  scene_change[tPyrLev - 1][1] = NO;
  scene_change[tPyrLev - 1][2] = NO;

  // mark end of sequence with a scene_change
  if ( Level_change == YES ) {
    for (j = eff_GOPsz; j <= 2; j++ ) {
      scene_change[tPyrLev - 1][j] = YES;
    }
    printf("scene_change[%d][eff_GOPsz=%d] = YES\n", tPyrLev-1, eff_GOPsz);
  }
  
  // scene_change from preceeding GOP 
  scene_change[tPyrLev - 1][0] = scene_change_help;

  if(scene_change[tPyrLev - 1][1] == NO && scene_change[tPyrLev - 1][2] == NO){
    printf("ME (bi) for frames %.3d / %.3d / %.3d\n", 
           curr, curr + dist, curr + 2 * dist);
    sc_ref1 = &scene_change[tPyrLev - 1][1];
    sc_ref2 = &scene_change[tPyrLev - 1][2];
    mv_ref1 = mv_ref_bigGOP[2] = yfmv[2];
    mv_ref2 = mv_ref_bigGOP[3] = yfmv[3];

	mv_ref3 = NULL;
	mv_ref4 = NULL;

    fr_cur  = &pyrTemp[tPyrLev - 1][1];
    fr_ref1 = &pyrTemp[tPyrLev - 1][0];
    fr_ref2 = &pyrTemp[tPyrLev - 1][2]; 
  } else if (scene_change[tPyrLev - 1][1] == NO) {
    printf("ME (left) for frames %.3d / %.3d / %.3d\n", 
           curr, curr + dist, curr + 2 * dist);
    sc_ref1 = &scene_change[tPyrLev - 1][1];
    sc_ref2 = NULL;
    mv_ref1 = mv_ref_bigGOP[2] = yfmv[2];
    mv_ref2 = mv_ref_bigGOP[3] = NULL;

	mv_ref3 = NULL;
	mv_ref4 = NULL;

    fr_cur  = &pyrTemp[tPyrLev - 1][1];
    fr_ref1 = &pyrTemp[tPyrLev - 1][0];
    fr_ref2 = NULL;
    assert(scene_change[tPyrLev - 1][2] == YES);
    mv_ref_bigGOP[3] = NULL; 
  } else if (scene_change[tPyrLev - 1][2] == NO) {

    printf("ME (right) for frames %.3d / %.3d / %.3d\n", 
           curr, curr + dist, curr + 2 * dist);
    sc_ref1 = &scene_change[tPyrLev - 1][2];
    sc_ref2 = NULL;
    mv_ref1 = mv_ref_bigGOP[3] = yfmv[3];
    mv_ref2 = mv_ref_bigGOP[2] = NULL;

	mv_ref3 = NULL;
	mv_ref4 = NULL;

    fr_cur  = &pyrTemp[tPyrLev - 1][1];
    fr_ref1 = &pyrTemp[tPyrLev - 1][2];
    fr_ref2 = NULL;
    assert(scene_change[tPyrLev - 1][1] == YES);
    mv_ref_bigGOP[2] = NULL; 
  } else {
    sc_ref1 = sc_ref2 = NULL;
    mv_ref1 = mv_ref_bigGOP[2] = NULL;
    mv_ref2 = mv_ref_bigGOP[3] = NULL;
    fr_cur = fr_ref1 = fr_ref2 = NULL;
  }
  // free_vector before motion estimation
  free_vector( yfmv[2], info );
  free_vector( yfmv[3], info );

  assert( (tPyrLev - 1) >= 2);

#ifndef VBR_DEC
  if( simul_skip == YES){
	buff_lambda = info.lambda[tPyrLev - 1];
	info.lambda[tPyrLev - 1] = orig_lambda[tPyrLev - 1];
  }
#endif
  
  if (sc_ref1 != NULL) {
    block_matching(mv_ref1, mv_ref2, mv_ref3, mv_ref4, fr_cur, fr_ref1, fr_ref2, 
		sc_ref1, sc_ref2, info, i, dist, info.subpel[i], upframe1, upframe2, upsamp_x);
  }

#ifndef VBR_DEC
  if( simul_skip == YES){
	info.lambda[tPyrLev - 1] = buff_lambda;
  }
#endif

#ifdef   MCTF_WITH_OBMC
  get_mv_side_information(info, 
	                      yfmv[2],				// fmv2 
                          yfmv[3],				// fmv3
		                  mv_ref_bigGOP[2],     // mv_ref2
                          mv_ref_bigGOP[3],     // mv_ref3
  						  frameMEinfo,  varblkarray, 
						  &total_varblk, tPyrLev - 1, 1,
						  scene_change[tPyrLev - 1][1], scene_change[tPyrLev - 1][2], 1);

  mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
	  // for the analysis of directional iblock 
	  // for convenience we separate from the function of temporal_analysis()
  directional_iblock_analysis_with_OBMC( &pyrFrs[0],				// L1
							   &pyrFrs[1],							// H1
							   &end_of_lastGOP, // dump in 1.GOP    // H0 &dump for j=0
                               &pyrTemp[tPyrLev - 1][0],		    // fr1  
                               &pyrTemp[tPyrLev - 1][1],            // fr2
                               &pyrTemp[tPyrLev - 1][2],            // fr3
                               tmp_yfmv, //yfmv[1] in 1.GOP         // fmv1
                               yfmv[2],                             // fmv2
							   yfmv[3],                             // fmv3
							   mv_ref_bigGOP[1],                    // mv_ref1
							   mv_ref_bigGOP[2],                    // mv_ref2
							   mv_ref_bigGOP[3],                    // mv_ref3
							   tPyrLev - 1, remaining_frs, info,
							   frameMEinfo, varblkarray,
							   scene_change[tPyrLev - 1][1], scene_change[tPyrLev - 1][2]);

#endif 

  temporal_analysis_with_OBMC( 
	                 &pyrFrs[0],							// L1
                     &pyrFrs[1],                            // H1
                     &end_of_lastGOP, // dump in 1.GOP      // H0
                     &pyrTemp[tPyrLev - 1][0],              // fr1
                     &pyrTemp[tPyrLev - 1][1],              // fr2
                     &pyrTemp[tPyrLev - 1][2],              // fr3
                     tmp_yfmv, //yfmv[1] in 1.GOP           // fmv1
                     yfmv[2],                               // fmv2
					 yfmv[3],                               // fmv3
                     mv_ref_bigGOP[1],                      // mv_ref1
                     mv_ref_bigGOP[2],                      // mv_ref2
                     mv_ref_bigGOP[3],                      // mv_ref3
                     tPyrLev - 1, remaining_frs, info,
					 frameMEinfo, varblkarray,
					 scene_change[tPyrLev - 1][1], scene_change[tPyrLev - 1][2]);

#else 

  get_mv_side_information(info, 
	                      yfmv[2],				// fmv2 
                          yfmv[3],				// fmv3
		                  mv_ref_bigGOP[2],     // mv_ref2
                          mv_ref_bigGOP[3],     // mv_ref3
  						  frameMEinfo,  varblkarray, 
						  &total_varblk, tPyrLev - 1, 1,
						  scene_change[tPyrLev - 1][1], scene_change[tPyrLev - 1][2], 1);

  mv_weight_info(info, frameMEinfo, varblkarray, total_varblk, i, j);

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
	  // for the analysis of directional iblock 
	  // for convenience we separate from the function of temporal_analysis()

  directional_iblock_analysis_with_OBMC( &pyrFrs[0],				// L1
							   &pyrFrs[1],							// H1
							   &end_of_lastGOP, // dump in 1.GOP    // H0 &dump for j=0
                               &pyrTemp[tPyrLev - 1][0],		    // fr1  
                               &pyrTemp[tPyrLev - 1][1],            // fr2
                               &pyrTemp[tPyrLev - 1][2],            // fr3
                               tmp_yfmv, //yfmv[1] in 1.GOP         // fmv1
                               yfmv[2],                             // fmv2
							   yfmv[3],                             // fmv3
							   mv_ref_bigGOP[1],                    // mv_ref1
							   mv_ref_bigGOP[2],                    // mv_ref2
							   mv_ref_bigGOP[3],                    // mv_ref3
							   tPyrLev - 1, remaining_frs, info,
							   frameMEinfo, varblkarray,
							   scene_change[tPyrLev - 1][1], scene_change[tPyrLev - 1][2]);

#endif 

  temporal_analysis( &pyrFrs[0],							// L1
                     &pyrFrs[1],                            // H1
                     &end_of_lastGOP, // dump in 1.GOP      // H0
                     &pyrTemp[tPyrLev - 1][0],              // fr1
                     &pyrTemp[tPyrLev - 1][1],              // fr2
                     &pyrTemp[tPyrLev - 1][2],              // fr3
                     tmp_yfmv, //yfmv[1] in 1.GOP           // fmv1
                     yfmv[2],                               // fmv2
					 yfmv[3],                               // fmv3
                     mv_ref_bigGOP[1],                      // mv_ref1
                     mv_ref_bigGOP[2],                      // mv_ref2
                     mv_ref_bigGOP[3],                      // mv_ref3
                     tPyrLev - 1, remaining_frs, info );

#endif 


  if (scene_change[tPyrLev - 1][0] == YES && scene_change[tPyrLev - 1][1] == YES) {
#ifdef COPYCOMPENSATION_WEIGHTING 
    wcopyframe(&pyrTemp[tPyrLev - 1][0], &pyrFrs[0], 
               LPW4[1] * copycomp_weight_low[tPyrLev - 1], info);
#else
    // propagate L frame with factor 1
    copyframe(&pyrTemp[tPyrLev - 1][0], &pyrFrs[0], info);
#endif
  }
  if (scene_change[tPyrLev - 1][1] == YES && scene_change[tPyrLev - 1][2] == YES) {
#ifdef COPYCOMPENSATION_WEIGHTING 
    wcopyframe(&pyrTemp[tPyrLev - 1][1], &pyrFrs[1], 
               HPW4[1] * copycomp_weight_high[tPyrLev - 1], info);
#else
    // propagate H frame with factor 1
    copyframe(&pyrTemp[tPyrLev - 1][1], &pyrFrs[1], info);
#endif
  }

//	if(simul_enc == NO){
	  // save last calculated lowpass frames of bigGOP
	  bigGOP = info.bigGOP;
	  for( i = 0; i < tPyrLev - 1; i++ ){
		copyframe( &pyrTemp[i + 1][bigGOP / 2 - 1],
				   &pyrFrs_first[i], info );
		bigGOP /= 2;
	  }
	  scene_change[tPyrLev - 1][0] = YES; // just to be sure
 
	  // save last calculated highpass frame of last level
	  copyframe( &pyrFrs[1], &end_of_lastGOP, info );
  
	  // save yfmv[3] for next GOP
	  free_vector ( tmp_yfmv, info );
	  mv_copy( yfmv[3], tmp_yfmv, info );
	  mv_ref_bigGOP[1] = mv_ref_bigGOP[3] ? tmp_yfmv : NULL;
  
	  // save scene_change for next GOP
	  scene_change_help = scene_change[tPyrLev - 1][2];
//  }// if not simulative encoding

  free_vector( yfmv[1], info );

  if(simul_enc == YES)
	  simul_alert = 1;
  else
	  simul_alert = 0;

}


/*
 *                                  mc_anal
 进行mctf  当前第几帧、当前第几个gop、info、是否为第一个gop、是否最后一个gop不足一个biggop、剩余帧。预估码率、
 */
long int
mc_anal( int curr, int GOP_counter, videoinfo info,
         enum FLAG first_GOP, enum FLAG Level_change,
         int remaining_frs, int simul_enc, float *upframe1, float *upframe2, float *upsamp_x )
{
  int GOPheader_bytes, mvbits; // GOP头的花费bit， mv的bit
  FILE *fpstat; // 写stst文件

  if( !( fpstat = fopen( info.statname, "at+" ) ) ) {
    printf( "Can not open %s\n", info.statname );
    exit( 1 );
  }
  fprintf( fpstat, " frame %.3d - %.3d (%.3d) ......\n", 
           curr, curr + info.GOPsz - 1, curr + info.eff_GOPsz -1 );
  fclose( fpstat );
  
  if( info.tPyrLev >= 1 ) { // 进行时域分解  对一个gop内的几个时域进行分解
    analscheme3( curr, info, first_GOP, Level_change, remaining_frs, simul_enc, upframe1, upframe2, upsamp_x );

    GOPheader_bytes = write_GOPheader( scene_change, info );  // 写sence_change进码流
    mvbits = mv_encoding( info, FrsRate, yfmv, GOP_counter, simul_enc, curr );
    free_mvs( yfmv, info );   // initialize MVs for next GOP

	printf("GOPheader_bytes = %d, mvbits = %d\n",GOPheader_bytes, mvbits);
 
  } 
  else // 不进行时域分解
  { // info.tPyrLev == 0 // necessary? KH, 2003-12-13 
    mvbits = 0;
    // no GOPheader output          
    copyframe( &pyrTemp[0][0], &pyrFrs[0], info );
  }

  // print_mvbits( info, FrsRate );        //pstatN.c

  return mvbits;
}


/*
 * denoise_mctf_anal_ezbc()
 *
 */
long int
denoise_mctf_anal_ezbc( int curr, int GOP_counter, videoinfo info,
                        enum FLAG first_GOP, enum FLAG Level_change,
                        int remaining_frs )
{
  int i, j, filterType = DENOISE_FILTER;
  YUVimage tempfr;
  long int output_GOP_bytes;
  char strtmp[256], rd_name[256];

  for( j = 0; j < info.GOPsz; j++ ) {
    Four_GOP[j] = pyrFrs[j];
  }
  for( i = 0; i < ((YUV420 == 1) ? 1 : 3); i++ ) {
    for( j = 0; j < info.GOPsz; j++ ) {
      Four_GOP[( i + 1 ) * info.GOPsz + j] = spatial_high[i][j];
    }
  }

  /************/
  /* spatial  */
  /* analysis */
  /************/
  info.ywidth  *= 2;
  info.yheight *= 2;
  info.cwidth  *= 2;
  info.cheight *= 2;

  for( i = 0; i < info.eff_GOPsz; i++ ) {

    frame_alloc( &tempfr, info );           
    read_frame( &tempfr, info, info.inname, curr + i, info.format );

    if( YUV420 == 1 ) { // down-conversion to 420 
      f444_420( info, &tempfr ); 
      spatial_anal( tempfr.Y, info.ywidth, info.yheight, pyrTemp[0][i].Y,
                    spatial_high[0][i].Y, spatial_high[0][i].U,
                    spatial_high[0][i].V, filterType );
      for( j = 0; j < info.cwidth * info.cheight / 4; j++ ) {
        pyrTemp[0][i].U[j] = tempfr.U[j] / 4;
        pyrTemp[0][i].V[j] = tempfr.V[j] / 4;
      }
    } else {
      spatial_anal_frame( &tempfr, info, &pyrTemp[0][i], &spatial_high[0][i],
                          &spatial_high[1][i], &spatial_high[2][i],
                          filterType );
    }

    free_frame( tempfr );

  }

  info.ywidth  /= 2;
  info.yheight /= 2;
  info.cwidth  /= 2;
  info.cheight /= 2; 

  /********/
  /* MCTF */
  /********/
  mc_anal( curr, GOP_counter, info, first_GOP, Level_change, remaining_frs, NO, NULL, NULL, NULL );  

  /********/
  /* EZBC */
  /********/
  output_GOP_bytes = ezbc3d_enc_GOP( curr, Four_GOP, info, GOP_counter );
  
  return output_GOP_bytes;
}


/*
 * mctf_anal_ezbc()
 * 编码一个gop
 */
long int
mctf_anal_ezbc( int curr, int GOP_counter, videoinfo info,
                enum FLAG first_GOP, enum FLAG Level_change,
                int remaining_frs, long int *sum_mv, int simul_enc, float *upframe1, float *upframe2, float *upsamp_x )
{
  int i,j;
  long int output_GOP_bytes;
  long target_bitrate;
  long subband_budget;

  long int mv_bits;

  int t_level,dist,nfrs;
  float bpp, pwr_factor;
  float H;

  printf("\nsimul_enc = %d\n",simul_enc);

  target_bitrate = (info.simul_rate * 1000 * info.GOPsz) / (8 * info.framerate);

  for( i = 0; i < info.eff_GOPsz; i++ ) {  // 读取eff_gop的帧
    read_frame( &pyrTemp[0][i], info, info.inname, curr + i, info.format );
  }

  //H = image_entropy((&pyrTemp[0][0])->Y,info.ywidth,info.yheight);
//  printf("image entropy for curr %d: %f\n",curr,H);

  /********/
  /* MCTF */
  /********/
  // 处理一个gop的信息
  *sum_mv = mc_anal( curr, GOP_counter, info, first_GOP, Level_change, remaining_frs, simul_enc, upframe1, upframe2, upsamp_x );
  *sum_mv = *sum_mv / 8;

  printf("\n GOP_counter = %d, sum_mv = %ld\n",GOP_counter,*sum_mv);

/*
  subband_budget = target_bitrate - *sum_mv;
  printf("\nsubband_budget = %d, target_bitrate = %d\n",subband_budget,target_bitrate);

  bpp = (float)subband_budget * 8 / (info.GOPsz * info.yheight * info.ywidth * 1.5);
  pwr_factor = pow(2, (-2)*bpp);
  printf("bpp = %f, pwr_factor = %f\n",bpp,pwr_factor);
*/

  /********/
  /* EZBC */
  /********/
#ifdef CNN_wavelet
  {
	  //// 标记目前编码了几个gop，来表示进度
	  //char name[256], data_file_name[256], start[10];
	  //strncpy(name, info.bitname, strlen(info.bitname) - 4);
	  //name[strlen(info.bitname) - 4] = '_';
	  //snprintf(start, 9, "%03d", info.start);
	  //for (int i = 0; i < 3; i++)
	  //{
		 // name[strlen(info.bitname) - 3 + i] = start[i];
	  //}
	  //name[strlen(info.bitname)] = '\0';               //
	  //sprintf(data_file_name, "%s.encode_number", name);
	  //std::ofstream file(data_file_name, std::ios::app);
	  //file << GOP_counter << std::endl;
	  //file.close();
  }
  {// 写文本文件
	  char name[256], data_file_name[256], start[10];
	  strncpy(name, info.bitname, strlen(info.bitname) - 4);
	  name[strlen(info.bitname) - 4] = '_';
	  snprintf(start, 9, "%03d", info.start);
	  for (int i = 0; i < 3; i++)
	  {
		  name[strlen(info.bitname) - 3 + i] = start[i];
	  }
	  name[strlen(info.bitname)] = '\0';               //
	  sprintf(data_file_name, "%s.data", name);
	  std::cout << data_file_name << std::endl;
	  std::ofstream myfile(data_file_name, std::ios::app);
	  for (int pic = 0; pic < info.GOPsz; pic++)
	  {
		  for (int i = 0; i < info.yheight; i++)
		  {
			  for (int j = 0; j < info.ywidth; j++)
			  {
				  //int temp = round(pyrFrs[pic].Y[i*info.ywidth + j]);
				  float temp = (pyrFrs[pic].Y[i*info.ywidth + j]);
				  myfile << temp << " ";
			  }
			  myfile << std::endl;
		  }
	  }
	  for (int pic = 0; pic < info.GOPsz; pic++)
	  {
		  for (int i = 0; i < info.cheight; i++)
		  {
			  for (int j = 0; j < info.cwidth; j++)
			  {
				  //int temp = round(pyrFrs[pic].U[i*info.cwidth + j]);
				  float temp = (pyrFrs[pic].U[i*info.cwidth + j]);
				  myfile << temp << " ";
			  }
			  myfile << std::endl;
		  }
	  }
	  for (int pic = 0; pic < info.GOPsz; pic++)
	  {
		  for (int i = 0; i < info.cheight; i++)
		  {
			  for (int j = 0; j < info.cwidth; j++)
			  {
				  //int temp = round(pyrFrs[pic].V[i*info.cwidth + j]);
				  float temp = (pyrFrs[pic].V[i*info.cwidth + j]);
				  myfile << temp << " ";
			  }
			  myfile << std::endl;
		  }
	  }
	  myfile.close();
  }
  {
	  //char name[256], data_file_name[256], start[10];
	  //strncpy(name, info.bitname, strlen(info.bitname) - 4);
	  //name[strlen(info.bitname) - 4] = '_';
	  //snprintf(start, 9, "%03d", info.start);
	  //for (int i = 0; i < 3; i++)
	  //{
		 // name[strlen(info.bitname) - 3 + i] = start[i];
	  //}
	  //name[strlen(info.bitname)] = '\0';               //
	  //sprintf(data_file_name, "%s.bindata", name);
	  //std::ofstream myfile(data_file_name, std::ios::binary | std::ios::app);
	  //uint8_t  *buf = new uint8_t[1];
	  //for (int pic = 0; pic < info.GOPsz; pic++)
	  //{
		 // for (int i = 0; i < info.yheight; i++)
		 // {
			//  for (int j = 0; j < info.ywidth; j++)
			//  {
			//	  int temp = round(pyrFrs[pic].Y[i*info.ywidth + j]);
			//	  short int temp1 = temp;
			//	  if (temp != temp1)
			//	  {
			//		  std::cout << temp << " " << temp1 << std::endl;
			//		  std::cout << "cuolw" << std::endl;
			//		  exit(0);
			//	  }
			//	  myfile.write((char*)(&temp1), sizeof(short int));
			//  }
		 // }
	  //}
	  //for (int pic = 0; pic < info.GOPsz; pic++)
	  //{
		 // for (int i = 0; i < info.cheight; i++)
		 // {
			//  for (int j = 0; j < info.cwidth; j++)
			//  {
			//	  int temp = round(pyrFrs[pic].U[i*info.cwidth + j]);
			//	  short int temp1 = temp;
			//	  if (temp != temp1)
			//	  {
			//		  std::cout << temp << " " << temp1 << std::endl;
			//		  std::cout << "cuolu" << std::endl;
			//		  exit(0);
			//	  }
			//	  myfile.write((char*)(&temp1), sizeof(short int));
			//  }
		 // }
	  //}
	  //for (int pic = 0; pic < info.GOPsz; pic++)
	  //{
		 // for (int i = 0; i < info.cheight; i++)
		 // {
			//  for (int j = 0; j < info.cwidth; j++)
			//  {
			//	  int temp = round(pyrFrs[pic].V[i*info.cwidth + j]);
			//	  short int temp1 = temp;
			//	  if (temp != temp1)
			//	  {
			//		  std::cout << temp << " " << temp1 << std::endl;
			//		  std::cout << "cuolv" << std::endl;
			//		  exit(0);
			//	  }
			//	  myfile.write((char*)(&temp1), sizeof(short int));
			//  }
		 // }
	  //}
	  //myfile.close();
  }
  output_GOP_bytes = ezbc3d_enc_GOP( curr, pyrFrs, info, GOP_counter );
  //output_GOP_bytes = 0;
  //exit(0);
#else
#ifdef NO_EZBC
	output_GOP_bytes = 0;
#else
	output_GOP_bytes = ezbc3d_enc_GOP(curr, pyrFrs, info, GOP_counter);
#endif
#endif
  printf("\n GOP_counter = %d, output_GOP_bytes = %ld\n",GOP_counter,output_GOP_bytes);

  return output_GOP_bytes;
}
