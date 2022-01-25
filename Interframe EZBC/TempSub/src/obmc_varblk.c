#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "structN.h"
#include "basic.h"
#include "iostream"
#define EXTERN extern
#include "coderN.h"
#include "util_filtering.h"


// 获得每个块的mv信息
void  get_child_mv_info(vector_ptr fmv2, vector_ptr fmv3, 
						int cx, int cy, int xblk, int yblk, int hor,
						int ver, int *vector_num, 
						ImageMEinfo *imagemeinfo, Varblkarrayinfo *imageblkarray, 
						int connect_info, 
						vector_ptr mv_ref2, vector_ptr mv_ref3, int t_level, videoinfo info) 
{
  int i, j, h_len, v_len, continue_sign=0;
  enum BiMode cur_bi_mode;
  enum spatialMODE  cur_iblock_spatial_mode; 

  int aff_mrg, skip_sign;

  float left_mvx, left_mvy, right_mvx, right_mvy; 
  
  float left_aff1_mvx,left_aff1_mvy,left_aff2_mvx,left_aff2_mvy,left_aff3_mvx,left_aff3_mvy;
  float right_aff1_mvx,right_aff1_mvy,right_aff2_mvx,right_aff2_mvy,right_aff3_mvx,right_aff3_mvy;

/*
  if(fmv2 == NULL)
	  printf("left NULL!\n");
  if(fmv3 == NULL)
	  printf("right NULL!\n");
*/

  switch (connect_info) 
  {
  case 0: 
	  if (fmv2->child && fmv3->child) // 可继续划分
	  {
		  get_child_mv_info(fmv2->child0,  fmv3->child0,  cx,  cy,        xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(fmv2->child1,  fmv3->child1,  cx+xblk/2, cy,  xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(fmv2->child2,  fmv3->child2,  cx,        cy+yblk/2, xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(fmv2->child3,  fmv3->child3,  cx+xblk/2, cy+yblk/2,  xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
	  }else
		  continue_sign = 1; 
	  break; 
  case 1: 
	  if (fmv2->child && mv_ref3==NULL)
	  {
		  get_child_mv_info(fmv2->child0,  NULL,  cx,  cy,        xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(fmv2->child1,  NULL,  cx+xblk/2, cy,  xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(fmv2->child2,  NULL,  cx,        cy+yblk/2, xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(fmv2->child3,  NULL,  cx+xblk/2, cy+yblk/2,  xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
	  }else
		  continue_sign = 1; 
	  break; 
  case 2: 
	  if (mv_ref2==NULL && fmv3->child)
	  {
		  get_child_mv_info(NULL,  fmv3->child0,  cx,  cy,        xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(NULL,  fmv3->child1,  cx+xblk/2, cy,  xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(NULL,  fmv3->child2,  cx,        cy+yblk/2, xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
		  get_child_mv_info(NULL,  fmv3->child3,  cx+xblk/2, cy+yblk/2,  xblk/2,  yblk/2,  hor,
							ver,  vector_num,  imagemeinfo,  imageblkarray, 
							connect_info,  mv_ref2,  mv_ref3, t_level, info);
	  }else
		  continue_sign = 1; 
	  break; 
  }

//  printf("continue_sign = %d\n",continue_sign);

  if (continue_sign) // 不可继续划分了。根据前面的划分模式，记录mv信息
  {
	  if (cx<hor && cy<ver){

		  h_len = (cx+xblk<=hor)?  xblk: hor-cx;
		  v_len = (cy+yblk<=ver)?  yblk: ver-cy; 
		  switch (connect_info)
		  {
		  case 0:  // bi-connected
			  cur_bi_mode = fmv2->bi_mode;
			  aff_mrg     = fmv2->aff_mrg;
			  skip_sign   = fmv2->skip_sign;

			  cur_iblock_spatial_mode = INVALID_SPATIAL_MODE;
			  if (info.bi_mv[t_level]) // bi-directional motion vectors
			  {
				  switch (cur_bi_mode)
				  {
				  case BI_CONNECTED:
				  case BI_PREDICTED:
				  case PARALLEL:
					  left_mvx  = fmv2->mvx;        left_mvy  = fmv2->mvy;
					  right_mvx = fmv3->mvx;        right_mvy = fmv3->mvy;
					  break; 
				  case BLOCK_MERGING:
					  if(fmv2->aff_mrg == YES){
						  printf("bingo BI!\n");
						  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
						  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					      left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

						  right_aff1_mvx = fmv3->aff1_mvx;right_aff1_mvy = fmv3->aff1_mvy;
					      right_aff2_mvx = fmv3->aff2_mvx;right_aff2_mvy = fmv3->aff2_mvy;
					      right_aff3_mvx = fmv3->aff3_mvx;right_aff3_mvy = fmv3->aff3_mvy;  

						  printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);
						  printf("left_aff1_mvx = %f, left_aff1_mvy = %f\nleft_aff2_mvx = %f, left_aff2_mvy = %f\nleft_aff3_mvx = %f, left_aff3_mvy = %f\n",
							 left_aff1_mvx, left_aff1_mvy, left_aff2_mvx, left_aff2_mvy, left_aff3_mvx, left_aff3_mvy );
						  printf("right_aff1_mvx = %f, right_aff1_mvy = %f\nright_aff2_mvx = %f, right_aff2_mvy = %f\nright_aff3_mvx = %f, right_aff3_mvy = %f\n\n",
							 right_aff1_mvx, right_aff1_mvy, right_aff2_mvx, right_aff2_mvy, right_aff3_mvx, right_aff3_mvy );
					  }else{
						  assert(fmv2->aff_mrg == NO);
						  left_mvx  = fmv2->mvx;        left_mvy  = fmv2->mvy;
						  right_mvx = fmv3->mvx;        right_mvy = fmv3->mvy;
					  }
					  break; 
				  case LEFT_CONNECTED:
					  left_mvx = fmv2->mvx;         left_mvy = fmv2->mvy;
					  right_mvx = (float)HUGE_VAL;  right_mvy = (float)HUGE_VAL;
					  break; 
				  case RIGHT_CONNECTED:
					  left_mvx = (float)HUGE_VAL;     left_mvy = (float)HUGE_VAL;
					  right_mvx = fmv3->mvx;          right_mvy = fmv3->mvy;
					  break; 
				  case DIRECTIONAL_IBLOCK:
					  left_mvx = left_mvy = right_mvx = right_mvy = (float)HUGE_VAL;
					  cur_iblock_spatial_mode = fmv2->iblock_spatial_mode;
					  assert(fmv2->iblock_spatial_mode == fmv3->iblock_spatial_mode ); 
					  break;
				//////////  Added by Yuan Liu  //////////
				  case LEFT_CONNECTED_AFF:
					  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
					  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					  left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

					  right_aff1_mvx = (float)HUGE_VAL;right_aff1_mvy = (float)HUGE_VAL;
					  right_aff2_mvx = (float)HUGE_VAL;right_aff2_mvy = (float)HUGE_VAL;
					  right_aff3_mvx = (float)HUGE_VAL;right_aff3_mvy = (float)HUGE_VAL;
					  break;

				case RIGHT_CONNECTED_AFF:
					  right_aff1_mvx = fmv3->aff1_mvx;right_aff1_mvy = fmv3->aff1_mvy;
					  right_aff2_mvx = fmv3->aff2_mvx;right_aff2_mvy = fmv3->aff2_mvy;
					  right_aff3_mvx = fmv3->aff3_mvx;right_aff3_mvy = fmv3->aff3_mvy;

					  left_aff1_mvx = (float)HUGE_VAL;left_aff1_mvy = (float)HUGE_VAL;
					  left_aff2_mvx = (float)HUGE_VAL;left_aff2_mvy = (float)HUGE_VAL;
					  left_aff3_mvx = (float)HUGE_VAL;left_aff3_mvy = (float)HUGE_VAL;
					  break;

				case BI_CONNECTED_AFF:
					  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
					  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					  left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

					  right_aff1_mvx = fmv3->aff1_mvx;right_aff1_mvy = fmv3->aff1_mvy;
					  right_aff2_mvx = fmv3->aff2_mvx;right_aff2_mvy = fmv3->aff2_mvy;
					  right_aff3_mvx = fmv3->aff3_mvx;right_aff3_mvy = fmv3->aff3_mvy;
					  break;
				/////////////////////////////////////////
				  case UNDEFINED:
					  assert(0); 
					  break;
				  default:
					  assert(0);
					  break; 
				  }
			  }

			  if ( !info.bi_mv[t_level] ) // bi-directional motion vectors
			  {
				  assert(0); // 不能进来
				  switch (cur_bi_mode)
				  {
				  case BI_CONNECTED:
				  case BI_PREDICTED:
				  case PARALLEL:
					  left_mvx  = fmv2->mvx;        left_mvy  = fmv2->mvy;
					  right_mvx = fmv3->mvx;        right_mvy = fmv3->mvy;
					  break; 
				  case BLOCK_MERGING:
					  if(fmv2->aff_mrg == YES){
						  printf("bingo!\n");
						  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
						  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					      left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

						  right_aff1_mvx = fmv3->aff1_mvx;right_aff1_mvy = fmv3->aff1_mvy;
					      right_aff2_mvx = fmv3->aff2_mvx;right_aff2_mvy = fmv3->aff2_mvy;
					      right_aff3_mvx = fmv3->aff3_mvx;right_aff3_mvy = fmv3->aff3_mvy;  
					  }else{
						  assert(fmv2->aff_mrg == NO);
						  left_mvx  = fmv2->mvx;        left_mvy  = fmv2->mvy;
						  right_mvx = fmv3->mvx;        right_mvy = fmv3->mvy;
					  }
					  break; 
				  case LEFT_CONNECTED:
					  left_mvx = fmv2->mvx;         left_mvy = fmv2->mvy;
					  right_mvx = (float)HUGE_VAL;  right_mvy = (float)HUGE_VAL;
					  break; 
				  case RIGHT_CONNECTED:   // RIGHT_CONNECTED mode is prevented 
					  assert(0); 
					  break; 
				  case DIRECTIONAL_IBLOCK:
					  left_mvx = left_mvy = right_mvx = right_mvy = (float)HUGE_VAL;
					  cur_iblock_spatial_mode = fmv2->iblock_spatial_mode;
					  assert(fmv2->iblock_spatial_mode == fmv3->iblock_spatial_mode ); 
					  break;

				//////////  Added by Yuan Liu  //////////
				  case LEFT_CONNECTED_AFF:
					  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
					  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					  left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

					  right_aff1_mvx = (float)HUGE_VAL;right_aff1_mvy = (float)HUGE_VAL;
					  right_aff2_mvx = (float)HUGE_VAL;right_aff2_mvy = (float)HUGE_VAL;
					  right_aff3_mvx = (float)HUGE_VAL;right_aff3_mvy = (float)HUGE_VAL;
					  break;

				case RIGHT_CONNECTED_AFF:
					  assert(0);
					  break;

				case BI_CONNECTED_AFF:
					  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
					  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					  left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

					  right_aff1_mvx = fmv3->aff1_mvx;right_aff1_mvy = fmv3->aff1_mvy;
					  right_aff2_mvx = fmv3->aff2_mvx;right_aff2_mvy = fmv3->aff2_mvy;
					  right_aff3_mvx = fmv3->aff3_mvx;right_aff3_mvy = fmv3->aff3_mvy;
					  break;
				  case UNDEFINED:
					  assert(0); 
					  break;
				  default:
					  assert(0);
					  break; 
				  }
			  }
			  break;
		  case 1:  // left connected
			  cur_bi_mode = fmv2->bi_mode;
			  aff_mrg     = fmv2->aff_mrg;
			  skip_sign   = fmv2->skip_sign;

			  cur_iblock_spatial_mode = INVALID_SPATIAL_MODE;
			  switch (cur_bi_mode)
			  {
			  case LEFT_CONNECTED:
			  case LEFT_PREDICTED:
				  left_mvx = fmv2->mvx;         left_mvy = fmv2->mvy;
				  right_mvx = (float)HUGE_VAL;  right_mvy = (float)HUGE_VAL;
				  break; 
			  case BLOCK_MERGING:
					  if(fmv2->aff_mrg == YES){
						  printf("bingo LEFT!\n");
						  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
						  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					      left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

						  right_aff1_mvx = (float)HUGE_VAL;right_aff1_mvy = (float)HUGE_VAL;
					      right_aff2_mvx = (float)HUGE_VAL;right_aff2_mvy = (float)HUGE_VAL;
					      right_aff3_mvx = (float)HUGE_VAL;right_aff3_mvy = (float)HUGE_VAL;  

						  printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);
						  printf("left_aff1_mvx = %f, left_aff1_mvy = %f\nleft_aff2_mvx = %f, left_aff2_mvy = %f\nleft_aff3_mvx = %f, left_aff3_mvy = %f\n",
							 left_aff1_mvx, left_aff1_mvy, left_aff2_mvx, left_aff2_mvy, left_aff3_mvx, left_aff3_mvy );
						  printf("right_aff1_mvx = %f, right_aff1_mvy = %f\nright_aff2_mvx = %f, right_aff2_mvy = %f\nright_aff3_mvx = %f, right_aff3_mvy = %f\n\n",
							 right_aff1_mvx, right_aff1_mvy, right_aff2_mvx, right_aff2_mvy, right_aff3_mvx, right_aff3_mvy );
					  }else{
						  assert(fmv2->aff_mrg == NO);
						  left_mvx  = fmv2->mvx;        left_mvy  = fmv2->mvy;
						  right_mvx = (float)HUGE_VAL;        right_mvy = (float)HUGE_VAL;

//						  printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);
//						  printf("left_mvx = %f, left_mvy = %f\n",left_mvx, left_mvy);
//						  printf("right_mvx = %f, right_mvy = %f\n\n",right_mvx, right_mvy);
					  }
					  break; 
			  case DIRECTIONAL_IBLOCK:
				  left_mvx = left_mvy = right_mvx = right_mvy = (float)HUGE_VAL;
				  cur_iblock_spatial_mode = fmv2->iblock_spatial_mode;
				  break;
			  //////////  Added by Yuan Liu  //////////
				  case LEFT_CONNECTED_AFF:
					  left_aff1_mvx = fmv2->aff1_mvx;left_aff1_mvy = fmv2->aff1_mvy;
					  left_aff2_mvx = fmv2->aff2_mvx;left_aff2_mvy = fmv2->aff2_mvy;
					  left_aff3_mvx = fmv2->aff3_mvx;left_aff3_mvy = fmv2->aff3_mvy;

					  right_aff1_mvx = (float)HUGE_VAL;right_aff1_mvy = (float)HUGE_VAL;
					  right_aff2_mvx = (float)HUGE_VAL;right_aff2_mvy = (float)HUGE_VAL;
					  right_aff3_mvx = (float)HUGE_VAL;right_aff3_mvy = (float)HUGE_VAL;
					  break;
			  /////////////////////////////////////////
			  case UNDEFINED:
				  assert(0); 
				  break;
			  default:
				  assert(0);
				  break; 
			  }
			  break;
		  case 2:  // right connected 
//			  assert(0);
			  cur_bi_mode = fmv3->bi_mode;
			  aff_mrg     = fmv3->aff_mrg;
			  skip_sign   = fmv3->skip_sign;

			  cur_iblock_spatial_mode = INVALID_SPATIAL_MODE;
			  switch (cur_bi_mode)
			  {
			  case LEFT_CONNECTED:
			  case LEFT_PREDICTED:
				  left_mvx = (float)HUGE_VAL;     left_mvy = (float)HUGE_VAL;
				  right_mvx = fmv3->mvx;          right_mvy = fmv3->mvy;
				  break; 
			  case BLOCK_MERGING:
					  if(fmv3->aff_mrg == YES){
						  printf("bingo RIGHT!\n");
						  left_aff1_mvx = (float)HUGE_VAL;left_aff1_mvy = (float)HUGE_VAL;
						  left_aff2_mvx = (float)HUGE_VAL;left_aff2_mvy = (float)HUGE_VAL;
					      left_aff3_mvx = (float)HUGE_VAL;left_aff3_mvy = (float)HUGE_VAL;

						  right_aff1_mvx = fmv3->aff1_mvx;right_aff1_mvy = fmv3->aff1_mvy;
					      right_aff2_mvx = fmv3->aff2_mvx;right_aff2_mvy = fmv3->aff2_mvy;
					      right_aff3_mvx = fmv3->aff3_mvx;right_aff3_mvy = fmv3->aff3_mvy;  

						  printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);
						  printf("left_aff1_mvx = %f, left_aff1_mvy = %f\nleft_aff2_mvx = %f, left_aff2_mvy = %f\nleft_aff3_mvx = %f, left_aff3_mvy = %f\n",
							 left_aff1_mvx, left_aff1_mvy, left_aff2_mvx, left_aff2_mvy, left_aff3_mvx, left_aff3_mvy );
						  printf("right_aff1_mvx = %f, right_aff1_mvy = %f\nright_aff2_mvx = %f, right_aff2_mvy = %f\nright_aff3_mvx = %f, right_aff3_mvy = %f\n\n",
							 right_aff1_mvx, right_aff1_mvy, right_aff2_mvx, right_aff2_mvy, right_aff3_mvx, right_aff3_mvy );
					  }else{
						  assert(fmv3->aff_mrg == NO);
						  left_mvx  = (float)HUGE_VAL;  left_mvy  = (float)HUGE_VAL;
						  right_mvx = fmv3->mvx;        right_mvy = fmv3->mvy;

//						  printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);
//						  printf("left_mvx = %f, left_mvy = %f\n",left_mvx, left_mvy);
//						  printf("right_mvx = %f, right_mvy = %f\n\n",right_mvx, right_mvy);
					  }
					  break;
			  case DIRECTIONAL_IBLOCK:
				  left_mvx = left_mvy = right_mvx = right_mvy = (float)HUGE_VAL;
				  cur_iblock_spatial_mode = fmv3->iblock_spatial_mode;
				  break;
			//////////  Added by Yuan Liu  //////////
			  case LEFT_CONNECTED_AFF:
					  right_aff1_mvx = fmv3->aff1_mvx;right_aff1_mvy = fmv3->aff1_mvy;
					  right_aff2_mvx = fmv3->aff2_mvx;right_aff2_mvy = fmv3->aff2_mvy;
					  right_aff3_mvx = fmv3->aff3_mvx;right_aff3_mvy = fmv3->aff3_mvy;

					  left_aff1_mvx = (float)HUGE_VAL;left_aff1_mvy = (float)HUGE_VAL;
					  left_aff2_mvx = (float)HUGE_VAL;left_aff2_mvy = (float)HUGE_VAL;
					  left_aff3_mvx = (float)HUGE_VAL;left_aff3_mvy = (float)HUGE_VAL;
			  break;
			  /////////////////////////////////////////
			  case UNDEFINED:
				  assert(0); 
				  break;
			  default:
				  assert(0); 
				  break; 
			  }
			  break; 
		  }

		//  /****************查看是否还有connect模式******************/
		//if(fmv2 != NULL)
		//  printf("%d ", fmv2->lifting_mode);
		//if(fmv3 != NULL)
		//  printf("%d ", fmv3->lifting_mode);
		///**********************************/

		  for (i=0; i<h_len; i++)
		  for (j=0; j<v_len; j++)
		  {
			imagemeinfo[(cy+j)*hor+(cx+i)].leftx     = cx;         // imagemeinfo可以作为
			imagemeinfo[(cy+j)*hor+(cx+i)].topy      = cy;
			imagemeinfo[(cy+j)*hor+(cx+i)].blksize   = xblk;	   
			imagemeinfo[(cy+j)*hor+(cx+i)].blocknum  = *vector_num; 
			imagemeinfo[(cy+j)*hor+(cx+i)].bi_mode   = cur_bi_mode;
			imagemeinfo[(cy+j)*hor+(cx+i)].iblock_spatial_mode = cur_iblock_spatial_mode;
			imagemeinfo[(cy+j)*hor+(cx+i)].skip_sign = skip_sign;

			if(cur_bi_mode == BLOCK_MERGING)
				imagemeinfo[(cy+j)*hor+(cx+i)].aff_mrg = aff_mrg;

			if(cur_bi_mode <= 6 || cur_bi_mode == 8 || (cur_bi_mode == 7 && aff_mrg == NO) ){
				imagemeinfo[(cy+j)*hor+(cx+i)].left_mvx  = left_mvx;
				imagemeinfo[(cy+j)*hor+(cx+i)].left_mvy  = left_mvy;
				imagemeinfo[(cy+j)*hor+(cx+i)].right_mvx = right_mvx;
				imagemeinfo[(cy+j)*hor+(cx+i)].right_mvy = right_mvy;
			}else if( (cur_bi_mode>=9 && cur_bi_mode<=11) || (cur_bi_mode == 7 && aff_mrg == YES) ){
				imagemeinfo[(cy+j)*hor+(cx+i)].left_mvx  = (left_aff2_mvx - left_aff1_mvx)*i/h_len + (left_aff3_mvx - left_aff1_mvx)*j/v_len + left_aff1_mvx;
				imagemeinfo[(cy+j)*hor+(cx+i)].left_mvy  = (left_aff2_mvy - left_aff1_mvy)*i/h_len + (left_aff3_mvy - left_aff1_mvy)*j/v_len + left_aff1_mvy;
				imagemeinfo[(cy+j)*hor+(cx+i)].right_mvx = (right_aff2_mvx - right_aff1_mvx)*i/h_len + (right_aff3_mvx - right_aff1_mvx)*j/v_len + right_aff1_mvx;
				imagemeinfo[(cy+j)*hor+(cx+i)].right_mvy = (right_aff2_mvy - right_aff1_mvy)*i/h_len + (right_aff3_mvy - right_aff1_mvy)*j/v_len + right_aff1_mvy;
			
				if( left_aff1_mvx == (float)HUGE_VAL || left_aff1_mvy == (float)HUGE_VAL ||
				  left_aff2_mvx == (float)HUGE_VAL || left_aff2_mvy == (float)HUGE_VAL ||
				  left_aff3_mvx == (float)HUGE_VAL || left_aff3_mvy == (float)HUGE_VAL){
					imagemeinfo[(cy+j)*hor+(cx+i)].left_mvx = (float)HUGE_VAL;
					imagemeinfo[(cy+j)*hor+(cx+i)].left_mvy = (float)HUGE_VAL;
				}
				if( right_aff1_mvx == (float)HUGE_VAL || right_aff1_mvy == (float)HUGE_VAL ||
				  right_aff2_mvx == (float)HUGE_VAL || right_aff2_mvy == (float)HUGE_VAL ||
				  right_aff3_mvx == (float)HUGE_VAL || right_aff3_mvy == (float)HUGE_VAL){
					imagemeinfo[(cy+j)*hor+(cx+i)].right_mvx = (float)HUGE_VAL;
					imagemeinfo[(cy+j)*hor+(cx+i)].right_mvy = (float)HUGE_VAL;
				}
			}

			if( (cur_bi_mode == 7 && aff_mrg == YES) || cur_bi_mode >= 9)
				imagemeinfo[(cy+j)*hor+(cx+i)].mk_sign = 1;
			else
				imagemeinfo[(cy+j)*hor+(cx+i)].mk_sign = 0;

			if( ( (cur_bi_mode == 7 && aff_mrg == YES) || cur_bi_mode >= 0) && (i == 0 || i == h_len-1 || j == 0 || j == v_len-1) )
				imagemeinfo[(cy+j)*hor+(cx+i)].mrg_blk_bdr = 1;
			else
				imagemeinfo[(cy+j)*hor+(cx+i)].mrg_blk_bdr = 0;


			if( (fmv2 != NULL && fmv2->two_comp_src == 1) || (fmv3 != NULL && fmv3->two_comp_src == 1) )
				imagemeinfo[(cy+j)*hor+(cx+i)].two_comp_src = 1;
			else
				imagemeinfo[(cy+j)*hor+(cx+i)].two_comp_src = 0;

			if( ( (fmv2 != NULL && fmv2->two_comp_src == 1) || (fmv3 != NULL && fmv3->two_comp_src == 1) ) && (i == 0 || i == h_len-1 || j == 0 || j == v_len-1) )
				imagemeinfo[(cy+j)*hor+(cx+i)].two_comp_bdr = 1;
			else
				imagemeinfo[(cy+j)*hor+(cx+i)].two_comp_bdr = 0;

		  }

 		 // side information for each block
		 imageblkarray[*vector_num].leftx        = cx;
         imageblkarray[*vector_num].topy         = cy;
		 imageblkarray[*vector_num].blksize      = xblk; // here we should use xblk
		 (*vector_num)++;  // the number of blocks
	}
  }
}



// get the block motion/mode information for each pixel and block 为每个像素和块获取块的运动/模式信息
void get_mv_side_information(videoinfo info, vector_ptr fmv2, vector_ptr fmv3, 
							 vector_ptr mv_ref2, vector_ptr mv_ref3, 
							 ImageMEinfo *imagemeinfo, Varblkarrayinfo *imageblkarray, 
							 int *total_blk, int t_level, int j, 
							 enum FLAG left_scene, enum FLAG right_scene, int encoder_sign)
{
  int  X, x, Y, y;
  int  vector_num, xnum, ynum, xblk, yblk;
  int  hor, ver, pos;
  int  connect_info;
  int s_level;
  FILE *pc;

//  printf("left_scene = %d, right_scene = %d\n",left_scene,right_scene);

//  pc = fopen("enc_block.txt","at");

  *total_blk = 0; 
  if (encoder_sign)  //in encoder use "left_scene" and "right_scene" for connection checking
  {
	  if (left_scene == NO  && right_scene == NO )  // 两侧场景都没有改变，则双向连接
		  connect_info = 0; // bi-connected 双向连接
	  if (left_scene == NO  && right_scene == YES)
		  connect_info = 1; // left connected
	  if (left_scene == YES && right_scene == NO )
		  connect_info = 2; // right connected
	  if (left_scene == YES && right_scene == YES)
	  {
		  printf("total_block = %d\n", 0);
//			  fprintf(pc,"%d\n",0);
		  return;           // intra frame
	  }
  }
  else   // in decoder use "mv_ref2" and "mv_ref3" for connection checking  
  {
	  if (mv_ref2 != NULL && mv_ref3 != NULL )
		  connect_info = 0; // bi-connected
	  if (mv_ref2 != NULL && mv_ref3 == NULL )
		  connect_info = 1; // left connected
	  if (mv_ref2 == NULL && mv_ref3 != NULL )
		  connect_info = 2; // right connected
	  if (mv_ref2 == NULL && mv_ref3 == NULL )
	  {
		  printf("total_block = %d\n", 0);
//			  fprintf(pc,"%d\n",0);
		  return;           // intra frame
	  }
  }

  // the spatial resolution reduction in decoder 
  s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));// 目前info.s_level = 0
  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];
  hor = info.ywidth  << s_level;
  ver = info.yheight << s_level;
  vector_num = 0;
  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {     
	  for( x = 0, X = 0; X < xnum; x += xblk, X++ ) { // 块为单位进行
		pos = Y * xnum + X;
		if (encoder_sign)
			get_child_mv_info( ( left_scene  == NO ) ? (&fmv2[pos]):NULL, 
							   ( right_scene == NO ) ? (&fmv3[pos]):NULL, 
							  x, y, xblk, yblk, hor,
							  ver, &vector_num, 
							  imagemeinfo, imageblkarray, connect_info, 
							  ( left_scene  == NO ) ? mv_ref2 : NULL,  
							  ( right_scene == NO ) ? mv_ref3 : NULL,  
							  t_level, info); 
		else
			get_child_mv_info( mv_ref2 ? (&fmv2[pos]):NULL, 
							   mv_ref3 ? (&fmv3[pos]):NULL, 
							  x, y, xblk, yblk, hor,
							  ver, &vector_num, 
							  imagemeinfo, imageblkarray, connect_info, 
							  mv_ref2, mv_ref3,  t_level, info); 
	}
  }

  *total_blk  = vector_num;
  printf("total_block = %d\n", *total_blk); // 返回的是这一帧的信息
  if(encoder_sign == YES && t_level == 0){
	  gop_block += *total_blk;
	  frame_cnt ++;
  }
//	fprintf(pc,"%d\n",*total_blk);

//  fclose(pc);
}

// 复制预定义的权重矩阵
void copy_weight_matrix(int blksize_cur, enum BiMode cur_bi_mode, 
			            float **self_weight_matrix, float **ver_weight_matrix, float **hor_weight_matrix)
{
	int i, j; 

	switch (blksize_cur) { 
	case 4:
		if (cur_bi_mode != DIRECTIONAL_IBLOCK)
		{
			for (i=0; i<blksize_cur; i++)
				for (j=0; j<blksize_cur; j++)
				{
					self_weight_matrix[i][j] = blk4_cur[i][j]; 
					if (i%2==0)
						ver_weight_matrix[i/2][j]  = blk4_ver[i/2][j];
					if (j%2==0)
						hor_weight_matrix[i][j/2]  = blk4_hor[i][j/2];
				}
		}else
		{
			for (i=0; i<blksize_cur; i++)
				for (j=0; j<blksize_cur; j++)
				{
					self_weight_matrix[i][j] = blk4_cur_iblk[i][j]; 
					if (i%2==0)
						ver_weight_matrix[i/2][j]  = blk4_ver_iblk[i/2][j];
					if (j%2==0)
						hor_weight_matrix[i][j/2]  = blk4_hor_iblk[i][j/2];
				}
		}
		break;
	case 8:
		if (cur_bi_mode != DIRECTIONAL_IBLOCK)
		{
			for (i=0; i<blksize_cur; i++)
				for (j=0; j<blksize_cur; j++)
				{
					self_weight_matrix[i][j] = blk8_cur[i][j]; 
					if (i%2==0)
						ver_weight_matrix[i/2][j]  = blk8_ver[i/2][j];
					if (j%2==0)
						hor_weight_matrix[i][j/2]  = blk8_hor[i][j/2];
				}
		}else
		{
			for (i=0; i<blksize_cur; i++)
				for (j=0; j<blksize_cur; j++)
				{
					self_weight_matrix[i][j] = blk8_cur_iblk[i][j]; 
					if (i%2==0)
						ver_weight_matrix[i/2][j]  = blk8_ver_iblk[i/2][j];
					if (j%2==0)
						hor_weight_matrix[i][j/2]  = blk8_hor_iblk[i][j/2];
				}
		}
		break;
	case 16: 
		for (i=0; i<blksize_cur; i++)
			for (j=0; j<blksize_cur; j++)
			{
				self_weight_matrix[i][j] = blk16_cur[i][j]; 
				if (i%2==0)
					ver_weight_matrix[i/2][j]  = blk16_ver[i/2][j];
				if (j%2==0)
					hor_weight_matrix[i][j/2]  = blk16_hor[i][j/2];
			}
		break;
	case 32:
		for (i=0; i<blksize_cur; i++)
			for (j=0; j<blksize_cur; j++)
			{
				self_weight_matrix[i][j] = blk32_cur[i][j]; 
				if (i%2==0)
					ver_weight_matrix[i/2][j]  = blk32_ver[i/2][j];
				if (j%2==0)
					hor_weight_matrix[i][j/2]  = blk32_hor[i][j/2];
			}
		break;
	case 64:
		for (i=0; i<blksize_cur; i++)
			for (j=0; j<blksize_cur; j++)
			{
				self_weight_matrix[i][j] = blk64_cur[i][j]; 
				if (i%2==0)
					ver_weight_matrix[i/2][j]  = blk64_ver[i/2][j];
				if (j%2==0)
					hor_weight_matrix[i][j/2]  = blk64_hor[i][j/2];
			}
		break;
	case 128:
		for (i=0; i<blksize_cur; i++)
			for (j=0; j<blksize_cur; j++)
			{
				self_weight_matrix[i][j] = blk128_cur[i][j]; 
				if (i%2==0)
					ver_weight_matrix[i/2][j]  = blk128_ver[i/2][j];
				if (j%2==0)
					hor_weight_matrix[i][j/2]  = blk128_hor[i][j/2];
			}
		break;
	case 256:
		for (i=0; i<blksize_cur; i++)
			for (j=0; j<blksize_cur; j++)
			{
				self_weight_matrix[i][j] = blk256_cur[i][j]; 
				if (i%2==0)
					ver_weight_matrix[i/2][j]  = blk256_ver[i/2][j];
				if (j%2==0)
					hor_weight_matrix[i][j/2]  = blk256_hor[i][j/2];
			}
		break;
	default:
		assert(0);
		printf("Error in mv_weight_info()!\n"); 
		break;
	}
}
       

// form the weights from self block, left neighbor and right neighbor 从自身块、左邻居、右邻居形成权重
void mv_weight_info(videoinfo info, ImageMEinfo *pixelmeinfo, Varblkarrayinfo *blockinfoarray,
					int total_varblk, int t_level, int tindex)
{
	int	  yhor, yver, blknum;
	int	  cx, cy, xblk, yblk, xblock, yblock, neighbor_effect;	
	int   row, col, i, j, blksize_cur, pos; 
	float top_weight, bottom_weight, left_weight, right_weight;
	enum BiMode cur_bi_mode; 
	float **self_weight_matrix, **ver_weight_matrix, **hor_weight_matrix;
	int  s_level; 

	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
	// we form the weighting coefficients for the full resolution frame ( original resolution ) 
	// for half, 1/4 resolution we just use the sub-sampled weighting coefficients
	yhor = info.ywidth<<s_level;   
	yver = info.yheight<<s_level;
	for (blknum=0; blknum<total_varblk; blknum++)
	{
		cx = blockinfoarray[blknum].leftx;    // starting position (cx, cy)
		cy = blockinfoarray[blknum].topy;
		xblk = blockinfoarray[blknum].blksize;
		yblk = blockinfoarray[blknum].blksize;
		xblock = (cx+xblk<=yhor)?  xblk : yhor-cx;  // the true block size xblock, yblock
		yblock = (cy+yblk<=yver)?  yblk : yver-cy;

		blksize_cur = pixelmeinfo[cy*yhor+cx].blksize;
		cur_bi_mode = pixelmeinfo[cy*yhor+cx].bi_mode;

		if(xblock<=0 || yblock<=0)	continue;

	  	self_weight_matrix = (float**) getarray(blksize_cur, sizeof(float*), "self_weight_matrix");
		for (i=0; i<blksize_cur; i++)
			self_weight_matrix[i]  = (float*) getarray(blksize_cur, sizeof(float), "self_weight_matrix[]");
		ver_weight_matrix  = (float**) getarray(blksize_cur/2, sizeof(float*), "ver_weight_matrix");
		for (i=0; i<blksize_cur/2; i++)
			ver_weight_matrix[i]   = (float*) getarray(blksize_cur, sizeof(float),  "ver_weight_matrix[]");
		hor_weight_matrix  = (float**) getarray(blksize_cur, sizeof(float*), "hor_weight_matrix");
		for (i=0; i<blksize_cur; i++)
			hor_weight_matrix[i]   = (float*) getarray(blksize_cur/2, sizeof(float), "hor_weight_matrix[]");


		copy_weight_matrix(blksize_cur, cur_bi_mode, self_weight_matrix, 
			ver_weight_matrix, hor_weight_matrix); // 复制预定义的权重到self_weight_matrix, ver_weight_matrix, hor_weight_matrix
		row = cy;
		for(i = 0; i < yblock; i++) {
			col = cx;
			for(j=0 ; j < xblock ; j++){

				// SHRINKING scheme between different block sizes
				// adjust the top_weight, bottom_weight, left_weight, and right_weight for different block sizes

				// top neighbor 
				pos = (cy-1)*yhor + col; // Note: we use col, not cx.
				if( (cy > 0) && pixelmeinfo[pos].bi_mode!=DIRECTIONAL_IBLOCK ) 
				{
					// the block size for top neighbor is equal to current neighbor's
					if (blksize_cur == pixelmeinfo[pos].blksize)
						top_weight = 1.0;
					//  the neighbor block size is smaller
					if (blksize_cur >  pixelmeinfo[pos].blksize)
					{
						// the effective range for the top neighbor
                        neighbor_effect = pixelmeinfo[pos].blksize / 2;
						if (i<neighbor_effect)
							top_weight = 1.0;
						else
							top_weight = 0.0;
					}
					// top neighbor block size is bigger
					if (blksize_cur < pixelmeinfo[pos].blksize)
						top_weight = 1.0;					
				}
				else // reflect the weight to self block
					top_weight = 0.0;
				
				// bottom neighbor 
				pos = (cy+yblk)*yhor + col; // use col instead of cx
				if((cy+yblk < yver) && pixelmeinfo[pos].bi_mode!=DIRECTIONAL_IBLOCK	) 
				{
					if (blksize_cur == pixelmeinfo[pos].blksize)
						bottom_weight = 1.0;
					if (blksize_cur >  pixelmeinfo[pos].blksize)
					{
                        neighbor_effect = pixelmeinfo[pos].blksize/2;
						if (yblk-(i+1)<neighbor_effect)
							bottom_weight = 1.0;
						else
							bottom_weight = 0.0;
					}
					if (blksize_cur < pixelmeinfo[pos].blksize)
						bottom_weight = 1.0;					
				}
				else 
					bottom_weight = 0.0;
				
				// left neighbor 
				pos = row*yhor + cx-1; // note we use row, not cy.
				if((cx > 0) && pixelmeinfo[pos].bi_mode!=DIRECTIONAL_IBLOCK	) 
				{
					if (blksize_cur == pixelmeinfo[pos].blksize)
						left_weight = 1.0;
					if (blksize_cur >  pixelmeinfo[pos].blksize)
					{
                        neighbor_effect = pixelmeinfo[pos].blksize/2;
						if (j<neighbor_effect)
							left_weight = 1.0;
						else
							left_weight = 0.0;
					}
					if (blksize_cur < pixelmeinfo[pos].blksize)
						left_weight = 1.0;					
				}
				else 
					left_weight = 0.0;
				
				// right neighbor 
				pos = row*yhor + cx+xblk; // note we use row, not cy.
				if((cx+xblk < yhor) && pixelmeinfo[pos].bi_mode!=DIRECTIONAL_IBLOCK	) 
				{
					if (blksize_cur == pixelmeinfo[pos].blksize)
						right_weight = 1.0;
					if (blksize_cur >  pixelmeinfo[pos].blksize)
					{
                        neighbor_effect = pixelmeinfo[pos].blksize/2;
						if (xblk-(j+1)<neighbor_effect)
							right_weight = 1.0;
						else
							right_weight = 0.0;
					}
					if (blksize_cur < pixelmeinfo[pos].blksize)
						right_weight = 1.0;					
				}
				else 
					right_weight = 0.0;
				

				// contribution from current block.
				pixelmeinfo[row*yhor+col].self_weight = self_weight_matrix[i][j];
				pixelmeinfo[row*yhor+col].v_weight   = 0;
				pixelmeinfo[row*yhor+col].h_weight   = 0;
				// weights from neighboring blocks.
				if (i < (yblk >> 1)) { 
					// top neighbor is effective, and also affected by shrinking scheme
					pixelmeinfo[row*yhor+col].v_weight     = ver_weight_matrix[i][j]*top_weight;
					// modify self weight (raise self_weight to compensate the top_weight in vertical direction)
					pixelmeinfo[row*yhor+col].self_weight += ver_weight_matrix[i][j]*(1-top_weight);
				}else if (i >= (yblk >> 1)) { 
					// bottom neighbor is effective, whose weights are symmetric to top neighbor's
					pixelmeinfo[row*yhor+col].v_weight     = ver_weight_matrix[yblk - i - 1][j]*bottom_weight;
					pixelmeinfo[row*yhor+col].self_weight += ver_weight_matrix[yblk - i - 1][j]*(1-bottom_weight);
				}
				if( (j < (xblk >> 1)) ) {   // left neighbor
					pixelmeinfo[row*yhor+col].h_weight     = hor_weight_matrix[i][j]*left_weight;
					pixelmeinfo[row*yhor+col].self_weight += hor_weight_matrix[i][j]*(1-left_weight);
				}else if( (j >= (xblk >> 1)) ) { // right
					pixelmeinfo[row*yhor+col].h_weight     = hor_weight_matrix[i][xblk - j - 1]*right_weight;
					pixelmeinfo[row*yhor+col].self_weight += hor_weight_matrix[i][xblk - j - 1]*(1-right_weight);
				}

				pixelmeinfo[row*yhor+col].self_weight /= (float)xblk; // normalize
				pixelmeinfo[row*yhor+col].h_weight    /= (float)xblk; // normalize
				pixelmeinfo[row*yhor+col].v_weight    /= (float)xblk; // normalize
				assert(fabs(pixelmeinfo[row*yhor+col].self_weight+pixelmeinfo[row*yhor+col].h_weight
							+pixelmeinfo[row*yhor+col].v_weight-1.0)< pow((float)10, int(-5) ) );
				col++;
			}
			row++;
		}
	  	
		for (i=0; i<blksize_cur; i++)
			free(self_weight_matrix[i]);
		for (i=0; i<blksize_cur/2; i++)
			free(ver_weight_matrix[i]);
		for (i=0; i<blksize_cur; i++)
			free(hor_weight_matrix[i]);
	  	free(self_weight_matrix); 
		free(ver_weight_matrix);
		free(hor_weight_matrix);

	}
}
