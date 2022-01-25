#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define EXTERN extern

#include "directional_iblock.h"
#include "mode_decision.h"
#include "bmeN.h"
#include "mvcodingN.h"
#include "basic.h"
#include "structN.h"
#include "util_filtering.h"
#include "bme_tools.h"
#include "coderN.h"
#include "mv_ec.h"
#include "miscN.h"

#define  DIRECT_LEN 1
#define  MERGE_LEN  1
#define  MERGE_DIR_UP_LEN 1
#define  MERGE_DIR_LEN  2
#define  AFF_IDX_OFFSET 3
#define  MRG_RANGE 10
#define  ITR_RANGE 6
#define  BI_RANGE  3
#define  PRL_RANGE 4

#define  SMALL_AFF_BLK_COEF  0.90
#define  AFF_BLK_COEF_32  0.94
#define  AFF_BLK_COEF_64  0.94


ModeInfo *tv11;
ModeInfo *tv12;
ModeInfo *tv13;

ModeInfo *tv21;
ModeInfo *tv22;
ModeInfo *tv23;

float    *block_buff1;
float    *block_buff2;

int get_idx_length( int aff_num ){
	int i,j;
	int getnum, val;
	int length = 0;

	if(aff_num >= 9){
		getnum = (aff_num + 1);
		while(getnum >= 4){
			getnum -= 4;
			length += 1;
		}
		length += 1;
		assert(getnum <= 3 && getnum >= 0);
		length += 2;
	}else if(aff_num >= 5 && aff_num <= 8){
		length = 3;
	}else if(aff_num >= 2 && aff_num <= 4){
		length = 2;
	}else if(aff_num == 1){
		length = 1;
	}else
		assert(0);

	return length;
}


//Added on 08.15.2016	////////////////////////////

int compare_mv(vector_ptr left1, vector_ptr right1, vector_ptr left2, vector_ptr right2){
	int eql;
	int diff_l = 1, diff_r = 1;

	assert(left1->bi_mode == right1->bi_mode && left2->bi_mode == right2->bi_mode);

	if(left1->bi_mode >= 9 || left2->bi_mode >= 9){
		diff_l = 1;
		diff_r = 1;
	}else{
		if(left1 != NULL && left2 != NULL){
			if(left1->mvx == left2->mvx && left1->mvy == left2->mvy)
				diff_l = 0;
		}else if(left1 == NULL && right1 == NULL){
			diff_l = 0;
		}

		if(right1 != NULL && right2 != NULL){
			if(right1->mvx == right2->mvx && right1->mvy == right2->mvy)
				diff_r = 0;
		}else if(right1 == NULL && right1 == NULL){
			diff_r = 0;
		}
	}

	return (diff_l + diff_r);
}

//cx,cy - the coordinate of current block's upper-left corner
//x_pos,y_pos the merge reference point
void get_aff_mrg_mv( vector_ptr fmv, int cx, int cy, int x_pos, int y_pos, int xblk, int yblk, float dx1, float dy1, 
					 float dx2, float dy2, float pmvx, float pmvy ){
	int i,j,k;
	float mvx1, mvy1, mvx2, mvy2, mvx3, mvy3;

	mvx1 = pmvx - (x_pos - cx) * dx1 - (y_pos - cy) * dx2;
	mvy1 = pmvy - (x_pos - cx) * dy1 - (y_pos - cy) * dy2;

	mvx2 = pmvx - (x_pos - cx - xblk) * dx1 - (y_pos - cy) * dx2;
	mvy2 = pmvy - (x_pos - cx - xblk) * dy1 - (y_pos - cy) * dy2;

	mvx3 = pmvx - (x_pos - cx) * dx1 - (y_pos - cy - yblk) * dx2;
	mvy3 = pmvy - (x_pos - cx) * dy1 - (y_pos - cy - yblk) * dy2;

	fmv->aff1_mvx = mvx1;
	fmv->aff1_mvy = mvy1;
	fmv->aff2_mvx = mvx2;
	fmv->aff2_mvy = mvy2;
	fmv->aff3_mvx = mvx3;
	fmv->aff3_mvy = mvy3;
}

/*
	CBD
	AZ
	E
*/
void get_merge_mv_info(vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv3_array, vector_ptr fmv4_array, vector_ptr *mrg_left, 
					   vector_ptr *mrg_right, int x_pos, int y_pos, int hor, int ver, int xblk, int yblk, videoinfo info, int t_level){
	int i,j,k;
	int get_xblk, get_yblk;
	vector_ptr left_mva = NULL, left_mvb = NULL, left_mvc = NULL, left_mvd = NULL, left_mv_tmp = NULL, left_mve = NULL;
	vector_ptr right_mva = NULL, right_mvb = NULL, right_mvc = NULL, right_mvd = NULL, right_mv_tmp = NULL, right_mve = NULL;

	float ula11,ula12,ulb11,ulb12,ulc11,ulc12,uld11,uld12,ult11,ult12,ule11,ule12;
	float ula21,ula22,ulb21,ulb22,ulc21,ulc22,uld21,uld22,ult21,ult22,ule21,ule22;

	float left_affx1, left_affy1, left_affx2, left_affy2, left_affx3, left_affy3, left_affx4, left_affy4;
	float right_affx1, right_affy1, right_affx2, right_affy2, right_affx3, right_affy3, right_affx4, right_affy4;
	float aff_dmvx1,aff_dmvx2,aff_dmvy1,aff_dmvy2;

	float med_l_px[4], med_l_py[4], med_r_px[4], med_r_py[4];
	float med_l_aff1x[4],med_l_aff1y[4],med_l_aff2x[4],med_l_aff2y[4],med_l_aff3x[4],med_l_aff3y[4];
	float med_r_aff1x[4],med_r_aff1y[4],med_r_aff2x[4],med_r_aff2y[4],med_r_aff3x[4],med_r_aff3y[4];

	int aff_blk; //Added on 03.14.2017

	for(i=0;i<=3;i++){
		med_l_px[i] = (float)HUGE_VAL;
		med_l_py[i] = (float)HUGE_VAL;
		med_l_aff1x[i] = (float)HUGE_VAL;
		med_l_aff1y[i] = (float)HUGE_VAL;
		med_l_aff2x[i] = (float)HUGE_VAL;
		med_l_aff2y[i] = (float)HUGE_VAL;
		med_l_aff3x[i] = (float)HUGE_VAL;
		med_l_aff3y[i] = (float)HUGE_VAL;

		med_r_px[i] = (float)HUGE_VAL;
		med_r_py[i] = (float)HUGE_VAL;
		med_r_aff1x[i] = (float)HUGE_VAL;
		med_r_aff1y[i] = (float)HUGE_VAL;
		med_r_aff2x[i] = (float)HUGE_VAL;
		med_r_aff2y[i] = (float)HUGE_VAL;
		med_r_aff3x[i] = (float)HUGE_VAL;
		med_r_aff3y[i] = (float)HUGE_VAL;
	}
	 
//A
	aff_blk = 0;

	left_mva = find_block(x_pos-1,y_pos+yblk-1,fmv1_array,info,t_level,&ula11,&ula12,&get_xblk,&get_yblk,0,0);

	find_block(x_pos-3,y_pos+yblk-3,fmv1_array,info,t_level,&left_affx1,&left_affy1,&get_xblk,&get_yblk,0,0);
	find_block(x_pos-2,y_pos+yblk-3,fmv1_array,info,t_level,&left_affx2,&left_affy2,&get_xblk,&get_yblk,0,0);
	find_block(x_pos-3,y_pos+yblk-2,fmv1_array,info,t_level,&left_affx3,&left_affy3,&get_xblk,&get_yblk,0,0);

	if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){

		assert(left_mva->bi_mode>=9 || (left_mva->bi_mode == 7 && left_mva->aff_mrg == YES) );

		aff_dmvx1 = (left_affx2 - left_affx1);
		aff_dmvy1 = (left_affy2 - left_affy1);

		aff_dmvx2 = (left_affx3 - left_affx1);
		aff_dmvy2 = (left_affy3 - left_affy1);

		get_aff_mrg_mv(mrg_left[0], x_pos, y_pos, x_pos-3, y_pos+yblk-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);

		aff_blk = 1;
	}

	if(fmv2_array != NULL){
		right_mva = find_block(x_pos-1,y_pos+yblk-1,fmv2_array,info,t_level,&ula21,&ula22,&get_xblk,&get_yblk,0,0);

		find_block(x_pos-3,y_pos+yblk-3,fmv2_array,info,t_level,&right_affx1,&right_affy1,&get_xblk,&get_yblk,0,0);
		find_block(x_pos-2,y_pos+yblk-3,fmv2_array,info,t_level,&right_affx2,&right_affy2,&get_xblk,&get_yblk,0,0);
		find_block(x_pos-3,y_pos+yblk-2,fmv2_array,info,t_level,&right_affx3,&right_affy3,&get_xblk,&get_yblk,0,0);

		if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){

			assert(right_mva->bi_mode>=9 || (right_mva->bi_mode == 7 && right_mva->aff_mrg == YES) );

			aff_dmvx1 = (right_affx2 - right_affx1);
			aff_dmvy1 = (right_affy2 - right_affy1);

			aff_dmvx2 = (right_affx3 - right_affx1);
			aff_dmvy2 = (right_affy3 - right_affy1);

			get_aff_mrg_mv(mrg_right[0], x_pos, y_pos, x_pos-3, y_pos+yblk-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);

			aff_blk = 1;
		}
	
		if(aff_blk == 1){
			if( left_mva != NULL && left_mva->is_predictor == YES ){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_aff_mrg_mv(mrg_left[0], x_pos, y_pos, x_pos-3, y_pos+yblk-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
			}

			if( right_mva != NULL && right_mva->is_predictor == YES ){
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_aff_mrg_mv(mrg_right[0], x_pos, y_pos, x_pos-3, y_pos+yblk-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
			}
		}
	}
//B
	aff_blk = 0;
	left_mvb = find_block(x_pos+xblk-1,y_pos-1,fmv1_array,info,t_level,&ulb11,&ulb12,&get_xblk,&get_yblk,0,0);

	find_block(x_pos+xblk-3,y_pos-3,fmv1_array,info,t_level,&left_affx1,&left_affy1,&get_xblk,&get_yblk,0,0);
	find_block(x_pos+xblk-2,y_pos-3,fmv1_array,info,t_level,&left_affx2,&left_affy2,&get_xblk,&get_yblk,0,0);
	find_block(x_pos+xblk-3,y_pos-2,fmv1_array,info,t_level,&left_affx3,&left_affy3,&get_xblk,&get_yblk,0,0);

	if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){

		assert( left_mvb->bi_mode>=9 || (left_mvb->bi_mode == 7 && left_mvb->aff_mrg == YES) );

		aff_dmvx1 = (left_affx2 - left_affx1);
		aff_dmvy1 = (left_affy2 - left_affy1);

		aff_dmvx2 = (left_affx3 - left_affx1);
		aff_dmvy2 = (left_affy3 - left_affy1);

		get_aff_mrg_mv(mrg_left[1], x_pos, y_pos, x_pos+xblk-3,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);

		aff_blk = 1;
	}

	if(fmv2_array != NULL){
		right_mvb = find_block(x_pos+xblk-1,y_pos-1,fmv2_array,info,t_level,&ulb21,&ulb22,&get_xblk,&get_yblk,0,0);
	
		find_block(x_pos+xblk-3,y_pos-3,fmv2_array,info,t_level,&right_affx1,&right_affy1,&get_xblk,&get_yblk,0,0);
		find_block(x_pos+xblk-2,y_pos-3,fmv2_array,info,t_level,&right_affx2,&right_affy2,&get_xblk,&get_yblk,0,0);
		find_block(x_pos+xblk-3,y_pos-2,fmv2_array,info,t_level,&right_affx3,&right_affy3,&get_xblk,&get_yblk,0,0);

		if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){

			assert(right_mvb->bi_mode>=9 || (right_mvb->bi_mode == 7 && right_mvb->aff_mrg == YES) );

			aff_dmvx1 = (right_affx2 - right_affx1);
			aff_dmvy1 = (right_affy2 - right_affy1);

			aff_dmvx2 = (right_affx3 - right_affx1);
			aff_dmvy2 = (right_affy3 - right_affy1);

			get_aff_mrg_mv(mrg_right[1], x_pos, y_pos,x_pos+xblk-3,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);

			aff_blk = 1;
		}

		if(aff_blk == 1){
			if( left_mvb != NULL && left_mvb->is_predictor == YES ){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_aff_mrg_mv(mrg_left[1], x_pos, y_pos, x_pos+xblk-3,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
			}
			if( right_mvb != NULL && right_mvb->is_predictor == YES ){
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_aff_mrg_mv(mrg_right[1], x_pos, y_pos,x_pos+xblk-3,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
			}
		}

	}
//C 
	aff_blk = 0;
	left_mvc = find_block(x_pos-1,y_pos-1,fmv1_array,info,t_level,&ulc11,&ulc12,&get_xblk,&get_yblk,0,0);

	find_block(x_pos-3,y_pos-3,fmv1_array,info,t_level,&left_affx1,&left_affy1,&get_xblk,&get_yblk,0,0);
	find_block(x_pos-2,y_pos-3,fmv1_array,info,t_level,&left_affx2,&left_affy2,&get_xblk,&get_yblk,0,0);
	find_block(x_pos-3,y_pos-2,fmv1_array,info,t_level,&left_affx3,&left_affy3,&get_xblk,&get_yblk,0,0);

	if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){

		assert(left_mvc->bi_mode>=9 || (left_mvc->bi_mode == 7 && left_mvc->aff_mrg == YES) );

		aff_dmvx1 = (left_affx2 - left_affx1);
		aff_dmvy1 = (left_affy2 - left_affy1);

		aff_dmvx2 = (left_affx3 - left_affx1);
		aff_dmvy2 = (left_affy3 - left_affy1);

		get_aff_mrg_mv(mrg_left[2], x_pos, y_pos, x_pos-3, y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);

		aff_blk = 1;
	}

	if(fmv2_array != NULL){
		right_mvc = find_block(x_pos-1,y_pos-1,fmv2_array,info,t_level,&ulc21,&ulc22,&get_xblk,&get_yblk,0,0);
	
		find_block(x_pos-3,y_pos-3,fmv2_array,info,t_level,&right_affx1,&right_affy1,&get_xblk,&get_yblk,0,0);
		find_block(x_pos-2,y_pos-3,fmv2_array,info,t_level,&right_affx2,&right_affy2,&get_xblk,&get_yblk,0,0);
		find_block(x_pos-3,y_pos-2,fmv2_array,info,t_level,&right_affx3,&right_affy3,&get_xblk,&get_yblk,0,0);

		if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){

			assert(right_mvc->bi_mode>=9 || (right_mvc->bi_mode == 7 && right_mvc->aff_mrg == YES) );

			aff_dmvx1 = (right_affx2 - right_affx1);
			aff_dmvy1 = (right_affy2 - right_affy1);

			aff_dmvx2 = (right_affx3 - right_affx1);
			aff_dmvy2 = (right_affy3 - right_affy1);

			get_aff_mrg_mv(mrg_right[2], x_pos, y_pos, x_pos-3, y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);

			aff_blk = 1;
		}

		if(aff_blk == 1){
			if(left_mvc != NULL && left_mvc->is_predictor == YES ){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_aff_mrg_mv(mrg_left[2], x_pos, y_pos, x_pos-3, y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
			}
			if(right_mvc != NULL && right_mvc->is_predictor == YES ){
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_aff_mrg_mv(mrg_right[2], x_pos, y_pos, x_pos-3, y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
			}
		}
	}

	aff_blk = 0;

	if( (left_mva != NULL || right_mva != NULL) && (left_mvb != NULL || right_mvb != NULL) && (left_mvc != NULL || right_mvc != NULL) );
	else{
//D
		left_mvd = find_block(x_pos+xblk,y_pos-1,fmv1_array,info,t_level,&uld11,&uld12,&get_xblk,&get_yblk,0,0);

		find_block(x_pos+xblk+1,y_pos-3,fmv1_array,info,t_level,&left_affx1,&left_affy1,&get_xblk,&get_yblk,0,0);
		find_block(x_pos+xblk+2,y_pos-3,fmv1_array,info,t_level,&left_affx2,&left_affy2,&get_xblk,&get_yblk,0,0);
		find_block(x_pos+xblk+1,y_pos-2,fmv1_array,info,t_level,&left_affx3,&left_affy3,&get_xblk,&get_yblk,0,0);

		if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){

			assert(left_mvd->bi_mode>=9 || (left_mvd->bi_mode == 7 && left_mvd->aff_mrg == YES) );

			aff_dmvx1 = (left_affx2 - left_affx1);
			aff_dmvy1 = (left_affy2 - left_affy1);

			aff_dmvx2 = (left_affx3 - left_affx1);
			aff_dmvy2 = (left_affy3 - left_affy1);

			get_aff_mrg_mv(mrg_left[3], x_pos, y_pos, x_pos+xblk+1,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
			aff_blk = 1;
		}

		if(fmv2_array != NULL){
			right_mvd = find_block(x_pos+xblk,y_pos-1,fmv2_array,info,t_level,&uld21,&uld22,&get_xblk,&get_yblk,0,0);
		    
			find_block(x_pos+xblk+1,y_pos-3,fmv2_array,info,t_level,&right_affx1,&right_affy1,&get_xblk,&get_yblk,0,0);
			find_block(x_pos+xblk+2,y_pos-3,fmv2_array,info,t_level,&right_affx2,&right_affy2,&get_xblk,&get_yblk,0,0);
			find_block(x_pos+xblk+1,y_pos-2,fmv2_array,info,t_level,&right_affx3,&right_affy3,&get_xblk,&get_yblk,0,0);

			if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){

				assert(right_mvd->bi_mode>=9 || (right_mvd->bi_mode == 7 && right_mvd->aff_mrg == YES) );

				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_aff_mrg_mv(mrg_right[3], x_pos, y_pos, x_pos+xblk+1,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
				aff_blk = 1;
			}
		
			if(aff_blk == 1){
				if(left_mvd != NULL && left_mvd->is_predictor == YES ){
					aff_dmvx1 = (left_affx2 - left_affx1);
					aff_dmvy1 = (left_affy2 - left_affy1);

					aff_dmvx2 = (left_affx3 - left_affx1);
					aff_dmvy2 = (left_affy3 - left_affy1);

					get_aff_mrg_mv(mrg_left[3], x_pos, y_pos, x_pos+xblk+1,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				}
				if(right_mvd != NULL && right_mvd->is_predictor == YES ){
					aff_dmvx1 = (right_affx2 - right_affx1);
					aff_dmvy1 = (right_affy2 - right_affy1);

					aff_dmvx2 = (right_affx3 - right_affx1);
					aff_dmvy2 = (right_affy3 - right_affy1);

					get_aff_mrg_mv(mrg_right[3], x_pos, y_pos, x_pos+xblk+1,y_pos-3,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
				}
			}
		}
	}

	if(left_mva != NULL){
		if( left_mva->bi_mode <= 6 || left_mva->bi_mode == 8 || (left_mva->bi_mode == 7 && left_mva->aff_mrg == NO) ){
			mrg_left[0]->mvx = left_mva->mvx;	mrg_left[0]->mvy = left_mva->mvy;	mrg_left[0]->lifting_mode = CONNECTED;
		}else{
			mrg_left[0]->mvx = ula11;	mrg_left[0]->mvy = ula12;	mrg_left[0]->lifting_mode = CONNECTED;
		}
	}
	if(right_mva != NULL){
		if( right_mva->bi_mode <= 6 || right_mva->bi_mode == 8 || (right_mva->bi_mode == 7 && right_mva->aff_mrg == NO) ){
			mrg_right[0]->mvx = right_mva->mvx;	mrg_right[0]->mvy = right_mva->mvy;	mrg_right[0]->lifting_mode = CONNECTED;
		}else{
			mrg_right[0]->mvx = ula21;	mrg_right[0]->mvy = ula22;	mrg_right[0]->lifting_mode = CONNECTED;
		}
	}
	if(left_mvb != NULL){
		if(left_mvb->bi_mode <= 6 || left_mvb->bi_mode == 8 || (left_mvb->bi_mode == 7 && left_mvb->aff_mrg == NO) ){
			mrg_left[1]->mvx = left_mvb->mvx;	mrg_left[1]->mvy = left_mvb->mvy;	mrg_left[1]->lifting_mode = CONNECTED;
		}else{
			mrg_left[1]->mvx = ulb11;	mrg_left[1]->mvy = ulb12;	mrg_left[1]->lifting_mode = CONNECTED;
		}
	}
	if(right_mvb != NULL){
		if( right_mvb->bi_mode <= 6 || right_mvb->bi_mode == 8 || (right_mvb->bi_mode == 7 && right_mvb->aff_mrg == NO) ){
			mrg_right[1]->mvx = right_mvb->mvx;	mrg_right[1]->mvy = right_mvb->mvy;	mrg_right[1]->lifting_mode = CONNECTED;
		}else{
			mrg_right[1]->mvx = ulb21;	mrg_right[1]->mvy = ulb22;	mrg_right[1]->lifting_mode = CONNECTED;
		}
	}
	if(left_mvc != NULL){
		if( left_mvc->bi_mode <= 6 || left_mvc->bi_mode == 8 || (left_mvc->bi_mode == 7 && left_mvc->aff_mrg == NO) ){
			mrg_left[2]->mvx = left_mvc->mvx;	mrg_left[2]->mvy = left_mvc->mvy;	mrg_left[2]->lifting_mode = CONNECTED;
		}else{
			mrg_left[2]->mvx = ulc11;	mrg_left[2]->mvy = ulc12;	mrg_left[2]->lifting_mode = CONNECTED;
		}
	}
	if(right_mvc != NULL){
		if( right_mvc->bi_mode <= 6 || right_mvc->bi_mode == 8 || (right_mvc->bi_mode == 7 && right_mvc->aff_mrg == NO) ){
			mrg_right[2]->mvx = right_mvc->mvx;	mrg_right[2]->mvy = right_mvc->mvy;	mrg_right[2]->lifting_mode = CONNECTED;
		}else{
			mrg_right[2]->mvx = ulc21;	mrg_right[2]->mvy = ulc22;	mrg_right[2]->lifting_mode = CONNECTED;
		}
	}
	if(left_mvd != NULL){
		if( left_mvd->bi_mode <= 6 || left_mvd->bi_mode == 8 || (left_mvd->bi_mode == 7 && left_mvd->aff_mrg == NO) ){
			mrg_left[3]->mvx = left_mvd->mvx;	mrg_left[3]->mvy = left_mvd->mvy;	mrg_left[3]->lifting_mode = CONNECTED;
		}else{
			mrg_left[3]->mvx = uld11;	mrg_left[3]->mvy = uld12;	mrg_left[3]->lifting_mode = CONNECTED;
		}
	}
	if(right_mvd != NULL){
		if( right_mvd->bi_mode <= 6 || right_mvd->bi_mode == 8 || (right_mvd->bi_mode == 7 && right_mvd->aff_mrg == NO) ){
			mrg_right[3]->mvx = right_mvd->mvx;	mrg_right[3]->mvy = right_mvd->mvy;	mrg_right[3]->lifting_mode = CONNECTED;
		}else{
			mrg_right[3]->mvx = uld21;	mrg_right[3]->mvy = uld22;	mrg_right[3]->lifting_mode = CONNECTED;
		}
	}

	////////////////////////
	for(i = 0;i <= 3;i ++){
	  if( (x_pos-(mrg_left[i]->mvx)<0||x_pos-(mrg_left[i]->mvx)+xblk>hor||y_pos-(mrg_left[i]->mvy)<0|| 
	       y_pos-(mrg_left[i]->mvy)+yblk>ver) && (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL ) ){
		mrg_left[i]->mvx = (float)HUGE_VAL;
        mrg_left[i]->mvy = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

		mrg_right[i]->mvx = (float)HUGE_VAL;
        mrg_right[i]->mvy = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
	  }
	  if( (x_pos-(mrg_right[i]->mvx)<0||x_pos-(mrg_right[i]->mvx)+xblk>hor||y_pos-(mrg_right[i]->mvy)<0||
		  y_pos-(mrg_right[i]->mvy)+yblk>ver)  && (mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL ) ){
		mrg_left[i]->mvx = (float)HUGE_VAL;
        mrg_left[i]->mvy = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

		mrg_right[i]->mvx = (float)HUGE_VAL;
        mrg_right[i]->mvy = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
	  }
	}

	//check for repeated candidates
	for(i = 3;i >= 1;i --){
		for(k = i-1;k >= 0;k --){
			if( (mrg_left[i]->mvx == mrg_left[k]->mvx && mrg_left[i]->mvy == mrg_left[k]->mvy) &&
			(mrg_right[i]->mvx == mrg_right[k]->mvx && mrg_right[i]->mvy == mrg_right[k]->mvy) ){
				mrg_left[i]->mvx = (float)HUGE_VAL;
				mrg_left[i]->mvy = (float)HUGE_VAL;
				mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
				mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
				mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

				mrg_right[i]->mvx = (float)HUGE_VAL;
				mrg_right[i]->mvy = (float)HUGE_VAL;
				mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
				mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
				mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
			}
		}
	}
	//check for HUGE_VALs
    k = 0;
    for(i = 0;i <= 3;i++){
	    if( (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL) ||
			(mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL) ){
		  med_l_px[k] = mrg_left[i]->mvx;
		  med_l_py[k] = mrg_left[i]->mvy;
		  med_l_aff1x[k] = mrg_left[i]->aff1_mvx;
		  med_l_aff1y[k] = mrg_left[i]->aff1_mvy;
		  med_l_aff2x[k] = mrg_left[i]->aff2_mvx;
		  med_l_aff2y[k] = mrg_left[i]->aff2_mvy;
		  med_l_aff3x[k] = mrg_left[i]->aff3_mvx;
		  med_l_aff3y[k] = mrg_left[i]->aff3_mvy;

		  med_r_px[k] = mrg_right[i]->mvx;
		  med_r_py[k] = mrg_right[i]->mvy;
		  med_r_aff1x[k] = mrg_right[i]->aff1_mvx;
		  med_r_aff1y[k] = mrg_right[i]->aff1_mvy;
		  med_r_aff2x[k] = mrg_right[i]->aff2_mvx;
		  med_r_aff2y[k] = mrg_right[i]->aff2_mvy;
		  med_r_aff3x[k] = mrg_right[i]->aff3_mvx;
		  med_r_aff3y[k] = mrg_right[i]->aff3_mvy;

		  k++;
	    }
    }
	assert(k <= 3);
	for(i = 0;i <= 3;i++){
		mrg_left[i]->mvx = med_l_px[i];
		mrg_left[i]->mvy = med_l_py[i];
		mrg_left[i]->aff1_mvx = med_l_aff1x[i];
		mrg_left[i]->aff1_mvy = med_l_aff1y[i];
		mrg_left[i]->aff2_mvx = med_l_aff2x[i];
		mrg_left[i]->aff2_mvy = med_l_aff2y[i];
		mrg_left[i]->aff3_mvx = med_l_aff3x[i];
		mrg_left[i]->aff3_mvy = med_l_aff3y[i];

		mrg_right[i]->mvx  = med_r_px[i];
		mrg_right[i]->mvy  = med_r_py[i];
		mrg_right[i]->aff1_mvx = med_r_aff1x[i];
		mrg_right[i]->aff1_mvy = med_r_aff1y[i];
		mrg_right[i]->aff2_mvx = med_r_aff2x[i];
		mrg_right[i]->aff2_mvy = med_r_aff2y[i];
		mrg_right[i]->aff3_mvx = med_r_aff3x[i];
		mrg_right[i]->aff3_mvy = med_r_aff3y[i];

	}
	for(i=0;i<=3;i++){
		med_l_px[i] = (float)HUGE_VAL;
		med_l_py[i] = (float)HUGE_VAL;
		med_l_aff1x[i] = (float)HUGE_VAL;
		med_l_aff1y[i] = (float)HUGE_VAL;
		med_l_aff2x[i] = (float)HUGE_VAL;
		med_l_aff2y[i] = (float)HUGE_VAL;
		med_l_aff3x[i] = (float)HUGE_VAL;
		med_l_aff3y[i] = (float)HUGE_VAL;

		med_r_px[i] = (float)HUGE_VAL;
		med_r_py[i] = (float)HUGE_VAL;
		med_r_aff1x[i] = (float)HUGE_VAL;
		med_r_aff1y[i] = (float)HUGE_VAL;
		med_r_aff2x[i] = (float)HUGE_VAL;
		med_r_aff2y[i] = (float)HUGE_VAL;
		med_r_aff3x[i] = (float)HUGE_VAL;
		med_r_aff3y[i] = (float)HUGE_VAL;
	}

	assert(mrg_left[3]->mvx == (float)HUGE_VAL && mrg_left[3]->mvy == (float)HUGE_VAL
		&& mrg_right[3]->mvx == (float)HUGE_VAL && mrg_right[3]->mvy == (float)HUGE_VAL);
	
//TMP merge candidate
	aff_blk = 0;

    for(i = 0; i <= 3; i++){
		if( (mrg_left[i]->mvx == (float)HUGE_VAL && mrg_left[i]->mvy == (float)HUGE_VAL) &&
			(mrg_right[i]->mvx == (float)HUGE_VAL && mrg_right[i]->mvy == (float)HUGE_VAL) ){
		  if(fmv3_array != NULL)
			left_mv_tmp = find_block(x_pos,y_pos,fmv3_array,info,t_level,&ult11,&ult12,&get_xblk,&get_yblk,0,0);

		  if(fmv2_array != NULL && fmv4_array != NULL)
			right_mv_tmp = find_block(x_pos,y_pos,fmv4_array,info,t_level,&ult21,&ult22,&get_xblk,&get_yblk,0,0);

		  if(left_mv_tmp != NULL){
			if( left_mv_tmp->bi_mode <= 6 || left_mv_tmp->bi_mode == 8 || (left_mv_tmp->bi_mode == 7 && left_mv_tmp->aff_mrg == NO) ){
				mrg_left[i]->mvx = left_mv_tmp->mvx;	mrg_left[i]->mvy = left_mv_tmp->mvy;	mrg_left[i]->lifting_mode = CONNECTED;
			}else{
				mrg_left[i]->mvx = ult11;	mrg_left[i]->mvy = ult12;	mrg_left[i]->lifting_mode = CONNECTED;
			}
		  }
		  if(right_mv_tmp != NULL){
			if( right_mv_tmp->bi_mode <= 6 || right_mv_tmp->bi_mode == 8 || (right_mv_tmp->bi_mode == 7 && right_mv_tmp->aff_mrg == NO) ){
				mrg_right[i]->mvx = right_mv_tmp->mvx;	mrg_right[i]->mvy = right_mv_tmp->mvy;	mrg_right[i]->lifting_mode = CONNECTED;
			}else{
				mrg_right[i]->mvx = ult21;	mrg_right[i]->mvy = ult22;	mrg_right[i]->lifting_mode = CONNECTED;
			}
		  }
//Additional TMP merge candidate
		  if( (mrg_left[i]->mvx == (float)HUGE_VAL && mrg_left[i]->mvy == (float)HUGE_VAL) &&
			(mrg_right[i]->mvx == (float)HUGE_VAL && mrg_right[i]->mvy == (float)HUGE_VAL) ){
			if(fmv3_array != NULL)
			  left_mv_tmp = find_block(x_pos-xblk,y_pos-yblk,fmv3_array,info,t_level,&ult11,&ult12,&get_xblk,&get_yblk,0,0);

		      if(fmv2_array != NULL && fmv4_array != NULL)
			    right_mv_tmp = find_block(x_pos-xblk,y_pos-yblk,fmv4_array,info,t_level,&ult21,&ult22,&get_xblk,&get_yblk,0,0);

		      if(left_mv_tmp != NULL){
				if( left_mv_tmp->bi_mode <= 6 || left_mv_tmp->bi_mode == 8 || (left_mv_tmp->bi_mode == 7 && left_mv_tmp->aff_mrg == NO) ){
					mrg_left[i]->mvx = left_mv_tmp->mvx;	mrg_left[i]->mvy = left_mv_tmp->mvy;	mrg_left[i]->lifting_mode = CONNECTED;
				}else{
					mrg_left[i]->mvx = ult11;	mrg_left[i]->mvy = ult12;	mrg_left[i]->lifting_mode = CONNECTED;
				}
			  }
			  if(right_mv_tmp != NULL){
				if( right_mv_tmp->bi_mode <= 6 || right_mv_tmp->bi_mode == 8 || (right_mv_tmp->bi_mode == 7 && right_mv_tmp->aff_mrg == NO) ){
					mrg_right[i]->mvx = right_mv_tmp->mvx;	mrg_right[i]->mvy = right_mv_tmp->mvy;	mrg_right[i]->lifting_mode = CONNECTED;
				}else{
					mrg_right[i]->mvx = ult21;	mrg_right[i]->mvy = ult22;	mrg_right[i]->lifting_mode = CONNECTED;
				}
			  }
		  }
//Additional TMP merge candidate
		  break;
		}
    }
	////////////////////////
	for(i = 0;i <= 3;i ++){
	  if( (x_pos-(mrg_left[i]->mvx)<0||x_pos-(mrg_left[i]->mvx)+xblk>hor||y_pos-(mrg_left[i]->mvy)<0|| 
	       y_pos-(mrg_left[i]->mvy)+yblk>ver) && (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL ) ){
		mrg_left[i]->mvx = (float)HUGE_VAL;
        mrg_left[i]->mvy = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

		mrg_right[i]->mvx = (float)HUGE_VAL;
        mrg_right[i]->mvy = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
	  }
	  if( (x_pos-(mrg_right[i]->mvx)<0||x_pos-(mrg_right[i]->mvx)+xblk>hor||y_pos-(mrg_right[i]->mvy)<0||
		  y_pos-(mrg_right[i]->mvy)+yblk>ver)  && (mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL ) ){
		mrg_left[i]->mvx = (float)HUGE_VAL;
        mrg_left[i]->mvy = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

		mrg_right[i]->mvx = (float)HUGE_VAL;
        mrg_right[i]->mvy = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
	  }
	}

	//check for repeated candidates
	for(i = 3;i >= 1;i --){
		for(k = i-1;k >= 0;k --){
			if( (mrg_left[i]->mvx == mrg_left[k]->mvx && mrg_left[i]->mvy == mrg_left[k]->mvy ) &&
			(mrg_right[i]->mvx == mrg_right[k]->mvx && mrg_right[i]->mvy == mrg_right[k]->mvy ) ){
				mrg_left[i]->mvx = (float)HUGE_VAL;
				mrg_left[i]->mvy = (float)HUGE_VAL;
				mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
				mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
				mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

				mrg_right[i]->mvx = (float)HUGE_VAL;
				mrg_right[i]->mvy = (float)HUGE_VAL;
				mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
				mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
				mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
			}
		}
	}
	//check for HUGE_VALs
    k = 0;
    for(i = 0;i <= 3;i++){
	    if( (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL) ||
			(mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL) ){
		  med_l_px[k] = mrg_left[i]->mvx;
		  med_l_py[k] = mrg_left[i]->mvy;
		  med_l_aff1x[k] = mrg_left[i]->aff1_mvx;
		  med_l_aff1y[k] = mrg_left[i]->aff1_mvy;
		  med_l_aff2x[k] = mrg_left[i]->aff2_mvx;
		  med_l_aff2y[k] = mrg_left[i]->aff2_mvy;
		  med_l_aff3x[k] = mrg_left[i]->aff3_mvx;
		  med_l_aff3y[k] = mrg_left[i]->aff3_mvy;

		  med_r_px[k] = mrg_right[i]->mvx;
		  med_r_py[k] = mrg_right[i]->mvy;
		  med_r_aff1x[k] = mrg_right[i]->aff1_mvx;
		  med_r_aff1y[k] = mrg_right[i]->aff1_mvy;
		  med_r_aff2x[k] = mrg_right[i]->aff2_mvx;
		  med_r_aff2y[k] = mrg_right[i]->aff2_mvy;
		  med_r_aff3x[k] = mrg_right[i]->aff3_mvx;
		  med_r_aff3y[k] = mrg_right[i]->aff3_mvy;
		  k++;
	    }
    }
	assert(k <= 4);
	for(i = 0;i <= 3;i++){
		mrg_left[i]->mvx = med_l_px[i];
		mrg_left[i]->mvy = med_l_py[i];
		mrg_left[i]->aff1_mvx = med_l_aff1x[i];
		mrg_left[i]->aff1_mvy = med_l_aff1y[i];
		mrg_left[i]->aff2_mvx = med_l_aff2x[i];
		mrg_left[i]->aff2_mvy = med_l_aff2y[i];
		mrg_left[i]->aff3_mvx = med_l_aff3x[i];
		mrg_left[i]->aff3_mvy = med_l_aff3y[i];

		mrg_right[i]->mvx  = med_r_px[i];
		mrg_right[i]->mvy  = med_r_py[i];
		mrg_right[i]->aff1_mvx = med_r_aff1x[i];
		mrg_right[i]->aff1_mvy = med_r_aff1y[i];
		mrg_right[i]->aff2_mvx = med_r_aff2x[i];
		mrg_right[i]->aff2_mvy = med_r_aff2y[i];
		mrg_right[i]->aff3_mvx = med_r_aff3x[i];
		mrg_right[i]->aff3_mvy = med_r_aff3y[i];

	}
	for(i=0;i<=3;i++){
		med_l_px[i] = (float)HUGE_VAL;
		med_l_py[i] = (float)HUGE_VAL;
		med_l_aff1x[i] = (float)HUGE_VAL;
		med_l_aff1y[i] = (float)HUGE_VAL;
		med_l_aff2x[i] = (float)HUGE_VAL;
		med_l_aff2y[i] = (float)HUGE_VAL;
		med_l_aff3x[i] = (float)HUGE_VAL;
		med_l_aff3y[i] = (float)HUGE_VAL;

		med_r_px[i] = (float)HUGE_VAL;
		med_r_py[i] = (float)HUGE_VAL;
		med_r_aff1x[i] = (float)HUGE_VAL;
		med_r_aff1y[i] = (float)HUGE_VAL;
		med_r_aff2x[i] = (float)HUGE_VAL;
		med_r_aff2y[i] = (float)HUGE_VAL;
		med_r_aff3x[i] = (float)HUGE_VAL;
		med_r_aff3y[i] = (float)HUGE_VAL;
	}
//E
	aff_blk = 0;

	for(i = 0; i <= 3; i++){
		if( (mrg_left[i]->mvx == (float)HUGE_VAL && mrg_left[i]->mvy == (float)HUGE_VAL) &&
			(mrg_right[i]->mvx == (float)HUGE_VAL && mrg_right[i]->mvy == (float)HUGE_VAL) ){

			left_mve = find_block(x_pos-1,y_pos+yblk,fmv1_array,info,t_level,&ule11,&ule12,&get_xblk,&get_yblk,0,0);

			find_block(x_pos-3,y_pos+yblk+1,fmv1_array,info,t_level,&left_affx1,&left_affy1,&get_xblk,&get_yblk,0,0);
			find_block(x_pos-2,y_pos+yblk+1,fmv1_array,info,t_level,&left_affx2,&left_affy2,&get_xblk,&get_yblk,0,0);
			find_block(x_pos-3,y_pos+yblk+2,fmv1_array,info,t_level,&left_affx3,&left_affy3,&get_xblk,&get_yblk,0,0);

			if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){

				assert(left_mve->bi_mode>=9  || (left_mve->bi_mode == 7 && left_mve->aff_mrg == YES) );

				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_aff_mrg_mv(mrg_left[i], x_pos, y_pos, x_pos-3,y_pos+yblk+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				aff_blk = 1;
			}

			if(fmv2_array != NULL){
				right_mve = find_block(x_pos-1,y_pos+yblk,fmv2_array,info,t_level,&ule21,&ule22,&get_xblk,&get_yblk,0,0);

				find_block(x_pos-3,y_pos+yblk+1,fmv2_array,info,t_level,&right_affx1,&right_affy1,&get_xblk,&get_yblk,0,0);
				find_block(x_pos-2,y_pos+yblk+1,fmv2_array,info,t_level,&right_affx2,&right_affy2,&get_xblk,&get_yblk,0,0);
				find_block(x_pos-3,y_pos+yblk+2,fmv2_array,info,t_level,&right_affx3,&right_affy3,&get_xblk,&get_yblk,0,0);

				if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){

					assert(right_mve->bi_mode>=9 || (right_mve->bi_mode == 7 && right_mve->aff_mrg == YES) );

					aff_dmvx1 = (right_affx2 - right_affx1);
					aff_dmvy1 = (right_affy2 - right_affy1);

					aff_dmvx2 = (right_affx3 - right_affx1);
					aff_dmvy2 = (right_affy3 - right_affy1);

					get_aff_mrg_mv(mrg_right[i], x_pos, y_pos, x_pos-3,y_pos+yblk+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
					aff_blk = 1;
				}
				
				if(aff_blk == 1){
					if(left_mve != NULL && left_mve->is_predictor == YES){
						aff_dmvx1 = (left_affx2 - left_affx1);
						aff_dmvy1 = (left_affy2 - left_affy1);

						aff_dmvx2 = (left_affx3 - left_affx1);
						aff_dmvy2 = (left_affy3 - left_affy1);

						get_aff_mrg_mv(mrg_left[i], x_pos, y_pos, x_pos-3,y_pos+yblk+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
					}
					if(right_mve != NULL && right_mve->is_predictor == YES){
						aff_dmvx1 = (right_affx2 - right_affx1);
						aff_dmvy1 = (right_affy2 - right_affy1);

						aff_dmvx2 = (right_affx3 - right_affx1);
						aff_dmvy2 = (right_affy3 - right_affy1);

						get_aff_mrg_mv(mrg_right[i], x_pos, y_pos, x_pos-3,y_pos+yblk+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
					}
				}
			
			}

			if(left_mve != NULL){
				if( left_mve->bi_mode <= 6 || left_mve->bi_mode == 8 || (left_mve->bi_mode == 7 && left_mve->aff_mrg == NO) ){
					mrg_left[i]->mvx = left_mve->mvx;	mrg_left[i]->mvy = left_mve->mvy;	mrg_left[i]->lifting_mode = CONNECTED;
				}else{
					mrg_left[i]->mvx = ule11;	mrg_left[i]->mvy = ule12;	mrg_left[i]->lifting_mode = CONNECTED;
				}
			}
			if(right_mve != NULL){
				if( right_mve->bi_mode <= 6 || right_mve->bi_mode == 8 || (right_mve->bi_mode == 7 && right_mve->aff_mrg == NO) ){
					mrg_right[i]->mvx = right_mve->mvx;	mrg_right[i]->mvy = right_mve->mvy;	mrg_right[i]->lifting_mode = CONNECTED;
				}else{
					mrg_right[i]->mvx = ule21;	mrg_right[i]->mvy = ule22;	mrg_right[i]->lifting_mode = CONNECTED;
				}
			}
			break;
		}
	}

	aff_blk = 0;
	////////////////////////
	for(i = 0;i <= 3;i ++){
	  if( (x_pos-(mrg_left[i]->mvx)<0||x_pos-(mrg_left[i]->mvx)+xblk>hor||y_pos-(mrg_left[i]->mvy)<0|| 
	       y_pos-(mrg_left[i]->mvy)+yblk>ver) && (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL ) ){
		mrg_left[i]->mvx = (float)HUGE_VAL;
        mrg_left[i]->mvy = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

		mrg_right[i]->mvx = (float)HUGE_VAL;
        mrg_right[i]->mvy = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
	  }
	  if( (x_pos-(mrg_right[i]->mvx)<0||x_pos-(mrg_right[i]->mvx)+xblk>hor||y_pos-(mrg_right[i]->mvy)<0||
		  y_pos-(mrg_right[i]->mvy)+yblk>ver)  && (mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL ) ){
		mrg_left[i]->mvx = (float)HUGE_VAL;
        mrg_left[i]->mvy = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

		mrg_right[i]->mvx = (float)HUGE_VAL;
        mrg_right[i]->mvy = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
		mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
	  }
	}

	//check for repeated candidates
	for(i = 3;i >= 1;i --){
		for(k = i-1;k >= 0;k --){
			if( (mrg_left[i]->mvx == mrg_left[k]->mvx && mrg_left[i]->mvy == mrg_left[k]->mvy ) &&
			(mrg_right[i]->mvx == mrg_right[k]->mvx && mrg_right[i]->mvy == mrg_right[k]->mvy ) ){
				mrg_left[i]->mvx = (float)HUGE_VAL;
				mrg_left[i]->mvy = (float)HUGE_VAL;
				mrg_left[i]->aff1_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff1_mvy = (float)HUGE_VAL;
				mrg_left[i]->aff2_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff2_mvy = (float)HUGE_VAL;
				mrg_left[i]->aff3_mvx = (float)HUGE_VAL;
				mrg_left[i]->aff3_mvy = (float)HUGE_VAL;

				mrg_right[i]->mvx = (float)HUGE_VAL;
				mrg_right[i]->mvy = (float)HUGE_VAL;
				mrg_right[i]->aff1_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff1_mvy = (float)HUGE_VAL;
				mrg_right[i]->aff2_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff2_mvy = (float)HUGE_VAL;
				mrg_right[i]->aff3_mvx = (float)HUGE_VAL;
				mrg_right[i]->aff3_mvy = (float)HUGE_VAL;
			}
		}
	}
	//check for HUGE_VALs
    k = 0;
    for(i = 0;i <= 3;i++){
	    if( (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL) ||
			(mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL) ){
		  med_l_px[k] = mrg_left[i]->mvx;
		  med_l_py[k] = mrg_left[i]->mvy;
		  med_l_aff1x[k] = mrg_left[i]->aff1_mvx;
		  med_l_aff1y[k] = mrg_left[i]->aff1_mvy;
		  med_l_aff2x[k] = mrg_left[i]->aff2_mvx;
		  med_l_aff2y[k] = mrg_left[i]->aff2_mvy;
		  med_l_aff3x[k] = mrg_left[i]->aff3_mvx;
		  med_l_aff3y[k] = mrg_left[i]->aff3_mvy;

		  med_r_px[k] = mrg_right[i]->mvx;
		  med_r_py[k] = mrg_right[i]->mvy;
		  med_r_aff1x[k] = mrg_right[i]->aff1_mvx;
		  med_r_aff1y[k] = mrg_right[i]->aff1_mvy;
		  med_r_aff2x[k] = mrg_right[i]->aff2_mvx;
		  med_r_aff2y[k] = mrg_right[i]->aff2_mvy;
		  med_r_aff3x[k] = mrg_right[i]->aff3_mvx;
		  med_r_aff3y[k] = mrg_right[i]->aff3_mvy;
		  k++;
	    }
    }
	assert(k <= 4);
	for(i = 0;i <= 3;i++){
		mrg_left[i]->mvx = med_l_px[i];
		mrg_left[i]->mvy = med_l_py[i];
		mrg_left[i]->aff1_mvx = med_l_aff1x[i];
		mrg_left[i]->aff1_mvy = med_l_aff1y[i];
		mrg_left[i]->aff2_mvx = med_l_aff2x[i];
		mrg_left[i]->aff2_mvy = med_l_aff2y[i];
		mrg_left[i]->aff3_mvx = med_l_aff3x[i];
		mrg_left[i]->aff3_mvy = med_l_aff3y[i];

		mrg_right[i]->mvx  = med_r_px[i];
		mrg_right[i]->mvy  = med_r_py[i];
		mrg_right[i]->aff1_mvx = med_r_aff1x[i];
		mrg_right[i]->aff1_mvy = med_r_aff1y[i];
		mrg_right[i]->aff2_mvx = med_r_aff2x[i];
		mrg_right[i]->aff2_mvy = med_r_aff2y[i];
		mrg_right[i]->aff3_mvx = med_r_aff3x[i];
		mrg_right[i]->aff3_mvy = med_r_aff3y[i];

	}
}
////////////////////////////////////////////////////

//Added on 01.17.2016
float med(float a, float b, float c)
{
  return ((a>b) ? ((a>c) ? ((b>c) ? b : c) : a) : 
          ((b>c) ? ((a>c) ? a : c) : b));
}

vector_ptr find_block(int x, int y, vector_ptr root, videoinfo info, int t_level, float *mvx, float *mvy,
					  int *get_xblk, int *get_yblk, int get_addx, int get_addy){ //Find the 
	int i,j,k,sum;
	int pos;
	int x0,y0;
	vector_ptr v1;
	int xblk = info.xblk[t_level],yblk = info.yblk[t_level];
	int xblk2, yblk2;
	
	int addx,addy;
	int int_mvx,int_mvy;
	float fmvx1,fmvy1,fmvx2,fmvy2,fmvx3,fmvy3,fmvx4,fmvy4;

	sum = 0;
	pos = 0;

	if(x < 0 || x >= info.ywidth || y < 0 || y >= info.yheight){
		v1 = NULL;
		*mvx = (float)HUGE_VAL;
		*mvy = (float)HUGE_VAL;
		return v1;
	}

	while(y >= sum){
		sum += info.yblk[t_level];
		pos += info.xnum[t_level];
	}
	pos -= info.xnum[t_level];

	sum = 0;
	while(x >= sum){
		sum += info.xblk[t_level];
		pos ++;
	}
	pos --;

	v1 = &root[pos];
	x0 = info.xblk[t_level] * ( pos % info.xnum[t_level]);
	y0 = info.yblk[t_level] * (int)(pos / info.xnum[t_level]);

//	printf("Original position: x0 = %d, y0 = %d\n",x0,y0);

	while(v1->child){
		if( (y0 + yblk/2) > y){
			if( (x0 + xblk/2) > x){
				v1 = v1->child0;
				xblk/=2;
				yblk/=2;
			}
			else{
				v1 = v1->child1;
				x0 = x0 + xblk/2;
				xblk/=2;
				yblk/=2;
			}
		}
		else{
			if( (x0 + xblk/2) > x){
				v1 = v1->child2;
				y0 = y0 + yblk/2;
				xblk/=2;
				yblk/=2;
			}
			else{
				v1 = v1->child3;
				x0 = x0 + xblk/2;
				y0 = y0 + yblk/2;
				xblk/=2;
				yblk/=2;
			}
		}	
	}

	//////////	Added on 04.03.2016		//////////
	xblk2 = ( x0 + xblk <= info.ywidth ) ? xblk : info.ywidth - x0;
	yblk2 = ( y0 + yblk <= info.yheight ) ? yblk : info.yheight - y0;
	//////////////////////////////////////////////

	if(get_addx == 0 && get_addy == 0){
		if( (x-x0+1)/xblk2 == 1 )
			addx = 1;
		else
			addx = 0;

		if( (y-y0+1)/yblk2 == 1 )
			addy = 1;
		else
			addy = 0;
	}else{
		addx = get_addx;
		addy = get_addy;
	}

	//////////	Modified on 01.05.2016	//////////
	if( (v1->bi_mode >= 0 && v1->bi_mode <= 6) || v1->bi_mode == 8 || (v1->bi_mode == 7 && v1->aff_mrg == NO) ){
		*mvx = v1->mvx;
		*mvy = v1->mvy;
	}
	else if( (v1->bi_mode >= 9 && v1->bi_mode <= 11) || (v1->bi_mode == 7 && v1->aff_mrg == YES) ){
//		*mvx = (v1->aff2_mvx - v1->aff1_mvx)*(x-x0 + addx)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-y0 + addy)/yblk2 + v1->aff1_mvx;
//		*mvy = (v1->aff2_mvy - v1->aff1_mvy)*(x-x0 + addx)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-y0 + addy)/yblk2 + v1->aff1_mvy;
//Added on 01.24.2019
		if(get_addx == 0 && get_addy == 0){
			*mvx = (v1->aff2_mvx - v1->aff1_mvx)*(x-x0 + addx)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-y0 + addy)/yblk2 + v1->aff1_mvx;
			*mvy = (v1->aff2_mvy - v1->aff1_mvy)*(x-x0 + addx)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-y0 + addy)/yblk2 + v1->aff1_mvy;
		}else if(get_addx == 0 && get_addy == 1){
			fmvx3 = (v1->aff2_mvx - v1->aff1_mvx)*(x-x0)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-1-y0)/yblk2 + v1->aff1_mvx;
			fmvy3 = (v1->aff2_mvy - v1->aff1_mvy)*(x-x0)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-1-y0)/yblk2 + v1->aff1_mvy;

			fmvx4 = (v1->aff2_mvx - v1->aff1_mvx)*(x-x0)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-2-y0)/yblk2 + v1->aff1_mvx;
			fmvy4 = (v1->aff2_mvy - v1->aff1_mvy)*(x-x0)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-2-y0)/yblk2 + v1->aff1_mvy;

			*mvx = fmvx3 + 2 * (fmvx3 - fmvx4);
			*mvy = fmvy3 + 2 * (fmvy3 - fmvy4);
		}else if(get_addx == 1 && get_addy == 0){
			fmvx1 = (v1->aff2_mvx - v1->aff1_mvx)*(x-1-x0)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-y0)/yblk2 + v1->aff1_mvx;
			fmvy1 = (v1->aff2_mvy - v1->aff1_mvy)*(x-1-x0)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-y0)/yblk2 + v1->aff1_mvy;

			fmvx2 = (v1->aff2_mvx - v1->aff1_mvx)*(x-2-x0)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-y0)/yblk2 + v1->aff1_mvx;
			fmvy2 = (v1->aff2_mvy - v1->aff1_mvy)*(x-2-x0)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-y0)/yblk2 + v1->aff1_mvy;

			*mvx = fmvx1 + 2 * (fmvx1 - fmvx2);
			*mvy = fmvy1 + 2 * (fmvy1 - fmvy2);
		}else if(get_addx == 1 && get_addy == 1){
			fmvx1 = (v1->aff2_mvx - v1->aff1_mvx)*(x-1-x0)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-1-y0)/yblk2 + v1->aff1_mvx;
			fmvy1 = (v1->aff2_mvy - v1->aff1_mvy)*(x-1-x0)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-1-y0)/yblk2 + v1->aff1_mvy;

			fmvx2 = (v1->aff2_mvx - v1->aff1_mvx)*(x-2-x0)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-1-y0)/yblk2 + v1->aff1_mvx;
			fmvy2 = (v1->aff2_mvy - v1->aff1_mvy)*(x-2-x0)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-1-y0)/yblk2 + v1->aff1_mvy;

			fmvx3 = (v1->aff2_mvx - v1->aff1_mvx)*(x-1-x0)/xblk2 + (v1->aff3_mvx - v1->aff1_mvx)*(y-2-y0)/yblk2 + v1->aff1_mvx;
			fmvy3 = (v1->aff2_mvy - v1->aff1_mvy)*(x-1-x0)/xblk2 + (v1->aff3_mvy - v1->aff1_mvy)*(y-2-y0)/yblk2 + v1->aff1_mvy;

			*mvx = fmvx1 + 2 * (fmvx1 - fmvx3) + 2 * (fmvx1 - fmvx2);
			*mvy = fmvy1 + 2 * (fmvy1 - fmvy3) + 2 * (fmvy1 - fmvy2);
		}else
			assert(0);
///////////////////////////////////////////////////////////////
	}

	if(get_addx == 1 || get_addy == 1){
		(*mvx) = (*mvx) * AFF_SUBPEL;
		(*mvy) = (*mvy) * AFF_SUBPEL;
		int_mvx = (int)(*mvx);
		int_mvy = (int)(*mvy);
		(*mvx) = (float)(int_mvx);
		(*mvy) = (float)(int_mvy);
		(*mvx) = (*mvx) / AFF_SUBPEL;
		(*mvy) = (*mvy) / AFF_SUBPEL;
	}
//////////////////////////////////////////////

	if(v1->bi_mode == UNDEFINED || v1->bi_mode == DIRECTIONAL_IBLOCK || v1->is_predictor == NO){
		v1 = NULL;
		*mvx = (float)HUGE_VAL;
		*mvy = (float)HUGE_VAL;
	}

//	printf("mvx = %f, mvy = %f\n",*mvx,*mvy);

	*get_xblk = xblk2;
	*get_yblk = yblk2;

	return v1;
}

float find_affine_SAD(float *frame_cur, float *frame_ref1, float *frame_ref2, float *upframe1, float *upframe2, ModeInfo v11, ModeInfo v12, ModeInfo v13,
					  ModeInfo v21, ModeInfo v22, ModeInfo v23, int xblk, int yblk, int cx, int cy, int hor, int ver, 
					  int best_mode, int subpel, int type){
//type 0: interpolation in both left and right
//type 1: interpolation in left only
//type 2: interpolation in right only

	float mse,fx1,fy1,fx2,fy2;
	float px1,py1,px2,py2;
	float ptemp,diff,sum;
	float mvx1,mvy1,mvx2,mvy2;
	float ptemp1,ptemp2;

	float dhorx1,dverx1,dhory1,dvery1,dhorx2,dverx2,dhory2,dvery2;	//Added by Yuan Liu on 04.23.2016
	float mvx1_int, mvx2_int, mvy1_int, mvy2_int;					//Added by Yuan Liu on 04.23.2016

	int x,y,m,m0;
	int uphor,upver;
	int upx11,upx12,upy11,upy12;//
	int upx21,upx22,upy21,upy22;
	int pos1,pos2;

	int scale = 1 << subpel;
	scale = scale << ADD_SUB; //Added on 03.13.2017

	float accu = 1/((float)scale);

	int x1,x2,y1,y2;
	float dX,dY;

//	printf("Enter affine SAD!\n");

	uphor = ( hor - 1 ) * scale + 1;
	upver = ( ver - 1 ) * scale + 1;

	dhorx1 = (v12.mvx - v11.mvx)/xblk;
	dhory1 = (v12.mvy - v11.mvy)/xblk;

	dverx1 = (v13.mvx - v11.mvx)/yblk;
	dvery1 = (v13.mvy - v11.mvy)/yblk;

	dhorx2 = (v22.mvx - v21.mvx)/xblk;
	dhory2 = (v22.mvy - v21.mvy)/xblk;

	dverx2 = (v23.mvx - v21.mvx)/yblk;
	dvery2 = (v23.mvy - v21.mvy)/yblk;

	m = cy * hor + cx;
	m0 = 0;
	sum = 0;

	mvx1 = v11.mvx;
	mvy1 = v11.mvy;

	mvx1_int = mvx1;
	mvy1_int = mvy1;

	mvx2 = v21.mvx;
	mvy2 = v21.mvy;

	mvx2_int = mvx2;
	mvy2_int = mvy2;

	for( y = 0; y < yblk; y++ ) {
		for( x = 0; x < xblk; x++ ) {

//			mvx1 = (v12.mvx - v11.mvx)*x/xblk + (v13.mvx - v11.mvx)*y/yblk + v11.mvx;
//			mvy1 = (v12.mvy - v11.mvy)*x/xblk + (v13.mvy - v11.mvy)*y/yblk + v11.mvy;
			
			px1 = (float)cx - mvx1;
			py1 = (float)cy - mvy1;
			
			fx1 = px1 + (float)x;
			fy1 = py1 + (float)y;
			
//			if( ( ( int )( fx1 * 16 ) ) % 16 == 0 )fx1 = floor(fx1);
//			if( ( ( int )( fy1 * 16 ) ) % 16 == 0 )fy1 = floor(fy1);
			
			if(best_mode == RIGHT_CONNECTED_AFF){//If RIGHT_CONNECTED

				x2 = ( int )ceil( fx1 );
				x1 = ( int )floor( fx1 );
				y2 = ( int )ceil( fy1 );
				y1 = ( int )floor( fy1 );

				if( x2 >= (hor-1) || x1 < 0 || y2 >= (ver-1) || y1 < 0 ){
//					assert(fx1 > (hor-1) || fx1 < 0 || fy1 > (ver-1) || fy1 < 0);
					ptemp = (float)HUGE_VAL;
				}
				else{
					upx11 = (int)(fx1/accu);
					upy11 = (int)(fy1/accu);

					pos1 = upy11*uphor + upx11;

					dX = fx1 * scale - upx11;
					dY = fy1 * scale - upy11;

					if(x1 == x2 && y1 == y2)
						ptemp = upframe2[pos1];
					else if(x1 == x2)
						ptemp = (1 - dY) * upframe2[pos1] + dY * upframe2[pos1 + uphor];
					else if(y1 == y2)
						ptemp = (1 - dX) * upframe2[pos1] + dX * upframe2[pos1 + 1];
					else
						ptemp = ( 1 - dY ) * ( ( 1 - dX ) * upframe2[pos1] + dX * upframe2[pos1 + 1] ) + dY * ( ( 1 - dX ) * upframe2[pos1 + uphor] + dX * upframe2[pos1 + uphor + 1] );
				}

				mvx1 += dhorx1;
				mvy1 += dhory1;

			}else if(best_mode == LEFT_CONNECTED_AFF){//If LEFT_CONNECTED
				x2 = ( int )ceil( fx1 );
				x1 = ( int )floor( fx1 );
				y2 = ( int )ceil( fy1 );
				y1 = ( int )floor( fy1 );

				if( x2 >= (hor-1) || x1 < 0 || y2 >= (ver-1) || y1 < 0 ){
					ptemp = (float)HUGE_VAL;
				}
				else{
					upx11 = (int)(fx1/accu);
					upy11 = (int)(fy1/accu);

					pos1 = upy11*uphor + upx11;

					dX = fx1 * scale - upx11;
					dY = fy1 * scale - upy11;

					if(x1 == x2 && y1 == y2)
						ptemp = upframe1[pos1];
					else if(x1 == x2)
						ptemp = (1 - dY) * upframe1[pos1] + dY * upframe1[pos1 + uphor];
					else if(y1 == y2)
						ptemp = (1 - dX) * upframe1[pos1] + dX * upframe1[pos1 + 1];
					else
						ptemp = ( 1 - dY ) * ( ( 1 - dX ) * upframe1[pos1] + dX * upframe1[pos1 + 1] ) + dY * ( ( 1 - dX ) * upframe1[pos1 + uphor] + dX * upframe1[pos1 + uphor + 1] );
				}

				mvx1 += dhorx1;
				mvy1 += dhory1;
			}
			else if(best_mode == BI_CONNECTED_AFF){//If BI_CONNECTED
//LEFT MV
			  if(type==0 || type==1 || type==3 || type == 4){
				x2 = ( int )ceil( fx1 );
				x1 = ( int )floor( fx1 );
				y2 = ( int )ceil( fy1 );
				y1 = ( int )floor( fy1 );

				if( x2 >= (hor-1) || x1 < 0 || y2 >= (ver-1) || y1 < 0 ){
//					assert(fx1 > (hor-1) || fx1 < 0 || fy1 > (ver-1) || fy1 < 0);
					ptemp1 = (float)HUGE_VAL;
				}
				else{
					upx11 = (int)(fx1/accu);
					upy11 = (int)(fy1/accu);

					pos1 = upy11*uphor + upx11;

					dX = fx1 * scale - upx11;
					dY = fy1 * scale - upy11;

					if(x1 == x2 && y1 == y2)
						ptemp1 = upframe1[pos1];
					else if(x1 == x2)
						ptemp1 = (1 - dY) * upframe1[pos1] + dY * upframe1[pos1 + uphor];
					else if(y1 == y2)
						ptemp1 = (1 - dX) * upframe1[pos1] + dX * upframe1[pos1 + 1];
					else
						ptemp1 = ( 1 - dY ) * ( ( 1 - dX ) * upframe1[pos1] + dX * upframe1[pos1 + 1] ) + dY * ( ( 1 - dX ) * upframe1[pos1 + uphor] + dX * upframe1[pos1 + uphor + 1] );

					if(type==0 || type==1)
						block_buff1[m0] = ptemp1;

					mvx1 += dhorx1;
					mvy1 += dhory1;
				}
			  }//if type = 0 or 1
			  else
				  ptemp1 = block_buff1[m0];

//RIGHT MV
			  if(type==0 || type==2 || type==3 || type == 4){
//				mvx2 = (v22.mvx - v21.mvx)*x/xblk + (v23.mvx - v21.mvx)*y/yblk + v21.mvx;
//				mvy2 = (v22.mvy - v21.mvy)*x/xblk + (v23.mvy - v21.mvy)*y/yblk + v21.mvy;

				px2 = (float)cx - mvx2;
				py2 = (float)cy - mvy2;

				fx2 = px2 + (float)x;
				fy2 = py2 + (float)y;

				x2 = ( int )ceil( fx2 );
				x1 = ( int )floor( fx2 );
				y2 = ( int )ceil( fy2 );
				y1 = ( int )floor( fy2 );

				if( x2 >= (hor-1) || x1 < 0 || y2 >= (ver-1) || y1 < 0 ){
//					assert(fx2 > (hor-1) || fx2 < 0 || fy2 > (ver-1) || fy2 < 0);
					ptemp2 = (float)HUGE_VAL;
				}
				else{
					upx21 = (int)(fx2/accu);
					upy21 = (int)(fy2/accu);

					pos2 = upy21*uphor + upx21;

					dX = fx2 * scale - upx21;
					dY = fy2 * scale - upy21;

					if(x1 == x2 && y1 == y2)
						ptemp2 = upframe2[pos2];
					else if(x1 == x2)
						ptemp2 = (1 - dY) * upframe2[pos2] + dY * upframe2[pos2 + uphor];
					else if(y1 == y2)
						ptemp2 = (1 - dX) * upframe2[pos2] + dX * upframe2[pos2 + 1];
					else
						ptemp2 = ( 1 - dY ) * ( ( 1 - dX ) * upframe2[pos2] + dX * upframe2[pos2 + 1] ) + dY * ( ( 1 - dX ) * upframe2[pos2 + uphor] + dX * upframe2[pos2 + uphor + 1] );
					
					if(type==0 || type==2)
						block_buff2[m0] = ptemp2;
				}
				mvx2 += dhorx2;
				mvy2 += dhory2;
			  }//type = 0 or 2
			  else
				  ptemp2 = block_buff2[m0];

				if(ptemp1 == (float)HUGE_VAL || ptemp2 == (float)HUGE_VAL)
					ptemp = (float)HUGE_VAL;
				else
					ptemp = 0.5f * ( ptemp1 + ptemp2 );

				m0++;
			}//if BI_AFF
			else
				assert(0);

			if(ptemp == (float)HUGE_VAL){
//				printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);
				return ptemp;
				break;
			}

			diff = frame_cur[m] - ptemp;
			sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
			m++; 

			if(type == 4)
				printf("%d\t",(int)fabs(diff) );

		}//  x	
		if(ptemp == (float)HUGE_VAL)break;
		m += hor - xblk;

		mvx1 = mvx1_int + dverx1;
		mvy1 = mvy1_int + dvery1;

		mvx1_int = mvx1;
		mvy1_int = mvy1;

		if(best_mode == BI_CONNECTED_AFF){
			mvx2 = mvx2_int + dverx2;
			mvy2 = mvy2_int + dvery2;

			mvx2_int = mvx2;
			mvy2_int = mvy2;
		}

		if(type == 4)
			printf("\n");

	}//  y                   

//	printf("affine SAD sum: %f\n",sum);
//	sum /= xblk * yblk;

//	printf("Exit affine SAD!\n");

/*	if(type == 4)
		printf("\n\n");
*/	
	return sum;
}

float find_affine_MSE(float *frame_cur, float *frame_ref1, float *frame_ref2, float *aff_ref_var, ModeInfo v11, ModeInfo v12, ModeInfo v13, ModeInfo v21, ModeInfo v22, ModeInfo v23, int xblk, int yblk, int cx, int cy, int hor, int ver, int best_mode){
	float mse,fx1,fy1,fx2,fy2;
	float px1,py1,px2,py2;
	float ptemp,diff,sum;
	float mvx1,mvy1,mvx2,mvy2;

	int x,y,m;

	int ref_x1,ref_y1,ref_x2,ref_y2;

	float aff_ref_mean1, aff_ref_mean2;

	m = cy * hor + cx;
	sum = 0;

	*aff_ref_var = 0;
	aff_ref_mean1 = 0;
	aff_ref_mean2 = 0;

	for( y = 0; y < yblk; y++ ) {
		for( x = 0; x < xblk; x++ ) {

			mvx1 = (v12.mvx - v11.mvx)*x/xblk + (v13.mvx - v11.mvx)*y/yblk + v11.mvx;
			mvy1 = (v12.mvy - v11.mvy)*x/xblk + (v13.mvy - v11.mvy)*y/yblk + v11.mvy;
				
			px1 = cx - mvx1;
			py1 = cy - mvy1;
				
			fx1 = px1 + x;
			fy1 = py1 + y;

			position( &ref_x1, &ref_y1, fx1, fy1, mvx1, mvy1, hor, ver );
			
			if(best_mode == RIGHT_CONNECTED_AFF){
				ptemp = interpolate( fx1, fy1, frame_ref2, hor, ver, TYPE );
				*aff_ref_var += frame_ref2[ref_y1 * hor + ref_x1] * frame_ref2[ref_y1 * hor + ref_x1];
				aff_ref_mean1 += frame_ref2[ref_y1 * hor + ref_x1];

			}else if(best_mode == LEFT_CONNECTED_AFF){
				ptemp = interpolate( fx1, fy1, frame_ref1, hor, ver, TYPE );
				*aff_ref_var += frame_ref1[ref_y1 * hor + ref_x1] * frame_ref1[ref_y1 * hor + ref_x1];
				aff_ref_mean1 += frame_ref1[ref_y1 * hor + ref_x1];
			}
			else if(best_mode == BI_CONNECTED_AFF){
				mvx2 = (v22.mvx - v21.mvx)*x/xblk + (v23.mvx - v21.mvx)*y/yblk + v21.mvx;
				mvy2 = (v22.mvy - v21.mvy)*x/xblk + (v23.mvy - v21.mvy)*y/yblk + v21.mvy;

				px2 = cx - mvx2;
				py2 = cy - mvy2;

				fx2 = px2 + x;
				fy2 = py2 + y;

				ptemp = 0.5f * ( interpolate( fx1, fy1, frame_ref1, hor, ver, TYPE ) + interpolate(fx2, fy2, frame_ref2, hor, ver, TYPE ) );

				position( &ref_x2, &ref_y2, fx2, fy2, mvx2, mvy2, hor, ver );

				*aff_ref_var += 0.5f * ( frame_ref1[ref_y1 * hor + ref_x1] * frame_ref1[ref_y1 * hor + ref_x1] + 
					frame_ref2[ref_y2 * hor + ref_x2] * frame_ref2[ref_y2 * hor + ref_x2] );
				aff_ref_mean1 += frame_ref1[ref_y1 * hor + ref_x1];
				aff_ref_mean2 += frame_ref2[ref_y2 * hor + ref_x2];
			}

			if(ptemp == (float)HUGE_VAL){
//				printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);
				return ptemp;
				break;
			}

			diff = frame_cur[m] - ptemp;
			sum += diff * diff;
			m++; 
			}
		if(ptemp == (float)HUGE_VAL)break;
		m += hor - xblk;
	}                   

//	printf("affine SSE sum: %f\n",sum);
	sum /= xblk * yblk;

	(*aff_ref_var) /= xblk * yblk;
	aff_ref_mean1 /= xblk * yblk;
	aff_ref_mean2 /= xblk * yblk;

	if(best_mode == LEFT_CONNECTED_AFF || best_mode == RIGHT_CONNECTED_AFF){
		(*aff_ref_var) -= aff_ref_mean1 * aff_ref_mean1;
	}else{
		assert(best_mode == BI_CONNECTED_AFF);
		(*aff_ref_var) -= 0.5f * (aff_ref_mean1 * aff_ref_mean1 + aff_ref_mean2 * aff_ref_mean2);
	}

//	printf("affine_ref_var = %f\n",*aff_ref_var);
	
	return sum;
}

float delta_v_search_bi(float *frame_cur, float *frame_ref1, float *frame_ref2,float *upframe1, float *upframe2, int ct1x, int ct1y, int ct2x, int ct2y, ModeInfo *v11,
					 ModeInfo *v12, ModeInfo *v13, ModeInfo *v21, ModeInfo *v22, ModeInfo *v23, ModeInfo *dv11, 
					 ModeInfo *dv12, ModeInfo *dv13, ModeInfo *dv21, ModeInfo *dv22, ModeInfo *dv23, int xblk, 
					 int yblk, int cx, int cy, int hor, int ver, int best_mode, int lambda, int subpel, int t_level, int type){
  // type 0: INTER
  // type 2: MERGE UP
  // type 3: MERGE INTER
  float i,j;
  int k = 0;
  float tx1,ty1,tx2,ty2,tx3,ty3;
  int max;
  float accu = 0.25, range;

  float min_j, get_sad, code_v11,code_v12,code_v13, code_v21,code_v22,code_v23, get_bit, get_mvx, get_mvy, best_sad;
  float min_j_lst;

  if(type == 2 || type == 3){
	max = 3;
    range = accu * BI_RANGE;
  }else{
	assert(type == 0);
	max = 3;
	range = accu * BI_RANGE;
  }

  range = range + t_level * 0.25;

  dv11->mvx = v11->mvx;
  dv11->mvy = v11->mvy;
  dv12->mvx = v12->mvx;
  dv12->mvy = v12->mvy;  
  dv13->mvx = v13->mvx;
  dv13->mvy = v13->mvy;

  dv21->mvx = v21->mvx;
  dv21->mvy = v21->mvy;
  dv22->mvx = v22->mvx;
  dv22->mvy = v22->mvy;  
  dv23->mvx = v23->mvx;
  dv23->mvy = v23->mvy;

  tv11->mvx = v11->mvx;
  tv11->mvy = v11->mvy;
  tv12->mvx = v12->mvx;
  tv12->mvy = v12->mvy;  
  tv13->mvx = v13->mvx;
  tv13->mvy = v13->mvy;

  tv21->mvx = v21->mvx;
  tv21->mvy = v21->mvy;
  tv22->mvx = v22->mvx;
  tv22->mvy = v22->mvy;  
  tv23->mvx = v23->mvx;
  tv23->mvy = v23->mvy;

//  code_v1 = get_bit_cost(lambda,dv1->mvx + 0.5,dv1->mvy,v1->mvx,v1->mvy,0,0,subpel);
//  printf("code v1 = %f\n",code_v1);

  code_v11 = 0;
  code_v12 = 0;
  code_v13 = 0;

  code_v21 = 0;
  code_v22 = 0;
  code_v23 = 0;

  min_j = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
  best_sad = min_j;

//  if(type == 0)
//	assert(min_j != (float)HUGE_VAL);

  //////////////	INITIALIZATION	/////////////
  if(type == 2){
	code_v12 = get_bit_cost(lambda,dv12->mvx,dv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel);
	code_v22 = get_bit_cost(lambda,dv22->mvx,dv22->mvy,v22->mvx,v22->mvy,ct2x,ct2y,subpel);
  }else if(type == 3){
	code_v13 = get_bit_cost(lambda,dv13->mvx,dv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel);
	code_v23 = get_bit_cost(lambda,dv23->mvx,dv23->mvy,v23->mvx,v23->mvy,ct2x,ct2y,subpel);
  }else{
	assert(type == 0);
	code_v11 = get_bit_cost(lambda,dv11->mvx,dv11->mvy,v11->mvx,v11->mvy,ct1x,ct1y,subpel);
	code_v12 = get_bit_cost(lambda,dv12->mvx,dv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel);
	code_v13 = get_bit_cost(lambda,dv13->mvx,dv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel);

	code_v21 = get_bit_cost(lambda,dv21->mvx,dv21->mvy,v21->mvx,v21->mvy,ct2x,ct2y,subpel);
	code_v22 = get_bit_cost(lambda,dv22->mvx,dv22->mvy,v22->mvx,v22->mvy,ct2x,ct2y,subpel);
	code_v23 = get_bit_cost(lambda,dv23->mvx,dv23->mvy,v23->mvx,v23->mvy,ct2x,ct2y,subpel);
  }
  min_j = min_j + code_v11 + code_v12 + code_v13 + code_v21 + code_v22 + code_v23;
  min_j_lst = min_j;
  ///////////////////////////////////////////////

  while(k < max){
//V11
    if(type == 0 || type == 1){
	  tx1 = cx - dv11->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty1 = cy - dv11->mvy;

//	  if(type == 0)
//		assert(tx1 >= 0 && tx1 < hor && ty1 >= 0 && ty1 < ver);

	  get_mvx = 0;
	  get_mvy = 0; 

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx1 - i >= 0 && tx1 - i < hor && ty1 - j >= 0 && ty1 - j < ver && (i!=0 || j!=0) ){
				tv11->mvx = dv11->mvx + i;
				tv11->mvy = dv11->mvy + j;

				if(i == (-1)*range && j == (-1)*range && type == 0 )
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*tv11,*dv12,*dv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
				else
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*tv11,*dv12,*dv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,1);
				
				get_bit = get_bit_cost(lambda,tv11->mvx,tv11->mvy,v11->mvx,v11->mvy,ct1x,ct1y,subpel);

				if( (get_sad + get_bit + code_v12 + code_v13 + code_v21 + code_v22 + code_v23) < min_j ){
			      min_j = get_sad + get_bit + code_v12 + code_v13 + code_v21 + code_v22 + code_v23;
				  code_v11 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv11->mvx - dv11->mvx;
				  get_mvy = tv11->mvy - dv11->mvy;
				}
			  }
		  }
	  }
	  dv11->mvx = dv11->mvx + get_mvx;
	  dv11->mvy = dv11->mvy + get_mvy;

	  assert( code_v11 == get_bit_cost(lambda,dv11->mvx,dv11->mvy,v11->mvx,v11->mvy,ct1x,ct1y,subpel) );

	}//type 0 or 1

//V12
	if(type == 0 || type == 2){
	  tx2 = cx + xblk - dv12->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty2 = cy - dv12->mvy;

//	  if(type == 0)
//		assert(tx2 >= 0 && tx2 < hor && ty2 >= 0 && ty2 < ver);

	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx2 - i >= 0 && tx2 - i < hor && ty2 - j >= 0 && ty2 - j < ver && (i!=0 || j!=0) ){
				tv12->mvx = dv12->mvx + i;
				tv12->mvy = dv12->mvy + j;

				if(i == (-1)*range && j == (-1)*range && type == 2 )
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*tv12,*dv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
				else
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*tv12,*dv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,1);

				get_bit = get_bit_cost(lambda,tv12->mvx,tv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel);

				if( (get_sad + code_v11 + get_bit + code_v13 + code_v21 + code_v22 + code_v23) < min_j ){
			      min_j = get_sad + code_v11 + get_bit + code_v13 + code_v21 + code_v22 + code_v23;
				  code_v12 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv12->mvx - dv12->mvx;
				  get_mvy = tv12->mvy - dv12->mvy;
				}
			  }
		  }
	  }
	  dv12->mvx = dv12->mvx + get_mvx;
	  dv12->mvy = dv12->mvy + get_mvy;

	  assert( code_v12 == get_bit_cost(lambda,dv12->mvx,dv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel) );

	}//type 0 or 2
  
//V13
	if(type == 0 || type == 3){
	  tx3 = cx - dv13->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty3 = cy + yblk - dv13->mvy;

//	  if(type == 0)
//		assert(tx3 >= 0 && tx3 < hor && ty3 >= 0 && ty3 < ver);

	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu ){
		  for(j = (-1)*range;j <= range; j+= accu ){
			  if(tx3 - i >= 0 && tx3 - i < hor && ty3 - j >= 0 && ty3 - j < ver && (i!=0 || j!=0) ){
				tv13->mvx = dv13->mvx + i;
				tv13->mvy = dv13->mvy + j;

				if(i == (-1)*range && j == (-1)*range && type == 3 )
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*tv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
				else
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*tv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,1);

				get_bit = get_bit_cost(lambda,tv13->mvx,tv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel);

				if( (get_sad + code_v11 + code_v12 + get_bit + code_v21 + code_v22 + code_v23) < min_j ){
			      min_j = get_sad + code_v11 + code_v12 + get_bit + code_v21 + code_v22 + code_v23;
				  code_v13 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv13->mvx - dv13->mvx;
				  get_mvy = tv13->mvy - dv13->mvy;
				}
			  }
		  }
	  }
	  dv13->mvx = dv13->mvx + get_mvx;
	  dv13->mvy = dv13->mvy + get_mvy;

	  assert( code_v13 == get_bit_cost(lambda,dv13->mvx,dv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel) );

	}//type 0 or 3

////////
//RIGHT SEARCH
//V21
    if(type == 0 || type == 1){
	  tx1 = cx - dv21->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty1 = cy - dv21->mvy;

//	  if(type == 0)
//		assert(tx1 >= 0 && tx1 < hor && ty1 >= 0 && ty1 < ver);

	  get_mvx = 0;
	  get_mvy = 0; 

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx1 - i >= 0 && tx1 - i < hor && ty1 - j >= 0 && ty1 - j < ver && (i!=0 || j!=0) ){
				tv21->mvx = dv21->mvx + i;
				tv21->mvy = dv21->mvy + j;

				if(i == (-1)*range && j == (-1)*range )
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*tv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
				else
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*tv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,2);

				get_bit = get_bit_cost(lambda,tv21->mvx,tv21->mvy,v21->mvx,v21->mvy,ct2x,ct2y,subpel);

				if( (get_sad + code_v11 + code_v12 + code_v13 + get_bit + code_v22 + code_v23) < min_j ){
			      min_j = get_sad + code_v11 + code_v12 + code_v13 + get_bit + code_v22 + code_v23;
				  code_v21 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv21->mvx - dv21->mvx;
				  get_mvy = tv21->mvy - dv21->mvy;
				}
			  }
		  }
	  }
	  dv21->mvx = dv21->mvx + get_mvx;
	  dv21->mvy = dv21->mvy + get_mvy;

	  assert( code_v21 == get_bit_cost(lambda,dv21->mvx,dv21->mvy,v21->mvx,v21->mvy,ct2x,ct2y,subpel) );

	}//type 0 or 1

//V22
	if(type == 0 || type == 2){
	  tx2 = cx + xblk - dv22->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty2 = cy - dv22->mvy;

//	  if(type == 0)
//		assert(tx2 >= 0 && tx2 < hor && ty2 >= 0 && ty2 < ver);

	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx2 - i >= 0 && tx2 - i < hor && ty2 - j >= 0 && ty2 - j < ver && (i!=0 || j!=0) ){
				tv22->mvx = dv22->mvx + i;
				tv22->mvy = dv22->mvy + j;

				if(i == (-1)*range && j == (-1)*range && type == 2 )
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*dv21,*tv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
				else
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*dv21,*tv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,2);

				get_bit = get_bit_cost(lambda,tv22->mvx,tv22->mvy,v22->mvx,v22->mvy,ct2x,ct2y,subpel);

				if( (get_sad + code_v11 + code_v12 + code_v13 + code_v21 + get_bit + code_v23) < min_j ){
			      min_j = get_sad + code_v11 + code_v12 + code_v13 + code_v21 + get_bit + code_v23;
				  code_v22 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv22->mvx - dv22->mvx;
				  get_mvy = tv22->mvy - dv22->mvy;
				}
			  }
		  }
	  }
	  dv22->mvx = dv22->mvx + get_mvx;
	  dv22->mvy = dv22->mvy + get_mvy;

	  assert( code_v22 == get_bit_cost(lambda,dv22->mvx,dv22->mvy,v22->mvx,v22->mvy,ct2x,ct2y,subpel) );

	}//type 0 or 2
  
//V23
	if(type == 0 || type == 3){
	  tx3 = cx - dv23->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty3 = cy + yblk - dv23->mvy;

//	  if(type == 0)
//		assert(tx3 >= 0 && tx3 < hor && ty3 >= 0 && ty3 < ver);

	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu ){
		  for(j = (-1)*range;j <= range; j+= accu ){
			  if(tx3 - i >= 0 && tx3 - i < hor && ty3 - j >= 0 && ty3 - j < ver && (i!=0 || j!=0) ){
				tv23->mvx = dv23->mvx + i;
				tv23->mvy = dv23->mvy + j;


				if(i == (-1)*range && j == (-1)*range && type == 3 )
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*dv21,*dv22,*tv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
				else
					get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*dv21,*dv22,*tv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,2);
				
				get_bit = get_bit_cost(lambda,tv23->mvx,tv23->mvy,v23->mvx,v23->mvy,ct2x,ct2y,subpel);

				if( (get_sad + code_v11 + code_v12 + code_v13 + code_v21 + code_v22 + get_bit) < min_j ){
			      min_j = get_sad + code_v11 + code_v12 + code_v13 + code_v21 + code_v22 + get_bit;
				  code_v23 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv23->mvx - dv23->mvx;
				  get_mvy = tv23->mvy - dv23->mvy;
				}
			  }
		  }
	  }
	  dv23->mvx = dv23->mvx + get_mvx;
	  dv23->mvy = dv23->mvy + get_mvy;

	  assert( code_v23 == get_bit_cost(lambda,dv23->mvx,dv23->mvy,v23->mvx,v23->mvy,ct2x,ct2y,subpel) );

	}//type 0 or 3

////////
	k++;
	assert(min_j <= min_j_lst);

	if( min_j == min_j_lst ){
		k = max + 1;
		break;
	}else{
		min_j_lst = min_j;
	}

  }//if k < max

  return best_sad;
}

float delta_v_search_prl(float *frame_cur, float *frame_ref1, float *frame_ref2, float *upframe1, float *upframe2, int ct1x, int ct1y, int ct2x, int ct2y, ModeInfo *v11,
					    ModeInfo *v12, ModeInfo *v13, ModeInfo *dv11, ModeInfo *dv12, ModeInfo *dv13, ModeInfo *dv21, ModeInfo *dv22, 
					    ModeInfo *dv23, int xblk, int yblk, int cx, int cy, int hor, int ver, int best_mode, int lambda, int subpel, 
					    int t_level, int type){
  // type 0: INTER
  // type 2: MERGE UP
  // type 3: MERGE INTER
  float i,j;
  int k = 0;
  float tx11,ty11,tx12,ty12,tx13,ty13,tx21,ty21,tx22,ty22,tx23,ty23;
  int max;
  float get, accu = 0.25, range;

  float min_j, get_sad, code_v11,code_v12,code_v13, code_v21,code_v22,code_v23, get_bit, get_mvx, get_mvy, best_sad;
  float min_j_lst;

  if(type == 2 || type == 3){
	  max = 3;
	  range = MRG_RANGE * accu;
  }else{
	assert(type == 0);
	max = 3;
	range = PRL_RANGE * accu;
  }

  range = range + t_level * 0.25;

  dv11->mvx = v11->mvx;
  dv11->mvy = v11->mvy;
  dv12->mvx = v12->mvx;
  dv12->mvy = v12->mvy;  
  dv13->mvx = v13->mvx;
  dv13->mvy = v13->mvy;

  dv21->mvx = (-1)*v11->mvx;
  dv21->mvy = (-1)*v11->mvy;
  dv22->mvx = (-1)*v12->mvx;
  dv22->mvy = (-1)*v12->mvy;  
  dv23->mvx = (-1)*v13->mvx;
  dv23->mvy = (-1)*v13->mvy;

  tv11->mvx = v11->mvx;
  tv11->mvy = v11->mvy;
  tv12->mvx = v12->mvx;
  tv12->mvy = v12->mvy;  
  tv13->mvx = v13->mvx;
  tv13->mvy = v13->mvy;

  code_v11 = 0;
  code_v12 = 0;
  code_v13 = 0;

  min_j = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*dv13,*dv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
  best_sad = min_j;

  //////////////	INITIALIZATION	/////////////
  if(type == 2){
	code_v12 = get_bit_cost(lambda,dv12->mvx,dv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel);
  }else if(type == 3){
	code_v13 = get_bit_cost(lambda,dv13->mvx,dv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel);
  }else{
	assert(type == 0);
	code_v11 = get_bit_cost(lambda,dv11->mvx,dv11->mvy,v11->mvx,v11->mvy,ct1x,ct1y,subpel);
	code_v12 = get_bit_cost(lambda,dv12->mvx,dv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel);
	code_v13 = get_bit_cost(lambda,dv13->mvx,dv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel);
  }
  min_j = min_j + code_v11 + code_v12 + code_v13;
  min_j_lst = min_j;
  ///////////////////////////////////////////////

  while(k < max){
//V11
    if(type == 0 || type == 1){
	  tx11 = cx - dv11->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty11 = cy - dv11->mvy;

	  tx21 = cx - dv21->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty21 = cy - dv21->mvy;

	  get_mvx = 0;
	  get_mvy = 0; 

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx11 - i >= 0 && tx11 - i < hor && ty11 - j >= 0 && ty11 - j < ver &&
				 tx21 + i >= 0 && tx21 + i < hor && ty21 + j >= 0 && ty21 + j < ver && (i!=0 || j!=0) ){

				tv11->mvx = dv11->mvx + i;
				tv11->mvy = dv11->mvy + j;

				tv21->mvx = (-1)*tv11->mvx;
				tv21->mvy = (-1)*tv11->mvy;


				get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*tv11,*dv12,*dv13,*tv21,*dv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);

				get_bit = get_bit_cost(lambda,tv11->mvx,tv11->mvy,v11->mvx,v11->mvy,ct1x,ct1y,subpel);

				if( (get_sad + get_bit + code_v12 + code_v13) < min_j ){
			      min_j = get_sad + get_bit + code_v12 + code_v13;
				  code_v11 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv11->mvx - dv11->mvx;
				  get_mvy = tv11->mvy - dv11->mvy;
				}
			  }
		  }
	  }
	  dv11->mvx = dv11->mvx + get_mvx;
	  dv11->mvy = dv11->mvy + get_mvy;

	  dv21->mvx = (-1)*dv11->mvx;
	  dv21->mvy = (-1)*dv11->mvy;

	  assert( code_v11 == get_bit_cost(lambda,dv11->mvx,dv11->mvy,v11->mvx,v11->mvy,ct1x,ct1y,subpel) );

	}//type 0 or 1

//V12
	if(type == 0 || type == 2){
	  tx12 = cx + xblk - dv12->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty12 = cy - dv12->mvy;
	  
	  tx22 = cx + xblk - dv22->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty22 = cy - dv22->mvy;

	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx12 - i >= 0 && tx12 - i < hor && ty12 - j >= 0 && ty12 - j < ver &&
				 tx22 + i >= 0 && tx22 + i < hor && ty22 + j >= 0 && ty22 + j < ver && (i!=0 || j!=0) ){

				tv12->mvx = dv12->mvx + i;
				tv12->mvy = dv12->mvy + j;
				
				tv22->mvx = (-1)*tv12->mvx;
				tv22->mvy = (-1)*tv12->mvy;

				get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*tv12,*dv13,*dv21,*tv22,*dv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);

				get_bit = get_bit_cost(lambda,tv12->mvx,tv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel);

				if( (get_sad + code_v11 + get_bit + code_v13) < min_j ){
			      min_j = get_sad + code_v11 + get_bit + code_v13;
				  code_v12 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv12->mvx - dv12->mvx;
				  get_mvy = tv12->mvy - dv12->mvy;
				}
			  }
		  }
	  }
	  dv12->mvx = dv12->mvx + get_mvx;
	  dv12->mvy = dv12->mvy + get_mvy;
	  
	  dv22->mvx = (-1)*dv12->mvx;
	  dv22->mvy = (-1)*dv12->mvy;

	  assert( code_v12 == get_bit_cost(lambda,dv12->mvx,dv12->mvy,v12->mvx,v12->mvy,ct1x,ct1y,subpel) );

	}//type 0 or 2
  
//V13
	if(type == 0 || type == 3){
	  tx13 = cx - dv13->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty13 = cy + yblk - dv13->mvy;

	  tx23 = cx - dv23->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty23 = cy + yblk - dv23->mvy;

//	  if(type == 0)
//		assert(tx3 >= 0 && tx3 < hor && ty3 >= 0 && ty3 < ver);

	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu ){
		  for(j = (-1)*range;j <= range; j+= accu ){
			  if(tx13 - i >= 0 && tx13 - i < hor && ty13 - j >= 0 && ty13 - j < ver &&
				 tx23 + i >= 0 && tx23 + i < hor && ty23 + j >= 0 && ty23 + j < ver && (i!=0 || j!=0) ){

				tv13->mvx = dv13->mvx + i;
				tv13->mvy = dv13->mvy + j;
				
				tv23->mvx = (-1)*tv13->mvx;
				tv23->mvy = (-1)*tv13->mvy;

				get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv11,*dv12,*tv13,*dv21,*dv22,*tv23,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);

				get_bit = get_bit_cost(lambda,tv13->mvx,tv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel);

				if( (get_sad + code_v11 + code_v12 + get_bit) < min_j ){
			      min_j = get_sad + code_v11 + code_v12 + get_bit;
				  code_v13 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv13->mvx - dv13->mvx;
				  get_mvy = tv13->mvy - dv13->mvy;
				}
			  }
		  }
	  }
	  dv13->mvx = dv13->mvx + get_mvx;
	  dv13->mvy = dv13->mvy + get_mvy;
	  
	  dv23->mvx = (-1)*dv13->mvx;
	  dv23->mvy = (-1)*dv13->mvy;

	  assert( code_v13 == get_bit_cost(lambda,dv13->mvx,dv13->mvy,v13->mvx,v13->mvy,ct1x,ct1y,subpel) );

	}//type 0 or 3

	k++;
	assert(min_j <= min_j_lst);

	if( min_j == min_j_lst ){
		k = max + 1;
		break;
	}
	else{
		min_j_lst = min_j;
	}

  }//if k < max

  return best_sad;

}

float delta_v_search(float *frame_cur, float *frame_ref1, float *frame_ref2, float *upframe1, float *upframe2, int ctx, int cty, ModeInfo *v11,
					 ModeInfo *v12, ModeInfo *v13, ModeInfo *dv1, ModeInfo *dv2, ModeInfo *dv3, int xblk, 
					 int yblk, int cx, int cy, int hor, int ver, int best_mode, int lambda, int subpel, int t_level, int type){
  // type 0: INTER
  // type 2: MERGE UP
  // type 3: MERGE INTER
  float i,j;
  int k = 0;
  float tx1,ty1,tx2,ty2,tx3,ty3;
  int max;
  float get, accu = 0.25, range;

  float min_j, get_sad, code_v1,code_v2,code_v3, get_bit, get_mvx, get_mvy, best_sad;
  float min_j_lst;

//  range = range * (t_level + 1);

  if(type == 2 || type == 3){
	  max = 3;
	  range = accu * MRG_RANGE;
  }else{
	assert(type == 0);
	max = 3;
	range = accu * ITR_RANGE;
  }

  range = range + t_level * 0.25;

  dv1->mvx = v11->mvx;
  dv1->mvy = v11->mvy;
  dv2->mvx = v12->mvx;
  dv2->mvy = v12->mvy;  
  dv3->mvx = v13->mvx;
  dv3->mvy = v13->mvy;

  tv11->mvx = v11->mvx;
  tv11->mvy = v11->mvy;
  tv12->mvx = v12->mvx;
  tv12->mvy = v12->mvy;  
  tv13->mvx = v13->mvx;
  tv13->mvy = v13->mvy;

  code_v1 = 0;
  code_v2 = 0;
  code_v3 = 0;

  min_j = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv1,*dv2,*dv3,*dv1,*dv2,*dv3,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
  best_sad = min_j;

//  if(type == 0)
//	assert(min_j != (float)HUGE_VAL);

  //////////////	INITIALIZATION	/////////////
  if(type == 2){
	code_v2 = get_bit_cost(lambda,dv2->mvx,dv2->mvy,v12->mvx,v12->mvy,ctx,cty,subpel);
  }else if(type == 3){
	code_v3 = get_bit_cost(lambda,dv3->mvx,dv3->mvy,v13->mvx,v13->mvy,ctx,cty,subpel);
  }else{
	assert(type == 0);
	code_v1 = get_bit_cost(lambda,dv1->mvx,dv1->mvy,v11->mvx,v11->mvy,ctx,cty,subpel);
	code_v2 = get_bit_cost(lambda,dv2->mvx,dv2->mvy,v12->mvx,v12->mvy,ctx,cty,subpel);
	code_v3 = get_bit_cost(lambda,dv3->mvx,dv3->mvy,v13->mvx,v13->mvy,ctx,cty,subpel);
  }
  min_j = min_j + code_v1 + code_v2 + code_v3;
  min_j_lst = min_j;
  ///////////////////////////////////////////////

  while(k < max){
//V1
    if(type == 0 || type == 1){
	  tx1 = cx - dv1->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty1 = cy - dv1->mvy;

//	  if(type == 0)
//		assert(tx1 >= 0 && tx1 < hor && ty1 >= 0 && ty1 < ver);

	  get_mvx = 0;
	  get_mvy = 0; 

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx1 - i >= 0 && tx1 - i < hor && ty1 - j >= 0 && ty1 - j < ver && (i!=0 || j!=0) ){
				tv11->mvx = dv1->mvx + i;
				tv11->mvy = dv1->mvy + j;

				get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*tv11,*dv2,*dv3,*tv11,*dv2,*dv3,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);

				get_bit = get_bit_cost(lambda,tv11->mvx,tv11->mvy,v11->mvx,v11->mvy,ctx,cty,subpel);

				if( (get_sad + get_bit + code_v2 + code_v3) < min_j ){
			      min_j = get_sad + get_bit + code_v2 + code_v3;
				  code_v1 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv11->mvx - dv1->mvx;
				  get_mvy = tv11->mvy - dv1->mvy;
				}
			  }
		  }
	  }
	  dv1->mvx = dv1->mvx + get_mvx;
	  dv1->mvy = dv1->mvy + get_mvy;

	  assert( code_v1 == get_bit_cost(lambda,dv1->mvx,dv1->mvy,v11->mvx,v11->mvy,ctx,cty,subpel) );

	}//type 0 or 1

//V2
	if(type == 0 || type == 2){
	  tx2 = cx + xblk - dv2->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty2 = cy - dv2->mvy;

//	  if(type == 0)
//		assert(tx2 >= 0 && tx2 < hor && ty2 >= 0 && ty2 < ver);

	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu){
		  for(j = (-1)*range;j <= range; j+= accu){
			  if(tx2 - i >= 0 && tx2 - i < hor && ty2 - j >= 0 && ty2 - j < ver && (i!=0 || j!=0) ){
				tv12->mvx = dv2->mvx + i;
				tv12->mvy = dv2->mvy + j;

				get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv1,*tv12,*dv3,*dv1,*tv12,*dv3,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);
				
				get_bit = get_bit_cost(lambda,tv12->mvx,tv12->mvy,v12->mvx,v12->mvy,ctx,cty,subpel);

				if( (get_sad + code_v1 + get_bit + code_v3) < min_j ){
			      min_j = get_sad + code_v1 + get_bit + code_v3;
				  code_v2 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv12->mvx - dv2->mvx;
				  get_mvy = tv12->mvy - dv2->mvy;
				}
			  }
		  }
	  }
	  dv2->mvx = dv2->mvx + get_mvx;
	  dv2->mvy = dv2->mvy + get_mvy;
	   
	  assert( code_v2 == get_bit_cost(lambda,dv2->mvx,dv2->mvy,v12->mvx,v12->mvy,ctx,cty,subpel) );

	}//type 0 or 2
  
//V3
	if(type == 0 || type == 3){
	  tx3 = cx - dv3->mvx;//This is employed to confirm that searched pixel position is within the frame boundary
	  ty3 = cy + yblk - dv3->mvy;

/*	  if(type == 0){
		  if(!(tx3 >= 0 && tx3 < hor && ty3 >= 0 && ty3 < ver)){
			printf("tx3 = %f, ty3 = %f\n",tx3,ty3);
			assert(0);
		  }
	  }
*/
	  get_mvx = 0;
	  get_mvy = 0;

	  for(i = (-1)*range;i <= range; i += accu ){
		  for(j = (-1)*range;j <= range; j+= accu ){
			  if(tx3 - i >= 0 && tx3 - i < hor && ty3 - j >= 0 && ty3 - j < ver && (i!=0 || j!=0) ){
				tv13->mvx = dv3->mvx + i;
				tv13->mvy = dv3->mvy + j;

				get_sad = find_affine_SAD(frame_cur, frame_ref1, frame_ref2,upframe1,upframe2,*dv1,*dv2,*tv13,*dv1,*dv2,*tv13,xblk,yblk,cx,cy,hor,ver,best_mode,subpel,0);

				get_bit = get_bit_cost(lambda,tv13->mvx,tv13->mvy,v13->mvx,v13->mvy,ctx,cty,subpel);

				if( (get_sad + code_v1 + code_v2 + get_bit) < min_j ){
			      min_j = get_sad + code_v1 + code_v2 + get_bit;
				  code_v3 = get_bit;
				  best_sad = get_sad;

				  get_mvx = tv13->mvx - dv3->mvx;
				  get_mvy = tv13->mvy - dv3->mvy;
				}
			  }
		  }
	  }
	  dv3->mvx = dv3->mvx + get_mvx;
	  dv3->mvy = dv3->mvy + get_mvy;

	  assert( code_v3 == get_bit_cost(lambda,dv3->mvx,dv3->mvy,v13->mvx,v13->mvy,ctx,cty,subpel) );

	}//type 0 or 3

////////
	k++;

	assert(min_j <= min_j_lst);

	if( min_j == min_j_lst ){
		k = max + 1;
		break;
	}
	else{
		min_j_lst = min_j;
	}

  }//if k < max

  return best_sad;
}

// 决策对于当前块，找到最好的块模式
void find_best_mode(vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3_array, vector_ptr fmv4_array,
                    float *fr_cur, float *fr_ref1, float *fr_ref2, 
                    float *upframe1, float *upframe2,
                    int x, int y, int xblk, int yblk, int maxx, int maxy, 
                    int hor, int ver, int t_level, videoinfo info,
                    int ctx1x, int ctx1y, int ctx2x, int ctx2y, 
                    float *pmv1x, float *pmv1y, float *pmv2x, float *pmv2y, int dec)
{
  vector tmv1, tmv2;
  float min_cost, mode_coding_cost;
  int i,j,k, mode_cnt, best_mode, bst_idx1, bst_idx2,bst_aff_idx, getnum;
  int pos1, pos2;
  float px1, py1, px2, py2;
  ModeInfo *v1, *v2;

//  ModeInfo mrg11,mrg12,mrg13,mrg21,mrg22,mrg23;

  //Added on 02.02.2016
  float min1_total_cost, min1_bit_cost, min1_sad_cost, min2_total_cost, min2_bit_cost, min2_sad_cost;
  float bst_mv1x = (float)HUGE_VAL, bst_mv1y = (float)HUGE_VAL,bst_mv2x = (float)HUGE_VAL, bst_mv2y = (float)HUGE_VAL;
  float getval, bst_mse,bst_sad,bst_aff_sad;
  float mvx_buff, mvy_buff;

  float pre_pmv1x = (float)HUGE_VAL, pre_pmv1y = (float)HUGE_VAL, pre_pmv2x = (float)HUGE_VAL, pre_pmv2y = (float)HUGE_VAL;
  int xblk2, yblk2,get_xblk,get_yblk;

  float aff_coef;  //added on 11.24.2016

  float best_trans_sad, two_comp_sad;  //added on 11.10.2016
  int   best_trans_mode;  //added on 11.10.2016
  float best_trans_bit_cost;
  ModeInfo trans_mvl, trans_mvr;
  float aff_sad_cache;

  vector two_comp_tmv;  //added on 11.13.2016
  float *D1, *D2;
  int *Dx;
  int range = 5;
  float get_mvx,get_mvy;
  float dis_x,dis_y;

  float t_coef = 1;

  int do_two_comp = 0;

  xblk2 = ( x + xblk <= hor ) ? xblk : hor - x;
  yblk2 = ( y + yblk <= ver ) ? yblk : ver - y;

  tv11 = new ModeInfo;
  tv12 = new ModeInfo;
  tv13 = new ModeInfo;

  tv21 = new ModeInfo;
  tv22 = new ModeInfo;
  tv23 = new ModeInfo;

  vector_ptr mrg_left[4], mrg_right[4];

  //Added by Yuan Liu on 01.23.2016

  float get_aff_bit_cost, get_aff_sad_cost, get_aff_total_cost, best_aff_bit_cost, best_aff_sad_cost, best_aff_total_cost;

  float mvx, mvy;
  
  //Added on 01.18.2016
  float aff_mse,aff_ref_var,aff_var;
  int med_cnt1, med_cnt2;

#ifdef AFFINE_APPLY
  //Added by Yuan Liu
  int size;
  float level;
  float *get1x = new float;
  float *get1y = new float;
  float *get2x = new float;
  float *get2y = new float;
  float *get3x = new float;
  float *get3y = new float;

  float do_affine, do_inter;
  float best_sad, bili_sad;

  ModeInfo merg_aff_v1,merg_aff_v2;
  ModeInfo merg_aff2_v1,merg_aff2_v2;
  ModeInfo get_v1[4],get_v2[4],get_v3[4];	//Save reference information for generating affine MV
  ModeInfo get2_v1[4],get2_v2[4],get2_v3[4]; //Save reference information for generating affine MV

  int aff_idx1, aff_idx2;
  int pred11 = -1,pred12 = -1,pred13 = -1;
  int pred21 = -1,pred22 = -1,pred23 = -1;
  int best_idx1, best_idx2;
  int count, round, mrg_cnt1, mrg_cnt2;

  float left_dmvx, left_dmvy, right_dmvx, right_dmvy;

  ModeInfo *dv11 = new ModeInfo;
  ModeInfo *dv12 = new ModeInfo;
  ModeInfo *dv13 = new ModeInfo;

  ModeInfo *dv21 = new ModeInfo;
  ModeInfo *dv22 = new ModeInfo;
  ModeInfo *dv23 = new ModeInfo;

  int aff_num1,aff_num2;

#endif

  if(xblk >= 64)
	aff_coef = AFF_BLK_COEF_64;
  else if(xblk == 32)
	aff_coef = AFF_BLK_COEF_32;
  else
	aff_coef = SMALL_AFF_BLK_COEF;

  ///////////////	Added on 01.18.2016  /////////////////
/*
  for(i = 0;i <= 3;i ++){
	  fmv1->mrg_mvx[i] = (float)HUGE_VAL;
	  fmv1->mrg_mvy[i] = (float)HUGE_VAL;
	  if(fmv2 != NULL){
		fmv2->mrg_mvx[i] = (float)HUGE_VAL;
		fmv2->mrg_mvy[i] = (float)HUGE_VAL;
	  }
  }
*/

  for(i = 0;i <= 3;i ++){
	  mrg_left[i] = NULL;
	  mrg_right[i] = NULL;
  }

  med_cnt1 = 0;
  for(i=0;i<=3;i++){
	  if(pmv1x[i] != (float)HUGE_VAL && pmv1y[i] != (float)HUGE_VAL)
		med_cnt1++;
  }
//  printf("med_cnt1 = %d\n",med_cnt1);
  
  switch(med_cnt1){
	  case 0:
		  pre_pmv1x = 0.0;
		  pre_pmv1y = 0.0;
		  break;

	  case 1:
	  case 2:
		  assert(pmv1x[3] == (float)HUGE_VAL && pmv1y[3] == (float)HUGE_VAL);
		  if(pmv1x[0] != (float)HUGE_VAL && pmv1y[0] != (float)HUGE_VAL){
			pre_pmv1x = pmv1x[0];
			pre_pmv1y = pmv1y[0];
		  }
		  else if(pmv1x[1] != (float)HUGE_VAL && pmv1y[1] != (float)HUGE_VAL){
			pre_pmv1x = pmv1x[1];
			pre_pmv1y = pmv1y[1];
		  }
		  else if(pmv1x[2] != (float)HUGE_VAL && pmv1y[2] != (float)HUGE_VAL){
			pre_pmv1x = pmv1x[2];
			pre_pmv1y = pmv1y[2];
		  }
		  break;

	  case 3:
		  assert(pmv1x[3] == (float)HUGE_VAL && pmv1y[3] == (float)HUGE_VAL);
		  assert(pmv1x[0] != (float)HUGE_VAL && pmv1y[0] != (float)HUGE_VAL 
			  && pmv1x[1] != (float)HUGE_VAL && pmv1y[1] != (float)HUGE_VAL 
			  && pmv1x[2] != (float)HUGE_VAL && pmv1y[2] != (float)HUGE_VAL);
		  pre_pmv1x = med(pmv1x[0], pmv1x[1], pmv1x[2]);
		  pre_pmv1y = med(pmv1y[0], pmv1y[1], pmv1y[2]);
		  break;

	  case 4:
		assert(pmv1x[3] != (float)HUGE_VAL && pmv1y[3] != (float)HUGE_VAL);
		pre_pmv1x = pmv1x[3];
		pre_pmv1y = pmv1y[3];
		break;

		default:
			  assert(0);
  }
  
  if(fmv2 != NULL){

	  med_cnt2 = 0;
	  for(i=0;i<=3;i++){
		  if(pmv2x[i] != (float)HUGE_VAL && pmv2y[i] != (float)HUGE_VAL)
			med_cnt2++;
	  }

//	  printf("med_cnt2 = %d\n",med_cnt1);
  
	  switch(med_cnt2){
		  case 0:
			  pre_pmv2x = 0.0;
			  pre_pmv2y = 0.0;
			  break;

		  case 1:
		  case 2:
			  assert(pmv2x[3] == (float)HUGE_VAL && pmv2y[3] == (float)HUGE_VAL);
			  if(pmv2x[0] != (float)HUGE_VAL && pmv2y[0] != (float)HUGE_VAL){
				pre_pmv2x = pmv2x[0];
				pre_pmv2y = pmv2y[0];
			  }
			  else if(pmv2x[1] != (float)HUGE_VAL && pmv2y[1] != (float)HUGE_VAL){
				pre_pmv2x = pmv2x[1];
				pre_pmv2y = pmv2y[1];
			  }
			  else if(pmv2x[2] != (float)HUGE_VAL && pmv2y[2] != (float)HUGE_VAL){
				pre_pmv2x = pmv2x[2];
				pre_pmv2y = pmv2y[2];
			  }
			  break;

		  case 3:
			  assert(pmv2x[3] == (float)HUGE_VAL && pmv2y[3] == (float)HUGE_VAL);
			  assert(pmv2x[0] != (float)HUGE_VAL && pmv2y[0] != (float)HUGE_VAL 
				  && pmv2x[1] != (float)HUGE_VAL && pmv2y[1] != (float)HUGE_VAL 
				  && pmv2x[2] != (float)HUGE_VAL && pmv2y[2] != (float)HUGE_VAL);
			  pre_pmv2x = med(pmv2x[0], pmv2x[1], pmv2x[2]);
			  pre_pmv2y = med(pmv2y[0], pmv2y[1], pmv2y[2]);
			  break;

		  case 4:
			  assert(pmv2x[3] != (float)HUGE_VAL && pmv2y[3] != (float)HUGE_VAL);
			  pre_pmv2x = pmv2x[3];
			  pre_pmv2y = pmv2y[3];
			  break;

		  default:
			  assert(0);
	  }
  }
  ///////////////////////////////////////////////////////

  /************************************************************************
   **************************** LEFT CONNECTED ****************************
   ************************************************************************/
  // mode has not been tested before
  assert(fmv1->mode_info[LEFT_CONNECTED].is_valid == NO &&
         (fmv2 == NULL || fmv2->mode_info[LEFT_CONNECTED].is_valid == NO));

  min1_sad_cost = (float)HUGE_VAL;
  min1_bit_cost = (float)HUGE_VAL;
  min1_total_cost = (float)HUGE_VAL;
  tmv1.med_idx = -1;
  bst_mv1x = (float)HUGE_VAL, bst_mv1y = (float)HUGE_VAL, bst_idx1 = -1;


  for(i=0;i<=3;i++){
    if(pmv1x[i] != (float)HUGE_VAL && pmv1y[i] != (float)HUGE_VAL ){
	  // MV search
	  full_search_fast(&tmv1.mvx, &tmv1.mvy, fr_cur, fr_ref1, NULL,
					   upframe1, NULL,
					   0., 0., x, y, xblk, yblk, maxx, maxy, hor, ver, 
					   info.subpel[t_level], info.lambda[t_level], 
					   pmv1x[i], pmv1y[i], ctx1x, ctx1y, 
					   &tmv1.sad_cost, &tmv1.bit_cost, &tmv1.total_cost, NO); // add parallel mode. mwi

	  tmv1.med_idx = i;

	  if(min1_total_cost > tmv1.total_cost){
		min1_total_cost = tmv1.total_cost;
		min1_bit_cost = tmv1.bit_cost;
		min1_sad_cost = tmv1.sad_cost;
		bst_idx1 = tmv1.med_idx;
		bst_mv1x = tmv1.mvx;
		bst_mv1y = tmv1.mvy;
	  }
	}
  }
  
  if(min1_total_cost != (float)HUGE_VAL){
	assert(bst_idx1 >= 0 && bst_idx1 <= 3);
	assert(tmv1.med_idx >= 0 && bst_mv1x!= (float)HUGE_VAL && bst_mv1y!= (float)HUGE_VAL);
	tmv1.mvx = bst_mv1x;
	tmv1.mvy = bst_mv1y;
	tmv1.total_cost = min1_total_cost;
	tmv1.bit_cost = min1_bit_cost;
	tmv1.sad_cost = min1_sad_cost;
	tmv1.med_idx = bst_idx1;
  }

  if(med_cnt1 == 0){
	assert(tmv1.med_idx == -1 && bst_mv1x == (float)HUGE_VAL && bst_mv1y == (float)HUGE_VAL);
	// MV search
	full_search_fast(&tmv1.mvx, &tmv1.mvy, fr_cur, fr_ref1, NULL,
					 upframe1, NULL,
					 0., 0., x, y, xblk, yblk, maxx, maxy, hor, ver, 
					 info.subpel[t_level], info.lambda[t_level], 
					 0, 0, ctx1x, ctx1y, 
					 &tmv1.sad_cost, &tmv1.bit_cost, &tmv1.total_cost, NO); // add parallel mode. mwi 
	assert(tmv1.total_cost < min1_total_cost);
  }
  // check for validity
  find_MSE(&tmv1, NULL, fr_cur, fr_ref1, NULL, x, y, xblk, yblk, hor, ver,
           t_level);

  for(i=0;i<=3;i++){
	  if(pmv1x[i] != (float)HUGE_VAL && pmv1y[i] != (float)HUGE_VAL){
	    getval = get_bit_cost(info.lambda[t_level], tmv1.mvx,tmv1.mvy,pmv1x[i],pmv1y[i], ctx1x, ctx1y, info.subpel[t_level]);
		if(getval < tmv1.bit_cost){
			tmv1.med_idx = i;
			tmv1.bit_cost = getval;
			tmv1.total_cost = tmv1.bit_cost + tmv1.sad_cost;
		}
	  }
  }

  if (tmv1.lifting_mode == CONNECTED) {

    v1 = &fmv1->mode_info[LEFT_CONNECTED];
    v1->is_valid = YES;
    v1->is_predictor = YES;
    v1->mvx = tmv1.mvx;
    v1->mvy = tmv1.mvy;
    v1->sad_cost = tmv1.sad_cost;
	v1->mse = tmv1.mse;
    v1->bit_cost = tmv1.bit_cost;
    v1->total_cost = tmv1.total_cost;
	v1->med_idx = tmv1.med_idx;
    
    if (fmv2 != NULL) {
      v2 = &fmv2->mode_info[LEFT_CONNECTED];
      v2->is_valid = YES;
      v2->is_predictor = NO;
      v2->sad_cost = v1->sad_cost;
	  v2->mse = v1->mse;
      v2->bit_cost = v1->bit_cost;
      v2->total_cost = v1->total_cost;
    }
  }

  /************************************************************************
   **************************** LEFT PREDICTED ****************************
   ************************************************************************/

  // mode has not been tested before
  assert(fmv1->mode_info[LEFT_PREDICTED].is_valid == NO &&
         (fmv2 == NULL || fmv2->mode_info[LEFT_PREDICTED].is_valid == NO));
  
  v1 = &fmv1->mode_info[LEFT_PREDICTED];
  v1->is_valid = YES;
  v1->is_predictor = YES;
  v1->mvx = tmv1.mvx;
  v1->mvy = tmv1.mvy;
  v1->sad_cost = tmv1.sad_cost;
  v1->mse = tmv1.mse;
  v1->bit_cost = tmv1.bit_cost;
  v1->total_cost = tmv1.total_cost;
  v1->med_idx = tmv1.med_idx;

  if (fmv2 != NULL) {
    v2 = &fmv2->mode_info[LEFT_PREDICTED];
    v2->is_valid = YES;
    v2->is_predictor = NO;
    v2->sad_cost = v1->sad_cost;
	v2->mse = tmv1.mse;
    v2->bit_cost = v1->bit_cost;
    v2->total_cost = v1->total_cost;
  }

  // the following modes are possible in bi-directional case (2 MVFs) only
  if (fmv2 != NULL) {

    /************************************************************************
     **************************** RIGHT CONNECTED ***************************
     ************************************************************************/
  
    // mode has not been tested before
    assert(fmv1->mode_info[RIGHT_CONNECTED].is_valid == NO &&
           fmv2->mode_info[RIGHT_CONNECTED].is_valid == NO);
    
  min2_sad_cost = (float)HUGE_VAL;
  min2_bit_cost = (float)HUGE_VAL;
  min2_total_cost = (float)HUGE_VAL;
  tmv2.med_idx = -1;
  bst_mv2x = (float)HUGE_VAL, bst_mv2y = (float)HUGE_VAL, bst_idx2 = -1;

  for(i=0;i<=3;i++){
    if(pmv2x[i] != (float)HUGE_VAL && pmv2y[i] != (float)HUGE_VAL ){
	  // MV search
	  full_search_fast(&tmv2.mvx, &tmv2.mvy, fr_cur, fr_ref2, NULL, 
                     upframe2, NULL,
                     0., 0., x, y, xblk, yblk, maxx, maxy, hor, ver, 
                     info.subpel[t_level], info.lambda[t_level], 
                     pmv2x[i], pmv2y[i], ctx2x, ctx2y, 
                     &tmv2.sad_cost, &tmv2.bit_cost, &tmv2.total_cost, NO); // add parallel mode. mwi 

	   tmv2.med_idx = i;

	  if(min2_total_cost > tmv2.total_cost){
		min2_total_cost = tmv2.total_cost;
		min2_bit_cost = tmv2.bit_cost;
		min2_sad_cost = tmv2.sad_cost;
		bst_idx2 = tmv2.med_idx;
		bst_mv2x = tmv2.mvx;
		bst_mv2y = tmv2.mvy;
	  }
	}
  }
  
  if(min2_total_cost != (float)HUGE_VAL){
	assert(bst_idx2 >= 0 && bst_idx2 <= 3);
	assert(tmv2.med_idx >= 0 && bst_mv2x!= (float)HUGE_VAL && bst_mv2y!= (float)HUGE_VAL);
	tmv2.mvx = bst_mv2x;
	tmv2.mvy = bst_mv2y;
	tmv2.total_cost = min2_total_cost;
	tmv2.bit_cost = min2_bit_cost;
	tmv2.sad_cost = min2_sad_cost;
	tmv2.med_idx = bst_idx2;
  }

  if(med_cnt2 == 0){
	assert(tmv2.med_idx == -1 && bst_mv2x == (float)HUGE_VAL && bst_mv2y == (float)HUGE_VAL);
	// MV search
	full_search_fast(&tmv2.mvx, &tmv2.mvy, fr_cur, fr_ref2, NULL, 
                     upframe2, NULL,
                     0., 0., x, y, xblk, yblk, maxx, maxy, hor, ver, 
                     info.subpel[t_level], info.lambda[t_level], 
                     0, 0, ctx2x, ctx2y, 
                     &tmv2.sad_cost, &tmv2.bit_cost, &tmv2.total_cost, NO); // add parallel mode. mwi 
	assert(tmv2.total_cost < min2_total_cost);
  }
    
    // check for validity
    find_MSE(&tmv2, NULL, fr_cur, fr_ref2, NULL, x, y, xblk, yblk, hor, ver,
             t_level);

	for(i=0;i<=3;i++){
	  if(pmv2x[i] != (float)HUGE_VAL && pmv2y[i] != (float)HUGE_VAL){
	    getval = get_bit_cost(info.lambda[t_level], tmv2.mvx,tmv2.mvy,pmv2x[i],pmv2y[i], ctx2x, ctx2y, info.subpel[t_level]);
		if(getval < tmv2.bit_cost){
			tmv2.med_idx = i;
			tmv2.bit_cost = getval;
			tmv2.total_cost = tmv2.bit_cost + tmv2.sad_cost;
		}
	  }
    }
    
    if (tmv2.lifting_mode == CONNECTED) {
      
      v2 = &fmv2->mode_info[RIGHT_CONNECTED];
      v2->is_valid = YES;
      v2->is_predictor = YES;
      v2->mvx = tmv2.mvx;
      v2->mvy = tmv2.mvy;
      v2->sad_cost = tmv2.sad_cost;
	  v2->mse = tmv2.mse;
      v2->bit_cost = tmv2.bit_cost;
      v2->total_cost = tmv2.total_cost;
	  v2->med_idx = tmv2.med_idx;
      
      v1 = &fmv1->mode_info[RIGHT_CONNECTED];
      v1->is_valid = YES;
      v1->is_predictor = NO;
      v1->sad_cost = v2->sad_cost;
	  v1->mse = tmv2.mse;
      v1->bit_cost = v2->bit_cost;
      v1->total_cost = v2->total_cost;
    }
    
    
    /************************************************************************
     *************************** RIGHT PREDICTED ****************************
     ************************************************************************/
    
    // mode has not been tested before
    assert(fmv1->mode_info[RIGHT_PREDICTED].is_valid == NO &&
           fmv2->mode_info[RIGHT_PREDICTED].is_valid == NO);
    
    v2 = &fmv2->mode_info[RIGHT_PREDICTED];
    v2->is_valid = YES;
    v2->is_predictor = YES;
    v2->mvx = tmv2.mvx;
    v2->mvy = tmv2.mvy;
    v2->sad_cost = tmv2.sad_cost;
	v2->mse = tmv2.mse;
    v2->bit_cost = tmv2.bit_cost;
    v2->total_cost = tmv2.total_cost;
	v2->med_idx = tmv2.med_idx;

    v1 = &fmv1->mode_info[RIGHT_PREDICTED];
    v1->is_valid = YES;
    v1->is_predictor = NO;
    v1->sad_cost = v2->sad_cost;
	v1->mse = tmv2.mse;
    v1->bit_cost = v2->bit_cost;
    v1->total_cost = v2->total_cost;

    
    /************************************************************************
     ****************************** BI CONNECTED ****************************
     ************************************************************************/

    // mode has not been tested before
    assert(fmv1->mode_info[BI_CONNECTED].is_valid == NO &&
           fmv2->mode_info[BI_CONNECTED].is_valid == NO);

  // left vector has already been estimated, so start with right side
/*
  min2_sad_cost = (float)HUGE_VAL;
  min2_bit_cost = (float)HUGE_VAL;
  min2_total_cost = (float)HUGE_VAL;
  tmv2.med_idx = -1;
  bst_mv2x = (float)HUGE_VAL, bst_mv2y = (float)HUGE_VAL;
*/

  assert(tmv1.mvx != (float)HUGE_VAL && tmv1.mvy != (float)HUGE_VAL && tmv2.mvx != (float)HUGE_VAL && tmv2.mvy != (float)HUGE_VAL);

  for(k=0;k<=3;k++){
    if(pmv2x[k] != (float)HUGE_VAL && pmv2y[k] != (float)HUGE_VAL ){
	  // MV search
	  full_search_fast(&tmv2.mvx, &tmv2.mvy, fr_cur, fr_ref2, fr_ref1, 
                     upframe2, upframe1,
                     tmv1.mvx, tmv1.mvy, x, y, xblk, yblk, maxx, maxy, 
                     hor, ver, info.subpel[t_level], info.lambda[t_level], 
                     pmv2x[k], pmv2y[k], ctx2x, ctx2y, 
                     &tmv2.sad_cost, &tmv2.bit_cost, &tmv2.total_cost, NO); // add parallel mode. mwi 
	  if(min2_total_cost > tmv2.total_cost){
		min2_total_cost = tmv2.total_cost;
		min2_bit_cost = tmv2.bit_cost;
		min2_sad_cost = tmv2.sad_cost;
		tmv2.med_idx = k;
		bst_mv2x = tmv2.mvx;
		bst_mv2y = tmv2.mvy;
	  }
	}
  }
  
  if(min2_total_cost != (float)HUGE_VAL){
	assert(tmv2.med_idx >= 0 && bst_mv2x!= (float)HUGE_VAL && bst_mv2y!= (float)HUGE_VAL);
	tmv2.mvx = bst_mv2x;
	tmv2.mvy = bst_mv2y;
	tmv2.total_cost = min2_total_cost;
	tmv2.bit_cost = min2_bit_cost;
	tmv2.sad_cost = min2_sad_cost;
  }

  if(med_cnt2 == 0){
	assert(tmv2.med_idx == -1 && bst_mv2x == (float)HUGE_VAL && bst_mv2y == (float)HUGE_VAL);
	// MV search
	full_search_fast(&tmv2.mvx, &tmv2.mvy, fr_cur, fr_ref2, fr_ref1, 
                     upframe2, upframe1,
                     tmv1.mvx, tmv1.mvy, x, y, xblk, yblk, maxx, maxy, 
                     hor, ver, info.subpel[t_level], info.lambda[t_level], 
                     0, 0, ctx2x, ctx2y, 
                     &tmv2.sad_cost, &tmv2.bit_cost, &tmv2.total_cost, NO); // add parallel mode. mwi 
	assert(tmv2.total_cost < min2_total_cost);
  } 
    
    for (i = 1; i < BIME_ITERATIONS + 1; i++) {

      if (i % 2) {

// re-estimate left vector
/*        min1_sad_cost = (float)HUGE_VAL;
		min1_bit_cost = (float)HUGE_VAL;
		min1_total_cost = (float)HUGE_VAL;
		tmv1.med_idx = -1;
		bst_mv1x = (float)HUGE_VAL, bst_mv1y = (float)HUGE_VAL;
*/

		for(k=0;k<=3;k++){
		  if(pmv1x[k] != (float)HUGE_VAL && pmv1y[k] != (float)HUGE_VAL ){
			// MV search
			full_search_fast(&tmv1.mvx, &tmv1.mvy, fr_cur, fr_ref1, fr_ref2, 
                         upframe1, upframe2,
                         tmv2.mvx, tmv2.mvy, x, y, xblk, yblk, maxx, maxy, 
                         hor, ver, info.subpel[t_level], info.lambda[t_level], 
                         pmv1x[k], pmv1y[k], ctx1x, ctx1y, 
                         &tmv1.sad_cost, &tmv1.bit_cost, &tmv1.total_cost, NO); // add parallel mode. mwi
			if(min1_total_cost > tmv1.total_cost){
			  min1_total_cost = tmv1.total_cost;
			  min1_bit_cost = tmv1.bit_cost;
			  min1_sad_cost = tmv1.sad_cost;
			  tmv1.med_idx = k;
			  bst_mv1x = tmv1.mvx;
			  bst_mv1y = tmv1.mvy;
			}
		  }
		}
  
		if(min1_total_cost != (float)HUGE_VAL){
		  assert(tmv1.med_idx >= 0 && bst_mv1x!= (float)HUGE_VAL && bst_mv1y!= (float)HUGE_VAL);
		  tmv1.mvx = bst_mv1x;
		  tmv1.mvy = bst_mv1y;
		  tmv1.total_cost = min1_total_cost;
		  tmv1.bit_cost = min1_bit_cost;
		  tmv1.sad_cost = min1_sad_cost;
		}

		if(med_cnt1 == 0){
		  assert(tmv1.med_idx == -1 && bst_mv1x == (float)HUGE_VAL && bst_mv1y == (float)HUGE_VAL);
		  // MV search
		  full_search_fast(&tmv1.mvx, &tmv1.mvy, fr_cur, fr_ref1, fr_ref2, 
                         upframe1, upframe2,
                         tmv2.mvx, tmv2.mvy, x, y, xblk, yblk, maxx, maxy, 
                         hor, ver, info.subpel[t_level], info.lambda[t_level], 
                         0, 0, ctx1x, ctx1y, 
                         &tmv1.sad_cost, &tmv1.bit_cost, &tmv1.total_cost, NO);// add parallel mode. mwi 
		  assert(tmv1.total_cost < min1_total_cost);
		}
      } else {

// re-estimate right vector
/*      min2_sad_cost = (float)HUGE_VAL;
		min2_bit_cost = (float)HUGE_VAL;
		min2_total_cost = (float)HUGE_VAL;
		tmv2.med_idx = -1;
		bst_mv2x = (float)HUGE_VAL, bst_mv2y = (float)HUGE_VAL;
*/

		for(k=0;k<=3;k++){
		  if(pmv2x[k] != (float)HUGE_VAL && pmv2y[k] != (float)HUGE_VAL ){
			// MV search
			full_search_fast(&tmv2.mvx, &tmv2.mvy, fr_cur, fr_ref2, fr_ref1, 
                         upframe2, upframe1,
                         tmv1.mvx, tmv1.mvy, x, y, xblk, yblk, maxx, maxy, 
                         hor, ver, info.subpel[t_level], info.lambda[t_level],
                         pmv2x[k], pmv2y[k], ctx2x, ctx2y, 
                         &tmv2.sad_cost, &tmv2.bit_cost, &tmv2.total_cost, NO); // add parallel mode. mwi
			if(min2_total_cost > tmv2.total_cost){
			  min2_total_cost = tmv2.total_cost;
			  min2_bit_cost = tmv2.bit_cost;
			  min2_sad_cost = tmv2.sad_cost;
			  tmv2.med_idx = k;
			  bst_mv2x = tmv2.mvx;
			  bst_mv2y = tmv2.mvy;
			}
		  }
		}
  
		if(min2_total_cost != (float)HUGE_VAL){
		  assert(tmv2.med_idx >= 0 && bst_mv2x!= (float)HUGE_VAL && bst_mv2y!= (float)HUGE_VAL);
		  tmv2.mvx = bst_mv2x;
		  tmv2.mvy = bst_mv2y;
		  tmv2.total_cost = min2_total_cost;
		  tmv2.bit_cost = min2_bit_cost;
		  tmv2.sad_cost = min2_sad_cost;
		}

		if(med_cnt2 == 0){
		  assert(tmv2.med_idx == -1 && bst_mv2x == (float)HUGE_VAL && bst_mv2y == (float)HUGE_VAL);
		  // MV search
		  full_search_fast(&tmv2.mvx, &tmv2.mvy, fr_cur, fr_ref2, fr_ref1, 
                         upframe2, upframe1,
                         tmv1.mvx, tmv1.mvy, x, y, xblk, yblk, maxx, maxy, 
                         hor, ver, info.subpel[t_level], info.lambda[t_level],
                         0, 0, ctx2x, ctx2y, 
                         &tmv2.sad_cost, &tmv2.bit_cost, &tmv2.total_cost, NO); // add parallel mode. mwi 
		  assert(tmv2.total_cost < min2_total_cost);
		} 
      }
    }

    // check for validity
    find_MSE(&tmv1, &tmv2, fr_cur, fr_ref1, fr_ref2, x, y, xblk, yblk, 
             hor, ver, t_level);

	/////////////	Added on 02.06.2016	 /////////////////////////
	for(i=0;i<=3;i++){
	  if(pmv1x[i] != (float)HUGE_VAL && pmv1y[i] != (float)HUGE_VAL){
	    getval = get_bit_cost(info.lambda[t_level], tmv1.mvx,tmv1.mvy,pmv1x[i],pmv1y[i], ctx1x, ctx1y, info.subpel[t_level]);
		if(getval < tmv1.bit_cost){
			tmv1.med_idx = i;
			tmv1.bit_cost = getval;
			tmv1.total_cost = tmv1.bit_cost + tmv1.sad_cost;
		}
	  }
    }

	for(i=0;i<=3;i++){
	  if(pmv2x[i] != (float)HUGE_VAL && pmv2y[i] != (float)HUGE_VAL){
	    getval = get_bit_cost(info.lambda[t_level], tmv2.mvx,tmv2.mvy,pmv2x[i],pmv2y[i], ctx2x, ctx2y, info.subpel[t_level]);
		if(getval < tmv2.bit_cost){
			tmv2.med_idx = i;
			tmv2.bit_cost = getval;
			tmv2.total_cost = tmv2.bit_cost + tmv2.sad_cost;
		}
	  }
    }

    assert(tmv1.lifting_mode == tmv2.lifting_mode);
    if (tmv1.lifting_mode == CONNECTED) {

      v1 = &fmv1->mode_info[BI_CONNECTED];
      v1->is_valid = YES;
      v1->is_predictor = YES;
      v1->mvx = tmv1.mvx;
      v1->mvy = tmv1.mvy;
      v1->sad_cost = tmv2.sad_cost; // right ME was last => get right side SAD
	  v1->mse = tmv1.mse;
      v1->bit_cost = tmv1.bit_cost + tmv2.bit_cost;
      v1->total_cost = v1->sad_cost + v1->bit_cost;
	  v1->med_idx = tmv1.med_idx;

      v2 = &fmv2->mode_info[BI_CONNECTED];
      v2->is_valid = YES;
      v2->is_predictor = YES;
      v2->mvx = tmv2.mvx;
      v2->mvy = tmv2.mvy;
      v2->sad_cost = v1->sad_cost;
	  v2->mse = tmv1.mse;
      v2->bit_cost = v1->bit_cost;
      v2->total_cost = v1->total_cost;
	  v2->med_idx = tmv2.med_idx;
    }


    /************************************************************************
     ****************************** BI PREDICTED ****************************
     ************************************************************************/

      v1 = &fmv1->mode_info[BI_PREDICTED];
      v1->is_valid = YES;
      v1->is_predictor = YES;
      v1->mvx = tmv1.mvx;
      v1->mvy = tmv1.mvy;
      v1->sad_cost = tmv2.sad_cost; // right ME was last => get right side SAD
	  v1->mse = tmv1.mse;
      v1->bit_cost = tmv1.bit_cost + tmv2.bit_cost;
      v1->total_cost = v1->sad_cost + v1->bit_cost;
	  v1->med_idx = tmv1.med_idx;

      v2 = &fmv2->mode_info[BI_PREDICTED];
      v2->is_valid = YES;
      v2->is_predictor = YES;
      v2->mvx = tmv2.mvx;
      v2->mvy = tmv2.mvy;
      v2->sad_cost = v1->sad_cost;
	  v2->mse = tmv1.mse;
      v2->bit_cost = v1->bit_cost;
      v2->total_cost = v1->total_cost;
	  v2->med_idx = tmv2.med_idx;

      /************************************************************************
       ******************************* PARALLEL *******************************
       ************************************************************************/

      // mode has not been tested before
      assert(fmv1->mode_info[PARALLEL].is_valid == NO &&
             fmv2->mode_info[PARALLEL].is_valid == NO);

      // MV search
      min1_sad_cost = (float)HUGE_VAL;
	  min1_bit_cost = (float)HUGE_VAL;
	  min1_total_cost = (float)HUGE_VAL;
	  tmv1.med_idx = -1;
	  bst_mv1x = (float)HUGE_VAL, bst_mv1y = (float)HUGE_VAL;


	  for(i=0;i<=3;i++){
	    if(pmv1x[i] != (float)HUGE_VAL && pmv1y[i] != (float)HUGE_VAL ){
			// MV search
			full_search_fast(&tmv1.mvx, &tmv1.mvy, fr_cur, fr_ref1, fr_ref2, 
                         upframe1, upframe2,
                         0., 0., x, y, xblk, yblk, maxx, maxy, 
                         hor, ver, info.subpel[t_level], info.lambda[t_level], 
                         pmv1x[i], pmv1y[i], ctx1x, ctx1y, 
                         &tmv1.sad_cost, &tmv1.bit_cost, &tmv1.total_cost, YES); // add parallel mode. mwi
			if(min1_total_cost > tmv1.total_cost){
			  min1_total_cost = tmv1.total_cost;
			  min1_bit_cost = tmv1.bit_cost;
			  min1_sad_cost = tmv1.sad_cost;
			  tmv1.med_idx = i;
			  bst_mv1x = tmv1.mvx;
			  bst_mv1y = tmv1.mvy;
			}
		}
	  }
  
	  if(min1_total_cost != (float)HUGE_VAL){
		  assert(tmv1.med_idx >= 0 && bst_mv1x!= (float)HUGE_VAL && bst_mv1y!= (float)HUGE_VAL);
		  tmv1.mvx = bst_mv1x;
		  tmv1.mvy = bst_mv1y;
		  tmv1.total_cost = min1_total_cost;
		  tmv1.bit_cost = min1_bit_cost;
		  tmv1.sad_cost = min1_sad_cost;
	  }

	  if(med_cnt1 == 0){
		assert(tmv1.med_idx == -1 && bst_mv1x == (float)HUGE_VAL && bst_mv1y == (float)HUGE_VAL);
		// MV search
		full_search_fast(&tmv1.mvx, &tmv1.mvy, fr_cur, fr_ref1, fr_ref2, 
                         upframe1, upframe2,
                         0., 0., x, y, xblk, yblk, maxx, maxy, 
                         hor, ver, info.subpel[t_level], info.lambda[t_level], 
                         0, 0, ctx1x, ctx1y, 
                         &tmv1.sad_cost, &tmv1.bit_cost, &tmv1.total_cost, YES);// add parallel mode. mwi 
		assert(tmv1.total_cost < min1_total_cost);
	  }

      // parallel mode: backward mv = -forward mv
      tmv2.mvx = -tmv1.mvx;
      tmv2.mvy = -tmv1.mvy;

      // check for validity
      find_MSE(&tmv1, &tmv2, fr_cur, fr_ref1, fr_ref2, x, y, xblk, yblk, 
               hor, ver, t_level);

	  for(i=0;i<=3;i++){
		if(pmv1x[i] != (float)HUGE_VAL && pmv1y[i] != (float)HUGE_VAL){
			getval = get_bit_cost(info.lambda[t_level], tmv1.mvx,tmv1.mvy,pmv1x[i],pmv1y[i], ctx1x, ctx1y, info.subpel[t_level]);
			if(getval < tmv1.bit_cost){
				tmv1.med_idx = i;
				tmv1.bit_cost = getval;
				tmv1.total_cost = tmv1.bit_cost + tmv1.sad_cost;
			}
		}
      }

      assert(tmv1.lifting_mode == tmv2.lifting_mode);
      if (tmv1.lifting_mode == CONNECTED) {

        v1 = &fmv1->mode_info[PARALLEL];
        v1->is_valid = YES;
        v1->is_predictor = YES;
        v1->mvx = tmv1.mvx;
        v1->mvy = tmv1.mvy;
        v1->sad_cost = tmv1.sad_cost;
		v1->mse = tmv1.mse;
        v1->bit_cost = tmv1.bit_cost;
        v1->total_cost = tmv1.total_cost;
		v1->med_idx = tmv1.med_idx;
        
        v2 = &fmv2->mode_info[PARALLEL];
        v2->is_valid = YES;
        v2->is_predictor = YES;
        v2->mvx = tmv2.mvx;
        v2->mvy = tmv2.mvy;
        v2->sad_cost = v1->sad_cost;
		v2->mse = tmv1.mse;
        v2->bit_cost = v1->bit_cost;
        v2->total_cost = v1->total_cost;
      }
	  //////////////////////////////////

  } // end of: the following modes are possible in bi-directional case (2 MVFs) only

      /************************************************************************
       **************************** BLOCK MERGING  ****************************
       ************************************************************************/

	  assert( (x >= 0) && (x <= (hor - xblk2)) && (y >= 0) && (y <= (ver - yblk2)) );

	  aff_var = variance(fr_cur + y * hor + x, xblk2, yblk2, hor);

// mode has not been tested before
//      assert(fmv1->mode_info[BLOCK_MERGING].is_valid == NO &&
//             (fmv2 == NULL || fmv2->mode_info[BLOCK_MERGING].is_valid == NO));

// no MV search
// spatial mode is bidirectional: mv = pmv

	  for(i=0;i<=3;i++){
			mrg_left[i] = new vector;
			mrg_right[i] = new vector;

			clean_mrg_mv(mrg_left[i]);
			clean_mrg_mv(mrg_right[i]);
	  }

	  get_merge_mv_info(fmv1_array,fmv2_array,fmv3_array,fmv4_array,mrg_left,mrg_right,x,y, hor, ver, xblk2,yblk2,info,t_level);

	  for(i=0;i<=3;i++){
		  if(mrg_left[i] != NULL){
			if( mrg_left[i]->mvx != (float)HUGE_VAL ){
				assert(mrg_left[i]->mvy != (float)HUGE_VAL);
				assert( (float(x)-mrg_left[i]->mvx >= 0) && (float(x)-mrg_left[i]->mvx<= (hor - xblk2))
					&& (float(y)-mrg_left[i]->mvy >= 0) && (float(y)-mrg_left[i]->mvy <= (ver - yblk2)) );
			}
/*
			fmv1->mrg_mvx[i] = mrg_left[i]->mvx;
			fmv1->mrg_mvy[i] = mrg_left[i]->mvy;

			fmv1->mrg_aff_mvx1[i] = mrg_left[i]->aff1_mvx;
			fmv1->mrg_aff_mvy1[i] = mrg_left[i]->aff1_mvy;
			fmv1->mrg_aff_mvx2[i] = mrg_left[i]->aff2_mvx;
			fmv1->mrg_aff_mvy2[i] = mrg_left[i]->aff2_mvy;
			fmv1->mrg_aff_mvx3[i] = mrg_left[i]->aff3_mvx;
			fmv1->mrg_aff_mvy3[i] = mrg_left[i]->aff3_mvy;
*/
		  }

		  if(mrg_right[i] != NULL){
			if( mrg_right[i]->mvx != (float)HUGE_VAL ){
				assert(mrg_right[i]->mvy != (float)HUGE_VAL);
				assert( (float(x)-mrg_right[i]->mvx >= 0) && (float(x)-mrg_right[i]->mvx<= (hor - xblk2))
					&& (float(y)-mrg_right[i]->mvy >= 0) && (float(y)-mrg_right[i]->mvy <= (ver - yblk2)) );
			}
			if(fmv2 != NULL){
/*
				fmv2->mrg_mvx[i] = mrg_right[i]->mvx;
				fmv2->mrg_mvy[i] = mrg_right[i]->mvy;

				fmv2->mrg_aff_mvx1[i] = mrg_right[i]->aff1_mvx;
				fmv2->mrg_aff_mvy1[i] = mrg_right[i]->aff1_mvy;
				fmv2->mrg_aff_mvx2[i] = mrg_right[i]->aff2_mvx;
				fmv2->mrg_aff_mvy2[i] = mrg_right[i]->aff2_mvy;
				fmv2->mrg_aff_mvx3[i] = mrg_right[i]->aff3_mvx;
				fmv2->mrg_aff_mvy3[i] = mrg_right[i]->aff3_mvy;
*/
			}
		  }
	  }

//////////////	BLK MRG		////////////////////
	  bst_sad = (float)HUGE_VAL;
	  bst_mse = (float)HUGE_VAL;
	  bst_aff_sad = (float)HUGE_VAL;

	  bst_idx1 = -1;
	  tmv1.med_idx = -1;
	  tmv2.med_idx = -1;
	  tmv1.aff_mrg = NO;
	  tmv2.aff_mrg = NO;

	  tmv1.mvx = (float)HUGE_VAL;
	  tmv1.mvy = (float)HUGE_VAL;
	  tmv2.mvx = (float)HUGE_VAL;
	  tmv2.mvy = (float)HUGE_VAL;

	  tmv1.aff1_mvx = (float)HUGE_VAL;tmv1.aff1_mvy = (float)HUGE_VAL;
      tmv1.aff2_mvx = (float)HUGE_VAL;tmv1.aff2_mvy = (float)HUGE_VAL;
	  tmv1.aff3_mvx = (float)HUGE_VAL;tmv1.aff3_mvy = (float)HUGE_VAL;
	  tmv2.aff1_mvx = (float)HUGE_VAL;tmv2.aff1_mvy = (float)HUGE_VAL;
	  tmv2.aff2_mvx = (float)HUGE_VAL;tmv2.aff2_mvy = (float)HUGE_VAL;
	  tmv2.aff3_mvx = (float)HUGE_VAL;tmv2.aff3_mvy = (float)HUGE_VAL;

	  for( i = 0; i <= 3; i++ ){

		  getval = (float)HUGE_VAL;

		  if( (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL) ||
			  (mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL)){
			tmv1.mvx = mrg_left[i]->mvx;
			tmv1.mvy = mrg_left[i]->mvy;
			tmv2.mvx = mrg_right[i]->mvx;
			tmv2.mvy = mrg_right[i]->mvy;

			tmv1.aff1_mvx = mrg_left[i]->aff1_mvx;tmv1.aff1_mvy = mrg_left[i]->aff1_mvy;
			tmv1.aff2_mvx = mrg_left[i]->aff2_mvx;tmv1.aff2_mvy = mrg_left[i]->aff2_mvy;
			tmv1.aff3_mvx = mrg_left[i]->aff3_mvx;tmv1.aff3_mvy = mrg_left[i]->aff3_mvy;
			tmv2.aff1_mvx = mrg_right[i]->aff1_mvx;tmv2.aff1_mvy = mrg_right[i]->aff1_mvy;
			tmv2.aff2_mvx = mrg_right[i]->aff2_mvx;tmv2.aff2_mvy = mrg_right[i]->aff2_mvy;
			tmv2.aff3_mvx = mrg_right[i]->aff3_mvx;tmv2.aff3_mvy = mrg_right[i]->aff3_mvy;

		    if( (tmv1.mvx != (float)HUGE_VAL && tmv1.mvy != (float)HUGE_VAL && tmv2.mvx != (float)HUGE_VAL && tmv2.mvy != (float)HUGE_VAL) ){
			    tmv1.sad_cost = Subpel_MCP_Error4(fr_cur, x, y, xblk, yblk, 
												fr_ref1, upframe1, tmv1.mvx, tmv1.mvy,
												fr_ref2, upframe2, tmv2.mvx, tmv2.mvy,
												hor, ver, info.subpel[t_level]);
				tmv2.sad_cost = tmv1.sad_cost;
			}else if( tmv1.mvx != (float)HUGE_VAL && tmv1.mvx != (float)HUGE_VAL ){
			    tmv1.sad_cost = Subpel_MCP_Error4(fr_cur, x, y, xblk, yblk, 
												fr_ref1, upframe1, tmv1.mvx, tmv1.mvy,
												NULL, NULL, tmv2.mvx, tmv2.mvy,
												hor, ver, info.subpel[t_level]);

				tmv2.sad_cost = tmv1.sad_cost;
			}else if( tmv2.mvx != (float)HUGE_VAL && tmv2.mvx != (float)HUGE_VAL ){
			    tmv2.sad_cost = Subpel_MCP_Error4(fr_cur, x, y, xblk, yblk, 
												fr_ref2, upframe2, tmv2.mvx, tmv2.mvy,
												NULL, NULL, tmv1.mvx, tmv1.mvx,
												hor, ver, info.subpel[t_level]);

				tmv1.sad_cost = tmv2.sad_cost;
			}else
			    assert(0);

			getval = (float)HUGE_VAL;

//			if(fmv2 != NULL){
				if( (mrg_left[i]->aff1_mvx != (float)HUGE_VAL && mrg_left[i]->aff1_mvy != (float)HUGE_VAL && mrg_left[i]->aff2_mvx != (float)HUGE_VAL &&
					mrg_left[i]->aff2_mvy != (float)HUGE_VAL && mrg_left[i]->aff3_mvx != (float)HUGE_VAL && mrg_left[i]->aff3_mvy != (float)HUGE_VAL) &&
					(mrg_right[i]->aff1_mvx != (float)HUGE_VAL && mrg_right[i]->aff1_mvy != (float)HUGE_VAL && mrg_right[i]->aff2_mvx != (float)HUGE_VAL &&
					mrg_right[i]->aff2_mvy != (float)HUGE_VAL && mrg_right[i]->aff3_mvx != (float)HUGE_VAL && mrg_right[i]->aff3_mvy != (float)HUGE_VAL) ){
					tv11->mvx = mrg_left[i]->aff1_mvx;
					tv11->mvy = mrg_left[i]->aff1_mvy;
					tv12->mvx = mrg_left[i]->aff2_mvx;
					tv12->mvy = mrg_left[i]->aff2_mvy;
					tv13->mvx = mrg_left[i]->aff3_mvx;
					tv13->mvy = mrg_left[i]->aff3_mvy;

					tv21->mvx = mrg_right[i]->aff1_mvx;
					tv21->mvy = mrg_right[i]->aff1_mvy;
					tv22->mvx = mrg_right[i]->aff2_mvx;
					tv22->mvy = mrg_right[i]->aff2_mvy;
					tv23->mvx = mrg_right[i]->aff3_mvx;
					tv23->mvy = mrg_right[i]->aff3_mvy;

					getval = find_affine_SAD(fr_cur,fr_ref1,fr_ref2,upframe1,upframe2,*tv11,*tv12,*tv13,*tv21,*tv22,*tv23,xblk2,yblk2,x,y,
						hor,ver,BI_CONNECTED_AFF,info.subpel[t_level],3);

				}else if(mrg_left[i]->aff1_mvx != (float)HUGE_VAL && mrg_left[i]->aff1_mvy != (float)HUGE_VAL && mrg_left[i]->aff2_mvx != (float)HUGE_VAL &&
					mrg_left[i]->aff2_mvy != (float)HUGE_VAL && mrg_left[i]->aff3_mvx != (float)HUGE_VAL && mrg_left[i]->aff3_mvy != (float)HUGE_VAL){
					tv11->mvx = mrg_left[i]->aff1_mvx;
					tv11->mvy = mrg_left[i]->aff1_mvy;
					tv12->mvx = mrg_left[i]->aff2_mvx;
					tv12->mvy = mrg_left[i]->aff2_mvy;
					tv13->mvx = mrg_left[i]->aff3_mvx;
					tv13->mvy = mrg_left[i]->aff3_mvy;

					getval = find_affine_SAD(fr_cur,fr_ref1,fr_ref2,upframe1,upframe2,*tv11,*tv12,*tv13,*tv11,*tv12,*tv13,xblk2,yblk2,x,y,
						hor,ver,LEFT_CONNECTED_AFF,info.subpel[t_level],3);
				}else if(mrg_right[i]->aff1_mvx != (float)HUGE_VAL && mrg_right[i]->aff1_mvy != (float)HUGE_VAL && mrg_right[i]->aff2_mvx != (float)HUGE_VAL &&
					mrg_right[i]->aff2_mvy != (float)HUGE_VAL && mrg_right[i]->aff3_mvx != (float)HUGE_VAL && mrg_right[i]->aff3_mvy != (float)HUGE_VAL){
					tv21->mvx = mrg_right[i]->aff1_mvx;
					tv21->mvy = mrg_right[i]->aff1_mvy;
					tv22->mvx = mrg_right[i]->aff2_mvx;
					tv22->mvy = mrg_right[i]->aff2_mvy;
					tv23->mvx = mrg_right[i]->aff3_mvx;
					tv23->mvy = mrg_right[i]->aff3_mvy;

					getval = find_affine_SAD(fr_cur,fr_ref1,fr_ref2,upframe1,upframe2,*tv21,*tv22,*tv23,*tv21,*tv22,*tv23,xblk2,yblk2,x,y,
						hor,ver,RIGHT_CONNECTED_AFF,info.subpel[t_level],3);
				}
//			}//if fmv2 != NULL

			if(getval < bst_aff_sad){
				bst_aff_sad = getval;
				bst_idx1 = i;
			}

			if(tmv1.sad_cost < bst_sad){
				bst_sad = tmv1.sad_cost;
				tmv1.med_idx = i;
				tmv2.med_idx = i;
			}

		  }//if effective MRG target
	  }// i

	  if(bst_aff_sad < (bst_sad * aff_coef) ){
		assert(bst_idx1 >= 0);
/*		printf("cx = %d, cy = %d, xblk = %d, yblk = %d, bst_aff_sad = %f, bst_sad = %f\n",x,y,xblk2,yblk2,bst_aff_sad,bst_sad);

		if(mrg_left[bst_idx1]->aff1_mvx != (float)HUGE_VAL)
			printf("LEFT\naff1_mvx = %f, aff1_mvy = %f\naff2_mvx = %f, aff2_mvy = %f\naff3_mvx = %f, aff3_mvy = %f\n\n",
			mrg_left[bst_idx1]->aff1_mvx,mrg_left[bst_idx1]->aff1_mvy,mrg_left[bst_idx1]->aff2_mvx,
			mrg_left[bst_idx1]->aff2_mvy,mrg_left[bst_idx1]->aff3_mvx,mrg_left[bst_idx1]->aff3_mvy);
				
		if(mrg_right[bst_idx1]->aff1_mvx != (float)HUGE_VAL)
			printf("RIGHT\naff1_mvx = %f, aff1_mvy = %f\naff2_mvx = %f, aff2_mvy = %f\naff3_mvx = %f, aff3_mvy = %f\n\n",
			mrg_right[bst_idx1]->aff1_mvx,mrg_right[bst_idx1]->aff1_mvy,mrg_right[bst_idx1]->aff2_mvx,
			mrg_right[bst_idx1]->aff2_mvy,mrg_right[bst_idx1]->aff3_mvx,mrg_right[bst_idx1]->aff3_mvy);
*/
		tmv1.med_idx = bst_idx1;
		tmv2.med_idx = bst_idx1;
		tmv1.aff_mrg = YES;
		tmv2.aff_mrg = YES;
		tmv1.sad_cost = bst_aff_sad;
		tmv2.sad_cost = bst_aff_sad;
	  }else{
		tmv1.aff_mrg = NO;
		tmv2.aff_mrg = NO;
	  }

	  if(tmv1.med_idx >= 0){

		if(tmv1.aff_mrg == YES){
			tmv1.mvx = mrg_left[tmv1.med_idx]->mvx;
			tmv1.mvy = mrg_left[tmv1.med_idx]->mvy;
			tmv2.mvx = mrg_right[tmv2.med_idx]->mvx;
			tmv2.mvy = mrg_right[tmv2.med_idx]->mvy;

			tmv1.aff1_mvx = mrg_left[tmv1.med_idx]->aff1_mvx;tmv1.aff1_mvy = mrg_left[tmv1.med_idx]->aff1_mvy;
			tmv1.aff2_mvx = mrg_left[tmv1.med_idx]->aff2_mvx;tmv1.aff2_mvy = mrg_left[tmv1.med_idx]->aff2_mvy;
			tmv1.aff3_mvx = mrg_left[tmv1.med_idx]->aff3_mvx;tmv1.aff3_mvy = mrg_left[tmv1.med_idx]->aff3_mvy;
			tmv2.aff1_mvx = mrg_right[tmv2.med_idx]->aff1_mvx;tmv2.aff1_mvy = mrg_right[tmv2.med_idx]->aff1_mvy;
			tmv2.aff2_mvx = mrg_right[tmv2.med_idx]->aff2_mvx;tmv2.aff2_mvy = mrg_right[tmv2.med_idx]->aff2_mvy;
			tmv2.aff3_mvx = mrg_right[tmv2.med_idx]->aff3_mvx;tmv2.aff3_mvy = mrg_right[tmv2.med_idx]->aff3_mvy;

			tv11->mvx = mrg_left[tmv1.med_idx]->aff1_mvx;	tv21->mvx = mrg_right[tmv1.med_idx]->aff1_mvx;
			tv11->mvy = mrg_left[tmv1.med_idx]->aff1_mvy;	tv21->mvy = mrg_right[tmv1.med_idx]->aff1_mvy;
			tv12->mvx = mrg_left[tmv1.med_idx]->aff2_mvx;	tv22->mvx = mrg_right[tmv1.med_idx]->aff2_mvx;
			tv12->mvy = mrg_left[tmv2.med_idx]->aff2_mvy;	tv22->mvy = mrg_right[tmv2.med_idx]->aff2_mvy;
			tv13->mvx = mrg_left[tmv2.med_idx]->aff3_mvx;	tv23->mvx = mrg_right[tmv2.med_idx]->aff3_mvx;
			tv13->mvy = mrg_left[tmv2.med_idx]->aff3_mvy;	tv23->mvy = mrg_right[tmv2.med_idx]->aff3_mvy;
		}else{
			tmv1.mvx = mrg_left[tmv1.med_idx]->mvx;
			tmv1.mvy = mrg_left[tmv1.med_idx]->mvy;
			tmv2.mvx = mrg_right[tmv2.med_idx]->mvx;
			tmv2.mvy = mrg_right[tmv2.med_idx]->mvy;

			tmv1.aff1_mvx = (float)HUGE_VAL;tmv1.aff1_mvy = (float)HUGE_VAL;
			tmv1.aff2_mvx = (float)HUGE_VAL;tmv1.aff2_mvy = (float)HUGE_VAL;
			tmv1.aff3_mvx = (float)HUGE_VAL;tmv1.aff3_mvy = (float)HUGE_VAL;
			tmv2.aff1_mvx = (float)HUGE_VAL;tmv2.aff1_mvy = (float)HUGE_VAL;
			tmv2.aff2_mvx = (float)HUGE_VAL;tmv2.aff2_mvy = (float)HUGE_VAL;
			tmv2.aff3_mvx = (float)HUGE_VAL;tmv2.aff3_mvy = (float)HUGE_VAL;
		}
	  }
////////////////////////////////////////////////

	  for(i=0;i<=3;i++){
		delete(mrg_left[i]);
		delete(mrg_right[i]);
      }
/////////////////////////////////////
        // check for validity
	  if(tmv1.aff_mrg == NO){
			if( (tmv1.mvx != (float)HUGE_VAL && tmv1.mvy != (float)HUGE_VAL) &&
			  (tmv2.mvx != (float)HUGE_VAL && tmv2.mvy != (float)HUGE_VAL) ){
				find_MSE(&tmv1, &tmv2, fr_cur, fr_ref1, fr_ref2, x, y, xblk, yblk, hor, ver, t_level);
				assert(tmv1.mse == tmv2.mse);
				tmv1.is_predictor = YES;
				tmv2.is_predictor = YES;
			}
			else if( tmv1.mvx != (float)HUGE_VAL && tmv1.mvy != (float)HUGE_VAL ){
				find_MSE(&tmv1, NULL, fr_cur, fr_ref1, NULL, x, y, xblk, yblk, hor, ver, t_level);
				tmv2.mse = tmv1.mse;
				tmv2.lifting_mode = IGNORED;
				tmv1.is_predictor = YES;
				tmv2.is_predictor = NO;
			}
			else if( tmv2.mvx != (float)HUGE_VAL && tmv2.mvy != (float)HUGE_VAL ){
				find_MSE(&tmv2, NULL, fr_cur, fr_ref2, NULL, x, y, xblk, yblk, hor, ver, t_level);
				tmv1.mse = tmv2.mse;
				tmv1.lifting_mode = IGNORED;
				tmv1.is_predictor = NO;
				tmv2.is_predictor = YES;
			}else{
				tmv1.lifting_mode = IGNORED;
				tmv2.lifting_mode = IGNORED;
				tmv1.is_predictor = NO;
				tmv2.is_predictor = NO;
			}
	  }else{
		  assert(tmv1.aff_mrg == YES);
		  tmv1.lifting_mode = IGNORED;
		  tmv2.lifting_mode = IGNORED;
		  tmv1.is_predictor = NO;
		  tmv2.is_predictor = NO;
	  }

/////////////////////////////////////////////////////////
		if( (tmv1.aff1_mvx != (float)HUGE_VAL && tmv1.aff1_mvy != (float)HUGE_VAL) &&
		  (tmv2.aff1_mvx != (float)HUGE_VAL && tmv2.aff1_mvy != (float)HUGE_VAL) ){

			assert(tmv1.lifting_mode == IGNORED && tmv2.lifting_mode == IGNORED && tmv1.aff_mrg == YES);

			aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*tv11,*tv12,*tv13,
				*tv21,*tv22,*tv23,xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF);
			
			if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
				aff_mse < IBLOCK_FACTOR * aff_ref_var )
			   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
				assert(tmv1.mse == tmv2.mse);
				tmv1.is_predictor = YES;
				tmv2.is_predictor = YES;
				tmv1.mse = aff_mse;
				tmv2.mse = aff_mse;
				tmv1.lifting_mode = CONNECTED;
				tmv2.lifting_mode = CONNECTED;
			}
		}
		else if( tmv1.aff1_mvx != (float)HUGE_VAL && tmv1.aff1_mvy != (float)HUGE_VAL ){

			assert(tmv1.lifting_mode == IGNORED && tmv2.lifting_mode == IGNORED && tmv1.aff_mrg == YES);

			aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*tv11,*tv12,*tv13,
				*tv11,*tv12,*tv13,xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF);

			if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
				aff_mse < IBLOCK_FACTOR * aff_ref_var )
			   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
				tmv2.lifting_mode = IGNORED;
				tmv1.is_predictor = YES;
				tmv2.is_predictor = NO;
				tmv1.mse = aff_mse;
				tmv2.mse = aff_mse;
				tmv1.lifting_mode = CONNECTED;
			}
		}
		else if( tmv2.aff1_mvx != (float)HUGE_VAL && tmv2.aff1_mvy != (float)HUGE_VAL ){

			assert(tmv1.lifting_mode == IGNORED && tmv2.lifting_mode == IGNORED && tmv2.aff_mrg == YES);

			aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*tv21,*tv22,*tv23,
				*tv21,*tv22,*tv23,xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF);

			if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
				aff_mse < IBLOCK_FACTOR * aff_ref_var )
			   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
				tmv1.lifting_mode = IGNORED;
				tmv1.is_predictor = NO;
				tmv2.is_predictor = YES;
				tmv1.mse = aff_mse;
				tmv2.mse = aff_mse;
				tmv2.lifting_mode = CONNECTED;
			}
		}
/////////////////////////////////////////////////////////

        if ( (tmv1.lifting_mode == CONNECTED || tmv2.lifting_mode == CONNECTED) ) {
		  
		  assert(tmv1.med_idx >= 0 && tmv1.med_idx == tmv2.med_idx);
          // generate SAD
		  
		  if( (tmv1.mvx != (float)HUGE_VAL && tmv1.mvy != (float)HUGE_VAL && tmv2.mvx != (float)HUGE_VAL && tmv2.mvy != (float)HUGE_VAL) )
			  tmv1.sad_cost = Subpel_MCP_Error4(fr_cur, x, y, xblk, yblk, 
												fr_ref1, upframe1, tmv1.mvx, tmv1.mvy,
												fr_ref2, upframe2, tmv2.mvx, tmv2.mvy,
												hor, ver, info.subpel[t_level]);
		  else if( tmv1.mvx != (float)HUGE_VAL && tmv1.mvx != (float)HUGE_VAL )
			  tmv1.sad_cost = Subpel_MCP_Error4(fr_cur, x, y, xblk, yblk, 
												fr_ref1, upframe1, tmv1.mvx, tmv1.mvy,
												NULL, NULL, tmv2.mvx, tmv2.mvy,
												hor, ver, info.subpel[t_level]);
		  else if( tmv2.mvx != (float)HUGE_VAL && tmv2.mvx != (float)HUGE_VAL )
			  tmv1.sad_cost = Subpel_MCP_Error4(fr_cur, x, y, xblk, yblk, 
												fr_ref2, upframe2, tmv2.mvx, tmv2.mvy,
												NULL, NULL, tmv1.mvx, tmv1.mvx,
												hor, ver, info.subpel[t_level]);
		  else
			  assert(tmv1.aff1_mvx != (float)HUGE_VAL || tmv2.aff1_mvx != (float)HUGE_VAL);


          v1 = &fmv1->mode_info[BLOCK_MERGING];
          v1->is_valid = YES;
		  v1->is_predictor = tmv1.is_predictor;
		  v1->lifting_mode = tmv1.lifting_mode;
          v1->mvx = tmv1.mvx;
          v1->mvy = tmv1.mvy;
          v1->sad_cost = tmv1.sad_cost;
		  v1->mse = tmv1.mse;
          v1->bit_cost = 0.;
          v1->total_cost = tmv1.sad_cost;
		  v1->med_idx = tmv1.med_idx;
		  
		  v1->aff_mrg = tmv1.aff_mrg;

		  if(v1->aff_mrg == YES){
			v1->aff1_mvx = tmv1.aff1_mvx;
			v1->aff1_mvy = tmv1.aff1_mvy;
			v1->aff2_mvx = tmv1.aff2_mvx;
			v1->aff2_mvy = tmv1.aff2_mvy;
			v1->aff3_mvx = tmv1.aff3_mvx;
			v1->aff3_mvy = tmv1.aff3_mvy;
		  }

		  if(fmv2 != NULL){
			  v2 = &fmv2->mode_info[BLOCK_MERGING];
			  v2->is_valid = YES;
			  v2->is_predictor = tmv2.is_predictor;
			  v2->lifting_mode = tmv2.lifting_mode;
			  v2->mvx = tmv2.mvx;
			  v2->mvy = tmv2.mvy;
			  v2->sad_cost = v1->sad_cost;
			  v2->mse = tmv1.mse;
			  v2->bit_cost = v1->bit_cost;
			  v2->total_cost = v1->total_cost;
			  v2->med_idx = tmv2.med_idx;

			  v2->aff_mrg = tmv2.aff_mrg;

			  if(v2->aff_mrg == YES){
				v2->aff1_mvx = tmv2.aff1_mvx;
				v2->aff1_mvy = tmv2.aff1_mvy;
				v2->aff2_mvx = tmv2.aff2_mvx;
				v2->aff2_mvy = tmv2.aff2_mvy;
				v2->aff3_mvx = tmv2.aff3_mvx;
				v2->aff3_mvy = tmv2.aff3_mvy;
			  }
		  }
        }

  /*************************************************************************
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DECISION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  *************************************************************************/

  min_cost = (float)HUGE_VAL;
  best_mode = -1;

  for (mode_cnt = 0; mode_cnt < NUMBER_OF_TRANS_MODES; mode_cnt++) {

    // prevent certain modes
#ifdef NO_CONNECT
	  if (get_mode_coding_cost((BiMode)mode_cnt, fmv2, info.bi_mv[t_level], t_level) > 0 && mode_cnt != BI_CONNECTED
		  && mode_cnt != LEFT_CONNECTED && mode_cnt != RIGHT_CONNECTED && mode_cnt != PARALLEL && mode_cnt != BLOCK_MERGING) {
		  //if (get_mode_coding_cost((BiMode)mode_cnt, fmv2, info.bi_mv[t_level], t_level) > 0) {
#else
		  if (get_mode_coding_cost((BiMode)mode_cnt, fmv2, info.bi_mv[t_level], t_level) > 0) {
#endif

      v1 = &fmv1->mode_info[mode_cnt];
      
      if (v1->is_valid == YES) {
        assert(fmv2 != NULL || (mode_cnt != BI_CONNECTED && 
                                mode_cnt != BI_PREDICTED ));

        // add mode coding cost
        mode_coding_cost = 
          info.lambda[t_level] * get_mode_coding_cost((BiMode)mode_cnt, 
										fmv2, info.bi_mv[t_level], t_level);

        assert(mode_coding_cost > 0);
        
        v1->bit_cost += mode_coding_cost;
        v1->total_cost += mode_coding_cost;
        
        if (fmv2 != NULL) {
          v2 = &fmv2->mode_info[mode_cnt];
          
          assert(v2->is_valid == YES);
          v2->bit_cost += mode_coding_cost;
          v2->total_cost += mode_coding_cost;
        }

        // what does it cost?
        if (v1->total_cost < min_cost) {
          //          assert(fabs(v1->sad_cost + v1->bit_cost - v1->total_cost) < 1.0);
          
          best_mode = mode_cnt;
          min_cost = v1->total_cost;
        }
      }
    }
  }

  assert(best_mode >= 0);

  // best mode was found
  v1 = &fmv1->mode_info[best_mode];
  fmv1->bi_mode = (BiMode)best_mode;
  if( fmv1->bi_mode == BLOCK_MERGING ){
	fmv1->lifting_mode = v1->lifting_mode;
	fmv1->is_predictor = v1->is_predictor;

	fmv1->aff_mrg = v1->aff_mrg;
    if(fmv1->aff_mrg == YES){
  	    fmv1->aff1_mvx = v1->aff1_mvx;	fmv1->aff1_mvy = v1->aff1_mvy;
	    fmv1->aff2_mvx = v1->aff2_mvx;	fmv1->aff2_mvy = v1->aff2_mvy;
	    fmv1->aff3_mvx = v1->aff3_mvx;	fmv1->aff3_mvy = v1->aff3_mvy;
    }
  }else{
	fmv1->lifting_mode = left_mode[best_mode];
	fmv1->is_predictor = v1->is_predictor;
  }
  fmv1->mvx = v1->mvx;
  fmv1->mvy = v1->mvy;
  fmv1->sad_cost = v1->sad_cost;
  fmv1->mse = v1->mse;
  fmv1->bit_cost = v1->bit_cost;
  fmv1->total_cost = v1->total_cost;

  fmv1->med_idx = v1->med_idx;	//Added by Yuan Liu on 01.23.2016

  assert(fmv1->total_cost == min_cost);
  
  if (fmv2 != NULL) {
    v2 = &fmv2->mode_info[best_mode];
    assert(v2->is_valid == YES);
    fmv2->bi_mode = (BiMode)best_mode;
    if( fmv2->bi_mode == BLOCK_MERGING ){
	  fmv2->lifting_mode = v2->lifting_mode;
	  fmv2->is_predictor = v2->is_predictor;

	  fmv2->aff_mrg = v2->aff_mrg;
      if(fmv2->aff_mrg == YES){
	      fmv2->aff1_mvx = v2->aff1_mvx;	fmv2->aff1_mvy = v2->aff1_mvy;
    	  fmv2->aff2_mvx = v2->aff2_mvx;	fmv2->aff2_mvy = v2->aff2_mvy;
    	  fmv2->aff3_mvx = v2->aff3_mvx;	fmv2->aff3_mvy = v2->aff3_mvy;
      }
    }else{
	  fmv2->lifting_mode = right_mode[best_mode];
	  fmv2->is_predictor = v2->is_predictor;
    }
    fmv2->mvx = v2->mvx;
    fmv2->mvy = v2->mvy;
    fmv2->sad_cost = v2->sad_cost;
	fmv2->mse = fmv1->mse;
    fmv2->bit_cost = v2->bit_cost;
    fmv2->total_cost = v2->total_cost;

	fmv2->med_idx = v2->med_idx;  //Added by Yuan Liu on 01.23.2016

  }

  /****************查看是否还有connect模式******************/
  //if (fmv2 != NULL)
	 // printf("%d ", fmv2->lifting_mode);
  //if (fmv1 != NULL)
	 // printf("%d ", fmv1->lifting_mode);
  /**********************************/

  /************************************************************************
  ****************************   AFFINE MODES   ***************************
  ************************************************************************/
#ifdef AFFINE_APPLY

assert(dec == 1 || dec == 0);

size = xblk2 * yblk2;
level = pow( (float)2, (int)t_level);

do_affine = 0;
do_inter  = 0;

//Decide if we should do affine search
if(dec == 0){
	if( (fmv1->sad_cost > 1*size*(1 + t_level*0.25)) && fmv1->sad_cost > 160 )
		do_affine = 1;

	if(do_affine == 1){
		if(xblk == 8){
			if(fmv1->sad_cost > 300)
				do_inter = 1;

		}else if(fmv1->sad_cost > 1.5*size*(1 + t_level*0.25))
			do_inter = 1;
	}

}else{
	assert(dec == 1);
	if( fmv1->sad_cost > (4*size*level) )
		do_affine = 1;

	if(do_affine == 1){
		if( fmv1->sad_cost > (6*size*level) )
			do_inter = 1;
	}
}
//Affine search decision

if( (x + xblk < hor && y + yblk < ver && xblk >= 16) && do_affine == 1 && 
(best_mode == BI_CONNECTED || best_mode == BLOCK_MERGING || best_mode == PARALLEL || best_mode == LEFT_CONNECTED || best_mode == RIGHT_CONNECTED) ){


	best_trans_sad = fmv1->sad_cost;
	best_trans_mode = fmv1->bi_mode;

	aff_idx1 = fmv1->med_idx;

	if(fmv2 != NULL)
		aff_idx2 = fmv2->med_idx;

	if(best_mode != BLOCK_MERGING){
		if(best_mode != RIGHT_CONNECTED){//for LEFT_CONNECTED, BI_CONNECTED, PARALLEL
			if(fmv1->med_idx == -1){
				left_dmvx = fmv1->mvx;
				left_dmvy = fmv1->mvy;
			}else{
				left_dmvx = fmv1->mvx - pmv1x[fmv1->med_idx];
				left_dmvy = fmv1->mvy - pmv1y[fmv1->med_idx];
			}
		}

		if(best_mode != LEFT_CONNECTED && best_mode != PARALLEL){//for RIGHT_CONNECTED, RIGHT_CONNECTED
			if(fmv2->med_idx == -1){
				right_dmvx = fmv2->mvx;
				right_dmvy = fmv2->mvy;
			}else{
				right_dmvx = fmv2->mvx - pmv2x[fmv2->med_idx];
				right_dmvy = fmv2->mvy - pmv2y[fmv2->med_idx];
			}
		}
	}else{
		left_dmvx = 0;
		left_dmvy = 0;

		right_dmvx = 0;
		right_dmvy = 0;
	}

	best_trans_bit_cost = fmv1->bit_cost - info.lambda[t_level] * get_mode_coding_cost((BiMode)best_trans_mode,fmv2,info.bi_mv[t_level], t_level);

	if(best_mode == BLOCK_MERGING){
		assert(best_trans_bit_cost == 0);
	}else if(best_mode == LEFT_CONNECTED || best_mode == PARALLEL){
		getval = get_bit_cost(info.lambda[t_level],left_dmvx,left_dmvy,0,0,ctx1x,ctx1y,info.subpel[t_level]);
		assert(best_trans_bit_cost == getval);
	}else if(best_mode == RIGHT_CONNECTED){
		getval = get_bit_cost(info.lambda[t_level],right_dmvx,right_dmvy,0,0,ctx2x,ctx2y,info.subpel[t_level]);
		assert(best_trans_bit_cost == getval);
	}else{
		getval = get_bit_cost(info.lambda[t_level],left_dmvx,left_dmvy,0,0,ctx1x,ctx1y,info.subpel[t_level]) + get_bit_cost(info.lambda[t_level],right_dmvx,right_dmvy,0,0,ctx2x,ctx2y,info.subpel[t_level]);
		assert(best_trans_bit_cost == getval);
	}

	trans_mvl.mvx = fmv1->mvx;
	trans_mvl.mvy = fmv1->mvy;

	if(fmv2 != NULL){
		trans_mvr.mvx = fmv2->mvx;
		trans_mvr.mvy = fmv2->mvy;
	}

	if(fmv2 != NULL)
		assert(fmv1->bi_mode == fmv2->bi_mode && fmv1->sad_cost == fmv2->sad_cost);


//	printf("best_mode = %d, fmv1x = %f, fmv1y = %f, fmv2x = %f, fmv2y = %f, cx = %d, cy = %d, xblk2 = %d, yblk2 = %d, Translational SAD: %f\n",best_mode,fmv1->mvx,fmv1->mvy,fmv2->mvx,fmv2->mvy,x,y,xblk,yblk,fmv2->sad_cost);
//	printf("Translational MV coding bits: %f\n",fmv1->bit_cost);

//	assert(fmv1->mode_info[AFFINE_MODE].is_valid == NO &&
//		(fmv2 == NULL || fmv2->mode_info[AFFINE_MODE].is_valid == NO));

	get_v1[0].mvx = fmv1->mvx;
	get_v1[0].mvy = fmv1->mvy;
	get_v2[0].mvx = fmv1->mvx;
	get_v2[0].mvy = fmv1->mvy;
	get_v3[0].mvx = fmv1->mvx;
	get_v3[0].mvy = fmv1->mvy;

	if(fmv2!=NULL){
		get2_v1[0].mvx = fmv2->mvx;
		get2_v1[0].mvy = fmv2->mvy;
		get2_v2[0].mvx = fmv2->mvx;
		get2_v2[0].mvy = fmv2->mvy;
		get2_v3[0].mvx = fmv2->mvx;
		get2_v3[0].mvy = fmv2->mvy;
	}

	if(fmv2!=NULL)
		assert(fmv1->bi_mode == fmv2->bi_mode);

	for(i=0;i<4;i++){
		get_v1[i].mvx = (float)HUGE_VAL;
		get_v1[i].mvy = (float)HUGE_VAL;
		if(i<=3){
			get_v2[i].mvx = (float)HUGE_VAL;
			get_v2[i].mvy = (float)HUGE_VAL;
			get_v3[i].mvx = (float)HUGE_VAL;
			get_v3[i].mvy = (float)HUGE_VAL;
		}
	}

///////////
//AFFINE V1
	i = 0;

//	printf("BLOCK B\n");
	find_block(x,y-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,1);if(*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL)i++;	//BLOCK B	
//	printf("BLOCK C\n");
	find_block(x-1,y,fmv1_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,1,0);if(*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL)i++;	//BLOCK C
//	printf("BLOCK A\n");
	find_block(x-1,y-1,fmv1_array,info,t_level,get3x,get3y,&get_xblk,&get_yblk,1,1);if(*get3x != (float)HUGE_VAL && *get3y != (float)HUGE_VAL)i++;	//BLOCK A

	if( (*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL) && ( fabs(*get2x - *get1x) < SMALL_DIFF && fabs(*get2y - *get1y) < SMALL_DIFF ) ){
		i--;
		*get2x = (float)HUGE_VAL;
		*get2y = (float)HUGE_VAL;
	}
	if( (*get3x != (float)HUGE_VAL && *get3y != (float)HUGE_VAL) && ( ( fabs(*get3x - *get1x) < SMALL_DIFF && fabs(*get3y - *get1y) < SMALL_DIFF) || ( fabs(*get3x - *get2x) < SMALL_DIFF && fabs(*get3y - *get2y) < SMALL_DIFF ) ) ){
		i--;
		*get3x = (float)HUGE_VAL;
		*get3y = (float)HUGE_VAL;
	}

	get_v1[0].mvx = *get1x;
	get_v1[0].mvy = *get1y;

	get_v1[1].mvx = *get2x;
	get_v1[1].mvy = *get2y;

	get_v1[2].mvx = *get3x;
	get_v1[2].mvy = *get3y;

//TEMP Predictor
	for(i = 0;i < 3; i ++){
		if( get_v1[i].mvx == (float)HUGE_VAL && get_v1[i].mvy == (float)HUGE_VAL && fmv3_array != NULL ){
			find_block(x,y,fmv3_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,0);

			for(j = 0; j < 3; j ++){
				if( (*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL) && ( fabs(*get1x - get_v1[j].mvx) < SMALL_DIFF && fabs(*get1y - get_v1[j].mvy) < SMALL_DIFF ) ){
					*get1x = (float)HUGE_VAL;
					*get1y = (float)HUGE_VAL;
				}
			}
			get_v1[i].mvx = *get1x;
			get_v1[i].mvy = *get1y;
			break;
		}
	}

				
//AFFINE V2
	i = 0;

//	printf("BLOCK D\n");
	find_block(x+xblk2-1,y-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);if(*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL)i++;//BLOCK D
//	printf("BLOCK E\n");
	find_block(x+xblk2,y-1,fmv1_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,0,1);if(*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL)i++;//BLOCK E

	if( (*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL) && ( fabs(*get2x - *get1x) < SMALL_DIFF && fabs(*get2y - *get1y) < SMALL_DIFF ) ){
		i--;
		*get2x = (float)HUGE_VAL;
		*get2y = (float)HUGE_VAL;
	}

	get_v2[0].mvx = *get1x;
	get_v2[0].mvy = *get1y;

	get_v2[1].mvx = *get2x;
	get_v2[1].mvy = *get2y;

	assert(i<=2);

	if(i >= 2){
		get_v2[2].mvx = ( ((*get1x!=(float)HUGE_VAL)?*get1x:0) + ((*get2x!=(float)HUGE_VAL)?*get2x:0) )/i;
		get_v2[2].mvy = ( ((*get1y!=(float)HUGE_VAL)?*get1y:0) + ((*get2y!=(float)HUGE_VAL)?*get2y:0) )/i;
	}

	if( (get_v2[2].mvx != (float)HUGE_VAL && get_v2[2].mvy != (float)HUGE_VAL) && ( ( fabs(get_v2[2].mvx - get_v2[0].mvx) < SMALL_DIFF && fabs(get_v2[2].mvy - get_v2[0].mvy) < SMALL_DIFF) || ( fabs(get_v2[2].mvx - get_v2[1].mvx) < SMALL_DIFF && fabs(get_v2[2].mvy - get_v2[1].mvy) < SMALL_DIFF ) ) ){
		get_v2[2].mvx = (float)HUGE_VAL;
		get_v2[2].mvy = (float)HUGE_VAL;
	}

//TEMP Predictor
	for(i = 0;i < 3; i ++){
		if( get_v2[i].mvx == (float)HUGE_VAL && get_v2[i].mvy == (float)HUGE_VAL && fmv3_array != NULL ){
			find_block(x+xblk2,y,fmv3_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,0);

			for(j = 0; j < 3; j ++){
				if( (*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL) && ( fabs(*get1x - get_v2[j].mvx) < SMALL_DIFF && fabs(*get1y - get_v2[j].mvy) < SMALL_DIFF ) ){
					*get1x = (float)HUGE_VAL;
					*get1y = (float)HUGE_VAL;
				}
			}
			get_v2[i].mvx = *get1x;
			get_v2[i].mvy = *get1y;
			break;
		}
	}

//AFFINE V3
	i = 0;

//	printf("BLOCK F\n");
	find_block(x-1,y+yblk2-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);if(*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL)i++;//BLOCK F
//	printf("BLOCK G\n");
	find_block(x-1,y+yblk2,fmv1_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,1,0);if(*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL)i++;//BLOCK G

	if( (*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL) && ( fabs(*get2x - *get1x) < SMALL_DIFF && fabs(*get2y - *get1y) < SMALL_DIFF) ){
		i--;
		*get2x = (float)HUGE_VAL;
		*get2y = (float)HUGE_VAL;
	}

	get_v3[0].mvx = *get1x;
	get_v3[0].mvy = *get1y;

	get_v3[1].mvx = *get2x;
	get_v3[1].mvy = *get2y;
			
	assert(i<=2);

	if(i >= 2){
		get_v3[2].mvx = ( ((*get1x!=(float)HUGE_VAL)?*get1x:0) + ((*get2x!=(float)HUGE_VAL)?*get2x:0) )/i;
		get_v3[2].mvy = ( ((*get1y!=(float)HUGE_VAL)?*get1y:0) + ((*get2y!=(float)HUGE_VAL)?*get2y:0) )/i;
	}

	if( (get_v3[2].mvx != (float)HUGE_VAL && get_v3[2].mvy != (float)HUGE_VAL) && ( ( fabs(get_v3[2].mvx - get_v3[0].mvx) < SMALL_DIFF && fabs(get_v3[2].mvy - get_v3[0].mvy) < SMALL_DIFF) || ( fabs(get_v3[2].mvx - get_v3[1].mvx) < SMALL_DIFF && fabs(get_v3[2].mvy - get_v3[1].mvy) < SMALL_DIFF ) ) ){
		get_v3[2].mvx = (float)HUGE_VAL;
		get_v3[2].mvy = (float)HUGE_VAL;
	}

//TEMP Predictor
	for(i = 0;i < 3; i ++){
		if( get_v3[i].mvx == (float)HUGE_VAL && get_v3[i].mvy == (float)HUGE_VAL && fmv3_array != NULL ){
			find_block(x,y+yblk2,fmv3_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,0);

			for(j = 0; j < 3; j ++){
				if( (*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL) && ( fabs(*get1x - get_v3[j].mvx) < SMALL_DIFF && fabs(*get1y - get_v3[j].mvy) < SMALL_DIFF) ){
					*get1x = (float)HUGE_VAL;
					*get1y = (float)HUGE_VAL;
				}
			}

			get_v3[i].mvx = *get1x;
			get_v3[i].mvy = *get1y;
			break;
		}
	}

	for(i=0;i<3;i++){
		fmv1->aff1_pred_mvx[i] = get_v1[i].mvx;
		fmv1->aff1_pred_mvy[i] = get_v1[i].mvy;
		fmv1->aff2_pred_mvx[i] = get_v2[i].mvx;
		fmv1->aff2_pred_mvy[i] = get_v2[i].mvy;
		fmv1->aff3_pred_mvx[i] = get_v3[i].mvx;
		fmv1->aff3_pred_mvy[i] = get_v3[i].mvy;
	}

////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// fmv2 //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

	for(i=0;i<4;i++){
		get2_v1[i].mvx = (float)HUGE_VAL;
		get2_v1[i].mvy = (float)HUGE_VAL;
		if(i<=3){
			get2_v2[i].mvx = (float)HUGE_VAL;
			get2_v2[i].mvy = (float)HUGE_VAL;
			get2_v3[i].mvx = (float)HUGE_VAL;
			get2_v3[i].mvy = (float)HUGE_VAL;
		}
	}

  if(fmv2 != NULL){
///////////
//AFFINE V1
	i = 0;

//	printf("BLOCK B\n");
	find_block(x,y-1,fmv2_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,1);if(*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL)i++;	//BLOCK B	
//	printf("BLOCK C\n");
	find_block(x-1,y,fmv2_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,1,0);if(*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL)i++;	//BLOCK C
//	printf("BLOCK A\n");
	find_block(x-1,y-1,fmv2_array,info,t_level,get3x,get3y,&get_xblk,&get_yblk,1,1);if(*get3x != (float)HUGE_VAL && *get3y != (float)HUGE_VAL)i++;	//BLOCK A

	if( (*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL) && ( fabs(*get2x - *get1x) < SMALL_DIFF && fabs(*get2y - *get1y) < SMALL_DIFF ) ){
		i--;
		*get2x = (float)HUGE_VAL;
		*get2y = (float)HUGE_VAL;
	}
	if( (*get3x != (float)HUGE_VAL && *get3y != (float)HUGE_VAL) && ( ( fabs(*get3x - *get1x) < SMALL_DIFF && fabs(*get3y - *get1y) < SMALL_DIFF ) || ( fabs(*get3x - *get2x) < SMALL_DIFF && fabs(*get3y - *get2y) < SMALL_DIFF) ) ){
		i--;
		*get3x = (float)HUGE_VAL;
		*get3y = (float)HUGE_VAL;
	}

	get2_v1[0].mvx = *get1x;
	get2_v1[0].mvy = *get1y;

	get2_v1[1].mvx = *get2x;
	get2_v1[1].mvy = *get2y;

	get2_v1[2].mvx = *get3x;
	get2_v1[2].mvy = *get3y;

//TEMP Predictor
	for(i = 0;i < 3; i ++){
		if( get2_v1[i].mvx == (float)HUGE_VAL && get2_v1[i].mvy == (float)HUGE_VAL && fmv4_array != NULL  ){
			find_block(x,y,fmv4_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,0);

			for(j = 0; j < 3; j ++){
				if( (*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL) && ( fabs(*get1x - get2_v1[j].mvx) < SMALL_DIFF && fabs(*get1y - get2_v1[j].mvy) < SMALL_DIFF) ){
					*get1x = (float)HUGE_VAL;
					*get1y = (float)HUGE_VAL;
				}
			}
			get2_v1[i].mvx = *get1x;
			get2_v1[i].mvy = *get1y;
			break;
		}
	}
				
//AFFINE V2
	i = 0;

//	printf("BLOCK D\n");
	find_block(x+xblk2-1,y-1,fmv2_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);if(*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL)i++;//BLOCK D
//	printf("BLOCK E\n");
	find_block(x+xblk2,y-1,fmv2_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,0,1);if(*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL)i++;//BLOCK E

	if( (*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL) && ( fabs(*get2x - *get1x) < SMALL_DIFF && fabs(*get2y - *get1y) < SMALL_DIFF) ){
		i--;
		*get2x = (float)HUGE_VAL;
		*get2y = (float)HUGE_VAL;
	}

	get2_v2[0].mvx = *get1x;
	get2_v2[0].mvy = *get1y;

	get2_v2[1].mvx = *get2x;
	get2_v2[1].mvy = *get2y;

	assert(i<=2);
			
	if(i >= 2){
		get2_v2[2].mvx = ( ((*get1x!=(float)HUGE_VAL)?*get1x:0) + ((*get2x!=(float)HUGE_VAL)?*get2x:0) )/i;
		get2_v2[2].mvy = ( ((*get1y!=(float)HUGE_VAL)?*get1y:0) + ((*get2y!=(float)HUGE_VAL)?*get2y:0) )/i;
	}

	if( (get2_v2[2].mvx != (float)HUGE_VAL && get2_v2[2].mvy != (float)HUGE_VAL) && ( ( fabs(get2_v2[2].mvx - get2_v2[0].mvx) < SMALL_DIFF && fabs(get2_v2[2].mvy - get2_v2[0].mvy) < SMALL_DIFF) || ( fabs(get2_v2[2].mvx - get2_v2[1].mvx) < SMALL_DIFF && fabs(get2_v2[2].mvy - get2_v2[1].mvy) < SMALL_DIFF ) ) ){
		get2_v2[2].mvx = (float)HUGE_VAL;
		get2_v2[2].mvy = (float)HUGE_VAL;
	}

//TEMP Predictor
	for(i = 0;i < 3; i ++){
		if( get2_v2[i].mvx == (float)HUGE_VAL && get2_v2[i].mvy == (float)HUGE_VAL && fmv4_array != NULL  ){
			find_block(x+xblk2,y,fmv4_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,0);

			for(j = 0; j < 3; j ++){
				if( (*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL) && ( fabs(*get1x - get2_v2[j].mvx) < SMALL_DIFF && fabs(*get1y - get2_v2[j].mvy) < SMALL_DIFF) ){
					*get1x = (float)HUGE_VAL;
					*get1y = (float)HUGE_VAL;
				}
			}
			get2_v2[i].mvx = *get1x;
			get2_v2[i].mvy = *get1y;
			break;
		}
	}

//AFFINE V3
	i = 0;

//	printf("BLOCK F\n");
	find_block(x-1,y+yblk2-1,fmv2_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);if(*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL)i++;//BLOCK F
//	printf("BLOCK G\n");
	find_block(x-1,y+yblk2,fmv2_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,1,0);if(*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL)i++;//BLOCK G

	if( (*get2x != (float)HUGE_VAL && *get2y != (float)HUGE_VAL) && ( fabs(*get2x - *get1x) < SMALL_DIFF && fabs(*get2y - *get1y) < SMALL_DIFF) ){
		i--;
		*get2x = (float)HUGE_VAL;
		*get2y = (float)HUGE_VAL;
	}

	get2_v3[0].mvx = *get1x;
	get2_v3[0].mvy = *get1y;

	get2_v3[1].mvx = *get2x;
	get2_v3[1].mvy = *get2y;

	assert(i<=2);

	if(i >= 2){
		get2_v3[2].mvx = ( ((*get1x!=(float)HUGE_VAL)?*get1x:0) + ((*get2x!=(float)HUGE_VAL)?*get2x:0) )/i;
		get2_v3[2].mvy = ( ((*get1y!=(float)HUGE_VAL)?*get1y:0) + ((*get2y!=(float)HUGE_VAL)?*get2y:0) )/i;
	}

	if( (get2_v3[2].mvx != (float)HUGE_VAL && get2_v3[2].mvy != (float)HUGE_VAL) && ( ( fabs(get2_v3[2].mvx - get2_v3[0].mvx) < SMALL_DIFF && fabs(get2_v3[2].mvy - get2_v3[0].mvy) < SMALL_DIFF) || ( fabs(get2_v3[2].mvx - get2_v3[1].mvx) < SMALL_DIFF && fabs(get2_v3[2].mvy - get2_v3[1].mvy) < SMALL_DIFF ) ) ){
		get2_v3[2].mvx = (float)HUGE_VAL;
		get2_v3[2].mvy = (float)HUGE_VAL;
	}

//TEMP Predictor
	for(i = 0;i < 3; i ++){
		if( get2_v3[i].mvx == (float)HUGE_VAL && get2_v3[i].mvy == (float)HUGE_VAL && fmv4_array != NULL ){
			find_block(x,y+yblk2,fmv4_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,0);

			for(j = 0; j < 3; j ++){
				if( (*get1x != (float)HUGE_VAL && *get1y != (float)HUGE_VAL) && ( fabs(*get1x - get2_v3[j].mvx) < SMALL_DIFF && fabs(*get1y - get2_v3[j].mvy) < SMALL_DIFF ) ){
					*get1x = (float)HUGE_VAL;
					*get1y = (float)HUGE_VAL;
				}
			}
			get2_v3[i].mvx = *get1x;
			get2_v3[i].mvy = *get1y;
			break;
		}
	}

	for(i=0;i<3;i++){
		fmv2->aff1_pred_mvx[i] = get2_v1[i].mvx;
		fmv2->aff1_pred_mvy[i] = get2_v1[i].mvy;
		fmv2->aff2_pred_mvx[i] = get2_v2[i].mvx;
		fmv2->aff2_pred_mvy[i] = get2_v2[i].mvy;
		fmv2->aff3_pred_mvx[i] = get2_v3[i].mvx;
		fmv2->aff3_pred_mvy[i] = get2_v3[i].mvy;
	}

  }//if fmv2 != NULL

  ////////////////////	Added on 05.01.2016	////////////////////

  aff_num1 = 0;

  for(i=0;i<3;i++){
	for(j=0;j<3;j++){
		for(k=0;k<3;k++){
			if( get_v1[i].mvx!=(float)HUGE_VAL && get_v2[j].mvx!=(float)HUGE_VAL && get_v3[k].mvx!=(float)HUGE_VAL ){
				assert(get_v1[i].mvy!=(float)HUGE_VAL && get_v2[j].mvy!=(float)HUGE_VAL && get_v3[k].mvy!=(float)HUGE_VAL);
				aff_num1 ++;
			}
		}
	}
  }

//  if(aff_num1 == 0)
//	  printf("zero pred!\n");

//  printf("aff_num1 = %d\n",aff_num1);

  aff_num2 = 0;

  for(i=0;i<3;i++){
	for(j=0;j<3;j++){
		for(k=0;k<3;k++){
			if( get2_v1[i].mvx!=(float)HUGE_VAL && get2_v2[j].mvx!=(float)HUGE_VAL && get2_v3[k].mvx!=(float)HUGE_VAL ){
				assert(get2_v1[i].mvy!=(float)HUGE_VAL && get2_v2[j].mvy!=(float)HUGE_VAL && get2_v3[k].mvy!=(float)HUGE_VAL);
				aff_num2 ++;
			}
		}
	}
  }

//  if(aff_num2 == 0)
//	  printf("zero pred!\n");

//  printf("aff_num2 = %d\n",aff_num2);

  assert(aff_num1 <= 27 && aff_num2 <= 27);

  if(fmv2 != NULL){
//	printf("cx = %d, cy = %d, xblk = %d, yblk = %d, aff_num1 = %d, aff_num2 = %d\n",x,y,xblk2,yblk2,aff_num1,aff_num2);
  }
  //////////////////////////////////////////////////////////////

	best_aff_bit_cost = (float)HUGE_VAL;
	best_aff_sad_cost = (float)HUGE_VAL;
	best_aff_total_cost = best_aff_sad_cost + best_aff_bit_cost;
	best_idx1 = -1;
	count = 0;
	pos1 = -1;

//LEFT SEARCH
//AFF_DIRECT MODE
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				if( get_v1[i].mvx!=(float)HUGE_VAL && get_v2[j].mvx!=(float)HUGE_VAL && get_v3[k].mvx!=(float)HUGE_VAL &&
					get_v1[i].mvy!=(float)HUGE_VAL && get_v2[j].mvy!=(float)HUGE_VAL && get_v3[k].mvy!=(float)HUGE_VAL  ){
					get_aff_sad_cost = find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2,get_v1[i],get_v2[j],get_v3[k],get2_v1[i],get2_v2[j],get2_v3[k],
						xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF,info.subpel[t_level],0);
					get_aff_bit_cost = info.lambda[t_level] * get_mode_coding_cost(LEFT_CONNECTED_AFF, fmv2,
						info.bi_mv[t_level], t_level) + info.lambda[t_level] * ( (int)(aff_num1/4) + AFF_IDX_OFFSET + DIRECT_LEN );

					get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;
					if( (get_aff_total_cost < best_aff_total_cost) ){
						assert(get_aff_total_cost != (float)HUGE_VAL);

						best_aff_bit_cost = get_aff_bit_cost;
						best_aff_sad_cost = get_aff_sad_cost;
						best_aff_total_cost = get_aff_total_cost;
						best_idx1 = i * 9 + j * 3 + k;
						pred11 = i;
						pred12 = j;
						pred13 = k;
						pos1 = count;
//						printf("count = %d, length = %d\n",count, (int)(count/4) );
					}
					count++;
				}
			}
		}
	}

	if( ( (best_aff_total_cost < (fmv1->total_cost + 3 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
		best_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
		(fmv1->bi_mode >= 9 && (best_aff_total_cost < fmv1->total_cost ) && best_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

			aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,get_v1[pred11],get_v2[pred12],get_v3[pred13],
				get2_v1[pred21],get2_v2[pred22],get2_v3[pred23],xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF);

//			printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

		if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
            aff_mse < IBLOCK_FACTOR * aff_ref_var )
           || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

			if(fmv2 != NULL)
				assert(fmv1->total_cost == fmv2->total_cost);

			assert(pos1 >= 0 && pos1 <= 26);

			assert(best_idx1 >= 0 && best_idx1 <= 26);
			fmv1->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

			fmv1->sad_cost = best_aff_sad_cost;
			fmv1->mse = aff_mse;
			fmv1->bit_cost = best_aff_bit_cost;
			fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
			fmv1->aff_idx = pos1;
			fmv1->direct_idx = DIRECT;
			fmv1->is_predictor = YES;
			fmv1->lifting_mode = CONNECTED;
			fmv1->aff1_mvx = get_v1[pred11].mvx;fmv1->aff1_mvy = get_v1[pred11].mvy;
			fmv1->aff2_mvx = get_v2[pred12].mvx;fmv1->aff2_mvy = get_v2[pred12].mvy;
			fmv1->aff3_mvx = get_v3[pred13].mvx;fmv1->aff3_mvy = get_v3[pred13].mvy;

			if(fmv2 != NULL){
				fmv2->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

				fmv2->sad_cost = fmv1->sad_cost;
				fmv2->mse = aff_mse;
				fmv2->bit_cost = fmv1->bit_cost;
				fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
				fmv2->aff_idx = -1;
				fmv2->direct_idx = DIRECT;
				fmv2->is_predictor = NO;
				fmv2->lifting_mode = IGNORED;
				fmv2->aff1_mvx = (float)HUGE_VAL;fmv2->aff1_mvy = (float)HUGE_VAL;
				fmv2->aff2_mvx = (float)HUGE_VAL;fmv2->aff2_mvy = (float)HUGE_VAL;
				fmv2->aff3_mvx = (float)HUGE_VAL;fmv2->aff3_mvy = (float)HUGE_VAL;
			}
		}//If MSE
	}

//AFF_INTER MODE
	if(pred11 >= 0 && pred12 >= 0 && pred13 >= 0 && do_inter == 1){
		assert(pos1 >= 0);

		best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx1x, ctx1y, &get_v1[pred11], &get_v2[pred12], &get_v3[pred13],dv11,dv12,dv13,xblk2,yblk2,x,y,hor,ver,
			LEFT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);

		get_aff_sad_cost = find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, *dv11, *dv12, *dv13, *dv21, *dv22,
			*dv23, xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF,info.subpel[t_level],0);

		get_aff_bit_cost = info.lambda[t_level] * ( (int)(pos1/4) + AFF_IDX_OFFSET + DIRECT_LEN + MERGE_LEN) + get_bit_cost(info.lambda[t_level],dv11->mvx,
					dv11->mvy,get_v1[pred11].mvx,get_v1[pred11].mvy,ctx1x,ctx1y,info.subpel[t_level]) + get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,
					get_v2[pred12].mvx,get_v2[pred12].mvy,ctx1x,ctx1y,info.subpel[t_level]) + get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,
					get_v3[pred13].mvx,get_v3[pred13].mvy,ctx1x,ctx1y,info.subpel[t_level]) + info.lambda[t_level] * get_mode_coding_cost(LEFT_CONNECTED_AFF, 
					fmv2,info.bi_mv[t_level], t_level);

		get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

		aff_sad_cache = get_aff_sad_cost;

		if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
			get_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
			(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

			aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*dv11, *dv12, *dv13, *dv21, *dv22,
				*dv23, xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF);

//			printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

			if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
				aff_mse < IBLOCK_FACTOR * aff_ref_var )
			   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

				fmv1->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

				fmv1->sad_cost = get_aff_sad_cost;
				fmv1->mse = aff_mse;
				fmv1->bit_cost = get_aff_bit_cost;
				fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
				fmv1->aff_idx = pos1;
				fmv1->direct_idx = INDIRECT;
				fmv1->merge_idx = INTER;
				fmv1->is_predictor = YES;
				fmv1->lifting_mode = CONNECTED;
				fmv1->aff1_mvx = dv11->mvx;		fmv1->aff1_mvy = dv11->mvy;
				fmv1->aff2_mvx = dv12->mvx;		fmv1->aff2_mvy = dv12->mvy;
				fmv1->aff3_mvx = dv13->mvx;		fmv1->aff3_mvy = dv13->mvy;
			
				fmv1->aff1_dmvx = dv11->mvx - get_v1[pred11].mvx;
				fmv1->aff1_dmvy = dv11->mvy - get_v1[pred11].mvy;
				fmv1->aff2_dmvx = dv12->mvx - get_v2[pred12].mvx;
				fmv1->aff2_dmvy = dv12->mvy - get_v2[pred12].mvy;
				fmv1->aff3_dmvx = dv13->mvx - get_v3[pred13].mvx;
				fmv1->aff3_dmvy = dv13->mvy - get_v3[pred13].mvy;

	//			printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n",fmv1->aff2_dmvx,fmv1->aff2_dmvy);
	
				if(fmv2 != NULL){
					fmv2->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

					fmv2->sad_cost = fmv1->sad_cost;
					fmv2->mse = aff_mse;
					fmv2->bit_cost = fmv1->bit_cost;
					fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
					fmv2->aff_idx = -1;
					fmv2->direct_idx = INDIRECT;
					fmv2->merge_idx = INTER;
					fmv2->is_predictor = NO;
					fmv2->lifting_mode = IGNORED;
					fmv2->aff1_mvx = (float)HUGE_VAL;fmv2->aff1_mvy = (float)HUGE_VAL;
					fmv2->aff2_mvx = (float)HUGE_VAL;fmv2->aff2_mvy = (float)HUGE_VAL;
					fmv2->aff3_mvx = (float)HUGE_VAL;fmv2->aff3_mvy = (float)HUGE_VAL;
				}
			}//If MSE
		}
	}
//INTER

//AFF_MERGE MODE
//UP MERGE
    //	printf("BLOCK B\n");
	find_block(x,y-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,1);//BLOCK B	
	merg_aff_v1.mvx = *get1x;
	merg_aff_v1.mvy = *get1y;
	//	printf("BLOCK D\n");
	find_block(x+xblk2-1,y-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);//BLOCK D
	merg_aff_v2.mvx = *get1x;
	merg_aff_v2.mvy = *get1y;

	if(merg_aff_v1.mvx != (float)HUGE_VAL && merg_aff_v1.mvy != (float)HUGE_VAL && 
		merg_aff_v2.mvx != (float)HUGE_VAL && merg_aff_v2.mvy != (float)HUGE_VAL){

		mrg_cnt1 = 0;

		for(i = 0;i < 2; i++){		
			if(get_v3[i].mvx != (float)HUGE_VAL && get_v3[i].mvy != (float)HUGE_VAL){
				best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx1x, ctx1y, &merg_aff_v1, &merg_aff_v2, &get_v3[i], dv11,dv12,dv13,xblk2,yblk2,x,y,hor,ver,
				LEFT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,3);

				get_aff_sad_cost = best_sad;
//					find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, merg_aff_v1, merg_aff_v2, *dv13, merg_aff_v1, merg_aff_v2, *dv13, xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF,info.subpel[t_level],0);
				get_aff_bit_cost = info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_UP_LEN + 1 ) + info.lambda[t_level] * get_mode_coding_cost(LEFT_CONNECTED_AFF,fmv2,info.bi_mv[t_level],
				t_level) + get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,get_v3[i].mvx,get_v3[i].mvy,ctx1x,ctx1y,info.subpel[t_level]);

				get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

				if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
					get_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
					(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

					aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var, merg_aff_v1, merg_aff_v2, *dv13, 
						merg_aff_v1, merg_aff_v2, *dv13, xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF);

//					printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

					if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
						aff_mse < IBLOCK_FACTOR * aff_ref_var )
					   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

						fmv1->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

						fmv1->sad_cost = get_aff_sad_cost;
						fmv1->mse = aff_mse;
						fmv1->bit_cost = get_aff_bit_cost;
						fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
						fmv1->aff_idx = mrg_cnt1;
						fmv1->direct_idx = INDIRECT;
						fmv1->merge_idx = MERGE;
						fmv1->merge_dir = UP;
						fmv1->is_predictor = YES;
						fmv1->lifting_mode = CONNECTED;
						fmv1->aff1_mvx = merg_aff_v1.mvx;fmv1->aff1_mvy = merg_aff_v1.mvy;
						fmv1->aff2_mvx = merg_aff_v2.mvx;fmv1->aff2_mvy = merg_aff_v2.mvy;
						fmv1->aff3_mvx = dv13->mvx;		 fmv1->aff3_mvy = dv13->mvy;

						fmv1->aff3_dmvx = dv13->mvx - get_v3[i].mvx;
						fmv1->aff3_dmvy = dv13->mvy - get_v3[i].mvy;

	//					printf("fmv1->aff3_dmvx = %f, fmv1->aff3_dmvy = %f\n",fmv1->aff3_dmvx,fmv1->aff3_dmvy);
	
						if(fmv2 != NULL){
							fmv2->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

							fmv2->sad_cost = fmv1->sad_cost;
							fmv2->mse = aff_mse;
							fmv2->bit_cost = fmv1->bit_cost;
							fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
							fmv2->aff_idx = -1;
							fmv2->direct_idx = INDIRECT;
							fmv2->merge_idx = MERGE;
							fmv2->merge_dir = UP;
							fmv2->is_predictor = NO;
							fmv2->lifting_mode = IGNORED;
							fmv2->aff1_mvx = (float)HUGE_VAL;fmv2->aff1_mvy = (float)HUGE_VAL;
							fmv2->aff2_mvx = (float)HUGE_VAL;fmv2->aff2_mvy = (float)HUGE_VAL;
							fmv2->aff3_mvx = (float)HUGE_VAL;fmv2->aff3_mvy = (float)HUGE_VAL;
						}
					}//If MSE
				}

				mrg_cnt1 ++;
			}//IF get != HUGE_VAL
		}
	}

//LEFT MERGE
	//	printf("BLOCK C\n");
	find_block(x-1,y,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,0);//BLOCK C	
	merg_aff_v1.mvx = *get1x;
	merg_aff_v1.mvy = *get1y;
	//	printf("BLOCK F\n");
	find_block(x-1,y+yblk2-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);//BLOCK F
	merg_aff_v2.mvx = *get1x;
	merg_aff_v2.mvy = *get1y;

	if(merg_aff_v1.mvx != (float)HUGE_VAL && merg_aff_v1.mvy != (float)HUGE_VAL && 
		merg_aff_v2.mvx != (float)HUGE_VAL && merg_aff_v2.mvy != (float)HUGE_VAL){

		mrg_cnt1 = 0;

		for(i = 0;i < 2; i++){		
			if(get_v2[i].mvx != (float)HUGE_VAL && get_v2[i].mvy != (float)HUGE_VAL){
				best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2, upframe1,upframe2,ctx1x, ctx1y, &merg_aff_v1, &get_v2[i], &merg_aff_v2, dv11,dv12,dv13,xblk2,yblk2,x,y,hor,ver,
				LEFT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,2);

				get_aff_sad_cost = best_sad;
//				find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, merg_aff_v1, *dv12, merg_aff_v2, merg_aff_v1, *dv12, merg_aff_v2, xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF,info.subpel[t_level],0);
				get_aff_bit_cost = info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + info.lambda[t_level] * get_mode_coding_cost(LEFT_CONNECTED_AFF, 
					fmv2,info.bi_mv[t_level], t_level) + get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,get_v2[i].mvx,get_v2[i].mvy,
					ctx1x,ctx1y,info.subpel[t_level]);

				get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

				if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
					get_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
					(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

					aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var, merg_aff_v1, *dv12, merg_aff_v2, 
						merg_aff_v1, *dv12, merg_aff_v2, xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF);

//					printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

					if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
						aff_mse < IBLOCK_FACTOR * aff_ref_var )
					   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

						fmv1->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

						fmv1->sad_cost = get_aff_sad_cost;
						fmv1->mse = aff_mse;
						fmv1->bit_cost = get_aff_bit_cost;
						fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
						fmv1->aff_idx = mrg_cnt1;
						fmv1->direct_idx = INDIRECT;
						fmv1->merge_idx = MERGE;
						fmv1->merge_dir = LEFT;
						fmv1->is_predictor = YES;
						fmv1->lifting_mode = CONNECTED;
						fmv1->aff1_mvx = merg_aff_v1.mvx;fmv1->aff1_mvy = merg_aff_v1.mvy;
						fmv1->aff2_mvx = dv12->mvx;		 fmv1->aff2_mvy = dv12->mvy;
						fmv1->aff3_mvx = merg_aff_v2.mvx;fmv1->aff3_mvy = merg_aff_v2.mvy;

						fmv1->aff2_dmvx = dv12->mvx - get_v2[i].mvx;
						fmv1->aff2_dmvy = dv12->mvy - get_v2[i].mvy;

	//					printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n",fmv1->aff2_dmvx,fmv1->aff2_dmvy);
	
						if(fmv2 != NULL){
							fmv2->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

							fmv2->sad_cost = fmv1->sad_cost;
							fmv2->mse = aff_mse;
							fmv2->bit_cost = fmv1->bit_cost;
							fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
							fmv2->aff_idx = -1;
							fmv2->direct_idx = INDIRECT;
							fmv2->merge_idx = MERGE;
							fmv2->merge_dir = LEFT;
							fmv2->is_predictor = NO;
							fmv2->lifting_mode = IGNORED;
							fmv2->aff1_mvx = (float)HUGE_VAL;fmv2->aff1_mvy = (float)HUGE_VAL;
							fmv2->aff2_mvx = (float)HUGE_VAL;fmv2->aff2_mvy = (float)HUGE_VAL;
							fmv2->aff3_mvx = (float)HUGE_VAL;fmv2->aff3_mvy = (float)HUGE_VAL;
						}
					}//If MSE
				}

				mrg_cnt1 ++;
			}
		}// i
	}

/* TEST	*/
	if(pred11 >= 0 && pred12 >= 0 && pred13 >= 0 && do_inter == 1){
		if( (best_trans_mode == LEFT_CONNECTED || best_trans_mode == LEFT_PREDICTED || 
			(best_trans_mode == BLOCK_MERGING && trans_mvl.mvx != (float)HUGE_VAL && trans_mvr.mvx == (float)HUGE_VAL) ) && fmv2 != NULL){
			if(best_trans_sad < best_aff_sad_cost){
//				printf("left trans better, bi_mode = %d, x = %d, y = %d, xblk = %d, yblk = %d\naff_sad = %f, trans_sad = %f\n\n",
//					best_trans_mode,x,y,xblk2,yblk2,best_aff_sad_cost,best_trans_sad);

				best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx1x, ctx1y, &trans_mvl, &trans_mvl, &trans_mvl,dv11,dv12,dv13,xblk2,yblk2,x,y,hor,ver,
					LEFT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);

				get_aff_sad_cost = best_sad;

				get_aff_bit_cost = get_bit_cost(info.lambda[t_level],dv11->mvx,dv11->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]) +
					get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]) +
					get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]);

				dv11->dmvx = dv11->mvx - trans_mvl.mvx;	dv11->dmvy = dv11->mvy - trans_mvl.mvy;
				dv12->dmvx = dv12->mvx - trans_mvl.mvx;	dv12->dmvy = dv12->mvy - trans_mvl.mvy;
				dv13->dmvx = dv13->mvx - trans_mvl.mvx;	dv13->dmvy = dv13->mvy - trans_mvl.mvy;

//Apply trans-based aff result if better
				if(best_trans_mode == BLOCK_MERGING){

					get_aff_bit_cost = get_aff_bit_cost + info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + 
						info.lambda[t_level] * get_mode_coding_cost(LEFT_CONNECTED_AFF, fmv2,info.bi_mv[t_level], t_level);

				}else if(best_trans_mode == LEFT_CONNECTED || best_trans_mode == LEFT_PREDICTED){
					get_aff_bit_cost = get_aff_bit_cost + get_bit_cost(info.lambda[t_level],left_dmvx,left_dmvy,0,0,ctx1x,ctx1y,info.subpel[t_level]);

					get_aff_bit_cost = get_aff_bit_cost + info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + 
						info.lambda[t_level] * get_mode_coding_cost(LEFT_CONNECTED_AFF, fmv2,info.bi_mv[t_level], t_level);
					
				}else
					assert(0);

				get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;	

//decision making
				if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
					get_aff_sad_cost < (fmv1->sad_cost * aff_coef * t_coef) ) ||
					(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){
						
					aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*dv11, *dv12, *dv13, *dv21, *dv22,
						*dv23, xblk2,yblk2,x,y,hor,ver,LEFT_CONNECTED_AFF);

					if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
						aff_mse < IBLOCK_FACTOR * aff_ref_var )
					   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

						fmv1->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

						fmv1->sad_cost = get_aff_sad_cost;
						fmv1->mse = aff_mse;
						fmv1->bit_cost = get_aff_bit_cost;
						fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;

						fmv1->aff_idx = -1;

						fmv1->direct_idx = INDIRECT;
						fmv1->merge_idx = MERGE;
						fmv1->merge_dir = TRAN_P;
						
						if(best_trans_mode == BLOCK_MERGING)
							fmv1->trans_pred_idx = DIR;
						else{
							assert(best_trans_mode == LEFT_CONNECTED);
							fmv1->trans_pred_idx = INDIR;
						}


						fmv1->is_predictor = YES;
						fmv1->lifting_mode = CONNECTED;

						fmv1->dmvx = left_dmvx;			fmv1->dmvy = left_dmvy;
						fmv1->aff1_mvx = dv11->mvx;		fmv1->aff1_mvy = dv11->mvy;
						fmv1->aff2_mvx = dv12->mvx;		fmv1->aff2_mvy = dv12->mvy;
						fmv1->aff3_mvx = dv13->mvx;		fmv1->aff3_mvy = dv13->mvy;
			
						fmv1->aff1_dmvx = dv11->dmvx;	fmv1->aff1_dmvy = dv11->dmvy;						
						fmv1->aff2_dmvx = dv12->dmvx;	fmv1->aff2_dmvy = dv12->dmvy;
						fmv1->aff3_dmvx = dv13->dmvx;	fmv1->aff3_dmvy = dv13->dmvy;
						
	
						if(fmv2 != NULL){
							fmv2->bi_mode = (BiMode)LEFT_CONNECTED_AFF;

							fmv2->sad_cost = fmv1->sad_cost;
							fmv2->mse = aff_mse;
							fmv2->bit_cost = fmv1->bit_cost;
							fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
							fmv2->aff_idx = -1;

							fmv2->direct_idx = INDIRECT;
							fmv2->merge_idx = MERGE;
							fmv2->merge_dir = TRAN_P;

							if(best_trans_mode == BLOCK_MERGING)
								fmv2->trans_pred_idx = DIR;
							else{
								assert(best_trans_mode == LEFT_CONNECTED);
								fmv2->trans_pred_idx = INDIR;
							}

							fmv2->is_predictor = NO;
							fmv2->lifting_mode = IGNORED;

							fmv2->dmvx = (float)HUGE_VAL;		fmv2->dmvy = (float)HUGE_VAL;
							fmv2->aff1_mvx = (float)HUGE_VAL;	fmv2->aff1_mvy = (float)HUGE_VAL;
							fmv2->aff2_mvx = (float)HUGE_VAL;	fmv2->aff2_mvy = (float)HUGE_VAL;
							fmv2->aff3_mvx = (float)HUGE_VAL;	fmv2->aff3_mvy = (float)HUGE_VAL;
						}
					}//If MSE
				}
			}
		}
	}//if pred
/*	END OF TEST	*/


//RIGHT SEARCH IF AVAILABLE
	if( fmv2 != NULL ){
	  block_buff1 = ( float * )getarray( xblk * yblk, sizeof( float ), "block_buff1" );
      block_buff2 = ( float * )getarray( xblk * yblk, sizeof( float ), "block_buff2" );

//DIRECT MODE
		best_aff_bit_cost = (float)HUGE_VAL;
		best_aff_sad_cost = (float)HUGE_VAL;
		best_aff_total_cost = best_aff_sad_cost + best_aff_bit_cost;
		best_idx2 = -1;
		count = 0;
		pos2 = -1;

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					if( get2_v1[i].mvx!=(float)HUGE_VAL && get2_v2[j].mvx!=(float)HUGE_VAL && get2_v3[k].mvx!=(float)HUGE_VAL &&
						get2_v1[i].mvy!=(float)HUGE_VAL && get2_v2[j].mvy!=(float)HUGE_VAL && get2_v3[k].mvy!=(float)HUGE_VAL ){
						get_aff_sad_cost = find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2,get2_v1[i],get2_v2[j],get2_v3[k],get_v1[i],get_v2[j],get_v3[k],
							xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF,info.subpel[t_level],0);
						get_aff_bit_cost = info.lambda[t_level] * get_mode_coding_cost(RIGHT_CONNECTED_AFF, fmv2,
							info.bi_mv[t_level], t_level) + info.lambda[t_level] * ( (int)(aff_num2/4) + AFF_IDX_OFFSET + DIRECT_LEN );
						get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;
						if( (get_aff_total_cost < best_aff_total_cost) ){
							assert(get_aff_total_cost != (float)HUGE_VAL);

							best_aff_bit_cost = get_aff_bit_cost;
							best_aff_sad_cost = get_aff_sad_cost;
							best_aff_total_cost = get_aff_total_cost;
							best_idx2 = i * 9 + j * 3 + k;
							pred21 = i;
							pred22 = j;
							pred23 = k;
							pos2 = count;
//							printf("count = %d, length = %d\n",count, (int)(count/4) );
						}
						count++;
					}
				}
			}
		}

		if( ( (best_aff_total_cost < (fmv1->total_cost + 3 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
			best_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
			(fmv1->bi_mode >= 9 && (best_aff_total_cost < fmv1->total_cost ) && best_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

			aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var, get2_v1[pred21],get2_v2[pred22],get2_v3[pred23],
				get_v1[pred11],get_v2[pred12],get_v3[pred13], xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF);

//			printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

			if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
				aff_mse < IBLOCK_FACTOR * aff_ref_var )
			  || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
				assert(fmv1->total_cost == fmv2->total_cost);
				assert(pos2 >= 0 && pos2 <= 26);
				assert(best_idx2 >= 0 && best_idx2 <= 26);
				fmv2->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;

				fmv2->sad_cost = best_aff_sad_cost;
				fmv2->mse = aff_mse;
				fmv2->bit_cost = best_aff_bit_cost;
				fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
				fmv2->aff_idx = pos2;
				fmv2->direct_idx = DIRECT;
				fmv2->is_predictor = YES;
				fmv2->lifting_mode = CONNECTED;
				fmv2->aff1_mvx = get2_v1[pred21].mvx;fmv2->aff1_mvy = get2_v1[pred21].mvy;
				fmv2->aff2_mvx = get2_v2[pred22].mvx;fmv2->aff2_mvy = get2_v2[pred22].mvy;
				fmv2->aff3_mvx = get2_v3[pred23].mvx;fmv2->aff3_mvy = get2_v3[pred23].mvy;

				fmv1->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;
		
				fmv1->sad_cost = fmv2->sad_cost;
				fmv1->mse = aff_mse;
				fmv1->bit_cost = fmv2->bit_cost;
				fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
				fmv1->aff_idx = -1;
				fmv1->direct_idx = DIRECT;
				fmv1->is_predictor = NO;
				fmv1->lifting_mode = IGNORED;
				fmv1->aff1_mvx = (float)HUGE_VAL;fmv1->aff1_mvy = (float)HUGE_VAL;
				fmv1->aff2_mvx = (float)HUGE_VAL;fmv1->aff2_mvy = (float)HUGE_VAL;
				fmv1->aff3_mvx = (float)HUGE_VAL;fmv1->aff3_mvy = (float)HUGE_VAL;
			}//If MSE
		}

//AFF_INTER MODE
	  if(pred21 >= 0 && pred22 >= 0 && pred23 >= 0 && do_inter == 1){
		assert(pos2 >= 0);

		best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2, upframe1,upframe2,ctx2x, ctx2y, &get2_v1[pred21], &get2_v2[pred22], &get2_v3[pred23],dv21,dv22,dv23,xblk2,yblk2,x,y,
			hor,ver,RIGHT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);

		get_aff_sad_cost = best_sad;
//		find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, *dv21, *dv22, *dv23, *dv11, *dv12,*dv13, xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF,info.subpel[t_level],0);

		get_aff_bit_cost = info.lambda[t_level] * ( (int)(pos2/4) + AFF_IDX_OFFSET + DIRECT_LEN + MERGE_LEN ) + get_bit_cost(info.lambda[t_level],dv21->mvx,dv21->mvy,
				get2_v1[pred21].mvx,get2_v1[pred21].mvy,ctx2x,ctx2y,info.subpel[t_level]) + get_bit_cost(info.lambda[t_level],dv22->mvx,dv22->mvy,get2_v2[pred22].mvx,
				get2_v2[pred22].mvy,ctx2x,ctx2y,info.subpel[t_level]) + get_bit_cost(info.lambda[t_level],dv23->mvx,dv23->mvy,get2_v3[pred23].mvx,
				get2_v3[pred23].mvy,ctx2x,ctx2y,info.subpel[t_level]) + info.lambda[t_level] * get_mode_coding_cost(RIGHT_CONNECTED_AFF, 
				fmv2,info.bi_mv[t_level], t_level);
	
		get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

		aff_sad_cache = get_aff_sad_cost;

		if( ( (get_aff_total_cost < (fmv2->total_cost + 1 * info.lambda[t_level]) ) && fmv2->bi_mode <= 8 &&
			get_aff_sad_cost < (fmv2->sad_cost * aff_coef) ) ||
			(fmv2->bi_mode >= 9 && (get_aff_total_cost < fmv2->total_cost ) && get_aff_sad_cost < fmv2->sad_cost) ){  //affine decision

			aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var, *dv21, *dv22, *dv23, 
				*dv11, *dv12,*dv13, xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF);

//			printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

			if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
				aff_mse < IBLOCK_FACTOR * aff_ref_var )
			  || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
				assert(fmv2->total_cost == fmv1->total_cost);

				fmv2->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;

				fmv2->sad_cost = get_aff_sad_cost;
				fmv2->mse = aff_mse;
				fmv2->bit_cost = get_aff_bit_cost;
				fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
				fmv2->aff_idx = pos2;
				fmv2->direct_idx = INDIRECT;
				fmv2->merge_idx = INTER;
				fmv2->is_predictor = YES;
				fmv2->lifting_mode = CONNECTED;
				fmv2->aff1_mvx = dv21->mvx;		fmv2->aff1_mvy = dv21->mvy;
				fmv2->aff2_mvx = dv22->mvx;		fmv2->aff2_mvy = dv22->mvy;
				fmv2->aff3_mvx = dv23->mvx;		fmv2->aff3_mvy = dv23->mvy;
			
				fmv2->aff1_dmvx = dv21->mvx - get2_v1[pred21].mvx;
				fmv2->aff1_dmvy = dv21->mvy - get2_v1[pred21].mvy;
				fmv2->aff2_dmvx = dv22->mvx - get2_v2[pred22].mvx;
				fmv2->aff2_dmvy = dv22->mvy - get2_v2[pred22].mvy;
				fmv2->aff3_dmvx = dv23->mvx - get2_v3[pred23].mvx;
				fmv2->aff3_dmvy = dv23->mvy - get2_v3[pred23].mvy;

//				printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n",fmv1->aff2_dmvx,fmv1->aff2_dmvy);

				fmv1->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;

				fmv1->sad_cost = fmv2->sad_cost;
				fmv1->mse = aff_mse;
				fmv1->bit_cost = fmv2->bit_cost;
				fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
				fmv1->aff_idx = -1;
				fmv1->direct_idx = INDIRECT;
				fmv1->merge_idx = INTER;
				fmv1->is_predictor = NO;
				fmv1->lifting_mode = IGNORED;
				fmv1->aff1_mvx = (float)HUGE_VAL;fmv1->aff1_mvy = (float)HUGE_VAL;
				fmv1->aff2_mvx = (float)HUGE_VAL;fmv1->aff2_mvy = (float)HUGE_VAL;
				fmv1->aff3_mvx = (float)HUGE_VAL;fmv1->aff3_mvy = (float)HUGE_VAL;
			}//If MSE
		}
	}
//INTER

//AFF_MERGE MODE
//UP MERGE
		//	printf("BLOCK B\n");
		find_block(x,y-1,fmv2_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,1);//BLOCK B	
		merg_aff_v1.mvx = *get1x;
		merg_aff_v1.mvy = *get1y;
		//	printf("BLOCK D\n");
		find_block(x+xblk2-1,y-1,fmv2_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);//BLOCK D
		merg_aff_v2.mvx = *get1x;
		merg_aff_v2.mvy = *get1y;

		if(merg_aff_v1.mvx != (float)HUGE_VAL && merg_aff_v1.mvy != (float)HUGE_VAL && 
			merg_aff_v2.mvx != (float)HUGE_VAL && merg_aff_v2.mvy != (float)HUGE_VAL){

			mrg_cnt1 = 0;

			for(i = 0;i < 2; i++){		
				if(get2_v3[i].mvx != (float)HUGE_VAL && get2_v3[i].mvy != (float)HUGE_VAL){
					best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx2x, ctx2y, &merg_aff_v1, &merg_aff_v2, &get2_v3[i], dv21,dv22,dv23,xblk2,
					yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,3);

					get_aff_sad_cost = best_sad;
//					find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2,merg_aff_v1, merg_aff_v2, *dv23, merg_aff_v1, merg_aff_v2,*dv23, xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF,info.subpel[t_level],0);
					
					get_aff_bit_cost = info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_UP_LEN + 1 ) + info.lambda[t_level] * get_mode_coding_cost(RIGHT_CONNECTED_AFF,fmv2,info.bi_mv[t_level],
					t_level) + get_bit_cost(info.lambda[t_level],dv23->mvx,dv23->mvy,get2_v3[i].mvx,get2_v3[i].mvy,ctx2x,ctx2y,info.subpel[t_level]);

					get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

					assert(fmv1->total_cost == fmv2->total_cost);

					if( ( (get_aff_total_cost < (fmv2->total_cost + 1 * info.lambda[t_level]) ) && fmv2->bi_mode <= 8 &&
					get_aff_sad_cost < (fmv2->sad_cost * aff_coef) ) ||
					(fmv2->bi_mode >= 9 && (get_aff_total_cost < fmv2->total_cost ) && get_aff_sad_cost < fmv2->sad_cost) ){  //affine decision


						aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var, merg_aff_v1, merg_aff_v2, *dv23,
							merg_aff_v1, merg_aff_v2,*dv23, xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF);

//						printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

						if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
							aff_mse < IBLOCK_FACTOR * aff_ref_var )
						    || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
							fmv2->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;

							fmv2->sad_cost = get_aff_sad_cost;
							fmv2->mse = aff_mse;
							fmv2->bit_cost = get_aff_bit_cost;
							fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
							fmv2->aff_idx = mrg_cnt1;
							fmv2->direct_idx = INDIRECT;
							fmv2->merge_idx = MERGE;
							fmv2->merge_dir = UP;
							fmv2->is_predictor = YES;
							fmv2->lifting_mode = CONNECTED;
							fmv2->aff1_mvx = merg_aff_v1.mvx;fmv2->aff1_mvy = merg_aff_v1.mvy;
							fmv2->aff2_mvx = merg_aff_v2.mvx;fmv2->aff2_mvy = merg_aff_v2.mvy;
							fmv2->aff3_mvx = dv23->mvx;		 fmv2->aff3_mvy = dv23->mvy;
	
							fmv2->aff3_dmvx = dv23->mvx - get2_v3[i].mvx;
							fmv2->aff3_dmvy = dv23->mvy - get2_v3[i].mvy;

	//						printf("fmv2->aff3_dmvx = %f, fmv2->aff3_dmvy = %f\n",fmv2->aff3_dmvx,fmv2->aff3_dmvy);

							fmv1->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;
		
							fmv1->sad_cost = fmv2->sad_cost;
							fmv1->mse = aff_mse;
							fmv1->bit_cost = fmv2->bit_cost;
							fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
							fmv1->aff_idx = -1;
							fmv1->direct_idx = INDIRECT;
							fmv1->merge_idx = MERGE;
							fmv1->merge_dir = UP;
							fmv1->is_predictor = NO;
							fmv1->lifting_mode = IGNORED;
							fmv1->aff1_mvx = (float)HUGE_VAL;fmv1->aff1_mvy = (float)HUGE_VAL;
							fmv1->aff2_mvx = (float)HUGE_VAL;fmv1->aff2_mvy = (float)HUGE_VAL;
							fmv1->aff3_mvx = (float)HUGE_VAL;fmv1->aff3_mvy = (float)HUGE_VAL;
						}//If MSE
					}

					mrg_cnt1 ++;
				}
			}
		}

//LEFT MERGE
		//	printf("BLOCK C\n");
		find_block(x-1,y,fmv2_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,0);//BLOCK C	
		merg_aff_v1.mvx = *get1x;
		merg_aff_v1.mvy = *get1y;
		//	printf("BLOCK F\n");
		find_block(x-1,y+yblk2-1,fmv2_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);//BLOCK F
		merg_aff_v2.mvx = *get1x;
		merg_aff_v2.mvy = *get1y;

		if(merg_aff_v1.mvx != (float)HUGE_VAL && merg_aff_v1.mvy != (float)HUGE_VAL && 
			merg_aff_v2.mvx != (float)HUGE_VAL && merg_aff_v2.mvy != (float)HUGE_VAL){

			mrg_cnt1 = 0;

			for(i = 0;i < 2; i++){		
				if(get2_v2[i].mvx != (float)HUGE_VAL && get2_v2[i].mvy != (float)HUGE_VAL){
					best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx2x, ctx2y, &merg_aff_v1, &get2_v2[i], &merg_aff_v2, dv21,dv22,dv23,xblk2,yblk2,
					x,y,hor,ver,RIGHT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,2);

					get_aff_sad_cost = best_sad;
//					find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, merg_aff_v1, *dv22, merg_aff_v2, merg_aff_v1, *dv22, merg_aff_v2, xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF,info.subpel[t_level],0);
					
					get_aff_bit_cost = info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + info.lambda[t_level] * get_mode_coding_cost(RIGHT_CONNECTED_AFF, 
						fmv2,info.bi_mv[t_level], t_level) + get_bit_cost(info.lambda[t_level],dv22->mvx,dv22->mvy,get2_v2[i].mvx,get2_v2[i].mvy,
						ctx2x,ctx2y,info.subpel[t_level]);

					assert(fmv1->total_cost == fmv2->total_cost);

					get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

					if( ( (get_aff_total_cost < (fmv2->total_cost + 1 * info.lambda[t_level]) ) && fmv2->bi_mode <= 8 &&
						get_aff_sad_cost < (fmv2->sad_cost * aff_coef) ) ||
						(fmv2->bi_mode >= 9 && (get_aff_total_cost < fmv2->total_cost ) && get_aff_sad_cost < fmv2->sad_cost) ){  //affine decision

						aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var, merg_aff_v1, *dv22, merg_aff_v2, 
							merg_aff_v1, *dv22, merg_aff_v2, xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF);

//						printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

						if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
							aff_mse < IBLOCK_FACTOR * aff_ref_var )
						    || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
							fmv2->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;

							fmv2->sad_cost = get_aff_sad_cost;
							fmv2->mse = aff_mse;
							fmv2->bit_cost = get_aff_bit_cost;
							fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
							fmv2->aff_idx = mrg_cnt1;
							fmv2->direct_idx = INDIRECT;
							fmv2->merge_idx = MERGE;
							fmv2->merge_dir = LEFT;
							fmv2->is_predictor = YES;
							fmv2->lifting_mode = CONNECTED;
							fmv2->aff1_mvx = merg_aff_v1.mvx;fmv2->aff1_mvy = merg_aff_v1.mvy;
							fmv2->aff2_mvx = dv22->mvx;		 fmv2->aff2_mvy = dv22->mvy;
							fmv2->aff3_mvx = merg_aff_v2.mvx;fmv2->aff3_mvy = merg_aff_v2.mvy;
	
							fmv2->aff2_dmvx = dv22->mvx - get2_v2[i].mvx;
							fmv2->aff2_dmvy = dv22->mvy - get2_v2[i].mvy;

	//						printf("fmv2->aff2_dmvx = %f, fmv2->aff2_dmvy = %f\n",fmv2->aff2_dmvx,fmv2->aff2_dmvy);

							fmv1->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;
		
							fmv1->sad_cost = fmv2->sad_cost;
							fmv1->mse = aff_mse;
							fmv1->bit_cost = fmv2->bit_cost;
							fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
							fmv1->aff_idx = -1;
							fmv1->direct_idx = INDIRECT;
							fmv1->merge_idx = MERGE;
							fmv1->merge_dir = LEFT;
							fmv1->is_predictor = NO;
							fmv1->lifting_mode = IGNORED;
							fmv1->aff1_mvx = (float)HUGE_VAL;fmv1->aff1_mvy = (float)HUGE_VAL;
							fmv1->aff2_mvx = (float)HUGE_VAL;fmv1->aff2_mvy = (float)HUGE_VAL;
							fmv1->aff3_mvx = (float)HUGE_VAL;fmv1->aff3_mvy = (float)HUGE_VAL;
						}//If MSE
					}

					mrg_cnt1 ++;
				}
			}
		}
//MERGE mode

/* TEST	*/
	if(pred21 >= 0 && pred22 >= 0 && pred23 >= 0 && do_inter == 1){
		if( (best_trans_mode == RIGHT_CONNECTED || best_trans_mode == RIGHT_PREDICTED ||
			(best_trans_mode == BLOCK_MERGING && trans_mvr.mvx != (float)HUGE_VAL && trans_mvl.mvx == (float)HUGE_VAL) ) && fmv2 != NULL){
			if(best_trans_sad < best_aff_sad_cost){
//				printf("right trans better, bi_mode = %d, x = %d, y = %d, xblk = %d, yblk = %d\naff_sad = %f, trans_sad = %f\n\n",
//					best_trans_mode,x,y,xblk2,yblk2,best_aff_sad_cost,best_trans_sad);

				best_sad = delta_v_search(fr_cur, fr_ref1, fr_ref2, upframe1,upframe2,ctx2x, ctx2y, &trans_mvr, &trans_mvr, &trans_mvr, dv21,dv22,dv23,xblk2,yblk2,x,y,
					hor,ver,RIGHT_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);

				get_aff_sad_cost = best_sad;

				get_aff_bit_cost = get_bit_cost(info.lambda[t_level],dv21->mvx,dv21->mvy,trans_mvr.mvx,trans_mvr.mvy,ctx2x,ctx2y,info.subpel[t_level]) +
					get_bit_cost(info.lambda[t_level],dv22->mvx,dv22->mvy,trans_mvr.mvx,trans_mvr.mvy,ctx2x,ctx2y,info.subpel[t_level]) +
					get_bit_cost(info.lambda[t_level],dv23->mvx,dv23->mvy,trans_mvr.mvx,trans_mvr.mvy,ctx2x,ctx2y,info.subpel[t_level]);

				dv21->dmvx = dv21->mvx - trans_mvr.mvx;	dv21->dmvy = dv21->mvy - trans_mvr.mvy;
				dv22->dmvx = dv22->mvx - trans_mvr.mvx;	dv22->dmvy = dv22->mvy - trans_mvr.mvy;
				dv23->dmvx = dv23->mvx - trans_mvr.mvx;	dv23->dmvy = dv23->mvy - trans_mvr.mvy;

//Apply trans-based aff result if better
				if(best_trans_mode == BLOCK_MERGING){

					get_aff_bit_cost = get_aff_bit_cost + info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + 
						info.lambda[t_level] * get_mode_coding_cost(RIGHT_CONNECTED_AFF, fmv2,info.bi_mv[t_level], t_level);

				}else if(best_trans_mode == RIGHT_CONNECTED || best_trans_mode == RIGHT_PREDICTED){
					get_aff_bit_cost = get_aff_bit_cost + get_bit_cost(info.lambda[t_level],right_dmvx,right_dmvy,0,0,ctx2x,ctx2y,info.subpel[t_level]);

					get_aff_bit_cost = get_aff_bit_cost + info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + 
						info.lambda[t_level] * get_mode_coding_cost(RIGHT_CONNECTED_AFF, fmv2,info.bi_mv[t_level], t_level);
					
				}else
					assert(0);

				get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;	

//decision making
				if( ( (get_aff_total_cost < (fmv2->total_cost + 1 * info.lambda[t_level]) ) && fmv2->bi_mode <= 8 &&
					get_aff_sad_cost < (fmv2->sad_cost * aff_coef * t_coef) ) ||
					(fmv2->bi_mode >= 9 && (get_aff_total_cost < fmv2->total_cost ) && get_aff_sad_cost < fmv2->sad_cost) ){
						
					aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*dv21, *dv22, *dv23, *dv11, *dv12,
						*dv13, xblk2,yblk2,x,y,hor,ver,RIGHT_CONNECTED_AFF);

					if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
						aff_mse < IBLOCK_FACTOR * aff_ref_var )
					   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

						fmv2->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;

						fmv2->sad_cost = get_aff_sad_cost;
						fmv2->mse = aff_mse;
						fmv2->bit_cost = get_aff_bit_cost;
						fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;

						fmv2->aff_idx = -1;

						fmv2->direct_idx = INDIRECT;
						fmv2->merge_idx = MERGE;
						fmv2->merge_dir = TRAN_P;
						
						if(best_trans_mode == BLOCK_MERGING)
							fmv2->trans_pred_idx = DIR;
						else{
							assert(best_trans_mode == RIGHT_CONNECTED);
							fmv2->trans_pred_idx = INDIR;
						}

						fmv2->is_predictor = YES;
						fmv2->lifting_mode = CONNECTED;

						fmv2->dmvx = right_dmvx;		fmv2->dmvy = right_dmvy;
						fmv2->aff1_mvx = dv21->mvx;		fmv2->aff1_mvy = dv21->mvy;
						fmv2->aff2_mvx = dv22->mvx;		fmv2->aff2_mvy = dv22->mvy;
						fmv2->aff3_mvx = dv23->mvx;		fmv2->aff3_mvy = dv23->mvy;
			
						fmv2->aff1_dmvx = dv21->dmvx;	fmv2->aff1_dmvy = dv21->dmvy;						
						fmv2->aff2_dmvx = dv22->dmvx;	fmv2->aff2_dmvy = dv22->dmvy;
						fmv2->aff3_dmvx = dv23->dmvx;	fmv2->aff3_dmvy = dv23->dmvy;


						fmv1->bi_mode = (BiMode)RIGHT_CONNECTED_AFF;

						fmv1->sad_cost = fmv2->sad_cost;
						fmv1->mse = aff_mse;
						fmv1->bit_cost = fmv2->bit_cost;
						fmv1->total_cost = fmv2->sad_cost + fmv2->bit_cost;
						fmv1->aff_idx = -1;

						fmv1->direct_idx = INDIRECT;
						fmv1->merge_idx = MERGE;
						fmv1->merge_dir = TRAN_P;

						if(best_trans_mode == BLOCK_MERGING)
							fmv1->trans_pred_idx = DIR;
						else{
							assert(best_trans_mode == RIGHT_CONNECTED);
							fmv1->trans_pred_idx = INDIR;
						}

						fmv1->is_predictor = NO;
						fmv1->lifting_mode = IGNORED;

						fmv1->dmvx = (float)HUGE_VAL;		fmv1->dmvy = (float)HUGE_VAL;
						fmv1->aff1_mvx = (float)HUGE_VAL;	fmv1->aff1_mvy = (float)HUGE_VAL;
						fmv1->aff2_mvx = (float)HUGE_VAL;	fmv1->aff2_mvy = (float)HUGE_VAL;
						fmv1->aff3_mvx = (float)HUGE_VAL;	fmv1->aff3_mvy = (float)HUGE_VAL;
					}//If MSE
				}
			}
		}
	}//if pred
/*	END OF TEST	*/


//BI-DIRECTIONAL SEARCH
//AFF_DIRECT MODE
		if(pos1 >= 0 && pos2 >= 0){
			best_aff_bit_cost = (float)HUGE_VAL;
			best_aff_sad_cost = (float)HUGE_VAL;
			best_aff_total_cost = best_aff_sad_cost + best_aff_bit_cost;

			for(round = 0;round < 3; round++){
				assert(get_v1[pred11].mvx != (float)HUGE_VAL && get_v1[pred11].mvy != (float)HUGE_VAL && get_v2[pred12].mvx != (float)HUGE_VAL && 
					get_v2[pred12].mvy != (float)HUGE_VAL && get_v3[pred13].mvx != (float)HUGE_VAL && get_v3[pred13].mvy != (float)HUGE_VAL);

				count = 0;

				//RIGHT SEARCH
				for(i=0;i<3;i++){
					for(j=0;j<3;j++){
						for(k=0;k<3;k++){
							if( get2_v1[i].mvx!=(float)HUGE_VAL && get2_v2[j].mvx!=(float)HUGE_VAL && get2_v3[k].mvx!=(float)HUGE_VAL &&
								get2_v1[i].mvy!=(float)HUGE_VAL && get2_v2[j].mvy!=(float)HUGE_VAL && get2_v3[k].mvy!=(float)HUGE_VAL ){
								get_aff_sad_cost = find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2,get_v1[pred11],get_v2[pred12],get_v3[pred13],get2_v1[i],
									get2_v2[j],get2_v3[k],xblk2,yblk2,x,y,hor,ver, BI_CONNECTED_AFF,info.subpel[t_level],0);
								get_aff_bit_cost = info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF, fmv2, info.bi_mv[t_level], t_level)
									+ info.lambda[t_level] * ((int)(aff_num1/4) + AFF_IDX_OFFSET) + info.lambda[t_level] * ((int)(aff_num2/4) + AFF_IDX_OFFSET) + info.lambda[t_level]*1;
								get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;
								if( (get_aff_total_cost < best_aff_total_cost) ){
									assert(get_aff_total_cost != (float)HUGE_VAL);

									best_aff_bit_cost = get_aff_bit_cost;
									best_aff_sad_cost = get_aff_sad_cost;
									best_aff_total_cost = get_aff_total_cost;
									best_idx2 = i * 9 + j * 3 + k;
									pred21 = i;
									pred22 = j;
									pred23 = k;
									pos2 = count;
//									printf("count = %d, length = %d\n",count, (int)(count/4) );
								}
								count++;
							}
						}
					}
				}

				assert(get2_v1[pred21].mvx != (float)HUGE_VAL && get2_v1[pred21].mvy != (float)HUGE_VAL && get2_v2[pred22].mvx != (float)HUGE_VAL && 
					get2_v2[pred22].mvy != (float)HUGE_VAL && get2_v3[pred23].mvx != (float)HUGE_VAL && get2_v3[pred23].mvy != (float)HUGE_VAL);

				count = 0;

				//LEFT SEARCH
				for(i=0;i<3;i++){
					for(j=0;j<3;j++){
						for(k=0;k<3;k++){
							if( get_v1[i].mvx!=(float)HUGE_VAL && get_v2[j].mvx!=(float)HUGE_VAL && get_v3[k].mvx!=(float)HUGE_VAL &&
								get_v1[i].mvy!=(float)HUGE_VAL && get_v2[j].mvy!=(float)HUGE_VAL && get_v3[k].mvy!=(float)HUGE_VAL ){
								get_aff_sad_cost = find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2,get_v1[i],get_v2[j],get_v3[k],get2_v1[pred21],
									get2_v2[pred22],get2_v3[pred23],xblk2,yblk2,x,y,hor,ver, BI_CONNECTED_AFF,info.subpel[t_level],0);

								get_aff_bit_cost = info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF, fmv2, info.bi_mv[t_level], t_level)
								+ info.lambda[t_level] * ((int)(aff_num1/4) + AFF_IDX_OFFSET) + info.lambda[t_level] * ((int)(aff_num2/4) + AFF_IDX_OFFSET) + info.lambda[t_level]*1;
								get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;
								if( (get_aff_total_cost < best_aff_total_cost) ){
									assert(get_aff_total_cost != (float)HUGE_VAL);

									best_aff_bit_cost = get_aff_bit_cost;
									best_aff_sad_cost = get_aff_sad_cost;
									best_aff_total_cost = get_aff_total_cost;
									best_idx1 = i * 9 + j * 3 + k;
									pred11 = i;
									pred12 = j;
									pred13 = k;
									pos1 = count;
//									printf("count = %d, length = %d\n",count, (int)(count/4) );
								}
								count++;
							}
						}
					}
				}
			}//round

			if( ( (best_aff_total_cost < (fmv1->total_cost + 3 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
				best_aff_sad_cost < (fmv1->sad_cost - 4 * info.lambda[t_level]) ) ||
				(fmv1->bi_mode >= 9 && (best_aff_total_cost < fmv1->total_cost ) && best_aff_sad_cost < fmv1->sad_cost) ){  //affine decision
				
				aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var, get_v1[pred11],get_v2[pred12],get_v3[pred13],get2_v1[pred21],
									get2_v2[pred22],get2_v3[pred23], xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF);

//				printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

				if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
					aff_mse < IBLOCK_FACTOR * aff_ref_var )
				    || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){
					assert(fmv1->total_cost == fmv2->total_cost);
					assert(pos1 >= 0 && pos1 <= 26);
					assert(pos2 >= 0 && pos2 <= 26);

					assert(best_idx2 >= 0 && best_idx2 <= 26);
					fmv2->bi_mode = (BiMode)BI_CONNECTED_AFF;

					fmv2->sad_cost = best_aff_sad_cost;
					fmv2->mse = aff_mse;
					fmv2->bit_cost = best_aff_bit_cost;
					fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
					fmv2->aff_idx = pos2;
					fmv2->direct_idx = DIRECT;
					fmv2->is_predictor = YES;
					fmv2->lifting_mode = CONNECTED;
					fmv2->aff1_mvx = get2_v1[pred21].mvx;fmv2->aff1_mvy = get2_v1[pred21].mvy;
					fmv2->aff2_mvx = get2_v2[pred22].mvx;fmv2->aff2_mvy = get2_v2[pred22].mvy;
					fmv2->aff3_mvx = get2_v3[pred23].mvx;fmv2->aff3_mvy = get2_v3[pred23].mvy;

					fmv1->bi_mode = (BiMode)BI_CONNECTED_AFF;
		
					fmv1->sad_cost = fmv2->sad_cost;
					fmv1->mse = aff_mse;
					fmv1->bit_cost = fmv2->bit_cost;
					fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
					fmv1->aff_idx = pos1;
					fmv1->direct_idx = DIRECT;
					fmv1->is_predictor = YES;
					fmv1->lifting_mode = CONNECTED;
					fmv1->aff1_mvx = get_v1[pred11].mvx;fmv1->aff1_mvy = get_v1[pred11].mvy;
					fmv1->aff2_mvx = get_v2[pred12].mvx;fmv1->aff2_mvy = get_v2[pred12].mvy;
					fmv1->aff3_mvx = get_v3[pred13].mvx;fmv1->aff3_mvy = get_v3[pred13].mvy;
				}//If MSE
			}

//BI-DIRECTIONAL INTER
		  if(do_inter == 1){
			best_sad = delta_v_search_bi(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx1x, ctx1y, ctx2x, ctx2y, &get_v1[pred11], &get_v2[pred12], &get_v3[pred13], 
				&get2_v1[pred21], &get2_v2[pred22], &get2_v3[pred23],dv11,dv12,dv13,dv21,dv22,dv23,xblk2,yblk2,x,y,
				hor,ver,BI_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);

			get_aff_sad_cost = best_sad;
			
//			find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, *dv11, *dv12, *dv13, *dv21, *dv22,*dv23, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF,info.subpel[t_level],0);

			get_aff_bit_cost = info.lambda[t_level] * ( (int)(pos1/4) + AFF_IDX_OFFSET ) + info.lambda[t_level] * ( (int)(pos2/4) + AFF_IDX_OFFSET )+
					get_bit_cost(info.lambda[t_level],dv11->mvx,dv11->mvy,get_v1[pred11].mvx,get_v1[pred11].mvy,ctx1x,ctx1y,info.subpel[t_level]) + 
					get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,get_v2[pred12].mvx,get_v2[pred12].mvy,ctx1x,ctx1y,info.subpel[t_level])
					+ get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,get_v3[pred13].mvx,get_v3[pred13].mvy,ctx1x,ctx1y,info.subpel[t_level])
					+ get_bit_cost(info.lambda[t_level],dv21->mvx,dv21->mvy,get2_v1[pred21].mvx,get2_v1[pred21].mvy,ctx2x,ctx2y,info.subpel[t_level]) 
					+ get_bit_cost(info.lambda[t_level],dv22->mvx,dv22->mvy,get2_v2[pred22].mvx,get2_v2[pred22].mvy,ctx2x,ctx2y,info.subpel[t_level])
					+ get_bit_cost(info.lambda[t_level],dv23->mvx,dv23->mvy,get2_v3[pred23].mvx,get2_v3[pred23].mvy,ctx2x,ctx2y,info.subpel[t_level])
					+ info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF,fmv2,info.bi_mv[t_level], t_level) + (DIRECT_LEN + MERGE_LEN) * info.lambda[t_level];

			get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

			aff_sad_cache = get_aff_sad_cost;

			if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
				get_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
				(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

				aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*dv11, *dv12, *dv13, 
							*dv21, *dv22, *dv23, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF);

//				printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

				if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
					aff_mse < IBLOCK_FACTOR * aff_ref_var )
				    || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

					assert(fmv2->total_cost == fmv1->total_cost);

					fmv2->bi_mode = (BiMode)BI_CONNECTED_AFF;

					fmv2->sad_cost = get_aff_sad_cost;
					fmv2->mse = aff_mse;
					fmv2->bit_cost = get_aff_bit_cost;
					fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
					fmv2->aff_idx = pos2;
					fmv2->direct_idx = INDIRECT;
					fmv2->merge_idx = INTER;
					fmv2->is_predictor = YES;
					fmv2->lifting_mode = CONNECTED;
					fmv2->aff1_mvx = dv21->mvx;		fmv2->aff1_mvy = dv21->mvy;
					fmv2->aff2_mvx = dv22->mvx;		fmv2->aff2_mvy = dv22->mvy;
					fmv2->aff3_mvx = dv23->mvx;		fmv2->aff3_mvy = dv23->mvy;
			
					fmv2->aff1_dmvx = dv21->mvx - get2_v1[pred21].mvx;
					fmv2->aff1_dmvy = dv21->mvy - get2_v1[pred21].mvy;
					fmv2->aff2_dmvx = dv22->mvx - get2_v2[pred22].mvx;
					fmv2->aff2_dmvy = dv22->mvy - get2_v2[pred22].mvy;
					fmv2->aff3_dmvx = dv23->mvx - get2_v3[pred23].mvx;
					fmv2->aff3_dmvy = dv23->mvy - get2_v3[pred23].mvy;

	//				printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n",fmv1->aff2_dmvx,fmv1->aff2_dmvy);

					fmv1->bi_mode = (BiMode)BI_CONNECTED_AFF;

					fmv1->sad_cost = fmv2->sad_cost;
					fmv1->mse = aff_mse;
					fmv1->bit_cost = fmv2->bit_cost;
					fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
					fmv1->aff_idx = pos1;
					fmv1->direct_idx = INDIRECT;
					fmv1->merge_idx = INTER;
					fmv1->is_predictor = YES;
					fmv1->lifting_mode = CONNECTED;
					fmv1->aff1_mvx = dv11->mvx;		fmv1->aff1_mvy = dv11->mvy;
					fmv1->aff2_mvx = dv12->mvx;		fmv1->aff2_mvy = dv12->mvy;
					fmv1->aff3_mvx = dv13->mvx;		fmv1->aff3_mvy = dv13->mvy;
			
					fmv1->aff1_dmvx = dv11->mvx - get_v1[pred11].mvx;
					fmv1->aff1_dmvy = dv11->mvy - get_v1[pred11].mvy;
					fmv1->aff2_dmvx = dv12->mvx - get_v2[pred12].mvx;
					fmv1->aff2_dmvy = dv12->mvy - get_v2[pred12].mvy;
					fmv1->aff3_dmvx = dv13->mvx - get_v3[pred13].mvx;
					fmv1->aff3_dmvy = dv13->mvy - get_v3[pred13].mvy;
				}//If MSE
			}
		  }
//INTER

/* TEST	*/
		if(fmv2 != NULL){
		  if( (best_trans_mode == BI_CONNECTED ||  
//			  best_trans_mode == PARALLEL || 
			(best_trans_mode == BLOCK_MERGING && trans_mvl.mvx != (float)HUGE_VAL && trans_mvr.mvx != (float)HUGE_VAL ) ) && do_inter == 1 ){
			  if(best_trans_sad < best_aff_sad_cost){
//			    printf("bi trans better, bi_mode = %d, x = %d, y = %d, xblk = %d, yblk = %d\naff_sad = %f, trans_sad = %f\n\n",
//			    		best_trans_mode,x,y,xblk2,yblk2,best_aff_sad_cost,best_trans_sad);
			    
				if(best_trans_mode == BI_CONNECTED || best_trans_mode == BLOCK_MERGING){
					best_sad = delta_v_search_bi(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx1x, ctx1y, ctx2x, ctx2y, &trans_mvl, &trans_mvl, 
						&trans_mvl, &trans_mvr, &trans_mvr, &trans_mvr, dv11,dv12,dv13,dv21,dv22,dv23,xblk2,yblk2,x,y,
						hor,ver,BI_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);
				}else{
					assert(best_trans_mode == PARALLEL);

					best_sad = delta_v_search_bi(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx1x, ctx1y, ctx2x, ctx2y, &trans_mvl, &trans_mvl, &trans_mvl, 
						&trans_mvr, &trans_mvr, &trans_mvr, dv11,dv12,dv13,dv21,dv22,dv23,xblk2,yblk2,x,y,
						hor,ver,BI_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);
				}

				get_aff_sad_cost = best_sad;

				if(1){
					get_aff_bit_cost = get_bit_cost(info.lambda[t_level],dv11->mvx,dv11->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]) +
						get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]) +
						get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]) + 
						get_bit_cost(info.lambda[t_level],dv21->mvx,dv21->mvy,trans_mvr.mvx,trans_mvr.mvy,ctx2x,ctx2y,info.subpel[t_level]) +
						get_bit_cost(info.lambda[t_level],dv22->mvx,dv22->mvy,trans_mvr.mvx,trans_mvr.mvy,ctx2x,ctx2y,info.subpel[t_level]) +
						get_bit_cost(info.lambda[t_level],dv23->mvx,dv23->mvy,trans_mvr.mvx,trans_mvr.mvy,ctx2x,ctx2y,info.subpel[t_level]);
				}else{
					assert(0);
					get_aff_bit_cost = get_bit_cost(info.lambda[t_level],dv11->mvx,dv11->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]) +
						get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]) +
						get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,trans_mvl.mvx,trans_mvl.mvy,ctx1x,ctx1y,info.subpel[t_level]);
				}

				dv11->dmvx = dv11->mvx - trans_mvl.mvx;	dv11->dmvy = dv11->mvy - trans_mvl.mvy;
				dv12->dmvx = dv12->mvx - trans_mvl.mvx;	dv12->dmvy = dv12->mvy - trans_mvl.mvy;
				dv13->dmvx = dv13->mvx - trans_mvl.mvx;	dv13->dmvy = dv13->mvy - trans_mvl.mvy;

				dv21->dmvx = dv21->mvx - trans_mvr.mvx;	dv21->dmvy = dv21->mvy - trans_mvr.mvy;
				dv22->dmvx = dv22->mvx - trans_mvr.mvx;	dv22->dmvy = dv22->mvy - trans_mvr.mvy;
				dv23->dmvx = dv23->mvx - trans_mvr.mvx;	dv23->dmvy = dv23->mvy - trans_mvr.mvy;

				if(best_trans_mode == BLOCK_MERGING){

					get_aff_bit_cost = get_aff_bit_cost + info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + 
						info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF, fmv2,info.bi_mv[t_level], t_level);

				}else if(best_trans_mode == BI_CONNECTED){
					get_aff_bit_cost = get_aff_bit_cost + get_bit_cost(info.lambda[t_level],right_dmvx,right_dmvy,0,0,ctx2x,ctx2y,info.subpel[t_level])
						 + get_bit_cost(info.lambda[t_level],left_dmvx,left_dmvy,0,0,ctx1x,ctx1y,info.subpel[t_level]);

					get_aff_bit_cost = get_aff_bit_cost + info.lambda[t_level] * ( DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 ) + 
						info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF, fmv2,info.bi_mv[t_level], t_level);
					
				}else
					assert(0);

				get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;	

//decision making
				if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
					get_aff_sad_cost < (fmv1->sad_cost * aff_coef * t_coef) ) ||
					(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){
						
					aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*dv11, *dv12, *dv13, *dv21, *dv22,
						*dv23, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF);

					if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
						aff_mse < IBLOCK_FACTOR * aff_ref_var )
					   || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

						fmv2->bi_mode = (BiMode)BI_CONNECTED_AFF;

						fmv2->sad_cost = get_aff_sad_cost;
						fmv2->mse = aff_mse;
						fmv2->bit_cost = get_aff_bit_cost;
						fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
						fmv2->aff_idx = -1;

						fmv2->direct_idx = INDIRECT;
						fmv2->merge_idx = MERGE;
						fmv2->merge_dir = TRAN_P;
						
						if(best_trans_mode == BLOCK_MERGING)
							fmv2->trans_pred_idx = DIR;
						else{
							assert(best_trans_mode == BI_CONNECTED || best_trans_mode == PARALLEL);
							fmv2->trans_pred_idx = INDIR;
						}

						fmv2->is_predictor = YES;
						fmv2->lifting_mode = CONNECTED;

						fmv2->dmvx = right_dmvx;		fmv2->dmvy = right_dmvy;
						fmv2->aff1_mvx = dv21->mvx;		fmv2->aff1_mvy = dv21->mvy;
						fmv2->aff2_mvx = dv22->mvx;		fmv2->aff2_mvy = dv22->mvy;
						fmv2->aff3_mvx = dv23->mvx;		fmv2->aff3_mvy = dv23->mvy;
			
						fmv2->aff1_dmvx = dv21->dmvx;	fmv2->aff1_dmvy = dv21->dmvy;						
						fmv2->aff2_dmvx = dv22->dmvx;	fmv2->aff2_dmvy = dv22->dmvy;
						fmv2->aff3_dmvx = dv23->dmvx;	fmv2->aff3_dmvy = dv23->dmvy;

						
	
						fmv1->bi_mode = (BiMode)BI_CONNECTED_AFF;

						fmv1->sad_cost = fmv2->sad_cost;
						fmv1->mse = aff_mse;
						fmv1->bit_cost = fmv2->bit_cost;
						fmv1->total_cost = fmv2->sad_cost + fmv2->bit_cost;
						fmv1->aff_idx = -1;

						fmv1->direct_idx = INDIRECT;
						fmv1->merge_idx = MERGE;
						fmv1->merge_dir = TRAN_P;

						if(best_trans_mode == BLOCK_MERGING)
							fmv1->trans_pred_idx = DIR;
						else{
							assert(best_trans_mode == BI_CONNECTED || best_trans_mode == PARALLEL);
							fmv1->trans_pred_idx = INDIR;
						}

						fmv1->is_predictor = YES;
						fmv1->lifting_mode = CONNECTED;

						fmv1->dmvx = left_dmvx;			fmv1->dmvy = left_dmvy;
						fmv1->aff1_mvx = dv11->mvx;		fmv1->aff1_mvy = dv11->mvy;
						fmv1->aff2_mvx = dv12->mvx;		fmv1->aff2_mvy = dv12->mvy;
						fmv1->aff3_mvx = dv13->mvx;		fmv1->aff3_mvy = dv13->mvy;
			
						fmv1->aff1_dmvx = dv11->dmvx;	fmv1->aff1_dmvy = dv11->dmvy;						
						fmv1->aff2_dmvx = dv12->dmvx;	fmv1->aff2_dmvy = dv12->dmvy;
						fmv1->aff3_dmvx = dv13->dmvx;	fmv1->aff3_dmvy = dv13->dmvy;
					}//If MSE
				}
			  }
		  }
		}//if fmv2   

/*	 END OF TEST	*/

		}// if(pos1 >= 0 && pos2 >= 0)

//BI-DIRECTIONAL MERGE
//BI-DIRECTIONAL PARALLEL MERGE
		if(pos1 >= 0){
			best_sad = delta_v_search_prl(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, ctx1x, ctx1y, ctx2x, ctx2y, &get_v1[pred11], &get_v2[pred12], &get_v3[pred13], 
				dv11,dv12,dv13,dv21,dv22,dv23,xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF,info.lambda[t_level],info.subpel[t_level],t_level,0);

			get_aff_sad_cost = best_sad;
			
//			find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, *dv11, *dv12, *dv13, *dv21, *dv22, *dv23, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF, info.subpel[t_level],0);

			get_aff_bit_cost = info.lambda[t_level] * ( (int)(pos1/4)+AFF_IDX_OFFSET ) + (DIRECT_LEN+MERGE_LEN+MERGE_DIR_LEN) * info.lambda[t_level]
					  + get_bit_cost(info.lambda[t_level],dv11->mvx,dv11->mvy,get_v1[pred11].mvx,get_v1[pred11].mvy,ctx1x,ctx1y,info.subpel[t_level]) + 
					    get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,get_v2[pred12].mvx,get_v2[pred12].mvy,ctx1x,ctx1y,info.subpel[t_level]) + 
					    get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,get_v3[pred13].mvx,get_v3[pred13].mvy,ctx1x,ctx1y,info.subpel[t_level]) +
					    info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF,fmv2,info.bi_mv[t_level], t_level);

			get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

			if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
				get_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
				(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

				aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,*dv11, *dv12, *dv13, 
							*dv21, *dv22, *dv23, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF);

//				printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

				if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
					aff_mse < IBLOCK_FACTOR * aff_ref_var )
				    || aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

					assert(dv11->mvx == (-1)*dv21->mvx && dv11->mvy == (-1)*dv21->mvy && 
						   dv12->mvx == (-1)*dv22->mvx && dv12->mvy == (-1)*dv22->mvy &&
						   dv13->mvx == (-1)*dv23->mvx && dv13->mvy == (-1)*dv23->mvy );

					assert(fmv2->total_cost == fmv1->total_cost);

					fmv2->bi_mode = (BiMode)BI_CONNECTED_AFF;

					fmv2->sad_cost = get_aff_sad_cost;
					fmv2->mse = aff_mse;
					fmv2->bit_cost = get_aff_bit_cost;
					fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
					fmv2->aff_idx = -1;
					fmv2->direct_idx = INDIRECT;
					fmv2->merge_idx = MERGE;
					fmv2->merge_dir = PAL_L;
					fmv2->is_predictor = YES;
					fmv2->lifting_mode = CONNECTED;
					fmv2->aff1_mvx = dv21->mvx;		fmv2->aff1_mvy = dv21->mvy;
					fmv2->aff2_mvx = dv22->mvx;		fmv2->aff2_mvy = dv22->mvy;
					fmv2->aff3_mvx = dv23->mvx;		fmv2->aff3_mvy = dv23->mvy;

	//				printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n",fmv1->aff2_dmvx,fmv1->aff2_dmvy);

					fmv1->bi_mode = (BiMode)BI_CONNECTED_AFF;

					fmv1->sad_cost = fmv2->sad_cost;
					fmv1->mse = aff_mse;
					fmv1->bit_cost = fmv2->bit_cost;
					fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
					fmv1->aff_idx = pos1;
					fmv1->direct_idx = INDIRECT;
					fmv1->merge_idx = MERGE;
					fmv1->merge_dir = PAL_L;
					fmv1->is_predictor = YES;
					fmv1->lifting_mode = CONNECTED;
					fmv1->aff1_mvx = dv11->mvx;		fmv1->aff1_mvy = dv11->mvy;
					fmv1->aff2_mvx = dv12->mvx;		fmv1->aff2_mvy = dv12->mvy;
					fmv1->aff3_mvx = dv13->mvx;		fmv1->aff3_mvy = dv13->mvy;
			
					fmv1->aff1_dmvx = dv11->mvx - get_v1[pred11].mvx;
					fmv1->aff1_dmvy = dv11->mvy - get_v1[pred11].mvy;
					fmv1->aff2_dmvx = dv12->mvx - get_v2[pred12].mvx;
					fmv1->aff2_dmvy = dv12->mvy - get_v2[pred12].mvy;
					fmv1->aff3_dmvx = dv13->mvx - get_v3[pred13].mvx;
					fmv1->aff3_dmvy = dv13->mvy - get_v3[pred13].mvy;
				}//If MSE
			}
		}
//BI-DIRECTIONAL PARALLEL MERGE

//UP MERGE
		//	printf("BLOCK B\n");
		find_block(x,y-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,0,1);//BLOCK B	
		merg_aff_v1.mvx = *get1x;
		merg_aff_v1.mvy = *get1y;
		//	printf("BLOCK D\n");
		find_block(x+xblk2-1,y-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);//BLOCK D
		merg_aff_v2.mvx = *get1x;
		merg_aff_v2.mvy = *get1y;

		//	printf("BLOCK B\n");
		find_block(x,y-1,fmv2_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,0,1);//BLOCK B	
		merg_aff2_v1.mvx = *get2x;
		merg_aff2_v1.mvy = *get2y;
		//	printf("BLOCK D\n");
		find_block(x+xblk2-1,y-1,fmv2_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,1,1);//BLOCK D
		merg_aff2_v2.mvx = *get2x;
		merg_aff2_v2.mvy = *get2y;

		if(merg_aff_v1.mvx != (float)HUGE_VAL && merg_aff_v1.mvy != (float)HUGE_VAL && merg_aff_v2.mvx != (float)HUGE_VAL && 
			merg_aff_v2.mvy != (float)HUGE_VAL && merg_aff2_v1.mvx != (float)HUGE_VAL && merg_aff2_v1.mvy != (float)HUGE_VAL && 
			merg_aff2_v2.mvx != (float)HUGE_VAL && merg_aff2_v2.mvy != (float)HUGE_VAL){

			mrg_cnt1 = 0;

			for(i = 0;i < 2; i++){		

				mrg_cnt2 = 0;

				for(j = 0;j < 2; j++){	
					if(get_v3[i].mvx != (float)HUGE_VAL && get_v3[i].mvy != (float)HUGE_VAL && 
					   get2_v3[j].mvx != (float)HUGE_VAL && get2_v3[j].mvy != (float)HUGE_VAL){

						best_sad = delta_v_search_bi(fr_cur,fr_ref1,fr_ref2,upframe1,upframe2, ctx1x, ctx1y, ctx2x, ctx2y, &merg_aff_v1, &merg_aff_v2, &get_v3[i], &merg_aff2_v1,
							&merg_aff2_v2, &get2_v3[j], dv11,dv12,dv13,dv21,dv22,dv23,xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF,info.lambda[t_level],
						    info.subpel[t_level],t_level,3);

						get_aff_sad_cost = best_sad;
						
//						find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, merg_aff_v1, merg_aff_v2, *dv13, merg_aff2_v1, merg_aff2_v2,*dv23, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF,info.subpel[t_level],0);
						
						get_aff_bit_cost = info.lambda[t_level] * (DIRECT_LEN + MERGE_LEN + MERGE_DIR_UP_LEN + 1 + 1) + info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF,
						 fmv2,info.bi_mv[t_level],t_level) + get_bit_cost(info.lambda[t_level],dv13->mvx,dv13->mvy,get_v3[i].mvx,get_v3[i].mvy,ctx1x,
						  ctx1y,info.subpel[t_level]) + get_bit_cost(info.lambda[t_level],dv23->mvx,dv23->mvy,get2_v3[j].mvx,get2_v3[j].mvy,ctx2x,
						  ctx2y,info.subpel[t_level]);

						assert(fmv1->total_cost == fmv2->total_cost);

						get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

						if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
						get_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
						(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

							aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,merg_aff_v1, merg_aff_v2, *dv13,
								merg_aff2_v1, merg_aff2_v2,*dv23, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF);

//							printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

							if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
								aff_mse < IBLOCK_FACTOR * aff_ref_var )
								|| aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

								fmv2->bi_mode = (BiMode)BI_CONNECTED_AFF;

								fmv2->sad_cost = get_aff_sad_cost;
								fmv2->mse = aff_mse;
								fmv2->bit_cost = get_aff_bit_cost;
								fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
								fmv2->aff_idx = mrg_cnt2;
								fmv2->direct_idx = INDIRECT;
								fmv2->merge_idx = MERGE;
								fmv2->merge_dir = UP;
								fmv2->is_predictor = YES;
								fmv2->lifting_mode = CONNECTED;
								fmv2->aff1_mvx = merg_aff2_v1.mvx;fmv2->aff1_mvy = merg_aff2_v1.mvy;
								fmv2->aff2_mvx = merg_aff2_v2.mvx;fmv2->aff2_mvy = merg_aff2_v2.mvy;
								fmv2->aff3_mvx = dv23->mvx;		 fmv2->aff3_mvy = dv23->mvy;
	
								fmv2->aff3_dmvx = dv23->mvx - get2_v3[j].mvx;
								fmv2->aff3_dmvy = dv23->mvy - get2_v3[j].mvy;

	//							printf("fmv2->aff3_dmvx = %f, fmv2->aff3_dmvy = %f\n",fmv2->aff3_dmvx,fmv2->aff3_dmvy);

								fmv1->bi_mode = (BiMode)BI_CONNECTED_AFF;
		
								fmv1->sad_cost = fmv2->sad_cost;
								fmv1->mse = aff_mse;
								fmv1->bit_cost = fmv2->bit_cost;
								fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
								fmv1->aff_idx = mrg_cnt1;
								fmv1->direct_idx = INDIRECT;
								fmv1->merge_idx = MERGE;
								fmv1->merge_dir = UP;
								fmv1->is_predictor = YES;
								fmv1->lifting_mode = CONNECTED;
								fmv1->aff1_mvx = merg_aff_v1.mvx;fmv1->aff1_mvy = merg_aff_v1.mvy;
								fmv1->aff2_mvx = merg_aff_v2.mvx;fmv1->aff2_mvy = merg_aff_v2.mvy;
								fmv1->aff3_mvx = dv13->mvx;		 fmv1->aff3_mvy = dv13->mvy;
	
								fmv1->aff3_dmvx = dv13->mvx - get_v3[i].mvx;
								fmv1->aff3_dmvy = dv13->mvy - get_v3[i].mvy;
							}//If MSE
						}

						mrg_cnt2 ++;
					}
				}//j
				
				if(get_v3[i].mvx != (float)HUGE_VAL && get_v3[i].mvy != (float)HUGE_VAL)
					mrg_cnt1 ++;
			}//i
		}
//LEFT MERGE
		//	printf("BLOCK C\n");
		find_block(x-1,y,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,0);//BLOCK C	
		merg_aff_v1.mvx = *get1x;
		merg_aff_v1.mvy = *get1y;
		//	printf("BLOCK F\n");
		find_block(x-1,y+yblk2-1,fmv1_array,info,t_level,get1x,get1y,&get_xblk,&get_yblk,1,1);//BLOCK F
		merg_aff_v2.mvx = *get1x;
		merg_aff_v2.mvy = *get1y;

		//	printf("BLOCK C\n");
		find_block(x-1,y,fmv2_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,1,0);//BLOCK C	
		merg_aff2_v1.mvx = *get2x;
		merg_aff2_v1.mvy = *get2y;
		//	printf("BLOCK F\n");
		find_block(x-1,y+yblk2-1,fmv2_array,info,t_level,get2x,get2y,&get_xblk,&get_yblk,1,1);//BLOCK F
		merg_aff2_v2.mvx = *get2x;
		merg_aff2_v2.mvy = *get2y;

		if(merg_aff_v1.mvx != (float)HUGE_VAL && merg_aff_v1.mvy != (float)HUGE_VAL && merg_aff_v2.mvx != (float)HUGE_VAL && 
			merg_aff_v2.mvy != (float)HUGE_VAL && merg_aff2_v1.mvx != (float)HUGE_VAL && merg_aff2_v1.mvy != (float)HUGE_VAL && 
			merg_aff2_v2.mvx != (float)HUGE_VAL && merg_aff2_v2.mvy != (float)HUGE_VAL){

			mrg_cnt1 = 0;

			for(i = 0;i < 2; i++){		

				mrg_cnt2 = 0;

				for(j = 0;j < 2; j++){	
					if(get_v2[i].mvx != (float)HUGE_VAL && get_v2[i].mvy != (float)HUGE_VAL && 
					   get2_v2[j].mvx != (float)HUGE_VAL && get2_v2[j].mvy != (float)HUGE_VAL){

						best_sad = delta_v_search_bi(fr_cur,fr_ref1,fr_ref2,upframe1,upframe2, ctx1x, ctx1y, ctx2x, ctx2y, &merg_aff_v1, &get_v2[i], &merg_aff_v2, &merg_aff2_v1,
							&get2_v2[j], &merg_aff2_v2, dv11,dv12,dv13,dv21,dv22,dv23,xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF,info.lambda[t_level],
						    info.subpel[t_level],t_level,2);
						 
						get_aff_sad_cost = best_sad;
						
//						find_affine_SAD(fr_cur, fr_ref1, fr_ref2,upframe1,upframe2, merg_aff_v1, *dv12, merg_aff_v2, merg_aff2_v1, *dv22,merg_aff2_v2,xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF,info.subpel[t_level],0);
						
						get_aff_bit_cost = info.lambda[t_level] * (DIRECT_LEN + MERGE_LEN + MERGE_DIR_LEN + 1 + 1) + info.lambda[t_level] * get_mode_coding_cost(BI_CONNECTED_AFF,
						 fmv2,info.bi_mv[t_level],t_level) + get_bit_cost(info.lambda[t_level],dv12->mvx,dv12->mvy,get_v2[i].mvx,get_v2[i].mvy,ctx1x,
						  ctx1y,info.subpel[t_level]) + get_bit_cost(info.lambda[t_level],dv22->mvx,dv22->mvy,get2_v2[j].mvx,get2_v2[j].mvy,ctx2x,
						  ctx2y,info.subpel[t_level]);

						assert(fmv1->total_cost == fmv2->total_cost);

						get_aff_total_cost = get_aff_sad_cost + get_aff_bit_cost;

						if( ( (get_aff_total_cost < (fmv1->total_cost + 1 * info.lambda[t_level]) ) && fmv1->bi_mode <= 8 &&
							get_aff_sad_cost < (fmv1->sad_cost * aff_coef) ) ||
							(fmv1->bi_mode >= 9 && (get_aff_total_cost < fmv1->total_cost ) && get_aff_sad_cost < fmv1->sad_cost) ){  //affine decision

							aff_mse = find_affine_MSE(fr_cur,fr_ref1,fr_ref2,&aff_ref_var,merg_aff_v1, *dv12, merg_aff_v2,
								merg_aff2_v1, *dv22, merg_aff2_v2, xblk2,yblk2,x,y,hor,ver,BI_CONNECTED_AFF);

//							printf("aff_mse = %f, aff_ref_var = %f, aff_var = %f\n",aff_mse,aff_ref_var,aff_var);

							if( ( aff_mse < IBLOCK_FACTOR * aff_var && 
								aff_mse < IBLOCK_FACTOR * aff_ref_var )
								|| aff_mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ){

								fmv2->bi_mode = (BiMode)BI_CONNECTED_AFF;

								fmv2->sad_cost = get_aff_sad_cost;
								fmv2->mse = aff_mse;
								fmv2->bit_cost = get_aff_bit_cost;
								fmv2->total_cost = fmv2->sad_cost + fmv2->bit_cost;
								fmv2->aff_idx = mrg_cnt2;
								fmv2->direct_idx = INDIRECT;
								fmv2->merge_idx = MERGE;
								fmv2->merge_dir = LEFT;
								fmv2->is_predictor = YES;
								fmv2->lifting_mode = CONNECTED;
								fmv2->aff1_mvx = merg_aff2_v1.mvx;fmv2->aff1_mvy = merg_aff2_v1.mvy;
								fmv2->aff2_mvx = dv22->mvx;		  fmv2->aff2_mvy = dv22->mvy;
								fmv2->aff3_mvx = merg_aff2_v2.mvx;fmv2->aff3_mvy = merg_aff2_v2.mvy;
							
	
								fmv2->aff2_dmvx = dv22->mvx - get2_v2[j].mvx;
								fmv2->aff2_dmvy = dv22->mvy - get2_v2[j].mvy;

	//							printf("fmv2->aff3_dmvx = %f, fmv2->aff3_dmvy = %f\n",fmv2->aff3_dmvx,fmv2->aff3_dmvy);

								fmv1->bi_mode = (BiMode)BI_CONNECTED_AFF;
		
								fmv1->sad_cost = fmv2->sad_cost;
								fmv1->mse = aff_mse;
								fmv1->bit_cost = fmv2->bit_cost;
								fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
								fmv1->aff_idx = mrg_cnt1;
								fmv1->direct_idx = INDIRECT;
								fmv1->merge_idx = MERGE;
								fmv1->merge_dir = LEFT;
								fmv1->is_predictor = YES;
								fmv1->lifting_mode = CONNECTED;
								fmv1->aff1_mvx = merg_aff_v1.mvx;fmv1->aff1_mvy = merg_aff_v1.mvy;
								fmv1->aff2_mvx = dv12->mvx;		 fmv1->aff2_mvy = dv12->mvy;
								fmv1->aff3_mvx = merg_aff_v2.mvx;fmv1->aff3_mvy = merg_aff_v2.mvy;
							
	
								fmv1->aff2_dmvx = dv12->mvx - get_v2[i].mvx;
								fmv1->aff2_dmvy = dv12->mvy - get_v2[i].mvy;
							}//If MSE
						}

						mrg_cnt2 ++;
					}
				}//j

				if(get_v2[i].mvx != (float)HUGE_VAL && get_v2[i].mvy != (float)HUGE_VAL)
					mrg_cnt1 ++;
			}//i
		}

//MERGE MODE
	  free(block_buff1);
      free(block_buff2);

	}//if fmv2 != NULL

}// if best mode
  
  delete dv11;
  delete dv12;
  delete dv13;

  delete dv21;
  delete dv22;
  delete dv23;

  delete get1x;
  delete get1y;
  delete get2x;
  delete get2y;
  delete get3x;
  delete get3y;

#endif

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
  if (xblk==8 || xblk==4)  // directional IBLOCK detection
  {
	  float fadditional_penalty = 0; // to restrict the drift effect from directional IBLOCK
	  int   mode_bits; 
	  tmv1.mvx = tmv1.mvy = (float)HUGE_VAL;
	  tmv1.total_cost = tmv1.sad_cost = tmv1.bit_cost = (float)HUGE_VAL;
	  // fmv1_array the root of the mv quad-tree 
	  directional_iblock_detection(info, fmv1_array, fr_cur,  x, y, xblk, yblk, hor, ver, t_level, 
		                           &tmv1, &fadditional_penalty);
	  // the total cost has been updated, i.e. there is a valid spatial mode 
	  // otherwise there does not exist a valid spatial mode 
	  if (tmv1.total_cost!= (float)HUGE_VAL)  
	  {
		 mode_bits = get_mode_coding_cost(DIRECTIONAL_IBLOCK, fmv2, info.bi_mv[t_level], t_level );
		 assert(mode_bits>0);
		 mode_coding_cost = (DIRECTIONAL_IBLOCK_BIAS+fadditional_penalty)*info.lambda[t_level] * mode_bits;
         tmv1.bit_cost   += mode_coding_cost;
         tmv1.total_cost += mode_coding_cost;

	     // update the information in fmv1 and fmv2
	     v1 = &fmv1->mode_info[DIRECTIONAL_IBLOCK];
         v1->is_valid = YES;
		 v1->is_predictor = NO;
		 v1->mvx = v1->mvy = (float)HUGE_VAL;
	     v1->iblock_spatial_mode = tmv1.iblock_spatial_mode;
		 v1->sad_cost = tmv1.sad_cost;
		 v1->bit_cost = tmv1.bit_cost;
		 v1->total_cost = tmv1.total_cost;

		 if (fmv2 != NULL) {
			v2 = &fmv2->mode_info[DIRECTIONAL_IBLOCK];
			v2->is_valid = YES;
			v2->is_predictor = NO;
			v2->mvx = v2->mvy = (float)HUGE_VAL;
			v2->iblock_spatial_mode =  tmv1.iblock_spatial_mode;
			v2->sad_cost = v1->sad_cost;
			v2->bit_cost = v1->bit_cost;
			v2->total_cost = v1->total_cost;
		 }
		 // best mode is replaced by DIRECTIONAL_IBLOCK
		 if ( fmv1->mode_info[best_mode].total_cost>tmv1.total_cost)
		 {
			fmv1->bi_mode = DIRECTIONAL_IBLOCK;
			fmv1->lifting_mode = left_mode[DIRECTIONAL_IBLOCK];
			fmv1->is_predictor = NO;
			fmv1->mvx = fmv1->mvy = (float)HUGE_VAL;
			fmv1->iblock_spatial_mode = v1->iblock_spatial_mode;
			fmv1->sad_cost = v1->sad_cost;
			fmv1->bit_cost = v1->bit_cost;
			fmv1->total_cost = v1->total_cost;
			fmv1->propagate_iblk = (fadditional_penalty>0) ? 1: 0; 
			if (fmv2 != NULL) {
				fmv2->bi_mode = DIRECTIONAL_IBLOCK;
				fmv2->lifting_mode = right_mode[DIRECTIONAL_IBLOCK];
				fmv2->is_predictor = NO;
				fmv2->mvx = fmv2->mvy = (float)HUGE_VAL;
				fmv2->iblock_spatial_mode = v2->iblock_spatial_mode;
				fmv2->sad_cost = v2->sad_cost;
				fmv2->bit_cost = v2->bit_cost;
				fmv2->total_cost = v2->total_cost;
				fmv1->propagate_iblk = (fadditional_penalty>0) ? 1: 0; 
			}
		}
	  }
  }
#endif 

  delete tv11;
  delete tv12;
  delete tv13;

  delete tv21;
  delete tv22;
  delete tv23;

}