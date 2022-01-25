#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "structN.h"
#include "basic.h"
#include "bmeN.h"
#include "mv_ec.h"
#include "mv_statistics.h"
#include "bme_tools.h"
#include "mode_decision.h"

#define EXTERN extern

#include "layer_mv.h"

#include "ioN.h"
#include "miscN.h"

//#define EC_TYPE GOLOMB
#define EC_TYPE  AR_BINARY
#define EC_USE_CONTEXTS             /* use contexts */

#define Yes      1
#define No       0

#define DIS1     1              /* frame distance for block matching */
#define DIS2     2
#define DIS4     4

//#define   DEBUG_SCALABLE_MV           // for debug use by Yongjun Wu
//#define   DEBUG_BLOCK_MODE_MV_INFO    // for debug use by Yongjun Wu
//#define DEBUG_LAYER_STRUCTURE         // for debug use by Yongjun Wu
#define LAYER_BLOCK_SIZE 32

float avg_blk_size[4];
int blk_num;

int debug_counter=0;
FILE *fpAGP_debug; 
FILE *fpsub_debug; 
char decoder_AGPdebug[80];
char decoder_subdebug[80];
char encoder_AGPdebug[80];

//Added on 02.05.2018
int SIMUL_START;// = 16
int SIMUL_POINT;// = 0
int SIMUL_LIMIT;// = 1
int SIMUL_LEVEL;// = 0

int SIMUL_ROUND;
int SIMUL_RANGE;


VLCtable blockmodeVLC[2] = {
  {0, 1},                       //CONNECTED huffman code = 0
  {1, 1}                        //PREDICTED huffman code = 1
};



FRAME_MOTION_FIELD  *frame_motion_field; 
FRAME_MOTION_FIELD  *frame_motion_field2;

SIMP_FRAME_MOTION_FIELD  *prev_frame_motion_field1;
SIMP_FRAME_MOTION_FIELD  *prev_frame_motion_field2;

SIMP_FRAME_MOTION_FIELD  *prev_frame_motion_field_left1;
SIMP_FRAME_MOTION_FIELD  *prev_frame_motion_field_left2;

SIMP_FRAME_MOTION_FIELD	*buffer_frame_motion_field1[G_LEVEL];
SIMP_FRAME_MOTION_FIELD	*buffer_frame_motion_field2[G_LEVEL];

SIMP_FRAME_MOTION_FIELD	*save_buffer_frame_motion_field1[G_LEVEL];
SIMP_FRAME_MOTION_FIELD	*save_buffer_frame_motion_field2[G_LEVEL];


//Added on 01.10.2018, for simul decoding
SIMP_FRAME_MOTION_FIELD	*buffer_frame_motion_field1_dec[G_LEVEL];
SIMP_FRAME_MOTION_FIELD	*buffer_frame_motion_field2_dec[G_LEVEL];

SIMP_FRAME_MOTION_FIELD	*save_buffer_frame_motion_field1_dec[G_LEVEL];
SIMP_FRAME_MOTION_FIELD	*save_buffer_frame_motion_field2_dec[G_LEVEL];

int aff_cum_idx[27];

int mv_res_bits;

int mv_bits;

int use_huff;

float calc_res_bits;
float calc_res_sad;

extern enum FLAG **scene_change;
extern enum FLAG **dec_scene_change;

int ctl_info;
int mv_info;
int idx_info;

int aff_mrg_blk;

int cnt_quad[4];
int cnt_tri[3];

int bi_mode_num012[12];
int bi_mode_num345[12];

float mean_mode_mse[NUMBER_OF_BI_MODES];
int   mean_mode_num[NUMBER_OF_BI_MODES];

void child_map_decode( vector_ptr fmv1, vector_ptr fmv2,
                       int meandepth, int x, int y, int xblk, int yblk,
                       int hor, int ver, int small, videoinfo info, int t_level );

int read_number_core( FILE * fp );
void write_number_core( int outputbyte, FILE * fp );


// the variables for the section of sub-symbol
unsigned char     store_splitted_bytes_array[LAYER_NUM][MAX_AGP_LEVEL][65535];
unsigned int      splitted_mv_byte_num_array[LAYER_NUM][MAX_AGP_LEVEL];    // the number of bytes for a sub-symbol
unsigned char     splitted_mv_byte_array[LAYER_NUM][MAX_AGP_LEVEL];   // a specific byte for some sub-symbol
         char     splitted_bit_num_array[LAYER_NUM][MAX_AGP_LEVEL];   // how many bits have been used

// the variables for additional sign
unsigned char     store_splitted_sign[LAYER_NUM][65535];
unsigned int      splitted_sign_byte_num[LAYER_NUM];
unsigned char     splitted_sign_byte[LAYER_NUM];
         char     splitted_sign_bit_num[LAYER_NUM];


int frame_mv_cnt;
float sum_mv2, sum_mv1;

/*
 *       putbits()
 * write rightmost n (0<=n<=32) bits of val to outfile
 * val : value needs to be output.
 * n   : rightmost n bits of val
 */
void
putbits( int val, int n )
{
  int i;
  unsigned int mask;

  assert(n >= 0);

  mask = 1 << ( n - 1 );        /* selects first (leftmost) bit */

  for( i = 0; i < n; i++ ) {
    if( val & mask ) {
      output_bit( 1 );
    } else {
      output_bit( 0 );
    }
    mask >>= 1;                 /* select next bit */
  }
}


int
getbits( int len )
{
  int i, bit, val = 0;

  for( i = 0; i < len; i++ ) {
    val <<= 1;
    input_bit( bit );
    val |= bit;
  }
  return val;
}

//cx,cy - the coordinate of current block's upper-left corner
//x_pos,y_pos the merge reference point
void get_field_aff_mrg_mv( vector_ptr fmv, int cx, int cy, int x_pos, int y_pos, int xblk, int yblk, float dx1, float dy1,
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

int field_compare_mv(vector_ptr left1, vector_ptr right1, vector_ptr left2, vector_ptr right2){
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

/////////  added by Yuan Liu  ////////

void clean_mrg_mv(vector_ptr fmv){

	fmv->lifting_mode = IGNORED;

	fmv->mvx = (float)HUGE_VAL;
	fmv->mvy = (float)HUGE_VAL;

	fmv->aff1_mvx = (float)HUGE_VAL;
	fmv->aff1_mvy = (float)HUGE_VAL;
	fmv->aff2_mvx = (float)HUGE_VAL;
	fmv->aff2_mvy = (float)HUGE_VAL;
	fmv->aff3_mvx = (float)HUGE_VAL;
	fmv->aff3_mvy = (float)HUGE_VAL;

}

/*
	CBD
	AZ
	E
*/
void get_field_merge_mv_info(vector_ptr *mrg_left, vector_ptr *mrg_right, int x_pos, int y_pos, int xblk, int yblk,
	int hor, int ver, videoinfo info, int t_level, FRAME_MOTION_FIELD  *cur_frame_motion_field1, FRAME_MOTION_FIELD *cur_frame_motion_field2,
	SIMP_FRAME_MOTION_FIELD  *prev_frame_motion_field1, SIMP_FRAME_MOTION_FIELD  *prev_frame_motion_field2, int type){
	int i,j,k;
	int x_dest, y_dest, x_dest2, y_dest2;
	int get_xblk, get_yblk;
	int is_predl_a, is_predl_b, is_predl_c, is_predl_d, is_predl_tmp, is_predl_e;
	int is_predr_a, is_predr_b, is_predr_c, is_predr_d, is_predr_tmp, is_predr_e;

	int aff_blk; //Added on 09.23.2016

	float ul1,ul2;
	float med_l_px[4], med_l_py[4], med_r_px[4], med_r_py[4];
	float med_l_aff1x[4],med_l_aff1y[4],med_l_aff2x[4],med_l_aff2y[4],med_l_aff3x[4],med_l_aff3y[4];
	float med_r_aff1x[4],med_r_aff1y[4],med_r_aff2x[4],med_r_aff2y[4],med_r_aff3x[4],med_r_aff3y[4];

	float left_affx1, left_affy1, left_affx2, left_affy2, left_affx3, left_affy3, left_affx4, left_affy4;
	float right_affx1, right_affy1, right_affx2, right_affy2, right_affx3, right_affy3, right_affx4, right_affy4;

	float aff_dmvx1, aff_dmvy1;//x-directional affine dmv
	float aff_dmvx2, aff_dmvy2;//y-directional affine dmv

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
	x_dest = x_pos - 1;
	y_dest = y_pos + yblk - 1;

	if(x_dest < 0 || y_dest >= ver){
		is_predl_a = 0;
		is_predr_a = 0;
	}else{
		is_predl_a = cur_frame_motion_field1[y_dest*hor+x_dest].available;

		if(type == 0)
			is_predr_a = cur_frame_motion_field2[y_dest*hor+x_dest].available;
		else
			is_predr_a = 0;

		if(is_predl_a == 1){
			mrg_left[0]->mvx = cur_frame_motion_field1[y_dest*hor+x_dest].mvx;
			mrg_left[0]->mvy = cur_frame_motion_field1[y_dest*hor+x_dest].mvy;

			left_affx1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-2)].mvx;
			left_affy1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-2)].mvy;

			left_affx2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-1)].mvx;
			left_affy2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-1)].mvy;

			left_affx3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest-2)].mvx;
			left_affy3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest-2)].mvy;

			if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_field_aff_mrg_mv(mrg_left[0], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
			
				aff_blk = 1;
			}

			mrg_left[0]->lifting_mode = CONNECTED;
		}

		if(is_predr_a == 1){
			mrg_right[0]->mvx = cur_frame_motion_field2[y_dest*hor+x_dest].mvx;
			mrg_right[0]->mvy = cur_frame_motion_field2[y_dest*hor+x_dest].mvy;

			right_affx1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-2)].mvx;
			right_affy1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-2)].mvy;

			right_affx2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-1)].mvx;
			right_affy2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-1)].mvy;

			right_affx3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest-2)].mvx;
			right_affy3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest-2)].mvy;

			if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_field_aff_mrg_mv(mrg_right[0], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
				
				aff_blk = 1;
			}

			mrg_right[0]->lifting_mode = CONNECTED;
		}

		if(is_predl_a == 1 && is_predr_a == 1 && aff_blk == 1){
			aff_dmvx1 = (left_affx2 - left_affx1);
			aff_dmvy1 = (left_affy2 - left_affy1);

			aff_dmvx2 = (left_affx3 - left_affx1);
			aff_dmvy2 = (left_affy3 - left_affy1);

			get_field_aff_mrg_mv(mrg_left[0], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
			aff_dmvx1 = (right_affx2 - right_affx1);
			aff_dmvy1 = (right_affy2 - right_affy1);

			aff_dmvx2 = (right_affx3 - right_affx1);
			aff_dmvy2 = (right_affy3 - right_affy1);

			get_field_aff_mrg_mv(mrg_right[0], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
		}

	}

//B
	aff_blk = 0;
	x_dest = x_pos + xblk - 1;
	y_dest = y_pos - 1;

	if(y_dest < 0 || x_dest >= hor){
		is_predl_b = 0;
		is_predr_b = 0;
	}else{
		is_predl_b = cur_frame_motion_field1[y_dest*hor+x_dest].available;

		if(type == 0)
			is_predr_b = cur_frame_motion_field2[y_dest*hor+x_dest].available;
		else
			is_predr_b = 0;

		if(is_predl_b == 1){
			mrg_left[1]->mvx = cur_frame_motion_field1[y_dest*hor+x_dest].mvx;
			mrg_left[1]->mvy = cur_frame_motion_field1[y_dest*hor+x_dest].mvy;
			
			left_affx1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-2)].mvx;
			left_affy1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-2)].mvy;

			left_affx2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-1)].mvx;
			left_affy2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-1)].mvy;

			left_affx3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest-2)].mvx;
			left_affy3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest-2)].mvy;

			if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_field_aff_mrg_mv(mrg_left[1], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
				aff_blk = 1;
			}

			mrg_left[1]->lifting_mode = CONNECTED;
		}

		if(is_predr_b == 1){
			mrg_right[1]->mvx = cur_frame_motion_field2[y_dest*hor+x_dest].mvx;
			mrg_right[1]->mvy = cur_frame_motion_field2[y_dest*hor+x_dest].mvy;
			
			right_affx1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-2)].mvx;
			right_affy1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-2)].mvy;

			right_affx2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-1)].mvx;
			right_affy2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-1)].mvy;

			right_affx3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest-2)].mvx;
			right_affy3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest-2)].mvy;

			if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_field_aff_mrg_mv(mrg_right[1], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
			
				aff_blk = 1;
			}

			mrg_right[1]->lifting_mode = CONNECTED;
		}

		if(is_predl_b == 1 && is_predr_b == 1 && aff_blk == 1){
			aff_dmvx1 = (left_affx2 - left_affx1);
			aff_dmvy1 = (left_affy2 - left_affy1);

			aff_dmvx2 = (left_affx3 - left_affx1);
			aff_dmvy2 = (left_affy3 - left_affy1);

			get_field_aff_mrg_mv(mrg_left[1], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
			aff_dmvx1 = (right_affx2 - right_affx1);
			aff_dmvy1 = (right_affy2 - right_affy1);

			aff_dmvx2 = (right_affx3 - right_affx1);
			aff_dmvy2 = (right_affy3 - right_affy1);

			get_field_aff_mrg_mv(mrg_right[1], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
		}

	}
//C
	aff_blk = 0;
	x_dest = x_pos - 1;
	y_dest = y_pos - 1;

	if(x_dest < 0 || y_dest < 0){
		is_predl_c = 0;
		is_predr_c = 0;
	}else{
		is_predl_c = cur_frame_motion_field1[y_dest*hor+x_dest].available;

		if(type == 0)
			is_predr_c = cur_frame_motion_field2[y_dest*hor+x_dest].available;
		else
			is_predr_c = 0;

		if(is_predl_c == 1){
			mrg_left[2]->mvx = cur_frame_motion_field1[y_dest*hor+x_dest].mvx;
			mrg_left[2]->mvy = cur_frame_motion_field1[y_dest*hor+x_dest].mvy;

			left_affx1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-2)].mvx;
			left_affy1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-2)].mvy;

			left_affx2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-1)].mvx;
			left_affy2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest-1)].mvy;

			left_affx3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest-2)].mvx;
			left_affy3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest-2)].mvy;

			if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_field_aff_mrg_mv(mrg_left[2], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
			
				aff_blk = 1;
			}

			mrg_left[2]->lifting_mode = CONNECTED;
		}

		if(is_predr_c == 1){
			mrg_right[2]->mvx = cur_frame_motion_field2[y_dest*hor+x_dest].mvx;
			mrg_right[2]->mvy = cur_frame_motion_field2[y_dest*hor+x_dest].mvy;

			right_affx1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-2)].mvx;
			right_affy1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-2)].mvy;

			right_affx2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-1)].mvx;
			right_affy2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest-1)].mvy;

			right_affx3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest-2)].mvx;
			right_affy3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest-2)].mvy;

			if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_field_aff_mrg_mv(mrg_right[2], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
			
				aff_blk = 1;
			}

			mrg_right[2]->lifting_mode = CONNECTED;
		}

		if(is_predl_c == 1 && is_predr_c == 1 && aff_blk == 1){
			aff_dmvx1 = (left_affx2 - left_affx1);
			aff_dmvy1 = (left_affy2 - left_affy1);

			aff_dmvx2 = (left_affx3 - left_affx1);
			aff_dmvy2 = (left_affy3 - left_affy1);

			get_field_aff_mrg_mv(mrg_left[2], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
			aff_dmvx1 = (right_affx2 - right_affx1);
			aff_dmvy1 = (right_affy2 - right_affy1);

			aff_dmvx2 = (right_affx3 - right_affx1);
			aff_dmvy2 = (right_affy3 - right_affy1);

			get_field_aff_mrg_mv(mrg_right[2], x_pos, y_pos, x_dest-2, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
		}
	}
//D
	aff_blk = 0;
	if( (is_predl_a == 1 || is_predr_a == 1) && (is_predl_b == 1 || is_predr_b == 1) && (is_predl_c == 1 || is_predr_c == 1) );
	else{
		x_dest = x_pos + xblk;
		y_dest = y_pos - 1;

		if(x_dest >= hor || y_dest < 0){
			is_predl_d = 0;
			is_predr_d = 0;
		}else{
			is_predl_d = cur_frame_motion_field1[y_dest*hor+x_dest].available;

			if(type == 0)
				is_predr_d = cur_frame_motion_field2[y_dest*hor+x_dest].available;
			else
				is_predr_d = 0;

			if(is_predl_d == 1){
				mrg_left[3]->mvx = cur_frame_motion_field1[y_dest*hor+x_dest].mvx;
				mrg_left[3]->mvy = cur_frame_motion_field1[y_dest*hor+x_dest].mvy;

				left_affx1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest+1)].mvx;
				left_affy1 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest+1)].mvy;

				left_affx2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest+2)].mvx;
				left_affy2 = cur_frame_motion_field1[(y_dest-2)*hor+(x_dest+2)].mvy;

				left_affx3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest+1)].mvx;
				left_affy3 = cur_frame_motion_field1[(y_dest-1)*hor+(x_dest+1)].mvy;

				if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){
					aff_dmvx1 = (left_affx2 - left_affx1);
					aff_dmvy1 = (left_affy2 - left_affy1);

					aff_dmvx2 = (left_affx3 - left_affx1);
					aff_dmvy2 = (left_affy3 - left_affy1);

					get_field_aff_mrg_mv(mrg_left[3], x_pos, y_pos, x_dest+1, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
					aff_blk = 1;
				}
				mrg_left[3]->lifting_mode = CONNECTED;
			}

			if(is_predr_d == 1){
				mrg_right[3]->mvx = cur_frame_motion_field2[y_dest*hor+x_dest].mvx;
				mrg_right[3]->mvy = cur_frame_motion_field2[y_dest*hor+x_dest].mvy;

				right_affx1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest+1)].mvx;
				right_affy1 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest+1)].mvy;

				right_affx2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest+2)].mvx;
				right_affy2 = cur_frame_motion_field2[(y_dest-2)*hor+(x_dest+2)].mvy;

				right_affx3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest+1)].mvx;
				right_affy3 = cur_frame_motion_field2[(y_dest-1)*hor+(x_dest+1)].mvy;

				if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){
					aff_dmvx1 = (right_affx2 - right_affx1);
					aff_dmvy1 = (right_affy2 - right_affy1);

					aff_dmvx2 = (right_affx3 - right_affx1);
					aff_dmvy2 = (right_affy3 - right_affy1);

					get_field_aff_mrg_mv(mrg_right[3], x_pos, y_pos, x_dest+1, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
				
					aff_blk = 1;
				}

				mrg_right[3]->lifting_mode = CONNECTED;
			}

			if(is_predl_d == 1 && is_predr_d == 1 && aff_blk == 1){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_field_aff_mrg_mv(mrg_left[3], x_pos, y_pos, x_dest+1, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_field_aff_mrg_mv(mrg_right[3], x_pos, y_pos, x_dest+1, y_dest-2,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
			}

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
    for(i = 0; i <= 3; i++){
		if( (mrg_left[i]->mvx == (float)HUGE_VAL && mrg_left[i]->mvy == (float)HUGE_VAL) &&
			(mrg_right[i]->mvx == (float)HUGE_VAL && mrg_right[i]->mvy == (float)HUGE_VAL) ){
		  x_dest = x_pos;
		  y_dest = y_pos;
		  if (x_dest < 0 || y_dest < 0) {
			assert(0);
		  } else {
			  is_predl_tmp = prev_frame_motion_field1[y_dest*hor+x_dest].available;

			  if(type == 0)
				is_predr_tmp = prev_frame_motion_field2[y_dest*hor+x_dest].available;
			  else
				is_predr_tmp = 0;

			  if(is_predl_tmp == 1){
				mrg_left[i]->mvx = prev_frame_motion_field1[y_dest*hor+x_dest].mvx;
				mrg_left[i]->mvy = prev_frame_motion_field1[y_dest*hor+x_dest].mvy;
				mrg_left[i]->lifting_mode = CONNECTED;
			  }
			  if(is_predr_tmp == 1){
				mrg_right[i]->mvx = prev_frame_motion_field2[y_dest*hor+x_dest].mvx;
				mrg_right[i]->mvy = prev_frame_motion_field2[y_dest*hor+x_dest].mvy;
				mrg_right[i]->lifting_mode = CONNECTED;
			  }
		  }
//Additional TMP merge candidate
		  if( is_predl_tmp == 0 && is_predr_tmp == 0 ){
			x_dest = x_pos - xblk;
			y_dest = y_pos - yblk;
			if (x_dest < 0 || y_dest < 0) {
			  is_predl_tmp = 0;
			  is_predr_tmp = 0;
		    } else {
			    is_predl_tmp = prev_frame_motion_field1[y_dest*hor+x_dest].available;

				if(type == 0)
					is_predr_tmp = prev_frame_motion_field2[y_dest*hor+x_dest].available;
				else
					is_predr_tmp = 0;
			    
				if(is_predl_tmp == 1){
				  mrg_left[i]->mvx = prev_frame_motion_field1[y_dest*hor+x_dest].mvx;
				  mrg_left[i]->mvy = prev_frame_motion_field1[y_dest*hor+x_dest].mvy;
				  mrg_left[i]->lifting_mode = CONNECTED;
			    }
			    if(is_predr_tmp == 1){
				  mrg_right[i]->mvx = prev_frame_motion_field2[y_dest*hor+x_dest].mvx;
				  mrg_right[i]->mvy = prev_frame_motion_field2[y_dest*hor+x_dest].mvy;
				  mrg_right[i]->lifting_mode = CONNECTED;
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
		  x_dest = x_pos - 1;
		  y_dest = y_pos + yblk;
		  if(x_dest < 0 || y_dest >= ver){
			is_predl_e = 0;
			is_predr_e = 0;
		  }else{
			is_predl_e = cur_frame_motion_field1[y_dest*hor+x_dest].available;

			if(type == 0)
				is_predr_e = cur_frame_motion_field2[y_dest*hor+x_dest].available;
			else
				is_predr_e = 0;

			if(is_predl_e == 1){
				mrg_left[i]->mvx = cur_frame_motion_field1[y_dest*hor+x_dest].mvx;
				mrg_left[i]->mvy = cur_frame_motion_field1[y_dest*hor+x_dest].mvy;

				left_affx1 = cur_frame_motion_field1[(y_dest+1)*hor+(x_dest-2)].mvx;
				left_affy1 = cur_frame_motion_field1[(y_dest+1)*hor+(x_dest-2)].mvy;

				left_affx2 = cur_frame_motion_field1[(y_dest+1)*hor+(x_dest-1)].mvx;
				left_affy2 = cur_frame_motion_field1[(y_dest+1)*hor+(x_dest-1)].mvy;

				left_affx3 = cur_frame_motion_field1[(y_dest+2)*hor+(x_dest-2)].mvx;
				left_affy3 = cur_frame_motion_field1[(y_dest+2)*hor+(x_dest-2)].mvy;

				if( (left_affx3 != left_affx1 || left_affy3 != left_affy1) || (left_affx2 != left_affx1 || left_affy2 != left_affy1) ){
					aff_dmvx1 = (left_affx2 - left_affx1);
					aff_dmvy1 = (left_affy2 - left_affy1);

					aff_dmvx2 = (left_affx3 - left_affx1);
					aff_dmvy2 = (left_affy3 - left_affy1);

					get_field_aff_mrg_mv(mrg_left[i], x_pos, y_pos, x_dest-2, y_dest+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
					aff_blk = 1;
				}
				mrg_left[i]->lifting_mode = CONNECTED;
			}

			if(is_predr_e == 1){
				mrg_right[i]->mvx = cur_frame_motion_field2[y_dest*hor+x_dest].mvx;
				mrg_right[i]->mvy = cur_frame_motion_field2[y_dest*hor+x_dest].mvy;

				right_affx1 = cur_frame_motion_field2[(y_dest+1)*hor+(x_dest-2)].mvx;
				right_affy1 = cur_frame_motion_field2[(y_dest+1)*hor+(x_dest-2)].mvy;

				right_affx2 = cur_frame_motion_field2[(y_dest+1)*hor+(x_dest-1)].mvx;
				right_affy2 = cur_frame_motion_field2[(y_dest+1)*hor+(x_dest-1)].mvy;

				right_affx3 = cur_frame_motion_field2[(y_dest+2)*hor+(x_dest-2)].mvx;
				right_affy3 = cur_frame_motion_field2[(y_dest+2)*hor+(x_dest-2)].mvy;

				if( (right_affx3 != right_affx1 || right_affy3 != right_affy1) || (right_affx2 != right_affx1 || right_affy2 != right_affy1) ){
					aff_dmvx1 = (right_affx2 - right_affx1);
					aff_dmvy1 = (right_affy2 - right_affy1);

					aff_dmvx2 = (right_affx3 - right_affx1);
					aff_dmvy2 = (right_affy3 - right_affy1);

					get_field_aff_mrg_mv(mrg_right[i], x_pos, y_pos, x_dest-2, y_dest+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
				
					aff_blk = 1;
				}
				mrg_right[i]->lifting_mode = CONNECTED;
			}

			if(is_predl_e == 1 && is_predr_e == 1 && aff_blk == 1){
				aff_dmvx1 = (left_affx2 - left_affx1);
				aff_dmvy1 = (left_affy2 - left_affy1);

				aff_dmvx2 = (left_affx3 - left_affx1);
				aff_dmvy2 = (left_affy3 - left_affy1);

				get_field_aff_mrg_mv(mrg_left[i], x_pos, y_pos, x_dest-2, y_dest+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,left_affx1,left_affy1);
				
				aff_dmvx1 = (right_affx2 - right_affx1);
				aff_dmvy1 = (right_affy2 - right_affy1);

				aff_dmvx2 = (right_affx3 - right_affx1);
				aff_dmvy2 = (right_affy3 - right_affy1);

				get_field_aff_mrg_mv(mrg_right[i], x_pos, y_pos, x_dest-2, y_dest+1,xblk,yblk,aff_dmvx1,aff_dmvy1,aff_dmvx2,aff_dmvy2,right_affx1,right_affy1);
			}
		  }
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

}


float medi(float a, float b, float c)
{
  return ((a>b) ? ((a>c) ? ((b>c) ? b : c) : a) : 
          ((b>c) ? ((a>c) ? a : c) : b));
}

int put_aff_mvs(vector_ptr fmv, int x, int y, int xblk, int yblk, int hor, int ver, FRAME_MOTION_FIELD* cur_frame_motion_field, SIMP_FRAME_MOTION_FIELD* prev_frame_motion_field, int enc_idx){
	int getnum, count = 0, length, aff_num = 0, left_num = 0, up_num = 0;
	int i,j,k;
	float aff1_mvx[3],aff1_mvy[3],aff2_mvx[3],aff2_mvy[3],aff3_mvx[3],aff3_mvy[3];
	float get1_mvx[3],get1_mvy[3],get2_mvx[3],get2_mvy[3],get3_mvx[3],get3_mvy[3];
	float tempx, tempy;
	int int_mvx, int_mvy;
	float sub = 0.25;

	length = 0;

	for(i=0;i<3;i++){
		aff1_mvx[i] = (float)HUGE_VAL;get1_mvx[i] = (float)HUGE_VAL;
		aff1_mvy[i] = (float)HUGE_VAL;get1_mvy[i] = (float)HUGE_VAL;
		aff2_mvx[i] = (float)HUGE_VAL;get2_mvx[i] = (float)HUGE_VAL;
		aff2_mvy[i] = (float)HUGE_VAL;get2_mvy[i] = (float)HUGE_VAL;
		aff3_mvx[i] = (float)HUGE_VAL;get3_mvx[i] = (float)HUGE_VAL;
		aff3_mvy[i] = (float)HUGE_VAL;get3_mvy[i] = (float)HUGE_VAL;
	}

// MV1
	if(y-1>=0){
		aff1_mvx[0] = cur_frame_motion_field[(y-1)*hor + x].mvx;
		aff1_mvy[0] = cur_frame_motion_field[(y-1)*hor + x].mvy;

		aff1_mvx[0] = cur_frame_motion_field[(y-2)*hor + x].mvx + 2 * (cur_frame_motion_field[(y-2)*hor + x].mvx - cur_frame_motion_field[(y-3)*hor + x].mvx);
		aff1_mvy[0] = cur_frame_motion_field[(y-2)*hor + x].mvy + 2 * (cur_frame_motion_field[(y-2)*hor + x].mvy - cur_frame_motion_field[(y-3)*hor + x].mvy);

		aff1_mvx[0] = aff1_mvx[0] * AFF_SUBPEL;
		aff1_mvy[0] = aff1_mvy[0] * AFF_SUBPEL;
		int_mvx = (int)(aff1_mvx[0]);
		int_mvy = (int)(aff1_mvy[0]);
		aff1_mvx[0] = (float)(int_mvx);
		aff1_mvy[0] = (float)(int_mvy);
		aff1_mvx[0] = aff1_mvx[0] / AFF_SUBPEL;
		aff1_mvy[0] = aff1_mvy[0] / AFF_SUBPEL;
		
		if(cur_frame_motion_field[(y-2)*hor + x].mvx == (float)HUGE_VAL){
			aff1_mvx[0] = (float)HUGE_VAL;
			aff1_mvy[0] = (float)HUGE_VAL;
		}

//		if(aff1_mvx[0]!= get1_mvx[0] || aff1_mvy[0]!= get1_mvy[0])
//			printf("aff1_mvx0 = %f, aff1_mvy0 = %f, get1_mvx0 = %f, get1_mvy0 = %f\n\n",aff1_mvx[0],aff1_mvy[0],get1_mvx[0],get1_mvy[0]);
	}

	if(x-1>=0){
		aff1_mvx[1] = cur_frame_motion_field[y*hor + (x-1)].mvx;
		aff1_mvy[1] = cur_frame_motion_field[y*hor + (x-1)].mvy;

		aff1_mvx[1] = cur_frame_motion_field[y*hor + (x-2)].mvx + 2 * (cur_frame_motion_field[y*hor + (x-2)].mvx - cur_frame_motion_field[y*hor + (x-3)].mvx);
		aff1_mvy[1] = cur_frame_motion_field[y*hor + (x-2)].mvy + 2 * (cur_frame_motion_field[y*hor + (x-2)].mvy - cur_frame_motion_field[y*hor + (x-3)].mvy);

		aff1_mvx[1] = aff1_mvx[1] * AFF_SUBPEL;
		aff1_mvy[1] = aff1_mvy[1] * AFF_SUBPEL;
		int_mvx = (int)(aff1_mvx[1]);
		int_mvy = (int)(aff1_mvy[1]);
		aff1_mvx[1] = (float)(int_mvx);
		aff1_mvy[1] = (float)(int_mvy);
		aff1_mvx[1] = aff1_mvx[1] / AFF_SUBPEL;
		aff1_mvy[1] = aff1_mvy[1] / AFF_SUBPEL;
		
		if(cur_frame_motion_field[y*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff1_mvx[1] = (float)HUGE_VAL;
			aff1_mvy[1] = (float)HUGE_VAL;
		}

//		if(aff1_mvx[1] != get1_mvx[1] || aff1_mvy[1] != get1_mvy[1])
//			printf("aff1_mvx1 = %f, aff1_mvy1 = %f, get1_mvx1 = %f, get1_mvy1 = %f\n\n",aff1_mvx[1],aff1_mvy[1],get1_mvx[1],get1_mvy[1]);
	}

	if(y-1>=0 && x-1>=0){
		aff1_mvx[2] = cur_frame_motion_field[(y-1)*hor + (x-1)].mvx;
		aff1_mvy[2] = cur_frame_motion_field[(y-1)*hor + (x-1)].mvy;

		aff1_mvx[2] = cur_frame_motion_field[(y-2)*hor + (x-2)].mvx + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y-3)*hor + (x-2)].mvx) + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y-2)*hor + (x-3)].mvx);
		aff1_mvy[2] = cur_frame_motion_field[(y-2)*hor + (x-2)].mvy + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y-3)*hor + (x-2)].mvy) + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y-2)*hor + (x-3)].mvy);
		
		aff1_mvx[2] = aff1_mvx[2] * AFF_SUBPEL;
		aff1_mvy[2] = aff1_mvy[2] * AFF_SUBPEL;
		int_mvx = (int)(aff1_mvx[2]);
		int_mvy = (int)(aff1_mvy[2]);
		aff1_mvx[2] = (float)(int_mvx);
		aff1_mvy[2] = (float)(int_mvy);
		aff1_mvx[2] = aff1_mvx[2] / AFF_SUBPEL;
		aff1_mvy[2] = aff1_mvy[2] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y-2)*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff1_mvx[2] = (float)HUGE_VAL;
			aff1_mvy[2] = (float)HUGE_VAL;
		}

//		if(aff1_mvx[2] != get1_mvx[2] || aff1_mvy[2] != get1_mvy[2])
//			printf("aff1_mvx2 = %f, aff1_mvy2 = %f, get1_mvx2 = %f, get1_mvy2 = %f\n\n",aff1_mvx[2],aff1_mvy[2],get1_mvx[2],get1_mvy[2]);
	}

	if( fabs(aff1_mvx[1] - aff1_mvx[0]) < SMALL_DIFF && fabs(aff1_mvy[1] - aff1_mvy[0]) < SMALL_DIFF && aff1_mvx[1] != (float)HUGE_VAL && aff1_mvy[1] != (float)HUGE_VAL){
		aff1_mvx[1] = (float)HUGE_VAL;
		aff1_mvy[1] = (float)HUGE_VAL;
	}
	if( ( ( fabs(aff1_mvx[2] - aff1_mvx[0]) < SMALL_DIFF && fabs(aff1_mvy[2] - aff1_mvy[0]) < SMALL_DIFF) || ( fabs(aff1_mvx[2] - aff1_mvx[1]) < SMALL_DIFF && fabs(aff1_mvy[2] - aff1_mvy[1]) < SMALL_DIFF) )
		&& (aff1_mvx[2] != (float)HUGE_VAL && aff1_mvy[2] != (float)HUGE_VAL) ){
		aff1_mvx[2] = (float)HUGE_VAL;
		aff1_mvy[2] = (float)HUGE_VAL;
	}

//TEMP
	for( i = 0; i < 3 ; i ++ ){
		if( aff1_mvx[i] == (float)HUGE_VAL && aff1_mvy[i] == (float)HUGE_VAL ){
			tempx = prev_frame_motion_field[y*hor + x].mvx;
			tempy = prev_frame_motion_field[y*hor + x].mvy;

			for( j = 0; j < 3 ; j ++ ){
				if( tempx != (float)HUGE_VAL && tempy != (float)HUGE_VAL && (fabs(tempx - aff1_mvx[j]) < SMALL_DIFF && fabs(tempy - aff1_mvy[j]) < SMALL_DIFF) ){
					tempx = (float)HUGE_VAL;
					tempy = (float)HUGE_VAL;
				}
			}
			aff1_mvx[i] = tempx;
			aff1_mvy[i] = tempy;
			break;
		}
	}

// MV2
	if(y-1>=0 && x+xblk-1<hor){
		aff2_mvx[0] = cur_frame_motion_field[(y-1)*hor + (x+xblk-1)].mvx;
		aff2_mvy[0] = cur_frame_motion_field[(y-1)*hor + (x+xblk-1)].mvy;

		aff2_mvx[0] = cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx - cur_frame_motion_field[(y-3)*hor + (x+xblk-2)].mvx) + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx - cur_frame_motion_field[(y-2)*hor + (x+xblk-3)].mvx);
		aff2_mvy[0] = cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy - cur_frame_motion_field[(y-3)*hor + (x+xblk-2)].mvy) + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy - cur_frame_motion_field[(y-2)*hor + (x+xblk-3)].mvy);
		
		aff2_mvx[0] = aff2_mvx[0] * AFF_SUBPEL;
		aff2_mvy[0] = aff2_mvy[0] * AFF_SUBPEL;
		int_mvx = (int)(aff2_mvx[0]);
		int_mvy = (int)(aff2_mvy[0]);
		aff2_mvx[0] = (float)(int_mvx);
		aff2_mvy[0] = (float)(int_mvy);
		aff2_mvx[0] = aff2_mvx[0] / AFF_SUBPEL;
		aff2_mvy[0] = aff2_mvy[0] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx == (float)HUGE_VAL){
			aff2_mvx[0] = (float)HUGE_VAL;
			aff2_mvy[0] = (float)HUGE_VAL;
		}

//		if(aff2_mvx[0]!=get2_mvx[0] || aff2_mvy[0]!=get2_mvy[0])
//			printf("aff2_mvx0 = %f, aff2_mvy0 = %f, get2_mvx0 = %f, get2_mvy0 = %f\n\n",aff2_mvx[0],aff2_mvy[0],get2_mvx[0],get2_mvy[0]);
	}
	if(y-1>=0 && x+xblk<hor){
		aff2_mvx[1] = cur_frame_motion_field[(y-1)*hor + (x+xblk)].mvx;
		aff2_mvy[1] = cur_frame_motion_field[(y-1)*hor + (x+xblk)].mvy;

		aff2_mvx[1] = cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvx + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvx - cur_frame_motion_field[(y-3)*hor + (x+xblk)].mvx);
		aff2_mvy[1] = cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvy + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvy - cur_frame_motion_field[(y-3)*hor + (x+xblk)].mvy);

		aff2_mvx[1] = aff2_mvx[1] * AFF_SUBPEL;
		aff2_mvy[1] = aff2_mvy[1] * AFF_SUBPEL;
		int_mvx = (int)(aff2_mvx[1]);
		int_mvy = (int)(aff2_mvy[1]);
		aff2_mvx[1] = (float)(int_mvx);
		aff2_mvy[1] = (float)(int_mvy);
		aff2_mvx[1] = aff2_mvx[1] / AFF_SUBPEL;
		aff2_mvy[1] = aff2_mvy[1] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvx == (float)HUGE_VAL){
			aff2_mvx[1] = (float)HUGE_VAL;
			aff2_mvy[1] = (float)HUGE_VAL;
		}

//		if(aff2_mvx[1]!=get2_mvx[1] || aff2_mvy[1]!=get2_mvy[1])
//			printf("aff2_mvx1 = %f, aff2_mvy1 = %f, get2_mvx1 = %f, get2_mvy1 = %f\n\n",aff2_mvx[1],aff2_mvy[1],get2_mvx[1],get2_mvy[1]);
	}

	if( fabs(aff2_mvx[1] - aff2_mvx[0]) < SMALL_DIFF && fabs(aff2_mvy[1] - aff2_mvy[0]) < SMALL_DIFF && aff2_mvx[1] != (float)HUGE_VAL && aff2_mvy[1] != (float)HUGE_VAL){
		aff2_mvx[1] = (float)HUGE_VAL;
		aff2_mvy[1] = (float)HUGE_VAL;
	}

	if(aff2_mvx[0] != (float)HUGE_VAL && aff2_mvy[0] != (float)HUGE_VAL && aff2_mvx[1] != (float)HUGE_VAL && aff2_mvy[1] != (float)HUGE_VAL){
		aff2_mvx[2] = (aff2_mvx[0] + aff2_mvx[1])/2;
		aff2_mvy[2] = (aff2_mvy[0] + aff2_mvy[1])/2;
	}

	if( ( ( fabs(aff2_mvx[2] - aff2_mvx[0]) < SMALL_DIFF && fabs(aff2_mvy[2] - aff2_mvy[0]) < SMALL_DIFF) || ( fabs(aff2_mvx[2] - aff2_mvx[1]) < SMALL_DIFF && fabs(aff2_mvy[2] - aff2_mvy[1]) < SMALL_DIFF) )
		&& (aff2_mvx[2] != (float)HUGE_VAL && aff2_mvy[2] != (float)HUGE_VAL) ){
		aff2_mvx[2] = (float)HUGE_VAL;
		aff2_mvy[2] = (float)HUGE_VAL;
	}

//TEMP
	for( i = 0; i < 3 ; i ++ ){
		if( aff2_mvx[i] == (float)HUGE_VAL && aff2_mvy[i] == (float)HUGE_VAL ){
			tempx = prev_frame_motion_field[y*hor + x + xblk].mvx;
			tempy = prev_frame_motion_field[y*hor + x + xblk].mvy;

			for( j = 0; j < 3 ; j ++ ){
				if( tempx != (float)HUGE_VAL && tempy != (float)HUGE_VAL && (fabs(tempx - aff2_mvx[j]) < SMALL_DIFF && fabs(tempy - aff2_mvy[j]) < SMALL_DIFF) ){
					tempx = (float)HUGE_VAL;
					tempy = (float)HUGE_VAL;
				}
			}
			aff2_mvx[i] = tempx;
			aff2_mvy[i] = tempy;
			break;
		}
	}

// MV3
	if(y+yblk-1<ver && x-1>=0){
		aff3_mvx[0] = cur_frame_motion_field[(y+yblk-1)*hor + (x-1)].mvx;
		aff3_mvy[0] = cur_frame_motion_field[(y+yblk-1)*hor + (x-1)].mvy;

		aff3_mvx[0] = cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk-3)*hor + (x-2)].mvx) + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk-2)*hor + (x-3)].mvx);
		aff3_mvy[0] = cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk-3)*hor + (x-2)].mvy) + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk-2)*hor + (x-3)].mvy);

		aff3_mvx[0] = aff3_mvx[0] * AFF_SUBPEL;
		aff3_mvy[0] = aff3_mvy[0] * AFF_SUBPEL;
		int_mvx = (int)(aff3_mvx[0]);
		int_mvy = (int)(aff3_mvy[0]);
		aff3_mvx[0] = (float)(int_mvx);
		aff3_mvy[0] = (float)(int_mvy);
		aff3_mvx[0] = aff3_mvx[0] / AFF_SUBPEL;
		aff3_mvy[0] = aff3_mvy[0] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff3_mvx[0] = (float)HUGE_VAL;
			aff3_mvy[0] = (float)HUGE_VAL;
		}

//		if(aff3_mvx[0]!=get3_mvx[0] || aff3_mvy[0]!=get3_mvy[0])
//			printf("aff3_mvx0 = %f, aff3_mvy0 = %f, get3_mvx0 = %f, get3_mvy0 = %f\n\n",aff3_mvx[0],aff3_mvy[0],get3_mvx[0],get3_mvy[0]);
	}
	if(y+yblk<ver && x-1>=0){
		aff3_mvx[1] = cur_frame_motion_field[(y+yblk)*hor + (x-1)].mvx;
		aff3_mvy[1] = cur_frame_motion_field[(y+yblk)*hor + (x-1)].mvy;

		aff3_mvx[1] = cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvx + 2*(cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk)*hor + (x-3)].mvx);
		aff3_mvy[1] = cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvy + 2*(cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk)*hor + (x-3)].mvy);
		
		aff3_mvx[1] = aff3_mvx[1] * AFF_SUBPEL;
		aff3_mvy[1] = aff3_mvy[1] * AFF_SUBPEL;
		int_mvx = (int)(aff3_mvx[1]);
		int_mvy = (int)(aff3_mvy[1]);
		aff3_mvx[1] = (float)(int_mvx);
		aff3_mvy[1] = (float)(int_mvy);
		aff3_mvx[1] = aff3_mvx[1] / AFF_SUBPEL;
		aff3_mvy[1] = aff3_mvy[1] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff3_mvx[1] = (float)HUGE_VAL;
			aff3_mvy[1] = (float)HUGE_VAL;
		}

//		if(aff3_mvx[1]!=get3_mvx[1] || aff3_mvy[1]!=get3_mvy[1])
//			printf("aff3_mvx1 = %f, aff3_mvy1 = %f, get3_mvx1 = %f, get3_mvy1 = %f\n\n",aff3_mvx[1],aff3_mvy[1],get3_mvx[1],get3_mvy[1]);
	}

	if( fabs(aff3_mvx[1] - aff3_mvx[0]) < SMALL_DIFF && fabs(aff3_mvy[1] - aff3_mvy[0]) < SMALL_DIFF && aff3_mvx[1] != (float)HUGE_VAL && aff3_mvy[1] != (float)HUGE_VAL){
		aff3_mvx[1] = (float)HUGE_VAL;
		aff3_mvy[1] = (float)HUGE_VAL;
	}

	if(aff3_mvx[0] != (float)HUGE_VAL && aff3_mvy[0] != (float)HUGE_VAL && aff3_mvx[1] != (float)HUGE_VAL && aff3_mvy[1] != (float)HUGE_VAL){
		aff3_mvx[2] = (aff3_mvx[0] + aff3_mvx[1])/2;
		aff3_mvy[2] = (aff3_mvy[0] + aff3_mvy[1])/2;
	}

	if( ( ( fabs(aff3_mvx[2] - aff3_mvx[0]) < SMALL_DIFF && fabs(aff3_mvy[2] - aff3_mvy[0]) < SMALL_DIFF) || ( fabs(aff3_mvx[2] - aff3_mvx[1]) < SMALL_DIFF && fabs(aff3_mvy[2] - aff3_mvy[1]) < SMALL_DIFF) )
		&& (aff3_mvx[2] != (float)HUGE_VAL && aff3_mvy[2] != (float)HUGE_VAL) ){
		aff3_mvx[2] = (float)HUGE_VAL;
		aff3_mvy[2] = (float)HUGE_VAL;
	}

//TEMP
	for( i = 0; i < 3 ; i ++ ){
		if( aff3_mvx[i] == (float)HUGE_VAL && aff3_mvy[i] == (float)HUGE_VAL ){
			tempx = prev_frame_motion_field[(y+yblk)*hor + x].mvx;
			tempy = prev_frame_motion_field[(y+yblk)*hor + x].mvy;

			for( j = 0; j < 3 ; j ++ ){
				if( tempx != (float)HUGE_VAL && tempy != (float)HUGE_VAL && (fabs(tempx - aff3_mvx[j]) < SMALL_DIFF && fabs(tempy - aff3_mvy[j]) < SMALL_DIFF) ){
					tempx = (float)HUGE_VAL;
					tempy = (float)HUGE_VAL;
				}
			}
			aff3_mvx[i] = tempx;
			aff3_mvy[i] = tempy;
			break;
		}
	}

/////////////////////////////

	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
		for(k=0;k<3;k++){
			if( aff1_mvx[i]!=(float)HUGE_VAL && aff2_mvx[j]!=(float)HUGE_VAL && aff3_mvx[k]!=(float)HUGE_VAL ){
				assert(aff1_mvy[i]!=(float)HUGE_VAL && aff2_mvy[j]!=(float)HUGE_VAL && aff3_mvy[k]!=(float)HUGE_VAL);
				aff_num ++;
			}
		}
	  }
    }

	for(i=0;i<2;i++){
		if( aff2_mvx[i]!=(float)HUGE_VAL){
			assert(aff2_mvy[i]!=(float)HUGE_VAL);
				left_num ++;
		}
	}
	for(i=0;i<2;i++){
		if( aff3_mvx[i]!=(float)HUGE_VAL){
			assert(aff3_mvy[i]!=(float)HUGE_VAL);
				up_num ++;
		}
	}

	if(enc_idx == 1){
		assert(left_num >= 1 && up_num >= 1 && aff_num >= 1);
		printf("aff_num = %d\n\n",aff_num);
	}

	printf("Predictors:\n");
	for(i = 0; i < 3; i ++ ){
		printf("aff1_mvx[%d] = %f, aff1_mvy[%d] = %f\naff2_mvx[%d] = %f, aff2_mvy[%d] = %f\naff3_mvx[%d] = %f, aff3_mvy[%d] = %f\n",
			i,aff1_mvx[i],i,aff1_mvy[i],i,aff2_mvx[i],i,aff2_mvy[i],i,aff3_mvx[i],i,aff3_mvy[i]);
	}
	printf("\n\n");

	for(i = 0; i < 3; i ++ ){
		if( !(aff1_mvx[i] == fmv->aff1_pred_mvx[i] && aff1_mvy[i] == fmv->aff1_pred_mvy[i] &&
			aff2_mvx[i] == fmv->aff2_pred_mvx[i] && aff2_mvy[i] == fmv->aff2_pred_mvy[i] &&
			aff3_mvx[i] == fmv->aff3_pred_mvx[i] && aff3_mvy[i] == fmv->aff3_pred_mvy[i] ) ){
			printf("Aff Pred Error!\nx = %d, y = %d, xblk = %d, yblk = %d\n",x,y,xblk,yblk);
			for(k = 0; k < 3; k ++ ){
				printf("aff1_mvx[%d] = %f, aff1_mvy[%d] = %f\naff2_mvx[%d] = %f, aff2_mvy[%d] = %f\naff3_mvx[%d] = %f, aff3_mvy[%d] = %f\n",
					k,fmv->aff1_pred_mvx[k],k,fmv->aff1_pred_mvy[k],k,fmv->aff2_pred_mvx[k],k,fmv->aff2_pred_mvy[k],k,fmv->aff3_pred_mvx[k],k,
					fmv->aff3_pred_mvy[k]);
			}
			assert(0);
		}

	}

	if(enc_idx == 1){
		if( fmv->direct_idx == DIRECT || (fmv->direct_idx == INDIRECT && fmv->merge_idx == INTER) ){
			
			if(aff_num >= 9){
				getnum = fmv->aff_idx;
				while(getnum >= 4){
					getnum -= 4;
					putbits(1,1);
					length += 1;
				}
				putbits(0,1);
				length += 1;
				assert(getnum <= 3 && getnum >= 0);
				putbits(getnum,2);
				length += 2;
			}else if(aff_num >= 5 && aff_num <= 8){
				putbits(fmv->aff_idx,3);
				length += 3;
			}else if(aff_num >= 3 && aff_num <= 4){
				putbits(fmv->aff_idx,2);
				length += 2;
			}else if(aff_num == 2){
				putbits(fmv->aff_idx,1);
				length += 1;
			}else
				assert(fmv->aff_idx == 0);

		}else{
			assert(fmv->direct_idx == INDIRECT && fmv->merge_idx == MERGE);

			if(fmv->merge_dir == UP){ //Merge to the left
				if(up_num == 2){
					putbits(fmv->aff_idx,1);  //1-bit index for affine matrix prediction
					length+=1;
				}else
					assert(fmv->aff_idx == 0);
			}
			else if(fmv->merge_dir == LEFT){ //Merge to the up
				if(left_num == 2){
					putbits(fmv->aff_idx,1);  //1-bit index for affine matrix prediction
					length+=1;
				}else
					assert(fmv->aff_idx == 0);
			}
			else if(fmv->merge_dir == PAL_L){ //Must be a parallel affine mode
				getnum = fmv->aff_idx;

				if(aff_num >= 9){
					while(getnum >= 4){
						getnum -= 4;
						putbits(1,1);
						length += 1;
					}
					putbits(0,1);
					length += 1;
					assert(getnum <= 3 && getnum >= 0);
					putbits(getnum,2);
					length += 2;
				}else if(aff_num >= 5 && aff_num <= 8){
					putbits(fmv->aff_idx,3);
					length += 3;
				}else if(aff_num >= 3 && aff_num <= 4){
					putbits(fmv->aff_idx,2);
					length += 2;
				}
				else if(aff_num == 2){
					putbits(fmv->aff_idx,1);
					length += 1;
				}else
					assert(fmv->aff_idx == 0);
			}else
				assert(fmv->merge_dir == TRAN_P);  //the translational predictor index has been transmitted instead

		}
	}//if enc_idx

	return length;
}


void get_aff_mvs(vector_ptr fmv, int x, int y, int xblk, int yblk, int hor, int ver, FRAME_MOTION_FIELD* cur_frame_motion_field, SIMP_FRAME_MOTION_FIELD* prev_frame_motion_field, int dec_idx){
	int getnum, count = 0, length, aff_num = 0, up_num = 0, left_num = 0;
	int i,j,k;
	float aff1_mvx[3],aff1_mvy[3],aff2_mvx[3],aff2_mvy[3],aff3_mvx[3],aff3_mvy[3];
	float get1_mvx[3],get1_mvy[3],get2_mvx[3],get2_mvy[3],get3_mvx[3],get3_mvy[3];
	float tempx,tempy;
	int int_mvx, int_mvy;

	int val;
	float sub = 0.25;

	int subpel = (int)(1/sub);

	for(i=0;i<3;i++){
		aff1_mvx[i] = (float)HUGE_VAL;get1_mvx[i] = (float)HUGE_VAL;
		aff1_mvy[i] = (float)HUGE_VAL;get1_mvy[i] = (float)HUGE_VAL;
		aff2_mvx[i] = (float)HUGE_VAL;get2_mvx[i] = (float)HUGE_VAL;
		aff2_mvy[i] = (float)HUGE_VAL;get2_mvy[i] = (float)HUGE_VAL;
		aff3_mvx[i] = (float)HUGE_VAL;get3_mvx[i] = (float)HUGE_VAL;
		aff3_mvy[i] = (float)HUGE_VAL;get3_mvy[i] = (float)HUGE_VAL;
	}

// MV1
	if(y-1>=0){
		aff1_mvx[0] = cur_frame_motion_field[(y-1)*hor + x].mvx;
		aff1_mvy[0] = cur_frame_motion_field[(y-1)*hor + x].mvy;

		aff1_mvx[0] = cur_frame_motion_field[(y-2)*hor + x].mvx + 2 * (cur_frame_motion_field[(y-2)*hor + x].mvx - cur_frame_motion_field[(y-3)*hor + x].mvx);
		aff1_mvy[0] = cur_frame_motion_field[(y-2)*hor + x].mvy + 2 * (cur_frame_motion_field[(y-2)*hor + x].mvy - cur_frame_motion_field[(y-3)*hor + x].mvy);

		aff1_mvx[0] = aff1_mvx[0] * AFF_SUBPEL;
		aff1_mvy[0] = aff1_mvy[0] * AFF_SUBPEL;
		int_mvx = (int)(aff1_mvx[0]);
		int_mvy = (int)(aff1_mvy[0]);
		aff1_mvx[0] = (float)(int_mvx);
		aff1_mvy[0] = (float)(int_mvy);
		aff1_mvx[0] = aff1_mvx[0] / AFF_SUBPEL;
		aff1_mvy[0] = aff1_mvy[0] / AFF_SUBPEL;
		
		if(cur_frame_motion_field[(y-2)*hor + x].mvx == (float)HUGE_VAL){
			aff1_mvx[0] = (float)HUGE_VAL;
			aff1_mvy[0] = (float)HUGE_VAL;
		}
	}

	if(x-1>=0){
		aff1_mvx[1] = cur_frame_motion_field[y*hor + (x-1)].mvx;
		aff1_mvy[1] = cur_frame_motion_field[y*hor + (x-1)].mvy;

		aff1_mvx[1] = cur_frame_motion_field[y*hor + (x-2)].mvx + 2 * (cur_frame_motion_field[y*hor + (x-2)].mvx - cur_frame_motion_field[y*hor + (x-3)].mvx);
		aff1_mvy[1] = cur_frame_motion_field[y*hor + (x-2)].mvy + 2 * (cur_frame_motion_field[y*hor + (x-2)].mvy - cur_frame_motion_field[y*hor + (x-3)].mvy);

		aff1_mvx[1] = aff1_mvx[1] * AFF_SUBPEL;
		aff1_mvy[1] = aff1_mvy[1] * AFF_SUBPEL;
		int_mvx = (int)(aff1_mvx[1]);
		int_mvy = (int)(aff1_mvy[1]);
		aff1_mvx[1] = (float)(int_mvx);
		aff1_mvy[1] = (float)(int_mvy);
		aff1_mvx[1] = aff1_mvx[1] / AFF_SUBPEL;
		aff1_mvy[1] = aff1_mvy[1] / AFF_SUBPEL;
		
		if(cur_frame_motion_field[y*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff1_mvx[1] = (float)HUGE_VAL;
			aff1_mvy[1] = (float)HUGE_VAL;
		}
	}

	if(y-1>=0 && x-1>=0){
		aff1_mvx[2] = cur_frame_motion_field[(y-1)*hor + (x-1)].mvx;
		aff1_mvy[2] = cur_frame_motion_field[(y-1)*hor + (x-1)].mvy;

		aff1_mvx[2] = cur_frame_motion_field[(y-2)*hor + (x-2)].mvx + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y-3)*hor + (x-2)].mvx) + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y-2)*hor + (x-3)].mvx);
		aff1_mvy[2] = cur_frame_motion_field[(y-2)*hor + (x-2)].mvy + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y-3)*hor + (x-2)].mvy) + 2*(cur_frame_motion_field[(y-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y-2)*hor + (x-3)].mvy);
		
		aff1_mvx[2] = aff1_mvx[2] * AFF_SUBPEL;
		aff1_mvy[2] = aff1_mvy[2] * AFF_SUBPEL;
		int_mvx = (int)(aff1_mvx[2]);
		int_mvy = (int)(aff1_mvy[2]);
		aff1_mvx[2] = (float)(int_mvx);
		aff1_mvy[2] = (float)(int_mvy);
		aff1_mvx[2] = aff1_mvx[2] / AFF_SUBPEL;
		aff1_mvy[2] = aff1_mvy[2] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y-2)*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff1_mvx[2] = (float)HUGE_VAL;
			aff1_mvy[2] = (float)HUGE_VAL;
		}
	}

	if( fabs(aff1_mvx[1] - aff1_mvx[0]) < SMALL_DIFF && fabs(aff1_mvy[1] - aff1_mvy[0]) < SMALL_DIFF && aff1_mvx[1] != (float)HUGE_VAL && aff1_mvy[1] != (float)HUGE_VAL){
		aff1_mvx[1] = (float)HUGE_VAL;
		aff1_mvy[1] = (float)HUGE_VAL;
	}
	if( ( ( fabs(aff1_mvx[2] - aff1_mvx[0]) < SMALL_DIFF && fabs(aff1_mvy[2] - aff1_mvy[0]) < SMALL_DIFF) || ( fabs(aff1_mvx[2] - aff1_mvx[1]) < SMALL_DIFF && fabs(aff1_mvy[2] - aff1_mvy[1]) < SMALL_DIFF) )
		&& (aff1_mvx[2] != (float)HUGE_VAL && aff1_mvy[2] != (float)HUGE_VAL) ){
		aff1_mvx[2] = (float)HUGE_VAL;
		aff1_mvy[2] = (float)HUGE_VAL;
	}

//TEMP
	for( i = 0; i < 3 ; i ++ ){
		if( aff1_mvx[i] == (float)HUGE_VAL && aff1_mvy[i] == (float)HUGE_VAL ){
			tempx = prev_frame_motion_field[y*hor + x].mvx;
			tempy = prev_frame_motion_field[y*hor + x].mvy;

			for( j = 0; j < 3 ; j ++ ){
				if( tempx != (float)HUGE_VAL && tempy != (float)HUGE_VAL && (fabs(tempx - aff1_mvx[j]) < SMALL_DIFF && fabs(tempy - aff1_mvy[j]) < SMALL_DIFF) ){
					tempx = (float)HUGE_VAL;
					tempy = (float)HUGE_VAL;
				}
			}
			aff1_mvx[i] = tempx;
			aff1_mvy[i] = tempy;
			break;
		}
	}

// MV2
	if(y-1>=0 && x+xblk-1<hor){
		aff2_mvx[0] = cur_frame_motion_field[(y-1)*hor + (x+xblk-1)].mvx;
		aff2_mvy[0] = cur_frame_motion_field[(y-1)*hor + (x+xblk-1)].mvy;

		aff2_mvx[0] = cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx - cur_frame_motion_field[(y-3)*hor + (x+xblk-2)].mvx) + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx - cur_frame_motion_field[(y-2)*hor + (x+xblk-3)].mvx);
		aff2_mvy[0] = cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy - cur_frame_motion_field[(y-3)*hor + (x+xblk-2)].mvy) + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy - cur_frame_motion_field[(y-2)*hor + (x+xblk-3)].mvy);
		
		aff2_mvx[0] = aff2_mvx[0] * AFF_SUBPEL;
		aff2_mvy[0] = aff2_mvy[0] * AFF_SUBPEL;
		int_mvx = (int)(aff2_mvx[0]);
		int_mvy = (int)(aff2_mvy[0]);
		aff2_mvx[0] = (float)(int_mvx);
		aff2_mvy[0] = (float)(int_mvy);
		aff2_mvx[0] = aff2_mvx[0] / AFF_SUBPEL;
		aff2_mvy[0] = aff2_mvy[0] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx == (float)HUGE_VAL){
			aff2_mvx[0] = (float)HUGE_VAL;
			aff2_mvy[0] = (float)HUGE_VAL;
		}
	}
	if(y-1>=0 && x+xblk<hor){
		aff2_mvx[1] = cur_frame_motion_field[(y-1)*hor + (x+xblk)].mvx;
		aff2_mvy[1] = cur_frame_motion_field[(y-1)*hor + (x+xblk)].mvy;

		aff2_mvx[1] = cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvx + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvx - cur_frame_motion_field[(y-3)*hor + (x+xblk)].mvx);
		aff2_mvy[1] = cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvy + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvy - cur_frame_motion_field[(y-3)*hor + (x+xblk)].mvy);

		aff2_mvx[1] = aff2_mvx[1] * AFF_SUBPEL;
		aff2_mvy[1] = aff2_mvy[1] * AFF_SUBPEL;
		int_mvx = (int)(aff2_mvx[1]);
		int_mvy = (int)(aff2_mvy[1]);
		aff2_mvx[1] = (float)(int_mvx);
		aff2_mvy[1] = (float)(int_mvy);
		aff2_mvx[1] = aff2_mvx[1] / AFF_SUBPEL;
		aff2_mvy[1] = aff2_mvy[1] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y-2)*hor + (x+xblk)].mvx == (float)HUGE_VAL){
			aff2_mvx[1] = (float)HUGE_VAL;
			aff2_mvy[1] = (float)HUGE_VAL;
		}
	}

	if( fabs(aff2_mvx[1] - aff2_mvx[0]) < SMALL_DIFF && fabs(aff2_mvy[1] - aff2_mvy[0]) < SMALL_DIFF && aff2_mvx[1] != (float)HUGE_VAL && aff2_mvy[1] != (float)HUGE_VAL){
		aff2_mvx[1] = (float)HUGE_VAL;
		aff2_mvy[1] = (float)HUGE_VAL;
	}

	if(aff2_mvx[0] != (float)HUGE_VAL && aff2_mvy[0] != (float)HUGE_VAL && aff2_mvx[1] != (float)HUGE_VAL && aff2_mvy[1] != (float)HUGE_VAL){
		aff2_mvx[2] = (aff2_mvx[0] + aff2_mvx[1])/2;
		aff2_mvy[2] = (aff2_mvy[0] + aff2_mvy[1])/2;
	}

	if( ( ( fabs(aff2_mvx[2] - aff2_mvx[0]) < SMALL_DIFF && fabs(aff2_mvy[2] - aff2_mvy[0]) < SMALL_DIFF) || ( fabs(aff2_mvx[2] - aff2_mvx[1]) < SMALL_DIFF && fabs(aff2_mvy[2] - aff2_mvy[1]) < SMALL_DIFF) )
		&& (aff2_mvx[2] != (float)HUGE_VAL && aff2_mvy[2] != (float)HUGE_VAL) ){
		aff2_mvx[2] = (float)HUGE_VAL;
		aff2_mvy[2] = (float)HUGE_VAL;
	}

//TEMP
	for( i = 0; i < 3 ; i ++ ){
		if( aff2_mvx[i] == (float)HUGE_VAL && aff2_mvy[i] == (float)HUGE_VAL ){
			tempx = prev_frame_motion_field[y*hor + x + xblk].mvx;
			tempy = prev_frame_motion_field[y*hor + x + xblk].mvy;

			for( j = 0; j < 3 ; j ++ ){
				if( tempx != (float)HUGE_VAL && tempy != (float)HUGE_VAL && (fabs(tempx - aff2_mvx[j]) < SMALL_DIFF && fabs(tempy - aff2_mvy[j]) < SMALL_DIFF) ){
					tempx = (float)HUGE_VAL;
					tempy = (float)HUGE_VAL;
				}
			}
			aff2_mvx[i] = tempx;
			aff2_mvy[i] = tempy;
			break;
		}
	}

// MV3
	if(y+yblk-1<ver && x-1>=0){
		aff3_mvx[0] = cur_frame_motion_field[(y+yblk-1)*hor + (x-1)].mvx;
		aff3_mvy[0] = cur_frame_motion_field[(y+yblk-1)*hor + (x-1)].mvy;

		aff3_mvx[0] = cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk-3)*hor + (x-2)].mvx) + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk-2)*hor + (x-3)].mvx);
		aff3_mvy[0] = cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk-3)*hor + (x-2)].mvy) + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk-2)*hor + (x-3)].mvy);

		aff3_mvx[0] = aff3_mvx[0] * AFF_SUBPEL;
		aff3_mvy[0] = aff3_mvy[0] * AFF_SUBPEL;
		int_mvx = (int)(aff3_mvx[0]);
		int_mvy = (int)(aff3_mvy[0]);
		aff3_mvx[0] = (float)(int_mvx);
		aff3_mvy[0] = (float)(int_mvy);
		aff3_mvx[0] = aff3_mvx[0] / AFF_SUBPEL;
		aff3_mvy[0] = aff3_mvy[0] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff3_mvx[0] = (float)HUGE_VAL;
			aff3_mvy[0] = (float)HUGE_VAL;
		}
	}
	if(y+yblk<ver && x-1>=0){
		aff3_mvx[1] = cur_frame_motion_field[(y+yblk)*hor + (x-1)].mvx;
		aff3_mvy[1] = cur_frame_motion_field[(y+yblk)*hor + (x-1)].mvy;

		aff3_mvx[1] = cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvx + 2*(cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk)*hor + (x-3)].mvx);
		aff3_mvy[1] = cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvy + 2*(cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk)*hor + (x-3)].mvy);
		
		aff3_mvx[1] = aff3_mvx[1] * AFF_SUBPEL;
		aff3_mvy[1] = aff3_mvy[1] * AFF_SUBPEL;
		int_mvx = (int)(aff3_mvx[1]);
		int_mvy = (int)(aff3_mvy[1]);
		aff3_mvx[1] = (float)(int_mvx);
		aff3_mvy[1] = (float)(int_mvy);
		aff3_mvx[1] = aff3_mvx[1] / AFF_SUBPEL;
		aff3_mvy[1] = aff3_mvy[1] / AFF_SUBPEL;

		if(cur_frame_motion_field[(y+yblk)*hor + (x-2)].mvx == (float)HUGE_VAL){
			aff3_mvx[1] = (float)HUGE_VAL;
			aff3_mvy[1] = (float)HUGE_VAL;
		}
	}

	if( fabs(aff3_mvx[1] - aff3_mvx[0]) < SMALL_DIFF && fabs(aff3_mvy[1] - aff3_mvy[0]) < SMALL_DIFF && aff3_mvx[1] != (float)HUGE_VAL && aff3_mvy[1] != (float)HUGE_VAL){
		aff3_mvx[1] = (float)HUGE_VAL;
		aff3_mvy[1] = (float)HUGE_VAL;
	}

	if(aff3_mvx[0] != (float)HUGE_VAL && aff3_mvy[0] != (float)HUGE_VAL && aff3_mvx[1] != (float)HUGE_VAL && aff3_mvy[1] != (float)HUGE_VAL){
		aff3_mvx[2] = (aff3_mvx[0] + aff3_mvx[1])/2;
		aff3_mvy[2] = (aff3_mvy[0] + aff3_mvy[1])/2;
	}

	if( ( ( fabs(aff3_mvx[2] - aff3_mvx[0]) < SMALL_DIFF && fabs(aff3_mvy[2] - aff3_mvy[0]) < SMALL_DIFF) || ( fabs(aff3_mvx[2] - aff3_mvx[1]) < SMALL_DIFF && fabs(aff3_mvy[2] - aff3_mvy[1]) < SMALL_DIFF) )
		&& (aff3_mvx[2] != (float)HUGE_VAL && aff3_mvy[2] != (float)HUGE_VAL) ){
		aff3_mvx[2] = (float)HUGE_VAL;
		aff3_mvy[2] = (float)HUGE_VAL;
	}

//TEMP
	for( i = 0; i < 3 ; i ++ ){
		if( aff3_mvx[i] == (float)HUGE_VAL && aff3_mvy[i] == (float)HUGE_VAL ){
			tempx = prev_frame_motion_field[(y+yblk)*hor + x].mvx;
			tempy = prev_frame_motion_field[(y+yblk)*hor + x].mvy;

			for( j = 0; j < 3 ; j ++ ){
				if( tempx != (float)HUGE_VAL && tempy != (float)HUGE_VAL && (fabs(tempx - aff3_mvx[j]) < SMALL_DIFF && fabs(tempy - aff3_mvy[j]) < SMALL_DIFF) ){
					tempx = (float)HUGE_VAL;
					tempy = (float)HUGE_VAL;
				}
			}
			aff3_mvx[i] = tempx;
			aff3_mvy[i] = tempy;
			break;
		}
	}

/////////////////////////////

	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
		for(k=0;k<3;k++){
			if( aff1_mvx[i]!=(float)HUGE_VAL && aff2_mvx[j]!=(float)HUGE_VAL && aff3_mvx[k]!=(float)HUGE_VAL ){
				assert(aff1_mvy[i]!=(float)HUGE_VAL && aff2_mvy[j]!=(float)HUGE_VAL && aff3_mvy[k]!=(float)HUGE_VAL);
				aff_num ++;
			}
		}
	  }
    }

	for(i=0;i<2;i++){
		if( aff2_mvx[i]!=(float)HUGE_VAL){
			assert(aff2_mvy[i]!=(float)HUGE_VAL);
				left_num ++;
		}
	}
	for(i=0;i<2;i++){
		if( aff3_mvx[i]!=(float)HUGE_VAL){
			assert(aff3_mvy[i]!=(float)HUGE_VAL);
				up_num ++;
		}
	}

	if(dec_idx == 1){

		if( !(fmv->direct_idx == INDIRECT && fmv->merge_idx == MERGE && fmv->merge_dir == TRAN_P) )
			assert(left_num >= 1 && up_num >= 1 && aff_num >= 1);
		
		printf("aff_num = %d\n\n",aff_num);
	}

	printf("Predictors:\n");
	for(i = 0; i < 3; i ++ ){
		printf("aff1_mvx[%d] = %f, aff1_mvy[%d] = %f\naff2_mvx[%d] = %f, aff2_mvy[%d] = %f\naff3_mvx[%d] = %f, aff3_mvy[%d] = %f\n",
			i,aff1_mvx[i],i,aff1_mvy[i],i,aff2_mvx[i],i,aff2_mvy[i],i,aff3_mvx[i],i,aff3_mvy[i]);
	}
	printf("\n\n");
	////////////	Added on 09.10.16	/////////////////
	if(dec_idx == 1){
		if(fmv->direct_idx == DIRECT || (fmv->direct_idx == INDIRECT && fmv->merge_idx == INTER) ){
			
			if(aff_num >= 9){
				val = 0;
				getnum = getbits(1);
				while(getnum == 1){
					val += 4;
					getnum = getbits(1);
				}

				assert(getnum == 0);
				getnum = getbits(2);
				val += getnum;
				fmv->aff_idx = val;

			}else if(aff_num >= 5 && aff_num <= 8){
				val = getbits(3);
				fmv->aff_idx = val;
			}else if(aff_num >= 3 && aff_num <= 4){
				val = getbits(2);
				fmv->aff_idx = val;
			}else if(aff_num == 2){
				val = getbits(1);
				fmv->aff_idx = val;
			}else
				fmv->aff_idx = 0;

		}else{
			assert(fmv->direct_idx == INDIRECT && fmv->merge_idx == MERGE);

			if(fmv->merge_dir == UP){
				if(up_num == 2){
					val = getbits(1);
					fmv->aff_idx = val;
				}else
					fmv->aff_idx = 0;
			}else if(fmv->merge_dir == LEFT){
				if(left_num == 2){
					val = getbits(1);
					fmv->aff_idx = val;
				}else
					fmv->aff_idx = 0;
			}else if(fmv->merge_dir == PAL_L){

				if(aff_num >= 9){
					val = 0;
					getnum = getbits(1);
					while(getnum == 1){
						val += 4;
						getnum = getbits(1);
					}

					assert(getnum == 0);
					getnum = getbits(2);
					val += getnum;
					fmv->aff_idx = val;

				}else if(aff_num >= 5 && aff_num <= 8){
					val = getbits(3);
					fmv->aff_idx = val;
				}else if(aff_num >= 3 && aff_num <= 4){
					val = getbits(2);
					fmv->aff_idx = val;
				}else if(aff_num == 2){
					val = getbits(1);
					fmv->aff_idx = val;
				}else
					fmv->aff_idx = 0;
			}else
				assert(fmv->merge_dir == TRAN_P);  //predictors are left/right translational MVs that have been decoded.
		}
	}//if dec_idx
	else{
		fmv->aff_idx = -1;
	}
	/////////////////////////////////////////////////////

	assert(fmv->aff_idx == -1 || (fmv->aff_idx >= 0 && fmv->aff_idx <= 26) );

	if(fmv->aff_idx >= 0)
		aff_cum_idx[fmv->aff_idx]++;

	if(fmv->direct_idx == DIRECT){

//		printf("fmv->aff_idx = %d\n",fmv->aff_idx);

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					if(aff1_mvx[i]!=(float)HUGE_VAL && aff1_mvy[i]!=(float)HUGE_VAL && aff2_mvx[j]!=(float)HUGE_VAL && aff2_mvy[j]!=(float)HUGE_VAL &&
						aff3_mvx[k]!=(float)HUGE_VAL && aff3_mvy[k]!=(float)HUGE_VAL){
						if(count == fmv->aff_idx){
							fmv->aff1_mvx = aff1_mvx[i];
							fmv->aff1_mvy = aff1_mvy[i];
							fmv->aff2_mvx = aff2_mvx[j];
							fmv->aff2_mvy = aff2_mvy[j];
							fmv->aff3_mvx = aff3_mvx[k];
							fmv->aff3_mvy = aff3_mvy[k];
							goto sign1;
						}
						count++;
					}
				}
			}
		}
sign1:		;
	}else if(fmv->direct_idx == INDIRECT){
		count = 0;

		if(fmv->merge_idx == MERGE){

			if(fmv->merge_dir == UP){
				assert(fmv->aff_idx == 0 || fmv->aff_idx == 1);
				fmv->aff1_mvx = cur_frame_motion_field[(y-1)*hor + x].mvx;
				fmv->aff1_mvy = cur_frame_motion_field[(y-1)*hor + x].mvy;

				fmv->aff1_mvx = cur_frame_motion_field[(y-2)*hor + x].mvx + 2*(cur_frame_motion_field[(y-2)*hor + x].mvx - cur_frame_motion_field[(y-3)*hor + x].mvx);
				fmv->aff1_mvy = cur_frame_motion_field[(y-2)*hor + x].mvy + 2*(cur_frame_motion_field[(y-2)*hor + x].mvy - cur_frame_motion_field[(y-3)*hor + x].mvy);

				fmv->aff1_mvx = fmv->aff1_mvx * AFF_SUBPEL;
				fmv->aff1_mvy = fmv->aff1_mvy * AFF_SUBPEL;
				int_mvx = (int)(fmv->aff1_mvx);
				int_mvy = (int)(fmv->aff1_mvy);
				fmv->aff1_mvx = (float)(int_mvx);
				fmv->aff1_mvy = (float)(int_mvy);
				fmv->aff1_mvx = fmv->aff1_mvx / AFF_SUBPEL;
				fmv->aff1_mvy = fmv->aff1_mvy / AFF_SUBPEL;

				if(cur_frame_motion_field[(y-2)*hor + x].mvx == (float)HUGE_VAL){
					fmv->aff1_mvx = (float)HUGE_VAL;
					fmv->aff1_mvy = (float)HUGE_VAL;
				}

				//////////////////////////////
				fmv->aff2_mvx = cur_frame_motion_field[(y-1)*hor + (x+xblk-1)].mvx;
				fmv->aff2_mvy = cur_frame_motion_field[(y-1)*hor + (x+xblk-1)].mvy;

				fmv->aff2_mvx = cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx - cur_frame_motion_field[(y-3)*hor + (x+xblk-2)].mvx) + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx - cur_frame_motion_field[(y-2)*hor + (x+xblk-3)].mvx);
				fmv->aff2_mvy = cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy - cur_frame_motion_field[(y-3)*hor + (x+xblk-2)].mvy) + 2*(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvy - cur_frame_motion_field[(y-2)*hor + (x+xblk-3)].mvy);
		
				fmv->aff2_mvx = fmv->aff2_mvx * AFF_SUBPEL;
				fmv->aff2_mvy = fmv->aff2_mvy * AFF_SUBPEL;
				int_mvx = (int)(fmv->aff2_mvx);
				int_mvy = (int)(fmv->aff2_mvy);
				fmv->aff2_mvx = (float)(int_mvx);
				fmv->aff2_mvy = (float)(int_mvy);
				fmv->aff2_mvx = fmv->aff2_mvx / AFF_SUBPEL;
				fmv->aff2_mvy = fmv->aff2_mvy / AFF_SUBPEL;

				if(cur_frame_motion_field[(y-2)*hor + (x+xblk-2)].mvx == (float)HUGE_VAL){
					fmv->aff2_mvx = (float)HUGE_VAL;
					fmv->aff2_mvy = (float)HUGE_VAL;
				}

				for(i=0;i<2;i++){
					if( aff3_mvx[i] != (float)HUGE_VAL ){
						assert( aff3_mvy[i] != (float)HUGE_VAL );

						if(fmv->aff_idx == count){
							fmv->aff3_mvx = aff3_mvx[i] + fmv->aff3_dmvx;
							fmv->aff3_mvy = aff3_mvy[i] + fmv->aff3_dmvy;
							break;
						}
						count++;
					}
				}

			}
			else if(fmv->merge_dir == LEFT){
				assert(fmv->aff_idx == 0 || fmv->aff_idx == 1);
				fmv->aff1_mvx = cur_frame_motion_field[y*hor + (x-1)].mvx;
				fmv->aff1_mvy = cur_frame_motion_field[y*hor + (x-1)].mvy;

				fmv->aff1_mvx = cur_frame_motion_field[y*hor + (x-2)].mvx + 2 * (cur_frame_motion_field[y*hor + (x-2)].mvx - cur_frame_motion_field[y*hor + (x-3)].mvx);
				fmv->aff1_mvy = cur_frame_motion_field[y*hor + (x-2)].mvy + 2 * (cur_frame_motion_field[y*hor + (x-2)].mvy - cur_frame_motion_field[y*hor + (x-3)].mvy);

				fmv->aff1_mvx = fmv->aff1_mvx * AFF_SUBPEL;
				fmv->aff1_mvy = fmv->aff1_mvy * AFF_SUBPEL;
				int_mvx = (int)(fmv->aff1_mvx);
				int_mvy = (int)(fmv->aff1_mvy);
				fmv->aff1_mvx = (float)(int_mvx);
				fmv->aff1_mvy = (float)(int_mvy);
				fmv->aff1_mvx = fmv->aff1_mvx / AFF_SUBPEL;
				fmv->aff1_mvy = fmv->aff1_mvy / AFF_SUBPEL;

				if(cur_frame_motion_field[y*hor + (x-2)].mvx == (float)HUGE_VAL){
					fmv->aff1_mvx = (float)HUGE_VAL;
					fmv->aff1_mvy = (float)HUGE_VAL;
				}

				for(i=0;i<2;i++){
					if( aff2_mvx[i] != (float)HUGE_VAL ){
						assert( aff2_mvy[i] != (float)HUGE_VAL );

						if(fmv->aff_idx == count){
							fmv->aff2_mvx = aff2_mvx[i] + fmv->aff2_dmvx;
							fmv->aff2_mvy = aff2_mvy[i] + fmv->aff2_dmvy;
							break;
						}
						count++;
					}
				}

				fmv->aff3_mvx = cur_frame_motion_field[(y+yblk-1)*hor + (x-1)].mvx;
				fmv->aff3_mvy = cur_frame_motion_field[(y+yblk-1)*hor + (x-1)].mvy;

				fmv->aff3_mvx = cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk-2)*hor + (x-3)].mvx) + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx - cur_frame_motion_field[(y+yblk-3)*hor + (x-2)].mvx);
				fmv->aff3_mvy = cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk-2)*hor + (x-3)].mvy) + 2*(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvy - cur_frame_motion_field[(y+yblk-3)*hor + (x-2)].mvy);

				fmv->aff3_mvx = fmv->aff3_mvx * AFF_SUBPEL;
				fmv->aff3_mvy = fmv->aff3_mvy * AFF_SUBPEL;
				int_mvx = (int)(fmv->aff3_mvx);
				int_mvy = (int)(fmv->aff3_mvy);
				fmv->aff3_mvx = (float)(int_mvx);
				fmv->aff3_mvy = (float)(int_mvy);
				fmv->aff3_mvx = fmv->aff3_mvx / AFF_SUBPEL;
				fmv->aff3_mvy = fmv->aff3_mvy / AFF_SUBPEL;

				if(cur_frame_motion_field[(y+yblk-2)*hor + (x-2)].mvx == (float)HUGE_VAL){
					fmv->aff3_mvx = (float)HUGE_VAL;
					fmv->aff3_mvy = (float)HUGE_VAL;
				}


			}else if(fmv->merge_dir == PAL_L){

				if(fmv->aff_idx >= 0){
					assert(fmv->aff1_dmvx != (float)HUGE_VAL && fmv->aff1_dmvy != (float)HUGE_VAL && fmv->aff2_dmvx != (float)HUGE_VAL &&
					fmv->aff2_dmvy != (float)HUGE_VAL && fmv->aff3_dmvx != (float)HUGE_VAL && fmv->aff3_dmvy != (float)HUGE_VAL);

					for(i=0;i<3;i++){
						for(j=0;j<3;j++){
							for(k=0;k<3;k++){
								if(aff1_mvx[i]!=(float)HUGE_VAL && aff1_mvy[i]!=(float)HUGE_VAL && aff2_mvx[j]!=(float)HUGE_VAL && aff2_mvy[j]!=(float)HUGE_VAL &&
									aff3_mvx[k]!=(float)HUGE_VAL && aff3_mvy[k]!=(float)HUGE_VAL){
									if(count == fmv->aff_idx){
										fmv->aff1_mvx = aff1_mvx[i] + fmv->aff1_dmvx;
										fmv->aff1_mvy = aff1_mvy[i] + fmv->aff1_dmvy;
										fmv->aff2_mvx = aff2_mvx[j] + fmv->aff2_dmvx;
										fmv->aff2_mvy = aff2_mvy[j] + fmv->aff2_dmvy;
										fmv->aff3_mvx = aff3_mvx[k] + fmv->aff3_dmvx;
										fmv->aff3_mvy = aff3_mvy[k] + fmv->aff3_dmvy;
										goto sign2;
									}
									count++;
								}
							}
						}
					}
sign2:				;
				}//if not parallel side
			}else{
				assert(fmv->merge_dir == TRAN_P);
				
				fmv->aff1_mvx = fmv->mvx + fmv->aff1_dmvx;
				fmv->aff1_mvy = fmv->mvy + fmv->aff1_dmvy;
				fmv->aff2_mvx = fmv->mvx + fmv->aff2_dmvx;
				fmv->aff2_mvy = fmv->mvy + fmv->aff2_dmvy;
				fmv->aff3_mvx = fmv->mvx + fmv->aff3_dmvx;
				fmv->aff3_mvy = fmv->mvy + fmv->aff3_dmvy;
			}
		}//if MERGE
		else{
			assert(fmv->merge_idx == INTER);
			//		printf("fmv->aff_idx = %d\n",fmv->aff_idx);

			assert(fmv->aff1_dmvx != (float)HUGE_VAL && fmv->aff1_dmvy != (float)HUGE_VAL && fmv->aff2_dmvx != (float)HUGE_VAL &&
				fmv->aff2_dmvy != (float)HUGE_VAL && fmv->aff3_dmvx != (float)HUGE_VAL && fmv->aff3_dmvy != (float)HUGE_VAL);

			for(i=0;i<3;i++){
				for(j=0;j<3;j++){
					for(k=0;k<3;k++){
						if(aff1_mvx[i]!=(float)HUGE_VAL && aff1_mvy[i]!=(float)HUGE_VAL && aff2_mvx[j]!=(float)HUGE_VAL && aff2_mvy[j]!=(float)HUGE_VAL &&
							aff3_mvx[k]!=(float)HUGE_VAL && aff3_mvy[k]!=(float)HUGE_VAL){
							if(count == fmv->aff_idx){
								fmv->aff1_mvx = aff1_mvx[i] + fmv->aff1_dmvx;
								fmv->aff1_mvy = aff1_mvy[i] + fmv->aff1_dmvy;
								fmv->aff2_mvx = aff2_mvx[j] + fmv->aff2_dmvx;
								fmv->aff2_mvy = aff2_mvy[j] + fmv->aff2_dmvy;
								fmv->aff3_mvx = aff3_mvx[k] + fmv->aff3_dmvx;
								fmv->aff3_mvy = aff3_mvy[k] + fmv->aff3_dmvy;
								goto sign3;
							}
							count++;
						}
					}
				}
			}
sign3:		;
		}
	}
}
////////////////////////////////

// the bit output for sub-symbols
void put_splitted_mvbits(char bit, char sub_index, int layer_index)
{
	(splitted_mv_byte_array[layer_index][sub_index]) <<=1;
	(splitted_mv_byte_array[layer_index][sub_index]) |= bit;
	if ( --(splitted_bit_num_array[layer_index][sub_index]) == 0 )
	{
		store_splitted_bytes_array[layer_index][sub_index][(splitted_mv_byte_num_array[layer_index][sub_index])++] 
			 = splitted_mv_byte_array[layer_index][sub_index];
		splitted_mv_byte_array[layer_index][sub_index] = 0;
		splitted_bit_num_array[layer_index][sub_index] = 8;
	}
}


// the bit output for addition sign
void put_splitted_signbits(char bit, int layer_index)
{
	(splitted_sign_byte[layer_index]) <<=1;
	splitted_sign_byte[layer_index] |= bit;
	if ( --(splitted_sign_bit_num[layer_index]) == 0 )
	{
		store_splitted_sign[layer_index][(splitted_sign_byte_num[layer_index])++] = 
						splitted_sign_byte[layer_index];
		splitted_sign_byte[layer_index]    = 0;
		splitted_sign_bit_num[layer_index] = 8;
	}
}

void output_huff_bits(vector_ptr fmv1, int num, int *cnt, videoinfo info, int *mapbits){
	int i,j,k;
	int max, min1, min2;
	int rank;

	assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
//Triple case
	if(num == 3){
		rank = 2;
		for(i = 0; i <= 2; i++){
			if(i != fmv1->med_idx){
				if(cnt[fmv1->med_idx] > cnt[i])
					rank --;
				else if(cnt[fmv1->med_idx] == cnt[i] && fmv1->med_idx < i)
					rank --;
			}
		}
		putbits(tripleHuffPred[rank].code, tripleHuffPred[rank].len);
		*mapbits += tripleHuffPred[rank].len;
		idx_info += tripleHuffPred[rank].len;
	}
//Quartet case
	else if(num == 4){
//Quartet huffman tree decision
		for(i = 0;i < num;i ++){
			rank = 3;
			for(j = 0;j < num;j ++){
				if(j != i){
					if(cnt[i] > cnt[j])
						rank --;
					else if(cnt[i] == cnt[j] && i < j)
						rank --;
				}
			}
			if(rank == 0)
				max = cnt[i];
			else if(rank == 2)
				min1 = cnt[i];
			else if(rank == 3)
				min2 = cnt[i];
		}
		if(max >  min1 + min2){  //Slanted version
			rank = 3;
			for(i = 0; i < num; i++){
				if(i != fmv1->med_idx){
					if(cnt[fmv1->med_idx] > cnt[i])
						rank --;
					else if(cnt[fmv1->med_idx] == cnt[i] && fmv1->med_idx < i)
						rank --;
				}
			}
			putbits(quartetHuffPredS[rank].code, quartetHuffPredS[rank].len);
			*mapbits += quartetHuffPredS[rank].len;
			idx_info += quartetHuffPredS[rank].len;
		}
		else{  //Flat version
			putbits(fmv1->med_idx,2);
			*mapbits += 2;
			idx_info += 2;
		}
	}
	else
		assert(0);

//Update distribution
	cnt[fmv1->med_idx]++;

}

void input_huff_bits(vector_ptr fmv1, int num, int *cnt, videoinfo info){
	int i,j,bit;
	int max, min1, min2, length;

	int val, getnum, rank;

//Triple case
	if(num == 3){
		rank = -1;
		val = 0;
		length = 0;

		while(rank == -1){
			val <<= 1;
			input_bit(bit);
			val |= bit;
			length++;

			for(i=0;i < num;i++){
				if(tripleHuffPred[i].code == val && tripleHuffPred[i].len == length){
					rank = i;
					break;
				}
			}
		}//while

		for(i=0;i< num;i++){
			getnum = 2;
			for(j=0;j< num;j++){
				if(j != i){
					if(cnt[i] > cnt[j])
						getnum --;
					else if(cnt[i] == cnt[j] && i < j)
						getnum --;
				}
			}
			if(getnum == rank){
				fmv1->med_idx = i;
				break;
			}
		}
	}
//Quartet case
	else if(num == 4){
//Quartet huffman tree decision
		for(i = 0;i <= 3;i ++){
			rank = 3;
			for(j = 0;j <= 3;j ++){
				if(j != i){
					if(cnt[i] > cnt[j])
						rank --;
					else if(cnt[i] == cnt[j] && i < j)
						rank --;
				}
			}
			if(rank == 0)
				max = cnt[i];
			else if(rank == 2)
				min1 = cnt[i];
			else if(rank == 3)
				min2 = cnt[i];
		}
		if(max > min1 + min2){  //Slanted version
			rank = -1;
			val = 0;
			length = 0;

			while(rank == -1){
				val <<= 1;
				input_bit(bit);
				val |= bit;
				length++;

				for(i=0;i < num;i++){
					if(quartetHuffPredS[i].code == val && quartetHuffPredS[i].len == length){
						rank = i;
						break;
					}
				}
			}//while

			for(i=0;i< num;i++){
				getnum = 3;
				for(j=0;j< num;j++){
					if(j != i){
						if(cnt[i] > cnt[j])
							getnum --;
						else if(cnt[i] == cnt[j] && i < j)
							getnum --;
					}
				}
				if(getnum == rank){
					fmv1->med_idx = i;
					break;
				}
			}
		}
		else{  //Flat version
			getnum = getbits(2);
			fmv1->med_idx = getnum;
		}
	}
	else
		assert(0);

//Update distribution
	cnt[fmv1->med_idx]++;
}

/*
 *  get_blockmode()
 * get block mode and if in PREDICTED, also get means
 */
void
get_blockmode( vector_ptr fmv, int meandepth, videoinfo info )
{
  int bit;

  input_bit( bit );
  if( bit == 0 )
    fmv->lifting_mode = CONNECTED;
  else {
    fmv->lifting_mode = PREDICTED;
  }
}

int get_mode_coding_cost(enum BiMode bi_mode, vector_ptr fmv2, int bi_sign, int t_level)
{
  if (fmv2 == NULL) {
    return uni_mode_VLC[bi_mode].len;
  } else {
	  if ( bi_sign )
	  {
		  if (t_level<=CORRELATED_TEMPORAL_LEVELS)
		  {
			  return bi_mode_VLC_1_012[bi_mode].len;
		  }else 
		  {
			  return bi_mode_VLC_1_345[bi_mode].len;
		  }
	  }
	  else
		  return bi_mode_VLC_0[bi_mode].len;
  }
}

int encode_block_mode(vector_ptr fmv1, vector_ptr fmv2, videoinfo info, int t_level)
{

  int length = 0;
  int getnum;

	if (fmv2 == NULL) {
		putbits(uni_mode_VLC[fmv1->bi_mode].code, 
		uni_mode_VLC[fmv1->bi_mode].len);
		length += uni_mode_VLC[fmv1->bi_mode].len;
	} else {
		// by Yongjun Wu
		if ( info.bi_mv[t_level] ) // bi-directional motion field, including RIGHT_CONNECTED mode 
		{
			if (t_level<=CORRELATED_TEMPORAL_LEVELS)
			{
				if(use_huff == 0){
					putbits(bi_mode_VLC_1_012[fmv1->bi_mode].code, 
							bi_mode_VLC_1_012[fmv1->bi_mode].len);
					length += bi_mode_VLC_1_012[fmv1->bi_mode].len;
				}else{
					assert(use_huff == 1);
					putbits(bi_mode_VLC_enc_1_012[fmv1->bi_mode].code, 
							bi_mode_VLC_enc_1_012[fmv1->bi_mode].len);
					length += bi_mode_VLC_enc_1_012[fmv1->bi_mode].len;
				}

			}else 
			{
				putbits(bi_mode_VLC_1_345[fmv1->bi_mode].code, 
						bi_mode_VLC_1_345[fmv1->bi_mode].len);
				length += bi_mode_VLC_1_345[fmv1->bi_mode].len;
			}
		}else
		{
			assert(0);
			putbits(bi_mode_VLC_0[fmv1->bi_mode].code, 
					bi_mode_VLC_0[fmv1->bi_mode].len);
			length += bi_mode_VLC_0[fmv1->bi_mode].len;
		}
	}

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
    //  code the spatial prediction mode for directional IBLOCK 
	if (fmv1->bi_mode==DIRECTIONAL_IBLOCK)
	{
	    putbits(spatialmodeVLC[fmv1->iblock_spatial_mode].code,
		        spatialmodeVLC[fmv1->iblock_spatial_mode].len);
		length += spatialmodeVLC[fmv1->iblock_spatial_mode].len;
	}
#endif

//////////////////	Added by Yuan Liu	////////////////////
	if(fmv1->bi_mode >= 9){

		if(fmv2 != NULL)
			assert(fmv1->bi_mode == fmv2->bi_mode);

		if(fmv1->direct_idx == DIRECT){  //This bit decides if the affine mode is DIRECT, 0 for DIRECT and 1 for non-DIRECT

			putbits(0,1);  //0 for direct affine mode
			length += 1;
		}
		else{	//INDIRECT case
			putbits(1,1);  //1 for indirect affine mode
			length += 1;

			if(fmv1->merge_idx == MERGE){ 
				putbits(0,1);  //0 for merge mode
				length += 1;

				if(fmv1->bi_mode == BI_CONNECTED_AFF){

					if( !(fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L || fmv1->merge_dir == TRAN_P) ) )
						assert(fmv1->aff_idx >= 0 && fmv2->aff_idx >= 0);

					if(fmv1->merge_dir == UP){ //Merge to the left
						putbits(fmv1->merge_dir,2);
						length += 2;
					}
					else if(fmv1->merge_dir == LEFT){ //Merge to the up
						putbits(fmv1->merge_dir,2);
						length += 2;
					}
					else if(fmv1->merge_dir == PAL_L){ //Must be a parallel affine mode
						putbits(fmv1->merge_dir,2);
						length += 2;
					}else{
						assert(fmv1->merge_dir = TRAN_P);
						putbits(fmv1->merge_dir,2);
						length += 2;

						//transmit trans MV residual or not
						putbits(fmv1->trans_pred_idx,1);
						length += 1;
					}
				}else if(fmv1->bi_mode == LEFT_CONNECTED_AFF){

					if(fmv1->merge_dir == UP){ //Merge to the left
						putbits(0,1);
						length += 1;
					}else if(fmv1->merge_dir == LEFT){ //Merge to the up, 10
						putbits(2,2);
						length += 2;
					}else{  //Must be trans pred, 11
						assert(fmv1->merge_dir == TRAN_P);
						putbits(3,2);
						length += 2;

						//transmit trans MV residual or not
						putbits(fmv1->trans_pred_idx,1);
						length += 1;
					}
				}else{
					assert(fmv1->bi_mode == RIGHT_CONNECTED_AFF);

					if(fmv1->merge_dir == UP){ //Merge to the left
						putbits(0,1);
						length += 1;
					}
					else if(fmv1->merge_dir == LEFT){ //Merge to the up, 10
						putbits(2,2);
						length += 2;
					}else{  //Must be trans pred, 11
						assert(fmv1->merge_dir == TRAN_P);
						putbits(3,2);
						length += 2;

						//transmit trans MV residual or not
						putbits(fmv1->trans_pred_idx,1);
						length += 1;
					}
				}
			}else{ //INTER affine mode
				assert(fmv1->merge_idx == INTER);

				putbits(1,1);  //1 for inter mode
				length += 1;
			}
		}
	}
/////////////////	Added by Yuan Liu	/////////////////////
  return length;
}

void decode_block_mode(vector_ptr fmv1, vector_ptr fmv2, int meandepth, videoinfo info, int t_level)
{
  int i, j;
  int bit, val;
  const VLCtable *bi_mode_VLC;
  int getnum;

  if (fmv2 == NULL) {

    fmv1->bi_mode = UNDEFINED;
	fmv1->mvx     = fmv1->mvy = (float)HUGE_VAL;

	fmv1->aff1_mvx = (float)HUGE_VAL; fmv1->aff1_mvy = (float)HUGE_VAL;
	fmv1->aff2_mvx = (float)HUGE_VAL; fmv1->aff2_mvy = (float)HUGE_VAL;
	fmv1->aff3_mvx = (float)HUGE_VAL; fmv1->aff3_mvy = (float)HUGE_VAL;

    i = 0;
    val = 0;
    while (fmv1->bi_mode == UNDEFINED) {
      i++;
      val <<= 1;
      input_bit(bit);
      val |= bit;
      for (j = 0; j < NUMBER_OF_BI_MODES; j++) {
        if ((uni_mode_VLC[j].len == i) && (uni_mode_VLC[j].code == val)) {
          fmv1->bi_mode = (BiMode)j;
          fmv1->lifting_mode = left_mode[j];
          break;
        }
      }
    }

  } else {

    fmv1->bi_mode = UNDEFINED;
	fmv1->mvx     = fmv1->mvy = (float)HUGE_VAL;
	fmv2->mvx     = fmv2->mvy = (float)HUGE_VAL;

	fmv1->aff1_mvx = (float)HUGE_VAL; fmv1->aff1_mvy = (float)HUGE_VAL;
	fmv1->aff2_mvx = (float)HUGE_VAL; fmv1->aff2_mvy = (float)HUGE_VAL;
	fmv1->aff3_mvx = (float)HUGE_VAL; fmv1->aff3_mvy = (float)HUGE_VAL;

	fmv2->aff1_mvx = (float)HUGE_VAL; fmv2->aff1_mvy = (float)HUGE_VAL;
	fmv2->aff2_mvx = (float)HUGE_VAL; fmv2->aff2_mvy = (float)HUGE_VAL;
	fmv2->aff3_mvx = (float)HUGE_VAL; fmv2->aff3_mvy = (float)HUGE_VAL;

    i = 0;
    val = 0;
	if (t_level<=CORRELATED_TEMPORAL_LEVELS)
	{
		if(use_huff == 0)
			bi_mode_VLC= bi_mode_VLC_1_012;
		else
			bi_mode_VLC= bi_mode_VLC_enc_1_012;

	}else{
		bi_mode_VLC= bi_mode_VLC_1_345;
	}

    while (fmv1->bi_mode == UNDEFINED) {
      i++;
      val <<= 1;
      input_bit(bit);
      val |= bit;
      for (j = 0; j < NUMBER_OF_BI_MODES; j++) {
		if ( info.bi_mv[t_level] )
		{
			if ((bi_mode_VLC[j].len == i) && (bi_mode_VLC[j].code == val)) {
			  fmv1->bi_mode = (BiMode)j;
			  fmv2->bi_mode = (BiMode)j;
			  fmv1->lifting_mode = left_mode[j];
			  fmv2->lifting_mode = right_mode[j];
			  break;
			}

		}else
		{
			if ((bi_mode_VLC_0[j].len == i) && (bi_mode_VLC_0[j].code == val)) {
			  fmv1->bi_mode = (BiMode)j;
			  fmv2->bi_mode = (BiMode)j;
			  fmv1->lifting_mode = left_mode[j];
			  fmv2->lifting_mode = right_mode[j];
			  break;
			}
		}
      }
    }
  }

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
    //  code the spatial prediction mode for directional IBLOCK 
	if (fmv1->bi_mode==DIRECTIONAL_IBLOCK)
	{
		fmv1->iblock_spatial_mode = INVALID_SPATIAL_MODE;
		i = 0;
		val = 0;
		while (fmv1->iblock_spatial_mode == INVALID_SPATIAL_MODE) {
			i++;
			val <<= 1;
			input_bit(bit);
			val |= bit;
			for (j = 0; j < PRE_MODE_NUM; j++) {
				if ((spatialmodeVLC[j].len == i) && (spatialmodeVLC[j].code == val)) {
					fmv1->iblock_spatial_mode = (spatialMODE)j;
					if (fmv2 != NULL)
						fmv2->iblock_spatial_mode = (spatialMODE)j;
					break;
				}
			}
		}
	}
#endif

	if(fmv1->bi_mode == BLOCK_MERGING){
		fmv1->aff_mrg = NO;

		if(fmv2 != NULL)
			fmv2->aff_mrg = NO;
	}

//////////////////	Added by Yuan Liu	////////////////////
	if(fmv1->bi_mode >= 9){

		val = getbits(1);

		if(val == 0){  //DIRECT affine mode
			fmv1->direct_idx = DIRECT;
			if(fmv2 != NULL)
				fmv2->direct_idx = DIRECT;
		}else if(val == 1){//INDIRECT affine mode

			fmv1->direct_idx = INDIRECT;
			if(fmv2 != NULL)
				fmv2->direct_idx = INDIRECT;

			val = getbits(1);
			if(val == 1){  //INTER mode

				fmv1->merge_idx = INTER;
				if(fmv2 != NULL)
					fmv2->merge_idx = INTER;
			}
			else if(val == 0){  //MERGE mode
				fmv1->merge_idx = MERGE;
				if(fmv2 != NULL)
					fmv2->merge_idx = MERGE;

				if(fmv1->bi_mode == LEFT_CONNECTED_AFF){
					val = getbits(1);//get MERGE direction

					if(val == 0){
						fmv1->merge_dir = UP;
						if(fmv2 != NULL)
							fmv2->merge_dir = UP;
					}else{
						assert(val == 1);
						val = getbits(1);

						if(val == 0){ //LEFT merge mode, code 10
							fmv1->merge_dir = LEFT;
							if(fmv2 != NULL)
								fmv2->merge_dir = LEFT;
						}else{ //TRAN_P merge mode, code 11
							assert(val == 1);

							fmv1->merge_dir = TRAN_P;
							if(fmv2 != NULL)
								fmv2->merge_dir = TRAN_P;

							val = getbits(1);
							fmv1->trans_pred_idx = val;
							if(fmv2 != NULL)
								fmv2->trans_pred_idx = val;
						}
					}
				}else if(fmv1->bi_mode == RIGHT_CONNECTED_AFF){
					val = getbits(1);//get MERGE direction
					
					if(val == 0){
						fmv1->merge_dir = UP;
						if(fmv2 != NULL)
							fmv2->merge_dir = UP;
					}else{
						assert(val == 1);
						val = getbits(1);

						if(val == 0){ //LEFT merge mode, code 10
							fmv1->merge_dir = LEFT;
							if(fmv2 != NULL)
								fmv2->merge_dir = LEFT;
						}else{ //TRAN_P merge mode, code 11
							assert(val == 1);

							fmv1->merge_dir = TRAN_P;
							if(fmv2 != NULL)
								fmv2->merge_dir = TRAN_P;

							val = getbits(1);
							fmv1->trans_pred_idx = val;
							if(fmv2 != NULL)
								fmv2->trans_pred_idx = val;
						}
					}
				}else{
					assert(fmv1->bi_mode == BI_CONNECTED_AFF);
					val = getbits(2);
					assert(val >= 0 && val <= 3);

					fmv1->merge_dir = val;
					fmv2->merge_dir = val;

					if( fmv1->merge_dir == PAL_L ){
						fmv1->aff_idx = 0;
						fmv2->aff_idx = -1;
					}else if(fmv1->merge_dir == TRAN_P){
						val = getbits(1);

						fmv1->trans_pred_idx = val;
						fmv2->trans_pred_idx = val;
					}

/*					if(val == 0){//UP merge
						fmv1->merge_dir = val;
						fmv2->merge_dir = val;
					}else{//LEFT or parallel merge
						assert(val == 1);
						val = getbits(1);
						val = val + 2;

						fmv1->merge_dir = val;
						fmv2->merge_dir = val;

						if( fmv1->merge_dir == PAL_L ){
							fmv1->aff_idx = 0;
							fmv2->aff_idx = -1;
						}
					}
*/
				}
			}
		}
	}
/////////////////	Added by Yuan Liu	/////////////////////
}


/****************************************************************************/
/*                              child_mv_encode_further()                   */
/****************************************************************************/
void
child_mv_encode_further( vector_ptr fmv1_array, vector_ptr fmv1, vector_ptr fmv2, 
                 float *pmvx, float *pmvy, int num_symbol, int subpel, 
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int encode_parallelmv, int t_level, int blk_thresh, int count)
{
  int dmvx, dmvy, cx, cy;
  int ctx_x, ctx_y;
  float mvpred_x, mvpred_y;
  float major_mvx, major_mvy, major_predx, major_predy, this_mvx, this_mvy;
  int   sub_symx, sub_symy, sub_sym_predx, sub_sym_predy; 
  int   AGP_scale, subpel_scale, sub_bit; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i,num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  printf("child_mv_encode_further\n");
  assert(0);

  if( fmv1->child && xblk>blk_thresh) {
    cx = x;
    cy = y;
    child_mv_encode_further(fmv1_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_encode_further(fmv1_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_encode_further(fmv1_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_encode_further(fmv1_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);
  } else {
    if( x >= hor || y >= ver )    return;

	if (!(fmv1->child))  // there are no children for this block with size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
	    this_mvx = fmv1->mvx;  this_mvy = fmv1->mvy;
	} else
	{
		// this is a merged (blk_thresh x blk_thresh ) block
		if ( !fmv1->mv_exist) return;   // no motion vector for this big block 

		// irregular motion area
		if ( fmv1->child0->child || fmv1->child1->child ||
			 fmv1->child2->child || fmv1->child3->child)
		{
			blk_thresh = 0; 
			child_mv_encode_further( fmv1_array, fmv1, fmv2, pmvx, pmvy,  num_symbol,  subpel, 
							x,  y,  xblk,  yblk,  hor,  ver, info,  encode_parallelmv,  t_level,  
							blk_thresh,  count);
			return;
		}
			 
		if ( fmv1->merge_sign  )
		{
			this_mvx = fmv1->sample_mvx;   
			this_mvy = fmv1->sample_mvy;	
		}else
		{
			blk_thresh = 0; 
			child_mv_encode_further( fmv1_array, fmv1, fmv2, pmvx, pmvy,  num_symbol,  subpel, 
							x,  y,  xblk,  yblk,  hor,  ver, info,  encode_parallelmv,  t_level,  
							blk_thresh,  count);
			return; 
		}
	}

    if (this_mvx == (float)HUGE_VAL || this_mvy == (float)HUGE_VAL){
	  if( (fmv1->bi_mode >= 9 && fmv1->bi_mode <= 11) || (fmv1->bi_mode == 7 && fmv1->aff_mrg == YES) ){
		  assert(fmv1->aff1_mvx != (float)HUGE_VAL && fmv1->aff1_mvy != (float)HUGE_VAL && fmv1->aff2_mvx != (float)HUGE_VAL &&
			  fmv1->aff2_mvy != (float)HUGE_VAL && fmv1->aff3_mvx != (float)HUGE_VAL && fmv1->aff3_mvy != (float)HUGE_VAL);
	  }
	  else{
		printf("error in mvcoding.c: mvx / mvy = (float)HUGE_VAL; " 
				 "x = %d, y = %d, xblk = %d, yblk = %d\n", x, y, xblk, yblk );
		exit ( 1 );
	  }
    }
	assert(0);
    get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, blk_thresh);
	/////////////////////////////////
	if( !((fmv1->bi_mode == PARALLEL) && !encode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}

		if(fmv1->med_idx == -1){
			assert(num == 0);
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
//			putbits(fmv1->med_idx,2);
			mvpred_x = mvpredstr_x[fmv1->med_idx];
			mvpred_y = mvpredstr_y[fmv1->med_idx];
		}
	}
	else if(fmv1->bi_mode == BLOCK_MERGING ){
	  num = 0;
	  for(i=0;i<=3;i++){
		  if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
			num++;
	  }
	  switch(num){
		  case 0:
			  mvpred_x = 0.0;
			  mvpred_y = 0.0;
			  break;

		  case 1:
		  case 2:
			  if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[0];
			    mvpred_y = mvpredstr_y[0];
			  }
			  else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[1];
			    mvpred_y = mvpredstr_y[1];
			  }
			  else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[2];
			    mvpred_y = mvpredstr_y[2];
			  }
			  else{
				assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
				mvpred_x = mvpredstr_x[3];
			    mvpred_y = mvpredstr_y[3];
			  }
			  break;

		  case 3:
			  if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
				mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
				mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
			  }
			  else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
				mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
				mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
			  }
			  else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
				mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
				mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
			  }
			  else{
				assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
				mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
				mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
			  }
			  break;
	  }
	}
	/////////////////////////////////

	if (!info.AGP_level[t_level]) // no AGP, layer structure( this function is always for layer structure)
	{
		fmv1->dmvx = this_mvx - mvpred_x;  
		fmv1->dmvy = this_mvy - mvpred_y;  
	}
	
	if (info.AGP_level[t_level])  // AGP, layer structure 
	{
		AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
		subpel_scale = 1<<subpel; 
		// major symbol
		major_mvx = (int)(this_mvx*AGP_scale)/(float)AGP_scale;
		// sub-symbol (already converted to integer)
		sub_symx = (char)(fabs( (this_mvx - major_mvx) * subpel_scale ));  
		major_mvy = (int)(this_mvy*AGP_scale)/(float)AGP_scale;
		sub_symy = (char)(fabs( (this_mvy - major_mvy) * subpel_scale ));  
		major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
		sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
		major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
		sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale));  
		// only major symbols are median predicted, sub-symbols are coded by binary sequence
		fmv1->dmvx = major_mvx - major_predx;
		fmv1->dmvy = major_mvy - major_predy;

#ifdef  DEBUG_SCALABLE_MV		
		fprintf(fpAGP_debug, "x=%03d  y=%03d  blk=%2d ", x, y, xblk ); 
		fprintf(fpAGP_debug, "major_mvx=%.1f  major_mvy=%.1f  major_predx=%.1f  major_predy=%.1f",
			    major_mvx,  major_mvy,  major_predx,  major_predy) ;
		fprintf(fpAGP_debug, " dmvx=%.1f  dmvy=%.1f  mvx=%.2f  mvy=%.2f\n", 
			                 fmv1->dmvx,  fmv1->dmvy, this_mvx, this_mvy) ;
#endif

		// code the sub-symbols as binary sequence 
		for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
		{
			put_splitted_mvbits( (sub_symx & (1<<sub_bit))>>sub_bit, sub_bit, 0);
			put_splitted_mvbits( (sub_symy & (1<<sub_bit))>>sub_bit, sub_bit, 0);
		}
		// additional sign bit
		if (major_mvx==0) // 0: positive, 1: negative
			put_splitted_signbits( (this_mvx<0), 0 );  
		if (major_mvy==0) // 0: positive, 1: negative
			put_splitted_signbits( (this_mvy<0), 0 );  

	}



#ifdef  DEBUG_LAYER_STRUCTURE
	char  base_file[80];
	FILE *fpbase; 
	// make base_file and enhance_file empty
	sprintf(base_file, "base%d.txt", count); 
	fpbase=fopen(base_file, "at");
	fprintf(fpbase, "x=%03d y=%03d blk=%d\t mvx=%.2f\t mvy=%.2f\t predx=%.2f predy=%.2f \n", 
		x, y, xblk, this_mvx, this_mvy, mvpred_x, mvpred_y); 
    fclose(fpbase);
#endif

	// prediction error after AGP 
    dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->dmvx);
    dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->dmvy);

    mvStat_setPos(x, y);
    mvStat_setDMV((float)dmvx, (float)dmvy);

	// get the context for coding this motion vector 
    ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, x, y, info, t_level, blk_thresh);
	ec_encode_word(dmvx, ctx_x);
	ec_update_model(dmvx, ctx_x);
	ec_encode_word(dmvy, ctx_y);
	ec_update_model(dmvy, ctx_y);
  	fmv1->is_predictor = YES;
	mvStat_writeDMVCTX();
	// update the motion feild in this frame 
	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
							  fmv1->dmvy, this_mvx, this_mvy); 

  }

}


/****************************************************************************/
/*                              child_mv_encode()                           */
/****************************************************************************/
void
child_mv_encode( vector_ptr fmv1_array, vector_ptr fmv1, vector_ptr fmv2, 
                 float *pmvx, float *pmvy, int num_symbol, int subpel, 
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int encode_parallelmv, int t_level, int blk_thresh, int count)
{
  int dmvx, dmvy, cx, cy;
  int ctx_x, ctx_y;
  int enc_trans;
  float mvpred_x, mvpred_y;
  int xblk2, yblk2; // debug
  float major_mvx, major_mvy, major_predx, major_predy, this_mvx, this_mvy;
  int   sub_symx, sub_symy, sub_sym_predx, sub_sym_predy; 
  int   AGP_scale, subpel_scale, sub_bit; 
  int aff1_dmvx, aff1_dmvy, aff2_dmvx, aff2_dmvy, aff3_dmvx, aff3_dmvy;

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0, num_dir = 0;
  int getnum;

  float getval;

  xblk2 = ( x + xblk <= hor ) ? xblk : hor - x;
  yblk2 = ( y + yblk <= ver ) ? yblk : ver - y;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk>blk_thresh) {
    cx = x;
    cy = y;
    child_mv_encode(fmv1_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_encode(fmv1_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_encode(fmv1_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_encode(fmv1_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);
  } else {
    if( x >= hor || y >= ver )    return;

	if (!fmv1->child)  // there are no children for this block with size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 

		this_mvx = fmv1->mvx;  
		this_mvy = fmv1->mvy;

	} else   //Never gonna happen, don't care!
	{
		printf("this is a merged (blk_thresh x blk_thresh ) block\n");
		// this is a merged (blk_thresh x blk_thresh ) block
		if ( !fmv1->mv_exist ) return;   // no motion vector for this big block 

		if ( fmv1->merge_sign  )
		{
			this_mvx = fmv1->sample_mvx;   
			this_mvy = fmv1->sample_mvy;	
		}else
		{
			child_mv_encode_further( fmv1_array, fmv1, fmv2, pmvx, pmvy,  num_symbol,  subpel, 
							x,  y,  xblk,  yblk,  hor,  ver, info,  encode_parallelmv,  t_level,  
							blk_thresh/2,  count);
			return; 
		}
	}

	assert(fmv1->bi_mode >= 0);	//Added ib 04.02.2016
	if(fmv2!=NULL)
		assert(fmv1->bi_mode == fmv2->bi_mode);

    if (this_mvx == (float)HUGE_VAL || this_mvy == (float)HUGE_VAL){
	  if( (fmv1->bi_mode >= 9 && fmv1->bi_mode <= 11) || (fmv1->bi_mode == 7 && fmv1->aff_mrg == YES)){
		  if(fmv1->aff1_mvx == (float)HUGE_VAL || fmv1->aff1_mvy == (float)HUGE_VAL || fmv1->aff2_mvx == (float)HUGE_VAL ||
			  fmv1->aff2_mvy == (float)HUGE_VAL || fmv1->aff3_mvx == (float)HUGE_VAL || fmv1->aff3_mvy == (float)HUGE_VAL){
			
			printf("\nbi_mode = %d, sad_cost = %f, med_idx = %d, x = %d, y = %d, xblk = %d, yblk = %d\nfmv1->mvx = %f, fmv1->mvy = %f, fmv1->dmvx = %f, fmv1->dmvy = %f\n"
				,fmv1->bi_mode, fmv1->sad_cost, fmv1->med_idx,x,y,xblk,yblk,fmv1->mvx,fmv1->mvy,fmv1->dmvx,fmv1->dmvy);

			if(fmv1->direct_idx == DIRECT){
				printf("aff_index = %d\nfmv1->aff1_mvx = %f, fmv1->aff1_mvy = %f\n fmv1->aff2_mvx = %f, fmv1->aff2_mvy = %f\n fmv1->aff3_mvx = %f, fmv1->aff3_mvy = %f\n"
				,fmv1->aff_idx,fmv1->aff1_mvx, fmv1->aff1_mvy, fmv1->aff2_mvx , fmv1->aff2_mvy, fmv1->aff3_mvx, fmv1->aff3_mvy);
			}
			else if(fmv1->direct_idx != DIRECT){
				printf("aff_index = %d\nfmv1->aff1_mvx = %f, fmv1->aff1_mvy = %f\n fmv1->aff2_mvx = %f, fmv1->aff2_mvy = %f\n fmv1->aff3_mvx = %f, fmv1->aff3_mvy = %f\nmerge_idx = %d, merge_dir = %d\n"
				,fmv1->aff_idx,fmv1->aff1_mvx, fmv1->aff1_mvy, fmv1->aff2_mvx , fmv1->aff2_mvy, fmv1->aff3_mvx, fmv1->aff3_mvy,fmv1->merge_idx,fmv1->merge_dir);
				if(fmv1->merge_idx == MERGE && fmv1->merge_dir == UP)
					printf("fmv1->aff3_dmvx = %f, fmv1->aff3_dmvy = %f\n", fmv1->aff3_dmvx, fmv1->aff3_dmvy);
				else if(fmv1->merge_idx == MERGE && fmv1->merge_dir == LEFT)
					printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n", fmv1->aff2_dmvx, fmv1->aff2_dmvy);
				else{
					assert(fmv1->merge_idx == INTER);
					printf("\nfmv1->aff1_dmvx = %f, fmv1->aff1_dmvy = %f\n", fmv1->aff1_dmvx, fmv1->aff1_dmvy);
					printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n", fmv1->aff2_dmvx, fmv1->aff2_dmvy);
					printf("fmv1->aff3_dmvx = %f, fmv1->aff3_dmvy = %f\n\n", fmv1->aff3_dmvx, fmv1->aff3_dmvy);
				}
			}
			assert(0);
		  }
	  }
	  else{
		printf("error in mvcoding.c: mvx / mvy = (float)HUGE_VAL; " 
				 "x = %d, y = %d, xblk = %d, yblk = %d\n", x, y, xblk, yblk );
		exit ( 1 );
	  }
    }

    get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
							x, y, xblk, yblk, info, t_level, blk_thresh);  //Get translational motion field, may be still useful in AFFINE MODEL

	/////////////////////////////////
	if( (fmv1->bi_mode <= 8 && !(fmv1->bi_mode == PARALLEL && !encode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING) ||
		( fmv1->bi_mode>=9 && fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == INDIR) ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL){
				num++;
				if(fmv1->bi_mode == BLOCK_MERGING){
				  if( (float(x)-mvpredstr_x[i] >= 0) && (float(x)-mvpredstr_x[i]<= (hor - xblk2)) && (float(y)-mvpredstr_y[i] >= 0)
					  && (float(y) - mvpredstr_y[i] <= (ver - yblk2)) )
						num_dir++;
				}
			}
		}

		if(num == 0){
			assert(fmv1->med_idx == -1);
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			if(fmv1->bi_mode == BLOCK_MERGING){
				if(num_dir > 0){//Encode the med_idx only if there is something fit for BLOCK_MERGING mode
					assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
					mvpred_x = mvpredstr_x[fmv1->med_idx];
					mvpred_y = mvpredstr_y[fmv1->med_idx];
				}
				else{
					assert(fmv1->med_idx == -1);
					mvpred_x = 0.0;
					mvpred_y = 0.0;
				}
			}
			else{
				assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
				mvpred_x = mvpredstr_x[fmv1->med_idx];
				mvpred_y = mvpredstr_y[fmv1->med_idx];
			}
		}
	}else if( (fmv1->bi_mode == PARALLEL && !encode_parallelmv) ){
			fmv1->med_idx = fmv2->med_idx;
	}

	if (!info.AGP_level[t_level])  // AGP, layer structure only or no layer
	{
		if( ((fmv1->bi_mode == PARALLEL) && !encode_parallelmv) || (fmv1->bi_mode == BLOCK_MERGING) ){
			fmv1->dmvx = 0.0;
			fmv1->dmvy = 0.0;
		}else if( fmv1->bi_mode >= 9 && fmv1->bi_mode <= 11 ){
			if(fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == INDIR){
				assert( fmv1->dmvx == (fmv1->mvx - mvpred_x) && fmv1->dmvy == (fmv1->mvy - mvpred_y) );
			}else{
				fmv1->dmvx = 0.0;
				fmv1->dmvy = 0.0;
			}

		}else{
			fmv1->dmvx = this_mvx - mvpred_x;  
			fmv1->dmvy = this_mvy - mvpred_y;
		}
		///////////////////////////////////////
		if( (fmv1->bi_mode == BLOCK_MERGING) && (fmv1->dmvx != 0 || fmv1->dmvy!= 0) ){
			assert(0);
		}
		///////////////////////////////////////
	}
	else if (info.AGP_level[t_level])  // AGP - Never happened, don't care!
	{
		AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
		subpel_scale = 1<<subpel; 
		// major symbol
		major_mvx = (int)(this_mvx*AGP_scale)/(float)AGP_scale;
		// sub-symbol (already converted to integer)
		sub_symx = (char)(fabs( (this_mvx - major_mvx) * subpel_scale ));  
		major_mvy = (int)(this_mvy*AGP_scale)/(float)AGP_scale;
		sub_symy = (char)(fabs( (this_mvy - major_mvy) * subpel_scale ));  
		major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
		sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
		major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
		sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale));  
		// only major symbols are median predicted, sub-symbols are coded by binary sequence
		fmv1->dmvx = major_mvx - major_predx;
		fmv1->dmvy = major_mvy - major_predy;

#ifdef  DEBUG_SCALABLE_MV		
		fprintf(fpAGP_debug, "x=%03d  y=%03d  blk=%2d ", x, y, xblk ); 
		fprintf(fpAGP_debug, "major_mvx=%.1f  major_mvy=%.1f  major_predx=%.1f  major_predy=%.1f",
			    major_mvx,  major_mvy,  major_predx,  major_predy) ;
		fprintf(fpAGP_debug, " dmvx=%.1f  dmvy=%.1f  mvx=%.2f  mvy=%.2f\n", 
			                 fmv1->dmvx,  fmv1->dmvy, this_mvx, this_mvy) ;
#endif

		if ( !info.layer_mv[t_level] ) // no layer structure: put sub-symbol
		{
			// motion vector in BLOCK_MERGING block should be perfectly predicted 
			assert((fmv1->bi_mode != BLOCK_MERGING)||
				   (sub_symx == sub_sym_predx && sub_symy==sub_sym_predy &&
					fmv1->dmvx ==0 && fmv1->dmvy == 0 ));
			if ( (fmv1->bi_mode != BLOCK_MERGING) && !( (fmv1->bi_mode == PARALLEL) && !encode_parallelmv ) ) {
				// code the sub-symbols as binary sequence 
				for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
				{
					put_splitted_mvbits( (sub_symx & (1<<sub_bit))>>sub_bit, sub_bit, 0);
					put_splitted_mvbits( (sub_symy & (1<<sub_bit))>>sub_bit, sub_bit, 0);
				}
				// additional sign bit
				if (major_mvx==0) // 0: positive, 1: negative
					put_splitted_signbits( (this_mvx<0), 0 );  
				if (major_mvy==0) // 0: positive, 1: negative
					put_splitted_signbits( (this_mvy<0), 0 );  
			}
		}else  // layer structure : put sub-symbol
		{
			// code the sub-symbols as binary sequence 
			for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
			{
				put_splitted_mvbits( (sub_symx & (1<<sub_bit))>>sub_bit, sub_bit, 0);
				put_splitted_mvbits( (sub_symy & (1<<sub_bit))>>sub_bit, sub_bit, 0);
			}
			// additional sign bit
			if (major_mvx==0) // 0: positive, 1: negative
				put_splitted_signbits( (this_mvx<0), 0 );  
			if (major_mvy==0) // 0: positive, 1: negative
				put_splitted_signbits( (this_mvy<0), 0 );  
		}
	}// AGP

#ifdef  DEBUG_LAYER_STRUCTURE
	char  base_file[80];
	FILE *fpbase; 
	// make base_file and enhance_file empty
	sprintf(base_file, "base%d.txt", count); 
	fpbase=fopen(base_file, "at");
	fprintf(fpbase, "x=%03d y=%03d blk=%d\t mvx=%.2f\t mvy=%.2f\t predx=%.2f predy=%.2f \n", 
		x, y, xblk, this_mvx, this_mvy, mvpred_x, mvpred_y); 
    fclose(fpbase);
#endif

    if (!info.layer_mv[t_level] ) // no layer structure, AGP or no AGP  
		assert( fmv1->bi_mode != BLOCK_MERGING || fmv1->dmvx==0 && fmv1->dmvy==0 ); 

	// prediction error after AGP or no AGP 
    dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->dmvx);
    dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->dmvy);

	if( fmv1->bi_mode >= 9 && fmv1->bi_mode <= 11 ){
		if(fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE){
			if(fmv1->merge_dir == LEFT){
				aff2_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff2_dmvx);
				aff2_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff2_dmvy);
	//			printf("aff2_dmvx = %d, aff2_dmvy = %d\n",aff2_dmvx,aff2_dmvy);
			}else if(fmv1->merge_dir == UP){
				aff3_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff3_dmvx);
				aff3_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff3_dmvy);
	//			printf("aff3_dmvx = %d, aff3_dmvy = %d\n",aff3_dmvx,aff3_dmvy);
			}else if( (fmv1->merge_dir == PAL_L && fmv1->aff_idx >= 0) || fmv1->merge_dir == TRAN_P ){
				aff1_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff1_dmvx);
				aff1_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff1_dmvy);
	//			printf("aff1_dmvx = %d, aff1_dmvy = %d\n",aff1_dmvx,aff1_dmvy);

				aff2_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff2_dmvx);
				aff2_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff2_dmvy);
	//			printf("aff2_dmvx = %d, aff2_dmvy = %d\n",aff2_dmvx,aff2_dmvy);

				aff3_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff3_dmvx);
				aff3_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff3_dmvy);
	//			printf("aff3_dmvx = %d, aff3_dmvy = %d\n",aff3_dmvx,aff3_dmvy);
			}
		}else if(fmv1->direct_idx == INDIRECT && fmv1->merge_idx == INTER){
			aff1_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff1_dmvx);
			aff1_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff1_dmvy);
//			printf("aff1_dmvx = %d, aff1_dmvy = %d\n",aff1_dmvx,aff1_dmvy);

			aff2_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff2_dmvx);
			aff2_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff2_dmvy);
//			printf("aff2_dmvx = %d, aff2_dmvy = %d\n",aff2_dmvx,aff2_dmvy);

			aff3_dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff3_dmvx);
			aff3_dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->aff3_dmvy);
//			printf("aff3_dmvx = %d, aff3_dmvy = %d\n",aff3_dmvx,aff3_dmvy);
		}
	}

    mvStat_setPos(x, y);
    mvStat_setDMV((float)dmvx, (float)dmvy);

    if (!info.layer_mv[t_level])  // no layer structure, AGP or no AGP
	{
		if ( (fmv1->bi_mode <=8 && fmv1->bi_mode != BLOCK_MERGING && !(fmv1->bi_mode == PARALLEL && !encode_parallelmv)) ||
			( (fmv1->bi_mode>=9 && fmv1->bi_mode<=11) && (fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == INDIR) ) ) // Translational model
		{
			ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field,  x, y, info, t_level, blk_thresh);

 		    if (EC_TYPE == AR_NARY)
			{
				if (dmvx > num_symbol / 2)
					dmvx -= num_symbol;
				else if (dmvx < -(num_symbol / 2))
					dmvx += num_symbol;
				if (dmvy > num_symbol / 2)
					dmvy -= num_symbol;
				else if (dmvy < -(num_symbol / 2))
					dmvy += num_symbol;
			 }
			 ec_encode_word(dmvx, ctx_x);
			 ec_update_model(dmvx, ctx_x);
			 ec_encode_word(dmvy, ctx_y);
			 ec_update_model(dmvy, ctx_y);
		}
		
		if( fmv1->bi_mode >=9 && fmv1->bi_mode <=11 ){  // Affine model employed

		  if(fmv1->direct_idx == INDIRECT){
//affine V1
			  if(fmv1->merge_idx == INTER || 
				(fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0 )
				|| (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
				ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field,  x, y, info, t_level, blk_thresh);

 				if (EC_TYPE == AR_NARY)
				{
					if (aff1_dmvx > num_symbol / 2)
						aff1_dmvx -= num_symbol;
					else if (aff1_dmvx < -(num_symbol / 2))
						aff1_dmvx += num_symbol;
					if (aff1_dmvy > num_symbol / 2)
						aff1_dmvy -= num_symbol;
					else if (aff1_dmvy < -(num_symbol / 2))
						aff1_dmvy += num_symbol;
				 }
				 ec_encode_word(aff1_dmvx, ctx_x);
				 ec_update_model(aff1_dmvx, ctx_x);
				 ec_encode_word(aff1_dmvy, ctx_y);
				 ec_update_model(aff1_dmvy, ctx_y);
			}
//affine V2
			if(fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && fmv1->merge_dir == LEFT) ||
			   (fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0)
			   || (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
				ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field,  x+xblk2-1, y, info, t_level, blk_thresh);

 				if (EC_TYPE == AR_NARY)
				{
					if (aff2_dmvx > num_symbol / 2)
						aff2_dmvx -= num_symbol;
					else if (aff2_dmvx < -(num_symbol / 2))
						aff2_dmvx += num_symbol;
					if (aff2_dmvy > num_symbol / 2)
						aff2_dmvy -= num_symbol;
					else if (aff2_dmvy < -(num_symbol / 2))
						aff2_dmvy += num_symbol;
				 }
				 ec_encode_word(aff2_dmvx, ctx_x);
				 ec_update_model(aff2_dmvx, ctx_x);
				 ec_encode_word(aff2_dmvy, ctx_y);
				 ec_update_model(aff2_dmvy, ctx_y);
			}
//affine V3
			if(fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && fmv1->merge_dir == UP) ||
			   (fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0) ||
			   (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
				ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, x, y+yblk2-1, info, t_level, blk_thresh);

 				if (EC_TYPE == AR_NARY)
				{
					if (aff3_dmvx > num_symbol / 2)
						aff3_dmvx -= num_symbol;
					else if (aff3_dmvx < -(num_symbol / 2))
						aff3_dmvx += num_symbol;
					if (aff3_dmvy > num_symbol / 2)
						aff3_dmvy -= num_symbol;
					else if (aff3_dmvy < -(num_symbol / 2))
						aff3_dmvy += num_symbol;
				 }
				 ec_encode_word(aff3_dmvx, ctx_x);
				 ec_update_model(aff3_dmvx, ctx_x);
				 ec_encode_word(aff3_dmvy, ctx_y);
				 ec_update_model(aff3_dmvy, ctx_y);
			}
		  }// if INDIRECT
		  else //DIRECT mode, don't code anything
		  {
			  assert(fmv1->direct_idx == DIRECT);
		  }
		}

		assert(fmv1->is_predictor == YES);
		mvStat_writeDMVCTX();
	}else  // layer structure, AGP or no AGP  
	{
		printf("layer structure detected!\n");
		assert(0);

	    ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field,
								     x, y, info, t_level, blk_thresh);
		ec_encode_word(dmvx, ctx_x);
		ec_update_model(dmvx, ctx_x);
		ec_encode_word(dmvy, ctx_y);
		ec_update_model(dmvy, ctx_y);
	  	fmv1->is_predictor = YES;
		mvStat_writeDMVCTX();
	}

	/////////////////////////
//	printf("bi_mode = %d, x = %d, y = %d, xblk = %d, yblk = %d, ctx_x = %d, ctx_y = %d\n fmv->dmvx = %f, fmv->dmvy = %f\n",
//		fmv1->bi_mode,x,y,xblk2,yblk2,ctx_x, ctx_y,fmv1->dmvx,fmv1->dmvy);
	/////////////////////////

	if(fmv1->bi_mode == 1 || fmv1->bi_mode == 2 || fmv1->bi_mode == 4 || fmv1->bi_mode == 5 || fmv1->bi_mode == 10 || fmv1->bi_mode == 11){
		calc_res_sad += fmv1->sad_cost;
	}
	else{
		calc_res_sad += fmv1->sad_cost/2;
	}

	// update the motion field in this frame 
	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
							  fmv1->dmvy, this_mvx, this_mvy); 

  }

}


/****************************************************************************/
/*                              child_mv_encode_enhance_sub()               */
/****************************************************************************/
void
child_mv_encode_enhance_sub( vector_ptr fmv1_array, vector_ptr fmv1, vector_ptr fmv2, 
                 float *pmvx, float *pmvy, int num_symbol, int subpel, 
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int encode_parallelmv, int t_level, int blk_thresh, 
				 int *first_available, int count, int merge_sign)
{
  int dmvx, dmvy, cx, cy;
  int ctx_x, ctx_y;
  float mvpred_x, mvpred_y;
  int xblk2, yblk2; // debug
  float major_mvx, major_mvy, major_predx, major_predy;
  int   sub_symx, sub_symy, sub_sym_predx, sub_sym_predy; 
  int   AGP_scale, subpel_scale, sub_bit; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  printf("child_mv_encode_enhance_sub\n");

  if( fmv1->child && xblk>blk_thresh) {
    cx = x;
    cy = y;
    child_mv_encode_enhance_sub(fmv1_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, first_available, 
					count, merge_sign);

    cx = x + xblk / 2;
    cy = y;
    child_mv_encode_enhance_sub(fmv1_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, first_available, 
					count, merge_sign);

    cx = x;
    cy = y + yblk / 2;
    child_mv_encode_enhance_sub(fmv1_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, first_available, 
					count, merge_sign);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_encode_enhance_sub(fmv1_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, first_available, 
					count, merge_sign);
  } else {
    if( x >= hor || y >= ver )      return;

	// no motion vector for this block on this side 
	if ((fmv1->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
	if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
	assert(0);
	// re-do the spatial prediction no matter it's first motion vector or not 
    get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
								      x, y, xblk, yblk, info, t_level, 0);
	/////////////////////////////////
	if( !((fmv1->bi_mode == PARALLEL) && !encode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}

		if(fmv1->med_idx == -1){
			assert(num == 0);
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
			putbits(fmv1->med_idx,2);
			mvpred_x = mvpredstr_x[fmv1->med_idx];
			mvpred_y = mvpredstr_y[fmv1->med_idx];
		}
	}
	else if(fmv1->bi_mode == BLOCK_MERGING ){
	  num = 0;
	  for(i=0;i<=3;i++){
		  if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
			num++;
	  }
	  switch(num){
		  case 0:
			  mvpred_x = 0.0;
			  mvpred_y = 0.0;
			  break;

		  case 1:
		  case 2:
			  if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[0];
			    mvpred_y = mvpredstr_y[0];
			  }
			  else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[1];
			    mvpred_y = mvpredstr_y[1];
			  }
			  else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[2];
			    mvpred_y = mvpredstr_y[2];
			  }
			  else{
				assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
				mvpred_x = mvpredstr_x[3];
			    mvpred_y = mvpredstr_y[3];
			  }
			  break;

		  case 3:
			  if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
				mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
				mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
			  }
			  else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
				mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
				mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
			  }
			  else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
				mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
				mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
			  }
			  else{
				assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
				mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
				mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
			  }
			  break;
	  }
	}
	/////////////////////////////////

	if ( !info.AGP_level[t_level] )
	{
		fmv1->dmvx =  fmv1->mvx - mvpred_x;
		fmv1->dmvy =  fmv1->mvy - mvpred_y;
	}

	if (info.AGP_level[t_level])  // AGP 
	{
		AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
		subpel_scale = 1<<subpel; 
		// major symbol
		major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
		// sub-symbol (already converted to integer)
		sub_symx = (char)(fabs( (fmv1->mvx - major_mvx) * subpel_scale ));  
		major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
		sub_symy = (char)(fabs( (fmv1->mvy - major_mvy) * subpel_scale ));  
		major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
		sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
		major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
		sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale)); 
		// only major symbols are median predicted, sub-symbols are coded by binary sequence
		fmv1->dmvx = major_mvx - major_predx;
		fmv1->dmvy = major_mvy - major_predy;
	}

	fmv1->is_predictor = YES;
	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
							  fmv1->dmvy, fmv1->mvx, fmv1->mvy); 

	if (!merge_sign)  	return; 

	if ( *first_available ) // the first motion vector is in base layer
	{
		*first_available = 0; 
		return;  
	}

	// because this is true motion field 
	assert((fmv1->bi_mode != BLOCK_MERGING)||(fmv1->dmvx==0 && fmv1->dmvy==0)); 

	// prediction error after AGP 
    dmvx = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->dmvx);
    dmvy = (int) ((1 << (subpel-info.AGP_level[t_level]) ) * fmv1->dmvy);

    mvStat_setPos(x, y);
    mvStat_setDMV((float)dmvx, (float)dmvy);

#ifdef  DEBUG_LAYER_STRUCTURE
	char  enhance_file[80];
	FILE *fpenhance; 
	// make base_file and enhance_file empty
	sprintf(enhance_file, "enhance%d.txt", count); 
	fpenhance=fopen(enhance_file, "at");
	fprintf(fpenhance, "x=%03d y=%03d blk=%d\t mvx=%.2f\t mvy=%.2f\t predx=%.2f predy=%.2f \n", 
		x, y, xblk, fmv1->mvx, fmv1->mvy, mvpred_x, mvpred_y); 
    fclose(fpenhance);
#endif

	if ( (fmv1->bi_mode != BLOCK_MERGING) && !((fmv1->bi_mode == PARALLEL) && !encode_parallelmv) ) 
	{
	    ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field,
								     x, y, info, t_level, blk_thresh);
  	    ec_encode_word(dmvx, ctx_x);
		ec_update_model(dmvx, ctx_x);
		ec_encode_word(dmvy, ctx_y);
		ec_update_model(dmvy, ctx_y);
	}

	if ( info.AGP_level[t_level] )
	{
		assert((fmv1->bi_mode != BLOCK_MERGING)||
			   (sub_symx == sub_sym_predx && sub_symy==sub_sym_predy &&
				fmv1->dmvx ==0 && fmv1->dmvy == 0 ));
		if (    (fmv1->bi_mode != BLOCK_MERGING) &&
			 !( (fmv1->bi_mode == PARALLEL) && !encode_parallelmv ) ) {
			// code the sub-symbols as binary sequence 
			for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
			{
				put_splitted_mvbits( (sub_symx & (1<<sub_bit))>>sub_bit, sub_bit, 1);
				put_splitted_mvbits( (sub_symy & (1<<sub_bit))>>sub_bit, sub_bit, 1);
			}
			// additional sign bit
			if (major_mvx==0) // 0: positive, 1: negative
				put_splitted_signbits( (fmv1->mvx<0), 1 );  
			if (major_mvy==0) // 0: positive, 1: negative
				put_splitted_signbits( (fmv1->mvy<0), 1 );  
		}
	}

	mvStat_writeDMVCTX();

	xblk2 = ( x + xblk <= hor ) ? xblk : hor - x;
	yblk2 = ( y + yblk <= ver ) ? yblk : ver - y;
	assert( (x + xblk2 - fmv1->mvx <= hor ) &&
			(y + yblk2 - fmv1->mvy <= ver ) &&
			(x - fmv1->mvx >= 0 ) &&
			(y - fmv1->mvy >= 0 ) );
	assert( (!(fmv1->bi_mode == PARALLEL)) || 
			((x + xblk2 - fmv2->mvx <= hor ) &&
			 (y + yblk2 - fmv2->mvy <= ver ) &&
			 (x - fmv2->mvx >= 0 ) &&
			 (y - fmv2->mvy >= 0 )) );
  }
}

/****************************************************************************/
/*                              child_mv_encode_enhance_further()          */
/****************************************************************************/
void
child_mv_encode_enhance_further( vector_ptr fmv1_array, vector_ptr fmv1, vector_ptr fmv2, 
                 float *pmvx, float *pmvy, int num_symbol, int subpel, 
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int encode_parallelmv, int t_level, int blk_thresh, int count)
{
  int cx, cy;
  float mvpred_x, mvpred_y;
  int   first_available; 
  float major_mvx, major_mvy, major_predx, major_predy;
  int   sub_symx, sub_symy, sub_sym_predx, sub_sym_predy; 
  int   AGP_scale, subpel_scale; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk>blk_thresh) {
    cx = x;
    cy = y;
    child_mv_encode_enhance_further(fmv1_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_encode_enhance_further(fmv1_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_encode_enhance_further(fmv1_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_encode_enhance_further(fmv1_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);
  } else {
    if( x >= hor || y >= ver )       return;

	if (!fmv1->child )  // if there are no children for this block with size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if (fmv1->lifting_mode == IGNORED)   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
		assert(0);
		// re-do the spatial prediction 
        get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
	 								      x, y, xblk, yblk, info, t_level, 0);

		/////////////////////////////////
		if( !((fmv1->bi_mode == PARALLEL) && !encode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
			for(i=0;i<=3;i++){
				if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
					num++;
			}

			if(fmv1->med_idx == -1){
				assert(num == 0);
				mvpred_x = 0.0;
				mvpred_y = 0.0;
			}
			else{
				assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
				putbits(fmv1->med_idx,2);
				mvpred_x = mvpredstr_x[fmv1->med_idx];
				mvpred_y = mvpredstr_y[fmv1->med_idx];
			}
		}
		else if(fmv1->bi_mode == BLOCK_MERGING ){
		  num = 0;
		  for(i=0;i<=3;i++){
			  if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		  }
		  switch(num){
			  case 0:
				  mvpred_x = 0.0;
				  mvpred_y = 0.0;
				  break;

			  case 1:
			  case 2:
				  if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[0];
					mvpred_y = mvpredstr_y[0];
				  }
				  else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[1];
					mvpred_y = mvpredstr_y[1];
				  }
				  else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[2];
					mvpred_y = mvpredstr_y[2];
				  }
				  else{
					assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
					mvpred_x = mvpredstr_x[3];
					mvpred_y = mvpredstr_y[3];
				  }
				  break;

			  case 3:
				  if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				  }
				  else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				  }
				  break;
		  }
		}
		/////////////////////////////////

		if ( !info.AGP_level[t_level] )
		{
			fmv1->dmvx =  fmv1->mvx - mvpred_x;
			fmv1->dmvy =  fmv1->mvy - mvpred_y;
		}else // with AGP 
		{
			AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
			subpel_scale = 1<<subpel; 
			// major symbol
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			// sub-symbol (already converted to integer)
			sub_symx = (char)(fabs( (fmv1->mvx - major_mvx) * subpel_scale ));  
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			sub_symy = (char)(fabs( (fmv1->mvy - major_mvy) * subpel_scale ));  
			major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
			sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
			major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
			sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale));  
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
		fmv1->is_predictor = YES;
     	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
						   		  fmv1->dmvy, fmv1->mvx, fmv1->mvy); 
		return;  
	} else
	{
		if ( !fmv1->mv_exist )  return; 

		// classified as irregular motion area
		if ( fmv1->child0->child || fmv1->child1->child ||
			 fmv1->child2->child || fmv1->child3->child)
			 fmv1->merge_sign = 0; 
			
		first_available = 1; 
		blk_thresh = 0;
		child_mv_encode_enhance_sub( fmv1_array, fmv1, fmv2, pmvx, pmvy, num_symbol,  subpel, 
									x, y,  xblk, yblk,  hor,  ver,  info,   encode_parallelmv,  
									t_level,  blk_thresh, &first_available, count, fmv1->merge_sign);
	}
  }
}


/****************************************************************************/
/*                              child_mv_encode_enhance()                    */
/****************************************************************************/
void
child_mv_encode_enhance( vector_ptr fmv1_array, vector_ptr fmv1, vector_ptr fmv2, 
                 float *pmvx, float *pmvy, int num_symbol, int subpel, 
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int encode_parallelmv, int t_level, int blk_thresh, int count)
{
  int cx, cy;
  float mvpred_x, mvpred_y;
  int   first_available; 
  float major_mvx, major_mvy, major_predx, major_predy;
  int   sub_symx, sub_symy, sub_sym_predx, sub_sym_predy; 
  int   AGP_scale, subpel_scale; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  printf("child_mv_encode_enhance\n");

  if( fmv1->child && xblk>blk_thresh) {
    cx = x;
    cy = y;
    child_mv_encode_enhance(fmv1_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_encode_enhance(fmv1_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_encode_enhance(fmv1_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_encode_enhance(fmv1_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);
  } else {
    if( x >= hor || y >= ver )       return;

	if (!fmv1->child )  // if there are no children for this block with size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if (fmv1->lifting_mode == IGNORED)   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
		assert(0);
		// re-do the spatial prediction 
        get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
	 								      x, y, xblk, yblk, info, t_level, 0);
		/////////////////////////////////
		if( !((fmv1->bi_mode == PARALLEL) && !encode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
			for(i=0;i<=3;i++){
				if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
					num++;
			}

			if(fmv1->med_idx == -1){
				assert(num == 0);
				mvpred_x = 0.0;
				mvpred_y = 0.0;
			}
			else{
				assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
				putbits(fmv1->med_idx,2);
				mvpred_x = mvpredstr_x[fmv1->med_idx];
				mvpred_y = mvpredstr_y[fmv1->med_idx];
			}
		}
		else if(fmv1->bi_mode == BLOCK_MERGING ){
		  num = 0;
		  for(i=0;i<=3;i++){
			  if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		  }
		  switch(num){
			  case 0:
				  mvpred_x = 0.0;
				  mvpred_y = 0.0;
				  break;

			  case 1:
			  case 2:
				  if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[0];
					mvpred_y = mvpredstr_y[0];
				  }
				  else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[1];
					mvpred_y = mvpredstr_y[1];
				  }
				  else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[2];
					mvpred_y = mvpredstr_y[2];
				  }
				  else{
					assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
					mvpred_x = mvpredstr_x[3];
					mvpred_y = mvpredstr_y[3];
				  }
				  break;

			  case 3:
				  if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				  }
				  else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				  }
				  break;
		  }
		}
		/////////////////////////////////

		if ( !info.AGP_level[t_level] )
		{
			fmv1->dmvx =  fmv1->mvx - mvpred_x;
			fmv1->dmvy =  fmv1->mvy - mvpred_y;
		}else
		{
			AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
			subpel_scale = 1<<subpel; 
			// major symbol
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			// sub-symbol (already converted to integer)
			sub_symx = (char)(fabs( (fmv1->mvx - major_mvx) * subpel_scale ));  
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			sub_symy = (char)(fabs( (fmv1->mvy - major_mvy) * subpel_scale ));  
			major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
			sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
			major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
			sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale));  
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
		fmv1->is_predictor = YES;
     	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
						   		  fmv1->dmvy, fmv1->mvx, fmv1->mvy); 
		return;  
	} else
	{
		if ( !fmv1->mv_exist ) return; 

		if ( fmv1->merge_sign )
		{
			first_available = 1; 
			blk_thresh = 0;
			child_mv_encode_enhance_sub( fmv1_array, fmv1, fmv2, pmvx, pmvy, num_symbol,  subpel, 
									x, y,  xblk, yblk,  hor,  ver,  info,   encode_parallelmv,  
									t_level,  blk_thresh, &first_available, count, fmv1->merge_sign);
		}else
		{
			child_mv_encode_enhance_further(  fmv1_array,  fmv1,  fmv2, 
							pmvx, pmvy,  num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
							info,  encode_parallelmv,  t_level,  blk_thresh/2,  count);
		}
	}
  }
}


/****************************************************************************/
/*                              child_map_pre_encode()                          */
/****************************************************************************/
void child_map_pre_encode( vector_ptr fmv1, vector_ptr fmv2, int x, int y,
                  int xblk, int yblk, int hor, int ver, int small, 
                  videoinfo info, int t_level )
{
	int cx, cy;

	assert(fmv2 == NULL || fmv1->child == fmv2->child);

	if( fmv1->child ) {

//	ctl_info ++;

    cx = x;
    cy = y;
    child_map_pre_encode( fmv1->child0, fmv2 ? fmv2->child0 : NULL, 
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );

    cx = x + xblk / 2;
    cy = y;
    child_map_pre_encode( fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );

    cx = x;
    cy = y + yblk / 2;
    child_map_pre_encode( fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_map_pre_encode( fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );
  } else {

	  if( x >= hor || y >= ver )     return;

	  if(fmv2 != NULL){
		/////////////////
		assert(fmv1->bi_mode >= 0 && fmv1->bi_mode <= 11);
		if (t_level<=CORRELATED_TEMPORAL_LEVELS)
			bi_mode_num012[fmv1->bi_mode] ++;
		else
			bi_mode_num345[fmv1->bi_mode] ++;
		/////////////////
	  }
  }

}

/****************************************************************************/
/*                              child_map_encode()                          */
/****************************************************************************/
void
child_map_encode( vector_ptr fmv1, vector_ptr fmv2, int *mapbit, int x, int y,
                  int xblk, int yblk, int hor, int ver, int small, 
                  videoinfo info, int t_level )
{
  int cx, cy;

  int middle;

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child ) {
    output_bit( 1 );
    ( *mapbit )++;
//	ctl_info ++;

    cx = x;
    cy = y;
    child_map_encode( fmv1->child0, fmv2 ? fmv2->child0 : NULL, mapbit, 
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );

    cx = x + xblk / 2;
    cy = y;
    child_map_encode( fmv1->child1, fmv2 ? fmv2->child1 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );

    cx = x;
    cy = y + yblk / 2;
    child_map_encode( fmv1->child2, fmv2 ? fmv2->child2 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_map_encode( fmv1->child3, fmv2 ? fmv2->child3 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level );
  } else {
	////////////////////
	avg_blk_size[t_level] += xblk;
	blk_num ++;
	////////////////////
    if( x >= hor || y >= ver )     return;
    if( xblk > small ) {
      output_bit( 0 );
      ( *mapbit )++;          /*for tree structure */
//	  ctl_info ++;
    }
    //output block mode
    middle = encode_block_mode(fmv1, fmv2, info, t_level);

	(*mapbit) += middle;
//	ctl_info += middle;

  }
}

// further selective mergence from 8x8->16x16
void
child_merge_encode_further( vector_ptr fmv1, vector_ptr fmv2, int *mapbit, int x, int y,
                  int xblk, int yblk, int hor, int ver, int small, 
                  videoinfo info, int t_level, int blk_thresh, int count )
{
  int cx, cy;

  assert(fmv2 == NULL || fmv1->child == fmv2->child);
  if( fmv1->child && xblk>blk_thresh) {

    cx = x;
    cy = y;
    child_merge_encode_further( fmv1->child0, fmv2 ? fmv2->child0 : NULL, mapbit, 
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );

    cx = x + xblk / 2;
    cy = y;
    child_merge_encode_further( fmv1->child1, fmv2 ? fmv2->child1 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );

    cx = x;
    cy = y + yblk / 2;
    child_merge_encode_further( fmv1->child2, fmv2 ? fmv2->child2 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_merge_encode_further( fmv1->child3, fmv2 ? fmv2->child3 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );
  } else {
    if( x >= hor || y >= ver )   return;

	if (!fmv1->child)     return; 

	if (!fmv1->mv_exist)  return; // no subsampled motion vector in this big block 

	// irregular motion area
	if ( fmv1->child0->child || fmv1->child1->child ||
		 fmv1->child2->child || fmv1->child3->child )
		 return;

	assert( fmv1->merge_sign==0 || fmv1->merge_sign==1);
	output_bit( fmv1->merge_sign );
	( *mapbit )++;           

	FILE *fmerge; 
	char merge_file[80];
	sprintf(merge_file, "merge_sign%d.txt", count); 
	fmerge = fopen(merge_file, "at"); 
	fprintf(fmerge, "%d\n", fmv1->merge_sign); 
	fclose(fmerge); 

  }
}


/****************************************************************************/
/*                              child_merge_encode()                        */
/****************************************************************************/
void
child_merge_encode( vector_ptr fmv1, vector_ptr fmv2, int *mapbit, int x, int y,
                  int xblk, int yblk, int hor, int ver, int small, 
                  videoinfo info, int t_level, int blk_thresh, int count )
{
  int cx, cy;

  assert(fmv2 == NULL || fmv1->child == fmv2->child);
  if( fmv1->child && xblk>blk_thresh) {

    cx = x;
    cy = y;
    child_merge_encode( fmv1->child0, fmv2 ? fmv2->child0 : NULL, mapbit, 
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );

    cx = x + xblk / 2;
    cy = y;
    child_merge_encode( fmv1->child1, fmv2 ? fmv2->child1 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );

    cx = x;
    cy = y + yblk / 2;
    child_merge_encode( fmv1->child2, fmv2 ? fmv2->child2 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_merge_encode( fmv1->child3, fmv2 ? fmv2->child3 : NULL, mapbit,
                      cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level, 
					  blk_thresh, count );
  } else {
    if( x >= hor || y >= ver )   return;

	if (!fmv1->child)     return; 

	if (!fmv1->mv_exist)  return; // no subsampled motion vector in this big block 

	printf("ain't returned from child_merge_encode");
	assert(0);

	FILE *fmerge; 
	char merge_file[80];
	sprintf(merge_file, "merge_sign%d.txt", count); 
	fmerge = fopen(merge_file, "at"); 
	fprintf(fmerge, "%d\n", fmv1->merge_sign); 
	fclose(fmerge); 
	
	assert( fmv1->merge_sign==0 || fmv1->merge_sign==1);
	output_bit( fmv1->merge_sign );
	( *mapbit )++;           
	if (! fmv1->merge_sign )
		// further selective mergence from 8x8->16x16
		child_merge_encode_further( fmv1, fmv2, mapbit,  x,  y,  xblk,  yblk,  
									hor,  ver,  small, info,  t_level,  blk_thresh/2,  count );
  }
}

/////////////////	Added by Yuan Liu on 01.30.2016   /////////////////////
/*****************************************************************************/
/*                            child_cand_encode()                            */
/*****************************************************************************/
void 
child_cand_encode( int *mapbit, vector_ptr fmv1_array, vector_ptr fmv1, vector_ptr fmv2_array, vector_ptr fmv2, 
                 float *pmvx, float *pmvy, int num_symbol, int subpel, 
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int encode_parallelmv, int t_level, int blk_thresh, int count )
{
	int xblk2,yblk2,cx,cy;
	float mvpredstr_x[4], mvpredstr_y[4];
	int i, num = 0, num_dir = 0, mrg_num = 0;
	int getnum;
	int ctx_x, ctx_y;
	int middle;

	FLAG do_fmv1 = YES, do_fmv2 = YES;
	int map_side;
	vector_ptr mrg_left[4],mrg_right[4];

	xblk2 = ( x + xblk <= hor ) ? xblk : hor - x;
	yblk2 = ( y + yblk <= ver ) ? yblk : ver - y;

    assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk>blk_thresh) {
	
	if(fmv2!=NULL)
		assert(fmv2->child);

    cx = x;
    cy = y;
    child_cand_encode(mapbit, fmv1_array, fmv1->child0, fmv2_array, fmv2 ? fmv2->child0 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_cand_encode(mapbit, fmv1_array, fmv1->child1, fmv2_array, fmv2 ? fmv2->child1 : NULL, 
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_cand_encode(mapbit, fmv1_array, fmv1->child2, fmv2_array, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_cand_encode(mapbit, fmv1_array, fmv1->child3, fmv2_array, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, encode_parallelmv, t_level, blk_thresh, count);
  } else {
    if( x >= hor || y >= ver )    return;

	if (!fmv1->child)  // there are no children for this block with size>=blk_thresh
	{
		if(fmv2 != NULL)
			assert(!fmv2->child);

		//We must decide lifting mode of MERGE blocks in ME stage, for subbands will be generated before MV encoding stage
		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   do_fmv1 = NO;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) do_fmv1 = NO; 
#endif 
		if (fmv1->lifting_mode == UNDECIDED)	assert(0);

		if(fmv2 != NULL){
			// no motion vector for this block on this side 
			if ((fmv2->lifting_mode == IGNORED))   do_fmv2 = NO;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
			if (fmv2->bi_mode==DIRECTIONAL_IBLOCK) do_fmv2 = NO; 
#endif 
			if (fmv2->lifting_mode == UNDECIDED)	assert(0);
		}
	}else
		assert(0);

	///////////////	MERGE detection	///////////////////////
	if( fmv1->bi_mode == BLOCK_MERGING || ( fmv1->bi_mode>=9 && fmv1->direct_idx == INDIRECT && 
		fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == DIR) ){

		assert(do_fmv1 == YES || do_fmv2 == YES);

		if(fmv2 != NULL)
			assert(fmv1->med_idx == fmv2->med_idx);

		for(i=0;i<=3;i++){
			mrg_left[i] = new vector;
			mrg_right[i] = new vector;

			clean_mrg_mv(mrg_left[i]);
			clean_mrg_mv(mrg_right[i]);
		}

		if(fmv1->aff_mrg == YES){
//			printf("Aff MRG blk: cx = %d, cy = %d, xblk = %d, yblk = %d\n",x,y,xblk2,yblk2);
		}

		if(fmv2 != NULL)
			get_field_merge_mv_info(mrg_left, mrg_right,x,y,xblk2,yblk2,hor,ver,info,t_level,frame_motion_field, frame_motion_field2, prev_frame_motion_field2, prev_frame_motion_field1,0);
		else
			get_field_merge_mv_info(mrg_left, mrg_right,x,y,xblk2,yblk2,hor,ver,info,t_level,frame_motion_field, frame_motion_field2, prev_frame_motion_field2, prev_frame_motion_field1,1);
		
		for(i=0;i<=3;i++){

			if( (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL) 
				|| (mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL) ){
				mrg_num ++;
			}
		}
//		printf("mrg_num = %d\n",mrg_num);
		if(mrg_num == 0)
			assert(0);

		///////////////////	MRG idx enc	//////////////////
		if(mrg_num == 0)
			assert(fmv1->med_idx == -1);
		else if(mrg_num == 1)
			assert(fmv1->med_idx == 0);
		else if( mrg_num == 2){
			putbits(fmv1->med_idx,1);
			*mapbit += 1;
			idx_info += 1;
		}else{
			assert(mrg_num == 3 || mrg_num == 4);
			if(mrg_num == 3)
				output_huff_bits(fmv1,mrg_num,cnt_tri,info,mapbit);
			else
				output_huff_bits(fmv1,mrg_num,cnt_quad,info,mapbit);
		}
		
		//////////////////////////////////////////////////
		if(fmv1->bi_mode == BLOCK_MERGING){
			if(mrg_left[fmv1->med_idx]->aff1_mvx != (float)HUGE_VAL || mrg_right[fmv1->med_idx]->aff1_mvx != (float)HUGE_VAL){
				assert( (mrg_left[fmv1->med_idx]->aff2_mvx != (float)HUGE_VAL && mrg_left[fmv1->med_idx]->aff3_mvx != (float)HUGE_VAL) ||
						(mrg_right[fmv1->med_idx]->aff2_mvx != (float)HUGE_VAL && mrg_right[fmv1->med_idx]->aff3_mvx != (float)HUGE_VAL) );
				aff_mrg_blk ++;

				assert(fmv1->aff_mrg == 0 || fmv1->aff_mrg == 1);
				putbits(fmv1->aff_mrg,1);
				*mapbit += 1;
				idx_info += 1;
			}
		}
		///////////////////////////////////////////////////

	}//if BLOCK_MERGING
	///////////////////////////////////////////////////////

  if(do_fmv1 == YES){
	////////////  Added by Yuan Liu on 01.23.2016  //////////////
	for(i=0;i<=3;i++){
		mvpredstr_x[i] = (float)HUGE_VAL;
	    mvpredstr_y[i] = (float)HUGE_VAL;
	}
	/////////////////////////////////////////////////////////////

	map_side = encode_parallelmv;
	assert(map_side == 1);

	if( ( fmv1->bi_mode <= 8 && !(fmv1->bi_mode == PARALLEL && !map_side) && fmv1->bi_mode != BLOCK_MERGING ) || 
		(fmv1->bi_mode >= 9 && fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == INDIR) ){
		get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, 0);
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL){
				num++;
				if(fmv1->bi_mode == BLOCK_MERGING){
				  if( (float(x)-mvpredstr_x[i] >= 0) && (float(x)-mvpredstr_x[i]<= (hor - xblk2)) && (float(y)-mvpredstr_y[i] >= 0)
					  && (float(y) - mvpredstr_y[i] <= (ver - yblk2)) )
						num_dir++;
				}
			}
		}

		if(num == 0)
			assert(fmv1->med_idx == -1);
		else if(num == 1){
				assert(fmv1->med_idx == 0);
		}
		else if(num == 2){
			assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 1);
			putbits(fmv1->med_idx,1);
			*mapbit += 1;
			idx_info += 1;

		}
		else{
			assert(num == 3 || num == 4);
			assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3);
			if(num == 3)
				output_huff_bits(fmv1,num,cnt_tri,info,mapbit);
			else
				output_huff_bits(fmv1,num,cnt_quad,info,mapbit);
			idx_info += 2;
		}

	}
	
	if( fmv1->bi_mode >= 9 ){
		middle = put_aff_mvs(fmv1,x,y,xblk2,yblk2,hor,ver,frame_motion_field, prev_frame_motion_field2,1);
		(*mapbit) += middle;
	}
	/////////////////////////////////////////////////////////////

	fmv1->is_predictor = YES;

//	printf("\nfmv->bi_mode=%d, fmv->med_idx=%d", fmv1->bi_mode,fmv1->med_idx);

	//////////////	Added by Yuan Liu	//////////////////
	if(fmv1->bi_mode >= 7 && fmv1->bi_mode <= 11){
		if( (fmv1->bi_mode == 7 && fmv1->aff_mrg == YES) || fmv1->bi_mode >= 9 )
			printf("\nbi_mode = %d, sad_cost = %f, skip = %d, med_idx = %d, x = %d, y = %d, xblk = %d, yblk = %d\nfmv1->mvx = %f, fmv1->mvy = %f, fmv1->dmvx = %f, fmv1->dmvy = %f\n"
			,fmv1->bi_mode, fmv1->sad_cost, fmv1->skip_sign, fmv1->med_idx,x,y,xblk2,yblk2,fmv1->mvx,fmv1->mvy,fmv1->dmvx,fmv1->dmvy);

		if(fmv1->bi_mode>=9  && fmv1->direct_idx == DIRECT){
			printf("aff_index = %d\nfmv1->aff1_mvx = %f, fmv1->aff1_mvy = %f\n fmv1->aff2_mvx = %f, fmv1->aff2_mvy = %f\n fmv1->aff3_mvx = %f, fmv1->aff3_mvy = %f\n"
			,fmv1->aff_idx,fmv1->aff1_mvx, fmv1->aff1_mvy, fmv1->aff2_mvx , fmv1->aff2_mvy, fmv1->aff3_mvx, fmv1->aff3_mvy);
		}
		else if(fmv1->bi_mode>=9  && fmv1->direct_idx != DIRECT){
			printf("aff_index = %d\nfmv1->aff1_mvx = %f, fmv1->aff1_mvy = %f\n fmv1->aff2_mvx = %f, fmv1->aff2_mvy = %f\n fmv1->aff3_mvx = %f, fmv1->aff3_mvy = %f\nmerge_idx = %d, merge_dir = %d\n"
			,fmv1->aff_idx,fmv1->aff1_mvx, fmv1->aff1_mvy, fmv1->aff2_mvx , fmv1->aff2_mvy, fmv1->aff3_mvx, fmv1->aff3_mvy,fmv1->merge_idx,fmv1->merge_dir);
			if(fmv1->merge_idx == MERGE && fmv1->merge_dir == UP)
				printf("fmv1->aff3_dmvx = %f, fmv1->aff3_dmvy = %f\n", fmv1->aff3_dmvx, fmv1->aff3_dmvy);
			else if(fmv1->merge_idx == MERGE && fmv1->merge_dir == LEFT)
				printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n", fmv1->aff2_dmvx, fmv1->aff2_dmvy);
			else if(fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && fmv1->merge_dir == PAL_L && fmv1->aff_idx >=0)
				|| (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
				printf("\nfmv1->aff1_dmvx = %f, fmv1->aff1_dmvy = %f\n", fmv1->aff1_dmvx, fmv1->aff1_dmvy);
				printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n", fmv1->aff2_dmvx, fmv1->aff2_dmvy);
				printf("fmv1->aff3_dmvx = %f, fmv1->aff3_dmvy = %f\n\n", fmv1->aff3_dmvx, fmv1->aff3_dmvy);
			}
		}
		if( (fmv1->bi_mode == 7 && fmv1->aff_mrg == YES) || fmv1->bi_mode >= 9 )
			printf("\n");
	}
	//////////////////////////////////////////////////////

	if(num == 2 && fmv1->bi_mode != BLOCK_MERGING && fmv1->bi_mode <= 8 )
		assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 1);

	if( fmv1->bi_mode <= 6 || fmv1->bi_mode == 8 || (fmv1->bi_mode == 7 && fmv1->aff_mrg == NO) ){
		assert( (x + xblk2 - fmv1->mvx <= hor ) && (y + yblk2 - fmv1->mvy <= ver ) &&
			(x - fmv1->mvx >= 0 ) && (y - fmv1->mvy >= 0 ) );
		assert( (!(fmv1->bi_mode == PARALLEL)) ||  
			((x + xblk2 - fmv2->mvx <= hor ) && (y + yblk2 - fmv2->mvy <= ver ) &&  
		(x - fmv2->mvx >= 0 ) &&  (y - fmv2->mvy >= 0 )) );
	}
	// update the motion field in this frame 
	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
							  fmv1->dmvy, fmv1->mvx, fmv1->mvy);
  }//if do_fmv1

/*************************************************************************/
/**************************	  fmv2	**************************************/
/*************************************************************************/
if(do_fmv2 == YES){
  if(fmv2 != NULL){
  ////////////  Added by Yuan Liu on 08.20.2016  //////////////
	for(i=0;i<=3;i++){
		mvpredstr_x[i] = (float)HUGE_VAL;
	    mvpredstr_y[i] = (float)HUGE_VAL;
	}
	/////////////////////////////////////////////////////////////

	map_side = (encode_parallelmv + 1) % 2;
	assert(map_side == 0);
	num = 0;
	num_dir = 0;

	///////////////// Added on 01.28.2016  //////////////////////
	if( (fmv2->bi_mode <= 8 && !(fmv2->bi_mode == PARALLEL && !map_side) && fmv2->bi_mode != BLOCK_MERGING) ||
		(fmv2->bi_mode >= 9 && fmv2->direct_idx == INDIRECT && fmv2->merge_idx == MERGE && fmv2->merge_dir == TRAN_P && fmv2->trans_pred_idx == INDIR) ){
		get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field2, prev_frame_motion_field1,
									  x, y, xblk, yblk, info, t_level, 0);
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL){
				num++;
				if(fmv2->bi_mode == BLOCK_MERGING){
				  if( (float(x)-mvpredstr_x[i] >= 0) && (float(x)-mvpredstr_x[i]<= (hor - xblk2)) && (float(y)-mvpredstr_y[i] >= 0)
					  && (float(y) - mvpredstr_y[i] <= (ver - yblk2)) )
						num_dir++;
				}
			}
		}
		
//		assert( fmv2->bi_mode >= 0 && fmv2->bi_mode <= 8 );

		if(num == 0)
			assert(fmv2->med_idx == -1);
		else if(num == 1){
				assert(fmv2->med_idx == 0);
		}
		else if(num == 2){
			assert(fmv2->med_idx >= 0 && fmv2->med_idx <= 1);
			putbits(fmv2->med_idx,1);
			*mapbit += 1;
			idx_info += 1;
		}
		else{
			assert(num == 3 || num == 4);
			assert(fmv2->med_idx >= 0 && fmv2->med_idx <= 3);
			if(num == 3)
				output_huff_bits(fmv2,num,cnt_tri,info,mapbit);
			else
				output_huff_bits(fmv2,num,cnt_quad,info,mapbit);
			idx_info += 2;
		}

	}else if( fmv2->bi_mode >= 9 ){
		if(fmv2->direct_idx == INDIRECT && fmv2->merge_idx == MERGE && fmv2->merge_dir == PAL_L)
			middle = put_aff_mvs(fmv2,x,y,xblk2,yblk2,hor,ver,frame_motion_field2, prev_frame_motion_field1,0);
		else
			middle = put_aff_mvs(fmv2,x,y,xblk2,yblk2,hor,ver,frame_motion_field2, prev_frame_motion_field1,1);

		(*mapbit) += middle;
	}
	/////////////////////////////////////////////////////////////

	fmv2->is_predictor = YES;

	//////////////	Added by Yuan Liu	//////////////////
	if(fmv2->bi_mode >= 7 && fmv2->bi_mode <= 11){
		if( (fmv2->bi_mode == 7 && fmv2->aff_mrg == YES) || fmv2->bi_mode >= 9 )
			printf("\nbi_mode = %d, sad_cost = %f, skip = %d, med_idx = %d, x = %d, y = %d, xblk = %d, yblk = %d\nfmv2->mvx = %f, fmv2->mvy = %f, fmv2->dmvx = %f, fmv2->dmvy = %f\n"
			,fmv2->bi_mode, fmv2->sad_cost, fmv2->skip_sign, fmv2->med_idx,x,y,xblk2,yblk2,fmv2->mvx,fmv2->mvy,fmv2->dmvx,fmv2->dmvy);

		if(fmv2->bi_mode>=9  && fmv2->direct_idx == DIRECT){
			printf("aff_index = %d\nfmv2->aff1_mvx = %f, fmv2->aff1_mvy = %f\n fmv2->aff2_mvx = %f, fmv2->aff2_mvy = %f\n fmv2->aff3_mvx = %f, fmv2->aff3_mvy = %f\n"
			,fmv2->aff_idx,fmv2->aff1_mvx, fmv2->aff1_mvy, fmv2->aff2_mvx , fmv2->aff2_mvy, fmv2->aff3_mvx, fmv2->aff3_mvy);
		}
		else if(fmv2->bi_mode>=9  && fmv2->direct_idx != DIRECT){
			printf("aff_index = %d\nfmv2->aff1_mvx = %f, fmv2->aff1_mvy = %f\n fmv2->aff2_mvx = %f, fmv2->aff2_mvy = %f\n fmv2->aff3_mvx = %f, fmv2->aff3_mvy = %f\nmerge_idx = %d, merge_dir = %d\n"
			,fmv2->aff_idx,fmv2->aff1_mvx, fmv2->aff1_mvy, fmv2->aff2_mvx , fmv2->aff2_mvy, fmv2->aff3_mvx, fmv2->aff3_mvy,fmv2->merge_idx,fmv2->merge_dir);
			if(fmv2->merge_idx == MERGE && fmv2->merge_dir == UP)
				printf("fmv2->aff3_dmvx = %f, fmv2->aff3_dmvy = %f\n", fmv2->aff3_dmvx, fmv2->aff3_dmvy);
			else if(fmv2->merge_idx == MERGE && fmv2->merge_dir == LEFT)
				printf("fmv2->aff2_dmvx = %f, fmv2->aff2_dmvy = %f\n", fmv2->aff2_dmvx, fmv2->aff2_dmvy);
			else if(fmv2->merge_idx == INTER || (fmv2->merge_idx == MERGE && fmv2->merge_dir == PAL_L && fmv2->aff_idx >=0)
				|| (fmv2->merge_idx == MERGE && fmv2->merge_dir == TRAN_P) ){
				printf("\nfmv2->aff1_dmvx = %f, fmv2->aff1_dmvy = %f\n", fmv2->aff1_dmvx, fmv2->aff1_dmvy);
				printf("fmv2->aff2_dmvx = %f, fmv2->aff2_dmvy = %f\n", fmv2->aff2_dmvx, fmv2->aff2_dmvy);
				printf("fmv2->aff3_dmvx = %f, fmv2->aff3_dmvy = %f\n\n", fmv2->aff3_dmvx, fmv2->aff3_dmvy);
			}
		}
		if( (fmv2->bi_mode == 7 && fmv2->aff_mrg == YES) || fmv2->bi_mode >= 9 )
			printf("\n\n");
	}
	//////////////////////////////////////////////////////

	if(num == 2 && fmv2->bi_mode != BLOCK_MERGING && fmv2->bi_mode <= 8 )
		assert(fmv2->med_idx >= 0 && fmv2->med_idx <= 1);

	if( fmv2->bi_mode <= 6 || fmv2->bi_mode == 8 || (fmv2->bi_mode == 7 && fmv2->aff_mrg == NO)  ){
		assert( (x + xblk2 - fmv2->mvx <= hor ) && (y + yblk2 - fmv2->mvy <= ver ) &&
			(x - fmv2->mvx >= 0 ) && (y - fmv2->mvy >= 0 ) );
		assert( (!(fmv2->bi_mode == PARALLEL)) ||  
			((x + xblk2 - fmv1->mvx <= hor ) && (y + yblk2 - fmv1->mvy <= ver ) &&  
		(x - fmv1->mvx >= 0 ) &&  (y - fmv1->mvy >= 0 )) );
	}
	// update the motion field in this frame 
	update_frame_motion_field(frame_motion_field2, x, y, xblk, yblk, info, fmv2, fmv2->dmvx, 
							  fmv2->dmvy, fmv2->mvx, fmv2->mvy);

  }//If fmv2 != NULL

}//if do_fmv2

  ///////////////	MERGE detection	//////////////////////
  if( fmv1->bi_mode == BLOCK_MERGING ){
	for(i=0;i<=3;i++){
		delete(mrg_left[i]);
		delete(mrg_right[i]);
    }
  }
  ///////////////////////////////////////////////////////

  assert(fmv1->bi_mode >= 0);
  mean_mode_mse[fmv1->bi_mode] += fmv1->mse;
  mean_mode_num[fmv1->bi_mode] ++;

}//child else

}
///////////////////////////////////////////////////////////////////////////

/*****************************************************************************/
/*                                mv_encode()                                */
/*****************************************************************************/
void
mv_encode( int *mvbit, int *mapbit, vector_ptr fmv1, vector_ptr fmv2,
           int encode_map, videoinfo info, int large, int t_level, int sub_mv[LAYER_NUM][MAX_AGP_LEVEL], 
		   int *sub_sign, int count, int GOP_counter )
{
  int partition, x, y, X, Y, xnum, ynum, xblk, yblk, hor, ver, pos, small, itemp, subpel;
  int num_symbol;
  float pmvx, pmvy;
  int sub_bit, blk_thresh, layer_num; 
  int i,j;

  if(fmv2 != NULL){
    for(i=0;i<=11;i++){
	  bi_mode_num012[i]=0;
	  bi_mode_num345[i]=0;
	}
	use_huff = 0;
  }

  // for motion field in a frame
  frame_motion_field = (FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(FRAME_MOTION_FIELD), "frame_motion_field"); 
  frame_motion_field2 = (FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(FRAME_MOTION_FIELD), "frame_motion_field2"); 

  blk_thresh = 0; 
  if (info.layer_mv[t_level])
	blk_thresh = LAYER_BLOCK_SIZE; 

#ifdef DEBUG_LAYER_STRUCTURE
  FILE *fpbase, *fpenhance; 
  char  base_file[80], enhance_file[80];

  // make base_file and enhance_file empty
  sprintf(base_file, "base%d.txt", count); 
  if (fpbase=fopen(base_file, "rt") )  
  {
	  fclose(fpbase);
	  fpbase=fopen(base_file, "wt");
	  fclose(fpbase); 
  }
  sprintf(enhance_file, "enhance%d.txt", count); 
  if (fpenhance=fopen(enhance_file, "rt") )
  {
	  fclose(fpenhance);
	  fpenhance=fopen(enhance_file, "wt");
	  fclose(fpenhance); 
  }
#endif

  // initialization for scalable motion vector coding: AGP and layer structure 
  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
	for (partition=0; partition<info.AGP_level[t_level]; partition++)
			sub_mv[layer_num][partition] = 0; 
  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
  {
	  for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)   // maximum 3 subsymbols
	  {
		  // for the subsymbol encoding
		  splitted_mv_byte_num_array[layer_num][sub_bit] = 0;
		  splitted_mv_byte_array[layer_num][sub_bit]     = 0;
		  splitted_bit_num_array[layer_num][sub_bit]     = 8;
	  }
	  // initialize for the additional sign
	 splitted_sign_byte_num[layer_num] = 0;
	 splitted_sign_byte[layer_num]     = 0;
	 splitted_sign_bit_num[layer_num]  = 8;
	 sub_sign[layer_num] = 0; 
  }

  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];
  subpel = info.subpel[t_level];
  hor = info.ywidth;
  ver = info.yheight;
  small = xblk;
  itemp = info.level[t_level];
  while( itemp != 1 ) {
    small /= 2;
    itemp--;
  }
  num_symbol = 0;
  
#ifdef  DEBUG_SCALABLE_MV  
  sprintf(encoder_AGPdebug, "encoder_AGP_debugfile_GOP%03d_count%03d.txt", GOP_counter, debug_counter);
  fpAGP_debug = fopen(encoder_AGPdebug, "wt");
  debug_counter++;
#endif 

  *mapbit = 0;                  /* map encoding */
#ifdef EC_USE_CONTEXTS
  ec_enc_init(info, EC_TYPE, num_symbol, 3, 1);
#else
  ec_enc_init(info, EC_TYPE, num_symbol, 1, 1);
#endif

//////////////////////////////////
  if( encode_map ) {
    for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {     /* coding loop */
      for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
        pos = Y * xnum + X;
        child_map_pre_encode( &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, 
                          x, y, xblk, yblk, hor, ver, small, info, t_level );
      }
    }
  }
//////////////////////////////////
  if(fmv2 != NULL){
	  if( bi_mode_num012[1] + bi_mode_num012[2] < bi_mode_num012[7] ){
		  printf("bi_mode_num012[1] = %d, bi_mode_num012[2] = %d, bi_mode_num012[7] = %d\n",bi_mode_num012[1],bi_mode_num012[2],bi_mode_num012[7]);
		  use_huff = 1;
	  }
	  putbits(use_huff,1);
	  ( *mapbit )++; 
  }

  if( encode_map ) {
    for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {     /* coding loop */
      for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
        pos = Y * xnum + X;
        child_map_encode( &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, 
                          mapbit, x, y, xblk, yblk, hor, ver, small, info, t_level );
      }
    }
  }

  // encode the merge_sign bit
  if( info.layer_mv[t_level] ) {
    for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {     /* coding loop */
      for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
        pos = Y * xnum + X;
        child_merge_encode( &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, 
                          mapbit, x, y, xblk, yblk, hor, ver, small, info, t_level, 
						  blk_thresh, count );
      }
    }
  }

  clear_frame_motion_field(frame_motion_field, info); 
  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {       /* mv coding */
    pmvx = 0.;
    pmvy = 0.;                  /* initialize to the zero per each row */

    for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
      pos = Y * xnum + X;
      child_mv_encode( fmv1, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, 
                       &pmvx, &pmvy, num_symbol, subpel, x, y, xblk, yblk,
                       hor, ver, info, encode_map, t_level, blk_thresh, count ); // add encode_map for mv coding in parallel mode. mwi 
    }
  }
  ec_enc_end1();

  if( (fmv2 != NULL && encode_map == 1) ){
	copy_frame_motion_field2(prev_frame_motion_field2, prev_frame_motion_field_left2, info);
	copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field_left1, info);
  }

  if( (fmv2 != NULL && encode_map == 0) ){
	copy_frame_motion_field2(prev_frame_motion_field_left2, prev_frame_motion_field2, info);
	copy_frame_motion_field2(prev_frame_motion_field_left1, prev_frame_motion_field1, info);
  }

  if( fmv2 != NULL && encode_map == 0 ){
	  putbits(0,14);
	  *mapbit += 14;
	  idx_info += 14;

	  if( fmv2 != NULL && encode_map == 0 ){

		clear_frame_motion_field(frame_motion_field, info);//This is to confirm that neighbor prediction candidates will be encoded in the same context order as
														 //MV residuals.
		clear_frame_motion_field(frame_motion_field2, info);

	    for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {       /* mv coding */
		  pmvx = 0.;
		  pmvy = 0.;                  /* initialize to the zero per each row */

		  for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
		    pos = Y * xnum + X;
		    child_cand_encode( mapbit, fmv2, &fmv2[pos], fmv1, fmv1 ? (&fmv1[pos]) : NULL, 
						   &pmvx, &pmvy, num_symbol, subpel, x, y, xblk, yblk,
						   hor, ver, info, 1, t_level, blk_thresh, count );
		  }
	    } 

		copy_frame_motion_field(frame_motion_field, prev_frame_motion_field2, info);
		copy_frame_motion_field(frame_motion_field2, prev_frame_motion_field1, info);

	  }

  }//(fmv2 != NULL && encode_map == 0) 
  else if(fmv2 == NULL){
	  putbits(0,14);
	  *mapbit += 14;
	  idx_info += 14;

	  clear_frame_motion_field(frame_motion_field2, info);

	  clear_frame_motion_field(frame_motion_field, info);//This is to confirm that neighbor prediction candidates will be encoded in the same context order as
														 //MV residuals.
	  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {       /* mv coding */
		pmvx = 0.;
		pmvy = 0.;                  /* initialize to the zero per each row */

		for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
		  pos = Y * xnum + X;
		  child_cand_encode( mapbit, fmv1, &fmv1[pos], fmv2, fmv2 ? (&fmv2[pos]) : NULL, 
						   &pmvx, &pmvy, num_symbol, subpel, x, y, xblk, yblk,
						   hor, ver, info, 1, t_level, blk_thresh, count );
		}
	  }

	  copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
	  copy_frame_motion_field(frame_motion_field, prev_frame_motion_field1, info);
  
  }//if fmv2 == NULL
  else{
	assert( fmv2 != NULL && encode_map == 1 );
	
	copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
	copy_frame_motion_field(frame_motion_field, prev_frame_motion_field1, info);
  }

  ec_enc_end2();

  mvbit[0] = 8 * outbyte - *mapbit;   // base layer motion vector bytes 

  if (info.AGP_level[t_level])
  {
	  // additional sign bits 
	  if (splitted_sign_bit_num[0] != 8) 
	  {
		  splitted_sign_byte[0] <<= splitted_sign_bit_num[0];
		  store_splitted_sign[0][splitted_sign_byte_num[0]++] = splitted_sign_byte[0];
	  }
	  done_splitted_bytes(store_splitted_sign[0], splitted_sign_byte_num[0], info, t_level);
	  // 2 bytes for the length of splitted additional sign
	  sub_sign[0] = 8* (splitted_sign_byte_num[0]+2); 

	  // save the sub-symbol bits in the order of significance 
	  for (sub_bit=info.AGP_level[t_level]-1; sub_bit>=0; sub_bit--)  
	  {
		  if (splitted_bit_num_array[0][sub_bit] != 8)
		  {
			  splitted_mv_byte_array[0][sub_bit] <<= splitted_bit_num_array[0][sub_bit];
			  store_splitted_bytes_array[0][sub_bit][splitted_mv_byte_num_array[0][sub_bit]++] = 
				                       splitted_mv_byte_array[0][sub_bit];
		  }
		  done_splitted_bytes(store_splitted_bytes_array[0][sub_bit], 
			                  splitted_mv_byte_num_array[0][sub_bit], 
			                  info, t_level);
		  // 2 bytes for the length of splitted motion vector bit-stream
		  sub_mv[0][sub_bit] = 8* (splitted_mv_byte_num_array[0][sub_bit]+2); 
	  }
  }


  if (info.layer_mv[t_level]) // layer structure: enhancement layer
  {
	  assert(0);
	  ec_enc_init(info, EC_TYPE, num_symbol, 3, 1);
	  clear_frame_motion_field(frame_motion_field, info); 
	  blk_thresh = LAYER_BLOCK_SIZE; 
	  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {      
		for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
		  pos = Y * xnum + X;
		  child_mv_encode_enhance( fmv1, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, 
						   &pmvx, &pmvy, num_symbol, subpel, x, y, xblk, yblk,
						   hor, ver, info, encode_map, t_level, blk_thresh, count ); 
		}
	  }
	  ec_enc_end1();
	  ec_enc_end2();
	 
	  mvbit[1] = 8 * outbyte;  // enhancement layer bytes
	  if (info.AGP_level[t_level])
	  {
		  // additional sign bits 
		  if (splitted_sign_bit_num[1] != 8) 
		  {
			  splitted_sign_byte[1] <<= splitted_sign_bit_num[1];
			  store_splitted_sign[1][splitted_sign_byte_num[1]++] = splitted_sign_byte[1];
		  }
		  done_splitted_bytes(store_splitted_sign[1], splitted_sign_byte_num[1], info, t_level);
		  // 2 bytes for the length of splitted additional sign
		  sub_sign[1] = 8* (splitted_sign_byte_num[1]+2); 
		  // save the sub-symbol bits in the order of significance 
		  for (sub_bit=info.AGP_level[t_level]-1; sub_bit>=0; sub_bit--)  
		  {
			  if (splitted_bit_num_array[1][sub_bit] != 8)
			  {
				  splitted_mv_byte_array[1][sub_bit] <<= splitted_bit_num_array[1][sub_bit];
				  store_splitted_bytes_array[1][sub_bit][splitted_mv_byte_num_array[1][sub_bit]++] = 
										   splitted_mv_byte_array[1][sub_bit];
			  }
			  done_splitted_bytes(store_splitted_bytes_array[1][sub_bit], 
								  splitted_mv_byte_num_array[1][sub_bit], 
								  info, t_level);
			  // 2 bytes for the length of splitted motion vector bit-stream
			  sub_mv[1][sub_bit] = 8* (splitted_mv_byte_num_array[1][sub_bit]+2); 
		  }
	  }

  }

#ifdef DEBUG_SCALABLE_MV
  fclose(fpAGP_debug); 
#endif   

  free(frame_motion_field); 
  free(frame_motion_field2); 

}

void write_mv_length(int AGP_level, int layer_mv, FILE *fp_mv, Rate FrsRate, 
					 int count, int sub_sign[LAYER_NUM], int itmp, int simul_enc)
{
	int partition, layer_num; 

	if (AGP_level && !layer_mv )  // AGP, no layer structure 
	{
		// major part: mapbits + mvbits (bytes)
		fprintf(fp_mv, "%d  ", (FrsRate.map[count] + FrsRate.mv[0][count])>>3); 
		// additional sign part
		fprintf(fp_mv, "%d  ", sub_sign[0]>>3);
		// sub-symbol part (bytes)
		for (partition=0; partition<AGP_level; partition++)
			fprintf(fp_mv, "%d  ",  FrsRate.submv[0][partition][count] >>3 );
	}

	if (!AGP_level && layer_mv ) // no AGP, layer structure 
	{
		// major part: mapbits + mvbits (bytes)
		fprintf(fp_mv, "%d  ", (FrsRate.map[count] + FrsRate.mv[0][count])>>3); 
		for (layer_num=1; layer_num<LAYER_NUM; layer_num++)
			fprintf(fp_mv, "%d  ", ( FrsRate.mv[layer_num][count])>>3); 
	}

	if (AGP_level && layer_mv ) // AGP,    layer structure 
	{
		// major part: mapbits + mvbits (bytes)
		fprintf(fp_mv, "%d  ", (FrsRate.map[count] + FrsRate.mv[0][count])>>3); 
		// additional sign part
		fprintf(fp_mv, "%d  ", sub_sign[0]>>3);
		// sub-symbol part (bytes)
		for (partition=0; partition<AGP_level; partition++)
			fprintf(fp_mv, "%d  ",  FrsRate.submv[0][partition][count] >>3 );
		for (layer_num=1; layer_num<LAYER_NUM; layer_num++)
		{
			// major part: mapbits + mvbits (bytes)
			fprintf(fp_mv, "%d  ", (FrsRate.mv[layer_num][count])>>3); 
			// additional sign part
			fprintf(fp_mv, "%d  ", sub_sign[layer_num]>>3);
			// sub-symbol part (bytes)
			for (partition=0; partition<AGP_level; partition++)
				fprintf(fp_mv, "%d  ",  FrsRate.submv[layer_num][partition][count] >>3 );
		}
	}

	// total sum for this set of motion vectors
//	if(simul_enc == NO)
		fprintf( fp_mv, "%d\n", itmp >> 3);
}



/*****************************************************************************/
/*                                mv_encoding()                              */
/*****************************************************************************/
int
mv_encoding( videoinfo info, Rate FrsRate, vector_ptr * yfmv, int GOP_counter, int simul_enc, int curr )
{
  int partition, i, j, sum, itmp, count, mv_count, dist, GOPsz;
  int eff_GOPsz[20]; // maximum number of temporal levels = 20;
  int mvrate1[LAYER_NUM], maprate1;
  int mvrate2[LAYER_NUM], maprate2;
  // for scalable motion vector coding 
  int sub_mv1[LAYER_NUM][MAX_AGP_LEVEL], sub_mv2[LAYER_NUM][MAX_AGP_LEVEL]; 
  int sub_sign1[LAYER_NUM], sub_sign2[LAYER_NUM]; 
  int layer_num; 

  int a,b;
  int res_frame, do_res = NO;
  int *level_res_frame;
  float getval;

  for(i = 0; i < NUMBER_OF_BI_MODES; i ++){
	  mean_mode_mse[i] = 0;
	  mean_mode_num[i] = 0;
  }

  level_res_frame = (int*)getarray(info.tPyrLev, sizeof(int), "level_res_frame");

  res_frame = info.act_last - curr + 1;
  if( res_frame > 0 && res_frame < info.GOPsz ){
	do_res = YES;
	for(i = 0; i < info.tPyrLev; i ++){
		level_res_frame[i] = res_frame;

		if(res_frame % 2 == 0)
			res_frame = res_frame / 2;
		else
			res_frame = (res_frame+1) / 2;
	}
  }else{
	res_frame = info.GOPsz;
	for(i = 0; i < info.tPyrLev; i ++)
		level_res_frame[i] = res_frame;
  }

  ctl_info = 0;
  mv_info = 0;
  idx_info = 0;

  aff_mrg_blk = 0;

  mv_res_bits = 0;
  calc_res_bits = 0;

  printf("GOP_counter = %d\n",GOP_counter);

//  test();

  for(i=0;i<=2;i++){
	cnt_tri[i]=0;
  }
  for(i=0;i<=3;i++){
	cnt_quad[i]=0;
  }
  // enum FLAG Level_change;  
  FILE *fpstat, *fp_mv;

  prev_frame_motion_field1 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field1"); 
  prev_frame_motion_field2 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field2"); 

  //Added on 08.20.2016
  prev_frame_motion_field_left1 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field_left1"); 
  prev_frame_motion_field_left2 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field_left2"); 
  //Added on 08.20.2016

  if(GOP_counter == 0 && simul_enc == NO){
	  for(i = 0;i <= info.tPyrLev - 2; i++){
		buffer_frame_motion_field1[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "buffer_frame_motion_field1");
		buffer_frame_motion_field2[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "buffer_frame_motion_field2");

		clear_frame_motion_field_simp(buffer_frame_motion_field1[i], info);
		clear_frame_motion_field_simp(buffer_frame_motion_field2[i], info);
//Added on 02.05.2018
		save_buffer_frame_motion_field1[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "save_buffer_frame_motion_field1");
		save_buffer_frame_motion_field2[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "save_buffer_frame_motion_field2");

		clear_frame_motion_field_simp(save_buffer_frame_motion_field1[i], info);
		clear_frame_motion_field_simp(save_buffer_frame_motion_field2[i], info);
	  }
  }

  if(simul_enc == YES){
	  for(i = 0; i <= info.tPyrLev - 2; i++){
		copy_frame_motion_field2(save_buffer_frame_motion_field1[i], buffer_frame_motion_field1[i], info);
		copy_frame_motion_field2(save_buffer_frame_motion_field2[i], buffer_frame_motion_field2[i], info);
	  }
  }

  if((info.GOPsz * GOP_counter) == SIMUL_POINT){
	  printf("\nMV_coding SIMUL_POINT = %d\n", SIMUL_POINT);

	  for(i = 0; i <= info.tPyrLev - 2; i++){
		copy_frame_motion_field2(buffer_frame_motion_field1[i], save_buffer_frame_motion_field1[i], info);
		copy_frame_motion_field2(buffer_frame_motion_field2[i], save_buffer_frame_motion_field2[i], info);
	  }
  }
 
  GOPsz = info.GOPsz;
  
  dist = ( int )pow( 2.0, ( double )( info.tPyrLev - 1 ) );
  
  fpstat = fopen( info.statname, "at+" );
  if( !( fp_mv = fopen( info.mvstatname, "at+" ) ) ) {
    printf( "can not open %s\n", info.mvstatname );
    exit( 1 );
  }
  
  // determine effective GOP size in level i
  for( i = 0; i < info.tPyrLev; i++ ) {
    eff_GOPsz[i] = GOPsz;  
    GOPsz /= 2;  
  }
  
  FrsRate.map[0] = 0;
  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
  {
	  FrsRate.mv[layer_num][0] = 0; 
	  // for scalable motion vector coding 
	  for (partition=0; partition<info.AGP_level[info.tPyrLev - 1]; partition++)
			FrsRate.submv[layer_num][partition][0]=0; 
  }

  mv_count = 0;
  count = 1;    // start with yfmv[2]
  sum = 0; 

  for( i = info.tPyrLev - 1; i >= 0; i-- ) {
    FrsRate.map[count] = 0;
    for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
	{
		FrsRate.mv[layer_num][count] = 0; 
		// for scalable motion vector coding: AGP 
		for (partition=0; partition<info.AGP_level[i]; partition++)
			FrsRate.submv[layer_num][partition][count]=0; 
	}

    fprintf( fpstat, "Fr%.2d: no motion information needed\n", count );  

//	if(simul_enc == NO)
		fprintf( fp_mv, "%d\n", 0);
	// the first set of motion vectors has already beend coded in previous GOP 
	// so skip it by count++,   Yongjun Wu
	// that means in level 3 motion vector sets: 1, 2 , 3          1  is skipped
	//            in level 2 motion vector sets: 4, 5 , 6, 7, 8    4  is skipped
	//            in level 1 motion vector sets: 9, 10,...,   17   9  is skipped
	//            in level 0 motion vector sets: 18,19,...,   34   18 is skipped 
    count++;  

    // scene_change[i][j = 0] is always YES!

	printf("LEVEL CHANGED!\n");
	clear_frame_motion_field_simp(prev_frame_motion_field1, info);
	clear_frame_motion_field_simp(prev_frame_motion_field2, info);

	//Added on 08.20.2016
	clear_frame_motion_field_simp(prev_frame_motion_field_left1, info);
	clear_frame_motion_field_simp(prev_frame_motion_field_left2, info);
	//Added on 08.20.2016

	if(GOP_counter >= 1 && i <= info.tPyrLev - 2){
		copy_frame_motion_field2(buffer_frame_motion_field1[i], prev_frame_motion_field1, info);
		copy_frame_motion_field2(buffer_frame_motion_field2[i], prev_frame_motion_field2, info);
	}

	avg_blk_size[i] = 0;
	blk_num = 0;

    for( j = 1; j <= eff_GOPsz[i]; j += 2 ) {

	  // initialization
	  maprate1 =  maprate2 =  0 ; 
	  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
	  {
		  mvrate1[layer_num] = sub_sign1[layer_num] =  0;
		  mvrate2[layer_num] = sub_sign2[layer_num] =  0;
		  for (partition=0; partition<info.AGP_level[i]; partition++)
			sub_mv1[layer_num][partition] = sub_mv2[layer_num][partition] = 0; 
	  }

	  calc_res_sad = 0;
		  
      if (j == eff_GOPsz[i]) { // single MVF
		printf("single MVF happened!\n");
		assert(0);
        if (scene_change[i][j] == NO) {
#ifdef  DEBUG_BLOCK_MODE_MV_INFO
		  write_block_mode_motion_vector("encoderleft", GOP_counter, count, 1, yfmv[count], info, i);
#endif 
		  // subsample motion vector according to block size
		  layer_structure_mv_subsample(info, i, j,  yfmv[count], 1); 
          mvStat_setFrame(GOP_counter, i, j);
          mv_encode(mvrate1, &maprate1, yfmv[count], NULL, 1, info, dist, i, sub_mv1, sub_sign1, 
					count, GOP_counter);
        } else {
		  maprate1 = 0;
		  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
		  {
			  mvrate1[layer_num] = sub_sign1[layer_num] =  0;
			  for (partition=0; partition<info.AGP_level[i]; partition++)
				sub_mv1[layer_num][partition] = 0; 
		  }
        }
      } else if (scene_change[i][j] == NO && scene_change[i][j + 1] == NO) {
		printf("\n\nBI MC, i = %d, j = %d, eff_GOPsz = %d, count = %d, count+1 = %d\n",i,j,eff_GOPsz[i], count,count+1);
//		printf("eff_GOPsz = %d\n",eff_GOPsz[i]);

#ifdef  DEBUG_BLOCK_MODE_MV_INFO
		write_block_mode_motion_vector("encoderleft", GOP_counter, count, 0, yfmv[count], info, i);
		write_block_mode_motion_vector("encoderright",GOP_counter, count+1, 0, yfmv[count+1], info, i);
#endif
 	    // subsample motion vector according to block size
	    layer_structure_mv_subsample(info, i, j,  yfmv[count], 1); 
        mvStat_setFrame(GOP_counter, i, j);
        mv_encode(mvrate1, &maprate1, yfmv[count], yfmv[count + 1], 1, 
                  info, dist, i, sub_mv1, sub_sign1, count, GOP_counter);
 	    // subsample motion vector according to block size
	    layer_structure_mv_subsample(info, i, j+1,  yfmv[count + 1], 1); 
        mvStat_setFrame(GOP_counter, i, j + 1);

		printf("\n\nBI MC, i = %d, j = %d, eff_GOPsz = %d, count+1 = %d, count = %d\n",i,j,eff_GOPsz[i], count+1,count);
//		printf("eff_GOPsz = %d\n",eff_GOPsz[i]);

        mv_encode(mvrate2, &maprate2, yfmv[count + 1], yfmv[count], 0, 
                  info, dist, i, sub_mv2, sub_sign2, count+1, GOP_counter);

      } else if (scene_change[i][j] == NO) {
		printf("\n\nLEFT MC, i = %d, j = %d, eff_GOPsz = %d, count = %d\n",i,j,eff_GOPsz[i], count);
//		printf("eff_GOPsz = %d\n",eff_GOPsz[i]);

#ifdef  DEBUG_BLOCK_MODE_MV_INFO
		write_block_mode_motion_vector("encoderleft", GOP_counter, count, 1, yfmv[count], info, i);
#endif 
 	    // subsample motion vector according to block size
	    layer_structure_mv_subsample(info, i, j,  yfmv[count], 1); 
        mvStat_setFrame(GOP_counter, i, j);

        mv_encode(mvrate1, &maprate1, yfmv[count], NULL, 1, 
                  info, dist, i, sub_mv1, sub_sign1, count, GOP_counter);
		maprate2 = 0;
		for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
		{
			mvrate2[layer_num] = sub_sign2[layer_num] =  0;
			for (partition=0; partition<info.AGP_level[i]; partition++)
				sub_mv2[layer_num][partition] = 0; 
		}
		/////////////// Added by Yuan Liu on 04.21.2016	///////////
//		if( GOP_counter != (get_GOP_num(info) - 1) ){
			copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
			clear_frame_motion_field_simp(prev_frame_motion_field1, info);
//		}
		///////////////////////////////////////////////////////////

      } else if (scene_change[i][j + 1] == NO) {

		/////////////// Added by Yuan Liu on 04.21.2016	///////////
//		if( GOP_counter != (get_GOP_num(info) - 1) ){
			copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
			clear_frame_motion_field_simp(prev_frame_motion_field1, info);
//		}
		///////////////////////////////////////////////////////////

		printf("\n\nRIGHT MC, i = %d, j = %d, eff_GOPsz = %d, count+1 = %d\n",i,j,eff_GOPsz[i], count+1);
//		printf("eff_GOPsz = %d\n",eff_GOPsz[i]);

#ifdef  DEBUG_BLOCK_MODE_MV_INFO
		write_block_mode_motion_vector("encoderright", GOP_counter, count+1, 2, yfmv[count+1], info, i);
#endif
	    // subsample motion vector according to block size
	    layer_structure_mv_subsample(info, i, j+1,  yfmv[count+1], 1); 
        mvStat_setFrame(GOP_counter, i, j + 1);
        mv_encode(mvrate2, &maprate2, yfmv[count + 1], NULL, 1, 
                  info, dist, i, sub_mv2, sub_sign2, count+1, GOP_counter);
	    maprate1 = 0;
		for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
		{
			mvrate1[layer_num] = sub_sign1[layer_num] =  0;
			for (partition=0; partition<info.AGP_level[i]; partition++)
			   sub_mv1[layer_num][partition] = 0; 
		}
      }else{   //Added by Yuan Liu on 04.21.2016
		assert(scene_change[i][j] == YES && scene_change[i][j + 1] == YES);
		printf("\n\nISOLATED FRAMES, i = %d, j = %d, eff_GOPsz = %d, count = %d, count+1 = %d\n",i,j,eff_GOPsz[i], count,count+1);
		printf("\n\nISOLATED FRAMES, i = %d, j = %d, eff_GOPsz = %d, count+1 = %d, count = %d\n",i,j,eff_GOPsz[i], count+1,count);
		clear_frame_motion_field_simp(prev_frame_motion_field1, info);
		clear_frame_motion_field_simp(prev_frame_motion_field2, info);
	  }

	  printf("calc_res_sad = %f\n\n",calc_res_sad);

      FrsRate.map[count] = maprate1;
	  itmp = maprate1;
	  ctl_info += maprate1;
	  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
	  {
		  FrsRate.mv[layer_num][count] = mvrate1[layer_num]; 
		  itmp += mvrate1[layer_num];
		  mv_info += mvrate1[layer_num];
		  // for scalable motion vector coding: AGP 
		  for (partition=0; partition<info.AGP_level[i]; partition++)
		  {
				FrsRate.submv[layer_num][partition][count]=sub_mv1[layer_num][partition]; 
				itmp += sub_mv1[layer_num][partition]; 
				mv_info += sub_mv1[layer_num][partition]; 
		  }
		  if (info.AGP_level[i]) {
			  assert(0);
			  itmp += sub_sign1[layer_num]; 
			  ctl_info += sub_sign1[layer_num];
		  }
	  }

      if (itmp > 0) {
		// mapbits & mvbits
        fprintf(fpstat, "Num%.2d: mapbits = %4d ", count, FrsRate.map[count]);
		for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
		{
			fprintf(fpstat, "major_mvbits%d = %5d  ", layer_num, FrsRate.mv[layer_num][count]);
			for (partition=0; partition<info.AGP_level[i]; partition++)
				fprintf(fpstat, "sub_mvbits%d[%d] = %5d  ", 
						layer_num, partition,  FrsRate.submv[layer_num][partition][count]);
			if (info.AGP_level[i]) 	
				fprintf(fpstat, "sub_sign%d = %5d  ", layer_num, sub_sign1[layer_num]);
		}
		// sum 
        fprintf(fpstat, "sum = %5d (%d)\n", itmp, itmp >> 3);
		write_mv_length(info.AGP_level[i], info.layer_mv[i], fp_mv, FrsRate, count, sub_sign1, itmp, simul_enc); 
        // the overall sum
		if(do_res == NO)
			sum += itmp;
		else{
			assert(do_res == YES);
			if(j < level_res_frame[i])
				sum += itmp;
			else
				printf("Current mv set skipped!\n");
		}

        mv_count++;
      } else {
        fprintf(fpstat, "Fr%.2d: no motion information needed\n", count); 

//		if(simul_enc == NO)
			fprintf( fp_mv, "%d\n", 0);
      }
      count++;

      if (j < eff_GOPsz[i]) {

		  FrsRate.map[count] = maprate2;
		  itmp = maprate2;
		  ctl_info += maprate2;
		  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
		  {
			  FrsRate.mv[layer_num][count] = mvrate2[layer_num]; 
			  itmp += mvrate2[layer_num];
			  mv_info += mvrate2[layer_num];
			  // for scalable motion vector coding: AGP 
			  for (partition=0; partition<info.AGP_level[i]; partition++)
			  {
					FrsRate.submv[layer_num][partition][count]=sub_mv2[layer_num][partition]; 
					itmp += sub_mv2[layer_num][partition]; 
					mv_info += sub_mv2[layer_num][partition]; 
			  }
			  if (info.AGP_level[i]) {
				  assert(0);
				  itmp += sub_sign2[layer_num]; 
				  ctl_info += sub_sign2[layer_num]; 
			  }
		  }

		  if (itmp > 0) {
			  // mapbits & mvbits
			  fprintf(fpstat, "Num%.2d: mapbits = %4d ", count, FrsRate.map[count]);
			  for (layer_num=0; layer_num<LAYER_NUM; layer_num++)
			  {
				  fprintf(fpstat, "major_mvbits%d = %5d  ", layer_num, FrsRate.mv[layer_num][count]);
				  for (partition=0; partition<info.AGP_level[i]; partition++)
						fprintf(fpstat, "sub_mvbits%d[%d] = %5d  ", 
								layer_num, partition,  FrsRate.submv[layer_num][partition][count]);
				  if (info.AGP_level[i]) 	
						fprintf(fpstat, "sub_sign%d = %5d  ", layer_num, sub_sign2[layer_num]);
			  }
	  		  // sum 
			  fprintf(fpstat, "sum = %5d (%d)\n", itmp, itmp >> 3);

			  write_mv_length(info.AGP_level[i], info.layer_mv[i], fp_mv, FrsRate, count, sub_sign2, itmp, simul_enc); 
			  // the overall sum
			  if(do_res == NO)
				sum += itmp;
			  else{
				assert(do_res == YES);
				if(j < level_res_frame[i])
					sum += itmp;
				else
					printf("Current mv set skipped!\n");
			  }

			  mv_count++;
		  } else {
			  fprintf(fpstat, "Fr%.2d: no motion information needed\n", count);
//			  if(simul_enc == NO)
				fprintf(fp_mv, "%d\n", 0);
		  }
          count++;
      }
    }// j
    dist /= 2;
	//////////////////////
	avg_blk_size[i] = avg_blk_size[i]/blk_num;

	printf("before copy!\n");

	if(i <= info.tPyrLev - 2){
		copy_frame_motion_field2(prev_frame_motion_field1, buffer_frame_motion_field1[i], info);
		copy_frame_motion_field2(prev_frame_motion_field2, buffer_frame_motion_field2[i], info);
	}
	//////////////////////
  }
  fprintf( fpstat, "---> %.2d MV sets -- total sum = %7d (%d bytes)\n", mv_count, sum, sum >> 3 );
  printf ( "---> %.2d MV sets -- total sum = %7d (%d bytes)\n\n", mv_count, sum, sum >> 3 );
  fclose( fpstat );
  fclose( fp_mv );

  free(prev_frame_motion_field1);
  free(prev_frame_motion_field2);

  free(prev_frame_motion_field_left1);
  free(prev_frame_motion_field_left2);

  free(level_res_frame);

  if(GOP_counter == (get_GOP_num(info) - 1) ){
	  for(i = 0;i < G_LEVEL; i++){
		free(buffer_frame_motion_field1[i]);
		free(buffer_frame_motion_field2[i]);
		free(save_buffer_frame_motion_field1[i]);
		free(save_buffer_frame_motion_field2[i]);
	  }
  }

  printf("\n");
  for(i = 0; i < NUMBER_OF_BI_MODES; i ++){
	  printf("mean_mode_mse[%d] = %f, num[%d] = %d\n", i, mean_mode_mse[i]/mean_mode_num[i], i, mean_mode_num[i] );
  }
  printf("\n");
  
  return ( sum );
}




void
alloc_vector_child( vector_ptr *fmv1, vector_ptr *fmv2, int meandepth,
                    int cx, int cy, int xblk, int yblk, int hor, int ver,
                    int small, videoinfo info, int t_level )
{
  *fmv1 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv1" );

  if (fmv2 != NULL) {
    *fmv2 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv2" );
  }

  if( cx < hor && cy < ver ) {  /* Feb23 */
    if( xblk > small ) {
      child_map_decode( *fmv1, fmv2 ? (*fmv2) : NULL, meandepth, cx, cy,
                        xblk, yblk, hor, ver, small, info, t_level );
    } else {

      (*fmv1)->child = 0;
      if (fmv2 != NULL) {
        (*fmv2)->child = 0;
      }

      /* !!! NOTE THAT BLOCK MODES ARE ALSO READ in child_mv_decode !!! */
      decode_block_mode(*fmv1, fmv2 ? (*fmv2) : NULL, meandepth, info, t_level);
#ifdef BLOCKMODE_STATISTICS
      record_blockmode_statistics(xblk, yblk, (*fmv1)->bi_mode, 
                                  fmv2 != NULL, hor, ver);
#endif
    }
  } else {
    (*fmv1)->child = 0;
    if (fmv2 != NULL) {
      (*fmv2)->child = 0;
    }
  }
}


/****************************************************************************/
/*                              child_map_decode()                          */
/* decode quature structure, block modes and means(Y,U,V) of PREDICTEDS    */
/****************************************************************************/
void
child_map_decode( vector_ptr fmv1, vector_ptr fmv2, int meandepth,
                  int x, int y, int xblk, int yblk, int hor, int ver,
                  int small, videoinfo info, int t_level )
{
  int bit, cx, cy;

  input_bit( bit );
  fmv1->child = bit;

  if (fmv2 != NULL) {
    fmv2->child = bit;
  }

  if( fmv1->child ) {
    cx = x;
    cy = y;
    alloc_vector_child(&(fmv1->child0), fmv2 ? (&(fmv2->child0)) : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level);
    
    cx = x + xblk / 2;
    cy = y;
    alloc_vector_child(&(fmv1->child1), fmv2 ? (&(fmv2->child1)) : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level);

    cx = x;
    cy = y + yblk / 2;
    alloc_vector_child(&(fmv1->child2), fmv2 ? (&(fmv2->child2)) : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    alloc_vector_child(&(fmv1->child3), fmv2 ? (&(fmv2->child3)) : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, info, t_level);
  
  } else {                      /* Feb23 */
    /* !!! NOTE THAT BLOCK MODES ARE ALSO READ in alloc_vector_child !!! */
    decode_block_mode(fmv1, fmv2, meandepth, info, t_level);
#ifdef BLOCKMODE_STATISTICS
    record_blockmode_statistics(xblk, yblk, fmv1->bi_mode, 
                                fmv2 != NULL, hor, ver);
#endif
  }
}


/****************************************************************************/
/*                              child_merge_decode()                        */
/****************************************************************************/
void
child_merge_decode_further( vector_ptr fmv1, vector_ptr fmv2, int meandepth,
                  int x, int y, int xblk, int yblk, int hor, int ver,
                  int small, videoinfo info, int t_level, int blk_thresh, int count )
{
  int bit, cx, cy;

  if( fmv1->child && xblk>blk_thresh ) {
    cx = x;
    cy = y;
    child_merge_decode_further(fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);
    
    cx = x + xblk / 2;
    cy = y;
    child_merge_decode_further(fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_merge_decode_further(fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_merge_decode_further(fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);
  
  } else {   

    if( x >= hor || y >= ver )   return;

	if (!fmv1->child)     return; 

	if (!fmv1->mv_exist)  return; // no subsampled motion vector in this big block 

	// irregular motion area
	if ( fmv1->child0->child || fmv1->child1->child ||
		 fmv1->child2->child || fmv1->child3->child)
		 return;

	input_bit(bit);
	fmv1->merge_sign = bit; 

	FILE *fmerge; 
	char merge_file[80];
	sprintf(merge_file, "decode_merge_sign%d.txt", count); 
	fmerge = fopen(merge_file, "at"); 
	fprintf(fmerge, "%d\n", bit); 
	fclose(fmerge); 
  }
}


/****************************************************************************/
/*                              child_merge_decode()                          */
/****************************************************************************/
void
child_merge_decode( vector_ptr fmv1, vector_ptr fmv2, int meandepth,
                  int x, int y, int xblk, int yblk, int hor, int ver,
                  int small, videoinfo info, int t_level, int blk_thresh, int count )
{
  int bit, cx, cy;

  if( fmv1->child && xblk>blk_thresh ) {
    cx = x;
    cy = y;
    child_merge_decode(fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);
    
    cx = x + xblk / 2;
    cy = y;
    child_merge_decode(fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_merge_decode(fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_merge_decode(fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                       meandepth, cx, cy, xblk / 2, yblk / 2, hor, ver, small, 
					   info, t_level, blk_thresh, count);
  
  } else {   

    if( x >= hor || y >= ver )   return;

	if (!fmv1->child)     return; 

	if (!fmv1->mv_exist)  return; // no subsampled motion vector in this big block 

	assert(0);

	input_bit(bit);
	fmv1->merge_sign = bit; 

	FILE *fmerge; 
	char merge_file[80];
	sprintf(merge_file, "decode_merge_sign%d.txt", count); 
	fmerge = fopen(merge_file, "at"); 
	fprintf(fmerge, "%d\n", bit); 
	fclose(fmerge); 

	if (!fmv1->merge_sign )
		child_merge_decode_further(  fmv1,  fmv2,  meandepth,
                   x,  y,  xblk,  yblk,  hor,  ver, small,  info,  t_level,  blk_thresh/2,  count );
  }
}



char get_splitted_mvbits(char whichsub, int layer_index)
{
	int value = 0;

	if ( --(splitted_bit_num_array[layer_index][whichsub]) <0 )
	{
		splitted_mv_byte_array[layer_index][whichsub] = 
			store_splitted_bytes_array[layer_index][whichsub][splitted_mv_byte_num_array[layer_index][whichsub]++];
		splitted_bit_num_array[layer_index][whichsub] = 7;
	}
	value = ((splitted_mv_byte_array[layer_index][whichsub] >> splitted_bit_num_array[layer_index][whichsub]) & 0x01);     

	return value;

}


char get_splitted_signbits(int layer_index)
{
	int value = 0;

	if ( --splitted_sign_bit_num[layer_index] <0 )
	{
		splitted_sign_byte[layer_index] = store_splitted_sign[layer_index][splitted_sign_byte_num[layer_index]++];
		splitted_sign_bit_num[layer_index] = 7;
	}
	value = ((splitted_sign_byte[layer_index] >> splitted_sign_bit_num[layer_index]) & 0x01);     
	return value;
}

/****************************************************************************/
/*                              child_mv_decode_enhance_sub()               */
/****************************************************************************/
void
child_mv_decode_enhance_sub( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int *first_available,  int count, int merge_sign )
{
  int dmvx, dmvy, cx, cy;
  int ctx_x, ctx_y;
  float mvpred_x, mvpred_y;
  char value; 
  int xblk2, yblk2; 
  float  major_mvx, major_mvy, sub_mvx, sub_mvy;
  float  major_predx, major_predy; 
  int    AGP_scale, subpel_scale, sub_symx, sub_symy, sub_bit, sub_signx, sub_signy; 
  char   sub_sym_predx, sub_sym_predy; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk > blk_thresh ) {     
    cx = x;
    cy = y;
    child_mv_decode_enhance_sub(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, 
					first_available, count, merge_sign);

    cx = x + xblk / 2;
    cy = y;
    child_mv_decode_enhance_sub(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, 
					first_available, count, merge_sign);

    cx = x;
    cy = y + yblk / 2;
    child_mv_decode_enhance_sub(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, 
					first_available, count, merge_sign);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_decode_enhance_sub(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, 
					first_available, count, merge_sign);
  } else {
    if( x >= hor || y >= ver )      return;

	// no motion vector for this block on this side 
	if ((fmv1->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
	if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 

	assert(0);

	if (!merge_sign)
	{
		// re-do the spatial prediction when it's first motion vector or not 
        get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, 0);
		/////////////////////////////////
		if( !((fmv1->bi_mode == PARALLEL) && !decode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
			for(i=0;i<=3;i++){
				if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
					num++;
			}

			if(num == 0){
				fmv1->med_idx = -1;
				mvpred_x = 0.0;
				mvpred_y = 0.0;
			}
			else{
				fmv1->med_idx = getbits(2);
				assert(fmv1->med_idx >= 0);
				mvpred_x = mvpredstr_x[fmv1->med_idx];
				mvpred_y = mvpredstr_y[fmv1->med_idx];
			}
		}
		else if(fmv1->bi_mode == BLOCK_MERGING ){
		  num = 0;
		  for(i=0;i<=3;i++){
			  if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		  }
		  switch(num){
			  case 0:
				  mvpred_x = 0.0;
				  mvpred_y = 0.0;
				  break;

			  case 1:
			  case 2:
				  if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[0];
					mvpred_y = mvpredstr_y[0];
				  }
				  else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[1];
					mvpred_y = mvpredstr_y[1];
				  }
				  else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[2];
					mvpred_y = mvpredstr_y[2];
				  }
				  else{
					assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
					mvpred_x = mvpredstr_x[3];
					mvpred_y = mvpredstr_y[3];
				  }
				  break;

			  case 3:
				  if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				  }
				  else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				  }
				  break;
		  }
		}
		/////////////////////////////////

		if ( !info.AGP_level[t_level] )
		{
			fmv1->dmvx =  fmv1->mvx - mvpred_x;
			fmv1->dmvy =  fmv1->mvy - mvpred_y;
		}else
		{
			AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
			subpel_scale = 1<<subpel; 
			// major symbol
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			// sub-symbol (already converted to integer)
			sub_symx = (char)(fabs( (fmv1->mvx - major_mvx) * subpel_scale ));  
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			sub_symy = (char)(fabs( (fmv1->mvy - major_mvy) * subpel_scale ));  
			major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
			sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
			major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
			sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale)); 
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
		fmv1->is_predictor = YES;
   	    update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
				fmv1->dmvy, fmv1->mvx, fmv1->mvy); 
		return; 
	}

	if ( *first_available ) // the first motion vector is skipped since it's already coded in base layer
	{
		// re-do the spatial prediction when it's first motion vector or not 
        get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, 0);
		/////////////////////////////////
		if( !((fmv1->bi_mode == PARALLEL) && !decode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
			for(i=0;i<=3;i++){
				if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
					num++;
			}

			if(num == 0){
				fmv1->med_idx = -1;
				mvpred_x = 0.0;
				mvpred_y = 0.0;
			}
			else{
				fmv1->med_idx = getbits(2);
				assert(fmv1->med_idx >= 0);
				mvpred_x = mvpredstr_x[fmv1->med_idx];
				mvpred_y = mvpredstr_y[fmv1->med_idx];
			}
		}
		else if(fmv1->bi_mode == BLOCK_MERGING ){
		  num = 0;
		  for(i=0;i<=3;i++){
			  if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		  }
		  switch(num){
			  case 0:
				  mvpred_x = 0.0;
				  mvpred_y = 0.0;
				  break;

			  case 1:
			  case 2:
				  if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[0];
					mvpred_y = mvpredstr_y[0];
				  }
				  else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[1];
					mvpred_y = mvpredstr_y[1];
				  }
				  else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
					mvpred_x = mvpredstr_x[2];
					mvpred_y = mvpredstr_y[2];
				  }
				  else{
					assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
					mvpred_x = mvpredstr_x[3];
					mvpred_y = mvpredstr_y[3];
				  }
				  break;

			  case 3:
				  if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				  }
				  else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				  }
				  else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				  }
				  break;
		  }
		}
		/////////////////////////////////
		if ( !info.AGP_level[t_level] )
		{
			fmv1->dmvx =  fmv1->mvx - mvpred_x;
			fmv1->dmvy =  fmv1->mvy - mvpred_y;
		}else
		{
			AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
			subpel_scale = 1<<subpel; 
			// major symbol
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			// sub-symbol (already converted to integer)
			sub_symx = (char)(fabs( (fmv1->mvx - major_mvx) * subpel_scale ));  
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			sub_symy = (char)(fabs( (fmv1->mvy - major_mvy) * subpel_scale ));  
			major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
			sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
			major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
			sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale)); 
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
		fmv1->is_predictor = YES;
   	    update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
				fmv1->dmvy, fmv1->mvx, fmv1->mvy); 
		*first_available = 0; 
		return;  
	}

	if (!bidir_exist) // the motion vectors on RIGHT side has been discarded
	{
		if ((fmv1->bi_mode == PARALLEL)&&(!decode_parallelmv)) {
			fmv1->mvx = -fmv2->mvx;
			fmv1->mvy = -fmv2->mvy;
		} 
		if (fmv1->bi_mode==RIGHT_CONNECTED)
			assert(0); // RIGHT_CONNECTED mode is prevented 
		return; 
	}

	if ( (fmv1->bi_mode != BLOCK_MERGING)&& !((fmv1->bi_mode == PARALLEL)&&(!decode_parallelmv)) ) {
	    ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, x, y, info, t_level, blk_thresh);

		dmvx = ec_decode_word(ctx_x);   // decode the mvx 
		ec_update_model(dmvx, ctx_x);
		if (EC_TYPE == AR_NARY) {
			if (dmvx > num_symbol / 2 - (int) ((1 << subpel) * *pmvx))
			  dmvx -= num_symbol;
			else if (dmvx < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvx))
			  dmvx += num_symbol;
		}
      
		dmvy = ec_decode_word(ctx_y);   // decode the mvy 
		ec_update_model(dmvy, ctx_y);
		if (EC_TYPE == AR_NARY) {
			if (dmvy > num_symbol / 2 - (int) ((1 << subpel) * *pmvy))
			  dmvy -= num_symbol;
			else if (dmvy < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvy))
			  dmvy += num_symbol;
		  }
	} else { // BLOCK_MERGING
	  dmvx = 0;
	  dmvy = 0;
	}

    get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, 0);
	/////////////////////////////////
	if( !((fmv1->bi_mode == PARALLEL) && !decode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}

		if(num == 0){
			fmv1->med_idx = -1;
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			fmv1->med_idx = getbits(2);
			assert(fmv1->med_idx >= 0);
			mvpred_x = mvpredstr_x[fmv1->med_idx];
			mvpred_y = mvpredstr_y[fmv1->med_idx];
		}
	}
	else if(fmv1->bi_mode == BLOCK_MERGING ){
		num = 0;
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}
		switch(num){
		  case 0:
			mvpred_x = 0.0;
			mvpred_y = 0.0;
			break;

		  case 1:
		  case 2:
			if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[0];
				mvpred_y = mvpredstr_y[0];
			}
			else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[1];
				mvpred_y = mvpredstr_y[1];
			}
			else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[2];
				mvpred_y = mvpredstr_y[2];
			}
			else{
				assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
				mvpred_x = mvpredstr_x[3];
				mvpred_y = mvpredstr_y[3];
			}
			break;

			case 3:
				if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				}
				else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				}
				break;
		}
	}
	/////////////////////////////////
	if ( !info.AGP_level[t_level] )
	{
		if ((fmv1->bi_mode == PARALLEL)&&(!decode_parallelmv)) {
			fmv1->mvx = -fmv2->mvx;
			fmv1->mvy = -fmv2->mvy;
		} else {
			fmv1->mvx = (float) dmvx / (1 << subpel) + mvpred_x;
			fmv1->mvy = (float) dmvy / (1 << subpel) + mvpred_y;
		}
	    fmv1->dmvx = fmv1->mvx - mvpred_x;
		fmv1->dmvy = fmv1->mvy - mvpred_y;
	}else
	{
		AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
		subpel_scale = 1<<subpel; 
		major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
		major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;

		if ((fmv1->bi_mode == PARALLEL)&&(!decode_parallelmv)) {
			fmv1->mvx = -fmv2->mvx;
			fmv1->mvy = -fmv2->mvy;
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
		} else if ( fmv1->bi_mode == BLOCK_MERGING )
		{
			fmv1->mvx = mvpred_x;
			fmv1->mvy = mvpred_y;
	        assert(dmvx ==0 && dmvy == 0 );
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
		}else		
		{
			major_mvx = (float) dmvx / AGP_scale  + major_predx;
			major_mvy = (float) dmvy / AGP_scale  + major_predy;
			// code the sub-symbols as binary sequence 
			sub_symx = sub_symy = 0; 
			for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
			{
				sub_symx <<= 1; 
				if ( sub_bit < info.AGP_exist[t_level] )
				{
					value = get_splitted_mvbits( sub_bit, 1 );
					sub_symx += value;
				}
				sub_symy <<= 1; 
				if ( sub_bit < info.AGP_exist[t_level] ) 
				{
					value = get_splitted_mvbits( sub_bit, 1 );
					sub_symy += value; 
				}
			}
			sub_mvx = (float)sub_symx/ (float)subpel_scale;  
			sub_mvy = (float)sub_symy/ (float)subpel_scale;  

			if (info.AGP_exist[t_level] && major_mvx==0) 
			{
				sub_signx = get_splitted_signbits(1);  // 0: positive, 1: negative
				if (sub_signx)
					sub_mvx *= -1;
			}else if (major_mvx<0)
				sub_mvx *= -1;
			if (info.AGP_exist[t_level] && major_mvy==0) 
			{
				sub_signy = get_splitted_signbits(1);  // 0: positive, 1: negative
				if (sub_signy )
					sub_mvy *= -1;
			}else if (major_mvy<0)
				sub_mvy *= -1;
			fmv1->mvx = major_mvx + sub_mvx;
			fmv1->mvy = major_mvy + sub_mvy; 
		}
		// only major symbols are median predicted, sub-symbols are coded by binary sequence
		fmv1->dmvx = major_mvx - major_predx;
		fmv1->dmvy = major_mvy - major_predy;

#ifdef  DEBUG_SCALABLE_MV
		fpAGP_debug = fopen(decoder_AGPdebug, "at");
		fprintf(fpAGP_debug, "x=%03d  y=%03d  blk=%2d ", x, y, xblk ); 
		fprintf(fpAGP_debug, "major_mvx=%.1f  major_mvy=%.1f  major_predx=%.1f  major_predy=%.1f",
			    major_mvx,  major_mvy,  major_predx,  major_predy) ;
		fprintf(fpAGP_debug, " dmvx=%.1f  dmvy=%.1f  mvx=%.2f  mvy=%.2f\n", 
			                 fmv1->dmvx,  fmv1->dmvy, fmv1->mvx, fmv1->mvy) ;
		fclose(fpAGP_debug);
#endif 

	}

	// double check 
    assert((fmv1->bi_mode != BLOCK_MERGING)||(fmv1->dmvx==0 && fmv1->dmvy==0)); // just to be shure. mwi
   
    mvStat_setPos(x, y);
    mvStat_setDMV((float)dmvx, (float)dmvy);
    mvStat_writeDMVCTX();
    fmv1->is_predictor = YES;
    update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
							  fmv1->dmvy, fmv1->mvx, fmv1->mvy); 


#ifdef  DEBUG_LAYER_STRUCTURE  // write the enhancement motion vectors 
    char  enhance_file[80]; 
    FILE *fpenhance; 
    // make base_file and enhance_file empty
    sprintf(enhance_file, "decode_enhance%d.txt", count); 
    fpenhance=fopen(enhance_file, "at");
    fprintf(fpenhance, "x=%03d y=%03d blk=%d\t mvx=%.2f\t mvy=%.2f\t predx=%.2f predy=%.2f \n", 
			x, y, xblk, fmv1->mvx, fmv1->mvy, mvpred_x, mvpred_y); 
    fclose(fpenhance);
#endif 
    
    xblk2 = ( x + xblk <= hor ) ? xblk : hor - x;
    yblk2 = ( y + yblk <= ver ) ? yblk : ver - y;
    
    assert( (x + xblk2 - fmv1->mvx <= hor ) &&
            (y + yblk2 - fmv1->mvy <= ver ) &&
            (x - fmv1->mvx >= 0 ) &&
            (y - fmv1->mvy >= 0 ) );

    if (!decode_parallelmv) {
      assert( (!(fmv1->bi_mode == PARALLEL)) || 
              ((x + xblk2 - fmv2->mvx <= hor ) &&
               (y + yblk2 - fmv2->mvy <= ver ) &&
               (x - fmv2->mvx >= 0 ) &&
               (y - fmv2->mvy >= 0 )) );
    }

  }
}

void
child_mv_decode_enhance_further( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count )
{
  int   cx, cy;
  float mvpred_x, mvpred_y;
  int    first_available; 
  float major_mvx, major_mvy, major_predx, major_predy;
  int   sub_symx, sub_symy;
  char  sub_sym_predx, sub_sym_predy;
  int   AGP_scale, subpel_scale; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk > blk_thresh ) {     
    cx = x;
    cy = y;
    child_mv_decode_enhance_further(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_decode_enhance_further(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_decode_enhance_further(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_decode_enhance_further(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);
  } else {

    if( x >= hor || y >= ver )      return;

	if (!(fmv1->child))  // if there are no children for this block with size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   return;

#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
		// re-do the spatial prediction 
		assert(0);
        get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, 0);
	/////////////////////////////////
	if( !((fmv1->bi_mode == PARALLEL) && !decode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}

		if(num == 0){
			fmv1->med_idx = -1;
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			fmv1->med_idx = getbits(2);
			assert(fmv1->med_idx >= 0);
			mvpred_x = mvpredstr_x[fmv1->med_idx];
			mvpred_y = mvpredstr_y[fmv1->med_idx];
		}
	}
	else if(fmv1->bi_mode == BLOCK_MERGING ){
		num = 0;
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}
		switch(num){
		  case 0:
			mvpred_x = 0.0;
			mvpred_y = 0.0;
			break;

		  case 1:
		  case 2:
			if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[0];
				mvpred_y = mvpredstr_y[0];
			}
			else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[1];
				mvpred_y = mvpredstr_y[1];
			}
			else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[2];
				mvpred_y = mvpredstr_y[2];
			}
			else{
				assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
				mvpred_x = mvpredstr_x[3];
				mvpred_y = mvpredstr_y[3];
			}
			break;

			case 3:
				if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				}
				else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				}
				break;
		}
	}
	/////////////////////////////////

		if ( !info.AGP_level[t_level] )
		{
			fmv1->dmvx =  fmv1->mvx - mvpred_x;
			fmv1->dmvy =  fmv1->mvy - mvpred_y;
		}else
		{
			AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
			subpel_scale = 1<<subpel; 
			// major symbol
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			// sub-symbol (already converted to integer)
			sub_symx = (char)(fabs( (fmv1->mvx - major_mvx) * subpel_scale ));  
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			sub_symy = (char)(fabs( (fmv1->mvy - major_mvy) * subpel_scale ));  
			major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
			sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
			major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
			sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale));  
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
		fmv1->is_predictor = YES;
   	    update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
				fmv1->dmvy, fmv1->mvx, fmv1->mvy); 
		return;  
	} else
	{
		if ( !fmv1->mv_exist )  return; 

		// irregular motion area
		if ( fmv1->child0->child || fmv1->child1->child ||
			 fmv1->child2->child || fmv1->child3->child)
			 fmv1->merge_sign = 0; 

		first_available = 1; 
		blk_thresh = 0;
		child_mv_decode_enhance_sub( fmv1_array, fmv2_array, fmv1,  fmv2, pmvx, pmvy, 
									 num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
									 info,  decode_parallelmv,  t_level,  bidir_exist, 
									 blk_thresh, &first_available,  count, fmv1->merge_sign );
	}

  }

}

void
child_mv_decode_enhance( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count )
{
  int   cx, cy;
  float mvpred_x, mvpred_y;
  int    first_available; 
  float major_mvx, major_mvy, major_predx, major_predy;
  int   sub_symx, sub_symy, sub_sym_predx, sub_sym_predy; 
  int   AGP_scale, subpel_scale; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk > blk_thresh ) {     
    cx = x;
    cy = y;
    child_mv_decode_enhance(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_decode_enhance(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_decode_enhance(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_decode_enhance(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);
  } else {

    if( x >= hor || y >= ver )      return;

	if (!(fmv1->child))  // if there are no children for this block with size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   return;

#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
		assert(0);
		// re-do the spatial prediction 
        get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, 0);
	  /////////////////////////////////
	if( !((fmv1->bi_mode == PARALLEL) && !decode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}

		if(num == 0){
			fmv1->med_idx = -1;
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			fmv1->med_idx = getbits(2);
			assert(fmv1->med_idx >= 0);
			mvpred_x = mvpredstr_x[fmv1->med_idx];
			mvpred_y = mvpredstr_y[fmv1->med_idx];
		}
	}
	else if(fmv1->bi_mode == BLOCK_MERGING ){
		num = 0;
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}
		switch(num){
		  case 0:
			mvpred_x = 0.0;
			mvpred_y = 0.0;
			break;

		  case 1:
		  case 2:
			if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[0];
				mvpred_y = mvpredstr_y[0];
			}
			else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[1];
				mvpred_y = mvpredstr_y[1];
			}
			else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[2];
				mvpred_y = mvpredstr_y[2];
			}
			else{
				assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
				mvpred_x = mvpredstr_x[3];
				mvpred_y = mvpredstr_y[3];
			}
			break;

			case 3:
				if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				}
				else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				}
				break;
		}
	}
	/////////////////////////////////

		if ( !info.AGP_level[t_level] )
		{
			fmv1->dmvx =  fmv1->mvx - mvpred_x;
			fmv1->dmvy =  fmv1->mvy - mvpred_y;
		}else
		{
			AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
			subpel_scale = 1<<subpel; 
			// major symbol
			major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
			// sub-symbol (already converted to integer)
			sub_symx = (char)(fabs( (fmv1->mvx - major_mvx) * subpel_scale ));  
			major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			sub_symy = (char)(fabs( (fmv1->mvy - major_mvy) * subpel_scale ));  
			major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
			sub_sym_predx = (char)(fabs( (mvpred_x - major_predx)*subpel_scale )); 
			major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
			sub_sym_predy = (char)(fabs( (mvpred_y - major_predy)*subpel_scale));  
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
		fmv1->is_predictor = YES;
   	    update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
				fmv1->dmvy, fmv1->mvx, fmv1->mvy); 
		return;  
	} else
	{
		if ( !fmv1->mv_exist ) return; 

		if ( fmv1->merge_sign )
		{
			first_available = 1; 
			blk_thresh = 0;
			child_mv_decode_enhance_sub( fmv1_array, fmv2_array, fmv1,  fmv2, pmvx, pmvy, 
										 num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
										 info,  decode_parallelmv,  t_level,  bidir_exist, 
										 blk_thresh, &first_available,  count, fmv1->merge_sign );
		}else
		{
			child_mv_decode_enhance_further(  fmv1_array,  fmv2_array,  fmv1,  fmv2,
                 pmvx, pmvy, num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
                  info,  decode_parallelmv,  t_level,  bidir_exist,  blk_thresh/2,  count );
		}
	}

  }

}

/****************************************************************************/
/*                              child_mv_decode_further()                           */
/****************************************************************************/
void
child_mv_decode_further( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count )
{
  int dmvx, dmvy, cx, cy;
  int ctx_x, ctx_y;
#ifdef MEDIAN_PREDICTION
  float mvpred_x, mvpred_y;
#endif
  char value; 
  float  major_mvx, major_mvy, sub_mvx, sub_mvy;
  float  major_predx, major_predy; 
  int    AGP_scale, subpel_scale, sub_symx, sub_symy, sub_bit, sub_signx, sub_signy; 

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk > blk_thresh ) {     
    cx = x;
    cy = y;
    child_mv_decode_further(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_decode_further(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_decode_further(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_decode_further(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);
  } else {

    if( x >= hor || y >= ver )      return;

	if ( !fmv1->child )  // there are no children for this block and size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
	} else
	{
		assert(0);
		if ( !fmv1->mv_exist )   return; // the whole (blk_thresh x blk_thresh) has no motion vector
		
		// irregular motion area
		if ( fmv1->child0->child || fmv1->child1->child ||
			 fmv1->child2->child || fmv1->child3->child)
		{
			blk_thresh = 0; 
			child_mv_decode_further( fmv1_array, fmv2_array, fmv1, fmv2, pmvx, pmvy,  
				             num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
                             info,  decode_parallelmv,  t_level,  bidir_exist, blk_thresh/2,  count );
			return;
		}
			 
		if ( !fmv1->merge_sign )
		{
			blk_thresh = 0; 
			child_mv_decode_further( fmv1_array, fmv2_array, fmv1, fmv2, pmvx, pmvy,  
				             num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
                             info,  decode_parallelmv,  t_level,  bidir_exist, blk_thresh/2,  count );
			return; 
		}
	}

	// the motion vectors on RIGHT side has been discarded for alternative reconstruction
	if (!bidir_exist) 
	{
		if ((fmv1->bi_mode == PARALLEL)&&(!decode_parallelmv)) {
			fmv1->mvx = -fmv2->mvx;
			fmv1->mvy = -fmv2->mvy;
		} 
		if (fmv1->bi_mode==RIGHT_CONNECTED)
			assert(0); // RIGHT_CONNECTED mode is prevented 
		return; 
	}

    ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, x, y, info, t_level, blk_thresh);
	dmvx = ec_decode_word(ctx_x);   // decode the mvx 
    ec_update_model(dmvx, ctx_x);
    dmvy = ec_decode_word(ctx_y);   // decode the mvy 
	ec_update_model(dmvy, ctx_y);

    get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, blk_thresh);
	/////////////////////////////////
	if( !((fmv1->bi_mode == PARALLEL) && !decode_parallelmv) && fmv1->bi_mode != BLOCK_MERGING ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}

		if(num == 0){
			fmv1->med_idx = -1;
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			fmv1->med_idx = getbits(2);
			assert(fmv1->med_idx >= 0);
			mvpred_x = mvpredstr_x[fmv1->med_idx];
			mvpred_y = mvpredstr_y[fmv1->med_idx];
		}
	}
	else if(fmv1->bi_mode == BLOCK_MERGING ){
		num = 0;
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL)
				num++;
		}
		switch(num){
		  case 0:
			mvpred_x = 0.0;
			mvpred_y = 0.0;
			break;

		  case 1:
		  case 2:
			if(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[0];
				mvpred_y = mvpredstr_y[0];
			}
			else if(mvpredstr_x[1] != (float)HUGE_VAL && mvpredstr_y[1] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[1];
				mvpred_y = mvpredstr_y[1];
			}
			else if(mvpredstr_x[2] != (float)HUGE_VAL && mvpredstr_y[2] != (float)HUGE_VAL){
				mvpred_x = mvpredstr_x[2];
				mvpred_y = mvpredstr_y[2];
			}
			else{
				assert(mvpredstr_x[3] != (float)HUGE_VAL && mvpredstr_y[3] != (float)HUGE_VAL);
				mvpred_x = mvpredstr_x[3];
				mvpred_y = mvpredstr_y[3];
			}
			break;

			case 3:
				if(mvpredstr_x[0] == (float)HUGE_VAL && mvpredstr_y[0] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[1],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[1],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[1] == (float)HUGE_VAL && mvpredstr_y[1] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[2],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[2],mvpredstr_y[3]);
				}
				else if(mvpredstr_x[2] == (float)HUGE_VAL && mvpredstr_y[2] == (float)HUGE_VAL){
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[3]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[3]);
				}
				else{
					assert(mvpredstr_x[3] == (float)HUGE_VAL && mvpredstr_y[3] == (float)HUGE_VAL);
					mvpred_x = medi(mvpredstr_x[0],mvpredstr_x[1],mvpredstr_x[2]);
					mvpred_y = medi(mvpredstr_y[0],mvpredstr_y[1],mvpredstr_y[2]);
				}
				break;
		}
	}
	/////////////////////////////////

	if (!info.AGP_level[t_level]) // no AGP, layer structure 
	{
		// subsampled motion vector 
		fmv1->sample_mvx = fmv1->mvx = (float) dmvx / (1 << subpel) + mvpred_x;
		fmv1->sample_mvy = fmv1->mvy = (float) dmvy / (1 << subpel) + mvpred_y;
	    fmv1->dmvx = fmv1->mvx - mvpred_x;  // update the prediction error
		fmv1->dmvy = fmv1->mvy - mvpred_y;
	}else  // AGP, layer structure 
	{
		AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
		subpel_scale = 1<<subpel; 
		major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
		major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;
		major_mvx = (float) dmvx / AGP_scale  + major_predx;
		major_mvy = (float) dmvy / AGP_scale  + major_predy;
		// decode the sub-symbols as binary sequence 
		sub_symx = sub_symy = 0; 
		for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
		{
			sub_symx <<= 1; 
			if ( sub_bit < info.AGP_exist[t_level] )
			{
				value = get_splitted_mvbits( sub_bit, 0 );
				sub_symx += value;
			}
			sub_symy <<= 1; 
			if ( sub_bit < info.AGP_exist[t_level] ) 
			{
				value = get_splitted_mvbits( sub_bit, 0 );
				sub_symy += value; 
			}
		}
		sub_mvx = (float)sub_symx/ (float)subpel_scale;  
		sub_mvy = (float)sub_symy/ (float)subpel_scale;  

		if (info.AGP_exist[t_level] && major_mvx==0) 
		{
			sub_signx = get_splitted_signbits(0);  // 0: positive, 1: negative
			if (sub_signx) sub_mvx *= -1;
		}else if (major_mvx<0)
			sub_mvx *= -1;
		if (info.AGP_exist[t_level] && major_mvy==0) 
		{
			sub_signy = get_splitted_signbits(0);  // 0: positive, 1: negative
			if (sub_signy ) sub_mvy *= -1;
		}else if (major_mvy<0)
			sub_mvy *= -1;
		fmv1->sample_mvx = fmv1->mvx = major_mvx + sub_mvx;
		fmv1->sample_mvy = fmv1->mvy = major_mvy + sub_mvy; 
		// only major symbols are median predicted, sub-symbols are coded by binary sequence
		fmv1->dmvx = major_mvx - major_predx;
		fmv1->dmvy = major_mvy - major_predy;

#ifdef  DEBUG_SCALABLE_MV
		fpAGP_debug = fopen(decoder_AGPdebug, "at");
		fprintf(fpAGP_debug, "x=%03d  y=%03d  blk=%2d ", x, y, xblk ); 
		fprintf(fpAGP_debug, "major_mvx=%.1f  major_mvy=%.1f  major_predx=%.1f  major_predy=%.1f",
			    major_mvx,  major_mvy,  major_predx,  major_predy) ;
		fprintf(fpAGP_debug, " dmvx=%.1f  dmvy=%.1f  mvx=%.2f  mvy=%.2f\n", 
			                 fmv1->dmvx,  fmv1->dmvy, fmv1->mvx, fmv1->mvy) ;
		fclose(fpAGP_debug);
#endif 
	}

    fmv1->is_predictor = YES;
	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
				fmv1->dmvy, fmv1->mvx, fmv1->mvy); 

    mvStat_setPos(x, y);
    mvStat_setDMV((float)dmvx, (float)dmvy);
    mvStat_writeDMVCTX();

#ifdef  DEBUG_LAYER_STRUCTURE  // for debug the layer structure coding 
	char  base_file[80]; 
    FILE *fpbase; 
    // make base_file and enhance_file empty
    sprintf(base_file, "decode_base%d.txt", count); 
    fpbase=fopen(base_file, "at");
    fprintf(fpbase, "x=%03d y=%03d blk=%d\t mvx=%.2f\t mvy=%.2f\t predx=%.2f predy=%.2f \n", 
			x, y, xblk, fmv1->mvx, fmv1->mvy, mvpred_x, mvpred_y); 
	fclose(fpbase);
#endif 
    
  }

}


/****************************************************************************/
/*                              child_mv_decode()                           */
/****************************************************************************/
void
child_mv_decode( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count )
{
  int dmvx, dmvy, cx, cy, xblk2, yblk2;
  int aff1_dmvx, aff1_dmvy, aff2_dmvx, aff2_dmvy, aff3_dmvx, aff3_dmvy;//Added by Yuan Liu on 03.20.2016
  int ctx_x, ctx_y;
#ifdef MEDIAN_PREDICTION
  float mvpred_x, mvpred_y;
#endif
  char value; 
  float  major_mvx, major_mvy, sub_mvx, sub_mvy;
  float  major_predx, major_predy; 
  int    AGP_scale, subpel_scale, sub_symx, sub_symy, sub_bit, sub_signx, sub_signy; 
  int	enc_trans = 0;
  int   getnum;

  ////////////  Added by Yuan Liu on 01.23.2016  //////////////
  float mvpredstr_x[4], mvpredstr_y[4];
  int i, num = 0;

  for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
  }
  /////////////////////////////////////////////////////////////

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  if( fmv1->child && xblk > blk_thresh ) {     
    cx = x;
    cy = y;
    child_mv_decode(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_mv_decode(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_mv_decode(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_mv_decode(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);
  } else {
    if( x >= hor || y >= ver )      return;

	if ( !fmv1->child )  // there are no children for this block and size>=blk_thresh
	{
		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
	} else  //This will never happen!
	{
		if ( !fmv1->mv_exist )   return; // the whole (blk_thresh x blk_thresh) has no motion vector

		if ( !fmv1->merge_sign )
		{
			child_mv_decode_further( fmv1_array, fmv2_array, fmv1, fmv2, pmvx, pmvy,  
				             num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
                             info,  decode_parallelmv,  t_level,  bidir_exist, blk_thresh/2,  count );
			return; 
		}
	}

	// the motion vectors on RIGHT side has been discarded for alternative reconstruction
	if (!bidir_exist) 
	{
		if ( (fmv1->bi_mode == PARALLEL)&&(!decode_parallelmv) ) {
			assert(0);
			fmv1->mvx = -fmv2->mvx;
			fmv1->mvy = -fmv2->mvy;
			printf("!bidir_exist parallel detected\n");
		} 
		if ( (fmv1->bi_mode == RIGHT_CONNECTED_AFF) || fmv1->bi_mode==RIGHT_CONNECTED)
			assert(0); // RIGHT_CONNECTED mode is prevented 
		return; 
	}

	assert(fmv1->bi_mode >= 0);	//Added ib 04.02.2016
	if(fmv2!=NULL)
		assert(fmv1->bi_mode == fmv2->bi_mode);
	
	xblk2 = ( x + xblk <= hor ) ? xblk : hor - x;
	yblk2 = ( y + yblk <= ver ) ? yblk : ver - y;

	if ( !info.layer_mv[t_level] )  // no layer structure, AGP or no AGP  
	{
		if ( ( fmv1->bi_mode<=8 && fmv1->bi_mode != BLOCK_MERGING && !(fmv1->bi_mode == PARALLEL && !decode_parallelmv) ) ||
	    ( fmv1->bi_mode>=9 && fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == INDIR) ){
		  ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, 
							x, y, info, t_level, blk_thresh);

		  dmvx = ec_decode_word(ctx_x);   // decode the mvx 
		  ec_update_model(dmvx, ctx_x);
		  if (EC_TYPE == AR_NARY) {
			if (dmvx > num_symbol / 2 - (int) ((1 << subpel) * *pmvx))
			  dmvx -= num_symbol;
			else if (dmvx < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvx))
			  dmvx += num_symbol;
		  }
      
		  dmvy = ec_decode_word(ctx_y);   // decode the mvy 
		  ec_update_model(dmvy, ctx_y);
		  if (EC_TYPE == AR_NARY) {
			if (dmvy > num_symbol / 2 - (int) ((1 << subpel) * *pmvy))
			  dmvy -= num_symbol;
			else if (dmvy < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvy))
			  dmvy += num_symbol;
		  }
		}else if( fmv1->bi_mode <=8 || ( fmv1->bi_mode>=9 && fmv1->direct_idx == INDIRECT 
		   && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == DIR) ){
		  dmvx = 0;
		  dmvy = 0;
		//////////////////////////  Added by Yuan Liu  /////////////////////////////
		}
		
		if( fmv1->bi_mode >=9 && fmv1->bi_mode<=11 ){  //gotta decode trans MVs if they are applied
		  
		    if( !(fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
				dmvx = 0;
				dmvy = 0;
			}

		  if(fmv1->direct_idx == INDIRECT){
//affine V1
			if(fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ||
			   (fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0) ){
//				printf("unusual encoding processes detected!\n");
//				assert(0);
				ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, 
				x, y, info, t_level, blk_thresh);

				aff1_dmvx = ec_decode_word(ctx_x);   // decode the mvx 
				ec_update_model(aff1_dmvx, ctx_x);
				if (EC_TYPE == AR_NARY) {
					if (aff1_dmvx > num_symbol / 2 - (int) ((1 << subpel) * *pmvx))
						aff1_dmvx -= num_symbol;
					else if (aff1_dmvx < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvx))
						aff1_dmvx += num_symbol;
				}
      
				aff1_dmvy = ec_decode_word(ctx_y);   // decode the mvy 
				ec_update_model(aff1_dmvy, ctx_y);
				if (EC_TYPE == AR_NARY) {
					if (aff1_dmvy > num_symbol / 2 - (int) ((1 << subpel) * *pmvy))
						aff1_dmvy -= num_symbol;
					else if (aff1_dmvy < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvy))
						aff1_dmvy += num_symbol;
				}
			}
//affine V2
			if(fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) || 
			 (fmv1->merge_idx == MERGE && fmv1->merge_dir == LEFT) || (fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0) ){
				ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, 
				x+xblk2-1, y, info, t_level, blk_thresh);

				aff2_dmvx = ec_decode_word(ctx_x);   // decode the mvx 
				ec_update_model(aff2_dmvx, ctx_x);
				if (EC_TYPE == AR_NARY) {
					if (aff2_dmvx > num_symbol / 2 - (int) ((1 << subpel) * *pmvx))
						aff2_dmvx -= num_symbol;
					else if (aff2_dmvx < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvx))
						aff2_dmvx += num_symbol;
				}
      
				aff2_dmvy = ec_decode_word(ctx_y);   // decode the mvy 
				ec_update_model(aff2_dmvy, ctx_y);
				if (EC_TYPE == AR_NARY) {
					if (aff2_dmvy > num_symbol / 2 - (int) ((1 << subpel) * *pmvy))
						aff2_dmvy -= num_symbol;
					else if (aff2_dmvy < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvy))
						aff2_dmvy += num_symbol;
				}
			}
//affine V3
			if(fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ||
			  (fmv1->merge_idx == MERGE && fmv1->merge_dir == UP) || (fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0) ){
				ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field, 
							x, y+yblk2-1, info, t_level, blk_thresh);

				aff3_dmvx = ec_decode_word(ctx_x);   // decode the mvx 
				ec_update_model(aff3_dmvx, ctx_x);
				if (EC_TYPE == AR_NARY) {
					if (aff3_dmvx > num_symbol / 2 - (int) ((1 << subpel) * *pmvx))
						aff3_dmvx -= num_symbol;
					else if (aff3_dmvx < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvx))
						aff3_dmvx += num_symbol;
				}
      
				aff3_dmvy = ec_decode_word(ctx_y);   // decode the mvy 
				ec_update_model(aff3_dmvy, ctx_y);
				if (EC_TYPE == AR_NARY) {
					if (aff3_dmvy > num_symbol / 2 - (int) ((1 << subpel) * *pmvy))
					aff3_dmvy -= num_symbol;
					else if (aff3_dmvy < -(num_symbol / 2) - (int) ((1 << subpel) * *pmvy))
					aff3_dmvy += num_symbol;
				}
			}
		  }
		  else if(fmv1->direct_idx == DIRECT){

		  }
		}
	}else // layer structure, AGP or no AGP 
	{
	    ec_get_contexts_motion_field(&ctx_x, &ctx_y, frame_motion_field,
								     x, y, info, t_level, blk_thresh);
		dmvx = ec_decode_word(ctx_x);   // decode the mvx 
  	    ec_update_model(dmvx, ctx_x);
 	    dmvy = ec_decode_word(ctx_y);   // decode the mvy 
		ec_update_model(dmvy, ctx_y);
	}
	
	if ( info.AGP_level[t_level] ) // AGP 
	{
		AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
		subpel_scale = 1<<subpel; 
		major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
		major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;

		if ( !info.layer_mv[t_level] ) // no layer
		{
			if ((fmv1->bi_mode == PARALLEL)&&(!decode_parallelmv)) {
				fmv1->mvx = -fmv2->mvx;
				fmv1->mvy = -fmv2->mvy;
				major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
				major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			} else if ( fmv1->bi_mode == BLOCK_MERGING )
			{
				fmv1->mvx = mvpred_x;
				fmv1->mvy = mvpred_y;
				assert(dmvx ==0 && dmvy == 0 );
				major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
				major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			}else		
			{
				major_mvx = (float) dmvx / AGP_scale  + major_predx;
				major_mvy = (float) dmvy / AGP_scale  + major_predy;
				// code the sub-symbols as binary sequence 
				sub_symx = sub_symy = 0; 
				for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
				{
					sub_symx <<= 1; 
					if ( sub_bit < info.AGP_exist[t_level] )
					{
						value = get_splitted_mvbits( sub_bit, 0 );
						sub_symx += value;
					}
					sub_symy <<= 1; 
					if ( sub_bit < info.AGP_exist[t_level] ) 
					{
						value = get_splitted_mvbits( sub_bit, 0 );
						sub_symy += value; 
					}
				}
				sub_mvx = (float)sub_symx/ (float)subpel_scale;  
				sub_mvy = (float)sub_symy/ (float)subpel_scale;  
				if (info.AGP_exist[t_level] && major_mvx==0) 
				{
					sub_signx = get_splitted_signbits(0);  // 0: positive, 1: negative
					if (sub_signx) sub_mvx *= -1;
				}else if (major_mvx<0)
					sub_mvx *= -1;
				if (info.AGP_exist[t_level] && major_mvy==0) 
				{
					sub_signy = get_splitted_signbits(0);  // 0: positive, 1: negative
					if (sub_signy ) sub_mvy *= -1;
				}else if (major_mvy<0)
					sub_mvy *= -1;
				fmv1->mvx = major_mvx + sub_mvx;
				fmv1->mvy = major_mvy + sub_mvy; 
			}
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;

#ifdef  DEBUG_SCALABLE_MV
			fpAGP_debug = fopen(decoder_AGPdebug, "at");
			fprintf(fpAGP_debug, "x=%03d  y=%03d  blk=%2d ", x, y, xblk ); 
			fprintf(fpAGP_debug, "major_mvx=%.1f  major_mvy=%.1f  major_predx=%.1f  major_predy=%.1f",
					major_mvx,  major_mvy,  major_predx,  major_predy) ;
			fprintf(fpAGP_debug, " dmvx=%.1f  dmvy=%.1f  mvx=%.2f  mvy=%.2f, ctx_x=%d\t ctx_y=%d\n", 
								 fmv1->dmvx,  fmv1->dmvy, fmv1->mvx, fmv1->mvy, ctx_x, ctx_x) ;
			fclose(fpAGP_debug);
#endif 
		}

		if ( info.layer_mv[t_level] ) // layer structure 
		{
			major_mvx = (float) dmvx / AGP_scale  + major_predx;
			major_mvy = (float) dmvy / AGP_scale  + major_predy;
			// code the sub-symbols as binary sequence 
			sub_symx = sub_symy = 0; 
			for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
			{
				sub_symx <<= 1; 
				if ( sub_bit < info.AGP_exist[t_level] )
				{
					value = get_splitted_mvbits( sub_bit, 0 );
					sub_symx += value;
				}
				sub_symy <<= 1; 
				if ( sub_bit < info.AGP_exist[t_level] ) 
				{
					value = get_splitted_mvbits( sub_bit, 0 );
					sub_symy += value; 
				}
			}
			sub_mvx = (float)sub_symx/ (float)subpel_scale;  
			sub_mvy = (float)sub_symy/ (float)subpel_scale;  
			if (info.AGP_exist[t_level] && major_mvx==0) 
			{
				sub_signx = get_splitted_signbits(0);  // 0: positive, 1: negative
				if (sub_signx) sub_mvx *= -1;
			}else if (major_mvx<0)
				sub_mvx *= -1;
			if (info.AGP_exist[t_level] && major_mvy==0) 
			{
				sub_signy = get_splitted_signbits(0);  // 0: positive, 1: negative
				if (sub_signy ) sub_mvy *= -1;
			}else if (major_mvy<0)
				sub_mvy *= -1;
			fmv1->sample_mvx = fmv1->mvx = major_mvx + sub_mvx;
			fmv1->sample_mvy = fmv1->mvy = major_mvy + sub_mvy; 
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
	} // if ( info.AGP_level[t_level] ) // AGP 
	else if(!info.AGP_level[t_level]){
		fmv1->dmvx = (float) dmvx / (1 << subpel);
		fmv1->dmvy = (float) dmvy / (1 << subpel);

		if(fmv1->bi_mode >= 9 && fmv1->bi_mode <= 11){
			if(fmv1->direct_idx == INDIRECT){
				if( fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && 
				  fmv1->merge_dir == PAL_L && fmv1->aff_idx >= 0) || (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
					fmv1->aff1_dmvx = (float)aff1_dmvx/ (1 << subpel);
					fmv1->aff1_dmvy = (float)aff1_dmvy/ (1 << subpel);
				}
				if( (fmv1->merge_idx == MERGE && fmv1->merge_dir == LEFT) || fmv1->merge_idx == INTER ||
					(fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0) 
					|| (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
					fmv1->aff2_dmvx = (float)aff2_dmvx/ (1 << subpel);
					fmv1->aff2_dmvy = (float)aff2_dmvy/ (1 << subpel);
				}
				if( (fmv1->merge_idx == MERGE && fmv1->merge_dir == UP) || fmv1->merge_idx == INTER ||
				  (fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx >= 0)
				    || (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
					fmv1->aff3_dmvx = (float)aff3_dmvx/ (1 << subpel);
					fmv1->aff3_dmvy = (float)aff3_dmvy/ (1 << subpel);
				}
			}
		}
	}

	/////////////////////////
//	printf("bi_mode = %d, x = %d, y = %d, xblk = %d, yblk = %d, ctx_x = %d, ctx_y = %d\n fmv->dmvx = %f, fmv->dmvy = %f\n",
//		fmv1->bi_mode,x,y,xblk2,yblk2,ctx_x, ctx_y,fmv1->dmvx,fmv1->dmvy);
	/////////////////////////

    fmv1->is_predictor = YES;
	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
				fmv1->dmvy, 0, 0); 
    mvStat_setPos(x, y);
    mvStat_setDMV((float)dmvx, (float)dmvy);
    mvStat_writeDMVCTX();

#ifdef  DEBUG_LAYER_STRUCTURE  // for debug the layer structure coding 
    char  base_file[80]; 
    FILE *fpbase; 
    // make base_file and enhance_file empty
    sprintf(base_file, "decode_base%d.txt", count); 
    fpbase=fopen(base_file, "at");
    fprintf(fpbase, "x=%03d y=%03d blk=%d\t mvx=%.2f\t mvy=%.2f\t predx=%.2f predy=%.2f \n", 
  		    x, y, xblk, fmv1->mvx, fmv1->mvy, mvpred_x, mvpred_y); 
	fclose(fpbase);
#endif 
    
  }

}

/////////////////	Added by Yuan Liu on 01.30.2016   /////////////////////
/*****************************************************************************/
/*                            child_cand_decode()                            */
/*****************************************************************************/
void 
child_cand_decode( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count )
{
  int dmvx, dmvy, cx, cy, xblk2, yblk2;
  int ctx_x, ctx_y;
#ifdef MEDIAN_PREDICTION
  float mvpred_x, mvpred_y;
#endif
  char value; 
  float  major_mvx, major_mvy, sub_mvx, sub_mvy;
  float  major_predx, major_predy; 
  int    AGP_scale, subpel_scale, sub_symx, sub_symy, sub_bit, sub_signx, sub_signy; 
  int	enc_trans = 0;
  int   getnum;

  float mvpredstr_x[4], mvpredstr_y[4];
  vector_ptr mrg_left[4], mrg_right[4];
  int i, num = 0, num_dir = 0, mrg_num = 0;

  int   map_side;

  int	do_fmv1 = YES, do_fmv2 = YES;

  assert(fmv2 == NULL || fmv1->child == fmv2->child);

  xblk2 = ( x + xblk <= hor ) ? xblk : hor - x;
  yblk2 = ( y + yblk <= ver ) ? yblk : ver - y;

  if( fmv1->child && xblk > blk_thresh ) {    
    
	if(fmv2 != NULL)
		assert(fmv2->child);

    cx = x;
    cy = y;
    child_cand_decode(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    child_cand_decode(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    child_cand_decode(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    child_cand_decode(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);
  } else {
    if( x >= hor || y >= ver )      return;

	if ( !fmv1->child )  // there are no children for this block and size>=blk_thresh
	{
		if(fmv2 != NULL)
			assert( !fmv2->child );

		// no motion vector for this block on this side 
		if ((fmv1->lifting_mode == IGNORED))   do_fmv1 = NO;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) do_fmv1 = NO;
#endif 

		if(fmv2 != NULL){
			// no motion vector for this block on this side 
			if ((fmv2->lifting_mode == IGNORED))   do_fmv2 = NO;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
			if (fmv2->bi_mode==DIRECTIONAL_IBLOCK) do_fmv2 = NO;
#endif 
		}

	}

	///////////////	MERGE detection	///////////////////////
	if( fmv1->bi_mode == BLOCK_MERGING ||
	  (fmv1->bi_mode >= 9 && fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == DIR) ){

		if(fmv1->bi_mode == BLOCK_MERGING)
			assert(do_fmv1 == YES && do_fmv2 == YES);

		for(i=0;i<=3;i++){
			mrg_left[i] = new vector;
			mrg_right[i] = new vector;

			clean_mrg_mv(mrg_left[i]);
			clean_mrg_mv(mrg_right[i]);
		}

		if(fmv2 != NULL)
			get_field_merge_mv_info(mrg_left, mrg_right,x,y,xblk2,yblk2,hor,ver,info,t_level,frame_motion_field, frame_motion_field2, prev_frame_motion_field2, prev_frame_motion_field1,0);
		else
			get_field_merge_mv_info(mrg_left, mrg_right,x,y,xblk2,yblk2,hor,ver,info,t_level,frame_motion_field, frame_motion_field2, prev_frame_motion_field2, prev_frame_motion_field1,1);
		

		///////////////////////////
/*		for(i = 0;i <= 3;i++){
			printf("mrg_left[%d].mvx = %f, mrg_left[%d].mvy = %f\n",i,mrg_left[i]->mvx,i,mrg_left[i]->mvy);

			if(fmv2 != NULL){
				printf("mrg_right[%d].mvx = %f, mrg_right[%d].mvy = %f\n",i,mrg_right[i]->mvx,i,mrg_right[i]->mvy);
			}
		}
*/		///////////////////////////

		for(i=0;i<=3;i++){	
			if( (mrg_left[i]->mvx != (float)HUGE_VAL && mrg_left[i]->mvy != (float)HUGE_VAL) 
				|| (mrg_right[i]->mvx != (float)HUGE_VAL && mrg_right[i]->mvy != (float)HUGE_VAL) ){
				mrg_num ++;
			}
		}
//		printf("mrg_num = %d\n",mrg_num);

		///////////////////	MRG idx dec	//////////////////
		if(mrg_num == 0)
			assert(0);
		else if(mrg_num == 1){
			fmv1->med_idx = 0;

			if(fmv2!=NULL)
				fmv2->med_idx = 0;
		}else if(mrg_num == 2){
			fmv1->med_idx = getbits(1);

			if(fmv2 != NULL)
				fmv2->med_idx = fmv1->med_idx;
		}else{
			assert(mrg_num == 3 || mrg_num == 4);
			if(mrg_num == 3)
				input_huff_bits(fmv1,mrg_num,cnt_tri,info);
			else
				input_huff_bits(fmv1,mrg_num,cnt_quad,info);

			assert(fmv1->med_idx >= 0 && fmv1->med_idx <= 3 );

			if(fmv2!=NULL)
				fmv2->med_idx = fmv1->med_idx;
		}

		//////////////////////////////////////////////////
		if(fmv1->bi_mode == BLOCK_MERGING){
			if(mrg_left[fmv1->med_idx]->aff1_mvx != (float)HUGE_VAL || mrg_right[fmv1->med_idx]->aff1_mvx != (float)HUGE_VAL){
				assert( (mrg_left[fmv1->med_idx]->aff2_mvx != (float)HUGE_VAL && mrg_left[fmv1->med_idx]->aff3_mvx != (float)HUGE_VAL) ||
						(mrg_right[fmv1->med_idx]->aff2_mvx != (float)HUGE_VAL && mrg_right[fmv1->med_idx]->aff3_mvx != (float)HUGE_VAL) );

				fmv1->aff_mrg = getbits(1);

				if(fmv2 != NULL)
					fmv2->aff_mrg = fmv1->aff_mrg;
			}else{
				fmv1->aff_mrg = NO;

				if(fmv2 != NULL)
					fmv2->aff_mrg = NO;
			}

			assert(fmv1->lifting_mode == UNDECIDED);
			if(fmv2 != NULL)
				assert(fmv2->lifting_mode == UNDECIDED);
		}
		////////////////////////////////////////////

		if(fmv2 != NULL)
			assert( (mrg_left[fmv1->med_idx]->mvx != (float)HUGE_VAL && mrg_left[fmv1->med_idx]->mvy != (float)HUGE_VAL) ||
			(mrg_right[fmv2->med_idx]->mvx != (float)HUGE_VAL && mrg_right[fmv2->med_idx]->mvy != (float)HUGE_VAL) );

		if(mrg_left[fmv1->med_idx]->mvx != (float)HUGE_VAL && mrg_left[fmv1->med_idx]->mvy != (float)HUGE_VAL){
			fmv1->is_predictor = YES;
			fmv1->lifting_mode = CONNECTED;
		}else{
			fmv1->is_predictor = NO;
			fmv1->lifting_mode = IGNORED;
			do_fmv1 = NO;
		}

		if(fmv2 != NULL){
			if(mrg_right[fmv2->med_idx]->mvx != (float)HUGE_VAL && mrg_right[fmv2->med_idx]->mvy != (float)HUGE_VAL){
				fmv2->is_predictor = YES;
				fmv2->lifting_mode = CONNECTED;
			}else{
				fmv2->is_predictor = NO;
				fmv2->lifting_mode = IGNORED;
				do_fmv2 = NO;
			}
		}else{
			do_fmv2 = NO;
		}
	}
	///////////////////////////////////////////////////////


  if(do_fmv1 == YES){
	////////////  Added by Yuan Liu on 01.23.2016  //////////////
    for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
    }
    /////////////////////////////////////////////////////////////

	map_side = decode_parallelmv;

    get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field, prev_frame_motion_field2,
									  x, y, xblk, yblk, info, t_level, blk_thresh);

	/////////////////////////////////
	if(fmv1->bi_mode == BLOCK_MERGING || (fmv1->bi_mode >= 9 && fmv1->direct_idx == INDIRECT && 
		fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == DIR) ){
		mvpred_x = mrg_left[fmv1->med_idx]->mvx;
		mvpred_y = mrg_left[fmv1->med_idx]->mvy;

		if( fmv1->bi_mode == BLOCK_MERGING ){
			fmv1->aff1_mvx = mrg_left[fmv1->med_idx]->aff1_mvx;fmv1->aff1_mvy = mrg_left[fmv1->med_idx]->aff1_mvy;
			fmv1->aff2_mvx = mrg_left[fmv1->med_idx]->aff2_mvx;fmv1->aff2_mvy = mrg_left[fmv1->med_idx]->aff2_mvy;
			fmv1->aff3_mvx = mrg_left[fmv1->med_idx]->aff3_mvx;fmv1->aff3_mvy = mrg_left[fmv1->med_idx]->aff3_mvy;
		}

		assert(mvpred_x != (float)HUGE_VAL && mvpred_y != (float)HUGE_VAL);
	}else if( ( fmv1->bi_mode <= 8 && !(fmv1->bi_mode == PARALLEL && !map_side) && fmv1->bi_mode != BLOCK_MERGING ) || 
	(fmv1->bi_mode >= 9 && fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P && fmv1->trans_pred_idx == INDIR) ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL){
				num++;
				if(fmv1->bi_mode == BLOCK_MERGING){
				  if( (float(x)-mvpredstr_x[i] >= 0) && (float(x)-mvpredstr_x[i]<= (hor - xblk2)) && (float(y)-mvpredstr_y[i] >= 0)
					  && (float(y) - mvpredstr_y[i] <= (ver - yblk2)) )
						num_dir++;
				}
			}
		}
		assert(num_dir <= num);

//		assert( (fmv1->bi_mode >= 0 && fmv1->bi_mode <= 8) ||  );

		if(num == 0)
			fmv1->med_idx = -1;
		else if(num == 1){
			assert(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL);
			fmv1->med_idx = 0;
		}
		else if(num == 2){
			getnum = getbits(1);
			fmv1->med_idx = getnum;
		}
		else if(num >= 3){
			assert(num == 3 || num == 4);
			if(num == 3)
				input_huff_bits(fmv1,num,cnt_tri,info);
			else
				input_huff_bits(fmv1,num,cnt_quad,info);
		}
		else
			assert(0);

		if(num == 0){
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}else{
			if( (fmv1->bi_mode == BLOCK_MERGING) && num_dir == 0){
				assert(fmv1->med_idx == -1);
				mvpred_x = 0.0;
				mvpred_y = 0.0;
			}
			else if( (fmv1->bi_mode >= 0 && fmv1->bi_mode <= 8) || (fmv1->bi_mode >= 9) ){
				assert(fmv1->med_idx >= 0);
				mvpred_x = mvpredstr_x[fmv1->med_idx];
				mvpred_y = mvpredstr_y[fmv1->med_idx];
			}
		}
	}
	else if( (fmv1->bi_mode == PARALLEL && !map_side) ){

		if( (fmv1->bi_mode == PARALLEL && !map_side) )
			fmv1->med_idx = fmv2->med_idx;
		else
			assert(0);
	}
	/////////////////////////////////

	if (!info.AGP_level[t_level])  // no AGP 
	{
		if (!info.layer_mv[t_level])  // no layer structure in this temporal level 
		{
			if ( ((fmv1->bi_mode <=8) && (fmv1->bi_mode == PARALLEL)&&(!map_side)) ) {
				fmv1->mvx = -fmv2->mvx;
				fmv1->mvy = -fmv2->mvy;
			} else if( fmv1->bi_mode <=8 || 
			  (fmv1->bi_mode >= 9 && fmv1->direct_idx == INDIRECT && fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
				fmv1->mvx = fmv1->dmvx + mvpred_x;
				fmv1->mvy = fmv1->dmvy + mvpred_y;
			}

			//////////////	Added by Yuan Liu	//////////////////
			if(fmv1->bi_mode >= 7 && fmv1->bi_mode <= 11){
				if( (fmv1->bi_mode >= 7 && fmv1->aff_mrg == YES) || fmv1->bi_mode >= 9 )
					printf("\nbi_mode = %d, med_idx = %d, x = %d, y = %d, xblk = %d, yblk = %d\nfmv1->mvx = %f, fmv1->mvy = %f, fmv1->dmvx = %f, fmv1->dmvy = %f\n"
					,fmv1->bi_mode, fmv1->med_idx,x,y,xblk2,yblk2,fmv1->mvx,fmv1->mvy, fmv1->dmvx, fmv1->dmvy);
				
				if(fmv1->bi_mode == 7 && fmv1->aff_mrg == YES){
					printf("AFF_MRG!\n");
					printf("fmv1->aff1_mvx = %f, fmv1->aff1_mvy = %f\n fmv1->aff2_mvx = %f, fmv1->aff2_mvy = %f\n fmv1->aff3_mvx = %f, fmv1->aff3_mvy = %f\n",
						fmv1->aff1_mvx, fmv1->aff1_mvy, fmv1->aff2_mvx , fmv1->aff2_mvy, fmv1->aff3_mvx, fmv1->aff3_mvy);
				}

				if(fmv1->bi_mode>=9  && fmv1->direct_idx == DIRECT){

					fmv1->aff_idx = -1;

					get_aff_mvs(fmv1,x,y,xblk2,yblk2,hor,ver,frame_motion_field, prev_frame_motion_field2,1);
					printf("aff_index = %d\nfmv1->aff1_mvx = %f, fmv1->aff1_mvy = %f\n fmv1->aff2_mvx = %f, fmv1->aff2_mvy = %f\n fmv1->aff3_mvx = %f, fmv1->aff3_mvy = %f\n"
					,fmv1->aff_idx,fmv1->aff1_mvx, fmv1->aff1_mvy, fmv1->aff2_mvx , fmv1->aff2_mvy, fmv1->aff3_mvx, fmv1->aff3_mvy);
				}
				else if(fmv1->bi_mode>=9  && fmv1->direct_idx != DIRECT){
					get_aff_mvs(fmv1,x,y,xblk2,yblk2,hor,ver,frame_motion_field, prev_frame_motion_field2,1);

					//////////////////////////////////////////
					if(fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L) && fmv1->aff_idx == -1 ){
						fmv1->aff1_mvx = (-1) * fmv2->aff1_mvx;
						fmv1->aff1_mvy = (-1) * fmv2->aff1_mvy;
						fmv1->aff2_mvx = (-1) * fmv2->aff2_mvx;
						fmv1->aff2_mvy = (-1) * fmv2->aff2_mvy;
						fmv1->aff3_mvx = (-1) * fmv2->aff3_mvx;
						fmv1->aff3_mvy = (-1) * fmv2->aff3_mvy;
					}
					//////////////////////////////////////////

					printf("aff_index = %d\nfmv1->aff1_mvx = %f, fmv1->aff1_mvy = %f\n fmv1->aff2_mvx = %f, fmv1->aff2_mvy = %f\n fmv1->aff3_mvx = %f, fmv1->aff3_mvy = %f\nmerge_idx = %d, merge_dir = %d\n"
					,fmv1->aff_idx,fmv1->aff1_mvx, fmv1->aff1_mvy, fmv1->aff2_mvx , fmv1->aff2_mvy, fmv1->aff3_mvx, fmv1->aff3_mvy,fmv1->merge_idx,fmv1->merge_dir);

					if(fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P)
						printf("fmv1->trans_pred_idx = %d\n",fmv1->trans_pred_idx);

					if(fmv1->merge_idx == MERGE && fmv1->merge_dir == UP)
						printf("fmv1->aff3_dmvx = %f, fmv1->aff3_dmvy = %f\n", fmv1->aff3_dmvx, fmv1->aff3_dmvy);
					else if(fmv1->merge_idx == MERGE && fmv1->merge_dir == LEFT)
						printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n", fmv1->aff2_dmvx, fmv1->aff2_dmvy);
					else if(fmv1->merge_idx == INTER || (fmv1->merge_idx == MERGE && (fmv1->merge_dir == PAL_L)
						&& fmv1->aff_idx >=0) || (fmv1->merge_idx == MERGE && fmv1->merge_dir == TRAN_P) ){
						printf("\nfmv1->aff1_dmvx = %f, fmv1->aff1_dmvy = %f\n", fmv1->aff1_dmvx, fmv1->aff1_dmvy);
						printf("fmv1->aff2_dmvx = %f, fmv1->aff2_dmvy = %f\n", fmv1->aff2_dmvx, fmv1->aff2_dmvy);
						printf("fmv1->aff3_dmvx = %f, fmv1->aff3_dmvy = %f\n\n", fmv1->aff3_dmvx, fmv1->aff3_dmvy);
					}
				}
			}
//Added on 08.08.2018
			frame_mv_cnt ++;
			sum_mv2 += (fmv1->mvx*fmv1->mvx + fmv1->mvy*fmv1->mvy);
			sum_mv1 += sqrt(fmv1->mvx*fmv1->mvx + fmv1->mvy*fmv1->mvy);
			//////////////////////////////////////////////////////

			if(fmv1->bi_mode <= 6 || fmv1->bi_mode == 8 || (fmv1->bi_mode == 7 && fmv1->aff_mrg == NO) ){// Added by Yuan Liu
				assert( (x + xblk2 - fmv1->mvx <= hor ) &&
						(y + yblk2 - fmv1->mvy <= ver ) &&
						(x - fmv1->mvx >= 0 ) &&
						(y - fmv1->mvy >= 0 ) );
				if (!map_side) {
				  assert( (!(fmv1->bi_mode == PARALLEL)) || 
						  ((x + xblk2 - fmv2->mvx <= hor ) &&
						   (y + yblk2 - fmv2->mvy <= ver ) &&
						   (x - fmv2->mvx >= 0 ) &&
						   (y - fmv2->mvy >= 0 )) );
				}
			}//Added by Yuan Liu

		}else  // layer structure 
		{   
			fmv1->sample_mvx = fmv1->mvx = fmv1->dmvx + mvpred_x;
			fmv1->sample_mvy = fmv1->mvy = fmv1->dmvy + mvpred_y;
		}
	    
		if(fmv1->bi_mode >= 0 && fmv1->bi_mode <= 8){
			fmv1->dmvx = fmv1->mvx - mvpred_x;  // update the prediction error
			fmv1->dmvy = fmv1->mvy - mvpred_y;
		}
	}
	
	if ( info.AGP_level[t_level] ) // AGP 
	{
		assert(0);
		AGP_scale    = 1<< (subpel-info.AGP_level[t_level]);
		subpel_scale = 1<<subpel; 
		major_predx = (int)(mvpred_x*AGP_scale)/(float)AGP_scale;
		major_predy = (int)(mvpred_y*AGP_scale)/(float)AGP_scale;

		if ( !info.layer_mv[t_level] ) // no layer
		{
			if ((fmv1->bi_mode == PARALLEL)&&(!map_side)) {
				fmv1->mvx = -fmv2->mvx;
				fmv1->mvy = -fmv2->mvy;
				major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
				major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			} else if ( fmv1->bi_mode == BLOCK_MERGING )
			{
				fmv1->mvx = mvpred_x;
				fmv1->mvy = mvpred_y;
				assert(dmvx ==0 && dmvy == 0 );
				major_mvx = (int)(fmv1->mvx*AGP_scale)/(float)AGP_scale;
				major_mvy = (int)(fmv1->mvy*AGP_scale)/(float)AGP_scale;
			}else		
			{
				major_mvx = (float) dmvx / AGP_scale  + major_predx;
				major_mvy = (float) dmvy / AGP_scale  + major_predy;
				// code the sub-symbols as binary sequence 
				sub_symx = sub_symy = 0; 
				for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
				{
					sub_symx <<= 1; 
					if ( sub_bit < info.AGP_exist[t_level] )
					{
						value = get_splitted_mvbits( sub_bit, 0 );
						sub_symx += value;
					}
					sub_symy <<= 1; 
					if ( sub_bit < info.AGP_exist[t_level] ) 
					{
						value = get_splitted_mvbits( sub_bit, 0 );
						sub_symy += value; 
					}
				}
				sub_mvx = (float)sub_symx/ (float)subpel_scale;  
				sub_mvy = (float)sub_symy/ (float)subpel_scale;  
				if (info.AGP_exist[t_level] && major_mvx==0) 
				{
					sub_signx = get_splitted_signbits(0);  // 0: positive, 1: negative
					if (sub_signx) sub_mvx *= -1;
				}else if (major_mvx<0)
					sub_mvx *= -1;
				if (info.AGP_exist[t_level] && major_mvy==0) 
				{
					sub_signy = get_splitted_signbits(0);  // 0: positive, 1: negative
					if (sub_signy ) sub_mvy *= -1;
				}else if (major_mvy<0)
					sub_mvy *= -1;
				fmv1->mvx = major_mvx + sub_mvx;
				fmv1->mvy = major_mvy + sub_mvy; 
			}
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;

#ifdef  DEBUG_SCALABLE_MV
			fpAGP_debug = fopen(decoder_AGPdebug, "at");
			fprintf(fpAGP_debug, "x=%03d  y=%03d  blk=%2d ", x, y, xblk ); 
			fprintf(fpAGP_debug, "major_mvx=%.1f  major_mvy=%.1f  major_predx=%.1f  major_predy=%.1f",
					major_mvx,  major_mvy,  major_predx,  major_predy) ;
			fprintf(fpAGP_debug, " dmvx=%.1f  dmvy=%.1f  mvx=%.2f  mvy=%.2f, ctx_x=%d\t ctx_y=%d\n", 
								 fmv1->dmvx,  fmv1->dmvy, fmv1->mvx, fmv1->mvy, ctx_x, ctx_x) ;
			fclose(fpAGP_debug);
#endif 
		}

		if ( info.layer_mv[t_level] ) // layer structure 
		{
			major_mvx = (float) dmvx / AGP_scale  + major_predx;
			major_mvy = (float) dmvy / AGP_scale  + major_predy;
			// code the sub-symbols as binary sequence 
			sub_symx = sub_symy = 0; 
			for (sub_bit=0; sub_bit<info.AGP_level[t_level]; sub_bit++)
			{
				sub_symx <<= 1; 
				if ( sub_bit < info.AGP_exist[t_level] )
				{
					value = get_splitted_mvbits( sub_bit, 0 );
					sub_symx += value;
				}
				sub_symy <<= 1; 
				if ( sub_bit < info.AGP_exist[t_level] ) 
				{
					value = get_splitted_mvbits( sub_bit, 0 );
					sub_symy += value; 
				}
			}
			sub_mvx = (float)sub_symx/ (float)subpel_scale;  
			sub_mvy = (float)sub_symy/ (float)subpel_scale;  
			if (info.AGP_exist[t_level] && major_mvx==0) 
			{
				sub_signx = get_splitted_signbits(0);  // 0: positive, 1: negative
				if (sub_signx) sub_mvx *= -1;
			}else if (major_mvx<0)
				sub_mvx *= -1;
			if (info.AGP_exist[t_level] && major_mvy==0) 
			{
				sub_signy = get_splitted_signbits(0);  // 0: positive, 1: negative
				if (sub_signy ) sub_mvy *= -1;
			}else if (major_mvy<0)
				sub_mvy *= -1;
			fmv1->sample_mvx = fmv1->mvx = major_mvx + sub_mvx;
			fmv1->sample_mvy = fmv1->mvy = major_mvy + sub_mvy; 
			// only major symbols are median predicted, sub-symbols are coded by binary sequence
			fmv1->dmvx = major_mvx - major_predx;
			fmv1->dmvy = major_mvy - major_predy;
		}
	} // if ( info.AGP_level[t_level] ) // AGP 

    fmv1->is_predictor = YES;
	update_frame_motion_field(frame_motion_field, x, y, xblk, yblk, info, fmv1, fmv1->dmvx, 
				fmv1->dmvy, fmv1->mvx, fmv1->mvy); 

#ifdef  DEBUG_LAYER_STRUCTURE  // for debug the layer structure coding 
    char  base_file[80]; 
    FILE *fpbase; 
    // make base_file and enhance_file empty
    sprintf(base_file, "decode_base%d.txt", count); 
    fpbase=fopen(base_file, "at");
    fprintf(fpbase, "x=%03d y=%03d blk=%d\t mvx=%.2f\t mvy=%.2f\t predx=%.2f predy=%.2f \n", 
  		    x, y, xblk, fmv1->mvx, fmv1->mvy, mvpred_x, mvpred_y); 
	fclose(fpbase);
#endif 
    
  }//if do_fmv1

/**********************************************************/
/*****************   fmv2	*******************************/
/**********************************************************/

if(do_fmv2 == YES){

  if(fmv2 != NULL){

	////////////  Added by Yuan Liu on 01.23.2016  //////////////
    for(i=0;i<=3;i++){
	  mvpredstr_x[i] = (float)HUGE_VAL;
	  mvpredstr_y[i] = (float)HUGE_VAL;
    }
    /////////////////////////////////////////////////////////////

	map_side = (decode_parallelmv + 1) % 2;
	assert(map_side == 0);

	num = 0;
	num_dir = 0;

    get_median_predictor_motion_field(mvpredstr_x, mvpredstr_y, frame_motion_field2,prev_frame_motion_field1,
									  x, y, xblk, yblk, info, t_level, blk_thresh);

	/////////////////////////////////
	if(fmv2->bi_mode == BLOCK_MERGING || (fmv2->bi_mode >= 9 && fmv2->direct_idx == INDIRECT && 
		fmv2->merge_idx == MERGE && fmv2->merge_dir == TRAN_P && fmv2->trans_pred_idx == DIR) ){
		mvpred_x = mrg_right[fmv2->med_idx]->mvx;
		mvpred_y = mrg_right[fmv2->med_idx]->mvy;

		if(fmv2->bi_mode == BLOCK_MERGING){
			fmv2->aff1_mvx = mrg_right[fmv2->med_idx]->aff1_mvx;fmv2->aff1_mvy = mrg_right[fmv2->med_idx]->aff1_mvy;
			fmv2->aff2_mvx = mrg_right[fmv2->med_idx]->aff2_mvx;fmv2->aff2_mvy = mrg_right[fmv2->med_idx]->aff2_mvy;
			fmv2->aff3_mvx = mrg_right[fmv2->med_idx]->aff3_mvx;fmv2->aff3_mvy = mrg_right[fmv2->med_idx]->aff3_mvy;
		}

		assert(mvpred_x != (float)HUGE_VAL && mvpred_y != (float)HUGE_VAL);
	}else if( (fmv2->bi_mode <= 8 && !(fmv2->bi_mode == PARALLEL && !map_side && fmv2->bi_mode != BLOCK_MERGING) ) ||
	  (fmv2->bi_mode >= 9 && fmv2->direct_idx == INDIRECT && fmv2->merge_idx == MERGE && fmv2->merge_dir == TRAN_P && fmv2->trans_pred_idx == INDIR) ){
		for(i=0;i<=3;i++){
			if(mvpredstr_x[i] != (float)HUGE_VAL && mvpredstr_y[i] != (float)HUGE_VAL){
				num++;
				if(fmv2->bi_mode == BLOCK_MERGING){
				  if( (float(x)-mvpredstr_x[i] >= 0) && (float(x)-mvpredstr_x[i]<= (hor - xblk2)) && (float(y)-mvpredstr_y[i] >= 0)
					  && (float(y) - mvpredstr_y[i] <= (ver - yblk2)) )
						num_dir++;
				}
			}
		}
		assert(num_dir <= num);

		if(num == 0)
			fmv2->med_idx = -1;
		else if(num == 1){
			assert(mvpredstr_x[0] != (float)HUGE_VAL && mvpredstr_y[0] != (float)HUGE_VAL);
			fmv2->med_idx = 0;
		}
		else if(num == 2){
			getnum = getbits(1);
			fmv2->med_idx = getnum;
		}
		else if(num >= 3){
			assert(num == 3 || num == 4);
			if(num == 3)
				input_huff_bits(fmv2,num,cnt_tri,info);
			else
				input_huff_bits(fmv2,num,cnt_quad,info);
		}
		else
			assert(0);

		if(num == 0){
			mvpred_x = 0.0;
			mvpred_y = 0.0;
		}
		else{
			if( (fmv2->bi_mode == BLOCK_MERGING) && num_dir == 0){
				assert(fmv2->med_idx == -1);
				mvpred_x = 0.0;
				mvpred_y = 0.0;
			}
			else if( (fmv2->bi_mode >= 0 && fmv2->bi_mode <= 8) || (fmv2->bi_mode >= 9) ){
				assert(fmv2->med_idx >= 0);
				mvpred_x = mvpredstr_x[fmv2->med_idx];
				mvpred_y = mvpredstr_y[fmv2->med_idx];
			}
		}
	}
	else if( (fmv2->bi_mode == PARALLEL && !map_side) ){

		if( (fmv2->bi_mode == PARALLEL && !map_side) )
			fmv2->med_idx = fmv2->med_idx;
		else
			assert(0);
	}
	/////////////////////////////////

	if (!info.AGP_level[t_level])  // no AGP 
	{
		if (!info.layer_mv[t_level])  // no layer structure in this temporal level 
		{
			if ( ((fmv2->bi_mode <=8) && (fmv2->bi_mode == PARALLEL)&&(!map_side)) ) {
				fmv2->mvx = -fmv1->mvx;
				fmv2->mvy = -fmv1->mvy;
			} else if( fmv2->bi_mode <=8 || (fmv2->bi_mode >= 9 && fmv2->direct_idx == INDIRECT &&
				fmv2->merge_idx == MERGE && fmv2->merge_dir == TRAN_P) ) {
				fmv2->mvx = fmv2->dmvx + mvpred_x;
				fmv2->mvy = fmv2->dmvy + mvpred_y;
			}

			//////////////	Added by Yuan Liu	//////////////////
			if(fmv2->bi_mode >= 7 && fmv2->bi_mode <= 11){
				if( (fmv2->bi_mode >= 7 && fmv2->aff_mrg == YES) || fmv2->bi_mode >= 9 )
					printf("\nbi_mode = %d, med_idx = %d, x = %d, y = %d, xblk = %d, yblk = %d\nfmv2->mvx = %f, fmv2->mvy = %f, fmv2->dmvx = %f, fmv2->dmvy = %f\n"
					,fmv2->bi_mode, fmv2->med_idx,x,y,xblk2,yblk2,fmv2->mvx,fmv2->mvy, fmv2->dmvx, fmv2->dmvy);

				if(fmv2->bi_mode == 7 && fmv2->aff_mrg == YES){
					printf("AFF_MRG!\n");
					printf("fmv2->aff1_mvx = %f, fmv2->aff1_mvy = %f\n fmv2->aff2_mvx = %f, fmv2->aff2_mvy = %f\n fmv2->aff3_mvx = %f, fmv2->aff3_mvy = %f\n",
						fmv2->aff1_mvx, fmv2->aff1_mvy, fmv2->aff2_mvx , fmv2->aff2_mvy, fmv2->aff3_mvx, fmv2->aff3_mvy);
				}

				if(fmv2->bi_mode>=9  && fmv2->direct_idx == DIRECT){
					
					fmv2->aff_idx = -1;

					get_aff_mvs(fmv2,x,y,xblk2,yblk2,hor,ver,frame_motion_field2,prev_frame_motion_field1,1);
					printf("aff_index = %d\nfmv2->aff1_mvx = %f, fmv2->aff1_mvy = %f\n fmv2->aff2_mvx = %f, fmv2->aff2_mvy = %f\n fmv2->aff3_mvx = %f, fmv2->aff3_mvy = %f\n"
					,fmv2->aff_idx,fmv2->aff1_mvx, fmv2->aff1_mvy, fmv2->aff2_mvx , fmv2->aff2_mvy, fmv2->aff3_mvx, fmv2->aff3_mvy);
				}
				else if(fmv2->bi_mode>=9  && fmv2->direct_idx != DIRECT){

					if(fmv2->direct_idx == INDIRECT && fmv2->merge_idx == MERGE && fmv2->merge_dir == PAL_L)
						get_aff_mvs(fmv2,x,y,xblk2,yblk2,hor,ver,frame_motion_field2,prev_frame_motion_field1,0);
					else
						get_aff_mvs(fmv2,x,y,xblk2,yblk2,hor,ver,frame_motion_field2,prev_frame_motion_field1,1);

					//////////////////////////////////////////
					if(fmv2->merge_idx == MERGE && fmv2->merge_dir == PAL_L && fmv2->aff_idx == -1 ){
						fmv2->aff1_mvx = (-1) * fmv1->aff1_mvx;
						fmv2->aff1_mvy = (-1) * fmv1->aff1_mvy;
						fmv2->aff2_mvx = (-1) * fmv1->aff2_mvx;
						fmv2->aff2_mvy = (-1) * fmv1->aff2_mvy;
						fmv2->aff3_mvx = (-1) * fmv1->aff3_mvx;
						fmv2->aff3_mvy = (-1) * fmv1->aff3_mvy;
					}
					//////////////////////////////////////////

					printf("aff_index = %d\nfmv2->aff1_mvx = %f, fmv2->aff1_mvy = %f\n fmv2->aff2_mvx = %f, fmv2->aff2_mvy = %f\n fmv2->aff3_mvx = %f, fmv2->aff3_mvy = %f\nmerge_idx = %d, merge_dir = %d\n"
					,fmv2->aff_idx,fmv2->aff1_mvx, fmv2->aff1_mvy, fmv2->aff2_mvx , fmv2->aff2_mvy, fmv2->aff3_mvx, fmv2->aff3_mvy,fmv2->merge_idx,fmv2->merge_dir);

					if(fmv2->merge_idx == MERGE && fmv2->merge_dir == TRAN_P)
						printf("fmv2->trans_pred_idx = %d\n",fmv2->trans_pred_idx);

					if(fmv2->merge_idx == MERGE && fmv2->merge_dir == UP)
						printf("fmv2->aff3_dmvx = %f, fmv2->aff3_dmvy = %f\n", fmv2->aff3_dmvx, fmv2->aff3_dmvy);
					else if(fmv2->merge_idx == MERGE && fmv2->merge_dir == LEFT)
						printf("fmv2->aff2_dmvx = %f, fmv2->aff2_dmvy = %f\n", fmv2->aff2_dmvx, fmv2->aff2_dmvy);
					else if(fmv2->merge_idx == INTER || (fmv2->merge_idx == MERGE && (fmv2->merge_dir == PAL_L)
						&& fmv2->aff_idx >=0) || (fmv2->merge_idx == MERGE && fmv2->merge_dir == TRAN_P) ){
						printf("\nfmv2->aff1_dmvx = %f, fmv2->aff1_dmvy = %f\n", fmv2->aff1_dmvx, fmv2->aff1_dmvy);
						printf("fmv2->aff2_dmvx = %f, fmv2->aff2_dmvy = %f\n", fmv2->aff2_dmvx, fmv2->aff2_dmvy);
						printf("fmv2->aff3_dmvx = %f, fmv2->aff3_dmvy = %f\n\n", fmv2->aff3_dmvx, fmv2->aff3_dmvy);
					}
				}
			}
			//////////////////////////////////////////////////////

			if(fmv2->bi_mode <= 6 || fmv2->bi_mode == 8 || (fmv2->bi_mode == 7 && fmv2->aff_mrg == NO) ){// Added by Yuan Liu
				assert( (x + xblk2 - fmv2->mvx <= hor ) &&
						(y + yblk2 - fmv2->mvy <= ver ) &&
						(x - fmv2->mvx >= 0 ) &&
						(y - fmv2->mvy >= 0 ) );
				if (!map_side) {
				  assert( (!(fmv2->bi_mode == PARALLEL)) || 
						  ((x + xblk2 - fmv1->mvx <= hor ) &&
						   (y + yblk2 - fmv1->mvy <= ver ) &&
						   (x - fmv1->mvx >= 0 ) &&
						   (y - fmv1->mvy >= 0 )) );
				}
			}//Added by Yuan Liu

		}else  // layer structure 
		{   
			fmv2->sample_mvx = fmv2->mvx = fmv2->dmvx + mvpred_x;
			fmv2->sample_mvy = fmv2->mvy = fmv2->dmvy + mvpred_y;
		}
	    
		if(fmv2->bi_mode >= 0 && fmv2->bi_mode <= 8){
			fmv2->dmvx = fmv2->mvx - mvpred_x;  // update the prediction error
			fmv2->dmvy = fmv2->mvy - mvpred_y;
		}
	}

  
	fmv2->is_predictor = YES;
	update_frame_motion_field(frame_motion_field2, x, y, xblk, yblk, info, fmv2, fmv2->dmvx, 
				fmv2->dmvy, fmv2->mvx, fmv2->mvy); 

  }//if fmv2 != NULL

}//if do_fmv2

if( fmv1->bi_mode == BLOCK_MERGING ){
	for(i=0;i<=3;i++){
		delete(mrg_left[i]);
		delete(mrg_right[i]);
	}
}

}//if child else

}

/****************************************************************************/
/*                              mv_decode()    解码mv 解码一帧                             */
/****************************************************************************/
long int
mv_decode( vector_ptr fmv1, vector_ptr fmv2, int decode_map, 
           int meandepth, videoinfo info, int large, int t_level, long int starting_pos,
		   int bidir_exist, int count, int GOP_counter)
{
  int x, y, X, Y, xnum, ynum, xblk, yblk, hor, ver, pos, small, itemp, subpel;
  int num_symbol;
  float pmvx, pmvy;
  long int mvBytes = 0, major_bytes, initial_pos, enhance_bytes;
  int sub_bit, blk_thresh; 
  int i,j;

  if(decode_map == 0 || fmv2 == NULL){
	  frame_mv_cnt = 0;
	  sum_mv2 = 0;
	  sum_mv1 = 0;
  }

  if(fmv2 != NULL){
    for(i=0;i<=11;i++){
	  bi_mode_num012[i]=0;
	  bi_mode_num345[i]=0;
	}
	use_huff = 0;
  }

  // for motion field in a frame 申请运动场
  frame_motion_field = (FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(FRAME_MOTION_FIELD), "frame_motion_field"); 
  frame_motion_field2 = (FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(FRAME_MOTION_FIELD), "frame_motion_field2"); 

  blk_thresh = 0; 
  if (info.layer_mv[t_level])
	blk_thresh = LAYER_BLOCK_SIZE; 

#ifdef DEBUG_LAYER_STRUCTURE
  char  base_file[80], enhance_file[80]; 
  FILE *fpbase, *fpenhance; 
  // make base_file and enhance_file empty
  sprintf(base_file, "decode_base%d.txt", count); 
  if (fpbase=fopen(base_file, "rt") )  
  {
	  fclose(fpbase);
	  fpbase=fopen(base_file, "wt");
	  fclose(fpbase); 
  }
  sprintf(enhance_file, "decode_enhance%d.txt", count); 
  if (fpenhance=fopen(enhance_file, "rt") )
  {
	  fclose(fpenhance);
	  fpenhance=fopen(enhance_file, "wt");
	  fclose(fpenhance); 
  }
#endif

  initial_pos = starting_pos; // for scalable motion vector coding 为可伸缩运动矢量编码

  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];
  subpel = info.subpel[t_level];
  hor = info.ywidth;
  ver = info.yheight;
  small = xblk;
  itemp = info.level[t_level];
  while( itemp != 1 ) {
    small /= 2;
    itemp--;
  }

#ifdef DEBUG_SCALABLE_MV  
  sprintf(decoder_AGPdebug, "decoder_AGP_debugfile_GOP%03d_count%03d.txt", GOP_counter, count);
  debug_counter++;
#endif

  num_symbol = 0;
  
  if (bidir_exist)  // the sign whethe the motion vectors on this side exist 标示在这边运动向量是否存在
  {
	  fseek(fpbit, initial_pos, SEEK_SET); 
	  major_bytes = ec_dec_preinit(fpbit);  // the number of bytes for major symbols in base layer 字节的数量为主要符号在基本层
	  starting_pos += major_bytes; 
	  if (info.AGP_exist[t_level])// 子像素存在多少
	  {
		  // sign bytes
		  get_splitted_bytes(store_splitted_sign[0], &(splitted_sign_byte_num[0]), info, t_level, fpbit, starting_pos);
		  starting_pos += splitted_sign_byte_num[0]+2;
		  splitted_sign_bit_num[0] = 0;
		  splitted_sign_byte_num[0] = 0;  // reset to starting position 0 
		  // get the sub-symbol bits in the order of significance 根据重要性获得子符号的bit
		  for (sub_bit=0; sub_bit<info.AGP_exist[t_level]; sub_bit++)  
		  {
			  // sub_bit=0 is the most significant sub-symbol 零是最重要的系数
			  get_splitted_bytes(store_splitted_bytes_array[0][sub_bit], 
								 &(splitted_mv_byte_num_array[0][sub_bit]), info, t_level, fpbit, starting_pos);
			  starting_pos += splitted_mv_byte_num_array[0][sub_bit]+2;
			  splitted_bit_num_array[0][sub_bit] = 0;
			  splitted_mv_byte_num_array[0][sub_bit] = 0; // reset to starting position 0 
		  }
	  }
	  fseek(fpbit, initial_pos+4, SEEK_SET); // reset the pointer to initial starting position for this mv set
  }

  if(fmv2 != NULL){
	use_huff = getbits(1);
  }

  if( decode_map ){
    for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {     
      for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
        pos = Y * xnum + X;
        child_map_decode( &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL,
                          meandepth, x, y, xblk, yblk, hor, ver,
                          small, info, t_level);
      }
    }
  }

//  printf("Exited from map!\n");

  if (info.layer_mv[t_level] )  // reconstruct the subsample procedure for layer structure coding 重建下采样结构为了层结构编码
	  layer_structure_mv_subsample(info, t_level, 0,  fmv1, 0);

  if (info.layer_mv[t_level])
  {
    for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {     
      for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
        pos = Y * xnum + X;
        child_merge_decode( &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL,
                          meandepth, x, y, xblk, yblk, hor, ver,
                          small, info, t_level, blk_thresh, count);
      }
    }
  }

  clear_frame_motion_field(frame_motion_field, info); // 初始化
  if (bidir_exist)
  {
#ifdef EC_USE_CONTEXTS
	 ec_dec_init(EC_TYPE, num_symbol, 3);
#else
	 ec_dec_init(EC_TYPE, num_symbol, 1);
#endif
  }

  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
    pmvx = 0.;
    pmvy = 0.;                  
    for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
      pos = Y * xnum + X;
      child_mv_decode( fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL,
                       &pmvx, &pmvy, num_symbol, subpel,
                       x, y, xblk, yblk, hor, ver, info, decode_map, t_level, bidir_exist, 
					   blk_thresh, count );  // add decode_map for mv coding in parallel mode. mwi 
    }
  }

  ec_dec_end1();

//  getbits(14);

  if( (fmv2 != NULL && decode_map == 1) ){
	copy_frame_motion_field2(prev_frame_motion_field2, prev_frame_motion_field_left2, info);
	copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field_left1, info);
  }

  if( (fmv2 != NULL && decode_map == 0) ){
	copy_frame_motion_field2(prev_frame_motion_field_left2, prev_frame_motion_field2, info);
	copy_frame_motion_field2(prev_frame_motion_field_left1, prev_frame_motion_field1, info);
  }

  if( fmv2 != NULL && decode_map == 0 ){

	  if( fmv2 != NULL && decode_map == 0 ){

		clear_frame_motion_field(frame_motion_field, info); 

		clear_frame_motion_field(frame_motion_field2, info); 

		for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
			pmvx = 0.;
			pmvy = 0.;                  
			for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
			  pos = Y * xnum + X;
			  child_cand_decode( fmv2, fmv1, &fmv2[pos], fmv1 ? (&fmv1[pos]) : NULL,
							   &pmvx, &pmvy, num_symbol, subpel,
							   x, y, xblk, yblk, hor, ver, info, 1, t_level, bidir_exist, 
							   blk_thresh, count );
			}
		}
	  }

	  copy_frame_motion_field(frame_motion_field, prev_frame_motion_field2, info);
	  copy_frame_motion_field(frame_motion_field2, prev_frame_motion_field1, info);

  }// if( fmv2 != NULL && decode_map == 0 )
  else if( fmv2 == NULL ){

	clear_frame_motion_field(frame_motion_field2, info); 

	clear_frame_motion_field(frame_motion_field, info); 

	for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
		pmvx = 0.;
		pmvy = 0.;                  
		for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
		  pos = Y * xnum + X;
		  child_cand_decode( fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL,
						   &pmvx, &pmvy, num_symbol, subpel,
						   x, y, xblk, yblk, hor, ver, info, 1, t_level, bidir_exist, 
						   blk_thresh, count );
		}
	}

	copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
    copy_frame_motion_field(frame_motion_field, prev_frame_motion_field1, info);

  }// if fmv2 == NULL
  else{
	assert( fmv2 != NULL && decode_map == 1 );
	
	copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
	copy_frame_motion_field(frame_motion_field, prev_frame_motion_field1, info);
  }

  // replace the sub-tree after subsampling
  if (info.layer_mv[t_level])
  {
	  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
		pmvx = 0.;
		pmvy = 0.;                  
		for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
			pos = Y * xnum + X;
		    replace_enhancement_mv( fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL,
                       &pmvx, &pmvy, num_symbol, subpel,
                       x, y, xblk, yblk, hor, ver, info, decode_map, t_level, bidir_exist, 
					   blk_thresh, count );  
		}
	  }
  }

  if (bidir_exist)
  {
	  mvBytes = major_bytes = ec_dec_end2();
 	  if (info.AGP_exist[t_level])
	  {
		  mvBytes += splitted_sign_byte_num[0]+2;
		  for (sub_bit=0; sub_bit<info.AGP_exist[t_level]; sub_bit++)
			  mvBytes += splitted_mv_byte_num_array[0][sub_bit]+2;	
	  }
  }else
	return 0; 

  if ( info.layer_exist[t_level] )
  {
	  // starting position for enhancement layer
	  initial_pos = starting_pos = initial_pos+mvBytes; 
	  fseek(fpbit, initial_pos, SEEK_SET);
	  enhance_bytes = ec_dec_preinit(fpbit);  
	  starting_pos += enhance_bytes;
	  if (info.AGP_exist[t_level])
	  {
		  // sign bytes
		  get_splitted_bytes(store_splitted_sign[1], &(splitted_sign_byte_num[1]), info, t_level, fpbit, starting_pos);
		  starting_pos += splitted_sign_byte_num[1]+2;
		  splitted_sign_bit_num[1] = 0;
		  splitted_sign_byte_num[1] = 0;  // reset to starting position 0 
		  // get the sub-symbol bits in the order of significance 
		  for (sub_bit=0; sub_bit<info.AGP_exist[t_level]; sub_bit++)  
		  {
			  // sub_bit=0 is the most significant sub-symbol
			  get_splitted_bytes(store_splitted_bytes_array[1][sub_bit], 
								 &(splitted_mv_byte_num_array[1][sub_bit]), info, t_level, fpbit, starting_pos);
			  starting_pos += splitted_mv_byte_num_array[1][sub_bit]+2;
			  splitted_bit_num_array[1][sub_bit] = 0;
			  splitted_mv_byte_num_array[1][sub_bit] = 0; // reset to starting position 0 
		  }
	  }
	  fseek(fpbit, initial_pos+4, SEEK_SET); // reset the pointer to initial starting position for this mv set

	  clear_frame_motion_field(frame_motion_field, info); 
	  ec_dec_init(EC_TYPE, num_symbol, 3);
	  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
		pmvx = 0.;
		pmvy = 0.;                  
		for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
		  pos = Y * xnum + X;
		  child_mv_decode_enhance( fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL,
						   &pmvx, &pmvy, num_symbol, subpel,
						   x, y, xblk, yblk, hor, ver, info, decode_map, t_level, bidir_exist, 
						   blk_thresh, count );  
		}
	  }
	  ec_dec_end1();
	  enhance_bytes = ec_dec_end2();
	  mvBytes += enhance_bytes; 
 	  if (info.AGP_exist[t_level])
	  {
		  mvBytes += splitted_sign_byte_num[1]+2;
		  for (sub_bit=0; sub_bit<info.AGP_exist[t_level]; sub_bit++)
			  mvBytes += splitted_mv_byte_num_array[1][sub_bit]+2;	
	  }
  }

  free(frame_motion_field); 
  free(frame_motion_field2); 

//Added on 08.08.2018
  if(decode_map == 0 || fmv2 == NULL){
	  if(frame_mv_cnt > 0){
		sum_mv2 /= frame_mv_cnt;
		sum_mv1 /= frame_mv_cnt;
		sum_mv2 = sum_mv2 - sum_mv1*sum_mv1;
		
		printf("VAR = %f, frame_mv_cnt = %d\n",sum_mv2, frame_mv_cnt);
	  }
  }

  return mvBytes;  // the number of motion vector bytes for this set of motion vectors 


}


/*****************************************************************************/
/*                              mv_decoding() 解码mv                              */
/*****************************************************************************/
long int
mv_decoding( videoinfo info, vector_ptr * yfmv, int GOP_counter, long int starting_pos, int simul_dec, int theo_dec )
{
  int i, j, count, dist, meandepth, GOPsz;
  int eff_GOPsz[20]; // maximum number of temporal levels = 20;
  long int mvBytes = 0;  // the total number of bytes for all sets of motion vectors
  long int cur_bytes;    // the number of bytes for current set of motion vectors

  for(i=0;i<=2;i++){
	cnt_tri[i]=0;
  }
  for(i=0;i<=3;i++){
	cnt_quad[i]=0;
  }

  for(i=0;i<27;i++)
	aff_cum_idx[i]=0;

  // starting_pos: the starting position for current set of motion vectors in the video bitstream

  prev_frame_motion_field1 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, // 先前帧的运动场
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field1"); 
  prev_frame_motion_field2 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field2");

  //Added on 08.20.2016
  prev_frame_motion_field_left1 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, // 先前帧的运动场左边
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field_left1"); 
  prev_frame_motion_field_left2 = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
							sizeof(SIMP_FRAME_MOTION_FIELD), "prev_frame_motion_field_left2"); 
  //Added on 08.20.2016

  if(GOP_counter == 0 && simul_dec == NO // 第一个gop。申请空间并初始化
//	  && theo_dec == YES
	  )
  {
	  for(i = 0;i < G_LEVEL; i++){
		buffer_frame_motion_field1_dec[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "buffer_frame_motion_field1_dec");
		buffer_frame_motion_field2_dec[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "buffer_frame_motion_field2_dec");

		clear_frame_motion_field_simp(buffer_frame_motion_field1_dec[i], info);
		clear_frame_motion_field_simp(buffer_frame_motion_field2_dec[i], info);

//Added on 02.05.2018
		save_buffer_frame_motion_field1_dec[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "buffer_frame_motion_field1_dec");
		save_buffer_frame_motion_field2_dec[i] = (SIMP_FRAME_MOTION_FIELD*)getarray(info.yheight*info.ywidth, 
								sizeof(SIMP_FRAME_MOTION_FIELD), "buffer_frame_motion_field2_dec");

		clear_frame_motion_field_simp(save_buffer_frame_motion_field1_dec[i], info);
		clear_frame_motion_field_simp(save_buffer_frame_motion_field2_dec[i], info);
	  }
  }

//Added on 02.08.2018
  if(simul_dec == YES){
	  for(i = 0; i <= info.tPyrLev - 2; i++){
		copy_frame_motion_field2(save_buffer_frame_motion_field1_dec[i], buffer_frame_motion_field1_dec[i], info);
		copy_frame_motion_field2(save_buffer_frame_motion_field2_dec[i], buffer_frame_motion_field2_dec[i], info);
	  }
  }

  if((info.GOPsz * GOP_counter) == SIMUL_POINT && theo_dec == NO)
  {
	  for(i = 0; i <= info.tPyrLev - 2; i++){
		copy_frame_motion_field2(buffer_frame_motion_field1_dec[i], save_buffer_frame_motion_field1_dec[i], info);
		copy_frame_motion_field2(buffer_frame_motion_field2_dec[i], save_buffer_frame_motion_field2_dec[i], info);
	  }
  }

  
  GOPsz = info.GOPsz;
  dist  = ( int )pow( 2.0, ( double )( info.tPyrLev - 1 ) );

  // determine effective GOP size in level i
  for( i = 0; i < info.tPyrLev; i++ ) {
    eff_GOPsz[i] = GOPsz;  
    GOPsz /= 2;
  }
  
  count = 1;
  for( i = info.tPyrLev - 1; i >= 0; i-- ) {// 逐时域层向下减
    meandepth = info.pixeldepth + ( i+1 + 1 ) / 2;  // i+1 is the tPyrLev  位深，越高的位深越大
    count++;

	printf("LEVEL CHANGED!\n");
	clear_frame_motion_field_simp(prev_frame_motion_field1, info);
	clear_frame_motion_field_simp(prev_frame_motion_field2, info);

	//Added on 08.20.2016
	clear_frame_motion_field_simp(prev_frame_motion_field_left1, info);
	clear_frame_motion_field_simp(prev_frame_motion_field_left2, info);
	//Added on 08.20.2016

	if(GOP_counter >= 1 && i <= info.tPyrLev - 2){
		copy_frame_motion_field2(buffer_frame_motion_field1_dec[i], prev_frame_motion_field1, info);
		copy_frame_motion_field2(buffer_frame_motion_field2_dec[i], prev_frame_motion_field2, info);
	}

    for( j = 1; j <= eff_GOPsz[i]; j += 2 ) // 一帧一帧去做
	{
      if (j == eff_GOPsz[i]) // 最后一帧
	  {  // single MVF
		printf("single MVF detected!\n");

        if( i >= info.t_level ) {
          if (dec_scene_change[i][j] == NO) {
            mvStat_setFrame(GOP_counter, i, j);
            cur_bytes = mv_decode(yfmv[count], NULL, 1,
                                 meandepth, info, dist, i, starting_pos, 1, count, GOP_counter );
			mvBytes += cur_bytes;   starting_pos += cur_bytes; 
			if ( info.layer_mv[i] && !info.layer_exist[i] ) // enhancement layer is discarded
				layer_structure_mv_trim(info, i, j, yfmv[count], count);
#ifdef  DEBUG_BLOCK_MODE_MV_INFO
		   write_block_mode_motion_vector("decoderleft", GOP_counter, count, 1, yfmv[count], info, i);
#endif 
          }
        }
        count++;
      } 
	  else 
	  {
        if( i >= info.t_level ) {
          if (dec_scene_change[i][j] == NO && dec_scene_change[i][j + 1] == NO) // 左右都可参考
		  {
			printf("\n\nBI MC, i = %d, j = %d, eff_GOPsz = %d, count = %d, count+1 = %d\n",i,j,eff_GOPsz[i], count,count+1);
            mvStat_setFrame(GOP_counter, i, j);
            cur_bytes = mv_decode(yfmv[count], yfmv[count + 1], 1,
                                 meandepth, info, dist, i, starting_pos, 1, count, GOP_counter );
			mvBytes += cur_bytes;   
			starting_pos += cur_bytes; 
            mvStat_setFrame(GOP_counter, i, j + 1);
			printf("\n\nBI MC, i = %d, j = %d, eff_GOPsz = %d, count+1 = %d, count = %d\n",i,j,eff_GOPsz[i], count+1,count);
            cur_bytes = mv_decode(yfmv[count + 1], yfmv[count], 0,
                                 meandepth, info, dist, i, starting_pos, info.bi_exist[i], count+1, GOP_counter );
			mvBytes += cur_bytes;   starting_pos += cur_bytes; 
			if ( info.layer_mv[i] && !info.layer_exist[i] ) // enhancement layer is discarded
			{
				layer_structure_mv_trim(info, i, j, yfmv[count], count);
				layer_structure_mv_trim(info, i, j+1, yfmv[count+1], count+1);
			}
#ifdef  DEBUG_BLOCK_MODE_MV_INFO
			write_block_mode_motion_vector("decoderleft", GOP_counter, count, 0, yfmv[count], info, i);
			write_block_mode_motion_vector("decoderright",GOP_counter, count+1, 0, yfmv[count+1], info, i);
#endif

          } 
		  else if (dec_scene_change[i][j] == NO) {
			printf("\n\nLEFT MC, i = %d, j = %d, eff_GOPsz = %d, count = %d\n",i,j,eff_GOPsz[i], count);

            mvStat_setFrame(GOP_counter, i, j);
            cur_bytes = mv_decode(yfmv[count], NULL, 1,
                                 meandepth, info, dist, i, starting_pos, 1, count, GOP_counter );
			mvBytes += cur_bytes;   starting_pos += cur_bytes; 
			if ( info.layer_mv[i] && !info.layer_exist[i] ) // enhancement layer is discarded
				layer_structure_mv_trim(info, i, j, yfmv[count], count);
#ifdef  DEBUG_BLOCK_MODE_MV_INFO
			write_block_mode_motion_vector("decoderleft", GOP_counter, count, 1, yfmv[count], info, i);
#endif 
			/////////////// Added by Yuan Liu on 04.21.2016	///////////
//			if(  simul_dec != MUDA ){
				copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
				clear_frame_motion_field_simp(prev_frame_motion_field1, info);
//			}
			///////////////////////////////////////////////////////////
          } 
		  else if (dec_scene_change[i][j + 1] == NO) {

			/////////////// Added by Yuan Liu on 04.21.2016	///////////
//			if( GOP_counter != (get_GOP_num(info) - 1) ){
				copy_frame_motion_field2(prev_frame_motion_field1, prev_frame_motion_field2, info);
				clear_frame_motion_field_simp(prev_frame_motion_field1, info);
//			}
			///////////////////////////////////////////////////////////
			printf("\n\nRIGHT MC, i = %d, j = %d, eff_GOPsz = %d, count+1 = %d\n",i,j,eff_GOPsz[i], count+1);

            mvStat_setFrame(GOP_counter, i, j + 1);
            cur_bytes = mv_decode(yfmv[count + 1], NULL, 1,
                                 meandepth, info, dist, i, starting_pos, 1, count+1, GOP_counter );
			mvBytes += cur_bytes;   starting_pos += cur_bytes; 
			if ( info.layer_mv[i] && !info.layer_exist[i] ) // enhancement layer is discarded
				layer_structure_mv_trim(info, i, j+1, yfmv[count+1], count+1);
#ifdef  DEBUG_BLOCK_MODE_MV_INFO
			write_block_mode_motion_vector("decoderright", GOP_counter, count+1, 2, yfmv[count+1], info, i);
#endif
          }
		  else{   //Added by Yuan Liu on 04.21.2016
			assert(dec_scene_change[i][j] == YES && dec_scene_change[i][j + 1] == YES);
			printf("\n\nISOLATED FRAMES, i = %d, j = %d, eff_GOPsz = %d, count = %d, count+1 = %d\n",i,j,eff_GOPsz[i], count,count+1);
			printf("\n\nISOLATED FRAMES, i = %d, j = %d, eff_GOPsz = %d, count+1 = %d, count = %d\n",i,j,eff_GOPsz[i], count+1,count);
			clear_frame_motion_field_simp(prev_frame_motion_field1, info);
			clear_frame_motion_field_simp(prev_frame_motion_field2, info);
		  }
        }
        count += 2;
      }
    }// j
	//////////////////////
	if(i <= info.tPyrLev - 2 && theo_dec == NO){
		copy_frame_motion_field2(prev_frame_motion_field1, buffer_frame_motion_field1_dec[i], info);
		copy_frame_motion_field2(prev_frame_motion_field2, buffer_frame_motion_field2_dec[i], info);
	}

	//////////////////////
	dist /= 2;
  }// i
 
  // free
  if(GOP_counter == (get_GOP_num(info) - 1) ){
	  for(i = 0;i < G_LEVEL; i++){
		free(buffer_frame_motion_field1_dec[i]);
		free(buffer_frame_motion_field2_dec[i]);
		free(save_buffer_frame_motion_field1_dec[i]);
		free(save_buffer_frame_motion_field2_dec[i]);
	  }
  }

  free(prev_frame_motion_field1);
  free(prev_frame_motion_field2);

  free(prev_frame_motion_field_left1);
  free(prev_frame_motion_field_left2);

/*
  printf("\n tmp lvl 0 - 2\n");
  for(i=0;i<=11;i++)
	  printf("bi_mode_num012[%d] = %d\n",i,bi_mode_num012[i]);

  printf("\n tmp lvl 3\n");
  for(i=0;i<=11;i++)
	  printf("bi_mode_num345[%d] = %d\n",i,bi_mode_num345[i]);

  printf("\n");
  for(i=0;i<27;i++)
	printf("aff_cum_idx[%d] = %d\n",i,aff_cum_idx[i]);
*/
  return mvBytes;
}


/*****************************************************************************/
/*                             child_mv_compare                              */
/*****************************************************************************/
void
child_mv_compare( vector_ptr fmv1, vector_ptr fmv2, videoinfo info )
{
  if( fmv1->child ) {
    child_mv_compare( fmv1->child0, fmv2->child0, info );
    child_mv_compare( fmv1->child1, fmv2->child1, info );
    child_mv_compare( fmv1->child2, fmv2->child2, info );
    child_mv_compare( fmv1->child3, fmv2->child3, info );
  } else {
    if( fmv1->mvx != fmv2->mvx ) {
      printf( "error x\n" );
      exit( 1 );
    }
    if( fmv1->mvy != fmv2->mvy ) {
      printf( "error y\n" );
      exit( 1 );
    }
  }
}



/****************************************************************************/
/*                              write_number()                              */
/****************************************************************************/
void
write_number(  )
{

  if( outbyte < 0 ) {
    printf( "error in write_number() outbyte = %d\n", outbyte );
    exit( 1 );
  }

  if( outbyte <= 0x7fffffff ) {
    write_number_core( outbyte, fpbit );
    //printf("in write_number() outbyte=%d\n", outbyte);
  } else {
    printf( "error in write_number() outbyte=%d\n", outbyte );
    exit( 1 );
  }

}


/****************************************************************************/
/*                                 read_number()                            */
/****************************************************************************/
long int
read_number( int with )
{

  if( with ) {
    inbyte = read_number_core( fpbit );
  } else
    inbyte = 65535;

  //printf("inbyte = %d ..............(arcodemv.c)\n", inbyte);
  return ( inbyte + 4 );        // 4 bytes for the value of inbytes
}


int
write_GOPheader( enum FLAG **scene_change, videoinfo info ) // 就是把场景是否改变写进码流，改变标记为1，未改变标记为零
{
  int i, j, eff_GOPsz, GOPsz, count = 0;
  enum FLAG Level_change;  // flag：YES/NO，用来标记是否是最后一个不满的gop

  GOPsz = info.GOPsz;
  
  if ( info.eff_GOPsz < info.GOPsz ) {
    Level_change = YES;
    eff_GOPsz = info.eff_GOPsz;
  }
  else {
    Level_change = NO;
    eff_GOPsz = info.GOPsz;
  }
  
  encode_init( info ); // 初始化一下编码

  for( i = 0; i < info.tPyrLev; i++ ) {
    for( j = 0; j <= GOPsz; j++ ) {
      if ( j != 0 ) { // scene_change[i][0] is always YES!

        if( scene_change[i][j] == YES ) {
          output_bit( 1 ); count++;
          //          printf("scene_change[%d][%d] = YES\n",i,j);
//		  printf("i = %d, j = %d, 1   ",i,j);
        } else {
          output_bit( 0 ); count++;
          //          printf("scene_change[%d][%d] = NO\n",i,j);
//		  printf("i = %d, j = %d, 0   ",i,j);
        }
        
      }
      else { 
        //    printf("scene_change[%d][%d] = YES --> not transmitted\n",i,j);
      }
    }  
    GOPsz /= 2;  
    
    if ( Level_change == YES )
      eff_GOPsz = (int) ceil ((double)info.eff_GOPsz / (pow ( (float)2, (int)i + 1)));  
    else
      eff_GOPsz = GOPsz;
  }
//  printf("\n\n");

  encode_end( 0, info );
  
  return outbyte;
}


// 标记每一帧是否发生场景改变，返回这个gop标记这些场景改变花费的bit
int
read_GOPheader( enum FLAG **dec_scene_change, videoinfo info )// 
{
  int i, j, bit, GOPsz, GOPheader_bytes, count = 0;
  
  GOPsz = info.GOPsz;

  decode_init(  );
  
  for( i = 0; i < info.tPyrLev; i++ ) { // 统计有多少帧
    for( j = 1; j <= GOPsz; j++ ) 
      count++;
    GOPsz /= 2;  
  } 
  
  GOPsz = info.GOPsz;

  GOPheader_bytes = ( count + 7 ) / 8;
  inbyte = GOPheader_bytes;

  for( i = 0; i < info.tPyrLev; i++ ) { // 逐时域层
    for( j = 0; j <= GOPsz; j++ ) {   // 逐gop帧
      if (j != 0) { // scene_change[i][0] is always YES!
        input_bit( bit );// 
        if( bit == 1 ) {
          dec_scene_change[i][j] = YES;
          //          printf("scene_change[%d][%d] = YES\n",i,j);
        } else {
          dec_scene_change[i][j] = NO;
          //          printf("scene_change[%d][%d] = NO\n",i,j);
        }
      }
      else { 
        dec_scene_change[i][j] = YES;
        //    printf("scene_change[%d][%d] = YES --> not transmitted\n",i,j);
      }
    }
    GOPsz /= 2;  
  }
  
  decode_end(  );

  return GOPheader_bytes;
}
