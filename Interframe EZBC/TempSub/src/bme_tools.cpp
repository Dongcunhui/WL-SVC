#include "bme_tools.h"
#include "mv_ec.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int enc_trans_info(int aff_idx){

	int enc_trans = 0;

	if(aff_idx == 3 || aff_idx == 7 || (aff_idx >= 11 && aff_idx <= 15) )
		enc_trans = 1;
	if(aff_idx == 19 || aff_idx == 23 || (aff_idx >= 27 && aff_idx <= 31) )
		enc_trans = 1;
	if(aff_idx == 35 || aff_idx == 39 || (aff_idx >= 43 && aff_idx <= 47) )
		enc_trans = 1;
	if(aff_idx >= 48 && aff_idx <= 63 )
		enc_trans = 1;

	return enc_trans;
}
////////////////////////////////

inline float median(float a, float b, float c)
{
  return ((a>b) ? ((a>c) ? ((b>c) ? b : c) : a) : 
          ((b>c) ? ((a>c) ? a : c) : b));
}

void get_dmv_from_tree(float *dmvx, float *dmvy, enum FLAG *is_predictor,
                       vector_ptr fmv, int x_dst, int y_dst,  
                       int x, int y, int xblk, int yblk, int hor, int ver, int blk_thresh)
{
  int cx, cy;
  int dx, dy;

  assert(!(x < 0 || x >= hor || y < 0 || y >= ver));
  
  dx = x_dst - x;
  dy = y_dst - y;

  assert(!(dx < 0 || dx >= xblk || dy < 0 || dy >= xblk));

  if (fmv->child  && xblk>blk_thresh)
  {
    if (dx < xblk / 2 && dy < yblk / 2) {
      cx = x;
      cy = y;
      get_dmv_from_tree(dmvx, dmvy, is_predictor, fmv->child0, x_dst, y_dst, 
                        cx, cy, xblk / 2, yblk / 2, hor, ver, blk_thresh);
    } else if (dx >= xblk / 2 && dy < yblk / 2) {
      cx = x + xblk / 2;
      cy = y;
      get_dmv_from_tree(dmvx, dmvy, is_predictor, fmv->child1, x_dst, y_dst,
                        cx, cy, xblk / 2, yblk / 2, hor, ver, blk_thresh);
    } else if (dx < xblk / 2 && dy >= yblk / 2) {
      cx = x;
      cy = y + yblk / 2;
      get_dmv_from_tree(dmvx, dmvy, is_predictor, fmv->child2, x_dst, y_dst,
                        cx, cy, xblk / 2, yblk / 2, hor, ver, blk_thresh);
    } else if (dx >= xblk / 2 && dy >= yblk / 2) {
      cx = x + xblk / 2;
      cy = y + yblk / 2;
      get_dmv_from_tree(dmvx, dmvy, is_predictor, fmv->child3, x_dst, y_dst,
                        cx, cy, xblk / 2, yblk / 2, hor, ver, blk_thresh);
    } else {
      printf("error in get_dmv_from_tree!\n");
      exit(1);
    }
  }
  else
  {
	  if ( !blk_thresh)  // no layered structure for motion vector coding 
	  {
		*dmvx = fmv->dmvx;
		*dmvy = fmv->dmvy;
		*is_predictor = fmv->is_predictor;
	  }else // for layered structure motion vector coding 
	  {
		  if ( ! fmv->child ) // all the information is availabe in this fmv
		  {
			  *dmvx = fmv->dmvx;
			  *dmvy = fmv->dmvy;
		      *is_predictor = fmv->is_predictor;
		  }else  // it's a subsampled (blk_thresh x blk_thresh ) block
		  {
			  if ( fmv->mv_exist )
			  {
				  *dmvx = fmv->dmvx;
			      *dmvy = fmv->dmvy;
		          *is_predictor = YES; 
			  }else
			  {
				  *dmvx = (float)HUGE_VAL;
			      *dmvy = (float)HUGE_VAL;
		          *is_predictor = NO; 
			  }
		  }
	  }
  }
}

void get_pred_from_tree(float *pmvx, float *pmvy, enum FLAG *is_predictor,
                        vector_ptr fmv, int x_dest, int y_dest, int x, int y, 
                        int xblk, int yblk, int hor, int ver, enum BiMode *block_mode,
						int *propagate_iblk, int blk_thresh)
{ //Modified by Yuan Liu on 01.04.2016
  int cx, cy;
  int dx, dy;

  int xblk2,yblk2, it_mvx,it_mvy;

  int subpel = 2;
	  
  float accu = 0.25;
  int addx,addy;
  float int_mvx,int_mvy;

  assert(!(x < 0 || x >= hor || y < 0 || y >= ver));
  
  dx = x_dest - x;
  dy = y_dest - y;

  assert(!(dx < 0 || dx >= xblk || dy < 0 || dy >= xblk));

  if (fmv->child && xblk>blk_thresh)
  {
    if (dx < xblk / 2 && dy < yblk / 2) {
      cx = x;
      cy = y;
      get_pred_from_tree(pmvx, pmvy, is_predictor, fmv->child0, x_dest, y_dest,
                         cx, cy, xblk / 2, yblk / 2, hor, ver, block_mode, propagate_iblk, blk_thresh);
    } else if (dx >= xblk / 2 && dy < yblk / 2) {
      cx = x + xblk / 2;
      cy = y;
      get_pred_from_tree(pmvx, pmvy, is_predictor, fmv->child1, x_dest, y_dest,
                         cx, cy, xblk / 2, yblk / 2, hor, ver, block_mode, propagate_iblk, blk_thresh);
    } else if (dx < xblk / 2 && dy >= yblk / 2) {
      cx = x;
      cy = y + yblk / 2;
      get_pred_from_tree(pmvx, pmvy, is_predictor, fmv->child2, x_dest, y_dest,
                         cx, cy, xblk / 2, yblk / 2, hor, ver, block_mode, propagate_iblk, blk_thresh);
    } else if (dx >= xblk / 2 && dy >= yblk / 2) {
      cx = x + xblk / 2;
      cy = y + yblk / 2;
      get_pred_from_tree(pmvx, pmvy, is_predictor, fmv->child3, x_dest, y_dest,
                         cx, cy, xblk / 2, yblk / 2, hor, ver, block_mode, propagate_iblk, blk_thresh);
    } else {
      printf("error in get_pred_from_tree!\n");
      exit(1);
    }
  }
  else
  {
    //////////////Added on 01.10.2016//////////////
	xblk2 = ( x + xblk <= hor) ? xblk : hor - x;
	yblk2 = ( y + yblk <= ver) ? yblk : ver - y;
	///////////////////////////////////////////////
	if( (dx+1)/xblk2 == 1 )
		addx = 1;
	else
		addx = 0;

	if( (dy+1)/yblk2 == 1 )
		addy = 1;
	else
		addy = 0;

	if (!blk_thresh) // not for layered structure motion vector coding 
	{
		*is_predictor = fmv->is_predictor;
		*block_mode   = fmv->bi_mode;   // by Yongjun Wu: get the block mode 
		*propagate_iblk = fmv->propagate_iblk; // by Yongjun Wu: get the propagation property for the iblock

		if( (fmv->bi_mode >= 0 && fmv->bi_mode <= 6) || fmv->bi_mode == 8 || (fmv->bi_mode == 7 && fmv->aff_mrg == NO) ){
			*pmvx = fmv->mvx;
			*pmvy = fmv->mvy;
		//////////	Added by Yuan Liu	//////////////////
		}else if( (fmv->bi_mode >= 9 && fmv->bi_mode <= 11) || (fmv->bi_mode == 7 && fmv->aff_mrg == YES) ){

			*pmvx = (fmv->aff2_mvx - fmv->aff1_mvx)*((float)dx + addx)/((float)xblk2) + (fmv->aff3_mvx - fmv->aff1_mvx)*((float)dy + addy)/((float)yblk2) + fmv->aff1_mvx;
			*pmvy = (fmv->aff2_mvy - fmv->aff1_mvy)*((float)dx + addx)/((float)xblk2) + (fmv->aff3_mvy - fmv->aff1_mvy)*((float)dy + addy)/((float)yblk2) + fmv->aff1_mvy;
		}
		else
			assert(fmv->is_predictor == NO);

		*pmvx = *pmvx * (1 << subpel);
		*pmvy = *pmvy * (1 << subpel);

		it_mvx = (int)(*pmvx);
		it_mvy = (int)(*pmvy);

		*pmvx = (float)(it_mvx);
		*pmvy = (float)(it_mvy);

		*pmvx = *pmvx / (1 << subpel);
		*pmvy = *pmvy / (1 << subpel);
		//////////////////////////////////////////////////

	}else // for layered structure motion vector coding 
	{
		printf("blk_thresh!!\n");
		assert(0);
		if ( !fmv->child ) // all the information is availabe in this fmv
		{
			*pmvx = fmv->mvx;
			*pmvy = fmv->mvy;
			*is_predictor   = fmv->is_predictor;
			*block_mode     = fmv->bi_mode;   // by Yongjun Wu: get the block mode 
			*propagate_iblk = fmv->propagate_iblk; // by Yongjun Wu: get the propagation property for the iblock
		}else  // it's a subsampled (blk_thresh x blk_thresh ) block
		{
			if ( fmv->mv_exist )
			{
				*pmvx = fmv->sample_mvx; 
				*pmvy = fmv->sample_mvy;
				*is_predictor = fmv->is_predictor;
				*block_mode   = UNDEFINED; 
				*propagate_iblk = -1;
			}else
			{
				*pmvx = (float)HUGE_VAL; 
				*pmvy = (float)HUGE_VAL;
				*is_predictor = NO; 
				*block_mode   = UNDEFINED; 
				*propagate_iblk = -1;
			}
		}
	}

  }
}

void get_dmv(float *dmvx, float *dmvy, enum FLAG *is_predictor,
             vector_ptr fmv, int x_dest, int y_dest, videoinfo info, int t_level, int blk_thresh)
{
  int xnum, ynum, xblk, yblk, hor, ver, X, Y, x, y, pos;

  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];
  hor = info.ywidth;
  ver = info.yheight;

  if (x_dest < 0 || x_dest >= hor || y_dest < 0 || y_dest >= ver)
  {
    *dmvx = 0;
    *dmvy = 0;
    *is_predictor = NO;

    return;
  }
    
  X = x_dest / xblk;
  Y = y_dest / yblk;
  pos = Y * xnum + X;
  x = X * xblk;
  y = Y * yblk;

  get_dmv_from_tree(dmvx, dmvy, is_predictor, &fmv[pos], x_dest, y_dest, x, y,
                    xblk, yblk, hor, ver, blk_thresh);
}

void get_predictor(float *pmvx, float *pmvy, enum FLAG *is_predictor,
                   vector_ptr fmv, int x_dest, int y_dest, videoinfo info, int t_level, 
				   enum BiMode *block_mode, int *propagate_iblk, int blk_thresh)
{
  int xnum, xblk, yblk, hor, ver;
  int X, Y, x, y, pos;

  // initialization
  xnum = info.xnum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];
  hor = info.ywidth;
  ver = info.yheight;

  if (x_dest < 0 || x_dest >= hor || y_dest < 0 || y_dest >= ver)
  {
    *pmvx = 0;
    *pmvy = 0;
    *is_predictor = NO;
	*block_mode   = UNDEFINED;

    return;
  }

  X = x_dest / xblk;
  Y = y_dest / yblk;
  pos = Y * xnum + X;
  x = X * xblk;
  y = Y * yblk;

  get_pred_from_tree(pmvx, pmvy, is_predictor, &fmv[pos], 
                     x_dest, y_dest, x, y, xblk, yblk, hor, ver, block_mode, propagate_iblk, blk_thresh);
}

void rec_clear_predictors(vector_ptr fmv)
{
  fmv->is_predictor = NO; // by Yongjun Wu: clear predictor sign for the whole quad-tree 
  if (fmv->child)
  {
    rec_clear_predictors(fmv->child0);
    rec_clear_predictors(fmv->child1);
    rec_clear_predictors(fmv->child2);
    rec_clear_predictors(fmv->child3);
  }
}

void clear_predictors(vector_ptr fmv, videoinfo info, int t_level)
{
  int X, Y, xnum, ynum, pos;
  
  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  
  for( Y = 0; Y < ynum; Y++ ) {
    for( X = 0; X < xnum; X++ ) {
      pos = Y * xnum + X;
      rec_clear_predictors(&fmv[pos]);
    }
  }
}

// get median predictor for block E, using blocks A, B, C, and D
// block scheme:  CBD
//                AN
//				  E
//Modified by Yuan Liu on 01.18.2016
void get_median_predictor(float *pmvx, float *pmvy, vector_ptr fmv, vector_ptr prev_fmv,
                          vector_ptr prev_fmv2, int x_pos, int y_pos, int xblock, int yblock,
                          videoinfo info, int t_level, int blk_thresh)
{
  int hor, ver, x_dest, y_dest, i, k;
  enum FLAG is_pred_a, is_pred_b, is_pred_c, is_pred_d, is_pred_e;
  enum FLAG is_pred_tmp;  //Added on 02.22.2016
  float mvx_a, mvy_a;
  float mvx_b, mvy_b;
  float mvx_c, mvy_c;
  float mvx_d, mvy_d;
  float mvx_e, mvy_e;

  float mvx_tmp, mvy_tmp;	//Added on 02.22.2016
  int blocks_avail, count;
  enum BiMode block_mode;
  int propagate_iblk=0;

  float med_px[4], med_py[4];

  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }

  // initialization
  hor = info.ywidth;
  ver = info.yheight;
  
  xblock = ( x_pos + xblock <= hor ) ? xblock : hor - x_pos;
  yblock = ( y_pos + yblock <= ver ) ? yblock : ver - y_pos;

  assert(x_pos >= 0 && x_pos < hor && y_pos >= 0 && y_pos < ver);

  // try to get predictors for blocks A, B, and C
  // A
  x_dest = x_pos - 1;
  y_dest = y_pos + yblock - 1;
  if (x_dest < 0 || y_dest >= ver) {
    is_pred_a = NO;
  } else {
    get_predictor(&mvx_a, &mvy_a, &is_pred_a, fmv, x_dest, y_dest, info, t_level, 
		&block_mode, &propagate_iblk, blk_thresh);
  }

  // B
  y_dest = y_pos - 1;
  x_dest = x_pos + xblock - 1;
  if (y_dest < 0 || x_dest >= hor) {
    is_pred_b = NO;
  } else {
    get_predictor(&mvx_b, &mvy_b, &is_pred_b, fmv, x_dest, y_dest, info, t_level, 
		&block_mode, &propagate_iblk, blk_thresh);
  }

  // C
  x_dest = x_pos - 1;
  y_dest = y_pos - 1;
  if (x_dest < 0 || y_dest < 0) {
    is_pred_c = NO;
  } else {
    get_predictor(&mvx_c, &mvy_c, &is_pred_c, fmv, x_dest, y_dest, info, t_level, 
		&block_mode, &propagate_iblk, blk_thresh);
  }

  // blocks A, B, and C available?
  if (is_pred_a == YES && is_pred_b == YES && is_pred_c == YES) {
	pmvx[0] = mvx_a;
	pmvy[0] = mvy_a;
	pmvx[1] = mvx_b;
	pmvy[1] = mvy_b;
	pmvx[2] = mvx_c;
	pmvy[2] = mvy_c;

  } else {
// try to get predictor for block D
// D
	x_dest = x_pos + xblock;
    y_dest = y_pos - 1;
    if (x_dest >= hor || y_dest < 0) {
      is_pred_d = NO;
    } else {
      get_predictor(&mvx_d, &mvy_d, &is_pred_d, fmv, x_dest, y_dest, info, t_level, 
		  &block_mode, &propagate_iblk, blk_thresh);
    }

////////////////////	Added by Yuan Liu	////////////////////
	if(is_pred_a == YES){
		pmvx[0] = mvx_a;
		pmvy[0] = mvy_a;
	}
	if(is_pred_b == YES){
		pmvx[1] = mvx_b;
		pmvy[1] = mvy_b;
	}
	if(is_pred_c == YES){
		pmvx[2] = mvx_c;
		pmvy[2] = mvy_c;
	}
	if(is_pred_d == YES){
		pmvx[3] = mvx_d;
		pmvy[3] = mvy_d;
	}
////////////////////////////////////////////////////////////////
  }

///////////// Modified by by Yuan Liu	//////////////////
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
		pmvx[i] = (float)HUGE_VAL;
		pmvy[i] = (float)HUGE_VAL;
    }
  }
//////////////////////////////////////////////////////////

//////////////  Added on 02.02.2016  /////////////////////////
  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 3);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }
  assert(pmvx[3] == (float)HUGE_VAL && pmvy[3] == (float)HUGE_VAL);
//////////////////////////////////////////////////////////////
  //Set median prediction candidate if available
  if(pmvx[0] != (float)HUGE_VAL && pmvy[0] != (float)HUGE_VAL && pmvx[1] != (float)HUGE_VAL && pmvy[1] != (float)HUGE_VAL
	  && pmvx[2] != (float)HUGE_VAL && pmvy[2] != (float)HUGE_VAL){
	pmvx[3] = median( pmvx[0],pmvx[1],pmvx[2]);
	pmvy[3] = median( pmvy[0],pmvy[1],pmvy[2]);
  }

  ///////////// Modified by by Yuan Liu	//////////////////
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
		pmvx[i] = (float)HUGE_VAL;
		pmvy[i] = (float)HUGE_VAL;
    }
  }
  //////////////////////////////////////////////////////////

  //////////////  Added on 02.21.2016  /////////////////////////
  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 4);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }
//////////////////////////////////////////////////////////////

//  if(pmvx[3]==(float)HUGE_VAL && pmvy[3]==(float)HUGE_VAL && pmvx[2]!=(float)HUGE_VAL && pmvy[2]!=(float)HUGE_VAL){
//	  assert( pmvx[0]!=(float)HUGE_VAL && pmvy[0]!=(float)HUGE_VAL && pmvx[1]!=(float)HUGE_VAL && pmvy[1]!=(float)HUGE_VAL );
  for( count = 0; count <= 3; count ++ ){
    if(pmvx[count] == (float)HUGE_VAL && pmvy[count] == (float)HUGE_VAL){
	  if(prev_fmv != NULL){
		x_dest = x_pos;
		y_dest = y_pos;
		assert(x_dest >= 0 && y_dest >= 0);
		get_predictor(&mvx_tmp, &mvy_tmp, &is_pred_tmp, prev_fmv, x_dest, y_dest, info, t_level, 
		  &block_mode, &propagate_iblk, blk_thresh);
		if(mvx_tmp != (float)HUGE_VAL && mvy_tmp != (float)HUGE_VAL && is_pred_tmp == YES){
			pmvx[count] = mvx_tmp;
			pmvy[count] = mvy_tmp;
		}
	  }
	  break;
	}
  }

//////////////////////////////////////////////////////////////
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
      pmvx[i] = (float)HUGE_VAL;
      pmvy[i] = (float)HUGE_VAL;
    }
  }
///////////////////////////////////////////////////////

  //////////////  Added on 02.27.2016  /////////////////////////
  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 4);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }

//block E
  for( count = 0; count <= 3; count ++ ){
    if(pmvx[count] == (float)HUGE_VAL && pmvy[count] == (float)HUGE_VAL){
		x_dest = x_pos - 1;
		y_dest = y_pos + yblock;
		if (x_dest < 0 || y_dest >= ver) {
		  is_pred_e = NO;
		} else {
		  get_predictor(&mvx_e, &mvy_e, &is_pred_e, fmv, x_dest, y_dest, info, t_level, 
			  &block_mode, &propagate_iblk, blk_thresh);
		  if(mvx_e != (float)HUGE_VAL && mvy_e != (float)HUGE_VAL && is_pred_e == YES){
			pmvx[count] = mvx_e;
			pmvy[count] = mvy_e;
		  }
		}
		break;
	}
  }

//////////////////////////////////////////////////////////////
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
      pmvx[i] = (float)HUGE_VAL;
      pmvy[i] = (float)HUGE_VAL;
    }
  }
///////////////////////////////////////////////////////

  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 4);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }

}

float get_bit_cost(float lambda, float mvx, float mvy, float pmvx, float pmvy,
                   int ctx_x, int ctx_y, int subpel)
{
  int dmvx, dmvy;
  float float_dmvx, float_dmvy;

  if (lambda > 0.) {

	dmvx = (int) ((1 << subpel) * (mvx - pmvx));
    dmvy = (int) ((1 << subpel) * (mvy - pmvy));

    float_dmvx = ((1 << subpel) * (mvx - pmvx));
    float_dmvy = ((1 << subpel) * (mvy - pmvy));

	if(float_dmvx != (int)float_dmvx || float_dmvy != (int)float_dmvy){
		return (float)HUGE_VAL;
	}else    
		return lambda * (ec_get_expected_length(dmvx, ctx_x) +
						 ec_get_expected_length(dmvy, ctx_y));
  } else {
    return 0.;
  }
}
