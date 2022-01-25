#include "hvsbm_fullres.h"
#include "mode_decision.h"
#include "mv_ec.h"
#include "basic.h"
#include "bmeN.h"
#include "bme_tools.h"
#include "util_filtering.h"
#include "memoryN.h"
#include "mvcodingN.h"
#include "miscN.h"
#define EXTERN extern

#define ANCHOR_LAMBDA 10



#include "coderN.h"
#include <assert.h>

float ratio0 = -1, ratio1 = -1;
float rate0 = -1, rate1 = -1, rate2 = -1;
float new_lambda[LAMBDA_ADPT_LEVEL];
int   revise0 = 0, revise1 = 0, revise2 = 0;


//Added on 09.10.2017
void lambda_revise(float *lambda, int t_level, float add_num){
	*lambda += add_num;

	if( *lambda <= 1 )
		*lambda -= add_num;
}



void print_info(vector_ptr fmv1, vector_ptr fmv2, int xblock, int yblock, int x, int y){

	if(fmv1->child){
		if(fmv2){
			print_info(fmv1->child0,fmv2->child0,xblock/2,yblock/2,x,y);
			print_info(fmv1->child1,fmv2->child1,xblock/2,yblock/2,x+xblock /2,y);
			print_info(fmv1->child2,fmv2->child2,xblock/2,yblock/2,x,y+yblock /2);
			print_info(fmv1->child3,fmv2->child3,xblock/2,yblock/2,x+xblock /2,y+yblock /2);
		}
		else{
			print_info(fmv1->child0,NULL,xblock/2,yblock/2,x,y);
			print_info(fmv1->child1,NULL,xblock/2,yblock/2,x+xblock /2,y);
			print_info(fmv1->child2,NULL,xblock/2,yblock/2,x,y+yblock /2);
			print_info(fmv1->child3,NULL,xblock/2,yblock/2,x+xblock /2,y+yblock /2);
		}
	}
	else{
		if(fmv1->bi_mode == DIRECTIONAL_IBLOCK)printf("amigo! ");
		printf("x = %d, y = %d, blkx = %d, blky = %d, fmv1_x = %1.2f, fmv1_y = %1.2f\n",x, y, xblock, yblock, fmv1->mvx,fmv1->mvy);
		if(fmv2!=NULL){
		if(fmv2->bi_mode == DIRECTIONAL_IBLOCK)printf("amigo! ");
			printf("x = %d, y = %d, blkx = %d, blky = %d, fmv2_x = %1.2f, fmv2_y = %1.2f\n",x, y, xblock, yblock, fmv2->mvx,fmv2->mvy);
		}
	}

}

void create_children(vector_ptr fmv)
{
  int i;

  assert(fmv->child == 0);

  fmv->child=1;

  fmv->child0 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child0" );
  fmv->child1 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child1" );
  fmv->child2 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child2" );
  fmv->child3 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child3" );

  fmv->child0->parent = fmv;
  fmv->child1->parent = fmv;
  fmv->child2->parent = fmv;
  fmv->child3->parent = fmv;

  fmv->child0->child = 0;
  fmv->child0->bi_mode = UNDEFINED;
  fmv->child0->lifting_mode = IGNORED;
  fmv->child0->is_predictor = NO;
  fmv->child0->propagate_iblk = 0; 

  fmv->child1->child = 0;
  fmv->child1->bi_mode = UNDEFINED;
  fmv->child1->lifting_mode = IGNORED;
  fmv->child1->is_predictor = NO;
  fmv->child1->propagate_iblk = 0; 
  
  fmv->child2->child = 0;
  fmv->child2->bi_mode = UNDEFINED;
  fmv->child2->lifting_mode = IGNORED;
  fmv->child2->is_predictor = NO;
  fmv->child2->propagate_iblk = 0; 

  fmv->child3->child = 0;
  fmv->child3->bi_mode = UNDEFINED;
  fmv->child3->lifting_mode = IGNORED;
  fmv->child3->is_predictor = NO;
  fmv->child3->propagate_iblk = 0; 

  for (i = 0; i < NUMBER_OF_BI_MODES; i++) {
    fmv->child0->mode_info[i] = invalid_mode_info;
    fmv->child1->mode_info[i] = invalid_mode_info;
    fmv->child2->mode_info[i] = invalid_mode_info;
    fmv->child3->mode_info[i] = invalid_mode_info;
  }
}

void kill_children(vector_ptr fmv)
{
  if(fmv->child) {
    kill_children(fmv->child0);
    kill_children(fmv->child1);
    kill_children(fmv->child2);
    kill_children(fmv->child3);

    free(fmv->child0);
    free(fmv->child1);
    free(fmv->child2);
    free(fmv->child3);

    fmv->child = 0;
  }
} 

// 向fmv1和fmv2写东西，来表示这个块的运动信息
void rec_hvsbm_fullres(vector_ptr fmv1_array, vector_ptr fmv2_array,
                       vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3_array, vector_ptr fmv4_array,
                       float *fr_cur, float *fr_ref1, float *fr_ref2,
                       float *upframe1, float *upframe2,
                       float *pmv1x, float *pmv1y, float *pmv2x, float *pmv2y,
                       int x, int y, int xblock, int yblock, 
                       int maxx, int maxy, int hor, int ver, 
                       int block_level, int t_level, videoinfo info, int dec)
{
  int i;
  float split_sad_cost, split_bit_cost, split_total_cost;
  int cx, cy,xblk,yblk;
  int ctx1x, ctx1y, ctx2x, ctx2y;
  float mvpred1x[4], mvpred1y[4], mvpred2x[4], mvpred2y[4];
#ifndef MEDIAN_PREDICTION
  int xblk, yblk;  
#endif

  if(block_level == 1 && (t_level == 0 || t_level == 1 || t_level == 2) && (x == 0 && y == 0) )
	  printf("lambda = %f\n",info.lambda[t_level]);

/////////////////// Added by Yuan Liu on 01.18.2016  /////////////////
  for(i=0;i<=3;i++){
	mvpred1x[i] = (float)HUGE_VAL;
	mvpred1y[i] = (float)HUGE_VAL;
	mvpred2x[i] = (float)HUGE_VAL;
	mvpred2y[i] = (float)HUGE_VAL;
  }
//////////////////////////////////////////////////////////////////////

  assert((fmv2_array == NULL && fmv2 == NULL) ||
         (fmv2_array != NULL && fmv2 != NULL));

  // initialization
  fmv1->is_predictor = NO;
  if (fmv2 != NULL) fmv2->is_predictor = NO;

//  printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",x,y,xblock,yblock);

  // block validity check
  if (x >= hor || y >= ver) {
    //  fmv1->mad = 0.;
    fmv1->sad_cost = 0.;
    fmv1->bit_cost = 0.;
    fmv1->total_cost = 0.;

    if (fmv2 != NULL) {
      //  fmv2->mad = 0.;
      fmv2->sad_cost = 0.;
      fmv2->bit_cost = 0.;
      fmv2->total_cost = 0.;
    }

    return;
  }

  // update contexts and predictors
#ifdef ECSIM_USE_CONTEXTS
  ec_get_contexts(&ctx1x, &ctx1y, fmv1_array, x, y, info, t_level, 0);
  if (fmv2_array != NULL) {
    ec_get_contexts(&ctx2x, &ctx2y, fmv2_array, x, y, info, t_level, 0);
  }
#else
  ctx1x = ctx1y = 0;
  ctx2x = ctx2y = 0;
#endif

#ifdef MEDIAN_PREDICTION
  get_median_predictor(mvpred1x, mvpred1y, fmv1_array, fmv3_array, fmv4_array, 
                       x, y, xblock, yblock, info, t_level, 0);
  if (fmv2 != NULL) {
    get_median_predictor(mvpred2x, mvpred2y, fmv2_array, fmv4_array, NULL,
                         x, y, xblock, yblock, info, t_level, 0);
  }
#else
  xblk = ( x + xblock <= hor ) ? xblock : hor - x;
  yblk = ( y + yblock <= ver ) ? yblock : ver - y;

  // store zig-zag predictors if pointing inside reference frame
  if (x - (*pmv1x) < 0 || x - (*pmv1x) + xblk > hor  || 
      y - (*pmv1y) < 0 || y - (*pmv1y) + yblk > ver) {
    mvpred1x = 0.;
    mvpred1y = 0.;
  } else {
    mvpred1x = *pmv1x; 
    mvpred1y = *pmv1y;
  }
  if (x - (*pmv2x) < 0 || x - (*pmv2x) + xblk > hor  || 
      y - (*pmv2y) < 0 || y - (*pmv2y) + yblk > ver) {
    mvpred2x = 0.;
    mvpred2y = 0.;
  } else {
    mvpred2x = *pmv2x;
    mvpred2y = *pmv2y;
  }
#endif

  // motion estimation and mode decision
  // note that block size decision starts with smallest block size
  find_best_mode(fmv1_array, fmv2_array, fmv1, fmv2, fmv3_array, fmv4_array, fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
                 x, y, xblock, yblock, maxx, maxy, hor, ver, t_level, info, 
                 ctx1x, ctx1y, ctx2x, ctx2y, 
                 mvpred1x, mvpred1y, mvpred2x, mvpred2y, dec);

#ifdef ECSIM_USE_CONTEXTS
  // update coding engine
  // note that in mvcoding.c, the MCFs are encoded separately!!!
  // solution: increase the number of contexts and use different contexts
  //           for fmv1 and fmv2, which may be problematic with mv_ec!
  assert(EC_SIM_TYPE != AR_NARY);
  if (fmv1->lifting_mode != IGNORED) {
	  fmv1->dmvx = fmv1->mvx - (fmv1->med_idx >= 0)? mvpred1x[fmv1->med_idx] : 0.0;
	  fmv1->dmvy = fmv1->mvy - (fmv1->med_idx >= 0)? mvpred1y[fmv1->med_idx] : 0.0;
      ec_update_model((int)((1 << info.subpel[t_level]) * fmv1->dmvx), ctx1x);
      ec_update_model((int)((1 << info.subpel[t_level]) * fmv1->dmvy), ctx1y);
  }
  if (fmv2 != NULL && fmv2->lifting_mode != IGNORED) {
      fmv2->dmvx = fmv2->mvx - (fmv2->med_idx >= 0)? mvpred2x[fmv2->med_idx] : 0.0;
      fmv2->dmvy = fmv2->mvy - (fmv2->med_idx >= 0)? mvpred2y[fmv2->med_idx] : 0.0;
      ec_update_model((int)((1 << info.subpel[t_level]) * fmv2->dmvx), ctx2x);
      ec_update_model((int)((1 << info.subpel[t_level]) * fmv2->dmvy), ctx2y);
  }
#endif

  // minimum block size (e.g. 4x4) still not reached?
  if (block_level < info.level[t_level]) {
    // create children
    create_children(fmv1);
    if (fmv2 != NULL) create_children(fmv2);

    // freeze coder state
    ec_freeze();

    // find block sizes, vectors and modes for children
    cx = x;
    cy = y;
    rec_hvsbm_fullres(fmv1_array, fmv2_array,
                      fmv1->child0, fmv2 ? fmv2->child0 : NULL, fmv3_array, fmv4_array,
                      fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
                      pmv1x, pmv1y, pmv2x, pmv2y,
                      cx, cy, xblock / 2, yblock / 2, maxx, maxy, 
                      hor, ver, block_level + 1, t_level, info, dec);
    cx = x + xblock / 2;
    cy = y;
    rec_hvsbm_fullres(fmv1_array, fmv2_array,
                      fmv1->child1, fmv2 ? fmv2->child1 : NULL, fmv3_array, fmv4_array,
                      fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
                      pmv1x, pmv1y, pmv2x, pmv2y,
                      cx, cy, xblock / 2, yblock / 2, maxx, maxy, 
                      hor, ver, block_level + 1, t_level, info, dec);
    cx = x;
    cy = y + yblock / 2;
    rec_hvsbm_fullres(fmv1_array, fmv2_array,
                      fmv1->child2, fmv2 ? fmv2->child2 : NULL, fmv3_array, fmv4_array,
                      fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
                      pmv1x, pmv1y, pmv2x, pmv2y,
                      cx, cy, xblock / 2, yblock / 2, maxx, maxy, 
                      hor, ver, block_level + 1, t_level, info, dec);
    cx = x + xblock / 2;
    cy = y + yblock / 2;
    rec_hvsbm_fullres(fmv1_array, fmv2_array,
                      fmv1->child3, fmv2 ? fmv2->child3 : NULL, fmv3_array, fmv4_array,
                      fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
                      pmv1x, pmv1y, pmv2x, pmv2y,
                      cx, cy, xblock / 2, yblock / 2, maxx, maxy, 
                      hor, ver, block_level + 1, t_level, info, dec);

    // compute splitting costs
    // fmv1 costs and fmv2 costs are synonym
    split_sad_cost = fmv1->child0->sad_cost + fmv1->child1->sad_cost + 
                     fmv1->child2->sad_cost + fmv1->child3->sad_cost;
    split_bit_cost = fmv1->child0->bit_cost + fmv1->child1->bit_cost +
                     fmv1->child2->bit_cost + fmv1->child3->bit_cost;
	// after directional IBLOCK is added in this assertion is not valid any more 
    // assert(split_bit_cost / info.lambda[t_level] == float(int(split_bit_cost / info.lambda[t_level])));

	if(fmv1->bi_mode >= 9 && fmv1->bi_mode <= 11){
//		printf("\nx = %d, y = %d, xblk = %d, yblk = %d",x,y,xblock,yblock);
//		printf("\nsplit_sad_cost = %f, split_bit_cost = %f\n",split_sad_cost,split_bit_cost);
//		printf("fmv1->sad_cost = %f, fmv1->bit_cost = %f\n",fmv1->sad_cost,fmv1->bit_cost);
	}

    split_total_cost = split_sad_cost + split_bit_cost;

    // are you sure you like your children?
    if (split_total_cost > fmv1->total_cost) {
      // if not -> kill them ;-)
      kill_children(fmv1);
      if (fmv2 != NULL) kill_children(fmv2);

      // get previous states (without splitting)
      ec_unfreeze(0);

#ifndef MEDIAN_PREDICTION
      // get new zig-zag predictor
      if (fmv1->is_predictor == YES) {
        *pmv1x = fmv1->mvx;
        *pmv1y = fmv1->mvy;
      } else {
        *pmv1x = mvpred1x;
        *pmv1y = mvpred1y;
      }
      if (fmv2 != NULL && fmv2->is_predictor == YES) {
        *pmv2x = fmv2->mvx;
        *pmv2y = fmv2->mvy;
      } else {
        *pmv2x = mvpred2x;
        *pmv2y = mvpred2y;
      }
#endif
    } else {
      // you like your children
      // current node will not be coded, update cost
      fmv1->is_predictor = NO;
      fmv1->sad_cost = split_sad_cost;
      fmv1->bit_cost = split_bit_cost;
	  // after directional IBLOCK is added in this assertion is not valid any more 
      // assert(fmv1->bit_cost / info.lambda[t_level] == float(int(fmv1->bit_cost / info.lambda[t_level])));
      if (fmv2 != NULL) {
        fmv2->is_predictor = NO;
        fmv2->sad_cost = split_sad_cost;
        fmv2->bit_cost = split_bit_cost;
      }

      // remove state from stack but keep current states (with splitting)
      ec_unfreeze(1);
    }
      
    // compute real cost for current block (consider splitting decision bit)
    fmv1->bit_cost += 1.0f * info.lambda[t_level];
	// after directional IBLOCK is added in this assertion is not valid any more 
    // assert(fmv1->bit_cost / info.lambda[t_level] == float(int(fmv1->bit_cost / info.lambda[t_level])));
    fmv1->total_cost = fmv1->sad_cost + fmv1->bit_cost;
    if (fmv2 != NULL) {
      fmv2->bit_cost = fmv1->bit_cost;
      fmv2->total_cost = fmv1->total_cost;
    }

  }
#ifndef MEDIAN_PREDICTION
 else {
    // leaf node: get zig-zag predictor
    if (fmv1->is_predictor == YES) {
      *pmv1x = fmv1->mvx;
      *pmv1y = fmv1->mvy;
    }
    if (fmv2 != NULL && fmv2->is_predictor == YES) {
      *pmv2x = fmv2->mvx;
      *pmv2y = fmv2->mvy;
    }
  }
#else
  assert((*pmv1x) == 0. && (*pmv1y) == 0. && (*pmv2x) == 0. && (*pmv2y) == 0.);
#endif
}

/*
 * hvsbm_fullres
 * hiearchical variable size block matching on full resolution  挑选出最好的mv，此函数内部会进行rd计算
 */
float hvsbm_fullres(vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3, vector_ptr fmv4,
                    float *fr_cur, float *fr_ref1, float *fr_ref2, 
                    float *upframe1, float *upframe2,
                    int t_level, int dist, videoinfo info, int dec)
{
  int hor, ver, xnum, ynum, xblk, yblk, maxx, maxy;
  int x, y, X, Y, pos;
  float pmv1x, pmv1y, pmv2x, pmv2y;
  float bit_sum, sad_sum, total_sum;
  float lambda_cost[4], best_cost;
  int best_add_count;

  float new_ratio, new_rate;

  FILE *p_r, *p_sad, *p_lambda, *p_ratio;
  int i,j;

  assert(fmv1 && fr_cur && fr_ref1);
  assert((fmv2 && fr_ref2) || (!fmv2 && !fr_ref2));

  // 全是注释
  if(fmv2 != NULL && t_level == 0){
//	  p_r = fopen("rate_cost0.txt","at");
//	  p_sad = fopen("sad_cost0.txt","at");

//	  p_ratio = fopen("RD_ratio0.txt","at");
//	  p_lambda = fopen("lambda_value0.txt","at");
/*
	  if(rate0 != -1){
		assert(rate0 > 0 && new_lambda[t_level] > 0 );
		info.lambda[t_level] = new_lambda[t_level];
	  }
*/
  }
  else if(fmv2 != NULL && t_level == 1)
  {
//	  p_r = fopen("rate_cost1.txt","at");
//	  p_sad = fopen("sad_cost1.txt","at");

//	  p_ratio = fopen("RD_ratio1.txt","at");
//	  p_lambda = fopen("lambda_value1.txt","at");
/*
	  if(rate1 != -1){
		assert(rate1 > 0 && new_lambda[t_level] > 0 );
		info.lambda[t_level] = new_lambda[t_level];
	  }
*/
  }
  else if(fmv2 != NULL && t_level == 2){
//	  p_r = fopen("rate_cost2.txt","at");
//	  p_sad = fopen("sad_cost2.txt","at");
/*
	  if(rate2 != -1){
		assert(rate2 > 0 && new_lambda[t_level] > 0 );
		info.lambda[t_level] = new_lambda[t_level];
	  }
*/
  }
  else if(fmv2 != NULL && t_level == 3)
  {
//	  p_r = fopen("rate_cost3.txt","at");
//	  p_sad = fopen("sad_cost3.txt","at");
  }

  // initialization
  hor = info.ywidth;
  ver = info.yheight;
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];
  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];

  maxx = maxy = get_searchrange(dist, info);
  printf("t_level = %d: block sizes = %dx%d - %dx%d, accuracy = %d, "
         "search range = %d, lambda = %f, xnum = %d, ynum = %d\n", t_level, xblk, yblk, 
         xblk >> (info.level[t_level] - 1), yblk >> (info.level[t_level] - 1), 
         info.subpel[t_level], maxx, info.lambda[t_level], xnum, ynum);

  // initialize entropy coding
#ifdef ECSIM_USE_CONTEXTS
  ec_enc_init(info, ECSIM_TYPE, 0, 3, 0);
#else
  ec_enc_init(info, ECSIM_TYPE, 0, 1, 0);
#endif
  
///////////////////////////////////////////////////
/*
    if(t_level == 0 && fmv2 != NULL){
		for(i = 0; i <= 3; i ++){
			lambda_revise(&info.lambda[t_level],t_level,1);

		  // perform ME + MD
		  bit_sum = 0.;
		  sad_sum = 0.;
		  total_sum = 0.;
		  for(y = 0, Y = 0; Y < ynum; y += yblk, Y++) {
			// initialize zig-zag predictors per row
			pmv1x = 0.;
			pmv1y = 0.;
			pmv2x = 0.;
			pmv2y = 0.;

			for(x = 0, X = 0; X < xnum; x += xblk, X++) {
			  pos = Y * xnum + X;
			  rec_hvsbm_fullres(fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, fmv3, fmv4,
								fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
								&pmv1x, &pmv1y, &pmv2x, &pmv2y,
								x, y, xblk, yblk, maxx, maxy, 
								hor, ver, 1, t_level, info, dec);
			  assert(fmv2 == NULL || (fmv1[pos].bit_cost == fmv2[pos].bit_cost &&
									  fmv1[pos].total_cost == fmv2[pos].total_cost));
			  bit_sum += fmv1[pos].bit_cost;
			  sad_sum += fmv1[pos].sad_cost;
			  total_sum += fmv1[pos].total_cost;
			}
		  }

		  lambda_cost[i] = sad_sum + bit_sum/info.lambda[t_level]*ANCHOR_LAMBDA ;
		  printf("total_sum = %f\n",lambda_cost[i]);

		  free_vector(fmv1, info);
		  free_vector(fmv2, info);
		}//for i

		best_cost = (float)HUGE_VAL;

		for(i = 0; i <= 3; i ++){
			if( lambda_cost[i] < best_cost ){
				best_cost = lambda_cost[i];
				best_add_count = i;
			}
		}

		lambda_revise(&info.lambda[t_level], t_level, best_add_count-3);
		printf("best_lambda = %f\n", info.lambda[t_level]);
	}
///////////////////////////////////////////////////
*/

  // perform ME + MD
  bit_sum = 0.;
  sad_sum = 0.;
  total_sum = 0.;
  for(y = 0, Y = 0; Y < ynum; y += yblk, Y++) { 
    // initialize zig-zag predictors per row
    pmv1x = 0.;
    pmv1y = 0.;
    pmv2x = 0.;
    pmv2y = 0.;

    for(x = 0, X = 0; X < xnum; x += xblk, X++) {  // 按块进行
      pos = Y * xnum + X;
      rec_hvsbm_fullres(fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, fmv3, fmv4,
                        fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
                        &pmv1x, &pmv1y, &pmv2x, &pmv2y,
                        x, y, xblk, yblk, maxx, maxy, 
                        hor, ver, 1, t_level, info, dec);
      assert(fmv2 == NULL || (fmv1[pos].bit_cost == fmv2[pos].bit_cost &&
                              fmv1[pos].total_cost == fmv2[pos].total_cost));
      bit_sum += fmv1[pos].bit_cost;
	  sad_sum += fmv1[pos].sad_cost;
      total_sum += fmv1[pos].total_cost;
    }
  }
  //  printf("total mv bits expected = %f\n", bit_sum / info.lambda[t_level]);


/*
  new_ratio = bit_sum / sad_sum;
  new_rate = bit_sum / info.lambda[t_level];

  if(fmv2 != NULL){
	  if(t_level == 0){
		  if(rate0 == -1){
			rate0 = bit_sum / info.lambda[t_level];
		  }else{
			  if( (new_rate/rate0) < THRES_RATIO ){
				lambda_revise(&info.lambda[t_level],t_level,2);

				free_vector(fmv1, info);
			    free_vector(fmv2, info);

				  bit_sum = 0.;
				  sad_sum = 0.;
				  total_sum = 0.;
				  for(y = 0, Y = 0; Y < ynum; y += yblk, Y++) {
					// initialize zig-zag predictors per row
					pmv1x = 0.;
					pmv1y = 0.;
					pmv2x = 0.;
					pmv2y = 0.;

					for(x = 0, X = 0; X < xnum; x += xblk, X++) {
					  pos = Y * xnum + X;
					  rec_hvsbm_fullres(fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, fmv3, fmv4,
										fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
										&pmv1x, &pmv1y, &pmv2x, &pmv2y,
										x, y, xblk, yblk, maxx, maxy, 
										hor, ver, 1, t_level, info, dec);
					  assert(fmv2 == NULL || (fmv1[pos].bit_cost == fmv2[pos].bit_cost &&
											  fmv1[pos].total_cost == fmv2[pos].total_cost));
					  bit_sum += fmv1[pos].bit_cost;
					  sad_sum += fmv1[pos].sad_cost;
					  total_sum += fmv1[pos].total_cost;
					}
				  }
			  }

		  }//else
		  new_lambda[t_level] = info.lambda[t_level];

	  }else if(t_level == 1){
		  if(rate1 == -1){
			rate1 = bit_sum / info.lambda[t_level];
		  }else{
			if( (new_rate/rate1) < THRES_RATIO ){
				lambda_revise(&info.lambda[t_level],t_level,3);

				free_vector(fmv1, info);
			    free_vector(fmv2, info);

				bit_sum = 0.;
				sad_sum = 0.;
				total_sum = 0.;
				for(y = 0, Y = 0; Y < ynum; y += yblk, Y++) {
					// initialize zig-zag predictors per row
					pmv1x = 0.;
					pmv1y = 0.;
					pmv2x = 0.;
					pmv2y = 0.;

					for(x = 0, X = 0; X < xnum; x += xblk, X++) {
					  pos = Y * xnum + X;
					  rec_hvsbm_fullres(fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, fmv3, fmv4,
										fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
										&pmv1x, &pmv1y, &pmv2x, &pmv2y,
										x, y, xblk, yblk, maxx, maxy, 
										hor, ver, 1, t_level, info, dec);
					  assert(fmv2 == NULL || (fmv1[pos].bit_cost == fmv2[pos].bit_cost &&
											  fmv1[pos].total_cost == fmv2[pos].total_cost));
					  bit_sum += fmv1[pos].bit_cost;
					  sad_sum += fmv1[pos].sad_cost;
					  total_sum += fmv1[pos].total_cost;
					}
				}
			}
		  }//else
		  new_lambda[t_level] = info.lambda[t_level];

	  }else if(t_level == 2){
		  if(rate2 == -1){
			rate2 = bit_sum / info.lambda[t_level];
		  }else{
			if( (new_rate/rate2) < THRES_RATIO ){
				lambda_revise(&info.lambda[t_level],t_level,4);

				free_vector(fmv1, info);
			    free_vector(fmv2, info);

				bit_sum = 0.;
				sad_sum = 0.;
				total_sum = 0.;
				for(y = 0, Y = 0; Y < ynum; y += yblk, Y++) {
					// initialize zig-zag predictors per row
					pmv1x = 0.;
					pmv1y = 0.;
					pmv2x = 0.;
					pmv2y = 0.;

					for(x = 0, X = 0; X < xnum; x += xblk, X++) {
					  pos = Y * xnum + X;
					  rec_hvsbm_fullres(fmv1, fmv2, &fmv1[pos], fmv2 ? (&fmv2[pos]) : NULL, fmv3, fmv4,
										fr_cur, fr_ref1, fr_ref2, upframe1, upframe2,
										&pmv1x, &pmv1y, &pmv2x, &pmv2y,
										x, y, xblk, yblk, maxx, maxy, 
										hor, ver, 1, t_level, info, dec);
					  assert(fmv2 == NULL || (fmv1[pos].bit_cost == fmv2[pos].bit_cost &&
											  fmv1[pos].total_cost == fmv2[pos].total_cost));
					  bit_sum += fmv1[pos].bit_cost;
					  sad_sum += fmv1[pos].sad_cost;
					  total_sum += fmv1[pos].total_cost;
					}
				}
			}
		  }
		  new_lambda[t_level] = info.lambda[t_level];
	  }
  }

*/

  // terminate entropy coding
  ec_enc_end1();
  ec_enc_end2();

/*  for(x=0 ; x <= info.tPyrLev - 1 ; x++){
	printf("layer mv = %d\n",info.layer_mv[x]);
	printf("AGP level = %d\n",info.AGP_level[x]);
  }*/

  if(fmv2 != NULL){

	if(t_level <= 2){
//		fclose(p_r);
//		fclose(p_sad);
	}
  }

  return total_sum;
}



void rdme(vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3, vector_ptr fmv4, 
          YUVimage *fr_cur, YUVimage *fr_ref1, YUVimage *fr_ref2, 
          enum FLAG *sc1, enum FLAG *sc2,
          int t_level, int dist, videoinfo info, float *upframe1, float *upframe2, float *upsamp_x )
{
  int hor, ver, xnum, ynum;
  float cost, min_cost;
  vector_ptr tmp_fmv1, tmp_fmv2;

  int dec;

  hor = info.ywidth;// 图像宽度
  ver = info.yheight; // 图像高度
  xnum = info.xnum[t_level]; // x方向块的数量
  ynum = info.ynum[t_level];// y 方向块的数量

  assert(fmv1 && fr_cur && fr_ref1 && sc1 && (*sc1 == NO));
  assert((fmv2 && fr_ref2 && sc2 && (*sc2 == NO)) || 
         (!fmv2 && !fr_ref2 && !sc2));
  
#ifdef PRE_INTERPOL
  // generate interpolated reference frames
  Interpolate_frame2(fr_ref1->Y, hor, ver, &upframe1, upsamp_x,
                     MY_MIN(PRE_INTERPOL, info.subpel[t_level]));
#else
  upframe1 = NULL;
#endif

  ////////////////////////
  if(fmv2 != NULL)
	  dec = 1;
  else
	  dec = 0;
  ////////////////////////

  tmp_fmv1 = (vector*)getarray(info.maxMBnum, sizeof(vector), "tmp_fmv1");
  alloc_vector(tmp_fmv1, info);

  // generate intra satd
  temporal_filter();
  min_cost = RDMD_INTRA_FACTOR * getFrameSATD(fr_cur->Y, hor, ver);
  assert(*sc1 == NO);
  *sc1 = YES;
  free_vector(fmv1, info);
  if (fmv2 != NULL) {
    assert(*sc2 == NO);
    *sc2 = YES;
    free_vector(fmv2, info);
  }

  // perform left search
  free_vector(tmp_fmv1, info);
  cost = hvsbm_fullres(tmp_fmv1, NULL, fmv3, NULL, fr_cur->Y, fr_ref1->Y, NULL, 
                       upframe1, NULL, t_level, dist, info, dec);
  //  printf("left pred cost = %f\n", cost);
  if (cost < min_cost) {
    min_cost = cost;
    free_vector(fmv1, info);
    mv_copy(tmp_fmv1, fmv1, info);
    *sc1 = NO;
    if (fmv2 != NULL) {
      free_vector(fmv2, info);
      *sc2 = YES;
    }
  }



  // bi-directional case?
  if (fmv2 != NULL
//	  && (simul_skip == NO || skip_frame == NO ) 
	  ) {

#ifdef PRE_INTERPOL
	  Interpolate_frame2(fr_ref2->Y, hor, ver, &upframe2, upsamp_x,
                       MY_MIN(PRE_INTERPOL, info.subpel[t_level]));
#else
    upframe2 = NULL;
#endif

    tmp_fmv2 = (vector*)getarray(info.maxMBnum, sizeof(vector), "tmp_fmv2");
    alloc_vector(tmp_fmv2, info);

    // right search
    free_vector(tmp_fmv2, info);
    cost = hvsbm_fullres(tmp_fmv2, NULL, fmv4, NULL, fr_cur->Y, fr_ref2->Y, NULL, 
                         upframe2, NULL, t_level, dist, info, dec);
    //    printf("right pred cost = %f\n", cost);
    if (cost < min_cost) {
      min_cost = cost;
      free_vector(fmv1, info);
      free_vector(fmv2, info);
      mv_copy(tmp_fmv2, fmv2, info);
      *sc1 = YES;
      *sc2 = NO;
    }
    
    // bi-directional search
    free_vector(tmp_fmv1, info);
    free_vector(tmp_fmv2, info);
    cost = hvsbm_fullres(tmp_fmv1, tmp_fmv2, fmv3, fmv4, fr_cur->Y, fr_ref1->Y, fr_ref2->Y,
                         upframe1, upframe2, t_level, dist, info, 0);
    //    printf("bi pred cost = %f\n", cost);
    if (cost < min_cost) {
      min_cost = cost;
      free_vector(fmv1, info);
      free_vector(fmv2, info);
      mv_copy(tmp_fmv1, fmv1, info);
      mv_copy(tmp_fmv2, fmv2, info);
      *sc1 = NO;
      *sc2 = NO;
    }

    // tidy up
    free_vector(tmp_fmv2, info);
    free(tmp_fmv2);
#ifdef PRE_INTERPOL
//    free(upframe2);
#endif
  }
  free_vector(tmp_fmv1, info);
  free(tmp_fmv1);
  
#ifdef PRE_INTERPOL
//  free(upframe1);
#endif

  printf("=> sc1 = %d, sc2 = %d\n", *sc1, sc2 ? *sc2 : 0);
}

int get_searchrange(int dist, videoinfo info)
{
  int sr;

  assert(info.searchrange >= 0);

  sr = dist * info.searchrange;

  return (sr < info.maxsearchrange) ? sr : info.maxsearchrange;
}
