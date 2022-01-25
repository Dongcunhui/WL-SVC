#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "iostream"
#define EXTERN extern
#include "basic.h"
#include "structN.h"
#include "coderN.h"
#include "bmeN.h"
#include "util_filtering.h"
#include "miscN.h"
#include "hvsbm_fullres.h"
#include "bme_tools.h"
#include "Choisubband.h"

#include "mv_ec.h"

#define   FixOL     0           /* Overlap mode for FSBM */
#define   VarOL     1           /* Overlap mode for HVSBM */

#define   bit_coef  1

void print_time( double sc );

float getFrameSATD(float *fr, int hor, int ver)
{
  int i;
  int lenx, leny;
  float *image;
  float sum;

  // initialization
  image = (float*)getarray(hor * ver, sizeof(float), "image");

  for (i = 0; i < hor * ver; i++) {
    image[i] = fr[i];
  }

  // wavelet decomposition
  lenx = hor;
  leny = ver;
  while (!((lenx % 2) || (leny % 2))) {
    // use normalized 9/7 filter
    analysis(image, 0, 0, lenx, leny, hor, ver, 7);
    
    lenx = (lenx >> 1);
    leny = (leny >> 1);
  }

  // generate SAD
  sum = 0.;
  for (i = 0; i < hor * ver; i++) {
    sum = (image[i] < 0) ? (sum - image[i]) : (sum + image[i]);
  }

  // tidy up
  free(image);
  
  return sum;
}


/*
 * calculate MCP_MSE
 * frame1 : current  
 * frame0 : reference 
 * (cx,cy): the coordinator of the upper left corner 
 * mvx    : horizontal motion vector
 * mvy    : vertical motion vector
 * xblk   : real block width
 * yblk   : real block height
 * hor    : frame width
 * ver    : frame height
 */
float
MCP_MSE( float *frame_cur, int cx, int cy, 
         float *frame_ref1, float mvx1, float mvy1,
         float *frame_ref2, float mvx2, float mvy2,
         int xblk, int yblk, int hor, int ver )
{
  int y, x, m;
  float px1, py1;
  float px2, py2;
  float diff, ptemp, mean, sum;

  m = cy * hor + cx;

  sum = 0.;
  mean = 0.;

  if (frame_ref2 == NULL) {
    px1 = cx - mvx1;
    py1 = cy - mvy1;

    for( y = 0; y < yblk; y++ ) {
      for( x = 0; x < xblk; x++ ) {
        ptemp = interpolate( px1 + x, py1 + y, frame_ref1, hor, ver, TYPE );
        diff = frame_cur[m] - ptemp;
        sum += diff * diff;
        mean += diff;
        m++;
      }                   /* x */
      m += hor - xblk;
    }                     /* y */
  } else {
    px1 = cx - mvx1;
    py1 = cy - mvy1;
    px2 = cx - mvx2;
    py2 = cy - mvy2;

    for( y = 0; y < yblk; y++ ) {
      for( x = 0; x < xblk; x++ ) {
        ptemp = interpolate( px1 + x, py1 + y, frame_ref1, hor, ver, TYPE ) +
          interpolate( px2 + x, py2 + y, frame_ref2, hor, ver, TYPE );
        diff = frame_cur[m] - 0.5f * ptemp;
        sum += diff * diff;
        mean += diff;
        m++;
      }                   /* x */
      m += hor - xblk;
    }                     /* y */
  }
  sum /= xblk * yblk;

  return sum;
}

/**************************************************************************
 *                             find_diff()                                 
 * calculate differential subband and variance of blocks at all tree nodes and decide mode.
 * fmv     : motion vector of a block at a tree node
 * frame1  : curren  
 * frame0  : reference 
 * (cx,cy) : the coordinator of the upper left corner
 * xblock  : block width
 * yblock  : block height
 * hor     : frame width
 * ver     : frame height
 * t_level : temporal decomposition level (begin with 0)
 *************************************************************************/
void
find_diff( vector_ptr fmv, float *frame1, float *frame2, float *D,
          int cx, int cy, int xblock, int yblock, int hor, int ver, 
          int t_level, int type )
{
  int xblk, yblk;
  float pfx, pfy, ptemp1, ptemp2, diff;
  float px1, py1, px2, py2;

  int x,y,m,m0;

  assert(fmv && frame1);

  m = cy * hor + cx;
  m0 = 0;

  xblk = ( cx + xblock <= hor ) ? xblock : hor - cx; // real block width
  yblk = ( cy + yblock <= ver ) ? yblock : ver - cy; // real block height

  if( xblk <= 0 || yblk <= 0 ) { /* if block size if null, then return */
    return;
  }

  px1 = cx - fmv->mvx;
  py1 = cy - fmv->mvy;
  px2 = cx + fmv->mvx;
  py2 = cy + fmv->mvy;

  if(type == 0){
	  for( y = 0; y < yblk; y++ ) {
		  for( x = 0; x < xblk; x++ ) {

			ptemp2 = interpolate( px2 + x, py2 + y, frame2, hor, ver, TYPE );
			diff = ptemp2 - frame1[m];
			D[m0] = diff;
			m++;
			m0++;
		  }                   /* x */
		  m += hor - xblk;
	  } 
  }else{
	assert(type == 1);
	
	for( y = 0; y < yblk; y++ ) {
		for( x = 0; x < xblk; x++ ) {
			ptemp1 = interpolate( px1 + x, py1 + y, frame1, hor, ver, TYPE );

			ptemp2 = interpolate( px2 + x, py2 + y, frame2, hor, ver, TYPE );

			diff = ptemp2 - ptemp1;
			D[m0] = diff;
			m++;
			m0++;
		}                   /* x */
		m += hor - xblk;
	} 
  }

  assert( m0 == (xblk*yblk) );

}


float
get_sad( float mvx, float mvy, float *D1, float *D2, int *Dx,
          int xblk, int yblk, int t_level, int type )
{
	int i,j,m;
	float pfx,pfy;
	int xpos,ypos;
	float ptemp, diff;

	float sum = 0;

	for(m = 0 ; m < (xblk * yblk) ; m ++){
		xpos = (m % xblk);
		ypos = (int)(m/yblk);

		pfx = xpos - mvx;
		pfy = ypos - mvy;

		if(xpos >= 0 && xpos < xblk && pfx >= 0 && pfx < xblk &&
			ypos >= 0 && ypos < yblk && pfy >= 0 && pfy < yblk ){
			ptemp = interpolate( pfx, pfy, D2, xblk, yblk, TYPE );

			diff = D1[ypos * xblk + xpos] - ptemp;
		}else{
			diff = D1[ypos * xblk + xpos];
		}
		
		sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);

		if(type == 1)
			Dx[ypos * xblk + xpos] = (int)fabs(diff);
	}

	return sum;
}


int
two_comp_est( float *mvx, float *mvy, vector_ptr fmv, float *frame_cur, float *frame_ref1, float *frame_ref2, float *SAD,
				int cx, int cy, int xblk, int yblk, int hor, int ver, int t_level ){

	float i,j;
	int m,n, a,b, real = 0;
	float pfx1, pfy1, ptemp1,pfx2, pfy2, ptemp2;
	float int_mvx, int_mvy, get_mvx, get_mvy, tpx, tpy;
	float best_sad,sum;
	int well_count, cum_count = 0;

	float diff;

	int range = 12;

	float thres = 2.5;

	printf("fmv1->mvx = %f, fmv1->mvy = %f\n",fmv->mvx,fmv->mvy);

	*SAD = (float)HUGE_VAL;

	for(i=(-1)*range ; i <= range; i ++){
		for(j=(-1)*range ; j <= range; j ++){

			tpx = fmv->mvx + i;
			tpy = fmv->mvy + j;
			well_count = 0;

			m = cy * hor + cx;
			best_sad = 0;

			for( a = 0; a < yblk; a ++ ){
				for( b = 0; b < xblk; b ++ ){
					pfx1 = (float)cx - tpx + b;
					pfy1 = (float)cy - tpy + a;
					pfx2 = (float)cx + tpx + b;
					pfy2 = (float)cy + tpy + a;

//					printf("pfx1 = %f, pfy1 = %f, pfx2 = %f, pfy2 = %f\n",pfx1,pfy1,pfx2,pfy2);

					if(  pfx1 >= hor || pfx1 < 0 || pfy1 >= ver || pfy1 < 0 || pfx2 >= hor || pfx2 < 0 || pfy2 >= ver || pfy2 < 0 ){
						printf("jump out!\n");
						goto sign_out1;
					}

					ptemp1 = interpolate( pfx1, pfy1, frame_ref1, hor, ver, TYPE );
					ptemp2 = interpolate( pfx2, pfy2, frame_ref2, hor, ver, TYPE );

//					printf("fr_cur = %f, ptemp1 = %f, ptemp2 = %f\n",frame_cur[m],ptemp1,ptemp2);

					diff = frame_cur[m] - 0.5f * (ptemp1 + ptemp2);

					best_sad = ( diff < 0 ) ? (best_sad - diff) : (best_sad + diff);

					if( fabs(diff) <= thres ){
						well_count ++;
//						printf("diff = %f\n",diff);
					}

					m++;
				}
				m += hor - xblk;
			}

			if( well_count > cum_count ){
				cum_count = well_count;
				*mvx = tpx;
				*mvy = tpy;
			}

			if(best_sad < *SAD)
				*SAD = best_sad;

sign_out1:		;

		} //j
	} //i

	int_mvx = *mvx;
	int_mvy = *mvy;

//sub-pixel search
	for(i= -0.75 ; i <= 0.75; i += 0.25){
		for(j= -0.75 ; j <= 0.75; j += 0.25){

			tpx = int_mvx + i;
			tpy = int_mvy + j;
			well_count = 0;

			m = cy * hor + cx;
			best_sad = 0;

			for( a = 0; a < yblk; a ++ ){
				for( b = 0; b < xblk; b ++ ){
					pfx1 = (float)cx - tpx + b;
					pfy1 = (float)cy - tpy + a;
					pfx2 = (float)cx + tpx + b;
					pfy2 = (float)cy + tpy + a;

//					printf("pfx1 = %f, pfy1 = %f, pfx2 = %f, pfy2 = %f\n",pfx1,pfy1,pfx2,pfy2);

					if(  pfx1 >= hor || pfx1 < 0 || pfy1 >= ver || pfy1 < 0 || pfx2 >= hor || pfx2 < 0 || pfy2 >= ver || pfy2 < 0 ){
						printf("jump out!\n");
						goto sign_out2;
					}

					ptemp1 = interpolate( pfx1, pfy1, frame_ref1, hor, ver, TYPE );
					ptemp2 = interpolate( pfx2, pfy2, frame_ref2, hor, ver, TYPE );

					diff = frame_cur[m] - 0.5f * (ptemp1 + ptemp2);

					best_sad = ( diff < 0 ) ? (best_sad - diff) : (best_sad + diff);

					if( fabs(diff) <= thres ){
						well_count ++;
//						printf("diff = %f\n",diff);
					}

					m++;
				}
				m += hor - xblk;
			}

			if( well_count > cum_count ){
				cum_count = well_count;
				*mvx = tpx;
				*mvy = tpy;
			}

			if(best_sad < *SAD)
				*SAD = best_sad;

sign_out2:		;

		} //j
	} //i

	printf("cx = %d, cy = %d, xblk = %d, yblk = %d, cum_count = %d, \nsearch sad = %f\n",cx,cy,xblk,yblk,cum_count,*SAD);

	if( cum_count > (xblk*yblk/3) ){
		real = 1;
		assert(*SAD != (float)HUGE_VAL);
	}

	return real;
}

/**************************************************************************
 *                             local_diff_search()                                 
 * calculate differential subband and variance of blocks at all tree nodes and decide mode.
 * fmv     : motion vector of a block at a tree node
 * frame1  : curren  
 * frame0  : reference 
 * (cx,cy) : the coordinator of the upper left corner
 * xblock  : block width
 * yblock  : block height
 * hor     : frame width
 * ver     : frame height
 * t_level : temporal decomposition level (begin with 0)
 *************************************************************************/
float
local_diff_search( float *mvx, float *mvy, float *D1, float *D2, int *Dx,
          int xblk, int yblk, int hor, int ver, int t_level )
{
  float px2, py2;
  float i,j;
  float pfx, pfy, ptemp;
  float px1, py1;
  float best_sad,sum;

  int mvx_int, mvy_int;

  int xpos,ypos, xpos_ref, ypos_ref;
  float diff;

  int x,y,m,m0;

  int range = 15;

  if( xblk <= 0 || yblk <= 0 ) { /* if block size if null, then return */
    return (float)HUGE_VAL;
  }

  best_sad = 0;

  for(m = 0 ; m < (xblk * yblk) ; m ++){
	ptemp = D1[m] - D2[m];
	best_sad = ( ptemp < 0 ) ? (best_sad - ptemp) : (best_sad + ptemp);
	Dx[m] = (int)fabs(ptemp);
  }

  printf("original sad = %f\n\n",best_sad);

  *mvx = 0;
  *mvy = 0;

//integer position search
  for(i=(-1)*range ; i <= range; i ++){
	for(j=(-1)*range ; j <= range; j ++){

		sum = 0;

		for(m = 0 ; m < (xblk * yblk) ; m ++){
			xpos = (m % xblk);
			ypos = (int)(m/yblk);

			xpos_ref = xpos - i;
			ypos_ref = ypos - j;

			if(xpos >= 0 && xpos < xblk && xpos_ref >= 0 && xpos_ref < xblk &&
				ypos >= 0 && ypos < yblk && ypos_ref >= 0 && ypos_ref < yblk ){
				diff = D1[ypos * xblk + xpos] - D2[ypos_ref * xblk + xpos_ref];
			}else{
				diff = D1[ypos * xblk + xpos];
			}
		
			sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
		}

		if(sum < best_sad){
			best_sad = sum;
			*mvx = i;
			*mvy = j;
			printf("best_sad = %f\n",best_sad);
		}
	}
  }

  mvx_int = *mvx;
  mvy_int = *mvy;

//quarter-pixel search
  printf("\n\nsub-pixel search\n\n");
  for(i= -0.75 ; i <= 0.75; i += 0.25){
	for(j= -0.75 ; j <= 0.75; j += 0.25){

		sum = 0;

		for(m = 0 ; m < (xblk * yblk) ; m ++){
			xpos = (m % xblk);
			ypos = (int)(m/yblk);

			pfx = xpos - mvx_int - i;
			pfy = ypos - mvy_int - j;

			if(xpos >= 0 && xpos < xblk && pfx >= 0 && pfx < xblk &&
				ypos >= 0 && ypos < yblk && pfy >= 0 && pfy < yblk ){
				ptemp = interpolate( pfx, pfy, D2, xblk, yblk, TYPE );

				diff = D1[ypos * xblk + xpos] - ptemp;
			}else{
				diff = D1[ypos * xblk + xpos];
			}
		
			sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
		}

		if(sum < best_sad){
			best_sad = sum;
			*mvx = mvx_int + i;
			*mvy = mvy_int + j;
			printf("best_sad = %f\n",best_sad);
		}
	}
  }

  printf("\nmvx = %f, mvy = %f\n\n",*mvx,*mvy);

//save difference block

  for(m = 0 ; m < (xblk * yblk) ; m ++){
	xpos = (m % xblk);
	ypos = (int)(m/yblk);

	pfx = xpos - *mvx;
	pfy = ypos - *mvy;

	if(xpos >= 0 && xpos < xblk && pfx >= 0 && pfx < xblk &&
		ypos >= 0 && ypos < yblk && pfy >= 0 && pfy < yblk ){

		ptemp = interpolate( pfx, pfy, D2, xblk, yblk, TYPE );

		diff = D1[ypos * xblk + xpos] - ptemp;
	}else{
		diff = D1[ypos * xblk + xpos];
	}

	Dx[m] = fabs(diff);
  }

  return best_sad;
}

/**************************************************************************
 *                             find_MSE()                                 
 * calculate MSE and variance of blocks at all tree nodes and decide mode.
 * fmv     : motion vector of a block at a tree node
 * frame1  : curren  
 * frame0  : reference 
 * (cx,cy) : the coordinator of the upper left corner
 * xblock  : block width
 * yblock  : block height
 * hor     : frame width
 * ver     : frame height
 * t_level : temporal decomposition level (begin with 0)
 *************************************************************************/

void
find_MSE( vector_ptr fmv1, vector_ptr fmv2,
          float *frame_cur, float *frame_ref1, float *frame_ref2,
          int cx, int cy, int xblock, int yblock, int hor, int ver, 
          int t_level )
{
  int xblk, yblk, px1, py1, px2, py2;
  float pfx, pfy, ref_var;

  assert(fmv1 && frame_ref1);
  assert((fmv2 && frame_ref2) || (!fmv2 && !frame_ref2));

  temporal_filter();

  xblk = ( cx + xblock <= hor ) ? xblock : hor - cx; // real block width
  yblk = ( cy + yblock <= ver ) ? yblock : ver - cy; // real block height

  if( xblk <= 0 || yblk <= 0 ) { /* if block size if null, then return */
    return;
  }

  if (fmv2 == NULL) {
    fmv1->mse = MCP_MSE( frame_cur, cx, cy, frame_ref1, fmv1->mvx, fmv1->mvy,
                         NULL, 0., 0., xblk, yblk, hor, ver );
    fmv1->var = variance( frame_cur + cy * hor + cx, xblk, yblk, hor );
    
    pfx = cx - fmv1->mvx;
    pfy = cy - fmv1->mvy;
    position( &px1, &py1, pfx, pfy, fmv1->mvx, fmv1->mvy, hor, ver );
    
    ref_var  = variance( frame_ref1 + py1 * hor + px1, xblk, yblk, hor );

//	printf("fmv1->mse = %f, fmv1->var = %f, ref_var = %f\n",fmv1->mse,fmv1->var,ref_var);
    
    if( ( fmv1->mse < IBLOCK_FACTOR * fmv1->var && 
          fmv1->mse < IBLOCK_FACTOR * ref_var )
        || fmv1->mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ) 
      {
        fmv1->lifting_mode = CONNECTED;
      } 
    else 
      {
        fmv1->lifting_mode = PREDICTED;
      }
  } else {
    fmv1->mse = fmv2->mse = MCP_MSE( frame_cur, cx, cy, 
                                     frame_ref1, fmv1->mvx, fmv1->mvy,
                                     frame_ref2, fmv2->mvx, fmv2->mvy,
                                     xblk, yblk, hor, ver );
    fmv1->var = variance( frame_cur + cy * hor + cx, xblk, yblk, hor );
    
    pfx = cx - fmv1->mvx;
    pfy = cy - fmv1->mvy;
    position( &px1, &py1, pfx, pfy, fmv1->mvx, fmv1->mvy, hor, ver );
    pfx = cx - fmv2->mvx;
    pfy = cy - fmv2->mvy;
    position( &px2, &py2, pfx, pfy, fmv2->mvx, fmv2->mvy, hor, ver );
    
    ref_var  =  0.5f *
      (variance( frame_ref1 + py1 * hor + px1, xblk, yblk, hor ) +
       variance( frame_ref2 + py2 * hor + px2, xblk, yblk, hor ));

//	printf("fmv1->mse = %f, fmv1->var = %f, ref_var = %f\n",fmv1->mse,fmv1->var,ref_var);
    
    if( ( fmv1->mse < IBLOCK_FACTOR * fmv1->var && 
          fmv1->mse < IBLOCK_FACTOR * ref_var )
        || fmv1->mse < NOISE_VAR * pow( LPW4[1], ( float )t_level ) ) 
      {
        fmv1->lifting_mode = CONNECTED;
        fmv2->lifting_mode = CONNECTED;
      } 
    else 
      {
        fmv1->lifting_mode = PREDICTED;
        fmv2->lifting_mode = PREDICTED;
      }
  }      
}


/*
 *                            block_matching()                              
 * block-based motion estimation 块为基础的运动补偿
 * fr1-- current  fr0-- reference                                      
 * dist: temporal distance between fr1 and fr0
 * level: temporal decomposition level (begin with 0)
 */
void
block_matching( vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3, vector_ptr fmv4, 
                YUVimage *fr_cur, YUVimage *fr_ref1, YUVimage *fr_ref2, 
                enum FLAG *sc1, enum FLAG *sc2,
                videoinfo info, int t_level, int dist, int subpel, float *upframe1, float *upframe2, float *upsamp_x )
{

  switch ( info.ME ) {
  case 3:
    rdme(fmv1, fmv2, fmv3, fmv4, fr_cur, fr_ref1, fr_ref2, sc1, sc2, t_level, dist, info, upframe1, upframe2, upsamp_x);// 运动估计的rd决策
    break;
  default:
    printf( "error in blockmatching()\n" );
    exit( 1 );
  }
}


static float first_pred[256*256];

static int *spiral_laplength;
static int **spiral_x, **spiral_y;

void
initialize_spiral_search(int search_range)
{
  int i, j, k;
  int len;

  spiral_laplength = (int*)getarray(search_range + 1, sizeof(int), 
                                    "spiral_laplength");
  spiral_x = (int**)getarray(search_range + 1, sizeof(int*), "spiral_x");
  spiral_y = (int**)getarray(search_range + 1, sizeof(int*), "spiral_y");

  spiral_x[0] = (int*)getarray(1, sizeof(int), "spiral_x[0]");
  spiral_y[0] = (int*)getarray(1, sizeof(int), "spiral_y[0]");

  spiral_laplength[0] = 1;
  spiral_x[0][0] = 0;
  spiral_y[0][0] = 0;

  for (i = 1, len = 8; i <= search_range; i++, len += 8) {
    spiral_laplength[i] = len;

    spiral_x[i] = (int*)getarray(len, sizeof(int), "spiral_x[i]");
    spiral_y[i] = (int*)getarray(len, sizeof(int), "spiral_y[i]");

    k = 0;

    spiral_x[i][k] =  0; spiral_y[i][k++] = -i;
    spiral_x[i][k] =  i; spiral_y[i][k++] =  0;
    spiral_x[i][k] =  0; spiral_y[i][k++] =  i;
    spiral_x[i][k] = -i; spiral_y[i][k++] =  0;

    for (j = 1; j < i; j++) {
      spiral_x[i][k] =  j; spiral_y[i][k++] = -i;
      spiral_x[i][k] =  i; spiral_y[i][k++] =  j;
      spiral_x[i][k] = -j; spiral_y[i][k++] =  i;
      spiral_x[i][k] = -i; spiral_y[i][k++] = -j;

      spiral_x[i][k] = -j; spiral_y[i][k++] = -i;
      spiral_x[i][k] =  i; spiral_y[i][k++] = -j;
      spiral_x[i][k] =  j; spiral_y[i][k++] =  i;
      spiral_x[i][k] = -i; spiral_y[i][k++] =  j;
    }

    spiral_x[i][k] =  i; spiral_y[i][k++] =  i;
    spiral_x[i][k] = -i; spiral_y[i][k++] =  i;
    spiral_x[i][k] =  i; spiral_y[i][k++] = -i;
    spiral_x[i][k] = -i; spiral_y[i][k++] = -i;

    assert(k == len);
  }
}

void
clean_up_spiral_search(int search_range)
{
  int i;

  for (i = 0; i <= search_range; i++) {
    free(spiral_x[i]);
    free(spiral_y[i]);
  }
  free(spiral_x);
  free(spiral_y);
  free(spiral_laplength);
}

void
full_search_fast( float *mvx1, float *mvy1, 
                  float *frame_cur, float *frame_ref1, float *frame_ref2,
                  float *upframe1, float *upframe2,
                  float mvx2, float mvy2, 
                  int cx, int cy, int xblock, int yblock, int maxx, int maxy,
                  int hor, int ver, int subpel, float lambda, 
                  float pmvx, float pmvy, int ctx_x, int ctx_y,
                  float *sad_cost, float *bit_cost, float *total_cost,
                  int do_parallel) // add parallel mode. mwi 
     /* frame1-- current  frame0-- reference */
{
  int i, scale, dx, dy, px, py, px2, py2, m, n, xblk, yblk;
  int n2, ppx, ppy; // add parallel mode. mwi 
  //  int upver, uphor;
  int hx, hy;
  float px1, py1;
  float ppx2, ppy2; // add parallel mode. mwi 
  int x, y;
#ifdef PRE_INTERPOL
  int inthor, intver, step, hx2, hy2;
#endif
  float sum, bits;
  float hmvx, hmvy, fmvx, fmvy;
  float *tmp_first_pred;
  int center_x, center_y, lap, pos, cont;

  assert(sad_cost != NULL && bit_cost != NULL && total_cost != NULL);

  // add parallel mode. mwi 
  if (do_parallel == YES) {
    assert(frame_ref2 != NULL);
    
    if (upframe1 != NULL) {
      assert(upframe2 != NULL);
    }
  }

#ifdef PRE_INTERPOL
  step = (1 << MY_MIN(PRE_INTERPOL, subpel)) ;
  step = step << ADD_SUB; //Added on 02.26.2017
  inthor = (hor - 1) * step + 1;
  intver = (ver - 1) * step + 1;
#endif

  xblk = ( cx + xblock <= hor ) ? xblock : hor - cx;
  yblk = ( cy + yblock <= ver ) ? yblock : ver - cy;

//  printf("cx = %d, cy = %d, xblk = %d, yblk = %d\n",cx,cy,xblk,yblk);

  assert(xblk > 0 && yblk > 0);

  /* calculate existing prediction */
  if ( (frame_ref2 != NULL) && (do_parallel == NO) ) {
    tmp_first_pred = first_pred;
    m = 0;
    n = cy * hor + cx;
    for (y = cy; y < cy + yblk; y++) {
      for (x = cx; x < cx + xblk; x++) { 
        px1 = float(x) - mvx2;
        py1 = float(y) - mvy2;
	
        if (px1 >= 0 && px1 <= (hor - 1) && py1 >= 0 && py1 <= (ver - 1)) {
          if (upframe2 && 
              floor(px1*step) == px1*step && floor(py1*step) == py1*step) {
            // pre-interpolated
            tmp_first_pred[m] = upframe2[int((py1 * inthor + px1) * step)];
          } else {
            // not yet interpolated
			assert(0);
            tmp_first_pred[m] = FIRinterpolate(px1, py1, frame_ref2, FIR_LEN, hor, ver);
          }
        } else {
          tmp_first_pred[m] = frame_cur[n];
        }

        m++;
        n++;
      }
      n += hor - xblk;
    }
    assert(m == xblk * yblk);
  } else {
    tmp_first_pred = NULL;
  }

  // changed by Yongjun Wu
  // search center
  center_x = (int)pmvx;  // (pmvx > 0) ? (int)(pmvx + 0.5f) : (int)(pmvx - 0.5f);
  center_y = (int)pmvy;  // (pmvy > 0) ? (int)(pmvy + 0.5f) : (int)(pmvy - 0.5f);

  // first candidate
  // position in reference frame
  px = cx - center_x;
  py = cy - center_y;
  assert(px >= 0 && px <= (hor - xblk) && py >= 0 && py <= (ver - yblk));

  // bit cost
  *bit_cost = bit_coef * get_bit_cost(lambda, (float)center_x, (float)center_y, 
                           pmvx, pmvy, ctx_x, ctx_y, subpel);

  // SAD
  m = cy * hor + cx;
  n = py * hor + px;
  if (do_parallel == NO) { // add parallel mode. mwi 
    *sad_cost = MCP_Error2(frame_cur, m, frame_ref1, n, xblk, yblk, hor,
                           tmp_first_pred, (float)HUGE_VAL,
                           NULL, 0); // add parallel mode. mwi 
  } else {
    // mirrored position in the second reference frame
    ppx = cx + center_x;
    ppy = cy + center_y;
    if (!(ppx >= 0 && ppx <= (hor - xblk) && ppy >= 0 && ppy <= (ver - yblk))){
      *sad_cost = (float)HUGE_VAL;
    } else {
      n2 = ppy * hor + ppx;
      *sad_cost = MCP_Error2(frame_cur, m, frame_ref1, n, xblk, yblk, hor,
                             tmp_first_pred, (float)HUGE_VAL,
                             frame_ref2, n2); // add parallel mode. mwi 
    }
  }    

  // get it
  if (*sad_cost != (float)HUGE_VAL) {
    *total_cost = *sad_cost + *bit_cost;
    *mvx1 = (float)center_x;
    *mvy1 = (float)center_y;
  } else {
    *total_cost = (float)HUGE_VAL;
    *mvx1 = (float)HUGE_VAL;
    *mvy1 = (float)HUGE_VAL;
  }

  // further candidates
  // spiral laps
  assert(maxx >= 1);
  for (lap = 1; lap <= maxx; lap++) {
    cont = 0;
    
    // lap positions
    for (pos = 0; pos < spiral_laplength[lap]; pos++) {
      // candidate vector
      dx = center_x + spiral_x[lap][pos];
      dy = center_y + spiral_y[lap][pos];

      // position in reference frame
      px = cx - dx;
      py = cy - dy;
      ppx = cx + dx;  // add parallel mode. mwi 
      ppy = cy + dy;  // add parallel mode. mwi 

      // boundary check // add parallel mode. mwi 
      if(px >= 0 && px <= (hor - xblk) && py >= 0 && py <= (ver - yblk)) { 
        if ((do_parallel == NO) || ((do_parallel == YES)&&
                                    (ppx >= 0 && ppx <= (hor - xblk) && ppy >= 0 && ppy <= (ver - yblk))) ) {
        
          // bit cost
          bits = bit_coef * get_bit_cost(lambda, (float)dx, (float)dy, pmvx, pmvy, 
                              ctx_x, ctx_y, subpel);
          if (bits < *total_cost) {
            cont = 1;
          
            // SAD
            m = cy * hor + cx;
            n = py * hor + px;
            if (do_parallel == NO) { // add parallel mode. mwi 
              sum = MCP_Error2(frame_cur, m, frame_ref1, n, xblk, yblk, hor,
                               tmp_first_pred, (*total_cost) - bits,
                               NULL, 0); // add parallel mode. mwi 
            } else {
              n2 = ppy * hor + ppx;
              sum = MCP_Error2(frame_cur, m, frame_ref1, n, xblk, yblk, hor,
                               tmp_first_pred, (*total_cost) - bits,
                               frame_ref2, n2); // add parallel mode. mwi 
              
            }
            if( sum + bits < *total_cost ) {
              *sad_cost = sum;
              *bit_cost = bits;
              *total_cost = sum + bits;
              *mvx1 = (float)dx;
              *mvy1 = (float)dy;
            }
          }
        }
      }
    }
    
    // found any matches?
    if (*total_cost == (float)HUGE_VAL) cont = 1;

    // vector cost exceeded minimum total cost in last lap?
    if (!cont) break;
  }

  // still no match?
  if (*total_cost == (float)HUGE_VAL) {
    *sad_cost = (float)HUGE_VAL;
    *bit_cost = (float)HUGE_VAL;
    *mvx1 = 0.0;
    *mvy1 = 0.0;
//	assert(0);
    return;
  }

  /* half-pixel search                */
  //  uphor = hor;
  //  upver = ver;
  m = cy * hor + cx;

  for( i = 0; i < subpel; i++ ) {
    hmvx = 0.;
    hmvy = 0.;
    px1 = cx - *mvx1;
    py1 = cy - *mvy1;
    ppx2 = cx + *mvx1; // add parallel mode. mwi 
    ppy2 = cy + *mvy1; // add parallel mode. mwi 
    scale = (1 << i);

#ifdef PRE_INTERPOL
    if (i >= MY_MIN(PRE_INTERPOL, subpel)) {
#endif

      for( hy = -1; hy <= 1; hy++ ) {
        for( hx = -1; hx <= 1; hx++ ) {
          
          if (!(hx == 0 && hy == 0)) {
            
            fmvx = float(hx / 2. / float(scale));
            fmvy = float(hy / 2. / float(scale));
            
            if (px1 + fmvx >= 0 && px1 + xblk + fmvx <= hor &&
                py1 + fmvy >= 0 && py1 + yblk + fmvy <= ver) {
              // add parallel mode. mwi 
              if ((do_parallel == NO) || ((do_parallel == YES)&&
                                          ppx2 - fmvx >= 0 && ppx2 + xblk - fmvx <= hor &&
                                          ppy2 - fmvy >= 0 && ppy2 + yblk - fmvy <= ver)) {
              
              bits = bit_coef * get_bit_cost(lambda, (*mvx1) - fmvx, (*mvy1) - fmvy,
                                  pmvx, pmvy, ctx_x, ctx_y, subpel);
              
              if (bits < *total_cost) {
                if (do_parallel == NO) { // add parallel mode. mwi 
                  sum = Subpel_MCP_Error( frame_cur, m, frame_ref1, 
                                          px1, fmvx, py1, fmvy, 
                                          xblk, yblk, hor, ver, 
                                          tmp_first_pred, 0, xblk,
                                          NULL, 0., 0.); // add parallel mode. mwi 
                } else {
                  sum = Subpel_MCP_Error( frame_cur, m, frame_ref1, 
                                          px1, fmvx, py1, fmvy, 
                                          xblk, yblk, hor, ver, 
                                          tmp_first_pred, 0, xblk,
                                          frame_ref2, ppx2, ppy2);
                }
                if( sum + bits < *total_cost ) {
                  *sad_cost = sum;
                  *bit_cost = bits;
                  *total_cost = sum + bits;
                  hmvx = fmvx;
                  hmvy = fmvy;
                }
              }
            }
            }
          }
        }                         /* hx */
      }                           /* hy */
#ifdef PRE_INTERPOL
    } else {

      px2 = ( int )( px1 * step );
      py2 = ( int )( py1 * step );  
      ppx = ( int )( ppx2 * step ); // add parallel mode. mwi 
      ppy = ( int )( ppy2 * step ); // add parallel mode. mwi   

      // take care: hx/hy is *negative* refinement to mvx1/mvy1
      for( hy = -1; hy <= 1; hy++ ) {
        for( hx = -1; hx <= 1; hx++ ) {

          //          printf("px1 = %f, px2 = %f, hx = %d, hy = %d
          if (!(hx == 0 && hy == 0)) {
            assert(MY_MIN(PRE_INTERPOL, subpel) - i - 1 >= 0);

            hx2 = hx << (MY_MIN(PRE_INTERPOL, subpel + ADD_SUB) - i - 1);//Modified on 03.01.2017
            hy2 = hy << (MY_MIN(PRE_INTERPOL, subpel + ADD_SUB) - i - 1);
            
            fmvx = float(hx / 2. / float(scale));
            fmvy = float(hy / 2. / float(scale));
              
            assert((px1 + xblk + fmvx <= hor) == (px2 + hx2 < (inthor - (xblk - 1) * step)) &&
                   (py1 + yblk + fmvy <= ver) == (py2 + hy2 < (intver - (yblk - 1) * step)));
            // parallel mode
            assert((ppx2 + xblk - fmvx <= hor) == (ppx - hx2 < (inthor - (xblk - 1) * step)) &&
                   (ppy2 + yblk - fmvy <= ver) == (ppy - hy2 < (intver - (yblk - 1) * step)));

            if (px1 + fmvx >= 0 && px1 + xblk + fmvx <= hor &&
                py1 + fmvy >= 0 && py1 + yblk + fmvy <= ver) {
              if ((do_parallel == NO) || 
                  ((do_parallel == YES)&&
                   (ppx2 - fmvx >= 0 && ppx2 + xblk - fmvx <= hor && 
                    ppy2 - fmvy >= 0 && ppy2 + yblk - fmvy <= ver))) {
                
                bits = bit_coef * get_bit_cost(lambda, (*mvx1) - fmvx, (*mvy1) - fmvy,
                                    pmvx, pmvy, ctx_x, ctx_y, subpel);
                
                if (do_parallel == NO) { // add parallel mode. mwi               
                  sum = Subpel_MCP_Error3(frame_cur, m, upframe1,
                                          px2 + hx2, py2 + hy2, 
                                          xblk, yblk, hor, ver,
                                          inthor, intver, step,
                                          tmp_first_pred, (*total_cost) - bits,
                                          NULL, 0, 0); // add parallel mode. mwi 

                } else {
                  sum = Subpel_MCP_Error3(frame_cur, m, upframe1,
                                          px2 + hx2, py2 + hy2, 
                                          xblk, yblk, hor, ver,
                                          inthor, intver, step,
                                          tmp_first_pred, (*total_cost) - bits,
                                          upframe2, ppx-hx2, ppy-hy2);
                }
                
                if( sum + bits < *total_cost ) {
                  
                  *sad_cost = sum;
                  *bit_cost = bits;
                  *total_cost = sum + bits;
                  hmvx = fmvx;
                  hmvy = fmvy;
                }
              }
            }
          }
        }                         /* hx */
      }                           /* hy */
    }
#endif

    (*mvx1) -= hmvx;
    (*mvy1) -= hmvy;
  }

  assert(!(fabs(*mvx1) > 1e-42 && fabs(*mvx1) < 1e-10));
  assert(!(fabs(*mvy1) > 1e-42 && fabs(*mvy1) < 1e-10));
}

/*****************************************************************************/
/*                          generate_child()                                 */
/*****************************************************************************/
void
generate_child( vector_ptr fmv, float *mvx, float *mvy, float *mad )
{
  /* allocate the memory to the child vector pointer */
  /* save the motion vector and estimation error */

  assert(fmv->child == 0);

  fmv->child = 1;
  fmv->child0 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child0" );
  fmv->child1 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child1" );
  fmv->child2 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child2" );
  fmv->child3 = ( vector_ptr ) getarray( 1, sizeof( vector ), "fmv->child3" );

  fmv->child0->parent = fmv;
  fmv->child1->parent = fmv;
  fmv->child2->parent = fmv;
  fmv->child3->parent = fmv;

  fmv->child0->child = 0;
  fmv->child0->lifting_mode = CONNECTED;
  fmv->child0->merge = YES;
  fmv->child0->mvx = mvx[0];
  fmv->child0->mvy = mvy[0];
  fmv->child0->mad = mad[0];

  fmv->child1->child = 0;
  fmv->child1->lifting_mode = CONNECTED;
  fmv->child1->merge = YES;
  fmv->child1->mvx = mvx[1];
  fmv->child1->mvy = mvy[1];
  fmv->child1->mad = mad[1];

  fmv->child2->child = 0;
  fmv->child2->lifting_mode = CONNECTED;
  fmv->child2->merge = YES;
  fmv->child2->mvx = mvx[2];
  fmv->child2->mvy = mvy[2];
  fmv->child2->mad = mad[2];

  fmv->child3->child = 0;
  fmv->child3->lifting_mode = CONNECTED;
  fmv->child3->merge = YES;
  fmv->child3->mvx = mvx[3];
  fmv->child3->mvy = mvy[3];
  fmv->child3->mad = mad[3];
}


/******************************************************************
                  Compute MCP error with half pixel accuracy

   m:        initial position in current frame-frame1
   (px, py): initial position in previous frame- frame0
   (hx, hy): sub-pixel search
   xblk: block width
   yblk: block height
   hor:   frame1 width
   ver:   frame1 height
   add parallel mode: 
   (ppx,ppy):initial position in second ref.frame- frame2
*******************************************************************/
float
Subpel_MCP_Error( float *frame1, int m, float *frame0, float px, float hx,
                  float py, float hy, int xblk, int yblk, int hor, int ver,
                  float *first_pred, int k, int xblk_pred,
                  float *frame2, float ppx, float ppy) // add parallel mode. mwi 
{
  int y, x;
  float diff, ptemp, sum;
  float ptemp2; // add parallel mode. mwi 

  sum = 0.;
  if (frame2 == NULL) {
  if (first_pred == NULL) {
    for( y = 0; y < yblk; y++ ) { /* calculate the error */
      for( x = 0; x < xblk; x++ ) {

        ptemp = interpolate(px + hx + x, py + hy + y, frame0, hor, ver, TYPE);

        diff = frame1[m] - ptemp;
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff); 

        m++;
      }     /* x */
      m += hor - xblk;
    }                             /* y */
  } else {
    for( y = 0; y < yblk; y++ ) { /* calculate the error */
      for( x = 0; x < xblk; x++ ) {

        ptemp = interpolate(px + hx + x, py + hy + y, frame0, hor, ver, TYPE);

        diff = frame1[m] - 0.5f * (first_pred[k] + ptemp);
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff); 
	
        k++;
        m++;
      }     /* x */
      k += xblk_pred - xblk;
      m += hor - xblk;
    }                             /* y */
  }
  } else {  
    // add parallel mode. mwi 
    for( y = 0; y < yblk; y++ ) { /* calculate the error */
      for( x = 0; x < xblk; x++ ) {

        ptemp  = interpolate(px + hx + x, py + hy + y, frame0, hor, ver, TYPE);
        ptemp2 = interpolate(ppx - hx + x, ppy - hy + y, frame2, hor, ver, TYPE);

        diff = frame1[m] - 0.5f * (ptemp + ptemp2);
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff); 
	
        m++;
      }     /* x */
      m += hor - xblk;
    }                             /* y */
  }    

  return sum;
}




/******************************************************************
                  Compute MCP error with half pixel accuracy

   m:        initial position in current frame-frame1
   (px, py): initial integer position in previous frame
   upframe:  upsampled previous frame. (generated by Interpolate_frame)
   (hx, hy): half-pixel search
   xblk: block width
   yblk: block height
   hor:   frame1 width
   ver:   frame1 height
   uphor: upsampled frame width
   upver: upsampled frame height
   step : step size in reference frame (2 for once upsampled image, 4 for twice upsampled image)
*******************************************************************/
float
Subpel_MCP_Error2( float *frame1, int m, float *upframe, int px, int hx,
                   int py, int hy, int xblk, int yblk, int hor, int ver,
                   int uphor, int upver, int step,
                   float *first_pred, int k, int xblk_pred )

{
  int y, x;
  float diff, sum;

  sum = 0.;
  if (first_pred == NULL) {
    for( y = 0; y < yblk; y++ ) { /* calculate the error */
      for( x = 0; x < xblk; x++ ) {

        diff =
          frame1[m] - upframe[( py * 2 + y * step + hy ) * uphor +
                              ( px * 2 + x * step ) + hx];
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);  

        m++;
      }                           /* x */
      m += hor - xblk;
    }                             /* y */
  } else {
    for( y = 0; y < yblk; y++ ) { /* calculate the error */
      for( x = 0; x < xblk; x++ ) {

        diff =
          frame1[m] - 0.5f * (upframe[( py * 2 + y * step + hy ) * uphor +
                                     ( px * 2 + x * step ) + hx] + 
                             first_pred[k]);
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);  

        k++;
        m++;
      }                           /* x */
      k += xblk_pred - xblk;
      m += hor - xblk;
    }                             /* y */
  }

  return sum;
}

float
Subpel_MCP_Error3( float *frame1, int m, float *upframe, int px, int py, 
                   int xblk, int yblk, int hor, int ver,
                   int uphor, int upver, int step, float *first_pred, 
                   float min_cost,
                   float *upframe2, int ppx, int ppy) // add parallel mode. mwi 

{
  int x, y, k;
  float diff, sum;

  sum = 0.;
  if (upframe2 == NULL) {
    if (first_pred == NULL) {
      for( y = 0; y < yblk; y++ ) { /* calculate the error */
        if (sum >= min_cost) break;
        
        for( x = 0; x < xblk; x++ ) {
          
          diff =
            frame1[m] - upframe[( py + y * step ) * uphor +
                                ( px + x * step ) ];
          sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);  
          
          m++;
        }                           /* x */
        
        m += hor - xblk;
      }                             /* y */
    } else {
      k = 0;
      
      for( y = 0; y < yblk; y++ ) { /* calculate the error */
        if (sum >= min_cost) break;
        
        for( x = 0; x < xblk; x++ ) {
          
          diff =
            frame1[m] - 0.5f * (upframe[( py + y * step ) * uphor +
                                       ( px + x * step ) ] + 
                               first_pred[k]);
          sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);  
          
          k++;
          m++;
        }                           /* x */
        
        m += hor - xblk;
      }                             /* y */
    }
  } else {
    // parallel mode
    for( y = 0; y < yblk; y++ ) { /* calculate the error */
      if (sum >= min_cost) break;
      
      for( x = 0; x < xblk; x++ ) {
        
        diff =
          frame1[m] - 0.5f * (upframe[( py + y * step ) * uphor +
                                     ( px + x * step ) ] + 
                             upframe2[( ppy + y * step ) * uphor +
                                     ( ppx + x * step ) ]);
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);  

        m++;
      }                           /* x */
      
      m += hor - xblk;
    }                             /* y */
  }
    
  return sum;
}

float
Subpel_MCP_Error4( float *frame_cur,  
                   int cx, int cy, int xblock, int yblock, 
                   float *frame_ref1, float *upframe_ref1, 
                   float mvx1, float mvy1,
                   float *frame_ref2, float *upframe_ref2, 
                   float mvx2, float mvy2,
                   int hor, int ver, int subpel )
{
  int x, y, m;
  int xblk, yblk;
  float px1, py1, px2, py2;
  float pred1, pred2;
  float diff, sum;
#ifdef PRE_INTERPOL
  int step, inthor, intver;

  step = (1 << MY_MIN(PRE_INTERPOL, subpel)) ;

  step = step << ADD_SUB; //Added on 02.26.2017

  inthor = (hor - 1) * step + 1;
  intver = (ver - 1) * step + 1;
  
  assert(frame_ref1 && upframe_ref1);
  assert(!frame_ref2 || upframe_ref2);
#endif

  xblk = ( cx + xblock <= hor ) ? xblock : hor - cx;
  yblk = ( cy + yblock <= ver ) ? yblock : ver - cy;

  assert(xblk > 0 && yblk > 0);

  sum = 0.;

  m = cy * hor + cx;
  for (y = cy; y < cy + yblk; y++) {
    for (x = cx; x < cx + xblk; x++) { 
      px1 = float(x) - mvx1;
      py1 = float(y) - mvy1;
      assert(px1 >= 0 && px1 <= (hor - 1) && py1 >= 0 && py1 <= (ver - 1));

#ifdef PRE_INTERPOL
      if (floor(px1*step) == px1*step && floor(py1*step) == py1*step) {
        // pre-interpolated
        pred1 = upframe_ref1[int((py1 * inthor + px1) * step)];
      } else {       
        // not yet interpolated
        pred1 = FIRinterpolate(px1, py1, frame_ref1, 12, hor, ver);
      }
#else
      // no pre-interpolation
      pred1 = FIRinterpolate(px1, py1, frame_ref1, 12, hor, ver);
#endif

      if (frame_ref2) {
        px2 = float(x) - mvx2;
        py2 = float(y) - mvy2;
	
        assert(px2 >= 0 && px2 <= (hor - 1) && py2 >= 0 && py2 <= (ver - 1));
        
#ifdef PRE_INTERPOL
        if (floor(px2*step) == px2*step && floor(py2*step) == py2*step) {
          // pre-interpolated
          pred2 = upframe_ref2[int((py2 * inthor + px2) * step)];
        } else {
          // not yet interpolated
          pred2 = FIRinterpolate(px2, py2, frame_ref2, 12, hor, ver);
        }
#else
        // no pre-interpolation
        pred2 = FIRinterpolate(px2, py2, frame_ref2, 12, hor, ver);
#endif
      } else {
        pred2 = pred1;
      }

      diff = frame_cur[m] - 0.5f * (pred1 + pred2);
      sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);  

      m++;
    }
    m += hor - xblk;
  }
    
  return sum;
}


/******************************************************************
                  Compute MCP error with integer accuracy

   m:     initial position in current frame-frame1
   n:     initial position in previous frame-frame0
   oxblk: block width
   oyblk: block height
   hor:   image width
*******************************************************************/
float
MCP_Error(float *frame1, int m, float *frame0, int n, int oxblk, int oyblk,
          int hor, float *first_pred, int k, int xblk_pred)		//Never applied in practice
{
  int y, x;
  float diff, sum;

  sum = 0.;
  if (first_pred == NULL) {
    for( y = 0; y < oyblk; y++ ) {
      for( x = 0; x < oxblk; x++ ) {
        diff = frame1[m] - frame0[n];
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
        /* mean absolute error */
        m++;
        n++;
      }
      m += hor - oxblk;
      n += hor - oxblk;
    }
  } else {
    for( y = 0; y < oyblk; y++ ) {
      for( x = 0; x < oxblk; x++ ) {
        diff = frame1[m] - 0.5f * (frame0[n] + first_pred[k]);
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
        /* mean absolute error */
        k++;
        m++;
        n++;
      }
      k += xblk_pred - oxblk;
      m += hor - oxblk;
      n += hor - oxblk;
    }
  }      

  return sum;
}

float
MCP_Error2(float *frame1, int m, float *frame0, int n, int oxblk, int oyblk,
           int hor, float *first_pred, float min_cost,
           float *frame2, int n2) // add parallel mode. mwi 
{
  int y, x, k;
  float diff, sum;

  sum = 0.;
  if (frame2 == NULL) { // add parallel mode. mwi 
    if (first_pred == NULL) {
      for( y = 0; y < oyblk; y++ ) {
        if (sum >= min_cost) break;
        
        for( x = 0; x < oxblk; x++ ) {
          diff = frame1[m] - frame0[n];
          sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
          /* mean absolute error */
          m++;
          n++;
        }
        
        m += hor - oxblk;
        n += hor - oxblk;
      }
    } else {
      k = 0;
      for( y = 0; y < oyblk; y++ ) {
        if (sum >= min_cost) break;
        
        for( x = 0; x < oxblk; x++ ) {
          diff = frame1[m] - 0.5f * (frame0[n] + first_pred[k]);
          sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
          /* mean absolute error */
          k++;
          m++;
          n++;
        }
        
        m += hor - oxblk;
        n += hor - oxblk;
      }
    }      
  } else { // add parallel mode. mwi 
    for( y = 0; y < oyblk; y++ ) {
      if (sum >= min_cost) break;
      
      for( x = 0; x < oxblk; x++ ) {
        diff = frame1[m] - 0.5f * (frame0[n] + frame2[n2]);
        sum = ( diff < 0 ) ? (sum - diff) : (sum + diff);
        /* mean absolute error */
        m++;
        n++;
        n2++;
      }

      m += hor - oxblk;
      n += hor - oxblk;
      n2+= hor - oxblk;
    }
    
  }

  return sum;
}



/*
 *                               position1D()                                  
 * find the matching 1D integer position for sub-pixel accurate temporal filtering 
 *  寻找pfx的整数
 */

int
position1D( float pfx, float mvx )
{
  int tmpx;
  if( ( int )pfx == pfx )
    tmpx = ( int )pfx;          // integer pixel 
  else 
	  if( ( ( int )( pfx * 2 ) ) / 2. == pfx )
		tmpx = ( mvx > 0. ) ? ( int )ceil( pfx ) : ( int )floor( pfx );     // half pixel 
	  else
		tmpx = nint( pfx );         // otherwise

  return tmpx;
}


/*
 *                               position()               返回亚像素准确时域滤波对应的整数位置                   
 * find the matching integer position for sub-pixel accurate temporal filtering 
 *
 */
void
position( int *px, int *py, float pfx, float pfy, float mvx, float mvy,
          int hor, int ver )
{
  int tmpx, tmpy;


  tmpx = position1D( pfx, mvx );// 返回整数
  tmpy = position1D( pfy, mvy ); // 返回整数

  // clipping  ???
  if( tmpx < 0 )
    tmpx = 0;
  else if( tmpx > hor - 1 )
    tmpx = hor - 1;

  if( tmpy < 0 )
    tmpy = 0;
  else if( tmpy > ver - 1 )
    tmpy = ver - 1;

  *px = tmpx;
  *py = tmpy;
}


/*****************************************************************************/
/*                                  inbound()                                */
/*****************************************************************************/
int
inbound( float x, float y, int hor, int ver )
{
  if( x >= 0. && x <= ( float )( hor - 1 ) && y >= 0.
      && y <= ( float )( ver - 1 ) )
    return 1;
  else
    return 0;
}


/*****************************************************************************/
/*                               get_cvector                                 */
/*****************************************************************************/
void
get_cvector( float *cmvx, float *cmvy, float *ymvx, float *ymvy, int yhor,
             int yver, int chor, int cver, videoinfo info, int t_level )
{
  int x, y, pos;
  float dx, dy;
  int t1, t2, scale;

  scale = 1 << info.subpel[t_level];

  /* get the motion vector for the chrominance component */
  /* from the luminance vectors */
  /* half */
  /* yfv    dx     t1    dx    cfv */
  /* +-0.5  +-.25    0     0.  0.  */
  /* +-1.0  +-.5     0   +-.5  +-0.5 */
  /* +-1.5  +-.75    0   +-.5  +-0.5 */
  /* +-2.0  +-1.0  +-1   0.    +-1.0 */
  /* +-2.5  +-1.25 +-1   0.    +-1.0 */
  /* +-3.0  +-1.50 +-1   +-.5  +-1.5 */

  /* quarter */
  /* yfv    dx     t1    dx    cfv */
  /* 0.25   0.125   0    0         */
  /* 0.5    0.25    0    0.25      */
  /* 0.75   0.375   0    0.25      */
  /* 1      0.5     0    0.5       */
  /* 1.25   0.625   0    0.5       */
  /* 1.5    0.75    0    0.75      */
  /* 1.75   0.875   0    0.75      */

  if( yver == 2 * cver && yhor == 2 * chor ) {
    for( y = 0; y < cver; y++ ) {
      for( x = 0; x < chor; x++ ) {

        dx = ymvx[( 2 * y ) * yhor + ( 2 * x )];

        if( dx == (float)HUGE_VAL ) 
          cmvx[y * chor + x] = dx;
        else {
          dx /= 2.0;
          t1 = ( int )dx;
          dx = ( float )( ( int )( dx * scale ) % scale ) / ( float )scale;
          cmvx[y * chor + x] = t1 + dx;
        }

        dy = ymvy[( 2 * y ) * yhor + ( 2 * x )];
        if( dy == (float)HUGE_VAL ) 
          cmvy[y * chor + x] = dy;
        else {
          dy /= 2.0;
          t2 = ( int )dy;
          dy = ( float )( ( int )( dy * scale ) % scale ) / ( float )scale;
          cmvy[y * chor + x] = t2 + dy;

        }
      }
    }
  }                             //420 
  else if( yver == cver && yhor == chor ) {     // 444
    for( y = 0; y < cver; y++ ) {
      for( x = 0; x < chor; x++ ) {
        pos = y * chor + x;
        cmvx[pos] = ymvx[pos];
        cmvy[pos] = ymvy[pos];
      }
    }
  } else {
    printf( "can not handle this case (mctfN.c)\n" );
    exit( 1 );
  }

}


/*****************************************************************************
 *                         block2pixel2()            从块为基础的运动向量变为像素为基础的向量 每个像素执行一次这个函数
 * get dense motion field from block-based motion vectors 
 * except blocks with CONNECTED mode.
 * 
 *****************************************************************************/
void
rec_block2pixel( float *mvx, float *mvy, vector_ptr fmv, int cx, int cy,
                 int xblk, int yblk, int hor, int ver,
                 enum LiftingMode block_mode )
{
  int i, j, xblock, yblock, pos;

  /* change the structure of motion vectors */
  /* from the block-based to the pixel-based */
  /* write the motion vector of the block recursively */

  if( fmv->child ) {
    rec_block2pixel( mvx, mvy, fmv->child0, cx, cy, xblk / 2, yblk / 2, hor,
                     ver, block_mode );
    rec_block2pixel( mvx, mvy, fmv->child1, cx + xblk / 2, cy, xblk / 2,
                     yblk / 2, hor, ver, block_mode );
    rec_block2pixel( mvx, mvy, fmv->child2, cx, cy + yblk / 2, xblk / 2,
                     yblk / 2, hor, ver, block_mode );
    rec_block2pixel( mvx, mvy, fmv->child3, cx + xblk / 2, cy + yblk / 2,
                     xblk / 2, yblk / 2, hor, ver, block_mode );
  } else {
    /* consider the small block around the boundaries */
    xblock = ( cx + xblk <= hor ) ? xblk : hor - cx;
    yblock = ( cy + yblk <= ver ) ? yblk : ver - cy;

    if( xblock <= 0 || yblock <= 0 ) {
      /*      printf("xblock<=0 || yblock<=0 in block2pixel2() !\n");*/
      return;
    }

	// 得到每一个像素点的mv
    if( fmv->lifting_mode == block_mode ) {
	  
      for( i = cy; i < cy + yblock; i++ ) {
        for( j = cx; j < cx + xblock; j++ ) {  //Modified by Yuan Liu
			if(fmv->bi_mode <= 6 || fmv->bi_mode == 8 || (fmv->bi_mode == 7 && fmv->aff_mrg == NO) ){
			  pos = i * hor + j;
			  mvx[pos] = fmv->mvx;
			  mvy[pos] = fmv->mvy;
			}
			else if( (fmv->bi_mode >= 9 && fmv->bi_mode <= 11) || (fmv->bi_mode == 7 && fmv->aff_mrg == YES) ){
			  
			  pos = i * hor + j;
			  mvx[pos] = (fmv->aff2_mvx - fmv->aff1_mvx)*((float)(j-cx))/((float)xblock) + (fmv->aff3_mvx - fmv->aff1_mvx)*((float)(i-cy))/((float)yblock) + fmv->aff1_mvx;
			  mvy[pos] = (fmv->aff2_mvy - fmv->aff1_mvy)*((float)(j-cx))/((float)xblock) + (fmv->aff3_mvy - fmv->aff1_mvy)*((float)(i-cy))/((float)yblock) + fmv->aff1_mvy;
//			  printf("%f\t",mvx[pos]);
			  if( fmv->aff1_mvx == (float)HUGE_VAL || fmv->aff1_mvy == (float)HUGE_VAL ||
				  fmv->aff2_mvx == (float)HUGE_VAL || fmv->aff2_mvy == (float)HUGE_VAL ||
				  fmv->aff3_mvx == (float)HUGE_VAL || fmv->aff3_mvy == (float)HUGE_VAL){

				mvx[pos] = (float)HUGE_VAL;
				mvy[pos] = (float)HUGE_VAL;
			  }

			}
        }
		if(fmv->bi_mode >= 9 && fmv->bi_mode <= 11);
//			printf("\n");
      }
	  if(fmv->bi_mode >= 9 && fmv->bi_mode <= 11);
//		printf("\n\n\n");
    } 
	else {
      for( i = cy; i < cy + yblock; i++ ) {
        for( j = cx; j < cx + xblock; j++ ) {
          pos = i * hor + j;
          mvx[pos] = (float)HUGE_VAL;  
          mvy[pos] = (float)HUGE_VAL;  
        }
      }
    }

  }
}


void
scale_down_mv_comp( float *mv, int hor, int ver, int scale, int precision )
{
  int x, y, pos;
  float dx;
  int t;

  for( y = 0; y < ver; y++ ) {
    for( x = 0; x < hor; x++ ) {

      pos = ( scale * y ) * ( hor * scale ) + ( scale * x );
      dx = mv[pos];

      if( dx == (float)HUGE_VAL )
        mv[y * hor + x] = dx;
      else {
        dx /= scale;
        t = ( int )dx;
        dx =
          ( float )( ( int )( dx * precision ) % precision ) /
          ( float )precision;
        mv[y * hor + x] = t + dx;
      }
    }
  }
}


void
scale_down_mv( float *ymvx, float *ymvy, float *cmvx, float *cmvy,
               videoinfo info, int t_level )
{
  int hor, ver, scale, precision, s_level;
  
  s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
  scale = 1 << s_level;
  precision = 1 << info.subpel[t_level];

  hor = info.ywidth >> s_level;
  ver = info.yheight >> s_level;
  scale_down_mv_comp( ymvx, hor, ver, scale, precision );
  scale_down_mv_comp( ymvy, hor, ver, scale, precision );

  hor = info.cwidth >> s_level;
  ver = info.cheight >> s_level;
  scale_down_mv_comp( cmvx, hor, ver, scale, precision );
  scale_down_mv_comp( cmvy, hor, ver, scale, precision );

}

// 块级mv转换为像素级mv，放进ymvx、ymvy、cmvx、cmvy中
void
blockmv2pixelmv( vector_ptr fmv, float **ymvx, float **ymvy, float **cmvx,
                 float **cmvy, enum LiftingMode block_mode, videoinfo info, int t_level )
{
  int x, y, X, Y, yhor, yver, chor, cver;
  int xnum, ynum, xblk, yblk, s_level;

  /* allocate the memory and initialization */
  s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
  if( s_level > 0 ) {
    info.ywidth  <<= s_level;
    info.yheight <<= s_level;
    info.cwidth  <<= s_level;
    info.cheight <<= s_level;
  }

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;
  xnum = info.xnum[t_level];
  ynum = info.ynum[t_level];
  xblk = info.xblk[t_level];
  yblk = info.yblk[t_level];

  *ymvx = ( float * )getarray( yhor * yver, sizeof( float ), "ymvx" );
  *ymvy = ( float * )getarray( yhor * yver, sizeof( float ), "ymvy" );
  *cmvx = ( float * )getarray( chor * cver, sizeof( float ), "cmvx" );
  *cmvy = ( float * )getarray( chor * cver, sizeof( float ), "cmvy" );

//  printf("s_level= %d\n",s_level);
  /* convert the block-based vectors to the pixel-based vectors 转化块为基础的向量为像素级向量 */
  for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {
    for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {// 
      rec_block2pixel( *ymvx, *ymvy, &fmv[Y * xnum + X], x, y, xblk, yblk,
                       yhor, yver, block_mode );
    }
  }

  /* get the chrominance motion vector from the luminance vector 从亮度向量中获得色度运动向量 */
  if( chor && cver ) {
    get_cvector( *cmvx, *cmvy, *ymvx, *ymvy, yhor, yver, chor, cver, info, t_level );
  }

  if( s_level > 0 ) //Never happen
    scale_down_mv( *ymvx, *ymvy, *cmvx, *cmvy, info, t_level );
}
