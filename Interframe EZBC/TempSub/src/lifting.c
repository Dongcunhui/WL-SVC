#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "basic.h"
#define EXTERN extern
#include "encoderN.h"
#include "bmeN.h"
#include "analsyn.h"
#include "structN.h"
#include "coderN.h"
#include "memoryN.h"
#include "util_filtering.h"
#include "iostream"
#include "fstream"
EXTERN float global_motion_active;


/*****************************************************************************/
/*                         determine_min_distance()                          */
/* determines minimal distance to nearest unconnected pixel   Hanke 03/03/12 */
/*****************************************************************************/

int 
determine_min_distance(unsigned char *pused0, unsigned char *pused1, int px,
                       int py, int hor, int ver, int SLTF_range, int *tpos)
{
  int tx, ty, t, xmin, ymin, xmax, ymax, ppos;
  
  *tpos = 0;
  ppos = py * hor + px;

  assert( pused0[ppos] == USED || pused1[ppos] == USED );
 
  for (t = 1; t <= SLTF_range; t++){
    // test boundaries
    xmin = (px - t <= 0) ? 0 : px - t;
    ymin = (py - t <= 0) ? 0 : py - t;
    xmax = (px + t >= hor - 1) ? hor - 1 : px + t;
    ymax = (py + t >= ver - 1) ? ver - 1 : py + t;
    
    for (ty = ymin; ty <= ymax; ty++){ 
      for (tx = xmin; tx <= xmax; tx++){
        *tpos = ty * hor + tx;
        if ((tx == xmin) || (tx == xmax) || (ty == ymin) || (ty == ymax)){
          if( pused0[*tpos] == UNUSED && pused1[*tpos] == UNUSED ){
            return t; // unconnected pixel found
          }
        }
      }
    }
  }
  
  // no match
  return (SLTF_range + 1);
}


/*****************************************************************************/
/*                        lowpass_transition_weight()                        */
/* weights connected pixels in the neighborhood of unconnected pixels (SLTF) */
/*****************************************************************************/

float 
lowpass_transition_weight(unsigned char *pused0, unsigned char *pused1, int px,
                          int py, int hor, int ver, int SLTF_range, int *tpos)
{
  int t = 0;
  float SLTF_weight;

  if( SLTF_range > 0 ){
    t = determine_min_distance (pused0, pused1, px, py, hor, ver, SLTF_range, tpos);
    assert (t > 0 && t <= SLTF_range + 1);
    SLTF_weight = (float)t / (float)(SLTF_range + 1);
  } else {
    SLTF_weight = 1.0;
  }
    
  return SLTF_weight;
}


/*****************************************************************************/
/*                              mc_analysis()                                */
/*                                                                           */
/*                               mv1  mv2   mv3                              */
/*                              <---  ---> <---                              */
/*                               |          |                                */
/*                      frp  fr0 | fr1   fr2| fr3                            */
/*                       . .   . |/ | \   | |/                               */
/*                       .  .  . /  |  \  | /                                */
/*                       .   . ./|  |   \ |/|                                */
/*                       .    H0 |  |    H1 |                                */
/*                       .    .\ |  |    /  |                                */
/*                       .   .  \|  |   /   |                                */
/*                       .  .    |  |  /    |                                */
/*                       . .     |\ | /     |                                */
/*                       ..      | \|/      |                                */
/*                       L0      |  L1      |                                */
/*                               |          |                                */
/*                                                                           */
/*  INPUT:  H0, fr1, fr2, fr3, mv1, mv2, mv3                                 */
/*  OUTPUT: L1, H1                                                           */
/*                                                                           */
/*****************************************************************************/
void
mc_analysis( float *L1, float *H1, float *H0, float *fr1, float *fr2, float *fr3,
             float *mvx1, float *mvy1, float *mvx2, float *mvy2, float *mvx3, float *mvy3,
             float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
             float *mvx3_int, float *mvy3_int, vector_ptr mv_ref1, vector_ptr mv_ref2,
             vector_ptr mv_ref3, int hor, int ver, int level, int remaining_frs,
             videoinfo info )
{
  int i;
  int cx, cy, px, py, cpos, ppos, case1, case2;
  float cfx, cfy, pfx0, pfy0, pfx1, pfy1;
  float *hweight, *lweight, *fweight, fweight_sum;
  float *ltmp0, *ltmp1, *ptmp0, *ptmp1, ptmp0_sum, ptmp1_sum;
  unsigned char *pused0, *pused1;

  //std::cout << level << ((mv_ref1==NULL)?0:1) <<" "<< ((mv_ref2==NULL)?0:1)<< " " << ((mv_ref3==NULL)?0:1)<< " " <<std::endl;

  int SLTF_range, SHTF_x_range, SHTF_y_range, SHTF_filter_size;
  int tx, ty, tpos, tcase, tx_min, ty_min, tx_max, ty_max;
  int fx, fy, fpos;
  float SLTF_weight = 1.0, *SHTF_filter;

  float sad=0;

  int case_biCon, case_leftCon, case_rightCon;
  int case_biPred, case_leftPred, case_rightPred, case_intra;
#ifdef COPYCOMPENSATION_WEIGHTING  
  float comp_weight;
#endif

  float max=0.0,min=0.0;

  pused0 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused0" );
  pused1 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused1" );
  ltmp0 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp0" );
  ltmp1 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp1" );


  for( i = 0; i < hor * ver; i++ ){
    pused0[i] = UNUSED;
    pused1[i] = UNUSED;
    ltmp0[i] = 0.0;
    ltmp1[i] = 0.0;
  }

  case_biCon = case_leftCon = case_rightCon = 0;
  case_biPred = case_leftPred = case_rightPred = case_intra = 0;

  SLTF_range   = info.SLTF_range; // 通过配置文件读入，默认为零
  SHTF_x_range = info.SHTF_range;
  SHTF_y_range = info.SHTF_range;
  //std::cout << SHTF_y_range << " " << SLTF_range << " " << SHTF_x_range << std::endl;

  assert( SLTF_range >= 0 && SLTF_range < 10 );
  assert( SHTF_x_range >= 0 && SHTF_x_range < 10 );
  assert( SHTF_y_range >= 0 && SHTF_y_range < 10 );
  
  // initialize SHTF_filter
  SHTF_filter_size = (2 * SHTF_x_range + 1) * (2 * SHTF_y_range + 1) ; // 默认算出为1

  SHTF_filter = ( float * )getarray( SHTF_filter_size, sizeof( float ), "SHTF_filter" );
  fweight = ( float * )getarray( SHTF_filter_size, sizeof( float ), "fweight" );
  ptmp0 = ( float * )getarray( SHTF_filter_size, sizeof( float ), "ptmp0" );
  ptmp1 = ( float * )getarray( SHTF_filter_size, sizeof( float ), "ptmp1" );
   
  for( fpos = 0; fpos < SHTF_filter_size; fpos++ ){
    SHTF_filter[fpos] = 1.0f / SHTF_filter_size;
  }

  /* generate temporal high subband H1 */
#ifdef COUT_MV_UPDATE
  int UVcom = 0;
  if (info.ywidth != hor) {
	  UVcom = 1;
  }
  char name[256], data_file_name[256], start[10];
  strncpy(name, info.bitname, strlen(info.bitname) - 4);
  name[strlen(info.bitname) - 4] = '\0';
  sprintf(data_file_name, "%s_enc_mv_update_test_yuv.txt", name);
  std::ofstream myfile(data_file_name, std::ios::app);
  if (myfile.is_open() == NULL)
  {
	  std::cout << "open file_y failed" << std::endl;
	  exit(0);
  }
  //if(UVcom==0)
	myfile << "new frame! It is UV " << UVcom << "\n";
  
#endif
  
  if( mv_ref2 != NULL || mv_ref3 != NULL ){ 
    for( cy = 0; cy < ver; cy++ ) {
      for( cx = 0; cx < hor; cx++ ) {  // 逐像素处理
        cpos = cy * hor + cx;  // 第几个像素点
        
#ifdef COPYCOMPENSATION_WEIGHTING 
       // assert(SHTF_filter_size == 1); // commented out by Yongjun Wu
        comp_weight = 1.0f;
#endif

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
		// when the pixel is in directional IBLOCK, continue
		if ( mv_ref2 != NULL && mv_ref3 != NULL )
		{
			if ( mvx2[cpos] == (float)HUGE_VAL  &&  mvy2[cpos] == (float)HUGE_VAL  &&
				 mvx3[cpos] == (float)HUGE_VAL  &&  mvy3[cpos] == (float)HUGE_VAL  &&
				 mvx2_int[cpos] == (float)HUGE_VAL && mvy2_int[cpos] == (float)HUGE_VAL &&
                 mvx3_int[cpos] == (float)HUGE_VAL && mvy3_int[cpos] == (float)HUGE_VAL )
			{
				 case_intra++;
				 continue;
			}
		} else if ( mv_ref2 != NULL && mv_ref3 == NULL)
		{
			if ( mvx2[cpos] == (float)HUGE_VAL  &&  mvy2[cpos] == (float)HUGE_VAL  &&
				 mvx2_int[cpos] == (float)HUGE_VAL && mvy2_int[cpos] == (float)HUGE_VAL )
			{
				 case_intra++;
				 continue;
			}
		} else if ( mv_ref2 == NULL && mv_ref3 != NULL )
		{
			if ( mvx3[cpos] == (float)HUGE_VAL  &&  mvy3[cpos] == (float)HUGE_VAL  &&
				 mvx3_int[cpos] == (float)HUGE_VAL && mvy3_int[cpos] == (float)HUGE_VAL )
			{
				 case_intra++;
				 continue;
			}
		}
#endif

        for( fpos = 0; fpos < SHTF_filter_size; fpos++ ){
          fweight[fpos] = 0.0;
          ptmp0[fpos] = 0.0;
          ptmp1[fpos] = 0.0;
        }
        fweight_sum = ptmp0_sum = ptmp1_sum = 0.0;
        pfx0 = pfy0 = pfx1 = pfy1 = 0.0;
        
        // get nearest inbound position  防止出界，但是默认情况下，不会出界
        tx_min = (cx - SHTF_x_range <= 0) ? 0 : cx - SHTF_x_range;
        ty_min = (cy - SHTF_y_range <= 0) ? 0 : cy - SHTF_y_range;
        tx_max = (cx + SHTF_x_range >= hor - 1) ? hor - 1 : cx + SHTF_x_range;
        ty_max = (cy + SHTF_y_range >= ver - 1) ? ver - 1 : cy + SHTF_y_range;  

        for (ty = ty_min, fy = 0; ty <= ty_max; ty++, fy++){ 
			for (tx = tx_min, fx = 0; tx <= tx_max; tx++, fx++){  //  SHTF_range大于0时才会多次执行，否则仅仅一次
            tpos = ty * hor + tx;                    // absolute position in MV / image arrays
            fpos = fy * (2 * SHTF_x_range + 1) + fx; // relative position in SHTF_filter array
            assert( fpos >= 0 && fpos < SHTF_filter_size );
#ifdef COUT_MV_UPDATE
			//if (UVcom == 0)
			{
				if (ty % 4 == 0 && tx % 4 == 0)
				{
					float left_mvx, left_mvy, right_mvx, right_mvy;
					// here the set of motion vectors are for Y component 
					left_mvx = (mvx2[tpos] == (float)HUGE_VAL) ? mvx2_int[tpos] : mvx2[tpos];
					left_mvy = (mvy2[tpos] == (float)HUGE_VAL) ? mvy2_int[tpos] : mvy2[tpos];
					right_mvx = (mvx3[tpos] == (float)HUGE_VAL) ? mvx3_int[tpos] : mvx3[tpos];
					right_mvy = (mvy3[tpos] == (float)HUGE_VAL) ? mvy3_int[tpos] : mvy3[tpos];
					if (left_mvx != (float)HUGE_VAL && left_mvy != (float)HUGE_VAL && right_mvx != (float)HUGE_VAL && right_mvy != (float)HUGE_VAL)  // 两个参考帧都可以参考
					{
						myfile << "0p " << " "<< ty/4 <<" "<<tx/4 << " " << left_mvx << " " << left_mvy << " " << right_mvx << " " << right_mvy << "\n";
						//myfile << "0p " << tpos << " " << left_mvx << " " << left_mvy << " " << right_mvx << " " << right_mvy << " ";
						//myfile << "0p " << tpos << " " << left_mvx << " " << left_mvy << " " << right_mvx << " " << right_mvy << "\n";

					}
					else if (right_mvx != (float)HUGE_VAL && right_mvy != (float)HUGE_VAL) // 右侧有参考帧
					{
						myfile << "1p " << " "<< ty/4 <<" "<<tx/4  << " " << right_mvx << " " << right_mvy << "\n";
						//myfile << "1p " << tpos << " " << right_mvx << " " << right_mvy << " ";
						//myfile << "1p " << tpos << " " << right_mvx << " " << right_mvy << "\n";
					}
					else if (left_mvx != (float)HUGE_VAL && left_mvy != (float)HUGE_VAL)
					{
						myfile << "2p " << " "<< ty/4 <<" "<<tx/4  << " " << left_mvx << " " << left_mvy << "\n";
						//myfile << "2p " << tpos << " " << left_mvx << " " << left_mvy << " ";
						//myfile << "2p " << tpos << " " << left_mvx << " " << left_mvy << "\n";
					}
				}
			}
			
#endif
            if( mvx2[tpos] != (float)HUGE_VAL && mvy2[tpos] != (float)HUGE_VAL &&
                mvx3[tpos] != (float)HUGE_VAL && mvy3[tpos] != (float)HUGE_VAL ){
              hweight = HPW1; // default mode (bi-directional)

              pfx0 = cx - mvx2[tpos];
              pfy0 = cy - mvy2[tpos];
              pfx1 = cx - mvx3[tpos];
              pfy1 = cy - mvy3[tpos];
              if( inbound( pfx0, pfy0, hor, ver ) && 
                  inbound( pfx1, pfy1, hor, ver ) ){
                ptmp0[fpos] = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
                ptmp1[fpos] = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
                assert( fweight[fpos] == 0.0 );
                fweight[fpos] = SHTF_filter[fpos]; // set filter postion
              }
              if ( tpos == cpos ) case_biCon++;
            }
            else if ( mvx3[tpos] != (float)HUGE_VAL && mvy3[tpos] != (float)HUGE_VAL ){
              hweight = HPW2; // forward mode
              pfx1 = cx - mvx3[tpos];
              pfy1 = cy - mvy3[tpos];
              if( inbound( pfx1, pfy1, hor, ver ) ){
                ptmp1[fpos] = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
                assert( fweight[fpos] == 0.0 );
                fweight[fpos] = SHTF_filter[fpos]; // set filter postion
              }
              if ( tpos == cpos ) case_rightCon++;
            }
            else if ( mvx2[tpos] != (float)HUGE_VAL && mvy2[tpos] != (float)HUGE_VAL ){
              hweight = HPW3; // backward mode
              pfx0 = cx - mvx2[tpos];
              pfy0 = cy - mvy2[tpos];
              if( inbound( pfx0, pfy0, hor, ver ) ){ 
                ptmp0[fpos] = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
                assert( fweight[fpos] == 0.0 );
                fweight[fpos] = SHTF_filter[fpos]; // set filter postion
              }
              if ( tpos == cpos ) case_leftCon++;
            }
            else{ 
              if( mvx2_int[tpos] != (float)HUGE_VAL && mvy2_int[tpos] != (float)HUGE_VAL &&
                  mvx3_int[tpos] != (float)HUGE_VAL && mvy3_int[tpos] != (float)HUGE_VAL ){
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW1_pred; // bi-predicted
                pfx0 = cx - mvx2_int[tpos];
                pfy0 = cy - mvy2_int[tpos];
                pfx1 = cx - mvx3_int[tpos];
                pfy1 = cy - mvy3_int[tpos];
                if( inbound( pfx0, pfy0, hor, ver ) && 
                    inbound( pfx1, pfy1, hor, ver ) ){
                  ptmp0[fpos] = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
                  ptmp1[fpos] = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
                  assert( fweight[fpos] == 0.0 );
                  fweight[fpos] = SHTF_filter[fpos]; // set filter postion
                }
                if ( tpos == cpos ) case_biPred++;
              }
              else if ( mvx3_int[tpos] != (float)HUGE_VAL && mvy3_int[tpos] != (float)HUGE_VAL ){
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW2_pred; // fwd-predicted
                pfx1 = cx - mvx3_int[tpos];
                pfy1 = cy - mvy3_int[tpos];
                if ( inbound( pfx1, pfy1, hor, ver ) ){
                  ptmp1[fpos] = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
                  assert( fweight[fpos] == 0.0 );
                  fweight[fpos] = SHTF_filter[fpos]; // set filter postion
                }
                if ( tpos == cpos ) case_rightPred++;
              }
              else if ( mvx2_int[tpos] != (float)HUGE_VAL && mvy2_int[tpos] != (float)HUGE_VAL ){
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW3_pred; // bwd-predicted
                pfx0 = cx - mvx2_int[tpos];
                pfy0 = cy - mvy2_int[tpos];
                if( inbound( pfx0, pfy0, hor, ver ) ){
                  ptmp0[fpos] = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
                  assert( fweight[fpos] == 0.0 );
                  fweight[fpos] = SHTF_filter[fpos]; // set filter postion
                }
                if ( tpos == cpos ) case_leftPred++;
              }
              else{ // intra mode (scene changes at both sides)
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW4_pred;
                
                if (tpos == cpos) case_intra++;
              }
            }
           
            assert( fweight[fpos] >= 0. && fweight[fpos] <= 1. ); 
            fweight_sum += fweight[fpos];
            ptmp0_sum += hweight[0] * fweight[fpos] * ptmp0[fpos];
            ptmp1_sum += hweight[2] * fweight[fpos] * ptmp1[fpos];

          }  // tx
        }   // ty

		////////////////////////////////////
		if(pfx0 < 0 || pfx0 >= hor)
			printf("error in pfx0 = %f\n",pfx0);
		if(pfx1 < 0 || pfx1 >= hor)
			printf("error in pfx1 = %f\n",pfx1);
		if(pfy0 < 0 || pfy0 >= ver)
			printf("error in pfy0 = %f\n",pfy0);
		if(pfy1 < 0 || pfy1 >= ver)
			printf("error in pfy1 = %f\n",pfy1);
		////////////////////////////////////

        
        // transition filtering
        ptmp0_sum = (fweight_sum == 0) ? 0.0f : (ptmp0_sum / fweight_sum);
        ptmp1_sum = (fweight_sum == 0) ? 0.0f : (ptmp1_sum / fweight_sum);

#ifdef COPYCOMPENSATION_WEIGHTING 
        H1[cpos] = comp_weight * ( ptmp0_sum + hweight[1] * fr2[cpos] + ptmp1_sum );
#ifdef COUT_MV_UPDATE
		//myfile << H1[cpos]<< " " << fr2[cpos]<< " "<<fr1[cpos]<<" "<<fr3[cpos] << "\n";
#endif
#else
        H1[cpos] = ( ptmp0_sum + hweight[1] * fr2[cpos] + ptmp1_sum );
#endif
//		if(cx >= 240 && cx < 256 && cy >= 192 && cy < 208 && hor == 352)
		if(hor == 352)
			sad += fabs(H1[cpos]);

//		printf("%f\t",fabs(H1[cpos]));

      }
//	  printf("\n");
    }

//	printf("max = %f, min = %f\n",max,min);

    if( SHTF_x_range > 0) 
		printf("SHTF (%d pixels)\n", SHTF_x_range ); 
        assert(case_biCon + case_leftCon + case_rightCon + case_biPred +
			   case_leftPred + case_rightPred + case_intra == hor * ver);
    
  }else{ // copy frame component
    for( cpos = 0; cpos < hor * ver; cpos++ ) {
#ifdef COPYCOMPENSATION_WEIGHTING 
      H1[cpos] = copycomp_weight_high[level] * HPW4[1] * fr2[cpos];
#else
      H1[cpos] = HPW4[1] * fr2[cpos];
#endif
    }
    case_intra = hor * ver;
  }

//  if(hor == 352)
//	  printf("sad = %f\n",sad);
  
  /* generate temporal low subband L1 */
  
  // step 1 - motion compensation of H0 and H1
  for( cy = 0; cy < ver; cy++ ) {
    for( cx = 0; cx < hor; cx++ ) { // 遍历H0
      cpos = cy * hor + cx;
      
      if ( mvx1[cpos] != (float)HUGE_VAL && mvy1[cpos] != (float)HUGE_VAL ){
        pfx0 = cx - mvx1[cpos]; //低频帧中对应的位置
        pfy0 = cy - mvy1[cpos];//低频帧中对应的位置
        position( &px, &py, pfx0, pfy0, mvx1[cpos], mvy1[cpos], hor, ver ); // 返回整数位置
        ppos = py * hor + px;
        assert((ppos >= 0) && ppos < (hor * ver));
        cfx = px + mvx1[cpos];
        cfy = py + mvy1[cpos];
        
        if( pused0[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ) {
          ltmp0[ppos] = interpolate( cfx, cfy, H0, hor, ver, TYPE );
          pused0[ppos] = USED;
#ifdef COUT_MV_UPDATE
		  //if (UVcom == 0)
		  {
			  //if (py % 4 == 0 && px % 4 == 0)   // 4倍下采样
			  //{
				 // float left_mvx, left_mvy, right_mvx, right_mvy;
				 // // here the set of motion vectors are for Y component 
				 // left_mvx = mvx1[cpos];
				 // left_mvy = mvy1[cpos];
				 // myfile << "1u " << py/4 <<" "<< px/4 << " " << left_mvx << " " << left_mvy << "\n";
				 // //myfile << "1u " <<ppos<< " " << left_mvx << " " << left_mvy << "\n";
			  //}

			  {  // 不下采样
				  float left_mvx, left_mvy, right_mvx, right_mvy;
				  // here the set of motion vectors are for Y component 
				  left_mvx = mvx1[cpos];
				  left_mvy = mvy1[cpos];
				  myfile << "1u " << py << " " << px << " " << left_mvx << " " << left_mvy << "\n";
				  //myfile << "1u " <<ppos<< " " << left_mvx << " " << left_mvy << "\n";
			  }
		  }
#endif
        }
      } 
    }
  }


  for( cy = 0; cy < ver; cy++ ) {
    for( cx = 0; cx < hor; cx++ ) { // 遍历H1
      cpos = cy * hor + cx;
      
      if ( mvx2[cpos] != (float)HUGE_VAL && mvy2[cpos] != (float)HUGE_VAL ){
        pfx1 = cx - mvx2[cpos];
        pfy1 = cy - mvy2[cpos];
        position( &px, &py, pfx1, pfy1, mvx2[cpos], mvy2[cpos], hor, ver );
        ppos = py * hor + px;  // 当前图像的整数位置
        cfx = px + mvx2[cpos]; // 参考图像的位置
        cfy = py + mvy2[cpos]; // 
        
        if( pused1[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ) {
          ltmp1[ppos] = interpolate( cfx, cfy, H1, hor, ver, TYPE );// 将参考图像的位置的像素插值后赋值到当前图像的位置
          pused1[ppos] = USED;
#ifdef COUT_MV_UPDATE
		  //if (UVcom == 0)
		  {
			  //if (py % 4 == 0 && px % 4 == 0) {// 4倍下采样
				 // float right_mvx, right_mvy;
				 // // here the set of motion vectors are for Y component 
				 // right_mvx = mvx2[cpos];
				 // right_mvy = mvy2[cpos];
				 // myfile << "2u " << py/4 <<" "<<px/4 << " " << right_mvx << " " << right_mvy << "\n";
			  //}

			  {  // 不下采样
				  float right_mvx, right_mvy;
				  // here the set of motion vectors are for Y component 
				  right_mvx = mvx2[cpos];
				  right_mvy = mvy2[cpos];
				  myfile << "2u " << py << " " << px << " " << right_mvx << " " << right_mvy << "\n";
			  }
		  }
#endif
        }
      }
    }
  }

  case1 = 0; 
  case2 = 0;
  tcase = 0;
  
  // step 2 - calculation of L1
  if( mv_ref1 != NULL || mv_ref2 != NULL ){
    for( py = 0; py < ver; py++ ) {
      for( px = 0; px < hor; px++ ) {
        ppos = py * hor + px;
        
        if( pused0[ppos] == USED && pused1[ppos] == USED ){
          lweight = LPW1; // default mode (bi-directional)
          case1++;
        }
        else if( pused0[ppos] == UNUSED && pused1[ppos] == USED ){
          lweight = LPW2; // forward mode
          case1++;
        }
        else if( pused0[ppos] == USED && pused1[ppos] == UNUSED ){
          lweight = LPW3; // backward mode 
          case1++;
        }
        else{
          lweight = LPW4; // unconnected pixels
          case2++;
        }
   
        // lowpass transition filtering
        if( pused0[ppos] == USED || pused1[ppos] == USED ){
          SLTF_weight = lowpass_transition_weight(pused0, pused1, px, py, hor, ver, SLTF_range, &tpos);
          assert (SLTF_weight > 0.0 && SLTF_weight <= 1.0);
          if( SLTF_weight != 1.0 ) tcase++;
        }
        
        assert( ppos >= 0 && ppos < hor * ver );
        L1[ppos] = SLTF_weight * ( lweight[0] * ltmp0[ppos] + lweight[2] * ltmp1[ppos] ) + lweight[1] * fr1[ppos];    
#ifdef COUT_MV_UPDATE
		//myfile <<"u "<<ppos<<" "<< L1[ppos] << " " << fr1[ppos]<<" "<<H0[ppos]<<" "<<H1[ppos] << "\n";
#endif
      }
    }
  
    if( SLTF_range > 0) 
      printf("SLTF (%d pixels) -- connected %d, transitions %d, unconnected %d\n", SLTF_range, case1-tcase, tcase, case2 ); 

    /* checking */
    if (case1 + case2 != hor * ver) {
      printf( "error in lifting.c -- connected: %d, unconnected: %d\n", case1, case2 );
      exit( 1 );
    }
  }
  else{ // copy frame component
    for( ppos = 0; ppos < hor * ver; ppos++ ) {
#ifdef COPYCOMPENSATION_WEIGHTING 
      L1[ppos] = copycomp_weight_low[level] * LPW4[1] * fr1[ppos];    
#else
      L1[ppos] = LPW4[1] * fr1[ppos];    
#endif
    }
  }
#ifdef COUT_MV_UPDATE
  myfile.close();
#endif  

  free( SHTF_filter );
  free( fweight );
  free( ptmp0 );
  free( ptmp1 );

  free( pused0 );
  free( pused1 );
  free( ltmp0 );
  free( ltmp1 );
    
}


// scale down the motion vectors by 2 for U V components 
void
scale_down_mv_by_scale( float *mvx, float *mvy, videoinfo info, int t_level, int shift )
{
  float dx, dy;
  int   t, precision, scale;

  if ( *mvx==(float)HUGE_VAL  && *mvy==(float)HUGE_VAL )
	  return; 

  scale = 1<<shift; 
  precision = 1<<info.subpel[t_level]; 

  dx = *mvx;
  assert(dx!=(float)HUGE_VAL); 
  dx /= scale;
  t = ( int )dx;
  dx = ( float )( ( int )( dx * precision ) % precision ) / ( float )precision;
  *mvx = t + dx;
  dy = *mvy; 
  assert(dy!=(float)HUGE_VAL);
  dy /= scale; 
  t = ( int )dy;
  dy = ( float )( ( int )( dy * precision ) % precision ) / ( float )precision;
  *mvy = t + dy;
}


// get interpolated value for a neighbor motion vector 分局邻居mv获得插值 为了OBMC
// pos and imagemeinfo are always in Y coordinate
// pos: neighbor position 
// fr1 and fr3 may be Y, may be U V componet, which is indicated by UVcom
// cx, cy, hor, ver may be in Y coordinate, may be in U V coordinate
void   get_value_for_neighbor_mv(int cpos, int pos, ImageMEinfo *imagemeinfo, float *fr1, float *fr3, 
								 float *neighbor_value, float *self_weight, float neighbor_weight,
								 int cx, int cy, videoinfo info, int hor, int ver, 
								 int  UVcom, int t_level)
{
//	Pay attention, this ver/hor has been halved for U/V components, but both cpos and pos belong to the full-size/luminance frame.
//	cpos - position of the current MC pixel
//  pos  - pixel of neighbor block; 
	float pfx0, pfy0, pfx1, pfy1; 
	float ptmp0, ptmp1;
	float *hweight;
	float left_mvx, left_mvy, right_mvx, right_mvy; 
	int   s_level; 

	int nx, ny, sx, sy; //x,y of neighbor/self blocks
	float dx1, dy1, dx2, dy2;

	// here the set of motion vectors are for Y component with full resolution 下面的if用来获得左右侧的mv，用来进行下面的邻居补偿
	if( 1 ){
		ny = pos / (hor<<UVcom);
		nx = pos - ny * (hor<<UVcom); // 当前像素在当前帧空间邻居的坐标
		sy = cpos / (hor<<UVcom);  
		sx = cpos - sy * (hor<<UVcom); // 当前像素的坐标

		if( !(nx == sx || ny == sy) )  //至少有一个相等
			printf("nx = %d, ny = %d, sx = %d, sy = %d\n",nx,ny,sx,sy);
		assert(nx == sx || ny == sy);

		if(nx == sx){ //AFFINE in Y dimension  说明是上下邻居
			if( imagemeinfo[pos].left_mvx != (float)HUGE_VAL ){//有左侧x mv
				assert( imagemeinfo[pos].left_mvy != (float)HUGE_VAL);// 断定一定有左侧 y的mv
				dx1 =  imagemeinfo[pos + (hor<<UVcom)].left_mvx - imagemeinfo[pos].left_mvx;
				dy1 =  imagemeinfo[pos + (hor<<UVcom)].left_mvy - imagemeinfo[pos].left_mvy;
				left_mvx = imagemeinfo[pos].left_mvx + dx1 * (sy - ny);
				left_mvy = imagemeinfo[pos].left_mvy + dy1 * (sy - ny);
			}else{
				assert( imagemeinfo[pos].left_mvx == (float)HUGE_VAL && imagemeinfo[pos].left_mvy == (float)HUGE_VAL );
				left_mvx = (float)HUGE_VAL;
				left_mvy = (float)HUGE_VAL;
			}

			if( imagemeinfo[pos].right_mvx != (float)HUGE_VAL ){
				assert( imagemeinfo[pos].right_mvy != (float)HUGE_VAL);
				dx2 =  imagemeinfo[pos + (hor<<UVcom)].right_mvx - imagemeinfo[pos].right_mvx;
				dy2 =  imagemeinfo[pos + (hor<<UVcom)].right_mvy - imagemeinfo[pos].right_mvy;
				right_mvx = imagemeinfo[pos].right_mvx + dx2 * (sy - ny);
				right_mvy = imagemeinfo[pos].right_mvy + dy2 * (sy - ny);
			}else{
				assert( imagemeinfo[pos].right_mvx == (float)HUGE_VAL && imagemeinfo[pos].right_mvy == (float)HUGE_VAL );
				right_mvx = (float)HUGE_VAL;
				right_mvy = (float)HUGE_VAL;
			}
		}
		else{ //AFFINE in X dimension  说明是左右邻居
			assert(ny == sy);
			if( imagemeinfo[pos].left_mvx != (float)HUGE_VAL ){
				assert( imagemeinfo[pos].left_mvy != (float)HUGE_VAL);
				dx1 =  imagemeinfo[pos + 1].left_mvx - imagemeinfo[pos].left_mvx;
				dy1 =  imagemeinfo[pos + 1].left_mvy - imagemeinfo[pos].left_mvy;
				left_mvx = imagemeinfo[pos].left_mvx + dx1 * (sx - nx);
				left_mvy = imagemeinfo[pos].left_mvy + dy1 * (sx - nx);
			}
			else{
				assert( imagemeinfo[pos].left_mvx == (float)HUGE_VAL && imagemeinfo[pos].left_mvy == (float)HUGE_VAL );
				left_mvx = (float)HUGE_VAL;
				left_mvy = (float)HUGE_VAL;
			}

			if( imagemeinfo[pos].right_mvx != (float)HUGE_VAL ){
				assert( imagemeinfo[pos].right_mvy != (float)HUGE_VAL);
				dx2 =  imagemeinfo[pos + 1].right_mvx - imagemeinfo[pos].right_mvx;
				dy2 =  imagemeinfo[pos + 1].right_mvy - imagemeinfo[pos].right_mvy;
				right_mvx = imagemeinfo[pos].right_mvx + dx2 * (sx - nx);
				right_mvy = imagemeinfo[pos].right_mvy + dy2 * (sx - nx);
			}else{
				assert( imagemeinfo[pos].right_mvx == (float)HUGE_VAL && imagemeinfo[pos].right_mvy == (float)HUGE_VAL );
				right_mvx = (float)HUGE_VAL;
				right_mvy = (float)HUGE_VAL;
			}
		}
	}
	else{
		assert(0);
		assert(imagemeinfo[pos].bi_mode <= 8);

		left_mvx   = imagemeinfo[pos].left_mvx;
		left_mvy   = imagemeinfo[pos].left_mvy;
		right_mvx  = imagemeinfo[pos].right_mvx;
		right_mvy  = imagemeinfo[pos].right_mvy; 
	}
	

	// scale down the neighbor motion vectors by corresponding resolution reduction 
	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
	if (s_level>0)
	{
		scale_down_mv_by_scale( &left_mvx,  &left_mvy, info, t_level,  s_level );
		scale_down_mv_by_scale( &right_mvx, &right_mvy, info, t_level, s_level );
	}

	ptmp0 = ptmp1 = (float)0.0; 
	if( left_mvx != (float)HUGE_VAL && 
		left_mvy != (float)HUGE_VAL &&
        right_mvx!= (float)HUGE_VAL && 
		right_mvy!= (float)HUGE_VAL) // 双向预测
	{         
		hweight = HPW1; // default mode (bi-directional)
		if (UVcom)
		{
			scale_down_mv_by_scale( &left_mvx, &left_mvy, info,   t_level,  1 );
			scale_down_mv_by_scale( &right_mvx, &right_mvy, info, t_level,  1 );
		}
        pfx0 = cx - left_mvx; // 在参考帧的位置
        pfy0 = cy - left_mvy;
        pfx1 = cx - right_mvx;
        pfy1 = cy - right_mvy;
        if( inbound( pfx0, pfy0, hor, ver ) && 
            inbound( pfx1, pfy1, hor, ver ) ){
			ptmp0 = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
            ptmp1 = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
			*neighbor_value += (ptmp0*hweight[0]+ptmp1*hweight[2])*neighbor_weight;
        }
		else  // raise self weight to compensate for invalid neighbor motion vector 提高自身的权重来补偿邻居的权重
			*self_weight += neighbor_weight; 
    } 
	else if (  right_mvx!= (float)HUGE_VAL && 
				 right_mvy!= (float)HUGE_VAL ) // 前向预测
	{
              hweight = HPW2; // forward mode
			  if (UVcom)
				scale_down_mv_by_scale( &right_mvx, &right_mvy, info, t_level, 1 );
              pfx1 = cx - right_mvx;
              pfy1 = cy - right_mvy;
              if( inbound( pfx1, pfy1, hor, ver ) ){
                ptmp1 = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
                *neighbor_value += ptmp1*hweight[2]*neighbor_weight;
              }else
				  *self_weight += neighbor_weight; 
     }else if ( left_mvx != (float)HUGE_VAL && 
				left_mvy != (float)HUGE_VAL){
              hweight = HPW3; // backward mode
			  if (UVcom)
				scale_down_mv_by_scale( &left_mvx, &left_mvy, info, t_level, 1 );
              pfx0 = cx - left_mvx;
              pfy0 = cy - left_mvy;
              if( inbound( pfx0, pfy0, hor, ver ) ){ 
                ptmp0 = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
				*neighbor_value += ptmp0*hweight[0]*neighbor_weight;
              }else
				  *self_weight += neighbor_weight; 
     }else
		 assert(0); 
}


// 通过自身mv获得预测值
void get_value_for_self_mv(int pos, ImageMEinfo *imagemeinfo, float *fr1, float *fr3, 
								 float *self_value, float self_weight, 
								 int cx, int cy, videoinfo info, int hor, int ver,
								 int *case_bi, int *case_left, int *case_right,
								int UVcom, int t_level)
{
	float *hweight; 
	float ptmp0, ptmp1;
	float pfx0, pfy0, pfx1, pfy1; 
	float left_mvx, left_mvy, right_mvx, right_mvy; 
	int   s_level; 

	// here the set of motion vectors are for Y component 
	left_mvx   = imagemeinfo[pos].left_mvx; // 左侧mv的x信息
	left_mvy   = imagemeinfo[pos].left_mvy;
	right_mvx  = imagemeinfo[pos].right_mvx;
	right_mvy  = imagemeinfo[pos].right_mvy; 
	// scale down the self motion vector by corresponding resolution reduction 
	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
	if (s_level>0)
	{
		scale_down_mv_by_scale( &left_mvx,  &left_mvy,  info, t_level,  s_level );
		scale_down_mv_by_scale( &right_mvx, &right_mvy, info, t_level,  s_level );
	}
	ptmp0 = ptmp1 = (float)0.0; 
	if( left_mvx != (float)HUGE_VAL && 
		left_mvy != (float)HUGE_VAL &&
        right_mvx!= (float)HUGE_VAL && 
		right_mvy!= (float)HUGE_VAL)  // 两个参考帧都可以参考
	{              
		hweight = HPW1; // default mode (bi-directional)
		if (UVcom)
		{
			scale_down_mv_by_scale( &left_mvx,  &left_mvy, info, t_level,  1 );
			scale_down_mv_by_scale( &right_mvx, &right_mvy, info, t_level, 1 );
		}
        pfx0 = cx - left_mvx;  // 获取参考像素的位置
        pfy0 = cy - left_mvy;
        pfx1 = cx - right_mvx;
        pfy1 = cy - right_mvy;  
		if ( ! ( inbound( pfx0, pfy0, hor, ver ) && 
                inbound( pfx1, pfy1, hor, ver ) ) )
				printf("Error here\n");
        assert( inbound( pfx0, pfy0, hor, ver ) && 
                inbound( pfx1, pfy1, hor, ver ) );
		ptmp0 = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
        ptmp1 = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
//    块级写出
		// information for current block in Y coordinate with full resolution 
//		if (!UVcom)
//		{
//			int leftx = imagemeinfo[pos].leftx;
//			int topy = imagemeinfo[pos].topy;  // 这个块上边界在图像中的位置
//			int xblk = imagemeinfo[pos].blksize;
//			int yblk = imagemeinfo[pos].blksize;
//			//printf("%d, %d, %d, %d\n", cx, cy, topy, leftx);
//			if (leftx == cx && topy == cy)
//			{
//				//std::ofstream myfile("interpolate fream.txt", std::ios::app);
//				for (int block_h = cy; block_h < topy + yblk; block_h++)
//				{
//					for (int block_w = cx; block_w < leftx + xblk; block_w++)
//					{
//						int count = (block_w - cx) + (block_h - cy) * yblk;
//						int cnn_pos = block_h * hor + block_w;
//						//printf("%d, %d, %d\n", cnn_pos, block_h, block_w);
//						float cnn_pfx0 = block_w - left_mvx;  // 获取参考像素的位置
//						float cnn_pfy0 = block_h - left_mvy;
//						float cnn_pfx1 = block_w - right_mvx;
//						float cnn_pfy1 = block_h - right_mvy;
//						if (!(inbound(cnn_pfx0, cnn_pfy0, hor, ver) &&
//								inbound(cnn_pfx1, cnn_pfy1, hor, ver)))
//								continue;
//						//if (!(inbound(cnn_pfx0, cnn_pfy0, hor, ver) &&
//						//	inbound(cnn_pfx1, cnn_pfy1, hor, ver)))
//						//	printf("Error here\n");
//						//assert(inbound(cnn_pfx0, cnn_pfy0, hor, ver) &&
//						//	inbound(cnn_pfx1, cnn_pfy1, hor, ver));
//						imagemeinfo[cnn_pos].cnn_ptmp0 = interpolate(cnn_pfx0, cnn_pfy0, fr1, hor, ver, TYPE);
//						imagemeinfo[cnn_pos].cnn_ptmp1 = interpolate(cnn_pfx1, cnn_pfy1, fr3, hor, ver, TYPE);
//						//myfile << cnn_pos << " " << imagemeinfo[cnn_pos].cnn_ptmp0 << " " << imagemeinfo[cnn_pos].cnn_ptmp1 << "\n";
//					}
//				}
//				//myfile.close();
//			}

		*self_value = (ptmp0*hweight[0] + ptmp1 * hweight[2])*self_weight; // 通过mv自身权重获得的补偿值
		(*case_bi)++;
#ifdef CNN_TEMPORAL_WAVELET
		/*******写出*****/
			//std::cout<< "mv*** "<<UVcom <<" 0 "<< pos << " " << ptmp0 << " " << ptmp1 <<"\n";
		//char name[256], data_file_name[256], start[10];
		//strncpy(name, info.bitname, strlen(info.bitname) - 9);
		//name[strlen(info.bitname) - 9] = '\0';
		//sprintf(data_file_name, "%s_enc_mv_org_yuv.txt", name);
		//std::ofstream myfile(data_file_name, std::ios::app);
		//if (myfile.is_open() == NULL)
		//{
		//	std::cout << "open file_y failed" << std::endl;
		//	exit(0);
		//}
		//myfile << "0" << " " << left_mvx << " " << left_mvy << " " << right_mvx << " " << right_mvy << " "<< pos << " " << ptmp0 << " " << ptmp1 <<" ";
		//myfile.close();

			//// 解码端
			//char name[256], data_file_name[256], start[10];
			//strncpy(name, info.bitname, strlen(info.bitname) - 14);
			//name[strlen(info.bitname) - 14] = '\0';
			//sprintf(data_file_name, "%s_dec_5_3_predict.txt", name);
			//std::ofstream myfile(data_file_name, std::ios::app);
			//if (myfile.is_open() == NULL)
			//{
			//	std::cout << "open file_y failed" << std::endl;
			//	exit(0);
			//}
			//myfile << pos << " " << ptmp0 << " " << ptmp1 << " ";
			//myfile.close();
		/*******写出*****/
#endif


    } 
	else if ( right_mvx!= (float)HUGE_VAL && 
		        right_mvy!= (float)HUGE_VAL){ // 右侧有参考帧
              hweight = HPW2; // forward mode
 			  if (UVcom)
				  scale_down_mv_by_scale( &right_mvx, &right_mvy, info, t_level, 1 );
              pfx1 = cx - right_mvx;
              pfy1 = cy - right_mvy;
              assert( inbound( pfx1, pfy1, hor, ver ) ); 
              ptmp1 = interpolate( pfx1, pfy1, fr3, hor, ver, TYPE );
              *self_value = ptmp1*hweight[2]*self_weight;
			  (*case_right)++; 
#ifdef CNN_TEMPORAL_WAVELET
			  /*******写出*****/
				  //std::cout << "mv*** "<<UVcom << " 1 " << pos << " " << ptmp1 << "\n";
			  //char name[256], data_file_name[256], start[10];
			  //strncpy(name, info.bitname, strlen(info.bitname) - 9);
			  //name[strlen(info.bitname) - 9] = '\0';
			  //sprintf(data_file_name, "%s_enc_mv_org_yuv.txt", name);
			  //std::ofstream myfile(data_file_name, std::ios::app);
			  //if (myfile.is_open() == NULL)
			  //{
				 // std::cout << "open file_y failed" << std::endl;
				 // exit(0);
			  //}
			  //myfile << "1" << " " << right_mvx << " " << right_mvy << " "<< pos << " " << ptmp1 << " ";
			  //myfile.close();

				  //// 解码端
				  //char name[256], data_file_name[256], start[10];
				  //strncpy(name, info.bitname, strlen(info.bitname) - 14);
				  //name[strlen(info.bitname) - 14] = '\0';
				  //sprintf(data_file_name, "%s_dec_5_3_predict.txt", name);
				  //std::ofstream myfile(data_file_name, std::ios::app);
				  //if (myfile.is_open() == NULL)
				  //{
				  //	std::cout << "open file_y failed" << std::endl;
				  //	exit(0);
				  //}
				  //myfile << pos << " " << ptmp0 << " " << ptmp1 << " ";
				  //myfile.close();
			  /*******写出*****/
#endif
     }
	else if ( left_mvx != (float)HUGE_VAL && 
				left_mvy != (float)HUGE_VAL ){
              hweight = HPW3; // backward mode
			  if (UVcom)
				  scale_down_mv_by_scale( &left_mvx,  &left_mvy, info, t_level, 1 );
              pfx0 = cx - left_mvx;
              pfy0 = cy - left_mvy;
              assert( inbound( pfx0, pfy0, hor, ver ) ); 
              ptmp0 = interpolate( pfx0, pfy0, fr1, hor, ver, TYPE );
			  *self_value = ptmp0*hweight[0]*self_weight;
			  (*case_left)++; 
#ifdef CNN_TEMPORAL_WAVELET
			  /*******写出*****/
				  //std::cout << "mv*** "<<UVcom << " 2 " << pos << " " << ptmp0 << "\n";
				  //char name[256], data_file_name[256], start[10];
				  //strncpy(name, info.bitname, strlen(info.bitname) - 9);
				  //name[strlen(info.bitname) - 9] = '\0';
				  //sprintf(data_file_name, "%s_enc_mv_org_yuv.txt", name);
				  //std::ofstream myfile(data_file_name, std::ios::app);
				  //if (myfile.is_open() == NULL)
				  //{
					 // std::cout << "open file_y failed" << std::endl;
					 // exit(0);
				  //}
				  //myfile <<"2"<<" "<< left_mvx<< " " <<left_mvy <<" "<< pos << " " << ptmp0 <<" ";
				  //myfile.close();

				  //// 解码端
				  //char name[256], data_file_name[256], start[10];
				  //strncpy(name, info.bitname, strlen(info.bitname) - 14);
				  //name[strlen(info.bitname) - 14] = '\0';
				  //sprintf(data_file_name, "%s_dec_5_3_predict.txt", name);
				  //std::ofstream myfile(data_file_name, std::ios::app);
				  //if (myfile.is_open() == NULL)
				  //{
				  //	std::cout << "open file_y failed" << std::endl;
				  //	exit(0);
				  //}
				  //myfile << pos << " " << ptmp0 << " " << ptmp1 << " ";
				  //myfile.close();

			  /*******写出*****/
#endif
     }
}


/*****************************************************************************/
/*                              mc_analysis()                                */
/*                                                                           */
/*                               mv1  mv2   mv3                              */
/*                              <---  ---> <---                              */
/*                               |          |                                */
/*                      frp  fr0 | fr1   fr2| fr3                            */
/*                       . .   . |/ | \   | |/                               */
/*                       .  .  . /  |  \  | /                                */
/*                       .   . ./|  |   \ |/|                                */
/*                       .    H0 |  |    H1 |                                */
/*                       .    .\ |  |    /  |                                */
/*                       .   .  \|  |   /   |                                */
/*                       .  .    |  |  /    |                                */
/*                       . .     |\ | /     |                                */
/*                       ..      | \|/      |                                */
/*                       L0      |  L1      |                                */
/*                               |          |                                */
/*                                                                           */
/*  INPUT:  H0, fr1, fr2, fr3, mv1, mv2, mv3                                 */
/*  OUTPUT: L1, H1                                                           */
/*               输出高频帧和低频帧                                                            */
/*****************************************************************************/

int count = 0;

void
mc_analysis_with_OBMC( float *L1, float *H1, float *H0, float *fr1, float *fr2, float *fr3,
             float *mvx1, float *mvy1, float *mvx2, float *mvy2, float *mvx3, float *mvy3,
             float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
             float *mvx3_int, float *mvy3_int, vector_ptr mv_ref1, vector_ptr mv_ref2,
             vector_ptr mv_ref3, int hor, int ver, int level, int remaining_frs,
             videoinfo info, ImageMEinfo *imagemeinfo, Varblkarrayinfo *varblkarray, int UVcom,
			 enum FLAG left_scene, enum FLAG right_scene, int type)
{
  int i;
  int   cx, cy, px, py, cpos, com_pos, ppos, case1, case2;
  float cfx, cfy, pfx0, pfy0, pfx1, pfy1;
  float *lweight;
  float *ltmp0, *ltmp1;
  unsigned char *pused0, *pused1;
  int    leftx, topy, xblk, yblk, disx, disy, col, row, pos;
  float  self_weight, ver_weight, hor_weight;
  float  self_value, neighbor_value; 
  int case_bi, case_left, case_right, case_intra; // 用于记录采用这四种预测模式的像素的个数
  float sad = 0;

  float getnum, diff;
  int pixel_count;
  FILE *pc;

  // UVcom: indication for U V component
  // for U V components we use the sub-sampled weighting coefficients 

#ifdef COPYCOMPENSATION_WEIGHTING  
  float comp_weight;
#endif
  //float comp_weight = 1;
  // the variables for update step in MCTF 用于update的变量 整帧处理
  pused0 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused0" );
  pused1 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused1" );
  ltmp0 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp0" );
  ltmp1 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp1" );
  for( i = 0; i < hor * ver; i++ )
  {
    pused0[i] = UNUSED;
    pused1[i] = UNUSED;
    ltmp0[i] = 0.0;
    ltmp1[i] = 0.0;
  }
#ifdef CNN_TEMPORAL_WAVELET
  char name[256], data_file_name[256], start[10];
  strncpy(name, info.bitname, strlen(info.bitname) - 9);
  name[strlen(info.bitname) - 9] = '\0';
  sprintf(data_file_name, "%s_enc_mv_org_yuv.txt", name);
  std::ofstream myfile(data_file_name, std::ios::app);
  if (myfile.is_open() == NULL)
  {
	  std::cout << "open file_y failed" << std::endl;
	  exit(0);
  }
	myfile << "new frame! It is UV "<<UVcom <<"\n";

#endif
  case_bi = case_left = case_right = case_intra= 0;

  /* generate temporal high subband H1 生成高频子带 */
  if( left_scene == NO || right_scene == NO )// 两边存在可以参考
  { 
    for( cy = 0; cy < ver; cy++ ) 
	{
      for( cx = 0; cx < hor; cx++ )  // 像素级滤波，为什么不进行块级滤波
	  { 
        cpos = (cy<<UVcom) * (hor<<UVcom) + (cx<<UVcom);  // the position in Y coordinate  当前位置在Y通道上的第几个像素点
		com_pos = cy*hor+cx;   // the position in this component, either Y or U V  当前位置在当前通道上的第几个像素点
#ifdef CNN_TEMPORAL_WAVELET

		float left_mvx, left_mvy, right_mvx, right_mvy;
		// here the set of motion vectors are for Y component 
		left_mvx = imagemeinfo[cpos].left_mvx; // 左侧mv的x信息
		left_mvy = imagemeinfo[cpos].left_mvy;
		right_mvx = imagemeinfo[cpos].right_mvx;
		right_mvy = imagemeinfo[cpos].right_mvy;
		if (left_mvx != (float)HUGE_VAL && left_mvy != (float)HUGE_VAL && right_mvx != (float)HUGE_VAL && right_mvy != (float)HUGE_VAL)  // 两个参考帧都可以参考
		{
			myfile << "0" << " " << left_mvx << " " << left_mvy << " " << right_mvx << " " << right_mvy << " ";

		}
		else if (right_mvx != (float)HUGE_VAL && right_mvy != (float)HUGE_VAL) // 右侧有参考帧
		{
			myfile << "1" << " " << right_mvx << " " << right_mvy << " ";
		}
		else if (left_mvx != (float)HUGE_VAL && left_mvy != (float)HUGE_VAL)
		{
			myfile << "2" << " " << left_mvx << " " << left_mvy << " ";
		}
		//myfile << cpos << "\n";

#endif

#ifdef COPYCOMPENSATION_WEIGHTING 
        comp_weight = 1.0f;
#endif

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
		// the pixel is in directional IBLOCK  in Y coordinate 
		if (imagemeinfo[cpos].bi_mode == DIRECTIONAL_IBLOCK )
		{
			case_intra++;// 使用帧内模式
			continue; 
		}
#endif

		assert(imagemeinfo[cpos].bi_mode != UNDEFINED);

		// information for current block in Y coordinate with full resolution 
		leftx     = imagemeinfo[cpos].leftx;   
		topy      = imagemeinfo[cpos].topy;  // 这个块上边界在图像中的位置
		xblk      = imagemeinfo[cpos].blksize; 
		yblk      = imagemeinfo[cpos].blksize;

//		if(UVcom == 0)
//			printf("cx = %d, cy = %d, topy = %d, leftx = %d\n",cx,cy,topy,leftx);

		// the distance from this pixel to the corner in Y coordinate 当前像素到当前块的左、上边沿的距离
		disx = (cx<<UVcom) - leftx;   
		disy = (cy<<UVcom) - topy;
		// weight coefficients for self motion vector and neighbor motion vectors 
		self_weight = imagemeinfo[cpos].self_weight; 
		hor_weight  = imagemeinfo[cpos].h_weight; 
		ver_weight  = imagemeinfo[cpos].v_weight;
		assert(self_weight+hor_weight+ver_weight==(float)1.0);
        self_value = neighbor_value = 0.0;
		// current position ( cx<<UVcom, cy<<UVcom) in Y coordinate  当前像素点在Y通道的第几行和第几列
        col = cx<<UVcom;   row = cy<<UVcom; 
#ifdef NO_OBMC
		get_value_for_self_mv(cpos, imagemeinfo, fr1, fr3, &self_value, 1.0,
			cx, cy, info, hor, ver, &case_bi, &case_left, &case_right,
			UVcom, level);  // 在这里面去获取相邻两帧的预测值
#else
		// contributions from vertical neighbor motion vector  来自竖直邻居mv的贡献
		if (disy < (yblk >> 1))  {    //  top neighbor is effective 当前像素到哪边更近，这里是到上面更近
			pos = (topy-3)*(hor<<UVcom) + col;  // 上面一个块的倒数第三行  邻居的位置
			if ( (topy > 0) && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK )
				get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1, fr3, &neighbor_value, 
					                      &self_weight, ver_weight, cx, cy, info, hor, ver, 
										  UVcom, level);
		} 
		else if (disy >= (yblk >> 1))  { // bottom neighbor is effective 
			pos = (topy+yblk+2)*(hor<<UVcom) + col; // 下面一个块的第二行
			if  ( (topy+yblk < (ver<<UVcom) )  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK ) 
				get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1, fr3, &neighbor_value, 
					                      &self_weight, ver_weight, cx, cy, info, hor, ver,
										  UVcom, level);
		}

        // contribution from horizontal neighbor motion vectors 来自水平邻居mv的贡献
		if (disx < (xblk >> 1))  { // left neighbor is effective 
			pos = row*(hor<<UVcom) + leftx-3; 
			if ( (leftx > 0)  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK ) 
				 get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1, fr3, &neighbor_value, 
					                       &self_weight, hor_weight, cx, cy, info, hor, ver,
										   UVcom, level);
		}
		else if (disx >= (xblk >> 1))  { // right neighbor is effective 
			pos = row*(hor<<UVcom) + leftx+xblk+2; 
			if ( (leftx+xblk < (hor<<UVcom) )  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK )
				 get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1, fr3, &neighbor_value, 
					                       &self_weight, hor_weight, cx, cy, info, hor, ver,
										   UVcom, level);
		}
		get_value_for_self_mv(cpos, imagemeinfo, fr1, fr3, &self_value, self_weight, 
							   cx, cy, info, hor, ver, &case_bi, &case_left, &case_right,
							   UVcom, level);  // 在这里面去获取相邻两帧的预测值
#endif
		///////////////////////////////////////////////////
		//MARK BLK MRG
/*		if( type == 1 && level == 0 ){
			if(imagemeinfo[cpos].bi_mode >= 9 || (imagemeinfo[cpos].bi_mode == 7 && ) )
				fr2[com_pos] = fr2[com_pos] / 0.4;
		}

		if(UVcom == 0 && level == 0){
			if(imagemeinfo[cpos].mrg_blk_bdr == 1)
				fr2[cpos] = 255;
		}
*/		////////////////////////////////////////////////////
/*		//MARK TWO COMP
		if( type == 1 && level == 0 ){
			if(imagemeinfo[cpos].two_comp_src == 1)
				fr2[com_pos] = fr2[com_pos] / 0.4;
		}

		if(UVcom == 0 && level == 0){
			if(imagemeinfo[cpos].two_comp_bdr == 1)
				fr2[cpos] = 255;
		}
*/		////////////////////////////////////////////////////
		if(imagemeinfo[cpos].skip_sign == YES){
			H1[com_pos] = 0;
		}else
			H1[com_pos] = comp_weight * ( HPW1[1] * fr2[com_pos] + self_value + neighbor_value );// 计算出来H1
#ifdef CNN_TEMPORAL_WAVELET
			//myfile << cpos <<"\n";
			myfile << com_pos << " " << H1[com_pos] << " " << fr1[com_pos] << " " << fr2[com_pos] << " " << fr3[com_pos] << "\n";
#endif
	  }
	}

}
  else  // if( left_scene == NO || right_scene == NO )
  { // copy frame component
    for( cpos = 0; cpos < hor * ver; cpos++ ) { // 逐个像素去做，使用的是帧内模式
#ifdef COPYCOMPENSATION_WEIGHTING 
      H1[cpos] = copycomp_weight_high[level] * HPW4[1] * fr2[cpos];
#else
      H1[cpos] = HPW4[1] * fr2[cpos];
#endif
    }
    case_intra = hor * ver;
  }
#ifdef CNN_TEMPORAL_WAVELET
  myfile.close();
#endif
  assert(case_intra + case_bi + case_left + case_right == hor*ver ); 
  /* generate temporal low subband L1 生成低频帧*/
  // step 1 - motion compensation of H0 and H1 运动补偿
  for( cy = 0; cy < ver; cy++ ) {
    for( cx = 0; cx < hor; cx++ ) {
      cpos = cy * hor + cx; // 当前位置
      if ( mvx1[cpos] != (float)HUGE_VAL && mvy1[cpos] != (float)HUGE_VAL ) // 是否接受左边高频帧的更新
	  {
        pfx0 = cx - mvx1[cpos];
        pfy0 = cy - mvy1[cpos]; // 在参考图像中的位置
        position( &px, &py, pfx0, pfy0, mvx1[cpos], mvy1[cpos], hor, ver );
        ppos = py * hor + px; // 取整之后的位置
        assert((ppos >= 0) && ppos < (hor * ver));
        cfx = px + mvx1[cpos];
        cfy = py + mvy1[cpos]; // 当前位置的浮点位置
        if( pused0[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ) {
          ltmp0[ppos] = interpolate( cfx, cfy, H0, hor, ver, TYPE );
          pused0[ppos] = USED;
        }
      } 
    }
  }


  for( cy = 0; cy < ver; cy++ ) {  // 右侧高频帧更新
    for( cx = 0; cx < hor; cx++ ) {
      cpos = cy * hor + cx;
      if ( mvx2[cpos] != (float)HUGE_VAL && mvy2[cpos] != (float)HUGE_VAL ){
        pfx1 = cx - mvx2[cpos];
        pfy1 = cy - mvy2[cpos];
        position( &px, &py, pfx1, pfy1, mvx2[cpos], mvy2[cpos], hor, ver );
        ppos = py * hor + px;
        cfx = px + mvx2[cpos]; 
        cfy = py + mvy2[cpos];
        if( pused1[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ) {
          ltmp1[ppos] = interpolate( cfx, cfy, H1, hor, ver, TYPE );
          pused1[ppos] = USED;
        }
      }
    }
  }

  
  case1 = 0; 
  case2 = 0;
  // step 2 - calculation of L1
  //mv_ref1 = NULL;
  //mv_ref2 = NULL;

  if (mv_ref1 != NULL || mv_ref2 != NULL) // 有参考
	  //if (0) // 有参考
  {
    for( py = 0; py < ver; py++ ) {
		for (px = 0; px < hor; px++) {
			ppos = py * hor + px;
        
        if( pused0[ppos] == USED && pused1[ppos] == USED ){
          lweight = LPW1; // default mode (bi-directional)
          case1++;
        }
        else if( pused0[ppos] == UNUSED && pused1[ppos] == USED ){
          lweight = LPW2; // forward mode
          case1++;
        }
        else if( pused0[ppos] == USED && pused1[ppos] == UNUSED ){
          lweight = LPW3; // backward mode 
          case1++;
        }
        else{
          lweight = LPW4; // unconnected pixels
          case2++;
        }
        assert( ppos >= 0 && ppos < hor * ver );
#ifdef NO_UPDATE
		//lweight[1] = LPW4[1] * copycomp_weight_low[level];
		//lweight[0] = 0;
		//lweight[2] = 0;
		//printf("%f, %f, %f\n", lweight[0], lweight[1], lweight[2]);
		//L1[ppos] = LPW4[1] * copycomp_weight_low[level] * fr1[ppos];    // 计算出L1
		lweight = LPW4;
		L1[ppos] = (lweight[0] * ltmp0[ppos] + lweight[2] * ltmp1[ppos]) + lweight[1] * fr1[ppos];    // 计算出L1
#else
        L1[ppos] = ( lweight[0] * ltmp0[ppos] + lweight[2] * ltmp1[ppos] ) + lweight[1] * fr1[ppos];    // 计算出L1
#endif

		/**************写出**********/
		//{
		//	char name[256], data_file_name[256], start[10];
		//	strncpy(name, info.bitname, strlen(info.bitname) - 9);
		//	name[strlen(info.bitname) - 9] = '\0';
		//	sprintf(data_file_name, "%s_L.txt", name);
		//	std::ofstream myfile(data_file_name, std::ios::app);
		//	if (myfile.is_open() == NULL)
		//	{
		//		std::cout << "open file_y failed" << std::endl;
		//		exit(0);
		//	}
		//	myfile << ppos << " " << L1[ppos] << "\n";
		//	myfile.close();
		//}
		/**************写出**********/
      }
    }
/*
	if(UVcom == 0 && level == 0){
		diff = 0;	pixel_count = 0;
//		printf("MV matrix:\n\n");
		for( cy = 0; cy < ver; cy++ ){
		  for( cx = 0; cx < hor; cx++ ){
			cpos = ( cy<<UVcom) * ( hor<<UVcom) + (cx<<UVcom);  
			getnum = (float)HUGE_VAL;

			if(imagemeinfo[cpos].left_mvx != (float)HUGE_VAL && imagemeinfo[cpos].right_mvx != (float)HUGE_VAL ){
				getnum = (imagemeinfo[cpos].left_mvx - imagemeinfo[cpos].right_mvx)/2;
			}else if(imagemeinfo[cpos].left_mvx != (float)HUGE_VAL){
				getnum = imagemeinfo[cpos].left_mvx;
			}else if(imagemeinfo[cpos].right_mvx != (float)HUGE_VAL){
				getnum = (-1)*imagemeinfo[cpos].right_mvx;
			}
			
			if( getnum != (float)HUGE_VAL ){
				if(buff_frameMEinfo[cpos].left_mvx != (float)HUGE_VAL && buff_frameMEinfo[cpos].right_mvx != (float)HUGE_VAL ){
					getnum = getnum - (buff_frameMEinfo[cpos].left_mvx - buff_frameMEinfo[cpos].right_mvx)/2;
				}else if(buff_frameMEinfo[cpos].left_mvx != (float)HUGE_VAL){
					getnum = getnum - buff_frameMEinfo[cpos].left_mvx;
				}else if(buff_frameMEinfo[cpos].right_mvx != (float)HUGE_VAL){
					getnum = getnum - (-1)*buff_frameMEinfo[cpos].right_mvx;
				}
			}

			if( getnum != (float)HUGE_VAL ){
				diff = diff + fabs(getnum);
				pixel_count ++;
			}

//copy
			buff_frameMEinfo[cpos].left_mvx = imagemeinfo[cpos].left_mvx;
			buff_frameMEinfo[cpos].left_mvy = imagemeinfo[cpos].left_mvy;
			buff_frameMEinfo[cpos].right_mvx = imagemeinfo[cpos].right_mvx;
			buff_frameMEinfo[cpos].right_mvy = imagemeinfo[cpos].right_mvy;

		  }
		}

		diff = diff / pixel_count;
//		printf("AVG diff = %f, pixel_count = %d\n",diff,pixel_count);
//		printf("\n\nMV matrix END:\n");

		assert( diff >= 0 );
		if(diff > global_motion_active)
			global_motion_active = diff;

		pc = fopen("enc_block.txt","at");
		fprintf(pc,"%f\n",diff);
		fclose(pc);
	}//if UVcom = 0
*/  
    // checking 
    if (case1 + case2 != hor * ver) {
      printf( "error in lifting.c -- connected: %d, unconnected: %d\n", case1, case2 );
      exit( 1 );
    }
  }
  else{ // copy frame component
    for( ppos = 0; ppos < hor * ver; ppos++ ) {
#ifdef COPYCOMPENSATION_WEIGHTING 
      L1[ppos] = copycomp_weight_low[level] * LPW4[1] * fr1[ppos];    // 计算出L1
#else
      L1[ppos] = LPW4[1] * fr1[ppos];    
#endif
    }
  }

  free( pused0 );
  free( pused1 );
  free( ltmp0 );
  free( ltmp1 );
    
}

/*****************************************************************************/
/*                              mc_analysis()                                */
/*                                                                           */
/*                               mv1  mv2   mv3                              */
/*                              <---  ---> <---                              */
/*                               |          |                                */
/*                      frp  fr0 | fr1   fr2| fr3                            */
/*                       . .   . |/ | \   | |/                               */
/*                       .  .  . /  |  \  | /                                */
/*                       .   . ./|  |   \ |/|                                */
/*                       .    H0 |  |    H1 |                                */
/*                       .    .\ |  |    /  |                                */
/*                       .   .  \|  |   /   |                                */
/*                       .  .    |  |  /    |                                */
/*                       . .     |\ | /     |                                */
/*                       ..      | \|/      |                                */
/*                       L0      |  L1      |                                */
/*                               |          |                                */
/*                                                                           */
/*  INPUT:  H0, fr1, fr2, fr3, mv1, mv2, mv3                                 */
/*  OUTPUT: L1, H1                                                           */
/*                                                                           */
/*****************************************************************************/


/*****************************************************************************/
/*                               mc_synthesis()                              */
/*                                                                           */
/*                        mv0   mv1   mv2                                    */
/*                        ---> <---   --->                                   */
/*                         |           |           |                         */
/*                     frp | fr0   fr1 | fr2   fr3 |                         */
/*                      . \|  |   / | \|  .   . . .|                         */
/*                      .  \  |  /  |  \  .  .  .  .                         */
/*                      .  |\ | /   |  |\ . .   .  |.                        */
/*                      .  |  H0    |  | H1     .  | H2                      */
/*                      .  |        |  |        .  |                         */
/*                      .  |        |  |        .  |                         */
/*                      .  |        |  |        .  |                         */
/*                      L0 |        L1 |        L2 |                         */
/*                         |           |           |                         */
/*                                                                           */
/*  INPUT:  H0, L1, H1, frp, mv0, mv1, mv2                                   */
/*  OUTPUT: fr0, fr1                                                         */
/*                                                                           */
/*****************************************************************************/ 
void
mc_synthesis( float *fr0, float *fr1, float *H0, float *L1, float *H1, float *frp,
              float *mvx0, float *mvy0, float *mvx1, float *mvy1, 
              float *mvx2, float *mvy2, float *mvx0_int, float *mvy0_int, 
              float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
              vector_ptr mv_ref0, vector_ptr mv_ref1, vector_ptr mv_ref2,
              int hor, int ver, int level, videoinfo info )
{
  int i;
  int cx, cy, px, py, cpos, ppos, case1, case2, case3;
  float cfx, cfy, pfx0, pfy0, pfx1, pfy1;
  float *hweight, *lweight, *fweight, fweight_sum;
  float *ltmp0, *ltmp1, *ptmp0, *ptmp1, ptmp0_sum, ptmp1_sum;
  unsigned char *pused0, *pused1;

  int SLTF_range, SHTF_x_range, SHTF_y_range, SHTF_filter_size;
  int tx, ty, tpos, tcase, tx_min, ty_min, tx_max, ty_max;
  int fx, fy, fpos;
  float SLTF_weight = 1.0, *SHTF_filter;
#ifdef COPYCOMPENSATION_WEIGHTING 
  float comp_weight;
#endif

  // step 0 - calculation of fr1 (A) / scene_change propagation
  if (level < info.t_level) 
  {
    if( mv_ref1 != NULL || mv_ref2 != NULL )
	{ 
      for( ppos = 0; ppos < hor * ver; ppos++ )
	  {
        fr1[ppos] = L1[ppos] / LPW4[1];   // 第一次合成
      }
    } 
	else 
	{
      for( ppos = 0; ppos < hor * ver; ppos++ )
	  {
#ifdef COPYCOMPENSATION_WEIGHTING
        fr1[ppos] = L1[ppos] / (LPW4[1] * copycomp_weight_low[level]);
#else
        fr1[ppos] = L1[ppos];
#endif
      }
    }
    return;
  }


  pused0 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused0" );
  pused1 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused1" );
  ltmp0 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp0" );
  ltmp1 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp1" );
  
  for( i = 0; i < hor * ver; i++ ){
    pused0[i] = UNUSED;
    pused1[i] = UNUSED;
    ltmp0[i] = 0.0;
    ltmp1[i] = 0.0;
  }
  
  SLTF_range   = info.SLTF_range;
  SHTF_x_range = info.SHTF_range;
  SHTF_y_range = info.SHTF_range;

  assert( SLTF_range >= 0 && SLTF_range < 10 );
  assert( SHTF_x_range >= 0 && SHTF_x_range < 10 );
  assert( SHTF_y_range >= 0 && SHTF_y_range < 10 );
  
  // initialize SHTF_filter
  SHTF_filter_size = (2 * SHTF_x_range + 1) * (2 * SHTF_y_range + 1) ;

  SHTF_filter = ( float * )getarray( SHTF_filter_size, sizeof( float ), "SHTF_filter" );
  fweight = ( float * )getarray( SHTF_filter_size, sizeof( float ), "fweight" );
  ptmp0 = ( float * )getarray( SHTF_filter_size, sizeof( float ), "ptmp0" );
  ptmp1 = ( float * )getarray( SHTF_filter_size, sizeof( float ), "ptmp1" );
  
  for( fpos = 0; fpos < SHTF_filter_size; fpos++ ){
    SHTF_filter[fpos] = 1.0f / SHTF_filter_size; 
  }
#ifdef COUT_MV_UPDATE

  int UVcom = 0;
  if (info.ywidth != hor) {
	  UVcom = 1;
  }
  bool new_file = 0;
  char name[256], data_file_name[256], start[10];
  strncpy(name, info.bitname, strlen(info.bitname) - 4);
  name[strlen(info.bitname) - 4] = '\0';
  sprintf(data_file_name, "%s_dec_mv_update_test_yuv.txt", name);
  
  { // 判断是否为新创建的文件
	std::ifstream fin(data_file_name);
	if (!fin) 
		new_file = 1;
	else
		fin.close();
  }
  
  std::ofstream myfile(data_file_name, std::ios::app);
  if (myfile.is_open() == NULL)
  {
	  std::cout << "open file_y failed" << std::endl;
	  exit(0);
  }
  if (new_file)
  {
	  myfile << "frame number: " << info.last - info.start + 1 << "\n";
	  myfile << "width height: " << info.ywidth <<" "<< info.yheight << "\n";
  }
  myfile << "new frame! It is UV " << UVcom << "\n";
#endif
  /* reconstruction of fr1 (A) */

  // step 1 - motion compensation of H0 and H1
  // MC of H0
  for( cy = 0; cy < ver; cy++ ){
    for( cx = 0; cx < hor; cx++ ){
      cpos = cy * hor + cx;
      
      if( mvx1[cpos] != (float)HUGE_VAL && mvy1[cpos] != (float)HUGE_VAL ) {
        pfx0 = cx - mvx1[cpos];
        pfy0 = cy - mvy1[cpos];
        position( &px, &py, pfx0, pfy0, mvx1[cpos], mvy1[cpos], hor, ver );
        ppos = py * hor + px;
        cfx = px + mvx1[cpos];
        cfy = py + mvy1[cpos];
        
        if( pused0[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ){
          pused0[ppos] = USED;
          ltmp0[ppos] = interpolate( cfx, cfy, H0, hor, ver, TYPE );
#ifdef COUT_MV_UPDATE
		  {
			  //if (py % 4 == 0 && px % 4 == 0) { // 4倍下采样
				 // float left_mvx, left_mvy, right_mvx, right_mvy;
				 // // here the set of motion vectors are for Y component 
				 // left_mvx = mvx1[cpos];
				 // left_mvy = mvy1[cpos];
				 // myfile << "1u " << py / 4 << " " << px / 4 << " " << left_mvx << " " << left_mvy << "\n";
				 // //myfile << "1u " <<ppos<< " " << left_mvx << " " << left_mvy << "\n";
			  //}
			  { // 不下采样
				  float left_mvx, left_mvy, right_mvx, right_mvy;
				  // here the set of motion vectors are for Y component 
				  left_mvx = mvx1[cpos];
				  left_mvy = mvy1[cpos];
				  myfile << "1u " << py << " " << px << " " << left_mvx << " " << left_mvy << "\n";
				  //myfile << "1u " <<ppos<< " " << left_mvx << " " << left_mvy << "\n";
			  }
		  }
#endif
        }
      }
    }
  }
  
  // MC of H1
  for( cy = 0; cy < ver; cy++ ){
    for( cx = 0; cx < hor; cx++ ){
      cpos = cy * hor + cx;
      
      if( mvx2[cpos] != (float)HUGE_VAL && mvy2[cpos] != (float)HUGE_VAL ) {
        pfx1 = cx - mvx2[cpos];
        pfy1 = cy - mvy2[cpos];
        
        position( &px, &py, pfx1, pfy1, mvx2[cpos], mvy2[cpos], hor, ver );
        ppos = py * hor + px;
        cfx = px + mvx2[cpos];
        cfy = py + mvy2[cpos];
        
        if( pused1[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ) {
          pused1[ppos] = USED;
          ltmp1[ppos] = interpolate( cfx, cfy, H1, hor, ver, TYPE );
#ifdef COUT_MV_UPDATE
		  {
			  //if (py % 4 == 0 && px % 4 == 0) { 4倍下采样
				 // float right_mvx, right_mvy;
				 // // here the set of motion vectors are for Y component 
				 // right_mvx = mvx2[cpos];
				 // right_mvy = mvy2[cpos];
				 // myfile << "2u " << py / 4 << " " << px / 4 << " " << right_mvx << " " << right_mvy << "\n";
			  //}

			  { // 不下采样
				  float right_mvx, right_mvy;
				  // here the set of motion vectors are for Y component 
				  right_mvx = mvx2[cpos];
				  right_mvy = mvy2[cpos];
				  myfile << "2u " << py << " " << px << " " << right_mvx << " " << right_mvy << "\n";
			  }
		  }
#endif
        }
      }
    }
  }

  case1 = case2 = tcase = 0;
  double sum_update = 0;
  double sum_update_resi = 0;
  double sum_update_L = 0;
  // step 2 - calculation of fr1 (A)
  if( mv_ref1 != NULL || mv_ref2 != NULL ){ 
    for( py = 0; py < ver; py++ ) {
      for( px = 0; px < hor; px++ ) {
        ppos = py * hor + px;
        
        if( pused0[ppos] == USED && pused1[ppos] == USED ){
          lweight = LPW1; // default mode (bi-directional)
          case1++;
        }
        else if( pused0[ppos] == UNUSED && pused1[ppos] == USED ){
          lweight = LPW2; // forward mode
          case1++;
        }
        else if( pused0[ppos] == USED && pused1[ppos] == UNUSED ){
          lweight = LPW3; // backward mode 
          case1++;
        }
        else{
          lweight = LPW4; // intra mode (unconnected pixels)
          case2++;
        }
   
        // lowpass transition filtering
        if( pused0[ppos] == USED || pused1[ppos] == USED ){
          SLTF_weight = lowpass_transition_weight(pused0, pused1, px, py, hor, ver, SLTF_range, &tpos);
          assert (SLTF_weight > 0.0 && SLTF_weight <= 1.0);
          if( SLTF_weight != 1.0 ) tcase++;
        }
        
        assert( ppos >= 0 && ppos < hor * ver );
        fr1[ppos] = 1 / lweight[1] * ( L1[ppos] - SLTF_weight * ( lweight[0] * ltmp0[ppos] + lweight[2] * ltmp1[ppos] )); 
		//std::cout << lweight[1] <<" "<<SLTF_weight<< std::endl;
		sum_update += (fr1[ppos]);
		sum_update_resi += ( (lweight[0] * ltmp0[ppos] + lweight[2] * ltmp1[ppos]) /lweight[1]);
		sum_update_L += ( double( L1[ppos])/ lweight[1]);
      }
    }
	//if (hor == 832)
	//{
	//	std::cout.precision(10);
	//	std::cout << "update " << sum_update << std::endl;
	//	std::cout << "update_resi " << sum_update_resi << std::endl;
	//	std::cout << "update_L " << sum_update_L << std::endl;
	//}

   
    if( SLTF_range > 0) 
      printf("SLTF (%d pixels) -- connected %d, transitions %d, unconnected %d\n", SLTF_range, case1-tcase, tcase, case2 ); 

    /* checking */
    if( case1 + case2 != hor * ver ){
      printf( "error in mctf() covered: %d, uncovered: %d\n", case1, case2 );
      exit( 1 );
    }
  }
  else{ // copy frame component
    for( ppos = 0; ppos < hor * ver; ppos++ ){
      // fr1[ppos] = 1 / LPW4[1] * L1[ppos];
#ifdef COPYCOMPENSATION_WEIGHTING
      fr1[ppos] = L1[ppos] / (LPW4[1] * copycomp_weight_low[level]);
#else
      fr1[ppos] = L1[ppos];
#endif
    }
  }

  case1 = case2 = case3 = 0;

  /* reconstruction of fr0 (B) */
  //double sum_predict = 0;
  //double sum_predict_H = 0;
  //double sum_predict_pre = 0;
  if( mv_ref0 != NULL || mv_ref1 != NULL ){ 
    for( cy = 0; cy < ver; cy++ ){
      for( cx = 0; cx < hor; cx++ ){
        cpos = cy * hor + cx;
        
#ifdef COPYCOMPENSATION_WEIGHTING 
       // assert(SHTF_filter_size == 1);  // commented by Yongjun Wu
        comp_weight = 1.0f;
#endif

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
		// when the pixel is in directional IBLOCK, continue
		if ( mv_ref0 != NULL && mv_ref1 != NULL )
		{
			if ( mvx0[cpos] == (float)HUGE_VAL  &&  mvy0[cpos] == (float)HUGE_VAL  &&
				 mvx1[cpos] == (float)HUGE_VAL  &&  mvy1[cpos] == (float)HUGE_VAL  &&
				 mvx0_int[cpos] == (float)HUGE_VAL && mvy0_int[cpos] == (float)HUGE_VAL &&
                 mvx1_int[cpos] == (float)HUGE_VAL && mvy1_int[cpos] == (float)HUGE_VAL )
			{
				case3++;
				continue;
			}
		} else if ( mv_ref0 != NULL && mv_ref1 == NULL)
		{
			if ( mvx0[cpos] == (float)HUGE_VAL  &&  mvy0[cpos] == (float)HUGE_VAL  &&
				 mvx0_int[cpos] == (float)HUGE_VAL && mvy0_int[cpos] == (float)HUGE_VAL )
			{
				case3++;
				continue;
			}

		} else if ( mv_ref0 == NULL && mv_ref1 != NULL )
		{
			if ( mvx1[cpos] == (float)HUGE_VAL  &&  mvy1[cpos] == (float)HUGE_VAL  &&
				 mvx1_int[cpos] == (float)HUGE_VAL && mvy1_int[cpos] == (float)HUGE_VAL )
			{
				case3++;
				continue;
			}
		}
#endif

        for( fpos = 0; fpos < SHTF_filter_size; fpos++ ){
          fweight[fpos] = 0.0;
          ptmp0[fpos] = 0.0;
          ptmp1[fpos] = 0.0;
        }
           
        // get nearest inbound position
        tx_min = (cx - SHTF_x_range <= 0) ? 0 : cx - SHTF_x_range;
        ty_min = (cy - SHTF_y_range <= 0) ? 0 : cy - SHTF_y_range;
        tx_max = (cx + SHTF_x_range >= hor - 1) ? hor - 1 : cx + SHTF_x_range;
        ty_max = (cy + SHTF_y_range >= ver - 1) ? ver - 1 : cy + SHTF_y_range;  

        fweight_sum = ptmp0_sum = ptmp1_sum = 0.0;
        pfx0 = pfy0 = pfx1 = pfy1 = 0.0;
        
        for (ty = ty_min, fy = 0; ty <= ty_max; ty++, fy++){ 
          for (tx = tx_min, fx = 0; tx <= tx_max; tx++, fx++){
            
            tpos = ty * hor + tx;                    // absolute position in MV / image arrays
            fpos = fy * (2 * SHTF_x_range + 1) + fx; // relative position in SHTF_filter array
            assert( fpos >= 0 && fpos < SHTF_filter_size );
#ifdef COUT_MV_UPDATE
			if (ty % 4 == 0 && tx % 4 == 0)
			{
				float left_mvx, left_mvy, right_mvx, right_mvy;
				// here the set of motion vectors are for Y component 
				left_mvx = (mvx0[tpos] == (float)HUGE_VAL) ? mvx0_int[tpos] : mvx0[tpos];
				left_mvy = (mvy0[tpos] == (float)HUGE_VAL) ? mvy0_int[tpos] : mvy0[tpos];
				right_mvx = (mvx1[tpos] == (float)HUGE_VAL) ? mvx1_int[tpos] : mvx1[tpos];
				right_mvy = (mvy1[tpos] == (float)HUGE_VAL) ? mvy1_int[tpos] : mvy1[tpos];
				if (left_mvx != (float)HUGE_VAL && left_mvy != (float)HUGE_VAL && right_mvx != (float)HUGE_VAL && right_mvy != (float)HUGE_VAL)  // 两个参考帧都可以参考
				{
					myfile << "0p " << " "<< ty/4 <<" "<<tx/4 << " " << left_mvx << " " << left_mvy << " " << right_mvx << " " << right_mvy << "\n";
					//myfile << "0p" << " " << left_mvx << " " << left_mvy << " " << right_mvx << " " << right_mvy << "\n";

				}
				else if (right_mvx != (float)HUGE_VAL && right_mvy != (float)HUGE_VAL) // 右侧有参考帧
				{
					myfile << "1p " << " "<< ty/4 <<" "<<tx/4  << " " << right_mvx << " " << right_mvy << "\n";
					//myfile << "1p" << " " << right_mvx << " " << right_mvy << "\n";
				}
				else if (left_mvx != (float)HUGE_VAL && left_mvy != (float)HUGE_VAL)
				{
					myfile << "2p " << " "<< ty/4 <<" "<<tx/4  << " " << left_mvx << " " << left_mvy << "\n";
					//myfile << "2p" << " " << left_mvx << " " << left_mvy << "\n";
				}
			}
#endif
            if( mvx0[tpos] != (float)HUGE_VAL && mvy0[tpos] != (float)HUGE_VAL &&
                mvx1[tpos] != (float)HUGE_VAL && mvy1[tpos] != (float)HUGE_VAL ){ 
              hweight = HPW1; // default mode (bi-directional)
              pfx0 = cx - mvx0[tpos];
              pfy0 = cy - mvy0[tpos];
              pfx1 = cx - mvx1[tpos];
              pfy1 = cy - mvy1[tpos];
              if( inbound( pfx0, pfy0, hor, ver ) && 
                  inbound( pfx1, pfy1, hor, ver ) ){
                ptmp0[fpos] = interpolate( pfx0, pfy0, frp, hor, ver, TYPE );
                ptmp1[fpos] = interpolate( pfx1, pfy1, fr1, hor, ver, TYPE );
                assert( fweight[fpos] == 0.0 );
                fweight[fpos] = SHTF_filter[fpos]; // set filter postion
              }
              if ( tpos == cpos ) case1++;
            }
            else if( mvx1[tpos] != (float)HUGE_VAL && mvy1[tpos] != (float)HUGE_VAL ){ 
              hweight = HPW2; // forward mode
              pfx1 = cx - mvx1[tpos];
              pfy1 = cy - mvy1[tpos];
              if( inbound( pfx1, pfy1, hor, ver ) ){
                ptmp1[fpos] = interpolate( pfx1, pfy1, fr1, hor, ver, TYPE );
                assert( fweight[fpos] == 0.0 );
                fweight[fpos] = SHTF_filter[fpos]; // set filter postion
              }
              if ( tpos == cpos ) case1++;
            }
            else if( mvx0[tpos] != (float)HUGE_VAL && mvy0[tpos] != (float)HUGE_VAL ){ 
              hweight = HPW3; // backward mode
              pfx0 = cx - mvx0[tpos];
              pfy0 = cy - mvy0[tpos];
              if( inbound( pfx0, pfy0, hor, ver ) ){
                ptmp0[fpos] = interpolate( pfx0, pfy0, frp, hor, ver, TYPE ); 
                assert( fweight[fpos] == 0.0 );
                fweight[fpos] = SHTF_filter[fpos]; // set filter postion
              }
              if ( tpos == cpos ) case1++;
            } 
            else{ // predicted modes
              if ( tpos == cpos ) case2++;
              
              if( mvx0_int[tpos] != (float)HUGE_VAL && mvy0_int[tpos] != (float)HUGE_VAL &&
                  mvx1_int[tpos] != (float)HUGE_VAL && mvy1_int[tpos] != (float)HUGE_VAL ){ 
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW1_pred; // bi-directional
                pfx0 = cx - mvx0_int[tpos];
                pfy0 = cy - mvy0_int[tpos];
                pfx1 = cx - mvx1_int[tpos];
                pfy1 = cy - mvy1_int[tpos];
                if( inbound( pfx0, pfy0, hor, ver ) && 
                    inbound( pfx1, pfy1, hor, ver ) ){
                  ptmp0[fpos] = interpolate( pfx0, pfy0, frp, hor, ver, TYPE );
                  ptmp1[fpos] = interpolate( pfx1, pfy1, fr1, hor, ver, TYPE );
                  assert( fweight[fpos] == 0.0 );
                  fweight[fpos] = SHTF_filter[fpos]; // set filter postion
                }
              }
              else if( mvx1_int[tpos] != (float)HUGE_VAL && mvy1_int[tpos] != (float)HUGE_VAL ){ 
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW2_pred; // predicted - forward mode
                pfx1 = cx - mvx1_int[tpos];
                pfy1 = cy - mvy1_int[tpos];
                if( inbound( pfx1, pfy1, hor, ver ) ){
                  ptmp1[fpos] = interpolate( pfx1, pfy1, fr1, hor, ver, TYPE );
                  assert( fweight[fpos] == 0.0 );
                  fweight[fpos] = SHTF_filter[fpos]; // set filter postion
                }
              }
              else if( mvx0_int[tpos] != (float)HUGE_VAL && mvy0_int[tpos] != (float)HUGE_VAL ){ 
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW3_pred; // predicted - backward mode
                pfx0 = cx - mvx0_int[tpos];
                pfy0 = cy - mvy0_int[tpos]; 
                if( inbound( pfx0, pfy0, hor, ver ) ){
                  ptmp0[fpos] = interpolate( pfx0, pfy0, frp, hor, ver, TYPE );
                  assert( fweight[fpos] == 0.0 );
                  fweight[fpos] = SHTF_filter[fpos]; // set filter postion
                } 
              }
              else{ // intra - intra mode (scene changes at both sides) 
#ifdef COPYCOMPENSATION_WEIGHTING_FORPREDICTED
                comp_weight = copycomp_weight_high[level];
#endif
                hweight = HPW4_pred;

              }
            }
            
            assert( fweight[fpos] >= 0. && fweight[fpos] <= 1. );  
            fweight_sum += fweight[fpos];
            ptmp0_sum += hweight[0] * fweight[fpos] * ptmp0[fpos];
            ptmp1_sum += hweight[2] * fweight[fpos] * ptmp1[fpos];      
          
          } // tx
        }  // ty
        
        // transition filtering
        ptmp0_sum = (fweight_sum == 0) ? 0.0f : (ptmp0_sum / fweight_sum);
        ptmp1_sum = (fweight_sum == 0) ? 0.0f : (ptmp1_sum / fweight_sum);

//		printf("%f\t",fabs(H0[cpos]));
         
#ifdef COPYCOMPENSATION_WEIGHTING 
        fr0[cpos] = 1 / hweight[1] * ( - ptmp0_sum + H0[cpos] / comp_weight - ptmp1_sum );
		//sum_predict += (fr0[cpos]);
		//sum_predict_H += (H0[cpos])/ hweight[1];
		//sum_predict_pre += ( - ptmp0_sum - ptmp1_sum)/ hweight[1];
#else
        fr0[cpos] = 1 / hweight[1] * ( - ptmp0_sum + H0[cpos] - ptmp1_sum );
#endif
      } 
//	  printf("\n");
    }
	//if (hor == 832)
	//{
	//	std::cout.precision(10);
	//	std::cout << "predict " << sum_predict << std::endl;
	//	std::cout << "predict_H " << sum_predict_H << std::endl;
	//	std::cout << "predict_pre " << sum_predict_pre << std::endl;
	//}

    if( SHTF_x_range > 0) 
      printf("SHTF (%d pixels)\n", SHTF_x_range ); 
    
    /* checking */
    if( case1 + case2 + case3 != hor * ver ){
      printf( "error in lifting.c -- default: %d, intra: %d\n", case1, case2 );
      exit( 1 );
    }
  }
  else{ // copy frame component
    for( cpos = 0; cpos < hor * ver; cpos++ ){
      // fr0[cpos] = 1 / HPW4[1] * H0[cpos];
#ifdef COPYCOMPENSATION_WEIGHTING 
      fr0[cpos] = H0[cpos] / (HPW4[1] * copycomp_weight_high[level]);
#else
      fr0[cpos] = H0[cpos];
#endif
    }
  } 
#ifdef COUT_MV_UPDATE
  myfile.close();
#endif  
  free( SHTF_filter );
  free( fweight );
  free( ptmp0 );
  free( ptmp1 );

  free( pused0 );
  free( pused1 );
  free( ltmp0 );
  free( ltmp1 );

}


/*****************************************************************************/
/*                               mc_synthesis()                              */
/*                                                                           */
/*                        mv0   mv1   mv2                                    */
/*                        ---> <---   --->                                   */
/*                         |           |           |                         */
/*                     frp | fr0   fr1 | fr2   fr3 |                         */
/*                      . \|  |   / | \|  .   . . .|                         */
/*                      .  \  |  /  |  \  .  .  .  .                         */
/*                      .  |\ | /   |  |\ . .   .  |.                        */
/*                      .  |  H0    |  | H1     .  | H2                      */
/*                      .  |        |  |        .  |                         */
/*                      .  |        |  |        .  |                         */
/*                      .  |        |  |        .  |                         */
/*                      L0 |        L1 |        L2 |                         */
/*                         |           |           |                         */
/*                                                                           */
/*  INPUT:  H0, L1, H1, frp, mv0, mv1, mv2                                   */
/*  OUTPUT: fr0, fr1                                                         */
/*                                                                           */
/*****************************************************************************/ 
// here hor, ver are the dimension for specific component and specific resolution 
// fr0[com_pos] = 1 / HPW1[1] * ( - self_value + H0[com_pos] / comp_weight - neighbor_value );
void
mc_synthesis_with_OBMC( float *fr0, float *fr1, float *H0, float *L1, float *H1, float *frp,
              float *mvx0, float *mvy0, float *mvx1, float *mvy1, 
              float *mvx2, float *mvy2, float *mvx0_int, float *mvy0_int, 
              float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
              vector_ptr mv_ref0, vector_ptr mv_ref1, vector_ptr mv_ref2,
              int hor, int ver, int level, videoinfo info,
			  ImageMEinfo *imagemeinfo, Varblkarrayinfo *varblkarray, int UVcom, int type )
{
	//count += 1;
	//printf("%d\n", count);
  int i, case1, case2;
  int cx, cy, px, py, cpos, ppos, com_pos;
  float cfx, cfy, pfx0, pfy0, pfx1, pfy1;
  float *lweight;
  float *ltmp0, *ltmp1; // 两个补偿完的帧
  unsigned char *pused0, *pused1;
  int    leftx, topy, xblk, yblk, disx, disy, col, row, pos;
  float  self_weight, ver_weight, hor_weight;
  float  self_value, neighbor_value; 
  int case_bi, case_left, case_right, case_intra ;
  int s_level;
  float sad;

  float getnum, getnum2, diff;
  int pixel_count;
  FILE *pc;
  float avg_mv;

  // the spatial resolution reduction in decoder 
  s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));

#ifdef COPYCOMPENSATION_WEIGHTING 
  float comp_weight;
#endif
  //float comp_weight = 1;

  // step 0 - calculation of fr1 (A) / scene_change propagation// 进不来
  if (level < info.t_level) 
  {
    if( mv_ref1 != NULL || mv_ref2 != NULL )
	{ 
      for( ppos = 0; ppos < hor * ver; ppos++ ){
        fr1[ppos] = L1[ppos] / LPW4[1]; 
      }
    } 
	else 
	{
      for( ppos = 0; ppos < hor * ver; ppos++ ){
#ifdef COPYCOMPENSATION_WEIGHTING
        fr1[ppos] = L1[ppos] / (LPW4[1] * copycomp_weight_low[level]);
#else
        fr1[ppos] = L1[ppos];
#endif
      }
    }
    return;
  }

  // variables for update step in MCTF 
  pused0 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused0" );
  pused1 = ( unsigned char * )getarray( hor * ver, sizeof( unsigned char ), "pused1" );
  ltmp0 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp0" );
  ltmp1 = ( float * )getarray( hor * ver, sizeof( float ), "ltmp1" );
  for( i = 0; i < hor * ver; i++ ){ // 初始化
    pused0[i] = UNUSED;
    pused1[i] = UNUSED;
    ltmp0[i] = 0.0;
    ltmp1[i] = 0.0;
  }
 
  // reconstruction of fr1 (A) 根据L1，使用update，重建fr1 
  // step 1 - motion compensation of H0 and H1 为H0 H1做运动补偿 但是只做了双向预测的，其他三种模式是什么情况，为啥没有做？？？
  // MC of H0  为H0做运动补偿
  for( cy = 0; cy < ver; cy++ ){
    for( cx = 0; cx < hor; cx++ ){
      cpos = cy * hor + cx;
      
      if( mvx1[cpos] != (float)HUGE_VAL && mvy1[cpos] != (float)HUGE_VAL ) { // 有参考mv
        pfx0 = cx - mvx1[cpos]; // 参考像素的x位置
        pfy0 = cy - mvy1[cpos]; // 参考像素的y位置
        position( &px, &py, pfx0, pfy0, mvx1[cpos], mvy1[cpos], hor, ver );// 返回整数位置，放在px， py
        ppos = py * hor + px; // 整数位置
        cfx = px + mvx1[cpos];
        cfy = py + mvy1[cpos];
        
        if( pused0[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ){ // 如果之前没有被使用
          pused0[ppos] = USED;
          ltmp0[ppos] = interpolate( cfx, cfy, H0, hor, ver, TYPE ); // 根据运动补偿插出帧
        }
      }
    }
  }
  
  // MC of H1  为H1做运动补偿
  for( cy = 0; cy < ver; cy++ ){
    for( cx = 0; cx < hor; cx++ ){
      cpos = cy * hor + cx;
      
      if( mvx2[cpos] != (float)HUGE_VAL && mvy2[cpos] != (float)HUGE_VAL ) {
        pfx1 = cx - mvx2[cpos];
        pfy1 = cy - mvy2[cpos];
        
        position( &px, &py, pfx1, pfy1, mvx2[cpos], mvy2[cpos], hor, ver );
        ppos = py * hor + px;
        cfx = px + mvx2[cpos];
        cfy = py + mvy2[cpos];
        
        if( pused1[ppos] == UNUSED && inbound( cfx, cfy, hor, ver ) ) {
          pused1[ppos] = USED;
          ltmp1[ppos] = interpolate( cfx, cfy, H1, hor, ver, TYPE );
        }
      }
    }
  }

  case1 = case2 =  0;
  // step 2 - calculation of fr1 (A) 重建fr1，一共有四种模式，双向预测、前向预测、后向预测、帧内预测
  //mv_ref1 = NULL;
  //mv_ref2 = NULL;
  //if (0) {
  //count++;

	if (mv_ref1 != NULL || mv_ref2 != NULL) {
    for( py = 0; py < ver; py++ ) {
      for( px = 0; px < hor; px++ ) {
        ppos = py * hor + px;

        if( pused0[ppos] == USED && pused1[ppos] == USED ){ // 表示是双向的预测
          lweight = LPW1; // default mode (bi-directional)
          case1++;
        }
        else if( pused0[ppos] == UNUSED && pused1[ppos] == USED ){
          lweight = LPW2; // forward mode
          case1++;
        }
        else if( pused0[ppos] == USED && pused1[ppos] == UNUSED ){
		  // for this case info.bi_mv[level] must be 1
		  // if info.bi_mv[level]==0, this case is impossible,
		  // because RIGHT_CONNECTED becomes LEFT_PREDICTED
          lweight = LPW3; // backward mode 
          case1++;
        }
        else{
          lweight = LPW4; // intra mode (unconnected pixels) intra 模式
          case2++;
        }
        assert( ppos >= 0 && ppos < hor * ver );
#ifdef NO_UPDATE
		//lweight = LPW4;
		//lweight[1] = LPW4[1] * copycomp_weight_low[level];
		//lweight[0] = 0;
		//lweight[2] = 0;
		//fr1[ppos] = 1 / ( LPW4[1] * copycomp_weight_low[level]) * (L1[ppos]); // 重建出的帧
		lweight = LPW4;
		fr1[ppos] = 1 / lweight[1] * (L1[ppos] - (lweight[0] * ltmp0[ppos] + lweight[2] * ltmp1[ppos])); // 重建出的帧
#else
        fr1[ppos] = 1 / lweight[1] * ( L1[ppos] -  ( lweight[0] * ltmp0[ppos] + lweight[2] * ltmp1[ppos] )); // 重建出的帧
#endif
      }
    }
  
    assert( case1 + case2 == hor * ver );
  }
  else{ // copy frame component
    for( ppos = 0; ppos < hor * ver; ppos++ ){
#ifdef COPYCOMPENSATION_WEIGHTING
      fr1[ppos] = L1[ppos] / (LPW4[1] * copycomp_weight_low[level]);
#else
      fr1[ppos] = L1[ppos];
#endif
    }
  }

  sad = 0;

  // the OBMC weighting coefficients for U V are sub-sampled version of those for Y
  // the OBMC weighting coefficients for low resolution frame are 
  // sub-sampled version of those in full resolution
  // reconstruction of fr0 (B) 重建fr0
  case_bi = case_left = case_right = case_intra = 0;
  if( mv_ref0 != NULL || mv_ref1 != NULL ){ // 
    for( cy = 0; cy < ver; cy++ ){
      for( cx = 0; cx < hor; cx++ ){
        cpos = ( (cy<<UVcom)<<s_level ) * ( (hor<<UVcom)<<s_level ) + ( (cx<<UVcom)<<s_level );  
		// in Y coordinate with full resolution
		com_pos = cy*hor+cx;   // the position in this component, either Y or U V and this resolution
        
#ifdef COPYCOMPENSATION_WEIGHTING 
        comp_weight = 1.0f;
#endif

#ifdef DIRECTIONAL_IBLOCK_EMPLOYED  // by Yongjun Wu
		// the pixel is in directional IBLOCK  in Y coordinate 
		if (imagemeinfo[cpos].bi_mode == DIRECTIONAL_IBLOCK ) // 帧内预测，跳出本个块
		{
			case_intra++;
			continue; 
		}
#endif

		// information for current block in Y coordinate with full resolution 当前块的信息
		leftx     = imagemeinfo[cpos].leftx;   
		topy      = imagemeinfo[cpos].topy;
		xblk      = imagemeinfo[cpos].blksize; 
		yblk      = imagemeinfo[cpos].blksize;
		// the distance from this pixel to the corner in Y coordinate with full resolution 
		disx = ( (cx<<UVcom)<<s_level ) - leftx;   
		disy = ( (cy<<UVcom)<<s_level ) - topy;
		// weight coefficients for self motion vector and neighbor motion vectors 
		self_weight = imagemeinfo[cpos].self_weight; 
		hor_weight  = imagemeinfo[cpos].h_weight; 
		ver_weight  = imagemeinfo[cpos].v_weight;
		assert(self_weight+hor_weight+ver_weight==(float)1.0);
        self_value = neighbor_value = 0.0;
		// current position ( (cx<<UVcom)<<s_level, (cy<<UVcom)<<s_level ) 
		// in Y coordinate with full resolution
        col = ( cx<<UVcom )<< s_level;   row = ( cy<<UVcom )<< s_level; 
		// contributions from vertical neighbor motion vector
		if (disy < (yblk >> 1))  {    //  top neighbor is effective 离上面近
			pos = (topy-3)*( (hor<<UVcom)<<s_level ) + col; // neighbor position 
			if ( (topy > 0) && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK ) // 不是最上面一行且上面的块不是帧内块
				get_value_for_neighbor_mv(cpos, pos, imagemeinfo, frp, fr1, &neighbor_value, 
					                      &self_weight, ver_weight, cx, cy, info, hor, ver, 
										  UVcom, level);
		} 
		else if (disy >= (yblk >> 1))  { // bottom neighbor is effective 离下面近
			pos = (topy+yblk+2)*( (hor<<UVcom)<<s_level) + col;  // bottom neighbor position
			if  ( (topy+yblk < ( (ver<<UVcom)<<s_level ) )  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK ) 
				get_value_for_neighbor_mv(cpos, pos, imagemeinfo, frp, fr1, &neighbor_value, 
					                      &self_weight, ver_weight, cx, cy, info, hor, ver,
										  UVcom, level);
		}

        // contribution from horizontal neighbor motion vector
		if (disx < (xblk >> 1))  { // left neighbor is effective 
			pos = row*( (hor<<UVcom)<<s_level ) + leftx-3; 
			if ( (leftx > 0)  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK ) 
				 get_value_for_neighbor_mv(cpos, pos, imagemeinfo, frp, fr1, &neighbor_value, 
					                       &self_weight, hor_weight, cx, cy, info, hor, ver,
										   UVcom, level);
		}
		else if (disx >= (xblk >> 1))  { // right neighbor is effective 
			pos = row*( (hor<<UVcom)<<s_level) + leftx+xblk+2;  // right neighbor in full resolution 
			if ( (leftx+xblk < ( (hor<<UVcom)<<s_level ) )  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK )
				 get_value_for_neighbor_mv(cpos, pos, imagemeinfo, frp, fr1, &neighbor_value, 
					                       &self_weight, hor_weight, cx, cy, info, hor, ver,
										   UVcom, level);
		}

		get_value_for_self_mv( cpos, imagemeinfo, frp, fr1, &self_value, self_weight, 
							   cx, cy, info, hor, ver, &case_bi, &case_left, &case_right,
							   UVcom, level);
		
		fr0[com_pos] = 1 / HPW1[1] * ( - self_value + H0[com_pos] / comp_weight - neighbor_value );
#ifdef CNN_TEMPORAL_WAVELET
		/*******写出*****/
		//if (!UVcom)
		//{
		//	char name[256], data_file_name[256], start[10];
		//	strncpy(name, info.bitname, strlen(info.bitname) - 9);
		//	name[strlen(info.bitname) - 9] = '\0';
		//	sprintf(data_file_name, "%s_enc_mv_org.txt", name);
		//	std::ofstream myfile(data_file_name, std::ios::app);
		//	if (myfile.is_open() == NULL)
		//	{
		//		std::cout << "open file_y failed" << std::endl;
		//		exit(0);
		//	}
		//	myfile << com_pos << " " << fr0[com_pos] << " " << frp[com_pos] << " " << fr1[com_pos] << "\n";
		//	myfile.close();
		//}
		/*******写出*****/
#endif

		///////////////////////////////////////////////////
		//MARK BLK MRG
/*		if( type == 1 && level == 0 ){
			if(imagemeinfo[cpos].bi_mode >= 9 || (imagemeinfo[cpos].bi_mode == 7 && imagemeinfo[cpos].aff_mrg == YES) )
				fr0[com_pos] = fr0[com_pos] / 0.4;
		}
*/
/*		if(UVcom == 0 && level == 0){
			if(imagemeinfo[cpos].mrg_blk_bdr == 1)
				fr0[cpos] = 255;
		}
*/		////////////////////////////////////////////////////
         
      }
    }
/*
	if(UVcom == 0 && level == 0){
		diff = 0;	pixel_count = 0;	avg_mv = 0;
//		printf("MV matrix:\n\n");
		for( cy = 0; cy < ver; cy++ ){
		  for( cx = 0; cx < hor; cx++ ){
			cpos = ( (cy<<UVcom)<<s_level ) * ( (hor<<UVcom)<<s_level ) + ( (cx<<UVcom)<<s_level );  
			getnum = (float)HUGE_VAL;	getnum2 = (float)HUGE_VAL;

			if(imagemeinfo[cpos].left_mvx != (float)HUGE_VAL && imagemeinfo[cpos].right_mvx != (float)HUGE_VAL ){
				getnum = (imagemeinfo[cpos].left_mvx - imagemeinfo[cpos].right_mvx)/2;
			}else if(imagemeinfo[cpos].left_mvx != (float)HUGE_VAL){
				getnum = imagemeinfo[cpos].left_mvx;
			}else if(imagemeinfo[cpos].right_mvx != (float)HUGE_VAL){
				getnum = (-1)*imagemeinfo[cpos].right_mvx;
			}
			
			if( getnum != (float)HUGE_VAL ){
				if(buff_frameMEinfo[cpos].left_mvx != (float)HUGE_VAL && buff_frameMEinfo[cpos].right_mvx != (float)HUGE_VAL ){
					getnum = getnum - (buff_frameMEinfo[cpos].left_mvx - buff_frameMEinfo[cpos].right_mvx)/2;
				}else if(buff_frameMEinfo[cpos].left_mvx != (float)HUGE_VAL){
					getnum = getnum - buff_frameMEinfo[cpos].left_mvx;
				}else if(buff_frameMEinfo[cpos].right_mvx != (float)HUGE_VAL){
					getnum = getnum - (-1)*buff_frameMEinfo[cpos].right_mvx;
				}
			}

			if(imagemeinfo[cpos].left_mvy != (float)HUGE_VAL && imagemeinfo[cpos].right_mvy != (float)HUGE_VAL ){
				getnum2 = (imagemeinfo[cpos].left_mvy - imagemeinfo[cpos].right_mvy)/2;
			}else if(imagemeinfo[cpos].left_mvy != (float)HUGE_VAL){
				getnum2 = imagemeinfo[cpos].left_mvy;
			}else if(imagemeinfo[cpos].right_mvy != (float)HUGE_VAL){
				getnum2 = (-1)*imagemeinfo[cpos].right_mvy;
			}
			
			if( getnum2 != (float)HUGE_VAL ){
				if(buff_frameMEinfo[cpos].left_mvy != (float)HUGE_VAL && buff_frameMEinfo[cpos].right_mvy != (float)HUGE_VAL ){
					getnum2 = getnum2 - (buff_frameMEinfo[cpos].left_mvy - buff_frameMEinfo[cpos].right_mvy)/2;
				}else if(buff_frameMEinfo[cpos].left_mvy != (float)HUGE_VAL){
					getnum2 = getnum2 - buff_frameMEinfo[cpos].left_mvy;
				}else if(buff_frameMEinfo[cpos].right_mvy != (float)HUGE_VAL){
					getnum2 = getnum2 - (-1)*buff_frameMEinfo[cpos].right_mvy;
				}
			}

			if( getnum != (float)HUGE_VAL ){
				diff = diff + fabs( sqrt(getnum*getnum + getnum2*getnum2) );
				pixel_count ++;
			}

			if(imagemeinfo[cpos].left_mvx != (float)HUGE_VAL && imagemeinfo[cpos].right_mvx != (float)HUGE_VAL ){
				avg_mv += ( sqrt(imagemeinfo[cpos].right_mvx*imagemeinfo[cpos].right_mvx + imagemeinfo[cpos].right_mvy*imagemeinfo[cpos].right_mvy) +
					sqrt(imagemeinfo[cpos].left_mvx*imagemeinfo[cpos].left_mvx + imagemeinfo[cpos].left_mvy*imagemeinfo[cpos].left_mvy) )/ 2;
			}else if(imagemeinfo[cpos].left_mvx != (float)HUGE_VAL){
				avg_mv += sqrt(imagemeinfo[cpos].left_mvx*imagemeinfo[cpos].left_mvx + imagemeinfo[cpos].left_mvy*imagemeinfo[cpos].left_mvy);
			}else if(imagemeinfo[cpos].right_mvx != (float)HUGE_VAL){
				avg_mv += sqrt(imagemeinfo[cpos].right_mvx*imagemeinfo[cpos].right_mvx + imagemeinfo[cpos].right_mvy*imagemeinfo[cpos].right_mvy);
			}

//copy
			buff_frameMEinfo[cpos].left_mvx = imagemeinfo[cpos].left_mvx;
			buff_frameMEinfo[cpos].left_mvy = imagemeinfo[cpos].left_mvy;
			buff_frameMEinfo[cpos].right_mvx = imagemeinfo[cpos].right_mvx;
			buff_frameMEinfo[cpos].right_mvy = imagemeinfo[cpos].right_mvy;

		  }
		}

		diff = diff / pixel_count;
		avg_mv = avg_mv / pixel_count;
//		printf("AVG diff = %f, pixel_count = %d\n",diff,pixel_count);
//		printf("\n\nMV matrix END:\n");
		pc = fopen("enc_block.txt","at");
		fprintf(pc,"%f\t%f\n",diff,avg_mv);

		assert( diff >= 0 );
		if(diff > global_motion_active)
			global_motion_active = diff;

		fclose(pc);
	}//if UVcom = 0
*/
    assert( case_bi + case_left + case_right + case_intra == hor * ver );
  }//if mv
  else{ // copy frame component  直接就是整帧复制
    for( cpos = 0; cpos < hor * ver; cpos++ ){
#ifdef COPYCOMPENSATION_WEIGHTING 
      fr0[cpos] = H0[cpos] / (HPW4[1] * copycomp_weight_high[level]);
#else
      fr0[cpos] = H0[cpos];
#endif
    }
  } 


  free( pused0 );
  free( pused1 );
  free( ltmp0 );
  free( ltmp1 );

}

