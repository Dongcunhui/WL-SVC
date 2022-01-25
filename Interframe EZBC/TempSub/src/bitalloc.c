#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "basic.h"
#include "structN.h"
#include "miscN.h"
#include "subband.h"
#include "general.h"

#include "layer_mv.h"

#include "zlib.h"
#include "dpx.h"

#define EXTERN  extern
#include "coderN.h"
#include "miscN.h"
#include "ioN.h"
#include "pstatN.h"
#include "basic.h"

enum STRING{
	MAIN, RES
};

/*
 * This function computes the bytes for each sub-bitstream in constant bit plane coding.
 * Since the bytes for subband header (mean, msb) are counted in the bytes of the highest bitplane, all subbands will be 
 * allocated bytes for at least that bitplane.
 * Here subband means sub-bitstream, since I group subbands from Y U and V at the same resolution level together.
 * num_subband: number of subbands of each frame, counting on subbands of Y U V 
 * subband_header_bytes[frame_index][subband_index]: bytes for the subband header
 * delta_r[frame_index][subband bitplane index]: increase of bytes from that subband bitplane
 * num_subband_bp[frame_index]: number of subband bitplanes of that frame
 * new_r[frame_index][subband bitplane index]: bytes upto that subband bitplane.
 *    : since subband header is included in the statistic of the first bitplane, we will sent the first bitplane of
 *      all subbands.
 */

// by Yongjun Wu 不是CBR也不是VBR的进入此函数
// change the bit allocation scheme back to the original version from Rensselaer Polytechnic Institute 
// Namely the bit section is read in and arranged in 3-D style instead of 2-D style from RWTH group 
int bit_alloc_VBR( videoinfo info, unsigned char *qp, long int sum_mv, long int overhead )
{
  enum FRAME_RES { SD2, SD, CIF, HD, INVALID_RES } frame_res;

  int i, j, k, m, counter, GOPsz, curr, last, ncomps, flag, exact_b = 0, b1,
    bit_plane, skipped_subband, log_my_scale = 0;
  int remaining_frs, frame_rate, max_msb = 0,  nbands[3]; // *GOP_header_bytes
  int num_of_GOP, seq_len, num_subband, num_subband_bp;
  int **subband_header_bytes, *total_subband_header_bytes; // all kinds of headers
  int eff_GOPsz[20], GOP_multiplier = 1, dist, s_level;
  long  ***delta_r, ***bytes_upto_bp, ***delta_new_r, ***new_r, **scan, *rate, **subband_bytes;
  double b, *sum_r, budget, sum_k, s_r, slope;
  char seq_name[250], bitplanename[250];
  char *subband_flag, **subband_mask;
  char bit_alloc_name[256];
  FILE *FID;  
  int qtree_depth = 0;
  int  *frame_max_bitplane, roll_level1, roll_level2; // by Yongjun Wu
  int  max_subplane=0, max_bitplane=0, *frame_max_subplane;
  int  exact_bp, exact_subband; 
  
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 ); 
  seq_name[strlen( info.bitname ) - 4] = '\0'; 

  if (info.org_yheight == 480 && info.org_ywidth == 832 )
	  frame_res = SD;
  else if (info.org_yheight == 288 && info.org_ywidth == 352 )
	  frame_res = CIF;
  else if (info.org_yheight == 1080 && info.org_ywidth == 1920)
	  frame_res = HD;
  else if (info.org_yheight == 576 && info.org_ywidth == 704)
	  frame_res = SD2;

  // now each frame is coded independently  --- by Yongjun Wu 
  // it's 2-D EZBC coder instead of 3-D EZBC coder 每一帧都被单独编码，使用的是2D而不是3D EZBC
  seq_len    = info.last - info.start + 1;
  frame_rate = info.framerate;
  num_of_GOP = get_GOP_num( info );

  if( info.denoise_flag == YES ) {
    GOP_multiplier = ( YUV420 == 1 ) ? 2 : 4;
    seq_len *= GOP_multiplier;
  }

  // maximum bitplane number in each frame by Yongjun Wu
  frame_max_bitplane = (int *)getarray(seq_len, sizeof(int), "frame_max_bitplane");  
  // maximum sub-bitplane number in each frame 
  frame_max_subplane = (int *)getarray(seq_len, sizeof(int), "frame_max_subplane");
  // delta_r[i][j][k]: each data item, [i]: frame number, [j]: subband number, [k]: sub-bitplane number 
  delta_r = ( long *** )getarray( seq_len, sizeof( long * ), "delta_r" );
  // accumulated data items 
  bytes_upto_bp = ( long *** )getarray( seq_len, sizeof( long * ), "bytes_upto_bp" );

  subband_header_bytes = ( int ** )getarray( seq_len, sizeof( int * ), "subband_header_bytes" );
  total_subband_header_bytes = ( int * )getarray( seq_len, sizeof( int ), "total_subband_header_bytes" );
  scan                 = (long **)getarray(seq_len, sizeof(long *), "scan");

  sprintf( bitplanename, "%s.rd_sample_dat", seq_name );
  bitplanename[strlen( seq_name ) + 15] = '\0';  //
  FID = fopen( bitplanename, "rb" );
  if( FID == NULL ) {
    printf( "can not open file bit_rd_sample.dat\n" );
    exit( 1 );
  }

  ncomps = 1;
  num_subband = 0;
  for( m = 0; m < ncomps; m++ ) { // 只进行一次
    // read in the number of substreams of each frame    
    fscanf( FID, "%d\n", &( nbands[m] ) );
#ifdef   	FREQUENCY_ROLL_OFF
	if (frame_res == HD ) //Added by Yuan Liu
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res == SD2)
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif

		roll_level1 = 5; roll_level2=4;  // the subband levels where we need frequency roll-off
	}
	else if (frame_res == SD )
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif

		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res==CIF)
	{
		num_subband += nbands[m]+2;    // the independent coded bitstreams 
		roll_level1=roll_level2=4;
	}
	else printf("Error in frequency roll-off!\n"); 

#else
    num_subband += nbands[m];      // the independent coded bitstreams 
#endif
  }
  fseek( FID, 0, SEEK_SET );

  for(i = 0; i < seq_len; i++){
	delta_r[i]       = (long **)getarray(num_subband, sizeof(long *), "delta_r[i]");
    bytes_upto_bp[i] = (long **)getarray(num_subband, sizeof(long *), "bytes_upto_bp[i]");
	subband_header_bytes[i] = (int *)getarray(num_subband, sizeof(int), "subband_header_bytes[i]");
  }

  // generate subband mask 生成子带mask
  subband_mask = ( char ** )getarray( seq_len, sizeof( char * ), "subband_mask(bitalloc.c)" );
  for( i = 0; i < seq_len; i++ ) {
    subband_mask[i] =  ( char * )getarray( num_subband, sizeof( char ), "subband_mask[i] (bitalloc.c)" );
    for( k = 0; k < num_subband; k++ ) {
      subband_mask[i][k] = 0;
    }
  }

  curr = info.start;
  last = info.last;
  counter = 0;
  skipped_subband = 0;
  // 设置mask等
  while( curr <= last ) { // 逐个gop处理
    remaining_frs = last - curr + 1;  
    // determine effective GOP size in level i
    GOPsz = info.GOPsz;
    for( i = 0; i <= info.tPyrLev; i++ ) {
      if ( remaining_frs < info.GOPsz )
        eff_GOPsz[i] = (int) ceil ((double) remaining_frs / (pow (2, i)));  
      else
        eff_GOPsz[i] = GOPsz;  
      GOPsz /= 2;  
    }
  
    // mask useless temporal high subband frames  --- now each frame is coded independently  把没用的时域高频子带盖住
    GOPsz = info.GOPsz;
    for( i = 0; i < info.t_level; i++ ) {
      for( j = eff_GOPsz[i+1]; j < eff_GOPsz[i]; j++ ) {
        for( k = 0; k < num_subband; k++ ) {
          subband_mask[counter + j][k] = 1;
          skipped_subband++;
        }
      }
    }

    //mask useless spatial high subband frames; if current GOP < default GOP size, GOPsz could be 0 把没用的空域高频子带盖住
    for( i = 0; i < MY_MAX( eff_GOPsz[0], 1 ); i++ ) {
      s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
      ncomps = 1; // since Y U V are grouped together
      m = 0;
      for( j = 0; j < ncomps; j++ ) {
        m += nbands[j];

#ifdef   	FREQUENCY_ROLL_OFF
		int subband_index;
		if (frame_res == HD )
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
	    else if (frame_res == SD2)
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
		else if (frame_res == SD )
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
		else if (frame_res==CIF)
			subband_index = nbands[j]+2;    // the independent coded bitstreams 
		else printf("Error in frequency roll-off!\n"); 
#endif

        for( k = m - s_level; k < m; k++ ) {

#ifdef   	FREQUENCY_ROLL_OFF	
			subband_mask[counter+i][subband_index-1] = 1;
			subband_index--; 
			skipped_subband++;

#ifdef ROLL_STRUCTURE_ONE
			if (k==roll_level1 || k==roll_level2)  // additional two bands 
			{
				for (int s=0; s<2; s++)
				{
					subband_mask[counter+i][subband_index-1] = 1;
					subband_index--; 
					skipped_subband++;
				}
			}
#endif

#ifdef ROLL_STRUCTURE_TWO
			if (k==roll_level1 || k==roll_level2)  
			{
				int band_num; 
				if (k==roll_level1) band_num = 11; // additional 11 bands 
				if (k==roll_level2) band_num = 2;  // additional 2  bands 
				for (int s=0; s<band_num; s++)
				{
					subband_mask[counter+i][subband_index-1] = 1;
					subband_index--; 
					skipped_subband++;
				}
			}
#endif

#else       // no frequency roll off 
			subband_mask[counter+i][k] = 1;
			skipped_subband++;
#endif

        }
      }
    }

    // mask useless denoise subbands
    if( info.denoise_flag == YES ){
      if( info.s_level > 0 ) {
        for( i = 1; i < GOP_multiplier; i++ ) {
          // mark all spatial highpass frames
          for( j = 0; j < eff_GOPsz[0]; j++ ) {
            for( k = 0; k < num_subband; k++ ) {
              subband_mask[counter + i*eff_GOPsz[0] + j][k] = 1;
              skipped_subband++;
            }
          }
        }
      } 
      else if( info.t_level > 0 ) {
        dist = 1 << info.t_level;
        for( i = 1; i < GOP_multiplier; i++ ) {
          // first mark all
          for( j = 0; j < eff_GOPsz[0]; j++ ) {
            for( k = 0; k < num_subband; k++ ) {
              subband_mask[counter + i*eff_GOPsz[0] + j][k] = 1;
              skipped_subband++;
            }
          }  
          // then select temporal subbands that remain in bit stream
          for( j = 0; j < eff_GOPsz[0]; j += dist ) {
            for( k = 0; k < num_subband; k++ ) {
              subband_mask[counter + i*eff_GOPsz[0] + j][k] = 0;
              skipped_subband--;
            }
          }
        }
      }
    }

    curr += info.GOPsz;
    counter += info.GOPsz * GOP_multiplier;

  }                             /* while */
  
  // *** overhead bytes for subband allocation table ***  为子带分配表的字节支出
  // sizeof(unsigned short int) is for header, 3 is the extra bytes at the end of each substream
  overhead += ( sizeof( unsigned short int ) + 3 ) * ( num_subband * seq_len - skipped_subband );
  for( i = 0; i < seq_len; i++ ) {
    num_subband = 0;
    for( m = 0; m < ncomps; m++ ) {
      // read in the number of subbands of each component          
      fscanf( FID, "%d\n", &( nbands[m] ) );
#ifdef   	FREQUENCY_ROLL_OFF
	  if (frame_res == HD)

#ifdef ROLL_STRUCTURE_ONE
	  {
		  num_subband += nbands[m] + 2 + 2;    // the independent coded bitstreams 
	  }

#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res == SD2)

#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;    // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res == SD )
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2;    // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res==CIF)
		num_subband += nbands[m]+2;      // the independent coded bitstreams 
	else printf("Error in frequency roll-off!\n"); 
#else
    num_subband += nbands[m];           // the independent coded bitstreams 
#endif

    }

	// read in the maximum bit plane of a frame 读一帧的最大比特平面 by Yongjun Wu
	fscanf(FID, "%d\n", &(frame_max_bitplane[i]));    
	if(frame_max_bitplane[i] > max_bitplane) max_bitplane = frame_max_bitplane[i];
    // read in the maximum sub-bitplane of a frame
    fscanf( FID, "%d\n", &( frame_max_subplane[i] ) );
    if( frame_max_subplane[i] > max_subplane ) max_subplane = frame_max_subplane[i];

	// scan bit planes of each subband
	for(j=0; j<num_subband; j++){
      delta_r[i][j] = (long *)getarray(frame_max_subplane[i]+1, sizeof(long), "delta_r[i][j]");
      bytes_upto_bp[i][j] = (long *)getarray(frame_max_subplane[i]+1, sizeof(long), "bytes_upto_bp[i][j]");	
	}

    for(j=0; j<num_subband; j++){
		subband_header_bytes[i][j] = 0;
	}
	total_subband_header_bytes[i] = 0;

	for(j=frame_max_subplane[i]; j>=0; j--){
	  for(m=0; m<num_subband; m++){
		//read in bitplane index
        fscanf(FID, "%d", &bit_plane);   
        if (j != bit_plane){
          printf(" j %d bit_plane %d bit plane index error !\n", j, bit_plane);
		  exit(1);
		}		 
	    //read in the size of each subband bitplane
		fscanf(FID, "%d %d\n", &(delta_r[i][m][j]), &(bytes_upto_bp[i][m][j]));    
			
		// subband head bytes are located in the first place
		if(subband_header_bytes[i][m]==0 && delta_r[i][m][j] !=0){
			subband_header_bytes[i][m] = delta_r[i][m][j];
		    total_subband_header_bytes[i] +=  subband_header_bytes[i][m];
		}
	  } //m
	}//j
  }   // loop of i 
  fclose( FID );


  
  // the set of variables for the idea of frequency roll-off 
  int  *shift_per_bitplane, max_per_shift =0;  
  int   max_bitplane_shift, max_subplane_shift;
  float  *subband_weighting  = (float *)getarray(num_subband, sizeof(float), "subband_weighting");
  // initialize all the varialbles to be 0
  memset(subband_weighting, 0, num_subband*sizeof(float));
  shift_per_bitplane = (int *)getarray(seq_len, sizeof(int), "shift_per_bitplane");
  memset(shift_per_bitplane, 0, seq_len*sizeof(int));
  max_subplane_shift = max_per_shift = max_bitplane_shift = 0; 

#ifdef FREQUENCY_ROLL_OFF
  FILE *froll;
  float roll_factors[15]; 
  int   factor_num; 

  for (i=0; i<seq_len; i++)
  {
	  shift_per_bitplane[i] = (frame_max_subplane[i]+1)/(frame_max_bitplane[i]+1); 
	  if (max_per_shift < shift_per_bitplane[i] )
			max_per_shift = shift_per_bitplane[i];
  }
  max_bitplane_shift = 5; 
  max_subplane_shift = max_bitplane_shift * max_per_shift; 

  // for SD sequence, both CIF and QCIF resolution need the idea of frequency roll-off 
  // the set of weighting coefficients are predefined 

  if ( (info.s_level==1 || info.s_level==2) && frame_res == HD )  
  {

#ifdef ROLL_STRUCTURE_ONE
	  // read the predefined scaling factors from a text file
	  froll=fopen("roll_factors_structure_one.txt", "rt"); 
	  for (factor_num=0; factor_num<6; factor_num++)
		  fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
	  if (info.s_level==1)
	  {
		  subband_weighting[7] = roll_factors[0]; //0.5;  
		  subband_weighting[8] = roll_factors[1]; //1.5; 
		  subband_weighting[9] = roll_factors[2]; //1.5; 
	  }

	  if (info.s_level==2)
	  {
		  subband_weighting[4] = roll_factors[3]; //0.5;
		  subband_weighting[5] = roll_factors[4]; //1.0; 
		  subband_weighting[6] = roll_factors[5]; //1.0; 
	  }
#endif


#ifdef ROLL_STRUCTURE_TWO
	  // read the predefined scaling factors from a text file
	  froll=fopen("roll_factors_structure_two.txt", "rt"); 
	  for (factor_num=0; factor_num<15; factor_num++)
		  fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
	  if (info.s_level==1)
		  for (factor_num=0; factor_num<12; factor_num++)
			  subband_weighting[7+factor_num] = roll_factors[factor_num];   
	  if (info.s_level==2)
	  {
		  subband_weighting[4] = roll_factors[12]; 
		  subband_weighting[5] = roll_factors[13]; 
		  subband_weighting[6] = roll_factors[14];  
	  }
#endif

  }

  if ( (info.s_level==1 || info.s_level==2) && (frame_res == SD || frame_res == SD2) )  
  {

#ifdef ROLL_STRUCTURE_ONE
	  // read the predefined scaling factors from a text file
	  froll=fopen("roll_factors_structure_one.txt", "rt"); 
	  for (factor_num=0; factor_num<6; factor_num++)
		  fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
	  if (info.s_level==1)
	  {
		  subband_weighting[7] = roll_factors[0]; //0.5;  
		  subband_weighting[8] = roll_factors[1]; //1.5; 
		  subband_weighting[9] = roll_factors[2]; //1.5; 
	  }

	  if (info.s_level==2)
	  {
		  subband_weighting[4] = roll_factors[3]; //0.5;
		  subband_weighting[5] = roll_factors[4]; //1.0; 
		  subband_weighting[6] = roll_factors[5]; //1.0; 
	  }
#endif


#ifdef ROLL_STRUCTURE_TWO
	  // read the predefined scaling factors from a text file
	  froll=fopen("roll_factors_structure_two.txt", "rt"); 
	  for (factor_num=0; factor_num<15; factor_num++)
		  fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
	  if (info.s_level==1)
		  for (factor_num=0; factor_num<12; factor_num++)
			  subband_weighting[7+factor_num] = roll_factors[factor_num];   
	  if (info.s_level==2)
	  {
		  subband_weighting[4] = roll_factors[12]; 
		  subband_weighting[5] = roll_factors[13]; 
		  subband_weighting[6] = roll_factors[14];  
	  }
#endif

  }
  
  // for CIF sequence, only CIF resolution needs the idea of frequency roll-off 
  if ( info.s_level==1 && frame_res == CIF ) 
  {
	  froll=fopen("roll_factors_structure_one.txt", "rt"); 
	  for (factor_num=0; factor_num<3; factor_num++)
		fscanf(froll, "%f", &(roll_factors[factor_num]));
	  fclose(froll); 
	  subband_weighting[4] = roll_factors[0]; //0.0;
	  subband_weighting[5] = roll_factors[1]; //0.0; 
	  subband_weighting[6] = roll_factors[2]; //0.5; 
  }


#endif

  // for interleaving use 
  printf("num_subband = %d, max_subplane = %d, max_subplane_shift = %d\n",num_subband, max_subplane, max_subplane_shift);
  num_subband_bp = num_subband * (max_subplane+1+max_subplane_shift); 
  for (i=0; i<seq_len; i++)
    scan[i] = (long *)getarray(num_subband_bp, sizeof(long), "scan[i]");

  // move data in delta_r to delta_new_r to make all substream look like having the same number of bitplanes  
  // in order to simplify the following bit allocation algorithm
  new_r = (long ***)getarray(seq_len, sizeof(long **), "new_r");
  delta_new_r = (long ***)getarray(seq_len, sizeof(long **), "delta_new_r");
  subband_flag = (char *)getarray(num_subband, sizeof(char), "subband_flag");
  // manipulate the bitstream shift in different subbands via subband_weighting & shift_per_bitplane
  for (i=0; i<seq_len; i++){
	new_r[i] = (long **)getarray(num_subband, sizeof(long *), "new_r[i]");
	delta_new_r[i] = (long **)getarray(num_subband, sizeof(long *), "delta_new_r[i]");
	for(j=0; j<num_subband; j++){
      new_r[i][j] = (long *)getarray(max_subplane+1+max_subplane_shift, sizeof(long), "new_r[i][j]");
      delta_new_r[i][j] = (long *)getarray(max_subplane+1+max_subplane_shift, sizeof(long), "delta_new_r[i][j]");
	  for(k=0; k<=max_subplane+max_subplane_shift; k++){
	    new_r[i][j][k]=0;
	    delta_new_r[i][j][k]=0;
	  }
	}

    for(m=0; m<num_subband; m++){
	  subband_flag[m]=0;
	}

	for(j=frame_max_subplane[i]; j>=0; j--){
	  for(m=0; m<num_subband; m++){
		if(subband_flag[m]==0 && delta_r[i][m][j] !=0){
		  delta_new_r[i][m][j+max_subplane_shift-(int)(subband_weighting[m]*shift_per_bitplane[i])]=0;
		  subband_flag[m]=1;
		}
		else
		  delta_new_r[i][m][j+max_subplane_shift-(int)(subband_weighting[m]*shift_per_bitplane[i])]=delta_r[i][m][j];
	  }
	}
	for(m=0; m<num_subband; m++){ 
	  delta_new_r[i][m][max_subplane+max_subplane_shift] = subband_header_bytes[i][m];
	}
  }
  free(subband_flag);


  // mask useless subbands according to the subband_mask
  for (i=0; i<seq_len; i++){
	for(m=0; m<num_subband; m++){
	  if(subband_mask[i][m] == 1){
		assert(0);
	    for(j=max_subplane+max_subplane_shift; j>=0; j--){
		  delta_new_r[i][m][j]=0;
		}
	  }
	}
  }

  // form the accumulated sum of bits for each subband and save the sum in new_r 
  for(i=0; i<seq_len; i++){	
	for(j=0; j<num_subband; j++){ 
  	  new_r[i][j][max_subplane+max_subplane_shift]=delta_new_r[i][j][max_subplane+max_subplane_shift];
	  for(k=max_subplane-1+max_subplane_shift; k>=0; k--){		
	    new_r[i][j][k] = new_r[i][j][k+1] + delta_new_r[i][j][k]; 
	  }
	}
  }


  for(i=0; i<seq_len; i++){	
	int count = num_subband_bp-1;
	for(k=max_subplane+max_subplane_shift; k>=0; k--){		
	  for(j=0; j<num_subband; j++){ 
	    if(count == num_subband_bp-1)
		  scan[i][count] = delta_new_r[i][j][max_subplane+max_subplane_shift];
	    else
	      scan[i][count] = scan[i][count+1] +  delta_new_r[i][j][k];
	    count--;
	  }
	}
  }

  // rate increasing curve is in sum_r
  sum_r = (double *)getarray(num_subband_bp, sizeof(double), "sum_r");
  for (j=0; j<num_subband_bp; j++){
    sum_r[j] = 0;
  }
  for (j=num_subband_bp-1; j>=0; j--){
    for (i=0; i<seq_len; i++){
      sum_r[j] = sum_r[j] + scan[i][j];
	}
  } // loop of i (no_stream)


  rate = ( long * )getarray( seq_len, sizeof( long ), "rate" );
  // allocate bytes to each frame
  budget = ( double )info.bitrate * 1000 / 8 * ( info.last - info.start + 1 ) / frame_rate; 
  budget = budget - overhead - sum_mv;

  printf("Overhead = %d, sum_mv = %ld, budget = %f\n", overhead, sum_mv,budget);
  
  if( sum_r[num_subband_bp - 1] > budget ) {   // bit plane max_msb
    printf( " bitalloc.c: can not handle this case: budget is too small ( %f ).\n", budget );
    // exit( 1 );
    return 0;
  } else {                      // else1
    if( sum_r[num_subband_bp - 1] == budget ) {
      exact_b = num_subband_bp - 1;
      flag = 0;
    } else {                    // else2
      if( sum_r[0] < budget ) {
        exact_b = 0;
        flag = 0;
        budget = sum_r[0];
        printf( " bitalloc.c: no more data in the bit stream for allocation.\n" );
      } else {                  //else3    // find the interval
        for( j = num_subband_bp - 2; j >= 0; j-- ) {   // the other bit planes    
          if( sum_r[j] > budget ) {
            b1 = j;             // record the index 
            flag = 1;
            break;
          } else {
            if( sum_r[j] == budget ) {
              exact_b = j;
              flag = 0;
              break;
            }
          }
        }                       // for(j)
      }                         // else3
    }                           // else2
  }                             // else1


  subband_bytes = ( long ** )getarray( seq_len, sizeof( long * ), "subband_bytes" );

  for( i = 0; i < seq_len; i++ ) {
    subband_bytes[i] =  ( long * )getarray( num_subband, sizeof( long ), "num_subband" );
    for( j = 0; j < num_subband; j++ ) {
      subband_bytes[i][j] = 0;
    }
  }


	if (flag == 0){  // perfect match
	  exact_bp = exact_b / num_subband; 
	  exact_subband = num_subband - 1 - exact_b%num_subband;
 	  for(i=0; i<seq_len; i++){			
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}
		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}
	  }

	}

	// get the interpolated bytes
	if (flag == 1){
      sum_k = 0.;
      s_r = 0.;
      for (j=0; j<seq_len; j++){
        sum_k = sum_k + scan[j][b1+1] - scan[j][b1];
        s_r = s_r + scan[j][b1];
      }
      b = (budget - s_r + b1 * sum_k)/sum_k;
   
      for(j=0; j<seq_len; j++){
        slope = scan[j][b1+1] - scan[j][b1];
        rate[j] = (long)(scan[j][b1] + slope*(b-b1));		
      }
		
	  exact_bp = b1 / num_subband; 
	  exact_subband = num_subband - 1 - b1%num_subband;
  	  for(i=0; i<seq_len; i++){		    
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}

		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}	 
		subband_bytes[i][exact_subband] -= scan[i][b1]-rate[i];
	  }
	} 


  for( i = 0; i < seq_len; i++ ) {
    for( j = 0; j < num_subband; j++ ) {
      // if( subband_bytes[i][j] > 0 ) {
		if( ! subband_mask[i][j] ) {
        subband_bytes[i][j] += 3;
      }
    }
  }



  /**********/
  /* OUTPUT */
  /**********/
  sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, info.bitrate );
  FID = fopen( bit_alloc_name, "wb" );

  fprintf( FID, "%d\n", num_subband );

  for( i = 0; i < seq_len; i++ ) {
    for( m = 0; m < num_subband; m++ ) {
      fprintf( FID, "%ld\n", subband_bytes[i][m] );
    }
  }
  fclose( FID );


  for(i=0; i<seq_len; i++){
	for(j=0; j<num_subband; j++){
	  free(delta_r[i][j]);
	  free(bytes_upto_bp[i][j]);
	  free(new_r[i][j]);
	  free(delta_new_r[i][j]);
	}

	free(subband_bytes[i]);
    free(subband_header_bytes[i]);

	free(delta_r[i]);
	free(new_r[i]);
	free(delta_new_r[i]);
	free(bytes_upto_bp[i]);

	free(scan[i]);
  }
  free(subband_bytes);
  free(delta_r);
  free(new_r);
  free(delta_new_r);
  free(bytes_upto_bp);
  free(rate);
  free(sum_r);
  free(frame_max_bitplane);
  free(frame_max_subplane);
  free(subband_header_bytes);
  free(total_subband_header_bytes);


  return 1;
}



/***************************************************/
/***************************************************/
/***************************************************/
/*
	CBR test
*/
/***************************************************/
/***************************************************/
/***************************************************/

int bit_alloc_CBR( videoinfo info, unsigned char *qp, long int overhead, unsigned long *gop_sum_mv, int curr_last )
{
  enum FRAME_RES { SD2, SD, CIF, HD, INVALID_RES } frame_res;

  long int GOP_overhead;

  int i, j, k, m, counter, GOPsz, curr, last, ncomps, flag, exact_b = 0, b1,
    bit_plane, skipped_subband, log_my_scale = 0;
  int remaining_frs, frame_rate, max_msb = 0,  nbands[3]; // *GOP_header_bytes
  int num_of_GOP, seq_len, num_subband, num_subband_bp;

  int *GOP_num_subband_bp;

  int **subband_header_bytes, *total_subband_header_bytes; // all kinds of headers
  int eff_GOPsz[20], GOP_multiplier = 1, dist, s_level;
  long  ***delta_r, ***bytes_upto_bp, ***delta_new_r, ***new_r, **scan, *rate, **subband_bytes;
  double b, *sum_r, budget, GOP_budget, sum_k, s_r, slope;
  char seq_name[250], bitplanename[250];
  char *subband_flag, **subband_mask;
  char bit_alloc_name[256];
  FILE *FID;  
  int qtree_depth = 0;
  int  *frame_max_bitplane, roll_level1, roll_level2; // by Yongjun Wu
  int  max_subplane=0, max_bitplane=0, *frame_max_subplane;
  int  exact_bp, exact_subband; 
  
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 ); 
  seq_name[strlen( info.bitname ) - 4] = '\0'; 

  if (info.org_yheight == 480 && info.org_ywidth == 832 )
	  frame_res = SD;
  else if (info.org_yheight == 288 && info.org_ywidth == 352 )
	  frame_res = CIF;
  else if (info.org_yheight == 1080 && info.org_ywidth == 1920)
	  frame_res = HD;
  else if (info.org_yheight == 576 && info.org_ywidth == 704)
	  frame_res = SD2;

  // now each frame is coded independently  --- by Yongjun Wu 
  // it's 2-D EZBC coder instead of 3-D EZBC coder 
  seq_len    = curr_last - info.start + 1;
  frame_rate = info.framerate;
  num_of_GOP = get_GOP_num( info );

  if( info.denoise_flag == YES ) {
    GOP_multiplier = ( YUV420 == 1 ) ? 2 : 4;
    seq_len *= GOP_multiplier;
  }

  printf("seq_len = %d\n",seq_len);

  // maximum bitplane number in each frame by Yongjun Wu
  frame_max_bitplane = (int *)getarray(seq_len, sizeof(int), "frame_max_bitplane");  
  // maximum sub-bitplane number in each frame 
  frame_max_subplane = (int *)getarray(seq_len, sizeof(int), "frame_max_subplane");
  // delta_r[i][j][k]: each data item, [i]: frame number, [j]: subband number, [k]: sub-bitplane number 
  delta_r = ( long *** )getarray( seq_len, sizeof( long * ), "delta_r" );
  // accumulated data items 
  bytes_upto_bp = ( long *** )getarray( seq_len, sizeof( long * ), "bytes_upto_bp" );

  subband_header_bytes = ( int ** )getarray( seq_len, sizeof( int * ), "subband_header_bytes" );
  total_subband_header_bytes = ( int * )getarray( seq_len, sizeof( int ), "total_subband_header_bytes" );
  scan                 = (long **)getarray(seq_len, sizeof(long *), "scan");

  GOP_num_subband_bp = (int *)getarray(num_of_GOP, sizeof( int * ), "GOP_num_subband_bp" );

  sprintf( bitplanename, "%s.rd_sample_dat", seq_name );
  bitplanename[strlen( seq_name ) + 15] = '\0';  //
  FID = fopen( bitplanename, "rb" );
  if( FID == NULL ) {
    printf( "can not open file bit_rd_sample.dat\n" );
    exit( 1 );
  }

  ncomps = 1;
  num_subband = 0;
  for( m = 0; m < ncomps; m++ ) {
    // read in the number of substreams of each frame    
    fscanf( FID, "%d\n", &( nbands[m] ) );
#ifdef   	FREQUENCY_ROLL_OFF
	if (frame_res == HD ) //Added by Yuan Liu
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res == SD2)
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif

		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res == SD )
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif

		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res==CIF)
	{
		num_subband += nbands[m]+2;    // the independent coded bitstreams 
		roll_level1=roll_level2=4;
	}
	else printf("Error in frequency roll-off!\n"); 

#else
    num_subband += nbands[m];      // the independent coded bitstreams 
#endif
  }
  fseek( FID, 0, SEEK_SET );

  for(i = 0; i < seq_len; i++){
	delta_r[i]       = (long **)getarray(num_subband, sizeof(long *), "delta_r[i]");
    bytes_upto_bp[i] = (long **)getarray(num_subband, sizeof(long *), "bytes_upto_bp[i]");
	subband_header_bytes[i] = (int *)getarray(num_subband, sizeof(int), "subband_header_bytes[i]");
  }

  // generate subband mask 
  subband_mask = ( char ** )getarray( seq_len, sizeof( char * ), "subband_mask(bitalloc.c)" );
  for( i = 0; i < seq_len; i++ ) {
    subband_mask[i] =  ( char * )getarray( num_subband, sizeof( char ), "subband_mask[i] (bitalloc.c)" );
    for( k = 0; k < num_subband; k++ ) {
      subband_mask[i][k] = 0;
    }
  }

  curr = info.start;
  last = curr_last;
  counter = 0;
  skipped_subband = 0;
  while( curr <= last ) {
    remaining_frs = last - curr + 1;  
    // determine effective GOP size in level i
    GOPsz = info.GOPsz;
    for( i = 0; i <= info.tPyrLev; i++ ) {
      if ( remaining_frs < info.GOPsz )
        eff_GOPsz[i] = (int) ceil ((double) remaining_frs / (pow (2, i)));  
      else
        eff_GOPsz[i] = GOPsz;  
      GOPsz /= 2;  
    }
  
    //mask useless temporal high subband frames  --- now each frame is coded independently
    GOPsz = info.GOPsz;
    for( i = 0; i < info.t_level; i++ ) {
      for( j = eff_GOPsz[i+1]; j < eff_GOPsz[i]; j++ ) {
        for( k = 0; k < num_subband; k++ ) {
		  assert(0);
          subband_mask[counter + j][k] = 1;
          skipped_subband++;
        }
      }
    }

    //mask useless spatial high subband frames; if current GOP < default GOP size, GOPsz could be 0
    for( i = 0; i < MY_MAX( eff_GOPsz[0], 1 ); i++ ) {
      s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
      ncomps = 1; // since Y U V are grouped together
      m = 0;
      for( j = 0; j < ncomps; j++ ) {
        m += nbands[j];

#ifdef   	FREQUENCY_ROLL_OFF
		int subband_index;
		if (frame_res == HD )
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
	    else if (frame_res == SD2)
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
		else if (frame_res == SD )
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
		else if (frame_res==CIF)
			subband_index = nbands[j]+2;    // the independent coded bitstreams 
		else printf("Error in frequency roll-off!\n"); 
#endif

      }
    }

    curr += info.GOPsz;
    counter += info.GOPsz * GOP_multiplier;
  }                             /* while */
  
  // *** overhead bytes for subband allocation table ***
  // sizeof(unsigned short int) is for header, 3 is the extra bytes at the end of each substream

//  printf("skipped_subband = %d\n",skipped_subband);

  overhead += ( sizeof( unsigned short int ) + 3 ) * ( num_subband * (info.last - info.start + 1) - skipped_subband );

  GOP_overhead = overhead / num_of_GOP;

  for( i = 0; i < seq_len; i++ ) {
    num_subband = 0;
    for( m = 0; m < ncomps; m++ ) {
      // read in the number of subbands of each component          
      fscanf( FID, "%d\n", &( nbands[m] ) );
#ifdef   	FREQUENCY_ROLL_OFF
	if (frame_res == HD )

#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;    // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res == SD2)

#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;    // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res == SD )
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2;    // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res==CIF)
		num_subband += nbands[m]+2;      // the independent coded bitstreams 
	else printf("Error in frequency roll-off!\n"); 
#else
    num_subband += nbands[m];           // the independent coded bitstreams 
#endif

    }

	// read in the maximum bit plane of a frame by Yongjun Wu
	fscanf(FID, "%d\n", &(frame_max_bitplane[i]));    
	if(frame_max_bitplane[i] > max_bitplane) max_bitplane = frame_max_bitplane[i];
    // read in the maximum sub-bitplane of a frame
    fscanf( FID, "%d\n", &( frame_max_subplane[i] ) );
    if( frame_max_subplane[i] > max_subplane ) max_subplane = frame_max_subplane[i];

	// scan bit planes of each subband
	for(j=0; j<num_subband; j++){
      delta_r[i][j] = (long *)getarray(frame_max_subplane[i]+1, sizeof(long), "delta_r[i][j]");
      bytes_upto_bp[i][j] = (long *)getarray(frame_max_subplane[i]+1, sizeof(long), "bytes_upto_bp[i][j]");	
	}

    for(j=0; j<num_subband; j++){
		subband_header_bytes[i][j] = 0;
	}
	total_subband_header_bytes[i] = 0;

	for(j=frame_max_subplane[i]; j>=0; j--){

	  for(m=0; m<num_subband; m++){
		//read in bitplane index
        fscanf(FID, "%d", &bit_plane);      
        if (j != bit_plane){
          printf(" j %d bit_plane %d bit plane index error !\n", j, bit_plane);
		  exit(1);
		}		 
	    //read in the size of each subband bitplane
		fscanf(FID, "%d %d\n", &(delta_r[i][m][j]), &(bytes_upto_bp[i][m][j]));    
			
		// subband head bytes are located in the first place
		if(subband_header_bytes[i][m]==0 && delta_r[i][m][j] !=0){
			subband_header_bytes[i][m] = delta_r[i][m][j];
		    total_subband_header_bytes[i] +=  subband_header_bytes[i][m];
		}
	  } //m
	}//j
  }   // loop of i 
  fclose( FID );


  
  // the set of variables for the idea of frequency roll-off 
  int  *shift_per_bitplane, max_per_shift =0;  
  int   max_bitplane_shift, max_subplane_shift;
  float  *subband_weighting  = (float *)getarray(num_subband, sizeof(float), "subband_weighting");
  // initialize all the varialbles to be 0
  memset(subband_weighting, 0, num_subband*sizeof(float));
  shift_per_bitplane = (int *)getarray(seq_len, sizeof(int), "shift_per_bitplane");
  memset(shift_per_bitplane, 0, seq_len*sizeof(int));
  max_subplane_shift = max_per_shift = max_bitplane_shift = 0; 

#ifdef FREQUENCY_ROLL_OFF
  FILE *froll;
  float roll_factors[15]; 
  int   factor_num; 

  for (i=0; i<seq_len; i++)
  {
	  shift_per_bitplane[i] = (frame_max_subplane[i]+1)/(frame_max_bitplane[i]+1); 

	  if (max_per_shift < shift_per_bitplane[i] )
		max_per_shift = shift_per_bitplane[i];
  }
  max_bitplane_shift = 5; 
  max_subplane_shift = max_bitplane_shift * max_per_shift; 

#endif

  // for interleaving use 
  printf("num_subband = %d, max_subplane = %d, max_subplane_shift = %d\n",num_subband, max_subplane, max_subplane_shift);
  num_subband_bp = num_subband * (max_subplane+1+max_subplane_shift); 
  for (i=0; i<seq_len; i++)
    scan[i] = (long *)getarray(num_subband_bp, sizeof(long), "scan[i]");

  // move data in delta_r to delta_new_r to make all substream look like having the same number of bitplanes  
  // in order to simplify the following bit allocation algorithm
  new_r = (long ***)getarray(seq_len, sizeof(long **), "new_r");
  delta_new_r = (long ***)getarray(seq_len, sizeof(long **), "delta_new_r");
  subband_flag = (char *)getarray(num_subband, sizeof(char), "subband_flag");
  // manipulate the bitstream shift in different subbands via subband_weighting & shift_per_bitplane
  for (i=0; i<seq_len; i++){
	new_r[i] = (long **)getarray(num_subband, sizeof(long *), "new_r[i]");
	delta_new_r[i] = (long **)getarray(num_subband, sizeof(long *), "delta_new_r[i]");
	for(j=0; j<num_subband; j++){
      new_r[i][j] = (long *)getarray(max_subplane+1+max_subplane_shift, sizeof(long), "new_r[i][j]");
      delta_new_r[i][j] = (long *)getarray(max_subplane+1+max_subplane_shift, sizeof(long), "delta_new_r[i][j]");
	  for(k=0; k<=max_subplane+max_subplane_shift; k++){
	    new_r[i][j][k]=0;
	    delta_new_r[i][j][k]=0;
	  }
	}

    for(m=0; m<num_subband; m++){
	  subband_flag[m]=0;
	}

	for(j=frame_max_subplane[i]; j>=0; j--){
	  for(m=0; m<num_subband; m++){
		if(subband_flag[m]==0 && delta_r[i][m][j] !=0){
		  delta_new_r[i][m][j+max_subplane_shift-(int)(subband_weighting[m]*shift_per_bitplane[i])]=0;
		  subband_flag[m]=1;
		}
		else
		  delta_new_r[i][m][j+max_subplane_shift-(int)(subband_weighting[m]*shift_per_bitplane[i])]=delta_r[i][m][j];
	  }
	}
	for(m=0; m<num_subband; m++){ 
	  delta_new_r[i][m][max_subplane+max_subplane_shift] = subband_header_bytes[i][m];
	}
  }
  free(subband_flag);


  // mask useless subbands according to the subband_mask
  for (i=0; i<seq_len; i++){ // NOT WORKING
	for(m=0; m<num_subband; m++){
	  if(subband_mask[i][m] == 1){
		assert(0);
	    for(j=max_subplane+max_subplane_shift; j>=0; j--){
		  delta_new_r[i][m][j]=0;
		}
	  }
	}
  } // NOT WORKING

  // form the accumulated sum of bits for each subband and save the sum in new_r 
  for(i=0; i<seq_len; i++){	
	for(j=0; j<num_subband; j++){ 
  	  new_r[i][j][max_subplane+max_subplane_shift]=delta_new_r[i][j][max_subplane+max_subplane_shift];
	  for(k=max_subplane-1+max_subplane_shift; k>=0; k--){		
	    new_r[i][j][k] = new_r[i][j][k+1] + delta_new_r[i][j][k]; 
	  }
	}
  }

  for(i=0; i<seq_len; i++){	
	int count = num_subband_bp-1;
	for(k=max_subplane+max_subplane_shift; k>=0; k--){		
	  for(j=0; j<num_subband; j++){ 
	    if(count == num_subband_bp-1){
//			assert(count == max_subplane+max_subplane_shift);
		  scan[i][count] = delta_new_r[i][j][max_subplane+max_subplane_shift];
		}else
	      scan[i][count] = scan[i][count+1] +  delta_new_r[i][j][k];
	    count--;
	  }
	}
  }

  // rate increasing curve is in sum_r
  sum_r = (double *)getarray(num_subband_bp, sizeof(double), "sum_r");


  for (j=0; j<num_subband_bp; j++){
    sum_r[j] = 0;
  }
  for (j=num_subband_bp-1; j>=0; j--){
    for (i=0; i<seq_len; i++){
      sum_r[j] = sum_r[j] + scan[i][j];
	}
  } // loop of i (no_stream)


  rate = ( long * )getarray( seq_len, sizeof( long ), "rate" );
  subband_bytes = ( long ** )getarray( seq_len, sizeof( long * ), "subband_bytes" );

  for( i = 0; i < seq_len; i++ ) {
    subband_bytes[i] =  ( long * )getarray( num_subband, sizeof( long ), "num_subband" );
    for( j = 0; j < num_subband; j++ ) {
      subband_bytes[i][j] = 0;
    }
  }

  ///////////////////////////////////	VBR   //////////////////////////////
  // allocate bytes to each frame
/*
  budget = ( double )info.bitrate * 1000 / 8 * ( info.last - info.start + 1 ) / frame_rate; 
  budget = budget - overhead - sum_mv;

  printf("num_subband_bp = %d, budget = %f\n\n",num_subband_bp,budget);
  
  if( sum_r[num_subband_bp - 1] > budget ) {   // bit plane max_msb
    printf( " bitalloc.c: can not handle this case: budget is too small ( %f ).\n", budget );
    // exit( 1 );
    return 0;
  } else {                      // else1
    if( sum_r[num_subband_bp - 1] == budget ) {
      exact_b = num_subband_bp - 1;
      flag = 0;
    } else {                    // else2
      if( sum_r[0] < budget ) {
        exact_b = 0;
        flag = 0;
        budget = sum_r[0];
        printf( " bitalloc.c: no more data in the bit stream for allocation.\n" );
      } else {                  //else3    // find the interval
        for( j = num_subband_bp - 2; j >= 0; j-- ) {   // the other bit planes    
          if( sum_r[j] > budget ) {
            b1 = j;             // record the index 
            flag = 1;
            break;
          } else {
            if( sum_r[j] == budget ) {
              exact_b = j;
              flag = 0;
              break;
            }
          }
        }                       // for(j)
      }                         // else3
    }                           // else2
  }                             // else1


	if (flag == 0){  // perfect match
	  exact_bp = exact_b / num_subband; 
	  exact_subband = num_subband - 1 - exact_b%num_subband;
 	  for(i=0; i<seq_len; i++){			
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}
		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}
	  }

	}

	// get the interpolated bytes
	if (flag == 1){
      sum_k = 0.;
      s_r = 0.;
      for (j=0; j<seq_len; j++){
        sum_k = sum_k + scan[j][b1+1] - scan[j][b1];
        s_r = s_r + scan[j][b1];
      }
      b = (budget - s_r + b1 * sum_k)/sum_k;
   
      for(j=0; j<seq_len; j++){
        slope = scan[j][b1+1] - scan[j][b1];
        rate[j] = (long)(scan[j][b1] + slope*(b-b1));		
      }
		
	  exact_bp = b1 / num_subband; 
	  exact_subband = num_subband - 1 - b1%num_subband;
  	  for(i=0; i<seq_len; i++){		    
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}

		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}	 
		subband_bytes[i][exact_subband] -= scan[i][b1]-rate[i];
	  }
	} 


  for( i = 0; i < seq_len; i++ ) {
    for( j = 0; j < num_subband; j++ ) {
      // if( subband_bytes[i][j] > 0 ) {
		if( ! subband_mask[i][j] ) {
        subband_bytes[i][j] += 3;
      }
    }
  }
*/
///////////////////////////////	CBR  //////////////////////////////////

  curr = 0;
  last = curr_last - info.start;

  for( i = 0; i < seq_len; i++ ) {
    for( j = 0; j < num_subband; j++ ) {
      subband_bytes[i][j] = 0;
    }
  }

  counter = 0;

  while( curr <= last ){

	GOP_budget = ( double )info.bitrate * 1000 / 8 * ( info.GOPsz ) / frame_rate; 
	GOP_budget = GOP_budget - GOP_overhead - gop_sum_mv[counter];
	printf("GOP_overhead = %d, gop_sum_mv = %ld, GOP_budget = %f\n", GOP_overhead, gop_sum_mv[counter],GOP_budget);
//Part 1
	for (j=0; j<num_subband_bp; j++){
		sum_r[j] = 0;
	}
	for (j=num_subband_bp-1; j>=0; j--){
		for (i=curr; i<curr+info.GOPsz; i++){
			sum_r[j] = sum_r[j] + scan[i][j];
		}
	} // loop of i (no_stream)

//Part 2
	if( sum_r[num_subband_bp - 1] > GOP_budget ) {   // bit plane max_msb
		printf( " bitalloc.c: can not handle this case: GOP_budget is too small ( %f ) for GOP %d.\n", GOP_budget, counter );
//		exit( 1 );
//		return 0;
//Revised on 06.16.2019
		GOP_budget = sum_r[num_subband_bp - 1];
		exact_b = num_subband_bp - 1;
		flag = 0;
	} else {                      // else1
		if( sum_r[num_subband_bp - 1] == GOP_budget ) {
		  exact_b = num_subband_bp - 1;
		  flag = 0;
		} else {                    // else2
		  if( sum_r[0] < GOP_budget ) {
			exact_b = 0;
			flag = 0;
			GOP_budget = sum_r[0];
			printf( " bitalloc.c: no more data in the bit stream for allocation.\n" );
		  } else {                  //else3    // find the interval
			for( j = num_subband_bp - 2; j >= 0; j-- ) {   // the other bit planes    
			  if( sum_r[j] > GOP_budget ) {
				b1 = j;             // record the index 
				flag = 1;
				break;
			  } else {
				if( sum_r[j] == GOP_budget ) {
				  exact_b = j;
				  flag = 0;
				  break;
				}
			  }
			}                       // for(j)
		  }                         // else3
		}                           // else2
	}                             // else1

//Part 3
	if (flag == 0){  // perfect match
//		assert(0);
	  printf("perfect match detected!\n");
	  exact_bp = exact_b / num_subband; 
	  exact_subband = num_subband - 1 - exact_b%num_subband;
 	  for(i=curr; i<curr+info.GOPsz; i++){			
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}
		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}
	  }// i
	}//if flag

//Part 4
	// get the interpolated bytes
	if (flag == 1){
      sum_k = 0.;
      s_r = 0.;
      for (j=curr; j<curr+info.GOPsz; j++){
        sum_k = sum_k + scan[j][b1+1] - scan[j][b1];
        s_r = s_r + scan[j][b1];
      }
      b = (GOP_budget - s_r + b1 * sum_k)/sum_k;
   
      for(j=curr; j<curr+info.GOPsz; j++){
        slope = scan[j][b1+1] - scan[j][b1];
        rate[j] = (long)(scan[j][b1] + slope*(b-b1));		
      }
		
	  exact_bp = b1 / num_subband; 
	  exact_subband = num_subband - 1 - b1%num_subband;
  	  for(i=curr; i<curr+info.GOPsz; i++){		    
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}

		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}	 
		subband_bytes[i][exact_subband] -= scan[i][b1]-rate[i];
	  }
	}

//Part 5
	for( i = curr; i < curr+info.GOPsz; i++ ) {
	  for( j = 0; j < num_subband; j++ ) {
		  if( ! subband_mask[i][j] ) {
          subband_bytes[i][j] += 3;
        }
      }
    }
	
	curr += info.GOPsz;
	counter ++;
  }//while curr <= last



  /**********/
  /* OUTPUT */
  /**********/
  sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, info.bitrate );
  FID = fopen( bit_alloc_name, "wb" );

  fprintf( FID, "%d\n", num_subband );

  for( i = 0; i < seq_len; i++ ) {
    for( m = 0; m < num_subband; m++ ) {
      fprintf( FID, "%ld\n", subband_bytes[i][m] );
    }
  }
  fclose( FID );


  for(i=0; i<seq_len; i++){
	for(j=0; j<num_subband; j++){
	  free(delta_r[i][j]);
	  free(bytes_upto_bp[i][j]);
	  free(new_r[i][j]);
	  free(delta_new_r[i][j]);
	}

	free(subband_bytes[i]);
    free(subband_header_bytes[i]);

	free(delta_r[i]);
	free(new_r[i]);
	free(delta_new_r[i]);
	free(bytes_upto_bp[i]);

	free(scan[i]);
  }
  free(subband_bytes);
  free(delta_r);
  free(new_r);
  free(delta_new_r);
  free(bytes_upto_bp);
  free(rate);
  free(sum_r);
  free(frame_max_bitplane);
  free(frame_max_subplane);
  free(subband_header_bytes);
  free(total_subband_header_bytes);

  free(GOP_num_subband_bp);

  return 1;
}



int bit_alloc_VBR_v2( videoinfo info, unsigned char *qp, long int overhead, unsigned long *gop_sum_mv, int curr_last, int vbr_bound )
{
  enum FRAME_RES { SD2, SD, CIF, HD, INVALID_RES } frame_res;

  long int GOP_overhead;
  long int sum_mv;

  int i, j, k, m, counter, GOPsz, curr, last, ncomps, flag, exact_b = 0, b1,
    bit_plane, skipped_subband, log_my_scale = 0;
  int remaining_frs, frame_rate, max_msb = 0,  nbands[3]; // *GOP_header_bytes
  int num_of_GOP, act_num_of_GOP, seq_len, num_subband, num_subband_bp, seq_len_vbr;

  int *GOP_num_subband_bp;

  int **subband_header_bytes, *total_subband_header_bytes; // all kinds of headers
  int eff_GOPsz[20], GOP_multiplier = 1, dist, s_level;
  long  ***delta_r, ***bytes_upto_bp, ***delta_new_r, ***new_r, **scan, *rate, **subband_bytes;
  double b, *sum_r, budget, GOP_budget, sum_k, s_r, slope;
  char seq_name[250], bitplanename[250];
  char *subband_flag, **subband_mask;
  char bit_alloc_name[256];
  FILE *FID;  
  int qtree_depth = 0;
  int  *frame_max_bitplane, roll_level1, roll_level2; // by Yongjun Wu
  int  max_subplane=0, max_bitplane=0, *frame_max_subplane;
  int  exact_bp, exact_subband; 

  int res_frame, getnum, res_exist, do_res = NO;
  int *level_res_frame, *seq_sort, res_level, level_frame;
  
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 ); 
  seq_name[strlen( info.bitname ) - 4] = '\0'; 

  if (info.org_yheight == 480 && info.org_ywidth == 832 )
	  frame_res = SD;
  else if (info.org_yheight == 288 && info.org_ywidth == 352 )
	  frame_res = CIF;
  else if (info.org_yheight == 1080 && info.org_ywidth == 1920)
	  frame_res = HD;
  else if (info.org_yheight == 576 && info.org_ywidth == 704)
	  frame_res = SD2;

  // now each frame is coded independently  --- by Yongjun Wu 
  // it's 2-D EZBC coder instead of 3-D EZBC coder 
  seq_len    = curr_last - info.start + 1;//For foreman, seq_len = 288

  level_res_frame = (int*)getarray(info.tPyrLev, sizeof(int), "level_res_frame");
  seq_sort = (int*)getarray(seq_len, sizeof(int), "seq_sort");
/*
  if(curr_last <= info.act_last)
	seq_len_vbr = seq_len - info.GOPsz;//For foreman, seq_len_vbr = 272
  else
	seq_len_vbr = info.act_last + 1;
*/  
  if( (info.act_last + 1) % info.GOPsz ){
	res_exist = YES;
  }else{
	res_exist = NO;
  }

  if( res_exist == YES ){
	  curr = (seq_len - info.GOPsz * 2 > 0) ? (seq_len - info.GOPsz * 2) : 0;

	  res_frame = info.act_last - curr + 1;
	  if( res_frame > 0 && res_frame < info.GOPsz ){
		do_res = YES;
		for(i =  0; i < info.tPyrLev ; i ++){
			getnum = res_frame / 2;
			level_res_frame[i] = getnum;
			res_frame = res_frame - getnum;
		}
	  }else{
		do_res = NO;
	  }
  }else{
	  assert(res_exist == NO);
	  do_res = NO;
	  curr = seq_len - info.GOPsz;
  }

  if(do_res == NO){
	  for(i = 0; i < seq_len; i ++){
		  if(i < curr )
			  seq_sort[i] = MAIN;
		  else
			  seq_sort[i] = RES;
	  }
  }else{
	assert(do_res == YES && res_exist == YES);
	for(i = 0; i < curr; i ++)
		seq_sort[i] = MAIN;
/////////////////////
	seq_sort[curr] = MAIN;
	res_frame = 1;
	res_level = info.tPyrLev - 1;
	for(i = curr + 1; i < curr + info.GOPsz; i ++){
		if( (i-curr) == res_frame*2 ){
			res_frame = res_frame*2;
			res_level --;
		}

		if( (i-curr-res_frame) < level_res_frame[res_level] )
			seq_sort[i] = MAIN;
		else
			seq_sort[i] = RES;
	}
////////////////////
	for(i = curr + info.GOPsz; i < seq_len; i ++)
		seq_sort[i] = RES;
  }

  getnum = 0;
  for(i = 0; i <  seq_len; i ++)
	  if(seq_sort[i] == MAIN)
		  getnum ++;
	  else
		assert(seq_sort[i] == RES);

  seq_len_vbr = getnum;

  printf("res_exist = %d, do_res = %d, main frame number: %d, res frame number: %d\n", res_exist, do_res, seq_len_vbr, seq_len-seq_len_vbr);

/////////////////////////////////////////////////////////////////////////////////
//  seq_len_vbr = seq_len;
  frame_rate = info.framerate;
  num_of_GOP = get_GOP_num( info );

  act_num_of_GOP = (long int) ( (seq_len_vbr - 1 - info.start + info.GOPsz) / info.GOPsz );

  if( info.denoise_flag == YES ) {
    GOP_multiplier = ( YUV420 == 1 ) ? 2 : 4;
    seq_len *= GOP_multiplier;
  }

  printf("seq_len = %d\n",seq_len);

  // maximum bitplane number in each frame by Yongjun Wu
  frame_max_bitplane = (int *)getarray(seq_len, sizeof(int), "frame_max_bitplane");  
  // maximum sub-bitplane number in each frame 
  frame_max_subplane = (int *)getarray(seq_len, sizeof(int), "frame_max_subplane");
  // delta_r[i][j][k]: each data item, [i]: frame number, [j]: subband number, [k]: sub-bitplane number 
  delta_r = ( long *** )getarray( seq_len, sizeof( long * ), "delta_r" );
  // accumulated data items 
  bytes_upto_bp = ( long *** )getarray( seq_len, sizeof( long * ), "bytes_upto_bp" );

  subband_header_bytes = ( int ** )getarray( seq_len, sizeof( int * ), "subband_header_bytes" );
  total_subband_header_bytes = ( int * )getarray( seq_len, sizeof( int ), "total_subband_header_bytes" );
  scan                 = (long **)getarray(seq_len, sizeof(long *), "scan");

  GOP_num_subband_bp = (int *)getarray(num_of_GOP, sizeof( int * ), "GOP_num_subband_bp" );

  sprintf( bitplanename, "%s.rd_sample_dat", seq_name );
  bitplanename[strlen( seq_name ) + 15] = '\0';  //
  FID = fopen( bitplanename, "rb" );
  if( FID == NULL ) {
    printf( "can not open file bit_rd_sample.dat\n" );
    exit( 1 );
  }

  ncomps = 1;
  num_subband = 0;
  for( m = 0; m < ncomps; m++ ) {
    // read in the number of substreams of each frame    
    fscanf( FID, "%d\n", &( nbands[m] ) );
#ifdef   	FREQUENCY_ROLL_OFF
	if (frame_res == HD ) //Added by Yuan Liu
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res == SD2)
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif

		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res == SD )
	{
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2;     // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif

		roll_level1 = 5; roll_level2=4;   // the subband levels where we need frequency roll-off
	}
	else if (frame_res==CIF)
	{
		num_subband += nbands[m]+2;    // the independent coded bitstreams 
		roll_level1=roll_level2=4;
	}
	else printf("Error in frequency roll-off!\n"); 

#else
    num_subband += nbands[m];      // the independent coded bitstreams 
#endif
  }
  fseek( FID, 0, SEEK_SET );

  for(i = 0; i < seq_len; i++){
	delta_r[i]       = (long **)getarray(num_subband, sizeof(long *), "delta_r[i]");
    bytes_upto_bp[i] = (long **)getarray(num_subband, sizeof(long *), "bytes_upto_bp[i]");
	subband_header_bytes[i] = (int *)getarray(num_subband, sizeof(int), "subband_header_bytes[i]");
  }

  // generate subband mask 
  subband_mask = ( char ** )getarray( seq_len, sizeof( char * ), "subband_mask(bitalloc.c)" );
  for( i = 0; i < seq_len; i++ ) {
    subband_mask[i] =  ( char * )getarray( num_subband, sizeof( char ), "subband_mask[i] (bitalloc.c)" );
    for( k = 0; k < num_subband; k++ ) {
      subband_mask[i][k] = 0;
    }
  }

  curr = info.start;
  last = curr_last;
  counter = 0;
  skipped_subband = 0;
  while( curr <= last ) {
    remaining_frs = last - curr + 1;  
    // determine effective GOP size in level i
    GOPsz = info.GOPsz;
    for( i = 0; i <= info.tPyrLev; i++ ) {
      if ( remaining_frs < info.GOPsz )
        eff_GOPsz[i] = (int) ceil ((double) remaining_frs / (pow (2, i)));  
      else
        eff_GOPsz[i] = GOPsz;  
      GOPsz /= 2;  
    }
  
    //mask useless temporal high subband frames  --- now each frame is coded independently
    GOPsz = info.GOPsz;
    for( i = 0; i < info.t_level; i++ ) {
      for( j = eff_GOPsz[i+1]; j < eff_GOPsz[i]; j++ ) {
        for( k = 0; k < num_subband; k++ ) {
		  assert(0);
          subband_mask[counter + j][k] = 1;
          skipped_subband++;
        }
      }
    }

    //mask useless spatial high subband frames; if current GOP < default GOP size, GOPsz could be 0
    for( i = 0; i < MY_MAX( eff_GOPsz[0], 1 ); i++ ) {
      s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
      ncomps = 1; // since Y U V are grouped together
      m = 0;
      for( j = 0; j < ncomps; j++ ) {
        m += nbands[j];

#ifdef   	FREQUENCY_ROLL_OFF
		int subband_index;
		if (frame_res == HD )
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
	    else if (frame_res == SD2)
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
		else if (frame_res == SD )
#ifdef ROLL_STRUCTURE_ONE
			subband_index = nbands[j]+2;   // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
			subband_index = nbands[j]+2+11;  // the independent coded bitstreams 
#endif
		else if (frame_res==CIF)
			subband_index = nbands[j]+2;    // the independent coded bitstreams 
		else printf("Error in frequency roll-off!\n"); 
#endif

      }
    }

    curr += info.GOPsz;
    counter += info.GOPsz * GOP_multiplier;
  }                             /* while */
  
  // *** overhead bytes for subband allocation table ***
  // sizeof(unsigned short int) is for header, 3 is the extra bytes at the end of each substream

//  printf("skipped_subband = %d\n",skipped_subband);

  overhead += ( sizeof( unsigned short int ) + 3 ) * ( num_subband * (info.last - info.start + 1) - skipped_subband );

  GOP_overhead = overhead / num_of_GOP;

  for( i = 0; i < seq_len; i++ ) {
    num_subband = 0;
    for( m = 0; m < ncomps; m++ ) {
      // read in the number of subbands of each component          
      fscanf( FID, "%d\n", &( nbands[m] ) );
#ifdef   	FREQUENCY_ROLL_OFF
	if (frame_res == HD )

#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;    // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res == SD2)

#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2+2;    // the independent coded bitstreams 
#endif

#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res == SD )
#ifdef ROLL_STRUCTURE_ONE
		num_subband += nbands[m]+2;    // the independent coded bitstreams 
#endif
#ifdef ROLL_STRUCTURE_TWO
		num_subband += nbands[m]+2+11;    // the independent coded bitstreams 
#endif
	else if (frame_res==CIF)
		num_subband += nbands[m]+2;      // the independent coded bitstreams 
	else printf("Error in frequency roll-off!\n"); 
#else
    num_subband += nbands[m];           // the independent coded bitstreams 
#endif

    }

	// read in the maximum bit plane of a frame by Yongjun Wu
	fscanf(FID, "%d\n", &(frame_max_bitplane[i]));    
	if(frame_max_bitplane[i] > max_bitplane) max_bitplane = frame_max_bitplane[i];
    // read in the maximum sub-bitplane of a frame
    fscanf( FID, "%d\n", &( frame_max_subplane[i] ) );
    if( frame_max_subplane[i] > max_subplane ) max_subplane = frame_max_subplane[i];

	// scan bit planes of each subband
	for(j=0; j<num_subband; j++){
      delta_r[i][j] = (long *)getarray(frame_max_subplane[i]+1, sizeof(long), "delta_r[i][j]");
      bytes_upto_bp[i][j] = (long *)getarray(frame_max_subplane[i]+1, sizeof(long), "bytes_upto_bp[i][j]");	
	}

    for(j=0; j<num_subband; j++){
		subband_header_bytes[i][j] = 0;
	}
	total_subband_header_bytes[i] = 0;

	for(j=frame_max_subplane[i]; j>=0; j--){

	  for(m=0; m<num_subband; m++){
		//read in bitplane index
        fscanf(FID, "%d", &bit_plane);      
        if (j != bit_plane){
          printf(" j %d bit_plane %d bit plane index error !\n", j, bit_plane);
		  exit(1);
		}		 
	    //read in the size of each subband bitplane
		fscanf(FID, "%d %d\n", &(delta_r[i][m][j]), &(bytes_upto_bp[i][m][j]));    
			
		// subband head bytes are located in the first place
		if(subband_header_bytes[i][m]==0 && delta_r[i][m][j] !=0){
			subband_header_bytes[i][m] = delta_r[i][m][j];
		    total_subband_header_bytes[i] +=  subband_header_bytes[i][m];
		}
	  } //m
	}//j
  }   // loop of i 
  fclose( FID );


  
  // the set of variables for the idea of frequency roll-off 
  int  *shift_per_bitplane, max_per_shift =0;  
  int   max_bitplane_shift, max_subplane_shift;
  float  *subband_weighting  = (float *)getarray(num_subband, sizeof(float), "subband_weighting");
  // initialize all the varialbles to be 0
  memset(subband_weighting, 0, num_subband*sizeof(float));
  shift_per_bitplane = (int *)getarray(seq_len, sizeof(int), "shift_per_bitplane");
  memset(shift_per_bitplane, 0, seq_len*sizeof(int));
  max_subplane_shift = max_per_shift = max_bitplane_shift = 0; 

#ifdef FREQUENCY_ROLL_OFF
  FILE *froll;
  float roll_factors[15]; 
  int   factor_num; 

  for (i=0; i<seq_len; i++)
  {
	  shift_per_bitplane[i] = (frame_max_subplane[i]+1)/(frame_max_bitplane[i]+1); 

	  if (max_per_shift < shift_per_bitplane[i] )
		max_per_shift = shift_per_bitplane[i];
  }
  max_bitplane_shift = 5; 
  max_subplane_shift = max_bitplane_shift * max_per_shift; 

#endif

  // for interleaving use 
  printf("num_subband = %d, max_subplane = %d, max_subplane_shift = %d\n",num_subband, max_subplane, max_subplane_shift);
  num_subband_bp = num_subband * (max_subplane+1+max_subplane_shift); 
  for (i=0; i<seq_len; i++)
    scan[i] = (long *)getarray(num_subband_bp, sizeof(long), "scan[i]");

  // move data in delta_r to delta_new_r to make all substream look like having the same number of bitplanes  
  // in order to simplify the following bit allocation algorithm
  new_r = (long ***)getarray(seq_len, sizeof(long **), "new_r");
  delta_new_r = (long ***)getarray(seq_len, sizeof(long **), "delta_new_r");
  subband_flag = (char *)getarray(num_subband, sizeof(char), "subband_flag");
  // manipulate the bitstream shift in different subbands via subband_weighting & shift_per_bitplane
  for (i=0; i<seq_len; i++){
	new_r[i] = (long **)getarray(num_subband, sizeof(long *), "new_r[i]");
	delta_new_r[i] = (long **)getarray(num_subband, sizeof(long *), "delta_new_r[i]");
	for(j=0; j<num_subband; j++){
      new_r[i][j] = (long *)getarray(max_subplane+1+max_subplane_shift, sizeof(long), "new_r[i][j]");
      delta_new_r[i][j] = (long *)getarray(max_subplane+1+max_subplane_shift, sizeof(long), "delta_new_r[i][j]");
	  for(k=0; k<=max_subplane+max_subplane_shift; k++){
	    new_r[i][j][k]=0;
	    delta_new_r[i][j][k]=0;
	  }
	}

    for(m=0; m<num_subband; m++){
	  subband_flag[m]=0;
	}

	for(j=frame_max_subplane[i]; j>=0; j--){
	  for(m=0; m<num_subband; m++){
		if(subband_flag[m]==0 && delta_r[i][m][j] !=0){
		  delta_new_r[i][m][j+max_subplane_shift-(int)(subband_weighting[m]*shift_per_bitplane[i])]=0;
		  subband_flag[m]=1;
		}
		else
		  delta_new_r[i][m][j+max_subplane_shift-(int)(subband_weighting[m]*shift_per_bitplane[i])]=delta_r[i][m][j];
	  }
	}
	for(m=0; m<num_subband; m++){ 
	  delta_new_r[i][m][max_subplane+max_subplane_shift] = subband_header_bytes[i][m];
	}
  }
  free(subband_flag);

  // mask useless subbands according to the subband_mask
  for (i=0; i<seq_len; i++){ // NOT WORKING
	for(m=0; m<num_subband; m++){
	  if(subband_mask[i][m] == 1){
		assert(0);
	    for(j=max_subplane+max_subplane_shift; j>=0; j--){
		  delta_new_r[i][m][j]=0;
		}
	  }
	}
  } // NOT WORKING

  // form the accumulated sum of bits for each subband and save the sum in new_r 
  for(i=0; i<seq_len; i++){	
	for(j=0; j<num_subband; j++){ 
  	  new_r[i][j][max_subplane+max_subplane_shift]=delta_new_r[i][j][max_subplane+max_subplane_shift];
	  for(k=max_subplane-1+max_subplane_shift; k>=0; k--){		
	    new_r[i][j][k] = new_r[i][j][k+1] + delta_new_r[i][j][k]; 
	  }
	}
  }

  for(i=0; i<seq_len; i++){	
	int count = num_subband_bp-1;
	for(k=max_subplane+max_subplane_shift; k>=0; k--){		
	  for(j=0; j<num_subband; j++){ 
	    if(count == num_subband_bp-1){
//			assert(count == max_subplane+max_subplane_shift);
		  scan[i][count] = delta_new_r[i][j][max_subplane+max_subplane_shift];
		}else
	      scan[i][count] = scan[i][count+1] +  delta_new_r[i][j][k];
	    count--;
	  }
	}
  }

// rate increasing curve is in sum_r
  sum_r = (double *)getarray(num_subband_bp, sizeof(double), "sum_r");

  for (j=0; j<num_subband_bp; j++){
    sum_r[j] = 0;
  }
  for (j=num_subband_bp-1; j>=0; j--){
    for (i=0; i<seq_len; i++){
      sum_r[j] = sum_r[j] + scan[i][j];
	}
  } // loop of i (no_stream)

  rate = ( long * )getarray( seq_len, sizeof( long ), "rate" );
  subband_bytes = ( long ** )getarray( seq_len, sizeof( long * ), "subband_bytes" );

  for( i = 0; i < seq_len; i++ ) {
    subband_bytes[i] =  ( long * )getarray( num_subband, sizeof( long ), "num_subband" );
    for( j = 0; j < num_subband; j++ ) {
      subband_bytes[i][j] = 0;
    }
  }

///////////////////////////////////	VBR   //////////////////////////////
// allocate bytes to each frame
//Step Main
  sum_mv = 0;
  for(i = 0; i < act_num_of_GOP; i ++){
	  sum_mv += gop_sum_mv[i];
  }

  budget = ( double )info.bitrate * 1000 / 8 * seq_len_vbr / frame_rate; 
  budget = budget - GOP_overhead * act_num_of_GOP - sum_mv;

  printf("Main overhead = %d, sum_mv = %d, budget = %f\n\n",GOP_overhead * act_num_of_GOP,sum_mv,budget);
//Part 1
	for (j=0; j<num_subband_bp; j++){
		sum_r[j] = 0;
	}
	for (j=num_subband_bp-1; j>=0; j--){
		for( i = 0; i < seq_len; i++ ) {
			if(seq_sort[i] == MAIN)
				sum_r[j] = sum_r[j] + scan[i][j];
		}
	} // loop of i (no_stream)

//Part 2
	if( sum_r[num_subband_bp - 1] > budget ) {   // bit plane max_msb
		printf( " bitalloc.c: can not handle this case: budget is too small ( %f ) for GOP %d.\n", budget, 0 );
		// exit( 1 );
		return 0;
	} else {                      // else1
		if( sum_r[num_subband_bp - 1] == budget ) {
		  exact_b = num_subband_bp - 1;
		  flag = 0;
		} else {                    // else2
		  if( sum_r[0] < budget ) {
			exact_b = 0;
			flag = 0;
			budget = sum_r[0];
			printf( " bitalloc.c: no more data in the bit stream for allocation.\n" );
		  } else {                  //else3    // find the interval
			for( j = num_subband_bp - 2; j >= 0; j-- ) {   // the other bit planes    
			  if( sum_r[j] > budget ) {
				b1 = j;             // record the index 
				flag = 1;
				break;
			  } else {
				if( sum_r[j] == budget ) {
				  exact_b = j;
				  flag = 0;
				  break;
				}
			  }
			}                       // for(j)
		  }                         // else3
		}                           // else2
	}                             // else1

//Part 3
	if (flag == 0){  // perfect match
//		assert(0);
	  printf("perfect match detected!\n");
	  exact_bp = exact_b / num_subband; 
	  exact_subband = num_subband - 1 - exact_b%num_subband;
 	  for( i = 0; i < seq_len; i++ ) {	
		  if(seq_sort[i] == MAIN){
			for(m=0; m<=exact_subband; m++){
			  subband_bytes[i][m] = new_r[i][m][exact_bp];
			}
			if(exact_bp < max_subplane){
			  for(m=exact_subband+1; m<num_subband; m++)
				subband_bytes[i][m] = new_r[i][m][exact_bp+1];
			}
		  }
	  }// i
	}//if flag
//	printf("out here!\n");
//Part 4
	// get the interpolated bytes
	if (flag == 1){
      sum_k = 0.;
      s_r = 0.;
      for( j = 0; j < seq_len; j++ ) {
		if(seq_sort[j] == MAIN){
			sum_k = sum_k + scan[j][b1+1] - scan[j][b1];
			s_r = s_r + scan[j][b1];
		}
      }
      b = (budget - s_r + b1 * sum_k)/sum_k;
   
      for( j = 0; j < seq_len; j++ ) {
		  if(seq_sort[j] == MAIN){
			slope = scan[j][b1+1] - scan[j][b1];
			rate[j] = (long)(scan[j][b1] + slope*(b-b1));	
		  }
      }
		
	  exact_bp = b1 / num_subband; 
	  exact_subband = num_subband - 1 - b1%num_subband;
  	  for( i = 0; i < seq_len; i++ ) {	
		  if(seq_sort[i] == MAIN){
			for(m=0; m<=exact_subband; m++){
			  subband_bytes[i][m] = new_r[i][m][exact_bp];
			}

			if(exact_bp < max_subplane){
			  for(m=exact_subband+1; m<num_subband; m++)
				subband_bytes[i][m] = new_r[i][m][exact_bp+1];
			}	 
			subband_bytes[i][exact_subband] -= scan[i][b1]-rate[i];
		  }
	  }
	}

//Part 5
	for( i = 0; i < seq_len; i++ ) {
		if(seq_sort[i] == MAIN){
		  for( j = 0; j < num_subband; j++ ) {
			  if( ! subband_mask[i][j] ) {
			  subband_bytes[i][j] += 3;
			}
		  }
		}
    }
	
//////////////////////////////////

//Step Residual
  sum_mv = gop_sum_mv[num_of_GOP-1];

  GOP_budget = ( double )info.bitrate * 1000 / 8 * (seq_len - seq_len_vbr) / frame_rate; 
  GOP_budget = GOP_budget - GOP_overhead - sum_mv;

  printf("Residual overhead = %d, sum_mv = %d, GOP_budget = %f\n\n",GOP_overhead,sum_mv,GOP_budget);
//Part 1
	for (j=0; j<num_subband_bp; j++){
		sum_r[j] = 0;
	}
	for (j=num_subband_bp-1; j>=0; j--){
		for( i = 0; i < seq_len; i++ ) {
			if(seq_sort[i] == RES)
				sum_r[j] = sum_r[j] + scan[i][j];
		}
	} // loop of i (no_stream)

//Part 2
	if( sum_r[num_subband_bp - 1] > GOP_budget ) {   // bit plane max_msb
		printf( " bitalloc.c: can not handle this case: GOP_budget is too small ( %f ) for GOP %d.\n", GOP_budget, counter );
		// exit( 1 );
		return 0;
	} else {                      // else1
		if( sum_r[num_subband_bp - 1] == GOP_budget ) {
		  exact_b = num_subband_bp - 1;
		  flag = 0;
		} else {                    // else2
		  if( sum_r[0] < GOP_budget ) {
			exact_b = 0;
			flag = 0;
			GOP_budget = sum_r[0];
			printf( " bitalloc.c: no more data in the bit stream for allocation.\n" );
		  } else {                  //else3    // find the interval
			for( j = num_subband_bp - 2; j >= 0; j-- ) {   // the other bit planes    
			  if( sum_r[j] > GOP_budget ) {
				b1 = j;             // record the index 
				flag = 1;
				break;
			  } else {
				if( sum_r[j] == GOP_budget ) {
				  exact_b = j;
				  flag = 0;
				  break;
				}
			  }
			}                       // for(j)
		  }                         // else3
		}                           // else2
	}                             // else1

//Part 3
	if (flag == 0){  // perfect match
//		assert(0);
	  printf("perfect match detected!\n");
	  exact_bp = exact_b / num_subband; 
	  exact_subband = num_subband - 1 - exact_b%num_subband;
 	  for( i = 0; i < seq_len; i++ ) {	
		  if(seq_sort[i] == RES){
			for(m=0; m<=exact_subband; m++){
			  subband_bytes[i][m] = new_r[i][m][exact_bp];
			}
			if(exact_bp < max_subplane){
			  for(m=exact_subband+1; m<num_subband; m++)
				subband_bytes[i][m] = new_r[i][m][exact_bp+1];
			}
		  }
	  }// i
	}//if flag

//Part 4
	// get the interpolated bytes
	if (flag == 1){
      sum_k = 0.;
      s_r = 0.;
      for( j = 0; j < seq_len; j++ ) {
		  if(seq_sort[j] == RES){
			sum_k = sum_k + scan[j][b1+1] - scan[j][b1];
			s_r = s_r + scan[j][b1];
		  }
      }
      b = (GOP_budget - s_r + b1 * sum_k)/sum_k;
   
      for( j = 0; j < seq_len; j++ ) {
		  if(seq_sort[j] == RES){
			slope = scan[j][b1+1] - scan[j][b1];
			rate[j] = (long)(scan[j][b1] + slope*(b-b1));	
		  }
      }
		
	  exact_bp = b1 / num_subband; 
	  exact_subband = num_subband - 1 - b1%num_subband;
  	  for( i = 0; i < seq_len; i++ ) {	  
		  if(seq_sort[i] == RES){
			for(m=0; m<=exact_subband; m++){
			  subband_bytes[i][m] = new_r[i][m][exact_bp];
			}

			if(exact_bp < max_subplane){
			  for(m=exact_subband+1; m<num_subband; m++)
				subband_bytes[i][m] = new_r[i][m][exact_bp+1];
			}	 
			subband_bytes[i][exact_subband] -= scan[i][b1]-rate[i];
		  }
	  }
	}

//Part 5
	for( i = 0; i < seq_len; i++ ) {
		if(seq_sort[i] == RES){
		  for( j = 0; j < num_subband; j++ ) {
			  if( ! subband_mask[i][j] ) {
			  subband_bytes[i][j] += 3;
			}
		  }
		}
    }

///////////////	FOR ENCODER	/////////////////////////////////////
/*
// allocate bytes to each frame
  budget = ( double )info.bitrate * 1000 / 8 * ( info.last - info.start + 1 ) / frame_rate; 
  budget = budget - overhead - sum_mv;

  printf("overhead = %d, sum_mv = %d, budget = %f\n\n",overhead,sum_mv,budget);
  
  if( sum_r[num_subband_bp - 1] > budget ) {   // bit plane max_msb
    printf( " bitalloc.c: can not handle this case: budget is too small ( %f ).\n", budget );
    // exit( 1 );
    return 0;
  } else {                      // else1
    if( sum_r[num_subband_bp - 1] == budget ) {
      exact_b = num_subband_bp - 1;
      flag = 0;
    } else {                    // else2
      if( sum_r[0] < budget ) {
        exact_b = 0;
        flag = 0;
        budget = sum_r[0];
        printf( " bitalloc.c: no more data in the bit stream for allocation.\n" );
      } else {                  //else3    // find the interval
        for( j = num_subband_bp - 2; j >= 0; j-- ) {   // the other bit planes    
          if( sum_r[j] > budget ) {
            b1 = j;             // record the index 
            flag = 1;
            break;
          } else {
            if( sum_r[j] == budget ) {
              exact_b = j;
              flag = 0;
              break;
            }
          }
        }                       // for(j)
      }                         // else3
    }                           // else2
  }                             // else1


	if (flag == 0){  // perfect match
	  exact_bp = exact_b / num_subband; 
	  exact_subband = num_subband - 1 - exact_b%num_subband;
 	  for(i=0; i<seq_len; i++){			
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}
		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}
	  }

	}

	// get the interpolated bytes
	if (flag == 1){
      sum_k = 0.;
      s_r = 0.;
      for (j=0; j<seq_len; j++){
        sum_k = sum_k + scan[j][b1+1] - scan[j][b1];
        s_r = s_r + scan[j][b1];
      }
      b = (budget - s_r + b1 * sum_k)/sum_k;
   
      for(j=0; j<seq_len; j++){
        slope = scan[j][b1+1] - scan[j][b1];
        rate[j] = (long)(scan[j][b1] + slope*(b-b1));		
      }
		
	  exact_bp = b1 / num_subband; 
	  exact_subband = num_subband - 1 - b1%num_subband;
  	  for(i=0; i<seq_len; i++){		    
        for(m=0; m<=exact_subband; m++){
          subband_bytes[i][m] = new_r[i][m][exact_bp];
		}

		if(exact_bp < max_subplane){
          for(m=exact_subband+1; m<num_subband; m++)
		    subband_bytes[i][m] = new_r[i][m][exact_bp+1];
		}	 
		subband_bytes[i][exact_subband] -= scan[i][b1]-rate[i];
	  }
	} 


  for( i = 0; i < seq_len; i++ ) {
    for( j = 0; j < num_subband; j++ ) {
      // if( subband_bytes[i][j] > 0 ) {
		if( ! subband_mask[i][j] ) {
        subband_bytes[i][j] += 3;
      }
    }
  }
*/

  /**********/
  /* OUTPUT */
  /**********/
  sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, info.bitrate );
  FID = fopen( bit_alloc_name, "wb" );

  fprintf( FID, "%d\n", num_subband );

  for( i = 0; i < seq_len; i++ ) {
    for( m = 0; m < num_subband; m++ ) {
      fprintf( FID, "%ld\n", subband_bytes[i][m] );
    }
  }
  fclose( FID );


  for(i=0; i<seq_len; i++){
	for(j=0; j<num_subband; j++){
	  free(delta_r[i][j]);
	  free(bytes_upto_bp[i][j]);
	  free(new_r[i][j]);
	  free(delta_new_r[i][j]);
	}

	free(subband_bytes[i]);
    free(subband_header_bytes[i]);

	free(delta_r[i]);
	free(new_r[i]);
	free(delta_new_r[i]);
	free(bytes_upto_bp[i]);

	free(scan[i]);
  }
  free(subband_bytes);
  free(delta_r);
  free(new_r);
  free(delta_new_r);
  free(bytes_upto_bp);
  free(rate);
  free(sum_r);
  free(frame_max_bitplane);
  free(frame_max_subplane);
  free(subband_header_bytes);
  free(total_subband_header_bytes);

  free(seq_sort);
  free(level_res_frame);
  free(GOP_num_subband_bp);

  return 1;
}