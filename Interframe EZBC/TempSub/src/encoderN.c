#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define  EXTERN 
#include "structN.h"
#include "basic.h"
#include "coderN.h"
#include "init_encN.h"
#include "memoryN.h"
#include "analsyn.h"
#include "ioN.h"
#include "pstatN.h"
#include "miscN.h"
#include "ezbc_enc_3d.h"
#include "mv_statistics.h"
#include "hvsbm_fullres.h"
#include "bmeN.h"

#include "layer_mv.h"

#define EXTERN extern
#define LAMBDA_ADAPT_RATIO_THRES 0.6
#define INIT_RANGE 7
#define LAMBDA_LOWER_LMT 3

#define THEO_LENGTH 48
#define THEO_COUNT 6
#define ADAPT_LEVEL 5
//#define THEO_MODEL

//#define LAMBDA_CTRL

//Added on 02.05.2018
EXTERN int SIMUL_START;
EXTERN int SIMUL_POINT;
EXTERN int SIMUL_LIMIT;
EXTERN int SIMUL_LEVEL;
EXTERN int SIMUL_RANGE;

EXTERN int SIMUL_ROUND;

unsigned long int *gop_mv;

float global_motion_active;

void read_command( int argc, char **argv, videoinfo * info );

int bit_alloc_CBR( videoinfo info, unsigned char *qp,
				  long int overhead, unsigned long *gop_sum_mv, int curr_last );


typedef struct lambda_rate{
	int bitrate;
	float lambda;
	float bst_psnr;
	char  bitfile[256];
	char  dec_video[256];
	long int total_bytes_past;
	long int total_bytes_past_buffer;
};

typedef struct lambda_buff{
	float buff_lambda[ADAPT_LEVEL];
	lambda_buff *next;
};

void linear_fitting(int maxn, float *a0, float *a1, double *x, double *y){

//    double x[maxn] = {-0.693,-0.511,-0.357,-0.223,-0.105};
//	double y[maxn] = {2.565, 2.485, 2.302, 2.197, 2.197};
    double xi = 0, x2 = 0, yi = 0, xy = 0;

    for(int i = 0; i < maxn; i++){
        xi += x[i], x2 += x[i] * x[i], yi += y[i], xy += x[i] * y[i];
    }
    *a0 = (yi * x2 - xy * xi) / (x2 * maxn - xi * xi);
    *a1 = (yi * xi - xy * maxn) / (xi * xi - x2 * maxn);
    printf("P(x) = %f%+fx\n", *a0, *a1);

}

void file_copy(char *dest_name, char *source_name){

	FILE * file1,*file2;  
     
    file1 = fopen(source_name,"rb");   
    file2 = fopen(dest_name,"wb");   
    if(!file1)  
    {  
        printf("File open failed！",source_name);  
        return;  
    }  
    char c;  
    long int index = 0;  
    fseek(file1,0,SEEK_END);          
    long int length = ftell(file1);      
    //printf("%d\n",length);         
    if(!length)  
        return;  
    while(!fseek(file1,index,SEEK_SET)) 
    {  
        fread(&c,1,1,file1);              
        fwrite(&c,1,1,file2);              
        if(index == length - 1)           
        {  
            break;  
        }  
        index++;                         
    }  
    fclose(file1);                       
    fclose(file2);                       
}

int compare_lambda_vector(lambda_buff *buff, float *lambda_vector){
	int i,j;
	int getnum = 0;
	int judge;
	lambda_buff *root, *create;

	assert(buff != NULL);
	root = buff;

	printf("\nCompare lambda vector:\n");

	while(root != NULL){
		judge = YES;

		for(i = 0; i < ADAPT_LEVEL; i ++)
			printf("%f\t",root->buff_lambda[i]);
		printf("\n");

		for(i = 0; i < ADAPT_LEVEL; i ++){
			if(root->buff_lambda[i] != lambda_vector[i]){
				judge = NO;
				break;
			}
		}
		if(judge == YES)//find a match in the queue!
			return 1;

		root = root->next;
	}

	assert(root == NULL);//if we arrive at this point, a new vector has been detected.

	root = buff;
	while(root->next != NULL){
		root = root->next;
	}

	assert(root->next == NULL);

	create = new lambda_buff;
	for(i = 0; i < ADAPT_LEVEL; i ++)
		create->buff_lambda[i] = lambda_vector[i];
	create->next = NULL;
	root->next = create;

	return 0;
}

void destroy_lambda_buff(lambda_buff *buff){
	lambda_buff *create;

	if(buff != NULL){
		do{
			create = buff;
			buff = create->next;
			delete(create);
		}while(buff->next != NULL);

		delete(buff);
	}
}


/*
 *                                main()                                    
 */
int
main( int argc, char **argv )
{
  int curr, curr_disp, last, i, j, remaining_frs, header_bytes; 
  long int num_of_GOP, *output_GOP_bytes, *dec_output_GOP_bytes;
  int GOP_counter, VBR_GOP_counter;
  long mark, elp;   // initial and elapsed time
  double duration;
  enum FLAG first_GOP;
  enum FLAG Level_change;
  FILE *fpio, *fplmbd;
  char mvstatname[512];
  char post_stream[256], seq_name[256], raw_buffer[256], enc_buff[256], dec_buff[256]; //Added on 01.14.2018

  char rd_name[256], rd_buff[256], mvst_buff[256];

  int estimated_overhead = 0, FAT_bytes, GOPheader_bytes = 2;
  long int sum_mv;

//for simul dec, added on 01.04.2018
  long int total_bytes_past = 0, total_bytes_past_buffer, vbr_total_bytes_past;
  long int *read_GOP_bytes;
//for simul dec


  float best_lambda[ADAPT_LEVEL];
  float curr_lambda[ADAPT_LEVEL];

  float buff_best_lambda[ADAPT_LEVEL];

  int simul_count = 0;
  int round = 0;
  int simul_enc = NO;

  int simul_alert = NO, in_simul = NO;
  int segment_init = YES;

  float simul_psnr;

  float best_simul_psnr[3];
  float theo_simul_psnr[THEO_COUNT];
  double theo_simul_psnr1[THEO_COUNT][INIT_RANGE], theo_simul_psnr2[THEO_COUNT][INIT_RANGE], theo_simul_psnr3[THEO_COUNT][INIT_RANGE+2];
  int getnum;

  enum FLAG VBR_first_GOP, VBR_Level_change;
  int VBR_remaining_frs, vbr_curr, vbr_remaining_frs;

  int BOTTOM_BITRATE;
  double theo_a[3], theo_b[3], theo_log_a[3], theo_log_b[3];//Added on 05.20.2019, lambda = a * R^b
  int maxn = 5;

//Added on 04.07.2019
  lambda_rate lr_record[6];

  lambda_buff *lmda_buff;
  lmda_buff = NULL;

  for(i = 0; i < THEO_COUNT; i ++){
	  lr_record[i].total_bytes_past = 0;
	  lr_record[i].total_bytes_past_buffer = 0;
  }

//Added on 03.26.2018
  int uphor, upver, scale;
  float *upframe1, *upframe2, *upsamp_x;

  int BUFF_SIMUL_POINT, BUFF_SIMUL_START;

  int cx, cy, cpos;
  float block_diff, block_ratio;  //Added on 06.29.2018
  float buff_gop_block, last_gop_diff;
  simul_skip = NO;
  skip_frame = NO;

  int gop_end = NO;

  for(i = 0;i < 3;i ++)
	  best_simul_psnr[i] = 0;

  float buffer_active = 0;
  int SIMUL_ADD;

  printf("Built on %s at %s.\n", __DATE__, __TIME__);

  videoinfo info;
  //  videoinfo info2;              // only used in info.denoise_flag == YES this case

  mark = clock(  );
  read_command( argc, argv, &info ); // 读命令行

  SIMUL_START = 4 * info.GOPsz;
  SIMUL_POINT = info.GOPsz;
  SIMUL_LIMIT = INIT_RANGE;
  SIMUL_LEVEL = 0;
  SIMUL_ROUND = 2;
  SIMUL_RANGE = SIMUL_LIMIT / 2;

  printf("Simulative bitrate = %d\n",info.simul_rate);

  /* open MV statistics file */
  strcpy(mvstatname, info.bitname);
  strcat(mvstatname, "_mvstatistics.log");
  mvStat_open(mvstatname);

  buff_frameMEinfo = (ImageMEinfo *)getarray( info.yheight * info.ywidth, sizeof( ImageMEinfo ), "buff_frameMEinfo" ); 

  for( cy = 0; cy < info.yheight; cy++ ){
	for( cx = 0; cx < info.ywidth; cx++ ){
		cpos = cy * info.ywidth + cx;  
		buff_frameMEinfo[cpos].left_mvx = 0;
		buff_frameMEinfo[cpos].left_mvy = 0;
		buff_frameMEinfo[cpos].right_mvx = 0;
		buff_frameMEinfo[cpos].right_mvy = 0;
	}
  }

//bit file buffer
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 );
  seq_name[strlen( info.bitname ) - 4] = '\0';
  sprintf( enc_buff, "%s_%d.bit", seq_name, 2 );

//subband sample file buffer
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 );  // Hanke, 16.09.02
  seq_name[strlen( info.bitname ) - 4] = '\0';               //
  sprintf( rd_name, "%s.rd_sample_dat", seq_name ); 
  sprintf( rd_buff, "%s_buff.rd_sample_dat", seq_name );

//mvby file buffer
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 );
  seq_name[strlen( info.bitname ) - 4] = '\0';
  sprintf( mvst_buff, "%s_buff.mvby", seq_name ); //info.mvstatname is the file!

  if( info.denoise_flag == YES ) {  
    info.ywidth  /= 2;
    info.yheight /= 2;
    info.cwidth  /= 2;
    info.cheight /= 2;
  }
  
  init_enc( &info );  // in init_encN.c     

  initialize_spiral_search(get_searchrange(1 << (info.tPyrLev-1), info));

  write_header( info.bitname, info );
  header_bytes = sizeof( videoheader );

  num_of_GOP = get_GOP_num( info );

  total_bytes_past += sizeof( videoheader );
  total_bytes_past += num_of_GOP * sizeof( long int );

  for(i = 0; i < THEO_COUNT; i ++){
	  lr_record[i].total_bytes_past += sizeof( videoheader );
	  lr_record[i].total_bytes_past += num_of_GOP * sizeof( long int );
  }

  read_GOP_bytes =
    ( long int * )getarray( num_of_GOP, sizeof( long int ),
                            "read_GOP_bytes" );

  //Added on 11.12.2017
  gop_mv = ( unsigned long * )getarray( num_of_GOP, sizeof( unsigned long ), "gop_mv" ); 

  output_GOP_bytes =
      ( long int * )getarray( num_of_GOP, sizeof( long int ),
			      "output_GOP_bytes" );
  for( i = 0; i < num_of_GOP; i++ )
      output_GOP_bytes[i] = 0;

  dec_output_GOP_bytes =
      ( long int * )getarray( num_of_GOP, sizeof( long int ),
			      "dec_output_GOP_bytes" );
  for( i = 0; i < num_of_GOP; i++ )
      dec_output_GOP_bytes[i] = 0;


  if( !( fpio = fopen( info.bitname, "a+b" ) ) ) {
      printf( "can not open: %s\n", info.bitname );
      exit( 1 );
  }

//  printf("Enc num_of_GOP = %d\n",num_of_GOP);

  // leave some space for keeping the GOP byte size 
  fwrite( output_GOP_bytes, sizeof( long int ), num_of_GOP, fpio );
  fclose( fpio );
  header_bytes += num_of_GOP * sizeof( long int );
/*
  orig_lambda[0] = info.lambda[0];
  orig_lambda[1] = info.lambda[1];
  orig_lambda[2] = info.lambda[2];

  best_lambda[0] = orig_lambda[0];
  best_lambda[1] = orig_lambda[1];
  best_lambda[2] = orig_lambda[2];

  curr_lambda[0] = orig_lambda[0];
  curr_lambda[1] = orig_lambda[1];
  curr_lambda[2] = orig_lambda[2];

  buff_best_lambda[0] = orig_lambda[0];
  buff_best_lambda[1] = orig_lambda[1];
  buff_best_lambda[2] = orig_lambda[2];
*/
  for(i = 0; i < ADAPT_LEVEL; i ++){
	orig_lambda[i] = info.lambda[i];
	best_lambda[i] = orig_lambda[i];
	curr_lambda[i] = orig_lambda[i];
	buff_best_lambda[i] = orig_lambda[i];
  }

  first_GOP = YES;
  curr = info.start;
  last = info.last;
  GOP_counter = 0;

//Added on 03.26.2018
  scale = 1 << info.subpel[0];
  assert(info.subpel[0] > 0 && info.subpel[0] == info.subpel[1]);
  scale = scale << ADD_SUB;
  uphor = ( info.ywidth - 1 ) * scale + 1;
  upver = ( info.yheight - 1 ) * scale + 1;

  upsamp_x = ( float * )getarray( uphor * info.yheight, sizeof( float ), "upsamp_x" );
  upframe1 = ( float * )getarray( uphor * upver, sizeof( float ), "upframe1" );
  upframe2 = ( float * )getarray( uphor * upver, sizeof( float ), "upframe2" );

#ifdef VBR_DEC
  vbr_total_bytes_past = total_bytes_past;
  printf("\nVBR_DEC, copy file to buffer!\n");
  file_copy(enc_buff, info.bitname);
  file_copy(rd_buff, rd_name);
  file_copy(mvst_buff, info.mvstatname);
  save_enc_status( &info );

  info.lambda_adapt_thres = 100000;
  SIMUL_POINT = 0;
#endif

  while( curr <= last ) {  //  逐GOP

//	fscanf(fplmbd,"%f %f %f\n",&info.lambda[0],&info.lambda[1],&info.lambda[2]);
        
    remaining_frs = last - curr + 1;       // 剩余帧 
	// 下面检测剩余帧是否还够一个gop  正常情况下，eff_GOPsz等于bigGOP，但是在最后不够一个GOP时，等于剩余帧
    if ( remaining_frs < info.bigGOP ){
      Level_change = YES;             // last bigGOP detected 
      info.eff_GOPsz = remaining_frs; // effective GOP size
      curr_disp = curr + info.GOPsz - 1;
      if (curr_disp > last) curr_disp = last;
    }else{
      Level_change = NO;              // full GOP 
      info.eff_GOPsz = info.bigGOP;
      curr_disp = curr + info.GOPsz - 1;
    }

	gop_block = 0;	frame_cnt = 0;

#ifdef INTRA_SUPPORT    
    switch ( info.intra ) {
      
    case YES:
      
      if( info.denoise_flag == YES ) {
        printf( "can not handle this case\n" );
        exit( 1 );
      }
      
#ifdef THREE_D
      output_GOP_bytes[GOP_counter] = three_D_anal( curr, GOP_counter, info );
#else
      output_GOP_bytes[GOP_counter] = intra_encode( curr, GOP_counter, info );
#endif
     
      break;
      
    case NO:
#endif

//Added new
#ifdef TEST_LAMBDA
	if( (info.GOPsz * GOP_counter) == SIMUL_START && simul_count < SIMUL_LIMIT ){ //simul cycle
		
		GOP_counter = GOP_counter - ( (SIMUL_START - SIMUL_POINT)/info.GOPsz );
		curr = info.start + SIMUL_POINT;
		if( curr == info.start )first_GOP = YES;
		curr_disp = curr + info.GOPsz - 1;
		remaining_frs = last - curr + 1; 

reselect_sign: ;
		if( simul_count < (SIMUL_LIMIT-1) ){ //Simul enc

			for(i = 2;i >= 0;i --){
				if(i != SIMUL_LEVEL){
					if(i <= 1){
						info.lambda[i] = best_lambda[i];
						assert(best_lambda[i] > 0);
					}else{
						assert(i == 2);
						for(j = ADAPT_LEVEL-1;j >= i;j --){
							info.lambda[j] = best_lambda[j];
							assert(best_lambda[j] > 0);
						}
					}
				}
			}

			if( simul_count == 0 ){
				if(round == 0)
					SIMUL_ADD = 1;
				else
					SIMUL_ADD = 1;

				info.lambda[SIMUL_LEVEL] = curr_lambda[SIMUL_LEVEL] - SIMUL_RANGE * SIMUL_ADD;
				if(SIMUL_LEVEL == 2){
					for(j = ADAPT_LEVEL-1;j > SIMUL_LEVEL;j --){
						info.lambda[j] = curr_lambda[j] - SIMUL_RANGE * SIMUL_ADD;
					}
				}

				if(info.lambda[SIMUL_LEVEL] < LAMBDA_LOWER_LMT){
					info.lambda[SIMUL_LEVEL] = LAMBDA_LOWER_LMT;
					if(SIMUL_LEVEL == 2){
						for(j = ADAPT_LEVEL-1;j > SIMUL_LEVEL;j --){
							info.lambda[j] = LAMBDA_LOWER_LMT;
						}
					}
				}
				if(SIMUL_LEVEL >= 1 && info.lambda[SIMUL_LEVEL] < info.lambda[SIMUL_LEVEL-1]){
#ifdef LAMBDA_CTRL
					info.lambda[SIMUL_LEVEL] = info.lambda[SIMUL_LEVEL-1];
#endif
				}

			}else{
				info.lambda[SIMUL_LEVEL] = curr_lambda[SIMUL_LEVEL] + SIMUL_ADD;
				if(SIMUL_LEVEL == 2){
					for(j = ADAPT_LEVEL-1;j > SIMUL_LEVEL;j --){
						info.lambda[j] = curr_lambda[j] + SIMUL_ADD;
					}
				}

				if( info.lambda[SIMUL_LEVEL] == buff_best_lambda[SIMUL_LEVEL] ){
					info.lambda[SIMUL_LEVEL] = info.lambda[SIMUL_LEVEL] + SIMUL_ADD;
					if(SIMUL_LEVEL == 2){
						for(j = ADAPT_LEVEL-1;j > SIMUL_LEVEL;j --){
							info.lambda[j] = info.lambda[j] + SIMUL_ADD;
						}
					}
				}

				for( i = SIMUL_LEVEL+1; i <= 2; i ++ ){
#ifdef LAMBDA_CTRL
					if(info.lambda[i] < info.lambda[i-1])
						info.lambda[i] = info.lambda[i-1];
#endif
				}
			}

			printf("\nlambda revised! lambda0 = %f, lambda1 = %f, lambda2 = %f, lambda3 = %f\n", info.lambda[0], info.lambda[1], info.lambda[2], info.lambda[3]);

			for(i = ADAPT_LEVEL-1;i >= 0;i --)
				curr_lambda[i] = info.lambda[i];

//Added on 05.04.2019
			if(lmda_buff != NULL){
				getnum = compare_lambda_vector(lmda_buff,info.lambda);
				if( getnum == 0 ){
					printf("\nLambda vector not tried, try it.\n");
				}else{
					assert(getnum == 1);
					printf("\nLambda vector already tried, skip it.\n");
					simul_count ++;
					goto reselect_sign;
				}
			}else{
				assert(lmda_buff == NULL);
				segment_init = NO;
				printf("\nLambda vector not tried, try it.\n");
				lmda_buff = new lambda_buff;
				lmda_buff->next = NULL;
				for(i = 0;i < ADAPT_LEVEL;i ++)
					lmda_buff->buff_lambda[i] = info.lambda[i];
			}
//Added on 05.04.2019
		}else{ //Encode with the best choice, when simul_count == (SIMUL_LIMIT-1)
			assert(simul_count == (SIMUL_LIMIT-1));
			for(i = 0; i <= 2; i ++){
				info.lambda[i] = best_lambda[i];
				assert(best_lambda[i] > 0);
				if(i == 2){
					for(j = ADAPT_LEVEL-1;j > i;j --){
						info.lambda[j] = best_lambda[j];
					}
				}
			}
			printf("\nfinal best lambda0 = %f, lambda1 = %f, lambda2 = %f,  lambda3 = %f\n", info.lambda[0], info.lambda[1], info.lambda[2], info.lambda[3]);
		}

		simul_count ++;

		if(simul_count == SIMUL_LIMIT && SIMUL_START < (info.last - 0 * info.GOPsz) ){
			if(SIMUL_LEVEL == 2){
				if( (round == SIMUL_ROUND) || 
					( (buff_best_lambda[0] == best_lambda[0]) &&
					(buff_best_lambda[1] == best_lambda[1]) &&
					(buff_best_lambda[2] == best_lambda[2]) ) ){//Has done three rounds

					printf("Overall end of current GOP simulation reached!\n");
					destroy_lambda_buff(lmda_buff);
					lmda_buff = NULL;

					if( (buff_best_lambda[0] == best_lambda[0]) &&
					(buff_best_lambda[1] == best_lambda[1]) &&
					(buff_best_lambda[2] == best_lambda[2]) )
						printf("Obtained same lambda group.\n");

					round = 0;
					simul_count = 0;

					SIMUL_START = BUFF_SIMUL_START;
					SIMUL_POINT = BUFF_SIMUL_POINT;
					in_simul = NO;//move forward to search for the next segment with similar feature;

					SIMUL_LEVEL = 0;

					gop_end = YES;
				}else{
					printf("Current round end reached!\n");
#ifdef THEO_MODEL
					if(round == 0){
						printf("SIMUL LEVEL 1:\n");
						for(i = 0; i < THEO_COUNT; i ++){
							for(j = 0; j < INIT_RANGE; j ++){
								printf("theo point = %d, lambda point = %d, YPSNR = %f\n",i,j,theo_simul_psnr1[i][j]);
							}
							printf("\n\n");
						}
						printf("SIMUL LEVEL 2:\n");
						for(i = 0; i < THEO_COUNT; i ++){
							for(j = 0; j < INIT_RANGE; j ++){
								printf("theo point = %d, lambda point = %d, YPSNR = %f\n",i,j,theo_simul_psnr2[i][j]);
							}
							printf("\n\n");
						}
						printf("SIMUL LEVEL 3:\n");
						for(i = 0; i < THEO_COUNT; i ++){
							for(j = 0; j < INIT_RANGE+2; j ++){
								printf("theo point = %d, lambda point = %d, YPSNR = %f\n",i,j,theo_simul_psnr3[i][j]);
							}
							printf("\n\n");
						}
					}
#endif
					round ++;
					simul_count = 0;
					SIMUL_LEVEL = 0;
//					best_simul_psnr = 0;
					for(i = 0;i < ADAPT_LEVEL;i ++){
						curr_lambda[i] = best_lambda[i];
						buff_best_lambda[i] = best_lambda[i];
					}
				}

				SIMUL_LIMIT = INIT_RANGE;
				SIMUL_RANGE = SIMUL_LIMIT / 2;

			}else{//Do it for higher/lower temporal level
				assert(SIMUL_LEVEL == 0 || SIMUL_LEVEL == 1);
				for(i = 0;i < ADAPT_LEVEL;i ++)
					curr_lambda[i] = best_lambda[i];

				SIMUL_LEVEL ++;
				printf("SIMUL_LEVEL changed! SIMUL_LEVEL = %d\n",SIMUL_LEVEL);
				simul_count = 0;

				if(SIMUL_LEVEL == 1)
					SIMUL_LIMIT += 0;
				else{
					assert(SIMUL_LEVEL == 2);
					SIMUL_LIMIT += 2;
				}

				SIMUL_RANGE = SIMUL_LIMIT / 2;
			}
		}

		if ( remaining_frs < info.bigGOP ){
		  Level_change = YES;             // last bigGOP detected 
		  info.eff_GOPsz = remaining_frs; // effective GOP size
		  curr_disp = curr + info.GOPsz - 1;
		  if (curr_disp > last) curr_disp = last;
		}else{
		  Level_change = NO;              // full GOP 
		  info.eff_GOPsz = info.bigGOP;
		  curr_disp = curr + info.GOPsz - 1;
		}

//copy related files back from buffer
		file_copy(info.bitname, enc_buff);
		file_copy(rd_name, rd_buff);
		file_copy(info.mvstatname, mvst_buff);

		printf( "******************************************************\n");
		printf( "*  frame %.3d ~ frame %.3d (%.3d) (remaining_frs: %.3d)  *\n", 
				curr, curr_disp, curr + info.eff_GOPsz-1, remaining_frs );
		printf( "******************************************************\n");

		resume_enc_status( &info );

		simul_skip = YES;

		output_GOP_bytes[GOP_counter] =
			  mctf_anal_ezbc( curr, GOP_counter, info, first_GOP,
			   Level_change, remaining_frs, &sum_mv, YES, upframe1, upframe2, upsamp_x);

		simul_skip = NO;

		if( gop_end == YES ){
			for(i = 0;i < ADAPT_LEVEL;i ++){
				orig_lambda[i] = best_lambda[i];
			}

			for(i = 0;i <= 2;i ++)
				best_simul_psnr[i] = 0;

			for(i = 0;i < ADAPT_LEVEL;i ++){
				curr_lambda[i] = orig_lambda[i];
				best_lambda[i] = orig_lambda[i];
				buff_best_lambda[i] = orig_lambda[i];
			}

			gop_end = NO;
		}

		simul_alert = YES;

		total_bytes_past = total_bytes_past_buffer;
		for(i = 0; i < THEO_COUNT; i ++){
			 lr_record[i].total_bytes_past = lr_record[i].total_bytes_past_buffer;
		}
	}else{ //if non-simul
#endif
//Added new
		printf( "******************************************************\n");
		printf( "*  frame %.3d ~ frame %.3d (%.3d) (remaining_frs: %.3d)  *\n", 
		        curr, curr_disp, curr + info.eff_GOPsz-1, remaining_frs );  // 当前gop要编码的帧的范围：curr起始，curr_disp结束
		printf( "******************************************************\n");

#ifdef TEST_LAMBDA
		if( (info.GOPsz * GOP_counter) == info.start + SIMUL_POINT){
//copy related files into buffer
#ifndef VBR_DEC
			printf("\ncopy file to buffer!\n");
			file_copy(enc_buff, info.bitname);
			file_copy(rd_buff, rd_name);
			file_copy(mvst_buff, info.mvstatname);
			save_enc_status( &info );
#else
			printf("\nVBR_DEC, file is ready in the buffer!\n");
#endif

			total_bytes_past_buffer = total_bytes_past;
			for(i = 0; i < THEO_COUNT; i ++){
			  lr_record[i].total_bytes_past_buffer = lr_record[i].total_bytes_past;
			}

			printf("\nlambda revised! lambda0 = %f, lambda1 = %f, lambda2 = %f, lambda3 = %f\n", info.lambda[0], info.lambda[1], info.lambda[2], info.lambda[3]);
		}
		gop_block = 0;	frame_cnt = 0;
		printf("simul_alert = %d\n",simul_alert);

//		if(curr == (SIMUL_START - info.GOPsz) )
//			simul_skip = YES;
#endif
		if( info.denoise_flag == YES ) {
        
			output_GOP_bytes[GOP_counter] =
				denoise_mctf_anal_ezbc( curr, GOP_counter, info, first_GOP,
									  Level_change, remaining_frs );
		} else {
			output_GOP_bytes[GOP_counter] =
				// 起始帧，第几个gop，视频信息，是否为第一个gop，是不是不足一个gop，剩余帧，mv和，NO，
			  mctf_anal_ezbc( curr, GOP_counter, info, first_GOP,
			   Level_change, remaining_frs, &sum_mv, NO, upframe1, upframe2, upsamp_x); // 真正编码GOP的入口
		}
#ifdef TEST_LAMBDA
		gop_block = gop_block / frame_cnt;
		printf("Current avg frame block number: %f\n",gop_block);

		if(first_GOP == YES){
			buff_gop_block = gop_block;
/*			SIMUL_POINT = info.GOPsz;
			SIMUL_START = 4 * info.GOPsz;
*/
			in_simul = NO;
			printf("Initial state! New buff_gop_block = %f\n",gop_block);
			printf("buff_gop_block = %f\n", buff_gop_block);
		}else{
			assert(first_GOP == NO);
			if( curr == (info.start + SIMUL_POINT) && in_simul == NO ){
				buff_gop_block = gop_block;
				printf("New buff_gop_block start! Buff_gop_block = %f\n",buff_gop_block);
			}

			if( curr == (info.start + SIMUL_START - 3 * info.GOPsz) && in_simul == NO ){
#ifndef VBR_DEC
				buff_gop_block = gop_block;
				printf("Updated buff_gop_block for new segment = %f\n", buff_gop_block);
#endif
			}

			if( curr == (info.start + SIMUL_START - 2 * info.GOPsz) && in_simul == NO ){
				if( fabs(gop_block - buff_gop_block) > info.lambda_adapt_thres || curr >= (info.last + 1 - 3 * info.GOPsz) ){
					printf("gop_block = %f, buff_gop_block = %f\n", gop_block, buff_gop_block);
					buff_gop_block = gop_block;
					BUFF_SIMUL_POINT = curr;
					BUFF_SIMUL_START = curr + 3 * info.GOPsz;
					in_simul = YES;
					printf("Feature change detected! SIMUL_START = %d\n",SIMUL_START);
					printf("Updated buff_gop_block = %f\n", buff_gop_block);
				}else{
					SIMUL_START += info.GOPsz;
					in_simul = NO;
					printf("Feature continued! SIMUL_POINT = %d, increased SIMUL_START = %d\n",SIMUL_POINT, SIMUL_START);
					printf("gop_block = %f, buff_gop_block = %f\n", gop_block, buff_gop_block);
				}
			}//when curr == (SIMUL_START - info.GOPsz)
		}//not first GOP

	}//if non-simul
#endif
#ifdef INTRA_SUPPORT
      break;
      
    default:
      
      printf( "error (encoderN.c)\n" );
      exit( 1 );
    
    }
#endif
	gop_mv[GOP_counter] = sum_mv;

#ifdef TEST_LAMBDA

//////////////////////////////////////////////////
#ifdef THEO_MODEL
	BOTTOM_BITRATE = info.simul_rate - 200;

	if( curr < (info.last - info.GOPsz) && round == 0 ){
		info.bitrate = BOTTOM_BITRATE;
		for(i = 0; i < THEO_COUNT; i ++){

			FAT_bytes = num_of_GOP * sizeof( long int );
			GOPheader_bytes = 2 * num_of_GOP;

			sprintf( post_stream, "%s_%d.bit", seq_name, info.bitrate );
			estimated_overhead = sizeof( videoheader ) + GOPheader_bytes + FAT_bytes;

#ifdef VBR_DEC
			test( info, &info.qp, estimated_overhead, gop_mv, curr + info.GOPsz - 1, 1 );
//simul decoding
//info.bitname 
			strcpy(raw_buffer,info.bitname);
			strcpy(info.bitname,post_stream);
/////////////////////////////////////////

			if( curr == (info.start + SIMUL_START - info.GOPsz) ){
				lr_record[i].total_bytes_past = vbr_total_bytes_past;
				VBR_GOP_counter = 0;
				VBR_first_GOP = YES;
				vbr_curr = info.start;
				VBR_Level_change = NO;

				while( VBR_GOP_counter <= GOP_counter ){
					vbr_remaining_frs = info.last - vbr_curr + 1;

					if( !( fpio = fopen( info.bitname, "rb" ) ) ) {
						printf( "can not open: %s\n", info.bitname );
						exit( 1 );
					}
					fseek( fpio, sizeof( videoheader ) + VBR_GOP_counter * sizeof( long int ), SEEK_SET );

					fread( &dec_output_GOP_bytes[VBR_GOP_counter], sizeof( long int ), 1, fpio );

					printf("dec_output_GOP_bytes[VBR_GOP_counter] = %d\n",dec_output_GOP_bytes[VBR_GOP_counter]);
					printf("Current total_bytes_past = %d\n",lr_record[i].total_bytes_past);

					fclose( fpio );

					info.GOPbytes = dec_output_GOP_bytes[VBR_GOP_counter];

					if( simul_alert == YES ){
						mctf_syn_ezbc( vbr_curr, VBR_GOP_counter, &lr_record[i].total_bytes_past, info,
										 VBR_first_GOP, Level_change, vbr_remaining_frs, YES, NO );

						simul_alert = NO;
					}else{
						mctf_syn_ezbc( vbr_curr, VBR_GOP_counter, &lr_record[i].total_bytes_past, info,
										 VBR_first_GOP, Level_change, vbr_remaining_frs, NO, NO );
					}

					lr_record[i].total_bytes_past += info.GOPbytes;
					VBR_GOP_counter ++;
					VBR_first_GOP = NO;
					vbr_curr += info.GOPsz;
				}//end while
			}//end if curr

			if( curr == (info.start + SIMUL_START - info.GOPsz) ){
				printf("\nCalculate PSNR info! Obtained by lambda0 = %f, lambda1 = %f, lambda2 = %f\n", curr_lambda[0],curr_lambda[1],curr_lambda[2]);

				if(curr - info.GOPsz < info.start)
					assert(0);
				else{
					if(SIMUL_LEVEL == 0){
						theo_simul_psnr1[i][simul_count] = calsnr(SIMUL_POINT, curr-1, info);//frame 12-24
						printf("Theo simul_count = %d, THEO_COUNT = %d, YPSNR = %f\n\n",simul_count, i, theo_simul_psnr1[i][simul_count]);
					}else if(SIMUL_LEVEL == 1){
						theo_simul_psnr2[i][simul_count] = calsnr(SIMUL_POINT, curr-1, info);//frame 8-22
						printf("Theo simul_count = %d, THEO_COUNT = %d, YPSNR = %f\n\n",simul_count, i, theo_simul_psnr2[i][simul_count]);
					}else{
						assert(SIMUL_LEVEL == 2);
						theo_simul_psnr3[i][simul_count] = calsnr(SIMUL_POINT, curr-1, info);//frame 4-20
						printf("Theo simul_count = %d, THEO_COUNT = %d, YPSNR = %f\n\n",simul_count, i, theo_simul_psnr3[i][simul_count]);
					}
				}
			}

//simul decoding ends for certain bitrate
			strcpy(info.decname, dec_buff);
			info.bitrate += 200;
			strcpy(info.bitname, raw_buffer);
#else
			test( info, &info.qp, estimated_overhead, gop_mv, curr + info.GOPsz - 1, 0 );
//simul decoding
//info.bitname 
			strcpy(raw_buffer,info.bitname);
			strcpy(info.bitname,post_stream);

			if( !( fpio = fopen( info.bitname, "rb" ) ) ) {
				printf( "can not open: %s\n", info.bitname );
				exit( 1 );
			}
			fseek( fpio, sizeof( videoheader ) + GOP_counter * sizeof( long int ), SEEK_SET );

			fread( &dec_output_GOP_bytes[GOP_counter], sizeof( long int ), 1, fpio );

			printf("dec_output_GOP_bytes[GOP_counter] = %d\n",dec_output_GOP_bytes[GOP_counter]);
			printf("Current theo total_bytes_past = %d\n",lr_record[i].total_bytes_past);

			fclose( fpio );

//decoding start
			strcpy(dec_buff, info.decname);
			sprintf( info.decname, "%s_%d.yuv", seq_name, info.bitrate );

			info.GOPbytes = dec_output_GOP_bytes[GOP_counter];

			mctf_syn_ezbc( curr, GOP_counter, &(lr_record[i].total_bytes_past), info,
				first_GOP, Level_change, remaining_frs, simul_alert, YES );

			lr_record[i].total_bytes_past += info.GOPbytes;

			if( curr == (info.start + SIMUL_START - info.GOPsz) ){
				printf("\nCalculate PSNR info! Obtained by lambda0 = %f, lambda1 = %f, lambda2 = %f\n", curr_lambda[0],curr_lambda[1],curr_lambda[2]);

				if(curr - info.GOPsz < info.start)
					assert(0);
				else{
					if(SIMUL_LEVEL == 0)
						theo_simul_psnr[i] = calsnr(SIMUL_POINT, curr, info);//frame 12-24
					else if(SIMUL_LEVEL == 1)
						theo_simul_psnr[i] = calsnr(SIMUL_POINT, curr, info);//frame 8-22
					else{
						assert(SIMUL_LEVEL == 2);
						theo_simul_psnr[i] = calsnr(SIMUL_POINT, curr, info);//frame 4-20
					}
				}

			}
//simul decoding ends for certain bitrate
			strcpy(info.decname, dec_buff);
			info.bitrate += 100;
			strcpy(info.bitname, raw_buffer);
#endif
		}//THEO_COUNT
/*
		printf("PSNR values for theoretical model:\n");
		getnum = BOTTOM_BITRATE;
		for(i = 0; i < THEO_COUNT; i ++){
			printf("Bitrate: %d, PSNR: %f\n", getnum,theo_simul_psnr[i]);
			getnum += 100;
		}
*/
	}//end if
#endif
//////////////////////////////////////////////////
	
	info.bitrate = info.simul_rate;
	FAT_bytes = num_of_GOP * sizeof( long int );
	GOPheader_bytes = 2 * num_of_GOP;

	sprintf( post_stream, "%s_%d.bit", seq_name, info.bitrate );

	estimated_overhead = sizeof( videoheader ) + GOPheader_bytes + FAT_bytes;

#ifdef VBR_DEC
	test( info, &info.qp, estimated_overhead, gop_mv, curr + info.GOPsz - 1, 1 );
#else
	test( info, &info.qp, estimated_overhead, gop_mv, curr + info.GOPsz - 1, 0 );
#endif
//simul decoding
//	info.bitname 
	strcpy(raw_buffer,info.bitname);
	strcpy(info.bitname,post_stream);

#ifdef VBR_DEC
	if( curr == (info.start + SIMUL_START - info.GOPsz) ){
		total_bytes_past = vbr_total_bytes_past;
		VBR_GOP_counter = 0;
		VBR_first_GOP = YES;
		vbr_curr = info.start;
		VBR_Level_change = NO;

		while( VBR_GOP_counter <= GOP_counter ){
			vbr_remaining_frs = info.last - vbr_curr + 1;

			if( !( fpio = fopen( info.bitname, "rb" ) ) ) {
				printf( "can not open: %s\n", info.bitname );
				exit( 1 );
			}
			fseek( fpio, sizeof( videoheader ) + VBR_GOP_counter * sizeof( long int ), SEEK_SET );

			fread( &dec_output_GOP_bytes[VBR_GOP_counter], sizeof( long int ), 1, fpio );

			printf("dec_output_GOP_bytes[VBR_GOP_counter] = %d\n",dec_output_GOP_bytes[VBR_GOP_counter]);
			printf("Current total_bytes_past = %d\n",total_bytes_past);

			fclose( fpio );

			info.GOPbytes = dec_output_GOP_bytes[VBR_GOP_counter];

			if( simul_alert == YES ){
				mctf_syn_ezbc( vbr_curr, VBR_GOP_counter, &total_bytes_past, info,
								 VBR_first_GOP, Level_change, vbr_remaining_frs, YES, NO );

				simul_alert = NO;
			}else{
				mctf_syn_ezbc( vbr_curr, VBR_GOP_counter, &total_bytes_past, info,
								 VBR_first_GOP, Level_change, vbr_remaining_frs, NO, NO );
			}

			total_bytes_past += info.GOPbytes;
			VBR_GOP_counter ++;
			VBR_first_GOP = NO;
			vbr_curr += info.GOPsz;
		}//end while
	}//end if
#else
	if( !( fpio = fopen( info.bitname, "rb" ) ) ) {
		printf( "can not open: %s\n", info.bitname );
		exit( 1 );
	}
	fseek( fpio, sizeof( videoheader ) + GOP_counter * sizeof( long int ), SEEK_SET );

	fread( &dec_output_GOP_bytes[GOP_counter], sizeof( long int ), 1, fpio );

	printf("dec_output_GOP_bytes[GOP_counter] = %d\n",dec_output_GOP_bytes[GOP_counter]);
	printf("Current total_bytes_past = %d\n",total_bytes_past);

    fclose( fpio );

	info.GOPbytes = dec_output_GOP_bytes[GOP_counter];

	if( simul_alert == YES ){
		mctf_syn_ezbc( curr, GOP_counter, &total_bytes_past, info,
						 first_GOP, Level_change, remaining_frs, YES, NO );

		simul_alert = NO;
	}else{
		mctf_syn_ezbc( curr, GOP_counter, &total_bytes_past, info,
						 first_GOP, Level_change, remaining_frs, NO, NO );
	}

	total_bytes_past += info.GOPbytes;
#endif
	strcpy(info.bitname,raw_buffer);

	if( curr == (info.start + SIMUL_START - info.GOPsz) ){
		printf("\nCalculate PSNR info! Obtained by lambda0 = %f, lambda1 = %f, lambda2 = %f, lambda3 = %f\n", curr_lambda[0],curr_lambda[1],curr_lambda[2],curr_lambda[3]);
#ifndef VBR_DEC
			if(curr - info.GOPsz < info.start)
				assert(0);
			else{
				if(SIMUL_LEVEL == 0)
					simul_psnr = calsnr(SIMUL_POINT, curr-1, info);//frame 12-24
				else if(SIMUL_LEVEL == 1)
					simul_psnr = calsnr(SIMUL_POINT, curr-1, info);//frame 8-22
				else{
					assert(SIMUL_LEVEL == 2);
					simul_psnr = calsnr(SIMUL_POINT, curr-1, info);//frame 4-20
				}
			}
#else
			if(curr - info.GOPsz < info.start)
				assert(0);
			else{
				if(SIMUL_LEVEL == 0)
					simul_psnr = calsnr(SIMUL_POINT, info.act_last, info);
				else if(SIMUL_LEVEL == 1)
					simul_psnr = calsnr(SIMUL_POINT, info.act_last, info);
				else{
					assert(SIMUL_LEVEL == 2);
					simul_psnr = calsnr(SIMUL_POINT, info.act_last, info);
				}
			}
#endif

		if( simul_psnr > (best_simul_psnr[SIMUL_LEVEL]) ){
			best_simul_psnr[SIMUL_LEVEL] = simul_psnr;
			for(i = 0;i < ADAPT_LEVEL;i ++)
				best_lambda[i] = curr_lambda[i];
		}
		printf("\nBest PSNR current level = %f, lambda0 = %f, lambda1 = %f, lambda2 = %f, lambda3 = %f\n",best_simul_psnr[SIMUL_LEVEL], best_lambda[0], best_lambda[1], best_lambda[2], best_lambda[3]);
		for(i = 0;i <= 2;i ++)
			printf("TEMP level %d, best_simul_psnr = %f\n",i,best_simul_psnr[i]);
	}
#endif

	first_GOP = NO; 
    GOP_counter++;
    curr += info.GOPsz;

  }                             /* while */

#ifdef TEST_LAMBDA

  // decode LAST GOP
  {
    info.eff_GOPsz = remaining_frs;
    remaining_frs = 0; // no bitstream extraction of last GOP, only temporal filtering
    
    printf( " decoding frame %d - %d (last_GOP) .....\n", curr - info.GOPsz + 1, last );

//	info.bitname 
	strcpy(raw_buffer,info.bitname);
	strcpy(info.bitname,post_stream);

	if( !( fpio = fopen( info.bitname, "rb" ) ) ) {
		printf( "can not open: %s\n", info.bitname );
		exit( 1 );
	}
	fseek( fpio, sizeof( videoheader ) + GOP_counter * sizeof( long int ), SEEK_SET );

	fread( &dec_output_GOP_bytes[GOP_counter], sizeof( long int ), 1, fpio );

	printf("dec_output_GOP_bytes[GOP_counter] = %d\n",dec_output_GOP_bytes[GOP_counter]);

    fclose( fpio );

	info.GOPbytes = dec_output_GOP_bytes[GOP_counter];

    mctf_syn_ezbc( curr, GOP_counter, &total_bytes_past, info,
                   first_GOP, Level_change, remaining_frs, NO, NO );
    
//    total_bytes_past += info.GOPbytes;

	strcpy(info.bitname,raw_buffer);
  }

#endif

  /* pstatN.c *//* cfr and pfr are just used as the memory in calsnr, and the data in cfr and pfr is not used in calsnr */

  if( !( fpio = fopen( info.bitname, "r+b" ) ) ) {
    printf( "can not open: %s\n", info.bitname );
    exit( 1 );
  }
  fseek( fpio, sizeof( videoheader ), SEEK_SET );
  fwrite( output_GOP_bytes, sizeof( long int ), num_of_GOP, fpio );
  fclose( fpio );

//  fclose(fplmbd);

  clean_up_spiral_search(get_searchrange(1 << (info.tPyrLev-1), info));

  elp = clock(  ) - mark;
  duration = ( double )elp / CLOCKS_PER_SEC;
  print_time( duration );

  mvStat_close();

  free(gop_mv); //Added on 11.12.2017
  free(dec_output_GOP_bytes);
  free(output_GOP_bytes);

  free(upsamp_x);
  free(upframe1);
  free(upframe2);

//  enc_destructor( info );
  free(buff_frameMEinfo);

  printf( "finished.\n" );
  return 0;
}


/*
 *              usage()                                    
 */
void
usage(  )
{
  printf( "3DSBCen parameter_file \n" );
}

/*
 *                                error_check()                              
 */
void
error_check( videoinfo info )
{
  int i;

  if( !strcmp( info.inname, "NULL" ) ) {
    printf( "error in read_command(): specify -inname input\n" );
    exit( 1 );
  }
  if( !strcmp( info.bitname, "NULL" ) ) {
    printf( "error in read_command(): specify -bitname bitfile\n" );
    exit( 1 );
  }
  if( !strcmp( info.statname, "NULL" ) ) {
    printf( "error in read_command(): specify -statname statusfile\n" );
    exit( 1 );
  }
  for (i = 0; i < MAX_TLEVELS; i++) {
    if( info.level[i] == -1 ) {
      printf( "error in read_command(): specify -MBlevel\n" );
      exit( 1 );
    }
  }
  if( info.start < 0 || info.start > info.last ) {
    printf( "error in read_command(): start %d last %d\n", info.start,
            info.last );
    printf( "0<= start, last <1000 and start<=last\n" );
    exit( 1 );
  }
  if( info.ywidth < 0 || info.yheight < 0 || info.cwidth < 0
      || info.cheight < 0 ) {
    printf( "error in read_command(): image size\n" );
    printf( " Y: %d x %d, C: %d x %d\n", info.ywidth, info.yheight,
            info.cwidth, info.cheight );
    exit( 1 );
  }
  if( info.start < 0 ) {
    printf( "can not handle this case (encoderN.c)\n" );
    exit( 1 );
  }

  if( info.ywidth % 8 ) {
    printf( "error: Y width %d should be multiple of 8.(encoderN.c)\n",
            info.ywidth );
    exit( 1 );
  }
  if( info.yheight % 8 ) {
    printf( "error: Y height %d should be multiple of 8.(encoderN.c)\n",
            info.yheight );
    exit( 1 );
  }
  if( info.cwidth % 8 ) {
    printf( "C width %d should be multiple of 8.\n", info.cwidth );     //exit(1);
  }
  if( info.cheight % 8 ) {
    printf( "C height %d should be multiple of 8.\n", info.cheight );   // exit(1);
  }
  if( info.bitrate < 0 ) {
    printf( "error in read_command(): rate %d\n", info.bitrate );
    exit( 1 );
  }
#ifdef INTRA_SUPPORT
  if( info.intra == YES && info.denoise_flag == YES ) {
    printf
      ( "error in read_command(): both info.intra == YES && info.denoise_flag == YES \n" );
    exit( 1 );
  }
#endif
  for (i = 0; i < MAX_TLEVELS; i++) {
    if (info.xblk[i] < (1 << (info.level[i] - 1)) || 
        info.yblk[i] < (1 << (info.level[i] - 1))) {
      printf("error in read_command(): too many block splitting levels (%d) "
             "for macroblock size (%dx%d) in t-level %d\n", 
             info.level[i], info.xblk[i], info.yblk[i], i);
      exit(1);
    }
  }
}

  /* quad tree structured should be cleaned before starting the next */
  /* so, clean the structure except 1st layer */
  /* fmv->child and fmv->lifting_mode should be initialized to 0 */

/*
 *                              read_command()                               
 */
void
read_command( int argc, char **argv, videoinfo * info )
{
  int i, j, k, argnum = 1, tpyr, range = 0;
  char iline[256], token[80], istring[80], strtmp[256];
  FILE *fppar;
  /********** setting the initial or default values ************/
  strcpy( info->inname, "NULL" );
  strcpy( info->bitname, "NULL" );
  strcpy( info->statname, "NULL" );
  strcpy( info->tmpname, "mvbit.tmp" );
  strcpy( info->jp2kname, "NULL" );
  strcpy( info->jp2k_decname, "NULL" );

  for (i = 0; i < MAX_TLEVELS; i++) {// 最大做MAX_TLEVELS次时域滤波
    info->level[i] = -1;
    info->lambda[i] = 24.;
    info->subpel[i] = 2;
	info->AGP_level[i]  = 0; 
	info->bi_mv[i]      = 1;   // indicate the existence of RIGHT_CONNECTED mode
	info->bi_exist[i]   = 1;   // indicate the existence of mv on RIGHT side 
  // bi_mv[]==1 RIGHT_CONNECTED mode exists: 
  // mv must be bi-directional and bi_exist[] must be 1
  // bi_mv[]==0 RIGHT_CONNECTED mode is prevented: 
  // mv on RIGHT side can be discarded and bi_exist[] can be 0 or 1. 
	
	info->layer_mv[i] = 0; // no layer structure for motion vector coding in each temporal level
	info->layer_exist[i] = 0; 
  }
  info->ME = 3;                 /* ???? not output to bitstream */
  info->start = -1;
  info->last = -1;
  info->act_last = -1;
  info->ywidth = 352;
  info->yheight = 240;
  info->cwidth = 176;
  info->cheight = 120;
  info->framerate = 24;
  info->bitrate = 0;
  info->verbose = 0;
  info->tPyrLev = 3;
  info->GOPsz = 0x1 << info->tPyrLev;
  info->bigGOP = ( 0x1 << ( info->tPyrLev + 1 ) ) - 1;// gop size * 2 - 1
  info->eff_GOPsz = info->bigGOP;
  info->GOPbytes = MAX_STD_INT;
  info->t_level = 0;
  info->s_level = 0;
  info->searchrange = -1;
  info->maxsearchrange = info->ywidth;
  info->format = YUV;
  info->pixeldepth = 8;
  info->maxMBnum = -1;

  info->denoise_flag = NO;

  info->SLTF_range = 0;
  info->SHTF_range = 0;

  
  
  /************* read the argument and set the value ********/

  if( argc == 1 ) {
    usage(  );
    exit( 1 );
  }

  for( i = 1; i < argc; i++ ) {  // 打开配置文件，放在fppar上 
    if( *( argv[i] ) == '-' ) {

      switch ( *( ++argv[i] ) ) {
      default:
        printf( "-%c such an option is not available\n", *( argv[i] ) );
        usage(  );
        exit( 1 );
      case 'h':
        usage(  );
        exit( 1 );
        break;
      }
    } else {
      switch ( argnum ) {
      default:
        printf( "more parameters are specified\n" );
        usage(  );
        exit( 1 );
      case 1:
        fppar = fopen( argv[i], "rt" );
        if( fppar == NULL ) {
          printf( "can not open file %s\n", argv[i] );
          usage(  );
          exit( 1 );
        }
        argnum++;
        break;
      }
    }
  }

  while( fgets( iline, 254, fppar ) ) { /* read one line */
    sscanf( iline, "%s", token );

    if( !strcmp( token, "-inname" ) ) {
      sscanf( iline, "%s%s", token, info->inname );
    }else if( !strcmp( token, "-jp2kname" ) ) {
		sscanf( iline, "%s%s", token, info->jp2kname );
    } else if( !strcmp( token, "-jp2k_decname" ) ) {
		sscanf( iline, "%s%s", token, info->jp2k_decname );
    } else if( !strcmp( token, "-decname" ) ) {
      sscanf( iline, "%s%s", token, info->decname );
    } else if( !strcmp( token, "-bitname" ) ) {
      sscanf( iline, "%s%s", token, info->bitname );
      strncpy( strtmp, info->bitname, strlen( info->bitname ) - 4 ); // Hanke, 16.09.02
      strtmp[strlen( info->bitname ) - 4] = '\0';                    //
      sprintf( info->mvstatname, "%s%s", strtmp, ".mvby" );          //
      sprintf( info->tsubname, "%s%s", strtmp, ".t" );               //
    } else if( !strcmp( token, "-statname" ) ) {
      sscanf( iline, "%s%s", token, info->statname );
    } else if( !strcmp( token, "-simulrate" ) ) {
		sscanf( iline, "%s%d", token, &( info->simul_rate) );
    } else if( !strcmp( token, "-lambdathres" ) ) {
		sscanf( iline, "%s%d", token, &( info->lambda_adapt_thres) );
    } else if( !strcmp( token, "-start" ) ) {
      sscanf( iline, "%s%d", token, &( info->start ) );
    } else if( !strcmp( token, "-last" ) ) {
      sscanf( iline, "%s%d", token, &( info->last ) );
    } else if( !strcmp( token, "-act_last" ) ) {
      sscanf( iline, "%s%d", token, &( info->act_last ) );
    }else if( !strcmp( token, "-size" ) ) {
      sscanf( iline, "%s%d%d%d%d", token, &( info->ywidth ),
              &( info->yheight ), &( info->cwidth ), &( info->cheight ) );
	  // save the original resolution for frequency roll-off 
	  info->org_yheight = info->yheight; 
	  info->org_ywidth  = info->ywidth; 

    } else if( !strcmp( token, "-framerate" ) ) {
      sscanf( iline, "%s%d", token, &( info->framerate ) );
    } else if( !strcmp( token, "-searchrange" ) ) {
      sscanf( iline, "%s%d", token, &( info->searchrange ) );
    } else if( !strcmp( token, "-maxsearchrange" ) ) {
      sscanf( iline, "%s%d", token, &( info->maxsearchrange ) );
    }
#ifdef SUPPORT_INTRA
    else if( !strcmp( token, "-intra" ) ) {
      sscanf( iline, "%s%s", token, istring );
      if( !strcmp( istring, "YES" ) ) {
        info->intra = YES;
      } else if( !strcmp( istring, "NO" ) ) {
        info->intra = NO;
      } else {
        printf( "-intra options error\n" );
        exit( 1 );
      }
    } 
#endif
    else if( !strcmp( token, "-denoise" ) ) {
      sscanf( iline, "%s%s", token, istring );
      if( !strcmp( istring, "YES" ) ) {
        info->denoise_flag = YES;
      } else if( !strcmp( istring, "NO" ) ) {
        info->denoise_flag = NO;
      } else {
        printf( "-denoise options error\n" );
        exit( 1 );
      }
    }

    else if( !strcmp( token, "-motion" ) ) {
      sscanf( iline, "%s%s", token, istring );
      if( (!strcmp(istring, "hvsbm")) || 
          (!strcmp(istring, "hvsbm_fullres")) ) {
        for (i = 0; i < MAX_TLEVELS; i++) {
          info->xblk[i] = 64;
          info->yblk[i] = 64;
          info->level[i] = 5;
        }
        info->ME = 3;
      } else {
        printf( "-motion options error (hvsbm)\n" );
        exit( 1 );
      }
    } else if( !strcmp( token, "-lambda" ) ) {
      sscanf(iline, "%s", token);
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%f", &(info->lambda[j]));
		// switch the ' ' back and forth to get the number for lambda  by Yongjun Wu
        for (; iline[k] != ' '; k++); 
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No lambda value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of lambda values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->lambda[k] = info->lambda[k-1];
      }
    } else if( !strcmp( token, "-afflambda" ) ) {
      sscanf(iline, "%s", token);
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%f", &(info->aff_lambda[j]));
		// switch the ' ' back and forth to get the number for lambda  by Yongjun Wu
        for (; iline[k] != ' '; k++); 
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No lambda value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of lambda values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->aff_lambda[k] = info->aff_lambda[k-1];
      }
    }else if( !strcmp( token, "-AGP_level" ) ) { // the sign for AGP of motion vector coding
      sscanf(iline, "%s", token);
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%d", &(info->AGP_level[j]));
		// switch the ' ' back and forth to get the number for lambda  by Yongjun Wu
        for (; iline[k] != ' '; k++); 
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No AGP_level value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of AGP_level values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->AGP_level[k] = info->AGP_level[k-1];
      }
    }  else if( !strcmp( token, "-bi_mv" ) ) {  // the sign for alternative Haar reconstruction in decoder
      sscanf(iline, "%s", token);
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%d", &(info->bi_mv[j]));
		// switch the ' ' back and forth to get the number for lambda  by Yongjun Wu
        for (; iline[k] != ' '; k++); 
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No bi_mv value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of bi_mv values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->bi_mv[k] = info->bi_mv[k-1];
      }
    } else if( !strcmp( token, "-layer_mv" ) ) { // the layerd structure for motion vector coding
      sscanf(iline, "%s", token);
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%d", &(info->layer_mv[j]));
		// switch the ' ' back and forth to get the number for lambda  by Yongjun Wu
        for (; iline[k] != ' '; k++); 
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No layer_mv value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of layer_mv values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->layer_mv[k] = info->layer_mv[k-1];
      }
    }else if( !strcmp( token, "-MVaccuracy" ) ) {
      sscanf(iline, "%s", token);
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%d", &(info->subpel[j]));
		// switch the ' ' back and forth to get the number for MVaccuracy  by Yongjun Wu
        for (; iline[k] != ' '; k++);
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No subpel value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of subpel values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->subpel[k] = info->subpel[k-1];
      }
    } else if( !strcmp( token, "-SLTF_range" ) ) {
      sscanf( iline, "%s%d", token, &range );
      if( range < 0 || range > 10 ) { 
        printf( "illegal SLTF transition range\n" );  
        exit( 1 );
      } else {
        info->SLTF_range = range;
      }
    } else if( !strcmp( token, "-SHTF_range" ) ) {
      sscanf( iline, "%s%d", token, &range );
      if( range < 0 || range > 10 ) { 
        printf( "illegal SHTF transition range\n" );  
        exit( 1 );
      } else {
        info->SHTF_range = range;
      }
    } else if( !strcmp( token, "-tPyrLev" ) ) {
      sscanf( iline, "%s%d", token, &tpyr ); // &( info->tPyrLev ) );
      info->tPyrLev = (short int) tpyr;
      
      if( info->tPyrLev < 1 || info->tPyrLev > MAX_TLEVELS ) { 
        printf( "illegal temporal pyr. layer number(1 - %d)\n",
                MAX_TLEVELS); 
        exit( 1 );
      } else {
        //          nFrsPyr = (int)pow(2.0, (double)info->tPyrLev);
        info->GOPsz = 0x1 << info->tPyrLev;
        info->bigGOP = ( 0x1 << ( info->tPyrLev + 1 ) ) - 1;
      }
    } else if( !strcmp( token, "-MBsize" ) ) {
      sscanf( iline, "%s", token );
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%d", &(info->xblk[j]));
        info->yblk[j] = info->xblk[j];
        for (; iline[k] != ' '; k++);
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No MBsize value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of MBsize values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->xblk[k] = info->xblk[k-1];
        info->yblk[k] = info->xblk[k];
      }
    } else if( !strcmp( token, "-MBlevel" ) ) {
      sscanf( iline, "%s", token );
      j = 0;
      for (k = strlen(token); iline[k] == ' '; k++);// 把空格读掉
      while (j < MAX_TLEVELS && k < (int)strlen(iline) && iline[k] != '#') {
        sscanf(&(iline[k]), "%d", &(info->level[j]));// 赋值到level中，有几个赋几个，剩下的用前面的赋值
        for (; iline[k] != ' '; k++);
        for (; iline[k] == ' '; k++);
        j++;
      }
      if (j == 0) {
        printf("ERROR: No MBlevel value given!\n");
        exit(1);
      } else if (j == MAX_TLEVELS) {
        printf("WARNING: Maxmimum number of MBlevel values reached!\n");
      }
      for (k = j; k < MAX_TLEVELS; k++) {
        info->level[k] = info->level[k-1];
      }
    } else if( token[0] == '-' ) {
      printf( "%s such a token is not available\n", token );
      usage(  );
      exit( 1 );
    } else {
    }
  }
  fclose( fppar );

  if(info->act_last == -1){
	  assert(info->last >= 0);
	  info->act_last = info->last;
  }

  // how many sub-symbols exist
  for (k = 0; k < MAX_TLEVELS; k++) {
        info->AGP_exist[k] = info->AGP_level[k];
		info->layer_exist[k] = info->layer_mv[k]; 
  }

  strncpy( strtmp, info->bitname, strlen( info->bitname ) - 4 ); // Hanke, 16.09.02 
  strtmp[strlen( info->bitname ) - 4] = '\0';                    // 
  sprintf( info->mvname, "%s.mv", strtmp );                      //

  error_check( *info );
}
