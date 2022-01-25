#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#define EXTERN
#include "basic.h"
#include "structN.h"
#include "coderN.h"
#include "initN.h"
#include "ioN.h"
#include "memoryN.h"
#include "miscN.h"
#include "init_decN.h"
#include "mvcodingN.h"
#include "analsyn.h"
#include "pstatN.h"
#include "chrom.h"
#include "mv_statistics.h"

void ezbc3d_dec_GOP( YUVimage * pyrFrs, videoinfo info, long total_bytes_past,
                     long int GOP_counter, int curr );

void read_command( int argc, char **argv, videoinfo * info );

void print_time( double sc );

float global_motion_active;

/*
 *                                main()                                    
 */
int
main( int argc, char **argv )
{
  int curr, last, remaining_frs, s_level; 
  long mark, elp;  // initial and elapsed time
  double duration;
  FILE *fp_stat, *fpio;
  videoinfo info;
  long int total_bytes_past = 0, total_bytes_past_buffer;
  long int *read_GOP_bytes;
  long int num_of_GOP, GOP_counter;
  enum FLAG first_GOP, Level_change;
  char mvstatname[512];

  int gop_psnr_start;
  int simul_count; //Added on 01.19.2018 估计的数量

  int cx, cy, cpos;
  FILE *pc;

  printf("Built on %s at %s.\n", __DATE__, __TIME__);

  mark = clock(  );

  read_command( argc, argv, &info );// 读命令行

  read_header( info.bitname, &info );
  total_bytes_past += sizeof( videoheader );

  /* open MV statistics file 打开MV统计文件，以便写入统计信息 */
  strcpy(mvstatname, info.statname);
  strcat(mvstatname, "_mvstatistics.log");
  mvStat_open(mvstatname);

  s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
  // save the original resolution for frequency roll-off by Yongjun Wu 
  info.org_yheight = info.yheight; 
  info.org_ywidth  = info.ywidth; 
  if( s_level > 0 ) {
    info.ywidth  >>= s_level;
    info.yheight >>= s_level;
    info.cwidth  >>= s_level;
    info.cheight >>= s_level;
  }

  switch ( info.format ) {
  case YUV:
  case RAS:
    info.pixeldepth = 8;
    break;
  case DPX:
    info.pixeldepth = 10;
    break;
  default:
    printf( "image format error format = %d(pstatN.c)\n", info.format );
    exit( 1 );
  }

  
  // 初始化
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

  num_of_GOP = get_GOP_num( info );
  // printf( "number of GOPs = %d\n", (int) num_of_GOP );

  read_GOP_bytes =
    ( long int * )getarray( num_of_GOP, sizeof( long int ),
                            "read_GOP_bytes" );
  if( !( fpio = fopen( info.bitname, "rb" ) ) ) {
    printf( "can not open: %s\n", info.bitname );
    exit( 1 );
  }
  fseek( fpio, total_bytes_past, SEEK_SET );
  fread( read_GOP_bytes, sizeof( long int ), num_of_GOP, fpio );

  fclose( fpio );
  total_bytes_past += num_of_GOP * sizeof( long int );  //printf("%d\n", total_bytes_past);

  if( !( fp_stat = fopen( info.statname, "wt" ) ) ) {
    printf( "Can not open %s\n", info.statname );
    exit( 1 );
  }
  fprintf( fp_stat, "\n\tDecoding %.3d-%.3d\n", info.start, info.last );
  fprintf( fp_stat, "\tY: %d x %d; C: %d x %d\n", info.ywidth,
           info.yheight, info.cwidth, info.cheight );
  fprintf( fp_stat, "\tOriginal video frame rate: %d fps\n", info.framerate );
  fprintf( fp_stat, "\tGOP size: %d frs\n", info.GOPsz );
  fprintf( fp_stat, "\tHeader size: %d bytes\n\n", sizeof( videoheader ) );
  fclose( fp_stat );

  init_dec( info );

  Level_change = NO; // full GOP 
  first_GOP = YES;
  curr = info.start;
  last = info.last;
  GOP_counter = 0;


  while( curr <= last ) {       /* MCTF decoding gop级别*/
    remaining_frs = last - curr + 1;               
 
    if ( remaining_frs < info.GOPsz ){
	  Level_change = YES;
	  info.eff_GOPsz = remaining_frs;	
	  // printf( "******* read last GOP (curr %d, eff_GOPsz %d, remaining_frs %d) *******\n", 
      //         curr, info.eff_GOPsz, remaining_frs );	
    } else { 
	  Level_change = NO;            // full GOP   
	  info.eff_GOPsz = info.GOPsz;  // effective GOP size
    }    
      
    if ( first_GOP == YES )  
	  printf( " decoding frame %d (first_GOP) .....\n", curr );
    else 
	  printf( " decoding frame %d - %d .....\n", curr - info.GOPsz + 1, curr );
    
#ifdef SUPPORT_INTRA
    switch ( info.intra ) {
    case YES:
        
      if( info.denoise_flag == YES ) {
        printf( "In intraframe coding mode, denoise = YES is not allowed\n" );
        exit( 1 );
      }

      info.GOPbytes = read_GOP_bytes[GOP_counter];

#ifdef THREE_D
      three_D_syn( curr, GOP_counter, &total_bytes_past, info );
#else
      intra_decode( curr, GOP_counter, &total_bytes_past, info );
#endif
      total_bytes_past += info.GOPbytes;
      // calsnr_seq(&cfr, &pfr, curr, curr+nFrsPyr-1, info); /* pstatN.c */
      /* cfr and pfr are just used as the memory in calsnr, and the data in cfr and pfr is not used in calsnr */

      break;

    case NO:
#endif
      info.GOPbytes = read_GOP_bytes[GOP_counter]; // 一个gop的字节数
	  printf("newfat = %d\n",info.GOPbytes);

	  simul_count = 0;

//	  gop_block = 0;	frame_cnt = 0;

      if( info.denoise_flag == YES ) {
        
        denoise_mctf_syn_ezbc( curr, GOP_counter, &total_bytes_past, info,
                               first_GOP, Level_change, remaining_frs);
      } else {
        
        mctf_syn_ezbc( curr, GOP_counter, &total_bytes_past, info,
                       first_GOP, Level_change, remaining_frs, NO, NO ); // 进行解码
      }
      
//	  printf("gop_block = %f\n",gop_block/frame_cnt);

      total_bytes_past += info.GOPbytes;
//	  printf("info.GOPbytes = %d\n",info.GOPbytes);
      // calsnr_seq(&cfr, &pfr, curr, curr+info.GOPsz-1, info); /* pstatN.c */
      /* cfr and pfr are just used as the memory in calsnr, and the data in cfr and pfr is not used in calsnr */
#ifdef SUPPORT_INTRA
      break;

	default:
      printf( "error in decoderN.c\n" );
      exit( 1 );
    }
#endif

    first_GOP = NO;
    GOP_counter++;
    curr += info.GOPsz;
    
  } // while
  
  // decode LAST GOP // 解码最后一个gop
  {
    info.eff_GOPsz = remaining_frs;
    remaining_frs = 0; // no bitstream extraction of last GOP, only temporal filtering
    
    printf( " decoding frame %d - %d (last_GOP) .....\n", curr - info.GOPsz + 1, last );
          
    info.GOPbytes = read_GOP_bytes[GOP_counter];

    if( info.denoise_flag == YES ) {
      denoise_mctf_syn_ezbc( curr, GOP_counter, &total_bytes_past, info,
                             first_GOP, Level_change, remaining_frs);
    } else {
      mctf_syn_ezbc( curr, GOP_counter, &total_bytes_past, info,
                     first_GOP, Level_change, remaining_frs, NO, NO );
    }
    //    total_bytes_past += info.GOPbytes;
  }

  gop_psnr_start = info.last - info.GOPsz + 1;

  if ( info.s_level == 0 ){ // Quarter == 0 
    if( info.denoise_flag == YES ){
      info.ywidth  *= 2;
      info.yheight *= 2;
      info.cwidth  *= 2;
      info.cheight *= 2;
      
	  calsnr( info.start, info.last, info );
      
      info.ywidth  /= 2;
      info.yheight /= 2;
      info.cwidth  /= 2;
      info.cheight /= 2;
    } 
	else 
	{
      calsnr( info.start, info.last, info );
    }
  }
  
  // }

  info.bitrate =
    ( int )( 8 * total_bytes_past *
             ( ( float )info.framerate / ( info.last - info.start + 1 ) ) );
  if( !( fp_stat = fopen( info.statname, "at+" ) ) ) {
    printf( "Can not open %s\n", info.statname );
    exit( 1 );
  }
  fprintf( fp_stat, "\n\t total bytes %d, rate %d bps \n",
           (int) total_bytes_past, info.bitrate );
  fclose( fp_stat );

  elp = clock(  ) - mark;
  duration = ( double )elp / CLOCKS_PER_SEC;
  print_time( duration );

  mvStat_close();

#ifdef BLOCKMODE_STATISTICS
  dump_blockmode_statistics(info.statname);
#endif

  printf( "finished.\n" );
  free(read_GOP_bytes);

  free(buff_frameMEinfo);
  return 0;
}



/*
 *                                usage()                                    
 */
void
usage(  )
{
  printf( "decoder: bitfile decname inname statname  \n" );
}

void
read_command( int argc, char **argv, videoinfo * info )
{
  int i, argnum = 1;
  for( i = 1; i < argc; i++ ) {
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
        strcpy( info->bitname, argv[i] );
        argnum++;
        break;
      case 2:
        strcpy( info->decname, argv[i] );
        argnum++;
        break;
      case 3:
        strcpy( info->inname, argv[i] );
        argnum++;
        break;
      case 4:
        strcpy( info->statname, argv[i] );
        argnum++;
        break;
	  case 5:
        strcpy( info->jp2k_decname, argv[i] );
        argnum++;
        break;
#ifdef CNN_wavelet
	  case 6:
		  strcpy(info->hpcbindata, argv[i]);
		  argnum++;
		  break;
#endif
      }
    }
  }
#ifdef CNN_wavelet
  if (argc != 7) {
#else
  if (argc != 6) {
#endif
    usage(  );
    exit( 1 );
  }

}
