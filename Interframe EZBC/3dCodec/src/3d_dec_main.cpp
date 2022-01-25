/* ========================================================================= */
/* Description: main() for decoding 3D subband coefficients                  */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */
#include "structN.h"
#include "basic.h"
#include "general.h"
#include "image_bw.h"
#include "subband.h"
#include "ar_code.h"
#include "dwt_bitplane_enc.h"
#include "video_utils.h"
#include "ezbc_dec_3d.h"
#include "miscN.h"
#include "util_filtering.h"

#include "ioN.h"

/*
 * ezbc3d_dec  获取一帧
 * pyrFrs: buffer for reconstructed frame, allocated inside this function 在这个函数中分配
 * info: sequence header
 * total_bytes_past: bytes needed to be skipped before decoding the bit stream 
 * GOP_counter: GOP index
 * return bytes_past: bytes for subbands including subband headers
 */
long int
ezbc3d_dec( int curr, YUVimage * pyrFrs, videoinfo info, long int total_bytes_past,
            long int GOP_counter )
{
  DECODER_TYPE decoder;
  long int bytes_past;
  // Chronometer code_time, total_time;
  Chronometer total_time;

  total_time.start(  );// 计时

  EzbcDec3d video_dec( info, &decoder );
  SubbandCodec::setup_luts(  );// 初始化上下文表

  video_dec.GOP_reset(  );// 重置清零
  video_dec.reset_pyrs(  );// 重建子带内的金字塔结构

  video_dec.initialization(  );
  //video_dec.set_byte_budget(info.GOPbytes );

  bytes_past = video_dec.subband_dec_pointer( info.bitname, total_bytes_past,
                                   info.s_level );

//  printf("\nout of subband_dec_pointer!\n");

  video_dec.decode_GOP(  );

  video_dec.GOP_mean_reset(  );

  video_dec.reconstruct_subbands(  );

  video_dec.free_lists(  );// 释放掉了

  video_dec.reconstruct_images( pyrFrs, curr );// 写进pyfrs

  SubbandCodec::delete_cxt_tables(  );

  //total_time.display("\n  Total execution time (I/O included) =");
  //puts(" ");
  return ( bytes_past );
}

/*
 * ezbc3d_dec
 * 一个GOP的帧重建
 */
void
ezbc3d_dec_GOP( YUVimage * pyrFrs, videoinfo info, long int total_bytes_past,
                long int GOP_counter, int curr )
{
  int i, j, k,  GOPsz, eff_GOPsz, count, dist, s_level;
  int temp_pos[1024]; // maximum number of temporal levels = 10; 时域索引值
  int band_idx;

  int start = 0;

  int jp2k_count;

  printf("ezbc3d_dec_GOP, t_level = %d, curr = %d\n",info.t_level,curr);
 
  GOPsz = info.GOPsz;  

  if ( info.eff_GOPsz < info.GOPsz )  
    eff_GOPsz = info.eff_GOPsz;
  else 
    eff_GOPsz = GOPsz;

  printf("eff_GOPsz = %d, GOPsz = %d\n",eff_GOPsz,GOPsz);
  
  // get temporal positions of pyrFrs + save
  temp_pos[0] = 0; // first lowpass frame
  temp_pos[1] = (int) pow (2.0, info.tPyrLev - 1); // first highpass frame 以3层为例，temp_pos[1] = 4
  count = 2; 
  for( i = info.tPyrLev - 1; i > 0; i-- ){ // i = 2 1 
    dist = (int) pow (2.0, i); // 4 2
    for (j = 0; j < (info.GOPsz / dist); j++){ // 0 1 、0 1 2 3
	  temp_pos[count] = dist / 2 + j * dist; // temp_pos[3] = 2 temp_pos[4] = 2+4 
	  //printf("%d\n", temp_pos[count]);  2 6 1 3 5 7
	  count++;
    } 
  }

  jp2k_count = 0;

  //skip temporal high subband frames; if current GOP < default GOP size, GOPsz could be 0
  GOPsz = info.GOPsz;     // save the GOPsz 记录
  info.GOPsz = 1;         // change info.GOPsz to be 1  gopsz从1逐渐向上
  s_level = info.s_level; // save s_level 记录
  info.s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));

  for( j = 0; j < GOPsz >> info.t_level; j++ ){// 进行了几次时域分解
    if ( temp_pos[j] < eff_GOPsz ){  // 取eff_GOPsz个帧

#ifndef JP2K_SUBBAND
      total_bytes_past += ezbc3d_dec( curr + j, &pyrFrs[j], info, total_bytes_past, GOP_counter ); // 解码一帧
#else
//	  printf("frame j = %d\n",j);
	  read_jp2k_frame( &pyrFrs[j], info, info.jp2k_decname, curr - info.start + jp2k_count, info.format );
	  jp2k_count ++;
#endif

#ifdef PREQUANT_WEIGHTING
      // post-quantization weighting 量化  放在pyrFrs
      assert(!(GOPsz % (1 << info.tPyrLev)));
      if( j < (GOPsz >> info.tPyrLev) ) // low frequency subband
	  {
        band_idx = info.tPyrLev - 1;
        wcopyframe(&pyrFrs[j], &pyrFrs[j], 1.0f / static_prequant_weight_low[band_idx],
                   info);
      } 
	  else // high frequency subband
	  {
        band_idx = int( floor( log2( GOPsz ) - log2( j + 1 ) ) );
        assert(band_idx >= 0 && band_idx < info.tPyrLev);
        wcopyframe(&pyrFrs[j], &pyrFrs[j], 1.0f / static_prequant_weight_high[band_idx],
                   info);
      }
#endif

      // memory test
      // alloc_pyr_frames( &pyrFrs[i * GOPsz + j], info );
    }
  } 
  info.s_level = s_level;

  // 去噪
  if ( info.denoise_flag == YES && info.s_level == 0){ // Quarter == 0
    for( i = 1; i < (( YUV420 == 1 ) ? 2 : 4); i++ ){
      for( j = 0; j < eff_GOPsz; j += (1 << info.t_level) ){
		assert(0);
#ifndef JP2K_SUBBAND
        total_bytes_past +=
          ezbc3d_dec( curr + j, &pyrFrs[i * GOPsz + j], info, total_bytes_past, GOP_counter );
#endif

      }
    }
  }
  
  info.GOPsz = GOPsz;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// end of file < CodeTree.C >
