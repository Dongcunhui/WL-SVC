/* ========================================================================= */
/* Description: main() for encoding 3D subband coefficients                  */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <iostream>
#include "general.h"
#include "image_bw.h"
#include "structN.h"
#include "subband.h"
#include "ar_code.h"
#include "dwt_bitplane_enc.h"
#include "video_utils.h"
#include "ezbc_enc_3d.h"
#include "util_filtering.h"
#include "miscN.h"

#include "ioN.h"

// #define DEBUG_MCTF_RESULTS   // for debug use by Yongjun Wu
/*
 *  3dezbc_enc()
 * curr: start frame index of this GOP
 * Frs : frame buffer
 * info: header
 * GOP_counter: GOP index
 * return total_bytes: total number of bytes generated, including subband header
 逐帧去编码
 */
long int
ezbc3d_enc( int curr, YUVimage * Frs, videoinfo info, long int GOP_counter, float *gop_mse )
{
//  printf("New Frame\n\n");
  ENCODER_TYPE encoder;         //  EZBC0a/dwt_bitplane_enc.h:typedef  Encoder  ENCODER_TYPE; Utils/ar_code.h:class Encoder 
  //int header_bytes;
  long int total_bytes = 0;

  // Chronometer code_time, total_time;
  Chronometer total_time;

  Encoder::object_count = 0;

  total_time.start(  );

  EzbcEnc3d video_enc( info, &encoder );        /*"class EzbcEnc3d: public EzbcCodec3d" in ezbc_enc_3d.h, class EzbcCodec3d defined in ezbc_codec_3d.h */

  SubbandCodec::setup_luts(  ); // defined in EZBC0a/dwt_bitplane_codec.h

  video_enc.load_images( Frs, info, gop_mse, curr ); // 导入加变换 transform EzbcEnc3d video_enc(info, &encoder);  load_images defined in ezbc_enc_3d.C

  //free_frames( Frs, info.GOPsz );

  //video_enc.set_byte_budget(info.GOPbytes);
  video_enc.GOP_transform_mean(  ); // 设置GOP的均值

  video_enc.initialization(  );

  video_enc.encode_GOP(  ); // 编码gop

  total_bytes = video_enc.write_files( curr, info.bitname );
  video_enc.free_lists(  );

  SubbandCodec::delete_cxt_tables(  );  //same as release_3dezbc_enc() below.

  return ( total_bytes );
}



void
release_3dezbc_enc(  )
{
  SubbandCodec::delete_cxt_tables(  );
}

/*
 *  3dezbc_enc()
 * curr: start frame index of this GOP
 * Frs : frame buffer
 * info: header
 * GOP_counter: GOP index
 * return total_bytes: total number of bytes generated, including subband headers
 将mc后的帧进行编码
 */
long int
ezbc3d_enc_GOP( int curr, YUVimage * pyrFrs, videoinfo info,
                long int GOP_counter )
{
  int i, j, k, GOPsz, eff_GOPsz, dist, count, GOP_multiplier;
  long int total_bytes = 0;
  int temp_pos[1024]; // maximum of 10 temporal levels
  int band_idx;

  float getval, gop_mse;

  if ( info.denoise_flag == YES ){
    GOP_multiplier = ( YUV420 == 1 ) ? 2 : 4;
  } else {
    GOP_multiplier = 1;
  }
  
  if ( info.eff_GOPsz < info.GOPsz ) // info.GOPsz here! 
    eff_GOPsz = info.eff_GOPsz; // = remaining_frs!
  else 
    eff_GOPsz = info.GOPsz;

  // get temporal positions of pyrFrs + save  获得帧的位置，暂时还没看懂
  temp_pos[0] = 0; // first lowpass frame
  temp_pos[1] = (int) pow (2.0, info.tPyrLev - 1); // first highpass frame
  count = 2; 
  for( i = info.tPyrLev - 1; i > 0; i-- ){
    dist = (int) pow (2.0, i);
    for (j = 0; j < (info.GOPsz / dist); j++){
	  temp_pos[count] = dist / 2 + j * dist;
	  count++;
    } 
  }
  GOPsz = info.GOPsz; // save the GOPsz
  info.GOPsz = 1;     // change info.GOPsz to be 1

  gop_mse = 1;
  
  // 逐帧进行ezbc编码
  for( j = 0; j < GOPsz; j++ ){

    if ( temp_pos[j] < eff_GOPsz ){

#ifdef DEBUG_MCTF_RESULTS
		// for debug use
		write_frame_into_file(pyrFrs[j], info, j, GOP_counter); 
#endif

#ifdef PREQUANT_WEIGHTING
      // pre-quantization weighting
      assert(!(GOPsz % (1 << info.tPyrLev)));
      if( j < (GOPsz >> info.tPyrLev) ) {
        // low frequency subband
        band_idx = info.tPyrLev - 1;

		printf("static_prequant_weight_low[%d] = %f\n\n",band_idx,static_prequant_weight_low[band_idx]);

        wcopyframe(&pyrFrs[j], &pyrFrs[j], static_prequant_weight_low[band_idx], info);
      } 
	  else {
        // high frequency subband

        band_idx = int( floor( log2( GOPsz ) - log2( j + 1 ) ) );

//		printf("j = %d, log2(j+1) = %f, log2(GOPsz) = %d, log2(GOPsz) - log2(j+1) = %f, floor(log2(GOPsz) - log2(j +1) ) = %f\n\n",j,log((float)j+1)/log(2),log2(GOPsz),
//			log2(GOPsz) - log((float)j+1)/log(2), floor(log2(GOPsz) - log((float)j+1)/log(2) ) );

        assert(band_idx >= 0 && band_idx < info.tPyrLev);
		
		printf("static_prequant_weight_high[%d] = %f\n\n",band_idx,static_prequant_weight_high[band_idx]);

        wcopyframe(&pyrFrs[j], &pyrFrs[j], static_prequant_weight_high[band_idx], info);
      }
#endif

#ifdef JP2K_SUBBAND
	write_jp2k_frame(pyrFrs[j],info,info.jp2kname,curr + j,info.format);
#else
      total_bytes += ezbc3d_enc( curr + j, &pyrFrs[j], info, GOP_counter, &gop_mse );
#endif
    }
  }

//  printf("\nGOP sum mse = %f\n\n",gop_mse);
  gop_mse = pow(gop_mse, (float)((float)1/(float)GOPsz));
//  printf("\nGOP avg mse = %f\n\n",gop_mse);

  // 下面是为了denoise模式
  for( i = 1; i < GOP_multiplier; i++ ){
    for( j = 0; j < eff_GOPsz; j++ ){

	  assert(0);

      // for sH subbands (denoise) no temporal decomposition is applied,
      // thus not pre-quantization weighting is necessary
#ifndef JP2K_SUBBAND
      total_bytes += ezbc3d_enc( curr + j, &pyrFrs[i * GOPsz + j], info, GOP_counter, &gop_mse );
#endif
    }
  }
  
  info.GOPsz = GOPsz;
  return ( total_bytes );
}


// end of file < CodeTree.C >
