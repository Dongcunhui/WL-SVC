
#ifndef __EZBC_CODEC_3D_H__
#define __EZBC_CODEC_3D_H__

#include "dwt_bitplane_codec.h"
#include "struct_sht.h"
#include "image_bw.h"

/*
 * ncomps: number of color components
 * nbands[3]: number of spatial subbands of each color component
 * num_Slev[3]: number of spatial pyramid levels of each color component
 * subband_ACcoder[color_component][frame][spatial subband]: which AC coder will be used
 * max_Slev:
 */


class EzbcCodec3d
{
protected:
  int GOPsz, ncomps, max_Slev;
  int **subband_ACcoder[3];
  int total_num_ACcoder;

  Image_Coord dim, cdim; // dim为亮度的大小， cdim为色度的大小
  int nbands[3], num_Slev[3];   //initialized in initialize() 空域小波分解：三个通道的空域子带个数， 三个通道的空域分解次数+1

  //long byte_budget;
  SUBBAND_TREE_TYPE **pyrY, **pyrU, **pyrV;
  Image_BW *imgY, *imgU, *imgV;
public:
   int s_level;  // by Yongjun Wu
   enum FRAME_RES  frame_res; 
   

    EzbcCodec3d( void );
    EzbcCodec3d( videoinfo & info );
    void initialize( void );
   ~EzbcCodec3d( void );
  //void set_byte_budget(long bytes){ byte_budget = bytes;}
};

#endif
