#ifndef __EZBC_DEC_3D_H__
#define __EZBC_DEC_3D_H__

#include "dwt_bitplane_dec.h"
#include "ezbc_codec_3d.h"

class EzbcDec3d:public EzbcCodec3d
{
protected:
  DECODER_TYPE * decoder;
  DEC_SUBBAND_TREE_TYPE *decY, *decU, *decV;
  int max_depth;
public:
  char **subband_mask[3];       // allocated in subband_dec_pointer and freed in ~EzbcDec3d
  long int byte_used;
  int *GOPmean[3], *GOPmean_shift[3];

    EzbcDec3d( void );
    EzbcDec3d( videoinfo & info, DECODER_TYPE * dec );
   ~EzbcDec3d( void );
  void initialization( void );
  long int subband_dec_pointer( char *file_name, long int total_bytes_past,
                                short int s_level );
  void reset_GOP_dec( void );
  void reset_pyrs( void );
  void reconstruct_images( YUVimage * Frs, int curr );
  void GOP_reset( void );
  void GOP_mean_reset( void );
  void decode_GOP_header( void );
  //  void decode_GOP_single_pass(void);
  void decode_GOP_RD_passes( void );
  void decode_GOP_RD_passes2( void );
  void end_of_decoding( int subband_bytes_used, int subband_byte_budget,
                        int m, int k, int i, char *subband_flag,
                        long int *byte_used );
  void decode_GOP( void );
  void reconstruct_subbands( void );
  void free_lists( void );
};



#endif
