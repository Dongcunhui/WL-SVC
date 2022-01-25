
#ifndef __EZBC_ENC_3D_H__
#define __EZBC_ENC_3D_H__

#include "dwt_bitplane_enc.h"
#include "ezbc_codec_3d.h"

long int ezbc3d_enc_GOP( int curr, YUVimage * pyrFrs, videoinfo info,
                         long int GOP_counter );

class EzbcEnc3d:public EzbcCodec3d
{

protected:

  int GOP_max_msb, max_depth, total_subplane;   // chen
  ENCODER_TYPE *encoder;

  ENC_SUBBAND_TREE_TYPE *encY, *encU, *encV;    // typedef EncSubbandTree ENC_SUBBAND_TREE_TYPE; in EZBC0a/dwt_bitplane_enc.h

public:

  int *GOPmean[3], *GOPmean_shift[3];

    EzbcEnc3d( void );
    EzbcEnc3d( videoinfo & info, ENCODER_TYPE * enc );

   ~EzbcEnc3d( void );
  void load_images( YUVimage * Frs, videoinfo &info, float *gop_mse, int curr );
  void initialization( void );
  void reset_GOP_enc( void );
  void GOP_transform_mean( void );
  //  void encode_GOP_single_pass(void);
  void encode_GOP_RD_passes( void );
  void encode_GOP_RD_passes2( void );
  void encode_GOP( void );
  void free_lists( void );
  long int get_byte_used( void )
  {
    return encoder->bytes_used(  );
  };
  long int write_files( long int GOP_counter, char *bitname );
};

#endif
