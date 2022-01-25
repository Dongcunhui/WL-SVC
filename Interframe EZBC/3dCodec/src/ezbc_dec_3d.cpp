/* in decode_GOP_RD_passes2(), I did some statistic work on used bytes for LL bands*/

#include "structN.h"
#include "ezbc_dec_3d.h"
#ifdef NDEBUG
#define myassert(seq) seq
#else
#define myassert(seq) assert(seq)
#endif

int read_substream_length( int *length, FILE * fp );


EzbcDec3d::EzbcDec3d( videoinfo & info, DECODER_TYPE * dec )
:EzbcCodec3d( info ), decoder( dec )
{
  decY = decU = decV = NULL;

  NEW_VECTOR( GOPmean[0], GOPsz, int, "GOPmean[0]" );
  NEW_VECTOR( GOPmean_shift[0], GOPsz, int, "GOPmean_shift[0]" );

  if( cdim.x ) {
    NEW_VECTOR( GOPmean[1], GOPsz, int, "GOPmean[1]" );
    NEW_VECTOR( GOPmean_shift[1], GOPsz, int, "GOPmean_shift[1]" );

    NEW_VECTOR( GOPmean[2], GOPsz, int, "GOPmean[2]" );
    NEW_VECTOR( GOPmean_shift[2], GOPsz, int, "GOPmean_shift[2]" );
  }

  NEW_VECTOR( decY, GOPsz, DEC_SUBBAND_TREE_TYPE, "decY" );
  if( cdim.x > 0 ) {
    NEW_VECTOR( decU, GOPsz, DEC_SUBBAND_TREE_TYPE, "decU" );
    NEW_VECTOR( decV, GOPsz, DEC_SUBBAND_TREE_TYPE, "decV" );
  } else {
    decU = decV = NULL;
  }

}

EzbcDec3d::~EzbcDec3d( void )
{
  if( decY )
    DELETE_VECTOR( decY );
  if( decU )
    DELETE_VECTOR( decU );
  if( decV )
    DELETE_VECTOR( decV );

  DELETE_VECTOR( GOPmean[0] );
  DELETE_VECTOR( GOPmean_shift[0] );

  if( cdim.x ) {
    DELETE_VECTOR( GOPmean[1] );
    DELETE_VECTOR( GOPmean_shift[1] );

    DELETE_VECTOR( GOPmean[2] );
    DELETE_VECTOR( GOPmean_shift[2] );

  }


  for( int m = 0; m < ncomps; m++ ) {
    for( int i = 0; i < GOPsz; i++ ) {
      DELETE_VECTOR( subband_mask[m][i] );
    }
    DELETE_VECTOR( subband_mask[m] );
  }

  if( decoder )
    delete[]decoder;
}

void
EzbcDec3d::reset_pyrs( void )
{
  int i, nbds = 1 + 3 * imgY[0].pyramid_levels(  );
  int *null_msb_layers;



  NEW_VECTOR( null_msb_layers, nbds, int, "null_msb_layers" );

  for( i = 0; i < GOPsz; i++ ) { //GOPsz = 1
    if( pyrY[i] )
		DELETE_OBJECT( pyrY[i] );
     myassert( pyrY[i] = new SUBBAND_TREE_TYPE( imgY[i], s_level, null_msb_layers ) );
     imgY[i].dispose(  );
    if( cdim.x ) {
      if( pyrU[i] )
        DELETE_OBJECT( pyrU[i] );
      myassert( pyrU[i] = new SUBBAND_TREE_TYPE( imgU[i], s_level, null_msb_layers ) );

      imgU[i].dispose(  );
      if( pyrV[i] )
        DELETE_OBJECT( pyrV[i] );

      myassert( pyrV[i] = new SUBBAND_TREE_TYPE( imgV[i], s_level, null_msb_layers ) );
      imgV[i].dispose(  );
    }
  }
  delete[]null_msb_layers;

}

void
EzbcDec3d::initialization( void )
{
  int i, j;
  DEC_SUBBAND_TREE_TYPE *decGOP[3] = { decY, decU, decV };
  SUBBAND_TREE_TYPE **GOPpyr[3] = { pyrY, pyrU, pyrV };

  EzbcCodec3d::initialize(  );


  for( i = 0; i < GOPsz; i++ ) {
    for( j = 0; j < ncomps; j++ ) {
      decGOP[j][i].SubbandTreeCodec::initialize( GOPpyr[j][i] );
    }
  }

  NEW_VECTOR( decoder, total_num_ACcoder, DECODER_TYPE, "decoder (dec_subband_tree.cpp)" );     //062602
  for( i = 0; i < GOPsz; i++ ) {
    for( j = 0; j < ncomps; j++ ) {
      decGOP[j][i].initialize( GOPpyr[j][i], subband_ACcoder[j][i], decoder );
    }
  }

}

void
EzbcDec3d::reset_GOP_dec( void )
{
  int i;

  for( i = GOPsz - 1; i >= 0; i-- ) {
    decY[i].reset_tree_dec( pyrY[i] );
    if( cdim.x > 0 ) {
      decU[i].reset_tree_dec( pyrU[i] );
      decV[i].reset_tree_dec( pyrV[i] );
    }
  }
}

void
EzbcDec3d::GOP_reset( void )
{
  int i;

  for( i = 0; i < GOPsz; i++ ) {
    imgY[i].reset( dim, 1, 0, 0, 0 );

    if( cdim.x ) {
      imgU[i].reset( cdim, 1, 0, 0, 0 );

      imgV[i].reset( cdim, 1, 0, 0, 0 );
    }
  }

}

void
EzbcDec3d::GOP_mean_reset( void )
{
  int i;

  for( i = 0; i < GOPsz; i++ ) {
    imgY[i].reset( GOPmean[0][i], GOPmean_shift[0][i], 0 );

    if( cdim.x ) {
      imgU[i].reset( GOPmean[1][i], GOPmean_shift[1][i], 0 );

      imgV[i].reset( GOPmean[2][i], GOPmean_shift[2][i], 0 );
    }
  }

}



void
EzbcDec3d::decode_GOP( void )
{
  decode_GOP_RD_passes2(  );
}



void
EzbcDec3d::reconstruct_subbands( void )
{
  int i;

  for( i = GOPsz - 1; i >= 0; i-- ) {
    decY[i].rec_subbands(  );
    if( cdim.x ) {
      decU[i].rec_subbands(  );
      decV[i].rec_subbands(  );
    }
  }
}

// 输入时域的buf， 当前帧  将img中的数据写进Frs
void
EzbcDec3d::reconstruct_images( YUVimage * Frs, int curr )
{
  int i, mean, mean_shift;

  if( 0 )
  {
	  for( i = GOPsz - 1; i >= 0; i-- ) {
		  imgY[i].high_band = YES;
		  imgU[i].high_band = YES;
		  imgV[i].high_band = YES;
	  }
  }
  else
  {
	  for( i = GOPsz - 1; i >= 0; i-- ) 
	  {
		  imgY[i].high_band = NO;
		  imgU[i].high_band = NO;
		  imgV[i].high_band = NO;
	  }
  }

  for( i = GOPsz - 1; i >= 0; i-- ) {
    mean = imgY[i].transform_mean(  );
    mean_shift = imgY[i].mean_shift(  );
    imgY[i].reset( dim, 1, mean, mean_shift, 0 );
    pyrY[i]->SubSetLayer2Image_BW( imgY[i] );

	imgY[i].s_level = s_level; // by Yongjun Wu
	imgY[i].frame_res = frame_res; 
    imgY[i].recover(  );
    imgY[i].write_float_image( Frs[i].Y );
    imgY[i].dispose(  );

    if( cdim.x ) {
      mean = imgU[i].transform_mean(  );
      mean_shift = imgU[i].mean_shift(  );
      imgU[i].reset( cdim, 1, mean, mean_shift, 0 );
      pyrU[i]->SubSetLayer2Image_BW( imgU[i] );

	  imgU[i].s_level = s_level;   // by Yongjun Wu
	  imgU[i].frame_res = frame_res; 
      imgU[i].recover(  );
      imgU[i].write_float_image( Frs[i].U );
      imgU[i].dispose(  );

      mean = imgV[i].transform_mean(  );
      mean_shift = imgV[i].mean_shift(  );
      imgV[i].reset( cdim, 1, mean, mean_shift, 0 );
      pyrV[i]->SubSetLayer2Image_BW( imgV[i] );

	  imgV[i].s_level = s_level;   // by Yongjun Wu
	  imgV[i].frame_res = frame_res; 
      imgV[i].recover(  );
      imgV[i].write_float_image( Frs[i].V );
      imgV[i].dispose(  );
    }
  }
}


/*
 * subband_dec_pointer
 * total_bytes_past: bytes needed to be skipped before decoding the bit stream
 * return bytes_past: bytes for all subbands including subband headers
 */
long int
EzbcDec3d::

subband_dec_pointer( char *file_name, long int total_bytes_past,
                     short int s_level )
{
  int i, m, k, past;
  long int bytes_past = 0;
  int subband_bytes = 0;
  DEC_SUBBAND_TREE_TYPE *decGOP[3] = { decY, decU, decV };
  DecSubband *dec_band;
  FILE *fpbit;
  int s; // by Yongjun Wu

  if( ( fpbit = fopen( file_name, "rb" ) ) == NULL ) {
    printf( "cannot open file %s\n", file_name );
    exit( 1 );
  }

  /*
   * generate subband mask
   */
  for( m = 0; m < ncomps; m++ ) {
    NEW_VECTOR( subband_mask[m], GOPsz, char *,
                "subband_mask[m](ezbc_dec_3d.cpp)" );
    for( i = 0; i < GOPsz; i++ ) {

      NEW_VECTOR( subband_mask[m][i], nbands[m], char,
                  "subband_mask[m][i] (ezbc_dec_3d.cpp)" );
      for( k = 0; k < nbands[m]; k++ )
        subband_mask[m][i][k] = 0;
      /*for(k = 0; k < nbands[m]-s_level*3; k++)
         subband_mask[m][i][k] = 0;
         for(k = nbands[m]-s_level*3; k < nbands[m]; k++)
         subband_mask[m][i][k] = 1; */
    }
  }

  /*
   * direct each AC decoder to corresponding subbitstream
   */
  int roll_level1, roll_level2; 
  if (frame_res == HD )   //Added by Yuan Liu
  {
	  roll_level1 = 5; roll_level2 = 4;  // for HD sequence both CIF and QCIF need frequency roll-off 
  }else if (frame_res == SD2)
  {
	  roll_level1 = 5; roll_level2 = 4;  // for SD sequence both CIF and QCIF need frequency roll-off 
  }else if (frame_res == SD )
  {
	  roll_level1 = 4; roll_level2 = 4;  // for SD sequence both CIF and QCIF need frequency roll-off 
  }else if (frame_res == CIF )
	  roll_level1 = roll_level2 = 4;     // for CIF sequence, only QCIF needs frequency roll-off 

  for( i = 0; i < GOPsz; i++ ) {
    for( m = 0; m < 1; m++ ) {  // since Y U V are grouped together

      for( k = 0; k < num_Slev[m]; k++ ) {
		  if( subband_mask[m][i][k * 3] == 0 ) {  // actually this "if" is not necessary

#ifdef    FREQUENCY_ROLL_OFF

#ifdef ROLL_STRUCTURE_ONE
			  if (k!=roll_level1 && k!=roll_level2)
			  {
				  dec_band = &(decGOP[m][i].dec_subs[k*3]);
				  if(fseek(fpbit, total_bytes_past, SEEK_SET) != 0){// chen
					  printf(" total_bytes_past = %d, fseek error(ezbc_dec_3d.cpp)\n", total_bytes_past);
					  exit(1);
				  }
				  past = read_substream_length(&subband_bytes, fpbit);
			      total_bytes_past += past;
			      bytes_past += past;  //subband header
//				  printf("open option 1!\n");
				  dec_band->sub_decoder->new_open_file(file_name, total_bytes_past, subband_bytes);
			      total_bytes_past += subband_bytes;  
   	 	          bytes_past += subband_bytes; // subband size (not include bytes for header)
			  }else
			  {
				  for (s=2; s>=0; s--)
				  {
					  dec_band = &(decGOP[m][i].dec_subs[k*3-s]);
					  if(fseek(fpbit, total_bytes_past, SEEK_SET) != 0){// chen
						printf(" total_bytes_past = %d, fseek error(ezbc_dec_3d.cpp)\n", total_bytes_past);
						exit(1);
					  }
					  past = read_substream_length(&subband_bytes, fpbit);
			          total_bytes_past += past;
			          bytes_past += past;  //subband header
//					  printf("open option 2!\n");
				      dec_band->sub_decoder->new_open_file(file_name, total_bytes_past, subband_bytes);
			          total_bytes_past += subband_bytes;  
   	 	              bytes_past += subband_bytes; // subband size (not include bytes for header)
				  }
			  }
#endif

#ifdef ROLL_STRUCTURE_TWO
			  if (k<roll_level2)
			  {
				  dec_band = &(decGOP[m][i].dec_subs[k*3]);
				  if(fseek(fpbit, total_bytes_past, SEEK_SET) != 0){// chen
					  printf(" total_bytes_past = %d, fseek error(ezbc_dec_3d.cpp)\n", total_bytes_past);
					  exit(1);
				  }
				  past = read_substream_length(&subband_bytes, fpbit);
			      total_bytes_past += past;
			      bytes_past += past;  //subband header
				  dec_band->sub_decoder->new_open_file(file_name, total_bytes_past, subband_bytes);
			      total_bytes_past += subband_bytes;  
   	 	          bytes_past += subband_bytes; // subband size (not include bytes for header)
			  }else if (k==roll_level2  || k==roll_level1)
			  {
				  int small_bands; 
				  if (k==roll_level1) small_bands = 12;  // 12 small subbands 
				  if (k==roll_level2) small_bands = 3;   // 3  small subbands 
				  for (s=0; s<small_bands; s++)
				  {
					  dec_band = &(decGOP[m][i].dec_subs[(k-1)*3+1+s]);
					  if(fseek(fpbit, total_bytes_past, SEEK_SET) != 0){// chen
						printf(" total_bytes_past = %d, fseek error(ezbc_dec_3d.cpp)\n", total_bytes_past);
						exit(1);
					  }
					  past = read_substream_length(&subband_bytes, fpbit);
			          total_bytes_past += past;
			          bytes_past += past;  //subband header
				      dec_band->sub_decoder->new_open_file(file_name, total_bytes_past, subband_bytes);
			          total_bytes_past += subband_bytes;  
   	 	              bytes_past += subband_bytes; // subband size (not include bytes for header)
				  }
			  }else
			  {
				  dec_band = &(decGOP[m][i].dec_subs[25]);
				  if(fseek(fpbit, total_bytes_past, SEEK_SET) != 0){// chen
					  printf(" total_bytes_past = %d, fseek error(ezbc_dec_3d.cpp)\n", total_bytes_past);
					  exit(1);
				  }
				  past = read_substream_length(&subband_bytes, fpbit);
			      total_bytes_past += past;
			      bytes_past += past;  //subband header
				  dec_band->sub_decoder->new_open_file(file_name, total_bytes_past, subband_bytes);
			      total_bytes_past += subband_bytes;  
   	 	          bytes_past += subband_bytes; // subband size (not include bytes for header)

			  }
#endif

#else
			  dec_band = &(decGOP[m][i].dec_subs[k*3]);
			  if(fseek(fpbit, total_bytes_past, SEEK_SET) != 0){// chen
				printf(" total_bytes_past = %d, fseek error(ezbc_dec_3d.cpp)\n", total_bytes_past);
				exit(1);
			  }
			  past = read_substream_length(&subband_bytes, fpbit);
			  total_bytes_past += past;
			  bytes_past += past;  //subband header
			  dec_band->sub_decoder->new_open_file(file_name, total_bytes_past, subband_bytes);
			  total_bytes_past += subband_bytes;  
		 	  bytes_past += subband_bytes; // subband size (not include bytes for header)
#endif

        }
      }
    }
  }

  fclose( fpbit );
  return bytes_past;
}


// the difference from encode_GOP_RD_passes is interleaving subbands of different frames
// subband_flag shows which subband needs to be skipped. This information comes from subband_mask and the process of decoding, we can not overwrite subband_mask, because it will be used to close files.
void
EzbcDec3d::decode_GOP_RD_passes2( void )
{
  int i, j, m, n, k, pass_idx, lev, msb, max_msb;
  int num_ACcoder[3];
  DecSubband *dec_band;
  DEC_SUBBAND_TREE_TYPE *decGOP[3] = { decY, decU, decV };
  int total_byte_budget = 0;
  char ***subband_flag[3], *substream_flag;
  int s; // by Yongjun Wu

  byte_used = 0;
  for( m = 0; m < ncomps; m++ ) {
    for( i = 1; i < GOPsz; i++ ) {
      if( nbands[m] != decGOP[m][i].subband_tree->get_nband(  ) ) {
        printf
          ( "can not handle this case in interleaving subbands in a GOP\n" );
        exit( 1 );
      }
    }
  }

  //allocate subband decoding flag
  NEW_VECTOR( substream_flag, total_num_ACcoder, char, "substream_flag" );

  for( m = 0; m < ncomps; m++ ) {
//	printf("nbands = %d\n",nbands[m]);
    NEW_VECTOR( subband_flag[m], GOPsz, char **,
                "subband_flag[m] (ezbc_dec_3d.cpp)" );
    for( i = 0; i < GOPsz; i++ ) {
      NEW_VECTOR( subband_flag[m][i], nbands[m], char *,
                  "subband_flag[m][i] (ezbc_dec_3d.cpp)" );
      for( k = 0; k < nbands[m]; k++ ) {
        subband_flag[m][i][k] = &substream_flag[subband_ACcoder[m][i][k]];      //several subbands share the same flag
        *subband_flag[m][i][k] = subband_mask[m][i][k];
      }
    }
  }

  // for frequency roll off by Yongjun Wu
  int roll_level1, roll_level2; 
  if (frame_res == HD ) //Added by Yuan Liu
  {
	  roll_level1 = 5; roll_level2 = 4;  // for HD sequence both CIF and QCIF need frequency roll-off 
  }else if (frame_res == SD2)
  {
	  roll_level1 = 5; roll_level2 = 4;  // for SD sequence both CIF and QCIF need frequency roll-off 
  }else if (frame_res == SD )
  {
	  roll_level1 = 4; roll_level2 = 4;  // for SD sequence both CIF and QCIF need frequency roll-off 
  }else if (frame_res == CIF )
	  roll_level1 = roll_level2 = 4;     // for CIF sequence, only QCIF needs frequency roll-off 

  //calculate total bytes of all substreams
  for( i = 0; i < GOPsz; i++ ) {
    m = 0;
    num_ACcoder[m] = decGOP[m][i].subband_tree->get_pyr_levels(  ) + 1;

    for( k = 0; k < num_ACcoder[m]; k++ ) {
#ifdef FREQUENCY_ROLL_OFF

#ifdef ROLL_STRUCTURE_ONE
	  if (k!=roll_level1  && k!=roll_level2)
	  {
		if(*subband_flag[m][i][k*3] == 0){ 
			dec_band = &(decGOP[m][i].dec_subs[k*3]);
			total_byte_budget += dec_band->sub_decoder->byte_budget;  
		}
	  }else
	  {
		  for (s=2; s>=0; s--)
		  {
			  if(*subband_flag[m][i][k*3-s] == 0){ 
				  dec_band = &(decGOP[m][i].dec_subs[k*3-s]);
				  total_byte_budget += dec_band->sub_decoder->byte_budget;  
			  }
		  }
	  }
#endif 

#ifdef ROLL_STRUCTURE_TWO
	  if (k<roll_level2)
	  {
		if(*subband_flag[m][i][k*3] == 0){ 
			dec_band = &(decGOP[m][i].dec_subs[k*3]);
			total_byte_budget += dec_band->sub_decoder->byte_budget;  
		}
	  }else if (k==roll_level2)
	  {
		  for (s=2; s>=0; s--)  // 3 independent bitstreams in this level 
		  {
			  if(*subband_flag[m][i][k*3-s] == 0){ 
				  dec_band = &(decGOP[m][i].dec_subs[k*3-s]);
				  total_byte_budget += dec_band->sub_decoder->byte_budget;  
			  }
		  }
	  }else if (k==roll_level1)
	  {
		  for (s=0; s<12; s++) // 12 independent bitstreams in this level 
		  {
			  if(*subband_flag[m][i][(k-1)*3+1+s] == 0){ 
				  dec_band = &(decGOP[m][i].dec_subs[(k-1)*3+1+s]);
				  total_byte_budget += dec_band->sub_decoder->byte_budget;  
			  }
		  }

	  }else
	  {
		if(*subband_flag[m][i][25] == 0){ 
			dec_band = &(decGOP[m][i].dec_subs[25]);
			total_byte_budget += dec_band->sub_decoder->byte_budget;  
		}
	  }
#endif 


#else
	  if(*subband_flag[m][i][k*3] == 0){ 
		  dec_band = &(decGOP[m][i].dec_subs[k*3]);
		  total_byte_budget += dec_band->sub_decoder->byte_budget;  
	  }
#endif

    }
  }


/**********
  decode 
mean values
***********/
  for( m = 0; m < ncomps; m++ ) {       // componenets
    k = 0;
    for( i = 0; i < GOPsz; i++ ) {
      if( *subband_flag[m][i][k] == 0 ) {
        GOPmean_shift[m][i] = decGOP[m][i].dec_subs[k].sub_decoder->decode_bits( 4 );   //printf("%d\n", mean_shift);
        if( decGOP[m][i].dec_subs[k].sub_decoder->decode_bits( 1 ) )
          GOPmean[m][i] =
            -( decGOP[m][i].dec_subs[k].sub_decoder->decode_bits( 10 ) );
        else
          GOPmean[m][i] =
            decGOP[m][i].dec_subs[k].sub_decoder->decode_bits( 10 );
      }
    }
  }


/**********
  decode 
subband msb
***********/
  for( m = 0; m < ncomps; m++ ) {       // componenets
    for( k = 0; k < nbands[m]; k++ ) {
      for( i = 0; i < GOPsz; i++ ) {
        if( *subband_flag[m][i][k] == 0 ) {
          decGOP[m][i].dec_subs[k].start_dec_subband(  );
        }
      }
    }
  }

  max_msb = 0;
  for( m = 0; m < ncomps; m++ ) {       // componenets
    for( i = 0; i < GOPsz; i++ ) {
      decGOP[m][i].subband_tree->reset_max_msb(  );
      msb = decGOP[m][i].subband_tree->get_max_msb(  ); //printf("msb = %d\n", msb);
      if( max_msb < msb )
        max_msb = msb;
    }
  }

  for( m = 0; m < ncomps; m++ ) {
    for( k = 0; k < nbands[m]; k++ ) {
      for( i = 0; i < GOPsz; i++ ) {
        if( *subband_flag[m][i][k] == 0 ) {
          dec_band = &( decGOP[m][i].dec_subs[k] );

          end_of_decoding( dec_band->sub_decoder->bytes_used(  ),
                           dec_band->sub_decoder->byte_budget, m, k, i,
                           subband_flag[m][i][k], &byte_used );
        }
      }
    }
  }


  max_depth = decY[GOPsz - 1].dec_subs[nbands[0] - 1].get_qtree_depth(  );

  // bitplane decoding
  for( n = max_msb + EXTRA_BIT; n >= EXTRA_BIT; n-- ) {
    /* 
       *subpass for LIP
     */
    pass_idx = 0;
    for( m = 0; m < ncomps; m++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        for( i = 0; i < GOPsz; i++ ) {

          if( *subband_flag[m][i][k] == 0 ) {

            dec_band = &( decGOP[m][i].dec_subs[k] );
            if( n < dec_band->max_bit_idx ) {

              dec_band->bit_idx = n;
              dec_band->cur_pass = pass_idx;
              dec_band->set_bit_idx_mask(  );
              dec_band->LSP_break_pt = dec_band->node_list.LSP_end + 1;
              dec_band->update_node_cxts(  );
              ( dec_band->*( dec_band->decode_LIP ) ) (  );

              end_of_decoding( dec_band->sub_decoder->bytes_used(  ),
                               dec_band->sub_decoder->byte_budget, m, k, i,
                               subband_flag[m][i][k], &byte_used );
            }
          }
        }                       // i
      }                         // k
    }                           // m

    /* 
       *subpass for LIS leaves
     */
    pass_idx++;
    for( m = 0; m < ncomps; m++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        for( i = 0; i < GOPsz; i++ ) {
          dec_band = &( decGOP[m][i].dec_subs[k] );

          if( *subband_flag[m][i][k] == 0 ) {

            if( n < dec_band->max_bit_idx ) {
              dec_band->cur_pass = pass_idx;
              ( dec_band->*( dec_band->decode_LIS_leaves ) ) (  );

              end_of_decoding( dec_band->sub_decoder->bytes_used(  ),
                               dec_band->sub_decoder->byte_budget, m, k, i,
                               subband_flag[m][i][k], &byte_used );
            }
          }
        }                       // i
      }                         // k
    }                           // m

    /* 
       *subpass for other LIS 
     */
    for( lev = 2; lev <= max_depth; lev++ ) {
      pass_idx++;
      for( m = 0; m < ncomps; m++ ) {
        for( k = 0; k < nbands[m]; k++ ) {
          for( i = 0; i < GOPsz; i++ ) {        // fr #     printf("frame= %d\n", i);

            if( *subband_flag[m][i][k] == 0 ) {

              dec_band = &( decGOP[m][i].dec_subs[k] );

              if( ( n == dec_band->max_bit_idx )
                  && ( lev == dec_band->qtree.depth ) ) {
                dec_band->bit_idx = n;
                dec_band->cur_pass = pass_idx;
                dec_band->set_bit_idx_mask(  );
                dec_band->LSP_break_pt = dec_band->node_list.LSP_end + 1;

#ifdef GET_PARENT_MODELS
                if( dec_band->par_cxt_qtree ) {

#ifdef INITIALIZE_JSIG_MODELS_FROM_PAR

                  MODEL_TYPE *par_jsig_models, *jsig_models;
                  int cxts, par_cxts;

                  cxts = par_cxts = 0;
                  for( j = dec_band->qtree.depth - 2; j >= 0; j-- ) {
                    //parent band with depth = qtree.depth - 1
                    cxts += dec_band->cxt_qtree.sig_cxts[j];
                    par_cxts += dec_band->par_cxt_qtree->sig_cxts[j];
                    assert( cxts == par_cxts );
                  }
                  jsig_models = dec_band->cxt_qtree.cxt_models +
                    dec_band->cxt_qtree.sig_offsets[0];
                  par_jsig_models = dec_band->par_cxt_qtree->cxt_models +
                    dec_band->par_cxt_qtree->sig_offsets[0];
                  for( j = cxts - 1; j >= 0; ) {
                    jsig_models[j].reset( par_jsig_models[j] );
                    jsig_models[j--].taub_scale(  );
                  }
#endif
                }               //par_cxt_qtree
#endif //GET_PARENT_MODELS


#ifdef LSP_BIT_IDX
                //empty lists in the top plane   

                for( j = dec_band->qtree.depth - 1; j >= 0; j-- )
                  dec_band->node_list.LSP_ids[n][j] =
                    dec_band->node_list.LSP_end;
#endif


                if( *subband_flag[m][i][k] == 0 ) {

                  ( dec_band->*( dec_band->decode_sig_node ) ) ( 0, lev - 1 );  //root node


                  end_of_decoding( dec_band->sub_decoder->bytes_used(  ),
                                   dec_band->sub_decoder->byte_budget, m, k,
                                   i, subband_flag[m][i][k], &byte_used );

                }

                if( *subband_flag[m][i][k] == 0 ) {

                  dec_band->decode_LIS_stack(  );

                  end_of_decoding( dec_band->sub_decoder->bytes_used(  ),
                                   dec_band->sub_decoder->byte_budget, m, k,
                                   i, subband_flag[m][i][k], &byte_used );
                }
              } else if( ( n < dec_band->max_bit_idx )
                         && ( lev < dec_band->qtree.depth ) ) {
                if( *subband_flag[m][i][k] == 0 ) {
                  dec_band->cur_pass = pass_idx;
                  dec_band->cur_lev = lev;

                  ( dec_band->*( dec_band->decode_cur_qtree_level ) ) (  );

                  end_of_decoding( dec_band->sub_decoder->bytes_used(  ),
                                   dec_band->sub_decoder->byte_budget, m, k,
                                   i, subband_flag[m][i][k], &byte_used );
                }
              }                 //else if
            }
          }                     // i
        }                       // for k      
      }                         // for m    
    }                           // for l ev

    /* 
       *subpass for LSP 
     */
    pass_idx++;
    for( m = 0; m < ncomps; m++ ) {
      for( k = nbands[m] - 1; k >= 0; k-- ) {
        for( i = 0; i < GOPsz; i++ ) {

          if( *subband_flag[m][i][k] == 0 ) {

            dec_band = &( decGOP[m][i].dec_subs[k] );
            if( n <= dec_band->max_bit_idx ) {
              dec_band->cur_pass = pass_idx;

              ( dec_band->*( dec_band->decode_LSP ) ) (  );
              end_of_decoding( dec_band->sub_decoder->bytes_used(  ),
                               dec_band->sub_decoder->byte_budget, m, k, i,
                               subband_flag[m][i][k], &byte_used );

              dec_band->reset_cxt_models(  );
            }                   //if(n <= dec_band->max_bit_idx)
          }
        }                       // i
      }                         //for k
    }                           // for m       
  }                             // for n


  if( byte_used != total_byte_budget ) {
    printf
      ( " \a byte_used %d total_byte_budget %d not all of bytes are decoded! (ezbc_dec_3d.cpp)\n",
        (int) byte_used, total_byte_budget );

    byte_used = total_byte_budget;

    // exit(1);
  }
  //close files
  for( i = 0; i < GOPsz; i++ ) {
    for( m = 0; m < 1; m++ ) {  // since Y U V are grouped together
      for( k = 0; k < num_ACcoder[m]; k++ ) {

#ifdef FREQUENCY_ROLL_OFF

#ifdef ROLL_STRUCTURE_ONE
		if (k!=roll_level1  && k!=roll_level2)
		{
			if(subband_mask[m][i][k*3] == 0){ 
				dec_band = &(decGOP[m][i].dec_subs[k*3]);
//				printf("Close option 1!\n");
				dec_band->sub_decoder->close_file();
			}
		}else
		{
			for (s=2; s>=0; s--)
			{
				if(subband_mask[m][i][k*3-s] == 0){ 
					dec_band = &(decGOP[m][i].dec_subs[k*3-s]);
//					printf("Close option 2!\n");
					dec_band->sub_decoder->close_file();
				}
			}
		}
#endif 

#ifdef ROLL_STRUCTURE_TWO
		if (k<roll_level2)
		{
			if(subband_mask[m][i][k*3] == 0){ 
				dec_band = &(decGOP[m][i].dec_subs[k*3]);
				dec_band->sub_decoder->close_file();
			}
		}else if (k==roll_level2)
		{
			for (s=2; s>=0; s--)  // 3 independent bitstreams in this level 
			{
				if(subband_mask[m][i][k*3-s] == 0){ 
					dec_band = &(decGOP[m][i].dec_subs[k*3-s]);
					dec_band->sub_decoder->close_file();
				}
			}
		}else  if (k==roll_level1)
		{
			for (s=0; s<12; s++)  // 12 independent bitstreams in this level 
			{
				if(subband_mask[m][i][(k-1)*3+1+s] == 0){ 
					dec_band = &(decGOP[m][i].dec_subs[(k-1)*3+1+s]);
					dec_band->sub_decoder->close_file();
				}
			}
		}else 
		{
			if(subband_mask[m][i][25] == 0){ 
				dec_band = &(decGOP[m][i].dec_subs[25]);
				dec_band->sub_decoder->close_file();
			}
		}
#endif

#else
	    if(subband_mask[m][i][k*3] == 0){ 
		  dec_band = &(decGOP[m][i].dec_subs[k*3]);
		  dec_band->sub_decoder->close_file();
		}
#endif

      }
    }
  }

  DELETE_VECTOR( substream_flag );
  for( m = 0; m < ncomps; m++ ) {
    for( i = 0; i < GOPsz; i++ ) {
      DELETE_VECTOR( subband_flag[m][i] );
    }
    DELETE_VECTOR( subband_flag[m] );
  }
}



void
EzbcDec3d::

end_of_decoding( int subband_bytes_used, int subband_byte_budget, int comp,
                 int band, int fr, char *subband_flag, long int *byte_used )
{
  if( subband_bytes_used > subband_byte_budget ) {

    *subband_flag = 1;

    *byte_used += subband_byte_budget;
    // printf("m %d k %d i %d used %d budget %d byte_used = %ld (ezbc_dec_3d.cpp)\n", comp, band, fr, subband_bytes_used, subband_byte_budget, *byte_used);
    /*if(byte_used >= total_byte_budget){
       return;
       } */
  }
}

void
EzbcDec3d::free_lists( void )
{
  int i, m, k;
  DEC_SUBBAND_TREE_TYPE *decGOP[3] = { decY, decU, decV };


  for( m = 0; m < ncomps; m++ ) {       // componenets
    for( i = 0; i < GOPsz; i++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        decGOP[m][i].dec_subs[k].clear_node_list(  );
        decGOP[m][i].dec_subs[k].delete_coding_stats(  );
      }
    }
  }

}
