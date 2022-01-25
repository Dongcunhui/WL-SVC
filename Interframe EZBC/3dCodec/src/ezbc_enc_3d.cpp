/* ========================================================================= */
/* Description: wavelet decomposotion structure                              */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* ========================================================================= */
#include "structN.h"
#include "ezbc_enc_3d.h"
#include <stdio.h>

float avg_mse = 1;

long int rd_count = 0;

int write_substream_length( int length, FILE * fp );

EzbcEnc3d::EzbcEnc3d( videoinfo & info, ENCODER_TYPE * enc )
:EzbcCodec3d( info ), encoder( enc )
{
  //encY = encU = encV = NULL; 


  NEW_VECTOR( GOPmean[0], GOPsz, int, "GOPmean[0]" );
  NEW_VECTOR( GOPmean_shift[0], GOPsz, int, "GOPmean_shift[0]" );

  if( cdim.x ) {
    NEW_VECTOR( GOPmean[1], GOPsz, int, "GOPmean[1]" );
    NEW_VECTOR( GOPmean_shift[1], GOPsz, int, "GOPmean_shift[1]" );

    NEW_VECTOR( GOPmean[2], GOPsz, int, "GOPmean[2]" );
    NEW_VECTOR( GOPmean_shift[2], GOPsz, int, "GOPmean_shift[2]" );
  }

  NEW_VECTOR( encY, GOPsz, ENC_SUBBAND_TREE_TYPE, "encY" );

  if( cdim.x > 0 ) {
    NEW_VECTOR( encU, GOPsz, ENC_SUBBAND_TREE_TYPE, "encU" );
    NEW_VECTOR( encV, GOPsz, ENC_SUBBAND_TREE_TYPE, "encV" );
  } else {
    encU = encV = NULL;
  }

}

EzbcEnc3d::~EzbcEnc3d( void )
{
  if( encY )
    DELETE_VECTOR( encY );
  if( encU )
    DELETE_VECTOR( encU );
  if( encV )
    DELETE_VECTOR( encV );

  DELETE_VECTOR( GOPmean[0] );
  DELETE_VECTOR( GOPmean_shift[0] );

  if( cdim.x ) {
    DELETE_VECTOR( GOPmean[1] );
    DELETE_VECTOR( GOPmean_shift[1] );

    DELETE_VECTOR( GOPmean[2] );
    DELETE_VECTOR( GOPmean_shift[2] );

  }
  //DELETE_VECTOR(byte_upto_bitplane); 062602

  if( encoder )
    delete[]encoder;

}

// 导入图像，并进行了小波变换
void
EzbcEnc3d::load_images( YUVimage *Frs, videoinfo &info, float *gop_mse, int curr )
{
  int i, GOP, band;

  float getnum;
  static int frame=-1;

  if( 0 )
  {
	  for( i = GOPsz - 1; i >= 0; i-- ) {
		  imgY[i].high_band = YES;
		  imgU[i].high_band = YES;
		  imgV[i].high_band = YES;
	  }
  }
  else{
	  for( i = GOPsz - 1; i >= 0; i-- ) {
		  imgY[i].high_band = NO;
		  imgU[i].high_band = NO;
		  imgV[i].high_band = NO;
	  }
  }
  frame++;
  GOP  = (int)(frame / ((int)(pow (2, info.tPyrLev))));
  band = (int)(frame % ((int)(pow (2, info.tPyrLev)))); 

  for( i = GOPsz - 1; i >= 0; i-- ) {
    //for(i = 0; i >= 0; i--){
    imgY[i].read_float_image( dim, Frs[i].Y ); // 从Frs读数据进imgY

    imgY[i].transform(  ); // 进行了二维小波变换（可选择是harr或其他小波），imgY是减掉均值的

	getnum = imgY[i].get_mse( );
//	printf("frame mse = %f\n",getnum);
	*gop_mse = *gop_mse * getnum;
//	printf("avg_mse = %f\n",*gop_mse);

    if( pyrY[i] )
      DELETE_OBJECT( pyrY[i] );


    //assert(pyrY[i] = new SUBBAND_TREE_TYPE(imgY[i]));
    pyrY[i] = new SUBBAND_TREE_TYPE( imgY[i], 0 );
    assert( pyrY[i] );
    imgY[i].dispose(  );

	// 下面进行U、V
    if( cdim.x ) {
      imgU[i].read_float_image( cdim, Frs[i].U );
      imgU[i].transform(  );
//        imgU[i].scale_down(); printf("scale_down in ezbc_enc_3d.cpp\n");

      if( pyrU[i] )
        DELETE_OBJECT( pyrU[i] );
      //assert(pyrU[i] = new SUBBAND_TREE_TYPE(imgU[i]));
      pyrU[i] = new SUBBAND_TREE_TYPE( imgU[i], 0 );
      assert( pyrU[i] );
      imgU[i].dispose(  );

      imgV[i].read_float_image( cdim, Frs[i].V );
      imgV[i].transform(  );

//        imgV[i].scale_down(); printf("scale_down in ezbc_enc_3d.cpp\n");
      if( pyrV[i] )
        DELETE_OBJECT( pyrV[i] );
      //assert(pyrV[i] = new SUBBAND_TREE_TYPE(imgV[i]));
      pyrV[i] = new SUBBAND_TREE_TYPE( imgV[i], 0 );
      assert( pyrV[i] );
      imgV[i].dispose(  );
    }
  } // i
}

/*
 * this function allocates AC coders to different spatial subband of 3 color components
 *
 * num_ACcoder[color_component][frame]: number of AC coders
 */ 
void
EzbcEnc3d::initialization( void )
{
  int i, j;
  ENC_SUBBAND_TREE_TYPE *encGOP[3] = { encY, encU, encV };
  SUBBAND_TREE_TYPE **GOPpyr[3] = { pyrY, pyrU, pyrV };

  EzbcCodec3d::initialize(  );


  for( i = 0; i < GOPsz; i++ ) {
    for( j = 0; j < ncomps; j++ ) {
      encGOP[j][i].SubbandTreeCodec::initialize( GOPpyr[j][i] );
    }
  }

  NEW_VECTOR( encoder, total_num_ACcoder, ENCODER_TYPE, "encoder (enc_subband_tree.cpp)" );     //062602
  for( i = 0; i < total_num_ACcoder; i++ ) {
    encoder[i].new_open_file(  );       //062602
  }

  for( i = 0; i < GOPsz; i++ ) {
    for( j = 0; j < ncomps; j++ ) {
      encGOP[j][i].initialize( GOPpyr[j][i], subband_ACcoder[j][i], encoder );
    }
  }

}



void
EzbcEnc3d::reset_GOP_enc( void )
{
  int i;

  for( i = GOPsz - 1; i >= 0; i-- ) {
    encY[i].reset_tree_enc( pyrY[i] );
    if( cdim.x > 0 ) {
      encU[i].reset_tree_enc( pyrU[i] );
      encV[i].reset_tree_enc( pyrV[i] );
    }
  }
}

void
EzbcEnc3d::GOP_transform_mean( void )
{
  int i;
  //int mean, mean_shift;

  for( i = 0; i < GOPsz; i++ ) {
    GOPmean[0][i] = imgY[i].transform_mean(  );
    GOPmean_shift[0][i] = imgY[i].mean_shift(  );

    if( cdim.x > 0 ) {
      GOPmean[1][i] = imgU[i].transform_mean(  );
      GOPmean_shift[1][i] = imgU[i].mean_shift(  );

      GOPmean[2][i] = imgV[i].transform_mean(  );
      GOPmean_shift[2][i] = imgV[i].mean_shift(  );
    }
  }
}


void
EzbcEnc3d::encode_GOP( void )
{

  encode_GOP_RD_passes2(  );
}


// the difference from encode_GOP_RD_passes is interleaving subbands of different frames
void
EzbcEnc3d::encode_GOP_RD_passes2( void )
{
  int i, j, m, n, k, pass_idx, lev, msb, max_msb, subplane_index;
  EncSubband *enc_band;         // class EncSubband defined in dwt_bitplane_enc.h
  ENC_SUBBAND_TREE_TYPE *encGOP[3] = { encY, encU, encV };
  long int total_byte_budget;

 /********
 initializ
 ********/
  //total_byte_budget = byte_budget;
  total_byte_budget = MAX_STD_INT;

  for( m = 0; m < ncomps; m++ ) {
    for( i = 1; i < GOPsz; i++ ) {
      if( nbands[m] != encGOP[m][i].subband_tree->get_nband(  ) ) {
        printf
          ( "can not handle this case in interleaving subbands in a GOP\n" );
        exit( 1 );
      }
    }
  }
  max_msb = 0;
  for( m = 0; m < ncomps; m++ ) {
    for( i = 0; i < GOPsz; i++ ) {
      msb = encGOP[m][i].subband_tree->get_max_msb(  );
      if( max_msb < msb )
        max_msb = msb;
    }
  }
  GOP_max_msb = max_msb;

  max_depth = encY[GOPsz - 1].enc_subs[nbands[0] - 1].get_qtree_depth(  );

  total_subplane = ( max_msb + 1 ) * ( 2 + max_depth );
  
   
  subplane_index = total_subplane;

  for( m = 0; m < ncomps; m++ ) {
    for( k = 0; k < nbands[m]; k++ ) {
      for( i = 0; i < GOPsz; i++ ) {
        enc_band = &( encGOP[m][i].enc_subs[k] );
        NEW_VECTOR( enc_band->sub_encoder->byte_upto_bitplane, subplane_index,
                    long int, "enc_band->sub_encoder->byte_upto_bitplane" );
      }
    }
  }

/**********
  encode 
mean values
***********/
  for( m = 0; m < ncomps; m++ ) {       // componenets  
    k = 0;
    for( i = 0; i < GOPsz; i++ ) {
      encGOP[m][i].enc_subs[k].sub_encoder->code_bits( 4, GOPmean_shift[m][i] );        //printf("%d\n", mean_shift);
      if( GOPmean[m][i] >= 0 ) {
        encGOP[m][i].enc_subs[k].sub_encoder->code_bits( 1, 0 );
        encGOP[m][i].enc_subs[k].sub_encoder->code_bits( 10, GOPmean[m][i] );
      } else {
        encGOP[m][i].enc_subs[k].sub_encoder->code_bits( 1, 1 );
        encGOP[m][i].enc_subs[k].sub_encoder->code_bits( 10, -GOPmean[m][i] );
      }
    }
  }

/**********
  encode 
subband msb
***********/
  for( m = 0; m < ncomps; m++ ) {
    for( k = 0; k < nbands[m]; k++ ) {
      for( i = 0; i < GOPsz; i++ ) {
        encGOP[m][i].enc_subs[k].start_enc_subband(  );
      }
    }
  }

  //printf("max_msb = %d Extra_bit = %d(ezbc_enc_3d.c)\n", max_msb, EXTRA_BIT); 
  subplane_index--;
  for( n = max_msb + EXTRA_BIT; n >= EXTRA_BIT; n-- ) { // EXTRA_BIT defined in Utils/subband.h   
    /* 
       *subpass for LIP
     */
    pass_idx = 0;
    for( m = 0; m < ncomps; m++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        for( i = 0; i < GOPsz; i++ ) {  // fr #     //printf("frame = %d\n", i);

          enc_band = &( encGOP[m][i].enc_subs[k] );

          if( n < enc_band->max_bit_idx ) {
            enc_band->bit_idx = n;
            enc_band->cur_pass = pass_idx;
            enc_band->set_mag_mask(  );
            enc_band->set_bit_idx_mask(  );
            enc_band->update_node_cxts(  );
            ( enc_band->*( enc_band->encode_LIP ) ) (  );       // void EncSubband::encode_LIP_cxt_AC() in dwt_bitplane_enc_cxt_AC.C
            if( encoder->bytes_used(  ) > total_byte_budget ) {
              return;
            }
          }
        }                       // i
      }                         // k
    }                           // m

    for( m = 0; m < ncomps; m++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        for( i = 0; i < GOPsz; i++ ) {
          enc_band = &( encGOP[m][i].enc_subs[k] );
          enc_band->sub_encoder->byte_upto_bitplane[subplane_index] =
            enc_band->sub_encoder->bytes_used(  );
        }
      }
    }
    subplane_index--;

    /* 
       *subpass for LIS leaves
     */
    pass_idx++;
    for( m = 0; m < ncomps; m++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        for( i = 0; i < GOPsz; i++ ) {
          enc_band = &( encGOP[m][i].enc_subs[k] );
          if( n < enc_band->max_bit_idx ) {
            enc_band->cur_pass = pass_idx;
            ( enc_band->*( enc_band->encode_LIS_leaves ) ) (  );        //void EncSubband::encode_LIS_leaves_cxt_AC() in dwt_bitplane_enc_cxt_AC.C

            if( encoder->bytes_used(  ) > total_byte_budget ) {
              return;
            }
          }
        }                       // i
      }                         // k
    }                           // m


    for( m = 0; m < ncomps; m++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        for( i = 0; i < GOPsz; i++ ) {
          enc_band = &( encGOP[m][i].enc_subs[k] );
          enc_band->sub_encoder->byte_upto_bitplane[subplane_index] =
            enc_band->sub_encoder->bytes_used(  );
        }
      }
    }
    subplane_index--;

    /* 
       *subpass for other LIS 
     */
    for( lev = 2; lev <= max_depth; lev++ ) {
      pass_idx++;
      for( m = 0; m < ncomps; m++ ) {
        for( k = 0; k < nbands[m]; k++ ) {
          for( i = 0; i < GOPsz; i++ ) {        // fr #     //printf("frame = %d\n", i);

            enc_band = &( encGOP[m][i].enc_subs[k] );

            if( ( n == enc_band->max_bit_idx )
                && ( lev == enc_band->qtree.depth ) ) {
              // the first bitplane
              enc_band->bit_idx = n;
              enc_band->cur_pass = pass_idx;
              enc_band->set_mag_mask(  );
              enc_band->set_bit_idx_mask(  );

#ifdef GET_PARENT_MODELS

              if( enc_band->par_cxt_qtree ) {


/*#ifdef INITIALIZE_SIGN_MODELS_FROM_PAR
      MODEL_TYPE *par_sign_models, *sign_models;
      sign_models =  enc_band->cxt_qtree.cxt_models +
	enc_band->cxt_qtree.sign_offset;
      par_sign_models = enc_band->par_cxt_qtree->cxt_models +
	enc_band->par_cxt_qtree->sign_offset;
      for(j = enc_band->cxt_qtree.sign_cxts - 1; j >= 0;){
	sign_models[j].reset(par_sign_models[j]);
	sign_models[j--].taub_scale();
      }
#endif*/

#ifdef INITIALIZE_JSIG_MODELS_FROM_PAR

                MODEL_TYPE *par_jsig_models, *jsig_models;
                int cxts, par_cxts;


                cxts = par_cxts = 0;
                for( j = enc_band->qtree.depth - 2; j >= 0; j-- ) {
                  //parent band with depth = qtree.depth - 1
                  cxts += enc_band->cxt_qtree.sig_cxts[j];
                  par_cxts += enc_band->par_cxt_qtree->sig_cxts[j];
                  assert( cxts == par_cxts );
                }
                jsig_models = enc_band->cxt_qtree.cxt_models +
                  enc_band->cxt_qtree.sig_offsets[0];
                par_jsig_models = enc_band->par_cxt_qtree->cxt_models +
                  enc_band->par_cxt_qtree->sig_offsets[0];
                for( j = cxts - 1; j >= 0; ) {
                  jsig_models[j].reset( par_jsig_models[j] );
                  jsig_models[j--].taub_scale(  );
                }

#endif

              }                 //par_cxt_qtree
#endif //GET_PARENT_MODELS

#ifdef LSP_BIT_IDX
              //empty lists in the top plane
              for( j = enc_band->qtree.depth - 1; j >= 0; j-- )
                enc_band->node_list.LSP_ids[n][j] =
                  enc_band->node_list.LSP_end;
#endif

              ( enc_band->*( enc_band->encode_sig_node ) ) ( 0, lev - 1 );      //root node
              enc_band->encode_LIS_stack(  );

              enc_band->coding_stats.bitplanes[n].passes[pass_idx].
                cumulative_bytes = enc_band->sub_encoder->bytes_used(  );

              if( encoder->bytes_used(  ) > total_byte_budget ) {
                return;
              }
            } else if( ( n < enc_band->max_bit_idx )
                       && ( lev < enc_band->qtree.depth ) ) {
              enc_band->cur_pass = pass_idx;
              enc_band->cur_lev = lev;
              ( enc_band->*( enc_band->encode_cur_qtree_level ) ) (  );

              if( encoder->bytes_used(  ) > total_byte_budget ) {
                return;
              }
            }
          }                     // for i
        }                       // for k
      }                         // for m

      for( m = 0; m < ncomps; m++ ) {
        for( k = 0; k < nbands[m]; k++ ) {
          for( i = 0; i < GOPsz; i++ ) {
            enc_band = &( encGOP[m][i].enc_subs[k] );
            enc_band->sub_encoder->byte_upto_bitplane[subplane_index] =
              enc_band->sub_encoder->bytes_used(  );
          }
        }
      }
      subplane_index--;
    }                           //lev

    /* 
       *subpass for LSP 
     */
    pass_idx++;
    for( m = 0; m < ncomps; m++ ) {
      for( k = nbands[m] - 1; k >= 0; k-- ) {
        for( i = 0; i < GOPsz; i++ ) {  // fr #     //printf("frame = %d\n", i);

          enc_band = &( encGOP[m][i].enc_subs[k] );
          if( n <= enc_band->max_bit_idx ) {
            enc_band->cur_pass = pass_idx;
            ( enc_band->*( enc_band->encode_LSP ) ) (  );       //void EncSubband::encode_LSP_cxt_AC()  in dwt_bitplane_enc_cxt_AC.C
            if( encoder->bytes_used(  ) > total_byte_budget ) {
              return;
            }
            enc_band->reset_cxt_models(  );
          }                     //if(n <= enc_band->max_bit_idx)
        }                       // for i
      }                         // for k
    }                           //for m

    for( m = 0; m < ncomps; m++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        for( i = 0; i < GOPsz; i++ ) {
          enc_band = &( encGOP[m][i].enc_subs[k] );
          enc_band->sub_encoder->byte_upto_bitplane[subplane_index] =
            enc_band->sub_encoder->bytes_used(  );
        }
      }
    }
    subplane_index--;

  }                             // for n


  /*
   * close files
   */
  // for frequency roll off by Yongjun Wu
  int roll_level1, roll_level2; 
  if (frame_res == HD )
  {
	  roll_level1 = 5; roll_level2 = 4; 
  }else if (frame_res == SD2)
  {
	  roll_level1 = 5; roll_level2 = 4; 
  }else if (frame_res == SD )
  {
	  roll_level1 = 4; roll_level2 = 4; 
  }else if (frame_res == CIF)
	  roll_level1 = roll_level2 = 4;

  m = 0;
  for( i = 0; i < GOPsz; i++ ) {
    for( int k = 0; k < num_Slev[m]; k++ ) {

#ifdef FREQUENCY_ROLL_OFF  // by Yongjun Wu 
		int s;

#ifdef  ROLL_STRUCTURE_ONE
	   if (k!=roll_level1  && k!=roll_level2)
	   {   // the three subbands are coded by one arithematic coder
		   enc_band = &(encGOP[m][i].enc_subs[k*3]);   
		   enc_band->sub_encoder->close_file();
	   }else
	   {
		   // the bitstream for these subbands are separate in order to do bitstream shift
		   for (s=0; s<3; s++)
		   {
			   enc_band = &(encGOP[m][i].enc_subs[k*3-s]);
			   enc_band->sub_encoder->close_file();
		   }
	   }
#endif 

#ifdef ROLL_STRUCTURE_TWO  // only for SD sequence 
	   if (k<roll_level2)
	   {   // the three subbands are coded by one arithematic coder
		   enc_band = &(encGOP[m][i].enc_subs[k*3]);
		   enc_band->sub_encoder->close_file();
	   }
	   if (k==roll_level2)  // in this level 3 subbands are coded separately 
		   for (s=0; s<3; s++)
		   {
			   enc_band = &(encGOP[m][i].enc_subs[(k-1)*3+1+s]);
			   enc_band->sub_encoder->close_file();
		   }
	   if (k==roll_level1)  // in this level 12 separate subband 
		   for (s=0; s<12; s++)
		   {
			   enc_band = &(encGOP[m][i].enc_subs[(k-1)*3+1+s]);
			   enc_band->sub_encoder->close_file();
		   }
	   // in this level the 3 subbands are coded by one arithematic coder
	   if (k==6)
	   {
		   enc_band = &(encGOP[m][i].enc_subs[25]);
		   enc_band->sub_encoder->close_file();
	   }
#endif 	   


#else	   
	   enc_band = &(encGOP[m][i].enc_subs[k*3]);
	   enc_band->sub_encoder->close_file();
#endif 
	
	}
  }
}


// 写文件
long int
EzbcEnc3d::write_files( long int GOP_counter, char *bitname )
{
  int i, m, k, bp;              //, **subband_flag
  long int bytes, total_bytes = 0;
  unsigned char *buffer;
  char name[256], command[256], strtmp[256];
  ENC_SUBBAND_TREE_TYPE *encGOP[3] = { encY, encU, encV };
  EncSubband *enc_band;         // class EncSubband defined in dwt_bitplane_enc.h
  FILE *fp, *bitstream, *substream;
  int s; // by Yongjun Wu  for frequency roll-off
#ifdef ROLL_STRUCTURE_TWO
  int small_bands; 
#endif 

  if( ( bitstream = fopen( bitname, "a+b" ) ) == NULL ) {
    printf( "can not open file %s in write_rd\n", bitname );

    exit( 1 );
  }


  /*
   * concantenate subbitstreams generated by different AC coders连接由不同AC编码器生成的子比特流
   */
  int roll_level1, roll_level2; 
  if (frame_res == HD )
  {
	  roll_level1 = 5; roll_level2 = 4; 
  }else if (frame_res == SD2)
  {
	  roll_level1 = 5; roll_level2 = 4; 
  }else if (frame_res == SD )
  {
	  roll_level1 = 4; roll_level2 = 4; 
  }else if (frame_res == CIF )
	  roll_level1 = roll_level2 = 4;

  m = 0;
  for( i = 0; i < GOPsz; i++ ) {

    for( k = 0; k < num_Slev[m]; k++ ) {

#ifdef 	FREQUENCY_ROLL_OFF // by Yongjun Wu: Structure one frequency roll off

#ifdef ROLL_STRUCTURE_ONE
	  if (k!=roll_level1  && k!=roll_level2)
	  {
		  enc_band = &(encGOP[m][i].enc_subs[k*3]);
		  bytes = enc_band->sub_encoder->bytes_used();
		  NEW_VECTOR(buffer, bytes, unsigned char, "buffer (ezbc_enc_3d.cpp)");
		  if((substream = fopen(enc_band->sub_encoder->temp_name, "rb")) == NULL){	  
			printf("can not open file %s\n", enc_band->sub_encoder->temp_name);	  
			exit(1);
		  }
		  if(fread(buffer, sizeof(unsigned char), bytes, substream) != (unsigned)bytes){
			 printf("read error (ezbc_enc_3d.cpp)\n");
			 exit(1);
		  }
		  fclose(substream);
		  total_bytes += write_substream_length(bytes, bitstream);
		  //write subband interior
		  if(fwrite(buffer, sizeof(unsigned char), bytes, bitstream) != (unsigned)bytes){
	  		printf("write error (ezbc_enc_3d.cpp)\n");
			exit(1);
		  }
		  total_bytes += bytes;
		  DELETE_VECTOR(buffer);
	  }else
	  {
		  for (s=2;s>=0; s--)
		  {   // the 3 subbands are coded by different arithematic coders
			  enc_band = &(encGOP[m][i].enc_subs[k*3-s]);
			  bytes = enc_band->sub_encoder->bytes_used();
			  NEW_VECTOR(buffer, bytes, unsigned char, "buffer (ezbc_enc_3d.cpp)");
			  if((substream = fopen(enc_band->sub_encoder->temp_name, "rb")) == NULL){	  
				  printf("can not open file %s\n", enc_band->sub_encoder->temp_name);	  
				  exit(1);
			  }
			  if(fread(buffer, sizeof(unsigned char), bytes, substream) != (unsigned)bytes){
				  printf("read error (ezbc_enc_3d.cpp)\n");
				  exit(1);
			  }
			  fclose(substream);
			  total_bytes += write_substream_length(bytes, bitstream);
			  //write subband interior
			  if(fwrite(buffer, sizeof(unsigned char), bytes, bitstream) != (unsigned)bytes){
				  printf("write error (ezbc_enc_3d.cpp)\n");
				  exit(1);
			  }
			  total_bytes += bytes;
			  DELETE_VECTOR(buffer);
		  }
	  }
#endif 

#ifdef 	ROLL_STRUCTURE_TWO
	  if (k<roll_level2)
	  {
		  enc_band = &(encGOP[m][i].enc_subs[k*3]);
		  bytes = enc_band->sub_encoder->bytes_used();
		  NEW_VECTOR(buffer, bytes, unsigned char, "buffer (ezbc_enc_3d.cpp)");
		  if((substream = fopen(enc_band->sub_encoder->temp_name, "rb")) == NULL){	  
			printf("can not open file %s\n", enc_band->sub_encoder->temp_name);	  
			exit(1);
		  }
		  if(fread(buffer, sizeof(unsigned char), bytes, substream) != (unsigned)bytes){
			 printf("read error (ezbc_enc_3d.cpp)\n");
			 exit(1);
		  }
		  fclose(substream);
		  total_bytes += write_substream_length(bytes, bitstream);
		  //write subband interior
		  if(fwrite(buffer, sizeof(unsigned char), bytes, bitstream) != (unsigned)bytes){
	  		printf("write error (ezbc_enc_3d.cpp)\n");
			exit(1);
		  }
		  total_bytes += bytes;
		  DELETE_VECTOR(buffer);
	  }else if (k==roll_level2  || k==roll_level1)
	  {
		  if (k==roll_level1) small_bands = 12;  // 12 small subbands 
		  if (k==roll_level2) small_bands = 3;   // 3  small subbands 
		  for (s=0;s<small_bands; s++)
		  {
			  enc_band = &(encGOP[m][i].enc_subs[(k-1)*3+1+s]);
			  bytes = enc_band->sub_encoder->bytes_used();
			  NEW_VECTOR(buffer, bytes, unsigned char, "buffer (ezbc_enc_3d.cpp)");
			  if((substream = fopen(enc_band->sub_encoder->temp_name, "rb")) == NULL){	  
				  printf("can not open file %s\n", enc_band->sub_encoder->temp_name);	  
				  exit(1);
			  }
			  if(fread(buffer, sizeof(unsigned char), bytes, substream) != (unsigned)bytes){
				  printf("read error (ezbc_enc_3d.cpp)\n");
				  exit(1);
			  }
			  fclose(substream);
			  total_bytes += write_substream_length(bytes, bitstream);
			  //write subband interior
			  if(fwrite(buffer, sizeof(unsigned char), bytes, bitstream) != (unsigned)bytes){
				  printf("write error (ezbc_enc_3d.cpp)\n");
				  exit(1);
			  }
			  total_bytes += bytes;
			  DELETE_VECTOR(buffer);
		  }
	  } else
	  {
		  enc_band = &(encGOP[m][i].enc_subs[25]);
		  bytes = enc_band->sub_encoder->bytes_used();
		  NEW_VECTOR(buffer, bytes, unsigned char, "buffer (ezbc_enc_3d.cpp)");
		  if((substream = fopen(enc_band->sub_encoder->temp_name, "rb")) == NULL){	  
			printf("can not open file %s\n", enc_band->sub_encoder->temp_name);	  
			exit(1);
		  }
		  if(fread(buffer, sizeof(unsigned char), bytes, substream) != (unsigned)bytes){
			 printf("read error (ezbc_enc_3d.cpp)\n");
			 exit(1);
		  }
		  fclose(substream);
		  total_bytes += write_substream_length(bytes, bitstream);
		  //write subband interior
		  if(fwrite(buffer, sizeof(unsigned char), bytes, bitstream) != (unsigned)bytes){
	  		printf("write error (ezbc_enc_3d.cpp)\n");
			exit(1);
		  }
		  total_bytes += bytes;
		  DELETE_VECTOR(buffer);
	  }
#endif		


#else		
	  enc_band = &(encGOP[m][i].enc_subs[k*3]);
	  bytes = enc_band->sub_encoder->bytes_used();
	  NEW_VECTOR(buffer, bytes, unsigned char, "buffer (ezbc_enc_3d.cpp)");
	  if((substream = fopen(enc_band->sub_encoder->temp_name, "rb")) == NULL){	  
		printf("can not open file %s\n", enc_band->sub_encoder->temp_name);	  
		exit(1);
	  }
	  if(fread(buffer, sizeof(unsigned char), bytes, substream) != (unsigned)bytes){
		printf("read error (ezbc_enc_3d.cpp)\n");
		exit(1);
	  }
	  fclose(substream);
	  total_bytes += write_substream_length(bytes, bitstream);
	  //write subband interior
	  if(fwrite(buffer, sizeof(unsigned char), bytes, bitstream) != (unsigned)bytes){
	  	printf("write error (ezbc_enc_3d.cpp)\n");
	    exit(1);
	  }
	  total_bytes += bytes;
	  DELETE_VECTOR(buffer);
#endif

    }
  }

  for( i = 0; i < Encoder::object_count; i++ ) {
    // Windows
    // sprintf( command, "del sub%03d.bit", i );
    // system( command );
    
    // Linux - Hanke, 16.09.02
    //    sprintf( command, "rm sub%03d.bit", i );
    //    system( command );

    sprintf( command, "sub%03d.bit", i );
    remove(command);
  }
  fclose( bitstream );


  /*
   * output subbitstream allocation table 输出子比特流分配表
   */
  strncpy( strtmp, bitname, strlen( bitname ) - 4 );  // Hanke, 16.09.02
  strtmp[strlen( bitname ) - 4] = '\0';               //
  sprintf( name, "%s.rd_sample_dat", strtmp );        //
  if( ( fp = fopen( name, "a+b" ) ) == NULL ) {
    printf( "can not open file in write_rd\n" );
    exit( 1 );
  }

  fprintf( fp, "%d\n", num_Slev[0] );
  rd_count ++;

  fprintf(fp, "%d\n", GOP_max_msb); // by Yongjun Wu in order to do frequency roll-off
  rd_count ++;
  for( i = 0; i < GOPsz; i++ ) {

    fprintf( fp, "%d\n", total_subplane - 1 );
	rd_count ++;

    for( bp = total_subplane - 1; bp >= 0; bp-- ) {
      m = 0;

      for( k = 0; k < num_Slev[m]; k++ ) {

#ifdef  FREQUENCY_ROLL_OFF  // Structure one frequency roll-off 

#ifdef  ROLL_STRUCTURE_ONE // 写入
		  if (k!=roll_level1  &&  k!=roll_level2)
		  {

			  enc_band = &(encGOP[m][i].enc_subs[k*3]);
			  if(bp==total_subplane-1)
				  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp], 
												enc_band->sub_encoder->byte_upto_bitplane[bp]);
	 		  else
				  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp]-
												enc_band->sub_encoder->byte_upto_bitplane[bp+1], 
												enc_band->sub_encoder->byte_upto_bitplane[bp]);

			  rd_count += 3;
		  }
		  else
		  {
			  for (s=2;s>=0; s--) // 3次
			  {
				  enc_band = &(encGOP[m][i].enc_subs[k*3-s]);
				  if(bp==total_subplane-1)
					  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp], 
													enc_band->sub_encoder->byte_upto_bitplane[bp]);
	 			  else
					  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp]-
										            enc_band->sub_encoder->byte_upto_bitplane[bp+1], 
													enc_band->sub_encoder->byte_upto_bitplane[bp]);

				  rd_count += 3;
			  }

		  }
#endif 

#ifdef  ROLL_STRUCTURE_TWO
		  if (k<roll_level2)
		  {
			  enc_band = &(encGOP[m][i].enc_subs[k*3]);
			  if(bp==total_subplane-1)
				  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp], 
												enc_band->sub_encoder->byte_upto_bitplane[bp]);
	 		  else
				  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp]-
												enc_band->sub_encoder->byte_upto_bitplane[bp+1], 
												enc_band->sub_encoder->byte_upto_bitplane[bp]);
		  }else if (k==roll_level2  || k==roll_level1)
		  {
			  if (k==roll_level1) small_bands = 12;  // 12 small subbands
			  if (k==roll_level2) small_bands = 3;   // 3  small subbands 
			  for (s=0;s<small_bands; s++)
			  {
				  enc_band = &(encGOP[m][i].enc_subs[(k-1)*3+1+s]);
				  if(bp==total_subplane-1)
					  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp], 
													enc_band->sub_encoder->byte_upto_bitplane[bp]);
	 			  else
					  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp]-
										            enc_band->sub_encoder->byte_upto_bitplane[bp+1], 
													enc_band->sub_encoder->byte_upto_bitplane[bp]);
			  }

		  }else
		  {
			  enc_band = &(encGOP[m][i].enc_subs[25]);
			  if(bp==total_subplane-1)
				  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp], 
												enc_band->sub_encoder->byte_upto_bitplane[bp]);
	 		  else
				  fprintf(fp, "%d %d %d\n", bp, enc_band->sub_encoder->byte_upto_bitplane[bp]-
												enc_band->sub_encoder->byte_upto_bitplane[bp+1], 
												enc_band->sub_encoder->byte_upto_bitplane[bp]);
		  }
#endif 


#else
		  enc_band = &(encGOP[m][i].enc_subs[k*3]);
		  if(bp==total_subplane-1)
		    fprintf(fp, "%d %d %d\n", bp, (int)enc_band->sub_encoder->byte_upto_bitplane[bp], 
										  (int)enc_band->sub_encoder->byte_upto_bitplane[bp]);
		  else
		    fprintf(fp, "%d %d %d\n", bp, (int)(enc_band->sub_encoder->byte_upto_bitplane[bp]-
			                              enc_band->sub_encoder->byte_upto_bitplane[bp+1]), 
										  enc_band->sub_encoder->byte_upto_bitplane[bp]);
#endif

      }
    }
  }

  fclose( fp );
  return ( total_bytes );
}






void
EzbcEnc3d::free_lists( void )
{
  int i, m, k;
  ENC_SUBBAND_TREE_TYPE *encGOP[3] = { encY, encU, encV };


  for( m = 0; m < ncomps; m++ ) {       // componenets
    for( i = 0; i < GOPsz; i++ ) {
      for( k = 0; k < nbands[m]; k++ ) {
        encGOP[m][i].enc_subs[k].clear_node_list(  );
        encGOP[m][i].enc_subs[k].delete_coding_stats(  );
      }
    }
  }
}
