#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "structN.h"

#include "zlib.h"
#include "dpx.h"

#define EXTERN  extern
#include "coderN.h"
#include "miscN.h"
#include "ioN.h"
#include "pstatN.h"
#include "basic.h"

#define INVALID_AGP_EXIST  10
#define INVALID_LAYER_EXIST  10
#define INVALID_BI_MV      -1
#define MV_PERCENT_THERSH  0.55  // the threshold for motion vector percentage 
                                 // above the threshold we discard some sub-symbols

int read_number_core( FILE * fp );

unsigned int read_number_sub( FILE * fp );

void write_number_sub(unsigned int outputbyte, FILE * fp );

int GOP_motion_bytes;


void write_number_core( int outputbyte, FILE * fp );
int read_GOPheader( enum FLAG **dec_scene_change, videoinfo info );
int read_substream_length( int *length, FILE * fp );
int write_substream_length( int length, FILE * fp );
int bit_alloc_VBR( videoinfo info, unsigned char *qp, 
                   long int sum_mv, long int overhead );
int bit_alloc_CBR( videoinfo info, unsigned char *qp,
				  long int overhead, unsigned long *gop_sum_mv, int curr_last );
int bit_alloc_VBR_v2( videoinfo info, unsigned char *qp, long int overhead,
					 unsigned long *gop_sum_mv, int curr_last, int vbr_bound );

int simulcast_pull( char *pre_stream, char *bit_alloc_name, 
                    int t_level, int s_level, int kbps, 
                    int qp );

// process major section for motion vector 
void process_mv_section(long int mvBytes, FILE *outfp, long int *out_bytes, long int *total_mv_bytes, 
						long int *in_bytes, FILE *fpbit, FILE *fptab, enum FLAG pulled_bitstream)
{
	unsigned char *data;

	write_number_core( mvBytes, outfp );
	*out_bytes      += 4;
	*total_mv_bytes += 4;
	data = ( unsigned char * )getarray( mvBytes, sizeof( unsigned char ), "data" );
	if( fread( data, sizeof( unsigned char ), mvBytes, fpbit ) != ( unsigned )mvBytes ) {
		printf( "fread error2 %d\n", mvBytes );
		exit( 1 );
	}
	*in_bytes += mvBytes;
	fwrite( data, sizeof( unsigned char ), mvBytes, outfp );
	*out_bytes      += mvBytes;
	*total_mv_bytes += mvBytes;
	if ( pulled_bitstream == NO ){
		fprintf( fptab, "%d \t %% motion_bytes\n", mvBytes + 4 ); 
		GOP_motion_bytes += mvBytes;
	}
	free( data );
}

// get one section of bits for AGP 
void  get_bytes_section(FILE *fpbit, FILE *outfp, FILE * fptab, int write_sign, long int *in_bytes, long int *out_bytes,
					  long int *total_mv_bytes, enum FLAG pulled_bitstream )
{
	unsigned int section_Bytes;
	unsigned char *data;

	section_Bytes = read_number_sub( fpbit ); // fpbit will move forward by 2 bytes 
	*in_bytes += 2; 
	if ( write_sign )
	{
		write_number_sub( section_Bytes, outfp );
		*out_bytes  += 2;
		*total_mv_bytes +=2; 
		data = ( unsigned char * )getarray( section_Bytes, sizeof( unsigned char ), "data" );
		if( fread( data, sizeof( unsigned char ), section_Bytes, fpbit ) != ( unsigned )section_Bytes ) {
			printf( "fread error2 %d\n", section_Bytes );
			exit( 1 );
		}
		*in_bytes += section_Bytes;
		fwrite( data, sizeof( unsigned char ), section_Bytes, outfp );
		*out_bytes      += section_Bytes;
		*total_mv_bytes += section_Bytes;
		if ( pulled_bitstream == NO ){
			fprintf( fptab, "%d \t %% motion_bytes\n", section_Bytes + 2 ); 
			GOP_motion_bytes += section_Bytes;
		}
		free( data );
	} else {
		fseek( fpbit, section_Bytes, SEEK_CUR );
		in_bytes += section_Bytes;
	}
}

// get the section for enhancement layer 
void  get_bytes_section_enhance(FILE *fpbit, FILE *outfp, FILE * fptab, int write_sign, long int *in_bytes, long int *out_bytes,
					  long int *total_mv_bytes, enum FLAG pulled_bitstream )
{
	unsigned int enhance_Bytes;
	unsigned char *data;

	enhance_Bytes = read_number_core( fpbit ); // fpbit will move forward by 2 bytes 
	*in_bytes += 4; 
	if ( write_sign )
	{
		write_number_core( enhance_Bytes, outfp );
		*out_bytes  += 4;
		*total_mv_bytes +=4; 
		data = ( unsigned char * )getarray( enhance_Bytes, sizeof( unsigned char ), "data" );
		if( fread( data, sizeof( unsigned char ), enhance_Bytes, fpbit ) != ( unsigned )enhance_Bytes ) {
			printf( "fread error2 %d\n", enhance_Bytes );
			exit( 1 );
		}
		*in_bytes += enhance_Bytes;
		fwrite( data, sizeof( unsigned char ), enhance_Bytes, outfp );
		*out_bytes      += enhance_Bytes;
		*total_mv_bytes += enhance_Bytes;
		if ( pulled_bitstream == NO ){
			fprintf( fptab, "%d \t %% motion_bytes\n", enhance_Bytes + 4 ); 
			GOP_motion_bytes += enhance_Bytes;
		}
		free( data );
	} else {
		fseek( fpbit, enhance_Bytes, SEEK_CUR );
		in_bytes += enhance_Bytes;
	}
}


void test( videoinfo info, unsigned char *qp, int estimated_overhead, unsigned long int *gop_mv, int curr_last, int type ){

  int i, j, k, nm, curr, last, num_of_GOP, GOPheader_bytes, Per_GOPheader_bytes, FAT_bytes;
  int count = 0, remaining_frs, num_subband, bytes, mvBytes;
  int kbps = 0, tmp, GOPsz, eff_GOPsz, is_pullable;
  int *subband_fat, *newsubband_fat;
  long int in_bytes = 0, out_bytes = 0, FAT_start = 0, overhead_start;
  long int overhead = 0, total_mv_bytes = 0;
  long int subband_bytes = 0, total_subband_bytes = 0, *newfat; 
  short int t_level, s_level, rel_s_level, target_rate;
  int counter, buffer;

  int getnum;
  enum FLAG CBR = YES;
  
  unsigned char *data;
  char pre_stream[256], post_stream[256], overhead_stream[256];
  char seq_name[256], bit_alloc_name[256], bit_alloc_table[256];

  unsigned char list_tlevel[25], list_slevel[25], list_qp[25], list_simulcast[25];
  unsigned short int list_rates[25], simul_bitstream = 0;
  unsigned long int sum_mv, exist_sum_mv[10], list_out_bytes[25], simul_bitstream_pos = 0;
                              
  enum FLAG pulled_bitstream; 
  FILE *outfp, *fpbytes, *fptab, *fphead;
  int   AGP_sign, sub_bit, AGP_exist_pull[MAX_TLEVELS], bi_exist[MAX_TLEVELS]; 
  long int major_mvBytes; 
  enum FLAG scene_change_left, scene_change_right; 
  int   layer_sign, layer_exist_pull[MAX_TLEVELS]; 

  int GOP_subband_bytes;

  pulled_bitstream = NO;

  in_bytes += sizeof( videoheader );
  num_of_GOP = (info.last - info.start + info.GOPsz) / info.GOPsz;
//  printf("num_of_GOP = %d\n",num_of_GOP);

  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 );
  seq_name[strlen( info.bitname ) - 4] = '\0';

  // GOPheader_bytes (cmp. mvcoding.c)
  GOPsz = info.GOPsz;
  for( i = 0; i < info.tPyrLev; i++ ) {
    for( j = 0; j < GOPsz; j++ ) 
      count++;
    GOPsz /= 2;  
  }
  GOPsz = info.GOPsz; 
  GOPheader_bytes = (( count + 7 ) / 8) * num_of_GOP; 

  Per_GOPheader_bytes = GOPheader_bytes / num_of_GOP; 
  
  // FAT bytes
  FAT_bytes = num_of_GOP * sizeof( long int ); // DENOISE ? -> only larger values

  if(type == 0)
	is_pullable = bit_alloc_CBR(info, qp, estimated_overhead, gop_mv, curr_last);
  else{
	assert(type == 1);

	sum_mv = 0;
	counter = 0;
	curr = 0;
    last = curr_last - info.start;

	while( curr <= last ){
		sum_mv += gop_mv[counter];
		printf("gop_counter = %d, gop_sum_mv = %d\n",counter,gop_mv[counter]);
		counter ++;
		curr += info.GOPsz;
	}

	buffer = info.last;
	info.last = curr_last;
	is_pullable = bit_alloc_VBR_v2(info, qp, estimated_overhead, gop_mv, curr_last, curr_last - info.GOPsz);
	info.last = buffer;
  }

  sprintf( bit_alloc_table, "%s_%d.alloc_table", seq_name, info.bitrate );
	if( !( fptab = fopen( bit_alloc_table, "wb" ) ) ) {
      printf( "can not open: %s\n", bit_alloc_table );
      exit( 1 );
    }

  if( !is_pullable ){
	printf("bitstream is not pullable!\n");
	exit(1);
  }
//
  if( !( fpbit = fopen( info.bitname, "rb" ) ) ) {
    printf( "can not open: %s\n", info.bitname );
    exit( 1 );
  }

  sprintf( post_stream, "%s_%d.bit", seq_name, info.bitrate ); 
  
  write_header( post_stream, info );
  out_bytes += sizeof( videoheader );
  
  if( !( outfp = fopen( post_stream, "r+b" ) ) ) {
    printf( "can not open: %s\n", post_stream );
    exit( 1 );
  }

  if ( pulled_bitstream == NO ){  
    fprintf( fptab, "%%---------------------------------------------------------------\n" );
    fprintf( fptab, "%s %% sequence_name\n", seq_name );
    fprintf( fptab, "%03d \t %% sequence_length\n", info.last - info.start + 1 );
    fprintf( fptab, "%03d \t %% start\n", info.start );
    fprintf( fptab, "%03d \t %% last\n", info.last  );
    fprintf( fptab, "%d   \t %% bit_rate\n", info.bitrate );
    fprintf( fptab, "%0.2f\t %% frame_rate\n",
             ((double)(info.framerate) / pow((float)2, info.t_level)));
    fprintf( fptab, "%d \t %% t_level\n", info.t_level );
    fprintf( fptab, "%d \t %% Y_size_x\n", info.ywidth  >> rel_s_level);
    fprintf( fptab, "%d \t %% Y_size_y\n", info.yheight >> rel_s_level);
    fprintf( fptab, "%d \t %% C_size_x\n", info.cwidth  >> rel_s_level);
    fprintf( fptab, "%d \t %% C_size_y\n", info.cheight >> rel_s_level);
    fprintf( fptab, "%d \t %% s_level\n", info.s_level );
    fprintf( fptab, "%d \t %% number_of_GOP\n", num_of_GOP );
    fprintf( fptab, "%d \t %% header/fat_size\n",
             sizeof( videoheader ) + FAT_bytes + GOPheader_bytes );
    fprintf( fptab, "%ld \t %% overhead_size\n", estimated_overhead
          - (sizeof( videoheader ) - FAT_bytes - GOPheader_bytes) );
    fprintf( fptab, "---------------------------------------------------------------\n" );
  }

    
  /*******/
  /* FAT */
  /*******/
  newfat = ( long int * )getarray( num_of_GOP, sizeof( long int ), "newfat" );
  
  //skip FAT
  in_bytes += num_of_GOP * sizeof( long int );
  fseek( fpbit, in_bytes, SEEK_SET );
  
  fseek( outfp, out_bytes, SEEK_SET );
  fwrite( newfat, sizeof( long int ), num_of_GOP, outfp );
  FAT_start = out_bytes;
  out_bytes += num_of_GOP * sizeof( long int );
  overhead  = out_bytes;
  
  simul_scene_change =
    ( enum FLAG ** )getarray( info.tPyrLev, sizeof( enum FLAG * ),
                              "simul_scene_change" );
  GOPsz = info.GOPsz;
  for( i = 0; i < info.tPyrLev; i++ ) {
    simul_scene_change[i] =
      ( enum FLAG * )getarray( GOPsz + 1, sizeof( enum FLAG ), "simul_scene_change" );
    GOPsz /= 2;
  }
  
  curr = info.start;
  last = curr_last;

  int GOPcounter = 0;

  sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, info.bitrate );
  if( CBR == NO || CBR == YES ) {
    if( !( fpbytes = fopen( bit_alloc_name, "rb" ) ) ) {
      printf( "can not open %s\n", bit_alloc_name );
      exit( 1 );
    }
  }

  fscanf( fpbytes, "%d\n", &num_subband );

  subband_fat =
    ( int * )getarray( num_subband, sizeof( int ), "subband_fat" );
  newsubband_fat =
    ( int * )getarray( num_subband, sizeof( int ), "newsubband_fat" );

  while( curr <= last ) {       /* MCTF decoding */
    
    remaining_frs = last - curr + 1;  //printf("remaining frames = %d\n", remaining_frs);                                  
     
    if ( remaining_frs < info.GOPsz ){
      // Level_change = YES;
      info.eff_GOPsz = remaining_frs;	
    } else { 
      // Level_change = NO;         // full GOP   
      info.eff_GOPsz = info.GOPsz;  // effective GOP size
    }    
      
//    printf( " pulling GOP %d: frame %d - %d\n", GOPcounter,
//            curr, curr + info.eff_GOPsz - 1 );
    if ( pulled_bitstream == NO ){
      fprintf( fptab, "%% GOP %2d %% frame %d - %d\n", GOPcounter, 
               curr, curr + info.eff_GOPsz - 1 );
    }

    GOPheader_bytes = 0;
    if( info.tPyrLev >= 1 ) {

      /*************/
      /* GOPheader */
      /*************/

      GOPheader_bytes = read_GOPheader( simul_scene_change, info );   //printf("GOPheader_bytes = %d\n", GOPheader_bytes);
      fseek( fpbit, -GOPheader_bytes, SEEK_CUR );
      
      data =
        ( unsigned char * )getarray( GOPheader_bytes, sizeof( unsigned char ),
                                     "data" );
      if( fread( data, sizeof( unsigned char ), GOPheader_bytes, fpbit ) !=
          ( unsigned )GOPheader_bytes ) {
        printf( "fread error1 \n" );
        exit( 1 );
      }
      in_bytes += GOPheader_bytes;
      
      fwrite( data, sizeof( unsigned char ), GOPheader_bytes, outfp );
      out_bytes += GOPheader_bytes;
      overhead  += GOPheader_bytes;
      if ( pulled_bitstream == NO ){
        fprintf( fptab, "%d \t %% GOPheader_bytes\n", GOPheader_bytes );
      }
      free( data );
      
      /******/
      /* MV */
      /******/

	  AGP_sign = 0;
	  layer_sign = 0;

      GOP_motion_bytes = 0; //Added on 10.06.2017

      GOPsz = 2;
      for( j = info.tPyrLev - 1; j >= 0; j-- ) {

		  if (info.bi_exist[j]) // bi-directional motion field 
		  {
			  for( k = 0; k <= GOPsz; k++ ) {
				if( simul_scene_change[j][k] == NO ) { 

					// scalable motion vector coding: AGP and layer structure 

                	if ( !AGP_sign && !layer_sign )  // no AGP, no layer structure 
					{
						mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
						in_bytes += 4;

//						printf("mvBytes = %d\n",mvBytes);

//						printf("mvBytes = %d\n",mvBytes);
						if( j >= info.t_level ) { // temporal scalability
						  process_mv_section(mvBytes, outfp, &out_bytes, &total_mv_bytes, &in_bytes, fpbit, 
							  fptab, pulled_bitstream);
						} else {
						  fseek( fpbit, mvBytes, SEEK_CUR );
						  in_bytes += mvBytes;
						}
					}
				 }
			  }
		  } // if (info.bi_mv[j])

		  GOPsz *= 2;
      } // for( j = info.tPyrLev - 1; j >= 0; j-- )

	  fprintf(fptab, "%ld \t %% GOP motion_bytes\n", GOP_motion_bytes);

    } 

    /*************/
    /* Bitplanes */
    /*************/

    newfat[GOPcounter] = 0;
  
    eff_GOPsz = info.eff_GOPsz;
    if( info.denoise_flag == YES ) {
      eff_GOPsz *= ( YUV420 == 1 ) ? 2 : 4;
    }
//	printf("num_subband = %d, eff_GOPsz = %d\n",num_subband, eff_GOPsz);

	GOP_subband_bytes = 0;
    
    for( j = 0; j < eff_GOPsz; j++ ) {
      subband_bytes = 0;
      for( nm = 0; nm < num_subband; nm++ ) {
        // read bit plane
        in_bytes += read_substream_length( &subband_fat[nm], fpbit );
        data =
          ( unsigned char * )getarray( subband_fat[nm],
                                       sizeof( unsigned char ), "data" );
        if( fread( data, sizeof( unsigned char ), subband_fat[nm], fpbit ) !=
            ( unsigned )subband_fat[nm] ) {
          printf( "read error3 %d %d\n", GOPcounter, subband_fat[nm] );
          exit( 1 );
        }
        in_bytes += sizeof( unsigned char ) * subband_fat[nm];  

        // write bit plane      
        fscanf( fpbytes, "%d\n", &tmp ); 
//		printf("tmp = %d\n",tmp);
        newsubband_fat[nm] = tmp;
        if( newsubband_fat[nm] != 0 ) {
          if( subband_fat[nm] < newsubband_fat[nm] ) {
            newsubband_fat[nm] = subband_fat[nm];
          }
          // write subband header
          bytes = write_substream_length( newsubband_fat[nm], outfp );
          out_bytes           += bytes;
          newfat[GOPcounter]  += bytes;
          subband_bytes       += bytes;
          total_subband_bytes += bytes;
      
          if( fwrite
              ( data, sizeof( unsigned char ), newsubband_fat[nm],
                outfp ) != ( unsigned )newsubband_fat[nm] ) {
            printf( "write error %d\n", newsubband_fat[nm] );
            exit( 1 );
          }
          out_bytes           += sizeof( unsigned char ) * newsubband_fat[nm];
          newfat[GOPcounter]  += sizeof( unsigned char ) * newsubband_fat[nm];
          subband_bytes       += sizeof( unsigned char ) * newsubband_fat[nm];
          total_subband_bytes += sizeof( unsigned char ) * newsubband_fat[nm];
        }
        free( data );
      }
      if ( pulled_bitstream == NO ){
        fprintf( fptab, "%ld \t %% subband_bytes\n", subband_bytes );
		GOP_subband_bytes += subband_bytes;
      }
    }
    
	fprintf(fptab, "%ld \t %% GOP subband_bytes\n", GOP_subband_bytes);

//	printf("newfat = %d\n",newfat[GOPcounter]);
	
    curr += info.GOPsz;
    GOPcounter++;
  }                           // while

//  assert((unsigned int)total_mv_bytes==sum_mv[info.t_level]); 
//  assert((unsigned int)total_mv_bytes<=exist_sum_mv[info.t_level]); 

  /*************/
  /* write FAT */
  /*************/
  fseek( outfp, FAT_start, SEEK_SET );
  
  if( fwrite( newfat, sizeof( long int ), num_of_GOP, outfp ) !=
      ( unsigned )num_of_GOP ) {
    printf( "write error \n" );
    exit( 1 );
  }
  
  for( i = 0; i < info.tPyrLev; i++ ) {
    free( simul_scene_change[i] );
  }
  free( simul_scene_change );

/*********************************************/
  fclose( fpbytes );
  fclose( fpbit );
  fclose( outfp );
  fclose( fptab );

  free( subband_fat );
  free( newsubband_fat );
  free( newfat );
}