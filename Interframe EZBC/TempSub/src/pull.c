/*
  This file pulls out bit stream from the pre-encoded bit stream

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "iostream"
#define EXTERN 
#include "structN.h"
#include "coderN.h"
#include "miscN.h"
#include "ioN.h"
#include "pstatN.h"
#include "basic.h"
#include "dpx.h"
#include "zlib.h" // zip-Library
#include "general.h"

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

void
pull_usage(  )
{
  printf( "pull pre_encoded_stream -t t_level -s s_level -r kbps [-C]\n" );
}

void
read_command( int argc, char **argv, short int *t_level, short int *s_level,
              int *kbps, int *AGP_exist, int *layer_exist, int *bi_exist, 
			  char *pre_stream, enum FLAG *CBR, enum FLAG *list, enum FLAG *VBR_ADAPT,
              char *list_name)
{
  int i, j, k, argnum = 1;

  *CBR = NO;
  *VBR_ADAPT = NO;
  *list = NO;
  *t_level = 0;
  *s_level = 0;
  list_name[0] = '\0';
  // how many sub-symbols exist
  for (j = 0; j < MAX_TLEVELS; j++)
  {
        AGP_exist[j]   = INVALID_AGP_EXIST;
		layer_exist[j] = INVALID_LAYER_EXIST;
		bi_exist[j]  = INVALID_BI_MV;  
  }
  
  for( i = 1; i < argc; i++ ) {
    if( *( argv[i] ) == '-' ) {

      switch ( *( ++argv[i] ) ) {
      default:
        printf( "-%c such an option is not available\n", *( argv[i] ) );
        pull_usage(  );
        exit( 1 );
      case 'h':
        pull_usage(  );
        exit( 1 );
        break;

      case 'y':   // the parameter for layer structure 
		j = 0;  
		while (j < MAX_TLEVELS ) {
			layer_exist[j] = atoi( argv[++i] );
			j++;
			if ( i==argc-1 )  // no more parameters 
				break;
			if ( *( argv[i+1] ) == '-' )  // next parameter
				break; 
		}
		if (j == MAX_TLEVELS) {
			printf("WARNING: Maxmimum number of layer_exist values reached!\n");
		}
		for (k = j; k < MAX_TLEVELS; k++) {
			layer_exist[k] = layer_exist[k-1];
		}
        break;

      case 'a':   // the parameter for AGP 
		j = 0;  
		while (j < MAX_TLEVELS ) {
			AGP_exist[j] = atoi( argv[++i] );
			j++;
			if ( i==argc-1 )  // no more parameters 
				break;
			if ( *( argv[i+1] ) == '-' )  // next parameter
				break; 
		}
		if (j == MAX_TLEVELS) {
			printf("WARNING: Maxmimum number of AGP_exist values reached!\n");
		}
		for (k = j; k < MAX_TLEVELS; k++) {
			AGP_exist[k] = AGP_exist[k-1];
		}
        break;
	  case 'b':  // the parameter for bi-directional motion field 
		j = 0;  
		while (j < MAX_TLEVELS ) {
			bi_exist[j] = atoi( argv[++i] );
			j++;
			if ( i==argc-1 )  // no more parameters 
				break;
			if ( *( argv[i+1] ) == '-' )  // next parameter
				break; 
		}
		if (j == MAX_TLEVELS) {
			printf("WARNING: Maxmimum number of bi_mv values reached!\n");
		}
		for (k = j; k < MAX_TLEVELS; k++) {
			bi_exist[k] = bi_exist[k-1];
		}
        break;

      case 't':
        *t_level = atoi( argv[++i] );
        break;
      case 's':
        *s_level = atoi( argv[++i] );
        break;
      case 'r':
        *kbps = atoi( argv[++i] );
        break;
      case 'C':
        *CBR = YES;
        break;
      case 'l': // name of list file for static output rates
        strcpy( list_name, argv[++i] );
        break;
      case 'n': // flag for display target rate list
        *list = YES;
        break;  
	  case 'V':
		  *VBR_ADAPT = YES;
		break;
      }
    } else {
      switch ( argnum ) {
      default:
        printf( "more parameters are specified\n" );
        pull_usage(  );
        exit( 1 );
      case 1:
        strcpy( pre_stream, argv[i] );
        argnum++;
        break;
      }
    }
  }
  
}


int
read_listfile( char *list_name, unsigned short int *list_rates, unsigned char *list_tlevel,
               unsigned char *list_slevel, unsigned char *list_simulcast, 
               char **list_simulname, videoinfo info )
{
  int i, counter = 0, rate = 0, tlevel = 0, slevel = 0, simulcast = 0;
  char iline[256], token[80];
  FILE *fplist;

  if( !( fplist = fopen( list_name, "rb" ) ) ) {
    printf( "can not open: %s\n", list_name );
    exit( 1 );
  }

  while( fgets( iline, 254, fplist ) ){
    sscanf( iline, "%s", token );
    
    if( !strcmp( token, "-listitem" ) ) {
      sscanf( iline, "%s%d%d%d%d%s", token, &tlevel, &slevel,
              &rate, &simulcast, list_simulname[counter] );
      list_tlevel[counter]    = tlevel;
      list_slevel[counter]    = slevel;
      list_rates[counter]     = rate;
      list_simulcast[counter] = simulcast;
      counter++;
    }      
  }
  if ( counter == 0 ){
    printf( "no listitems in: %s\n", list_name );
    exit( 1 );
  }
  
  for ( i = 1; i < counter; i++ ){
    if ( list_tlevel[i] < info.t_level ){
      printf( "target t_level too low in: %s\n", list_name );
      exit( 1 );
    } 
    if ( list_slevel[i] < info.s_level ){
      printf( "target s_level too low in: %s\n", list_name );
      exit( 1 );
    } 
    if ( list_rates[i] >= info.bitrate ){
      printf( "target rate too high in: %s\n", list_name );
      exit( 1 );
    } 
  }

  fclose( fplist );
  
  return counter;
}


// 在后缀为.mvby的文件中读取mv的每一帧bit数
unsigned long int
read_mv_bytes( char *seq_name, videoinfo info, int AGP_sign, int layer_sign,
			  unsigned long int *exist_sum_mv )
{
  int j, k, remaining_frs, curr, last;
  int mv, mv1, mv2, GOPsz, eff_GOPsz[20];
  unsigned long int sum_mv;
  char mvstatname[250];
  FILE *FID;
  int  major_mv, sign_bytes, sub_bit, sub_mv, total_mv ; 
  
  GOPsz = info.GOPsz;

  // determine effective GOP size in level j 决定有效的gop的大小
  for( j = 0; j < info.tPyrLev; j++ ) {
    eff_GOPsz[j] = GOPsz;  
    GOPsz /= 2;
  }
  
  curr = info.start;
  last = info.last;
  // *exist_sum_mv: the number of existing motion vector bytes up to this temporal level
  sum_mv = *exist_sum_mv = 0; 
  sprintf( mvstatname, "%s.mvby", seq_name );    // Hanke, 16.09.02
  FID = fopen( mvstatname, "rb" );
  if( FID == NULL ) {
    printf( "can not open file mvby\n" );
    exit( 1 );
  }
  
  while( curr <= last ) {     // 处理完所有帧
    remaining_frs = last - curr + 1; 

    for( j = info.tPyrLev - 1; j >= 0; j-- ) { // 处理一个gop
	  if ( info.bi_exist[j] )  // bi-directional motion field 双向运动场
	  {
		  for( k = 0; k <= eff_GOPsz[j]; k++ ) {
            // scalable motion vector coding 
//			  printf("eff_GOPsz[j] = %d\n",eff_GOPsz[j]);
			if ( !AGP_sign  && !layer_sign) // no AGP, no layer structure
			{
				fscanf( FID, "%d\n", &mv );
				if( j >= info.t_level )  
				{
					sum_mv += mv;
					*exist_sum_mv += mv; // here *exist_sum_mv is the same as sum_mv
				}
			}

			if ( AGP_sign  && !layer_sign ) // AGP,    no layer structure 
			{
				fscanf( FID, "%d ", &major_mv );
				if (major_mv==0) continue;  // no motion vector, i.e. scene change 
				if( j >= info.t_level ) sum_mv += major_mv;   // major_mv is always necessary
				if (info.AGP_level[j])  // AGP in this temporal level 
				{
					fscanf( FID, "%d ", &sign_bytes); 
					if (info.AGP_exist[j]  && j >= info.t_level ) sum_mv += sign_bytes;
					for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
					{
						fscanf( FID, "%d ", &sub_mv); 
						if (sub_bit<info.AGP_exist[j] && j >= info.t_level ) sum_mv += sub_mv;
					}
					fscanf( FID, "%d ", &total_mv); // total bytes for this set of motion vectors 
					if (j>=info.t_level) *exist_sum_mv += total_mv ; 
				}else
					if (j>=info.t_level) *exist_sum_mv += major_mv ;
			}

			if ( !AGP_sign  && layer_sign ) // no AGP,    layer structure 
			{
				fscanf( FID, "%d ", &major_mv );
				if (major_mv==0) continue;  // no motion vector, i.e. scene change 
				if( j >= info.t_level ) sum_mv += major_mv;   // major_mv is always necessary
				if (info.layer_mv[j])  // layer structure in this temporal level 
				{
					fscanf( FID, "%d ", &sub_mv); 
					if (info.layer_exist[j] && j >= info.t_level ) sum_mv += sub_mv;
					fscanf( FID, "%d ", &total_mv); // total bytes for this set of motion vectors 
					if (j>=info.t_level) *exist_sum_mv += total_mv ; 
				}else
					if (j>=info.t_level) *exist_sum_mv += major_mv ;
			}

			if ( AGP_sign  && layer_sign ) // AGP,     layer structure 
			{
				fscanf( FID, "%d ", &major_mv );
				if (major_mv==0)		continue;  // no motion vector, i.e. scene change 
				if( j >= info.t_level ) sum_mv += major_mv;   // major_mv is always necessary
				if (info.AGP_level[j])  // AGP in this temporal level 
				{
					fscanf( FID, "%d ", &sign_bytes); 
					if (info.AGP_exist[j]  && j >= info.t_level ) sum_mv += sign_bytes;
					for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
					{
						fscanf( FID, "%d ", &sub_mv); 
						if (sub_bit<info.AGP_exist[j] && j >= info.t_level ) sum_mv += sub_mv;
					}
				}

				if ( info.layer_mv[j])
				{
					fscanf( FID, "%d ", &sub_mv); 
					if ( info.layer_exist[j] && j >= info.t_level ) sum_mv += sub_mv;
					if ( info.AGP_level[j])
					{
						fscanf( FID, "%d ", &sign_bytes); 
						if (info.AGP_exist[j]  && j >= info.t_level ) sum_mv += sign_bytes;
						for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
						{
							fscanf( FID, "%d ", &sub_mv); 
							if (sub_bit<info.AGP_exist[j] && j >= info.t_level ) sum_mv += sub_mv;
						}
					}

				}

				if ( info.layer_mv[j] || info.AGP_level[j] )
				{
					fscanf( FID, "%d ", &total_mv); // total bytes for this set of motion vectors 
					if (j>=info.t_level) *exist_sum_mv += total_mv ; 
				}else
					if (j>=info.t_level) *exist_sum_mv += major_mv ;
			}

		  }
	  } // if ( info.bi_exist[j] )

	  // alternative Haar reconstruction 
	  if ( ! info.bi_exist[j] )  // not bi-directional motion field, alternative reconstruction with Harr 
	  {
		  assert(0);
		  k=0; 
		  while ( k <= eff_GOPsz[j] )
		  {
			  if (AGP_sign==0)  // without AGP ( non-scalable motion vectors )
			  {
				  if (k==0)  // the first set of motion vectors are always 0
				  {
					  fscanf( FID, "%d\n", &mv );
					  assert(mv==0);
					  k++; 
				  }else
				  {
					  fscanf( FID, "%d\n", &mv1 );  // left set mv
					  k++;
					  fscanf( FID, "%d\n", &mv2 );  // right set mv
					  k++;
					  if( j >= info.t_level ) 
					  {
						  // bi-directional motion vector or left set motion vector
						  if ( ( mv1 && mv2 ) || ( mv1 && !mv2) )
						  {
							  sum_mv += mv1;  // only left set is reserved finally
							  *exist_sum_mv += mv1+mv2;
						  }
						  
						  if (!mv1 && mv2)  // right set motion vector
						  {
							  sum_mv += mv2;  
							  *exist_sum_mv += mv1+mv2;
						  }
					  }
				  }
			  }else  // with AGP
			  {

			  } // with AGP
		  }
	  }  // if ( ! info.bi_mv[j] )

    }
    
    curr += info.GOPsz;
  }
  fclose( FID );
  
  return sum_mv;
}


void
gop_read_mv_bytes( char *seq_name, videoinfo info, int AGP_sign, int layer_sign,
			  unsigned long int *gop_exist_sum_mv, int num_of_GOP, unsigned long int *gop_sum_mv )
{
  int i, j, k, remaining_frs, curr, last;
  int mv, mv1, mv2, GOPsz, eff_GOPsz[20];
  unsigned long int sum_mv;
  char mvstatname[250];
  FILE *FID;
  int  major_mv, sign_bytes, sub_bit, sub_mv, total_mv ; 
  
  GOPsz = info.GOPsz;

  // determine effective GOP size in level j
  for( j = 0; j < info.tPyrLev; j++ ) {
    eff_GOPsz[j] = GOPsz;  
    GOPsz /= 2;
  }
  
  curr = info.start;
  last = info.last;

  assert( (last - curr + 1) == num_of_GOP*info.GOPsz );



  // *exist_sum_mv: the number of existing motion vector bytes up to this temporal level
  sum_mv = *gop_exist_sum_mv = 0; 
  sprintf( mvstatname, "%s.mvby", seq_name );    // Hanke, 16.09.02
  FID = fopen( mvstatname, "rb" );
  if( FID == NULL ) {
    printf( "can not open file mvby\n" );
    exit( 1 );
  }
  
  for(i = 0; i < num_of_GOP; i ++){
	gop_sum_mv[i] = 0;
  }

  i = 0;

  while( curr <= last ) {     
    remaining_frs = last - curr + 1; 
    
    for( j = info.tPyrLev - 1; j >= 0; j-- ) {

	  if ( info.bi_exist[j] )  // bi-directional motion field 
	  {
		  for( k = 0; k <= eff_GOPsz[j]; k++ ) {

            // scalable motion vector coding 

			if ( !AGP_sign  && !layer_sign) // no AGP, no layer structure
			{
				fscanf( FID, "%d\n", &mv );
				if( j >= info.t_level )  
				{
					sum_mv += mv;
					gop_sum_mv[i] += mv; // here *exist_sum_mv is the same as sum_mv
				}
			}

		  }
	  } // if ( info.bi_exist[j] )

    }
    
    curr += info.GOPsz;
	i ++;
  }
  fclose( FID );

}

// automatically allocate the bits between motion vector and subband data
// according to some empirically bit-allcoation model 
void automatic_bit_allocation_with_AGP(unsigned long int *sum_mv, videoinfo *info,
									   char *seq_name, int AGP_sign, int layer_sign)
{
	double budget; 
	double mv_percent; 
	videoinfo info_bak;
	unsigned long int exist_sum_mv;
	int i; 

	budget = ( double )info->bitrate * 1000 / 8 * ( info->last - info->start + 1 ) / info->framerate;
	mv_percent = *sum_mv/budget;  // the percentage of motion vector bytes in the available budget
	if (mv_percent<=MV_PERCENT_THERSH)
		return;   // tolerable percentage according empirical bit allocation model 
	else
	{
		// discard motion vector sub-symbols starting from the top temporal level 
		for (i=info->t_level; i<=info->tPyrLev-1; i++)
		{
			while (info->AGP_exist[i]>0 && mv_percent>MV_PERCENT_THERSH)
			{
				(info->AGP_exist[i])--;
				info_bak = *info; 
				*sum_mv = read_mv_bytes( seq_name, info_bak, AGP_sign, layer_sign, &exist_sum_mv );
				mv_percent = *sum_mv/budget;
			}
			if (mv_percent<=MV_PERCENT_THERSH)
				return;
		}
		for (i=info->t_level; i<=info->tPyrLev-1; i++)
			assert( info->AGP_exist[i]==0 );
	}
}

/*
unsigned long int 
write_zip( char *bit_alloc_name, FILE *fphead, enum FLAG compress_flag )
{
  int err;
  unsigned long int input_size = 0, in_len, out_len = 1000000000;
  unsigned char *input_data , *compressed;
  FILE *fpbytes;

  // read input_file
  if( !( fpbytes = fopen( bit_alloc_name, "rb" ) ) ) {
    printf( "can not open %s\n", bit_alloc_name );
    exit( 1 );
  }
  fseek( fpbytes, 0, SEEK_END );
  input_size = ftell( fpbytes );
  fseek( fpbytes, 0, SEEK_SET );

  input_data =
    ( unsigned char * )getarray( input_size, sizeof( unsigned char ),
                                 "input_data" );
  if( fread( input_data, sizeof( unsigned char ), input_size, fpbytes ) !=
      ( unsigned ) input_size ) {
    printf( "fread error \n" );
    exit( 1 );
  }
  fclose ( fpbytes );

  in_len = input_size * sizeof( unsigned char );
 
  compressed =
      ( unsigned char * )getarray( in_len, sizeof( unsigned char ),
                                   "compressed" );   
  
  if( compress_flag == YES ){
    err = compress(compressed, &out_len, (Bytef*)input_data, in_len);
    if (err) {
      printf("error during zip compression.\n");
      exit(1);
    }
    assert ( out_len <= in_len );
  
    // write zip to overhead stream
    fwrite( compressed, sizeof( unsigned char ), out_len, fphead );
  } else {
    // write uncompressed input_data to overhead stream
    fwrite( input_data, sizeof( unsigned char ), in_len, fphead );
  }
  
  free( input_data );
  free( compressed );
  return ( compress_flag == YES ) ? out_len : in_len;
}
*/
/*
void
extract_zip( char *bit_alloc_name, unsigned char *compressed, 
             unsigned long int in_len, enum FLAG decompress_flag )
{
  int err;
  unsigned long int out_len = 20000000;
  unsigned char *decompressed;
  FILE *fpbytes;

  decompressed =
    ( unsigned char * )getarray( out_len, sizeof( unsigned char ),
                                 "decompressed" );
  
  if( decompress_flag == YES ){
    // extract zip archive
    err = uncompress((Bytef*)decompressed, &out_len, compressed, in_len);
    if (err) {
      printf("error during zip uncompression.\n");
      exit(1);
    }
  }
  
  // write data to file
  if( !( fpbytes = fopen( bit_alloc_name, "wb" ) ) ) {
    printf( "can not open %s\n", bit_alloc_name );
    exit( 1 );
  }

  if( decompress_flag == YES ){
    fwrite( decompressed, sizeof( unsigned char ), out_len, fpbytes );
  } else {
    fwrite( compressed, sizeof( unsigned char ), in_len, fpbytes );
  }
 
  fclose ( fpbytes );
  free ( decompressed );
}
*/

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



int
main( int argc, char **argv )
{
  int i, j, k, nm, curr, last, num_of_GOP, GOPheader_bytes, Per_GOPheader_bytes, FAT_bytes;
  int count = 0, remaining_frs, num_subband, bytes, mvBytes;
  int kbps = 0, tmp, GOPsz, eff_GOPsz, is_pullable;
  int *subband_fat, *newsubband_fat;
  long int in_bytes = 0, out_bytes = 0, FAT_start = 0, overhead_start;
  long int overhead = 0, estimated_overhead = 0, GOP_estimated_overhead = 0, total_mv_bytes = 0;
  long int subband_bytes = 0, total_subband_bytes = 0, *newfat; 
  short int t_level, s_level, rel_s_level, target_rate;

  int getnum;
  
  unsigned char *data;
  char pre_stream[256], post_stream[256], overhead_stream[256];
  char seq_name[256], seq_name2[256], bit_alloc_name[256], bit_alloc_table[256];
  char list_name[256], simulcast_name[256], **list_simulname;
  
  int list_size = 0, list_counter = 0, max_tlevel = 0;
  unsigned char list_tlevel[25], list_slevel[25], list_qp[25], list_simulcast[25];
  unsigned short int list_rates[25], simul_bitstream = 0;
  unsigned long int sum_mv[10], exist_sum_mv[10], list_out_bytes[25], simul_bitstream_pos = 0;

  unsigned long int *gop_sum_mv[10], *gop_exist_sum_mv[10];
                              
  enum FLAG CBR, display_list, pulled_bitstream, VBR_ADAPT; 
  FILE *outfp, *fpbytes, *fptab, *fphead;
  int   AGP_sign, sub_bit, AGP_exist_pull[MAX_TLEVELS], bi_exist[MAX_TLEVELS]; 
  long int major_mvBytes; 
  enum FLAG scene_change_left, scene_change_right; 
  int   layer_sign, layer_exist_pull[MAX_TLEVELS]; 

  int GOP_subband_bytes;
  
  videoinfo info;
  
  printf("Built on %s at %s.\n", __DATE__, __TIME__);

  read_command( argc, argv, &t_level, &s_level, &kbps, AGP_exist_pull, layer_exist_pull,
	            bi_exist, pre_stream, &CBR, &display_list, &VBR_ADAPT, list_name );

  if( kbps <= 0 ) {
    printf( "kbps parameter is not right\n" );
    pull_usage(  );
    exit( 1 );
  }

  /**********/
  /* header */
  /**********/
  
  read_header( pre_stream, &info );// 读头信息，并将头信息转换为info
  in_bytes += sizeof( videoheader );

  num_of_GOP = get_GOP_num( info );
  printf("num_of_GOP = %d\n",num_of_GOP);

  if(CBR == YES || VBR_ADAPT == YES){
	  for(i = 0; i < 10; i ++){
		  gop_sum_mv[i] = ( unsigned long * )getarray( num_of_GOP, sizeof( unsigned long ), "gop_sum_mv" );
		  gop_exist_sum_mv[i] = ( unsigned long * )getarray( num_of_GOP, sizeof( unsigned long ), "gop_sum_mv" );
	  }
  }

  // scalable motion vector coding 可伸缩运动向量编码,以及可伸缩分辨率编码
  AGP_sign = layer_sign = 0; 
  for (i = 0; i < MAX_TLEVELS; i++) 
  {
	  if ( info.AGP_level[i] ) 	  
		  AGP_sign = 1; 
	  // no value for AGP_exist_pull, then we use default value for it
	  if ( AGP_exist_pull[i]== INVALID_AGP_EXIST ) AGP_exist_pull[i] = info.AGP_level[i];
	  // for resolution scalability, when we reduce one level of resolution, we discard one least 
	  // significant sub-symbol until no sub-symbol can be discarded
	  // hence, AGP_exist_pull can not be larger than info.AGP_level[i]-s_level
	  if ( AGP_exist_pull[i]>info.AGP_level[i]-s_level ) 
		  AGP_exist_pull[i]= MY_MAX (0, info.AGP_level[i]-s_level); 
	  info.AGP_exist[i] = AGP_exist_pull[i]; 

	  if ( info.layer_mv[i]  )     layer_sign = 1;
      // no value for layer_exist_pull, the we use the default values 
	  if ( layer_exist_pull[i]== INVALID_LAYER_EXIST ) layer_exist_pull[i] = 1;
	  if ( layer_exist_pull[i]>info.layer_mv[i] ) 
		  layer_exist_pull[i]= info.layer_mv[i]; 
	  info.layer_exist[i] = layer_exist_pull[i]; 

	  // alternative Haar reconstruction 
	  // only when info.bi_mv[i] == 0, mv on RIGHT side can be discarded 
	  assert( bi_exist[i] == INVALID_BI_MV || bi_exist[i] == 0 || bi_exist[i] == 1 );
	  assert( info.bi_exist[i] == 1 || info.bi_exist[i] == 0 ); 
	  assert( info.bi_mv[i]== 1 || info.bi_mv[i] == 0 ); 
	  if ( bi_exist[i] != INVALID_BI_MV && info.bi_mv[i]==0)
		  info.bi_exist[i] = bi_exist[i]; 
  }
  
  pulled_bitstream = ( info.bitrate != 0 ) ? YES : NO;

  if ( pulled_bitstream == NO ){
    info.bitrate   = kbps;
    info.t_level   = max_tlevel = t_level;
    info.s_level   = s_level;

  }
  
  rel_s_level    = info.s_level - (info.denoise_flag == YES);// 

  if( strlen(list_name) > 0 ){  //NOT WORKING
    if( pulled_bitstream == YES ){
      printf("bitstream already pulled before - please use -n option to display pullable rates.\n");
      exit( 1 );
    }
    // get array for simulcast_names
    list_simulname = ( char ** ) malloc( 25 * sizeof ( char * ) );
    for( i = 0; i < 25; i++ ) { 
      list_simulname[i] = ( char * ) malloc( 256 * sizeof ( char ) );
    }

    // read in list file for output rates etc.
    list_counter = read_listfile( list_name, list_rates, list_tlevel, list_slevel, 
                                  list_simulcast, list_simulname, info );
    if( list_counter > 0 ){ 
      for ( i = 0; i < list_counter; i++ ){
        if ( list_tlevel[i] > max_tlevel){
          max_tlevel = list_tlevel[i];
        } 
      }
    }
  } //NOT WORKING

  strcpy( info.bitname, pre_stream );
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 );
  seq_name[strlen( info.bitname ) - 4] = '\0';

  // GOPheader_bytes (cmp. mvcoding.c)
  GOPsz = info.GOPsz;
  for( i = 0; i < info.tPyrLev; i++ ) {
    for( j = 0; j < GOPsz; j++ ) 
      count++;// 记录有多少个时域子带
    GOPsz /= 2;  
  }
  GOPsz = info.GOPsz; 
  GOPheader_bytes = (( count + 7 ) / 8) * num_of_GOP; // 每个时域子带一个bit。记录整个序列需要多少字节

  Per_GOPheader_bytes = GOPheader_bytes / num_of_GOP; 

  printf("num_of_GOP = %d, GOPsz = %d\n",num_of_GOP,GOPsz);
  
  // FAT bytes
  FAT_bytes = num_of_GOP * sizeof( long int ); // DENOISE ? -> only larger values 每个gop占一个long int

  if( pulled_bitstream == NO ){
    // read MV bytes
     
	info.org_yheight = info.yheight;
	info.org_ywidth  = info.ywidth; 

	//printf("max_tlevel = %d\n",max_tlevel);

	for ( i = 0; i <= max_tlevel; i++ ){
	  getnum = 0;
		
	  info.t_level = i;
	  sum_mv[i] = read_mv_bytes( seq_name, info, AGP_sign, layer_sign, &exist_sum_mv[i] ); // 读取每一帧mv的字节数
	  printf("%d \t %d \n", sum_mv[i], exist_sum_mv[i]);
	  if(CBR == YES || VBR_ADAPT == YES){
		  gop_read_mv_bytes( seq_name, info, AGP_sign, layer_sign,
				  gop_exist_sum_mv[i], num_of_GOP, gop_sum_mv[i] );

		  for(j = 0; j < num_of_GOP; j ++)
			getnum += gop_sum_mv[i][j];

		  printf("getnum = %d, sum_mv = %d\n",getnum, sum_mv[i]);
	  }
	}
	info.t_level = t_level;

//	if (AGP_sign && layer_sign) // automatic bit allocation according to empirical model 
//		automatic_bit_allocation_with_AGP(&sum_mv[t_level], &info, seq_name, AGP_sign, layer_sign ); 
    
    sprintf( bit_alloc_table, "%s_%d.alloc_table", seq_name, info.bitrate );
    if( !( fptab = fopen( bit_alloc_table, "wb" ) ) ) {
      printf( "can not open: %s\n", bit_alloc_table );
      exit( 1 );
    }
/*    
    if( list_counter > 0 ){ //NOT WORKING
		list_size = sizeof( unsigned char ) + 2 * sizeof( unsigned long int ) 
        + list_counter * ( sizeof( unsigned short int ) + 4 * sizeof( unsigned char )
                           + sizeof( unsigned long int ) ) ; // List + overhead FAT
      estimated_overhead += list_size;
      
      printf( "---------------------------------------------------------------\n" );
      printf( " List file : %s\n", list_name);
      printf( " List size : %d bytes\n", list_size ); 
      fprintf( fptab, "%%---------------------------------------------------------------\n" );
      fprintf( fptab, "%s %% list_file\n", list_name );
      fprintf( fptab, "%d \t %% list_size\n", list_size );
      
      // open overhead substream
      sprintf( overhead_stream, "%s.overhead", seq_name );
      if( !( fphead = fopen( overhead_stream , "wb" ) ) ) {
        printf( "can not open: %s\n", overhead_stream );
        exit( 1 );
      }

	  printf("list_counter = %d\n",list_counter);
 
      for ( i = 0; i < list_counter; i++ ){
        if( list_simulcast[i] ){
          // SIMULCAST --> determine filesize for overhead
		  assert(0);
          strncpy( seq_name2, list_simulname[i], strlen( list_simulname[i] ) - 4 );
          seq_name2[strlen( list_simulname[i] ) - 4] = '\0';
         
          is_pullable = simulcast_pull( list_simulname[i], seq_name2, list_tlevel[i],
                                        list_slevel[i], list_rates[i], 0 );
          if ( is_pullable == -1 ){
            printf("Simulcast-File %s not pullable\n", list_simulname[i]);
            exit(1);
          }
          list_qp[i] = is_pullable;
          
          // write file(s)
          if( list_simulcast[i] == 1){
            sprintf( simulcast_name, "%s_%d.bit", seq_name2, list_rates[i] );
             list_out_bytes[i] = write_zip( simulcast_name, fphead, NO );
            estimated_overhead += list_out_bytes[i];
          } else if( list_simulcast[i] == 2) {
            sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name2, list_rates[i] );
             list_out_bytes[i] = write_zip( bit_alloc_name, fphead, YES );
            estimated_overhead += list_out_bytes[i];
          } else {
            printf("Unknown Simulcast-Mode!\n");
            exit (1);
          }
          
          printf( " Target %2d : t_level %d, s_level %d, rate %5d, qp %2d, simulcast %d --> %6ld bytes\n", 
                  i, list_tlevel[i], list_slevel[i], list_rates[i], list_qp[i], list_simulcast[i], 
                  list_out_bytes[i] );
          fprintf( fptab, "%ld \t %% target %d : t_level %d, s_level %d, rate %5d, qp %2d, simulcast %d\n", 
                   list_out_bytes[i], i, list_tlevel[i], list_slevel[i], list_rates[i], list_qp[i],
                   list_simulcast[i]);
        } else {
          info.t_level = list_tlevel[i];
          info.s_level = list_slevel[i];
          info.bitrate = list_rates[i];
          overhead     = sizeof( videoheader ) + GOPheader_bytes + FAT_bytes;
          
          // bit_alloc file (bytes_per_GOP) calculate for every target rate
          if( CBR == NO ){          
            is_pullable = bit_alloc_VBR( info, &info.qp, sum_mv[info.t_level], overhead );
            
            if (!is_pullable){
              list_simulcast[i] = 1;
              i--;
            } else {
              sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, info.bitrate );
              list_out_bytes[i] = write_zip( bit_alloc_name, fphead, YES );
              estimated_overhead += list_out_bytes[i];
              list_qp[i] = info.qp; 
              printf( " Target %2d : t_level %d, s_level %d, rate %5d, qp %2d, simulcast %d --> %6ld bytes\n", 
                      i, list_tlevel[i], list_slevel[i], list_rates[i], list_qp[i], list_simulcast[i],
                      list_out_bytes[i] );
              fprintf( fptab, "%ld \t %% target %d : t_level %d, s_level %d, rate %5d, qp %2d, simulcast %d\n", 
                       list_out_bytes[i], i, list_tlevel[i], list_slevel[i], list_rates[i], list_qp[i],
                       list_simulcast[i]);
            }
          }
        }
      }
      fclose ( fphead );
        
      // reset values to main output bit stream
      info.t_level = t_level;
      info.s_level = s_level;
      info.bitrate = kbps;
    } //NOT WORKING
*/
  
    estimated_overhead += sizeof( videoheader ) + GOPheader_bytes + FAT_bytes;

	printf("GOPheader_bytes = %d\n",GOPheader_bytes);

	GOP_estimated_overhead = estimated_overhead / num_of_GOP;// 每个gop的overhead
    
    if( CBR == NO ){
	  if(VBR_ADAPT == NO)	
		is_pullable = bit_alloc_VBR( info, &info.qp, sum_mv[info.t_level], estimated_overhead );// 全是NO
	  else
		is_pullable = bit_alloc_VBR_v2(info, &info.qp, estimated_overhead, gop_sum_mv[info.t_level], info.last, info.last - info.GOPsz);
	
	}else if(CBR == YES){
		is_pullable = bit_alloc_CBR( info, &info.qp, estimated_overhead, gop_sum_mv[info.t_level], info.last );
	}
    
    if( !is_pullable ){
      printf("bitstream is not pullable!\n");
      exit(1);
    }
  }   // end of pulled_bitstream == NO
  
  if( !( fpbit = fopen( pre_stream, "rb" ) ) ) {
    printf( "can not open: %s\n", pre_stream );
    exit( 1 );
  }
/*       
  // get and extract data from pulled bitstream
  if( pulled_bitstream == YES ){ //NOT WORKING
    FAT_start = in_bytes;
    fseek( fpbit, -sizeof( unsigned long int ), SEEK_END );
    fread( &overhead_start, sizeof( unsigned long int ), 1, fpbit ); 
    estimated_overhead = ftell( fpbit ) - overhead_start;
    fseek( fpbit, overhead_start, SEEK_SET );
    in_bytes = overhead_start;
    
    // read target rate list 
    fread( &list_counter, sizeof( unsigned char ), 1, fpbit );
    in_bytes += sizeof( unsigned char );
    for ( i = 0; i < list_counter; i++ ){
      fread( &list_rates[i], sizeof( unsigned short int ), 1, fpbit );
      fread( &list_tlevel[i], sizeof( unsigned char ), 1, fpbit );
      fread( &list_slevel[i], sizeof( unsigned char ), 1, fpbit );
      fread( &list_qp[i], sizeof( unsigned char ), 1, fpbit );
      fread( &list_simulcast[i], sizeof( unsigned char ), 1, fpbit );
      fread( &list_out_bytes[i], sizeof( unsigned long int ), 1, fpbit );
      in_bytes += sizeof( unsigned short int )+ 4 * sizeof( unsigned char )
        + sizeof( unsigned long int );
    }
    
    unsigned long overhead_size = 0;
    unsigned char *input_data;
 
    // read in total overhead size
    fread( &overhead_size, sizeof( unsigned long int ), 1, fpbit );
    in_bytes += sizeof( unsigned long int );
    
    target_rate = -1;
    for ( i = 0; i < list_counter; i++ ){
      if( list_rates[i] == kbps ){
        target_rate = i;
        info.bitrate = list_rates[i];
        info.t_level = list_tlevel[i];
        info.s_level = list_slevel[i];
        info.qp = list_qp[i];
        rel_s_level  = info.s_level - (info.denoise_flag == YES);
        printf( " Target %2d : t_level %d, s_level %d, rate %5d, qp %2d\n", 
                i, list_tlevel[i], list_slevel[i], list_rates[i], list_qp[i] );
        break;
      } else {
        if ( list_simulcast[i] == 1 ){
          simul_bitstream = i;
          simul_bitstream_pos = in_bytes;
        }
        in_bytes += list_out_bytes[i];
        fseek( fpbit, in_bytes, SEEK_SET );
      }
    }
    if ( target_rate < 0 || display_list == YES ){
      printf( "possible target output configurations: \n" );
      printf( " Mainstream : t_level %d, s_level %d, rate %5d  ( decodeable without pull )\n", 
              info.t_level, info.s_level, info.bitrate);
      for ( i = 0; i < list_counter; i++ ){
        printf( " Target %2d  : t_level %d, s_level %d, rate %5d, qp %2d, simulcast %d --> %6ld bytes\n", 
                i, list_tlevel[i], list_slevel[i], list_rates[i], list_qp[i], list_simulcast[i],
                list_out_bytes[i] );
      }
      exit (1);
    }
    
    // read in zipped overhead information file
    fseek( fpbit, in_bytes, SEEK_SET );
        
    input_data =
      ( unsigned char * )getarray( list_out_bytes[target_rate], sizeof( unsigned char ),
                                   "input_data" );
    if( fread( input_data, sizeof( unsigned char ), list_out_bytes[target_rate], fpbit ) !=
        ( unsigned ) list_out_bytes[target_rate] ) {
      printf( "fread error \n" );
      exit( 1 );
    }
     
    // extract and save overhead file
    if( list_simulcast[i] ){
		assert(0);
      if( list_simulcast[i] == 2 ){ // extract bytes_per_GOP file
        sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, list_rates[target_rate] );
        extract_zip( bit_alloc_name, input_data, list_out_bytes[target_rate], YES ); 
        
        free( input_data );
        // and extract simulcast main bitstream
        fseek( fpbit, simul_bitstream_pos, SEEK_SET );
        
        input_data =
          ( unsigned char * )getarray( list_out_bytes[simul_bitstream], sizeof( unsigned char ),
                                       "input_data" );
        if( fread( input_data, sizeof( unsigned char ), list_out_bytes[simul_bitstream], fpbit ) !=
            ( unsigned ) list_out_bytes[simul_bitstream] ) {
          printf( "fread error \n" );
          exit( 1 );
        }
        sprintf( simulcast_name, "%s_%d.bit", seq_name, list_rates[simul_bitstream] );
        extract_zip( simulcast_name, input_data, list_out_bytes[simul_bitstream], NO );
             
        // proceed with pulling simulcast-files for target rate
        is_pullable = simulcast_pull( simulcast_name, seq_name, list_tlevel[target_rate],
                                      list_slevel[target_rate], list_rates[target_rate], 
                                      list_qp[target_rate] );
        if ( is_pullable == -1 ){
          printf("Simulcast-File %s not pullable\n", simulcast_name);
          exit(1);
        }
      } else { // only extract simulcast main bitstream
        sprintf( simulcast_name, "%s_%d.bit", seq_name, list_rates[target_rate] );
        extract_zip( simulcast_name, input_data, list_out_bytes[target_rate], NO );
        printf("Finished.\n");
      }
      exit(1);
    } else{
      sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, list_rates[target_rate] );
      extract_zip( bit_alloc_name, input_data, list_out_bytes[target_rate], YES );
    }

    // go to start of bitstream
    in_bytes = FAT_start;
    fseek( fpbit, FAT_start, SEEK_SET );
    
    free( input_data );
  } //NOT WORKING
*/
  sprintf( post_stream, "%s_%d.bit", seq_name, info.bitrate ); 
  
  write_header( post_stream, info );
  out_bytes += sizeof( videoheader );
  
  if( !( outfp = fopen( post_stream, "r+b" ) ) ) {
    printf( "can not open: %s\n", post_stream );
    exit( 1 );
  }
  
  printf( "---------------------------------------------------------------\n" );
  printf( " Sequence name     : %s\n", seq_name );
  printf( " Sequence length   : %03d frames ( %03d - %03d )\n",
          info.last - info.start + 1, info.start, info.last );
  printf( " Output bit rate   : %d kbit/sec\n", info.bitrate );
  printf( " Output frame rate : %0.2f frames/sec (t_level %d)\n",
          ((double)(info.framerate) / pow((float)2, info.t_level)), info.t_level );
  printf( " Output frame size : (Y) %d x %d (C) %d x %d (s_level %d)\n",
          info.ywidth >> rel_s_level, info.yheight >> rel_s_level,
          info.cwidth >> rel_s_level, info.cheight >> rel_s_level,
          info.s_level );
  printf( " Number of GOPs    : %d\n", num_of_GOP );
  printf( " Header/FAT bytes  : %d bytes\n", 
          sizeof( videoheader ) + FAT_bytes + GOPheader_bytes );
  // printf( " Header size       : %d bytes\n", sizeof( videoheader ) );
  // printf( " GOP FAT bytes     : %d bytes\n", FAT_bytes ); 
  // printf( " GOP header bytes  : %d bytes\n", GOPheader_bytes );
  printf( " Overhead bytes    : %ld bytes\n", estimated_overhead
          - (sizeof( videoheader ) - FAT_bytes - GOPheader_bytes) 
          * ( pulled_bitstream == NO ));
#ifdef MY_SCALE
  printf( " MY_SCALE          : %d (%d extra bits)\n", MY_SCALE, (int)log2(MY_SCALE));
#endif
  printf( "---------------------------------------------------------------\n" );
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
  
  scene_change =
    ( enum FLAG ** )getarray( info.tPyrLev, sizeof( enum FLAG * ),
                              "scene_change" );
  GOPsz = info.GOPsz;
  for( i = 0; i < info.tPyrLev; i++ ) {
    scene_change[i] =
      ( enum FLAG * )getarray( GOPsz + 1, sizeof( enum FLAG ), "scene_change" );
    GOPsz /= 2;
  }
  
  curr = info.start;
  last = info.last;

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
      
    printf( " pulling GOP %d: frame %d - %d\n", GOPcounter,
            curr, curr + info.eff_GOPsz - 1 );
    if ( pulled_bitstream == NO ){
      fprintf( fptab, "%% GOP %2d %% frame %d - %d\n", GOPcounter, 
               curr, curr + info.eff_GOPsz - 1 );
    }
    
    GOPheader_bytes = 0;
    if( info.tPyrLev >= 1 ) {

      /*************/
      /* GOPheader */
      /*************/

      GOPheader_bytes = read_GOPheader( scene_change, info );   //printf("%d\n", GOPheader_bytes);
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

      GOP_motion_bytes = 0; //Added on 10.06.2017

      GOPsz = 2;
      for( j = info.tPyrLev - 1; j >= 0; j-- ) {

		  if (info.bi_exist[j]) // bi-directional motion field 
		  {
			  for( k = 0; k <= GOPsz; k++ ) {
				if( scene_change[j][k] == NO ) { 

					// scalable motion vector coding: AGP and layer structure 

                	if ( !AGP_sign && !layer_sign )  // no AGP, no layer structure 
					{
						mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
						in_bytes += 4;

						printf("mvBytes = %d\n",mvBytes);
						if( j >= info.t_level ) { // temporal scalability
						  process_mv_section(mvBytes, outfp, &out_bytes, &total_mv_bytes, &in_bytes, fpbit, 
							  fptab, pulled_bitstream);
						} else {
						  fseek( fpbit, mvBytes, SEEK_CUR );
						  in_bytes += mvBytes;
						}
					}

					if ( AGP_sign && !layer_sign )  // AGP, no layer structure 
					{
						major_mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
						in_bytes += 4;
						if( j >= info.t_level ) { // temporal scalability
							process_mv_section(major_mvBytes, outfp, &out_bytes, &total_mv_bytes, &in_bytes, fpbit, 
								  fptab, pulled_bitstream);
							if ( info.AGP_level[j] )
							{
								if (info.AGP_exist[j]) // sign bytes
									get_bytes_section(fpbit, outfp, fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
								else
									get_bytes_section(fpbit, outfp, fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
								for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
									if ( sub_bit < info.AGP_exist[j])
										  get_bytes_section(fpbit, outfp,  fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
															pulled_bitstream);
									else
										  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
															pulled_bitstream);
							}
						} else {
						  fseek( fpbit, major_mvBytes, SEEK_CUR );
						  in_bytes += major_mvBytes;
						  if ( info.AGP_level[j] )
						  {
							  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
												pulled_bitstream);
							  for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
								  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													pulled_bitstream);
						  }
						}
					}

					if ( !AGP_sign && layer_sign )   // no AGP, layer structure 
					{
						major_mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
						in_bytes += 4;
						if( j >= info.t_level ) { // temporal scalability
							process_mv_section(major_mvBytes, outfp, &out_bytes, &total_mv_bytes, &in_bytes, fpbit, 
								  fptab, pulled_bitstream);
							if ( info.layer_mv[j] )
							{
								if (info.layer_exist[j]) // enhancement bytes
									get_bytes_section_enhance(fpbit, outfp, fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
								else
									get_bytes_section_enhance(fpbit, outfp, fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
							}
						} else {
						  fseek( fpbit, major_mvBytes, SEEK_CUR );
						  in_bytes += major_mvBytes;
						  if ( info.layer_mv[j] )
						  {
							  get_bytes_section_enhance(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
												pulled_bitstream);
						  }
						}
					}

					if ( AGP_sign && layer_sign )  // AGP, layer structure 
					{
						major_mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
						in_bytes += 4;
						if( j >= info.t_level ) { // temporal scalability
							process_mv_section(major_mvBytes, outfp, &out_bytes, &total_mv_bytes, &in_bytes, fpbit, 
								  fptab, pulled_bitstream);
							if ( info.AGP_level[j] )  // AGP for base layer
							{
								if (info.AGP_exist[j]) // sign bytes
									get_bytes_section(fpbit, outfp, fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
								else
									get_bytes_section(fpbit, outfp, fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
								for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
									if ( sub_bit < info.AGP_exist[j])
										  get_bytes_section(fpbit, outfp,  fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
															pulled_bitstream);
									else
										  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
															pulled_bitstream);
							}
							if ( info.layer_mv[j] ) 
							{
								if (info.layer_exist[j]) // enhancement bytes
								{
									get_bytes_section_enhance(fpbit, outfp, fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
									if ( info.AGP_level[j] )  // AGP for enhancement layer
									{
										if (info.AGP_exist[j]) // sign bytes
											get_bytes_section(fpbit, outfp, fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
															  pulled_bitstream);
										else
											get_bytes_section(fpbit, outfp, fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
															  pulled_bitstream);
										for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
											if ( sub_bit < info.AGP_exist[j])
												  get_bytes_section(fpbit, outfp,  fptab, 1, &in_bytes, &out_bytes, &total_mv_bytes,
																	pulled_bitstream);
											else
												  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
																	pulled_bitstream);
									}
								}
								else
								{
									get_bytes_section_enhance(fpbit, outfp, fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													  pulled_bitstream);
									if ( info.AGP_level[j] )
									{
										get_bytes_section(fpbit, outfp, fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
														  pulled_bitstream);
										for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
											  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
																pulled_bitstream);

									}
								}
							}
						} else {
						  fseek( fpbit, major_mvBytes, SEEK_CUR );
						  in_bytes += major_mvBytes;
						  if ( info.AGP_level[j] )
						  {
							  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
												pulled_bitstream);
							  for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
								  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													pulled_bitstream);
						  }
						  if ( info.layer_mv[j] )
						  {
							  get_bytes_section_enhance(fpbit, outfp, fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													pulled_bitstream);
							  if ( info.AGP_level[j] )
							  {
								  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
													pulled_bitstream);
								  for (sub_bit=0; sub_bit<info.AGP_level[j]; sub_bit++)
									  get_bytes_section(fpbit, outfp,  fptab, 0, &in_bytes, &out_bytes, &total_mv_bytes,
														pulled_bitstream);
							  }
						  }
						}
					}



				 }
			  }
		  } // if (info.bi_mv[j])

		  // alternative Haar reconstruction 
		  if ( ! info.bi_exist[j] ) // not bi-directional motion field 
		  {
			  assert(0);
			  k = 0; 
			  while ( k <= GOPsz )
			  {
				  if ( k== 0)
				  {
					  // the first set motion vector is always not existing
					  assert( scene_change[j][k] == YES ); 
					  k++; 
				  }
				  else{
					  scene_change_left = scene_change[j][k];
					  k++;
					  scene_change_right = scene_change[j][k];
					  k++;
					  // we always need one set of motion vectors if there exist motion vectors 
					  if ( (scene_change_left==NO || scene_change_right==NO) )
						  if ( AGP_sign==0 )
						  {
							  mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
							  in_bytes += 4;
							  if( j >= info.t_level ) { // temporal scalability
								  process_mv_section(mvBytes, outfp, &out_bytes, &total_mv_bytes, &in_bytes, fpbit, 
									fptab, pulled_bitstream);
							  } else {
								fseek( fpbit, mvBytes, SEEK_CUR );
								in_bytes += mvBytes;
							  }
						  }else
						  {

						  }

					  if ( scene_change_left==NO && scene_change_right==NO )
						  if ( AGP_sign==0 )
						  {
							  mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
							  in_bytes += 4;
	  						  fseek( fpbit, mvBytes, SEEK_CUR );
							  in_bytes += mvBytes;
						  }else {

						  }
				  }
			  }  // while ( k <= GOPsz )
		  } // if ( ! info.bi_mv[j]) 

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
	
    curr += info.GOPsz;
    GOPcounter++;
  }                           // while

  assert((unsigned int)total_mv_bytes==sum_mv[info.t_level]); 
  assert((unsigned int)total_mv_bytes<=exist_sum_mv[info.t_level]); 

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
    free( scene_change[i] );
  }
  free( scene_change );

  /******************/
  /* write overhead */
  /******************/
  if( list_counter > 0 && pulled_bitstream == NO ){ // NOT WORKING
	assert(0);
    fseek( outfp, out_bytes, SEEK_SET );
    overhead_start = out_bytes;
    
    // write list 
    fwrite( &list_counter, sizeof( unsigned char ), 1, outfp );
    out_bytes += sizeof( unsigned char );
    for ( i = 0; i < list_counter; i++ ){
      fwrite( &list_rates[i], sizeof( unsigned short int ), 1, outfp );
      fwrite( &list_tlevel[i], sizeof( unsigned char ), 1, outfp );
      fwrite( &list_slevel[i], sizeof( unsigned char ), 1, outfp );
      fwrite( &list_qp[i], sizeof( unsigned char ), 1, outfp );
      fwrite( &list_simulcast[i], sizeof( unsigned char ), 1, outfp );
      fwrite( &list_out_bytes[i], sizeof( unsigned long int ), 1, outfp );
      out_bytes += sizeof( unsigned short int )+ 4 * sizeof( unsigned char )
        + sizeof( unsigned long int );
    }
    
    unsigned long input_size = 0, in_len;
    unsigned char *input_data;

    // load overhead file 
    if( !( fphead = fopen( overhead_stream , "rb" ) ) ) {
      printf( "can not open %s\n", overhead_stream );
      exit( 1 );
    }
    fseek( fphead, 0, SEEK_END );
    input_size = ftell( fphead );
    fseek( fphead, 0, SEEK_SET );
    
    input_data =
      ( unsigned char * )getarray( input_size, sizeof( unsigned char ),
                                   "input_data" );
    if( fread( input_data, sizeof( unsigned char ), input_size, fphead ) !=
        ( unsigned ) input_size ) {
      printf( "fread error \n" );
      exit( 1 );
    }
    fclose( fphead );
  
    // write overhead file to bitstream 
    in_len = input_size * sizeof( unsigned char );
    fwrite( &in_len, sizeof( unsigned long int ), 1, outfp );
    fwrite( input_data, sizeof( unsigned char ), in_len, outfp );
    out_bytes += sizeof( unsigned long int ) + in_len * sizeof( unsigned char );
    
    // write overhead position to bitstream   
    fwrite( &overhead_start, sizeof( unsigned long int ), 1, outfp );
    out_bytes += sizeof( unsigned long int );
     
    overhead += (out_bytes - overhead_start);
    free( input_data );
  }// NOT WORKING

  if ( pulled_bitstream == NO )
    assert ( overhead == estimated_overhead );
    
  printf( "---------------------------------------------------------------\n" );
  printf( " total overhead    : %ld bytes\n", overhead);
  printf( " + motion data     : %ld bytes out of existing motion data %ld bytes \n", 
			total_mv_bytes,  exist_sum_mv[info.t_level] );
  printf( " + subband data    : %ld bytes\n", total_subband_bytes );
  printf( " == total output   : %ld bytes (%ld bytes read)\n", 
          out_bytes, in_bytes );
  printf( "---------------------------------------------------------------\n" );
  
  if ( pulled_bitstream == NO ){
    fprintf( fptab, "%%---------------------------------------------------------------\n" );
    fprintf( fptab, "%ld \t %% == total_input\n", in_bytes );
    fprintf( fptab, "%%---------------------------------------------------------------\n" );
    fprintf( fptab, "%ld \t %% total_overhead \n", overhead );
    fprintf( fptab, "%ld \t %% + motion_data  \n", total_mv_bytes );
    fprintf( fptab, "%ld \t %% + subband_data \n", total_subband_bytes );
    fprintf( fptab, "%ld \t %% == total_output\n", out_bytes );
    fprintf( fptab, "%%---------------------------------------------------------------\n" );
    fclose( fptab );
  }
  
  fclose( fpbytes );
  fclose( fpbit );
  fclose( outfp );

  free( subband_fat );
  free( newsubband_fat );
  free( newfat );

  if(CBR == YES || VBR_ADAPT == YES){
	  for(i = 0; i < 10; i ++){
		  free(gop_sum_mv[i]);
		  free(gop_exist_sum_mv[i]);
	  }
  }
  
  return 0;
}

        
int
simulcast_pull( char *pre_stream, char *output_name, int t_level, 
                int s_level, int kbps, int qp ) 
{
  int i, j, k, nm, curr, last, num_of_GOP, GOPheader_bytes, FAT_bytes;
  int count = 0, remaining_frs, num_subband, bytes, mvBytes;
  int tmp, GOPsz, eff_GOPsz, is_pullable;
  int *subband_fat, *newsubband_fat;
  long int in_bytes = 0, out_bytes = 0, FAT_start = 0, *newfat; 
  long int overhead = 0, estimated_overhead = 0, total_mv_bytes = 0;
  long int subband_bytes = 0, total_subband_bytes = 0; 
  short int rel_s_level;
  
  unsigned char *data;
  char post_stream[256], seq_name[256], bit_alloc_name[256];
  int AGP_sign = 0, layer_sign=0; 
  
  unsigned long int sum_mv, exist_sum_mv;
                              
  enum FLAG CBR = NO, pulled_bitstream; 
  FILE *outfp, *fpbytes;
   
  videoinfo info;
  

  /**********/
  /* header */
  /**********/
  
  read_header( pre_stream, &info ); 
  in_bytes += sizeof( videoheader );

  pulled_bitstream = ( info.bitrate != 0 ) ? YES : NO;
  int display_operations = (pulled_bitstream == YES); // display pull operations
 
  info.bitrate = kbps;     
  info.t_level = t_level;
  info.s_level = s_level;
  rel_s_level  = s_level - (info.denoise_flag == YES);
  info.qp      = qp;
  
  strcpy( info.bitname, pre_stream );
  strncpy( seq_name, info.bitname, strlen( info.bitname ) - 4 );
  seq_name[strlen( info.bitname ) - 4] = '\0';

  num_of_GOP = get_GOP_num( info );

  // GOPheader_bytes (cmp. mvcoding.c)
  GOPsz = info.GOPsz;
  for( i = 0; i < info.tPyrLev; i++ ) {
    for( j = 0; j < GOPsz; j++ ) 
      count++;
    GOPsz /= 2;  
  }
  GOPsz = info.GOPsz; 
  GOPheader_bytes = (( count + 7 ) / 8) * num_of_GOP; 
  
  // FAT bytes
  FAT_bytes = num_of_GOP * sizeof( long int );

  if( pulled_bitstream == NO ){
  
    // read MV bytes
    sum_mv = read_mv_bytes( seq_name, info, AGP_sign, layer_sign,  &exist_sum_mv ); 
    
    // pull if target rate was found in list, then unzip data
    estimated_overhead += sizeof( videoheader ) + GOPheader_bytes + FAT_bytes;
    
    if( CBR == NO )
      is_pullable = bit_alloc_VBR( info, &info.qp, sum_mv, estimated_overhead );
    
    if( !is_pullable ){
      return -1;
    }
  }
  
  if( !( fpbit = fopen( pre_stream, "rb" ) ) ) { // total bitstream
    printf( "can not open: %s\n", pre_stream );
    exit( 1 );
  }
    
  if (pulled_bitstream == NO){
    sprintf( post_stream, "%s_%d.bit", seq_name, info.bitrate ); 
  } else {
    sprintf( post_stream, "%s_%d.bit", output_name, info.bitrate ); 
  }
  
  write_header( post_stream, info );
  out_bytes += sizeof( videoheader );
  
  if( !( outfp = fopen( post_stream, "r+b" ) ) ) {
    printf( "can not open: %s\n", post_stream );
    exit( 1 );
  }
    
  if( display_operations ){
    printf( "---------------------------------------------------------------\n" );
    printf( " Sequence name     : %s\n", seq_name );
    printf( " Sequence length   : %03d frames ( %03d - %03d )\n",
            info.last - info.start + 1, info.start, info.last );
    printf( " Output bit rate   : %d kbit/sec\n", info.bitrate );
    printf( " Output frame rate : %0.2f frames/sec (t_level %d)\n",
            ((double)(info.framerate) / pow((float)2, info.t_level)), info.t_level );
    printf( " Output frame size : (Y) %d x %d (C) %d x %d (s_level %d)\n",
            info.ywidth << rel_s_level, info.yheight << rel_s_level,
            info.cwidth << rel_s_level, info.cheight << rel_s_level,
          info.s_level );
    printf( " Number of GOPs    : %d\n", num_of_GOP );
    // printf( " Header size       : %d bytes\n", sizeof( videoheader ) );
    // printf( " GOP FAT bytes     : %d bytes\n", FAT_bytes ); 
    // printf( " GOP header bytes  : %d bytes\n", GOPheader_bytes );
    printf( " Total overhead    : %ld bytes\n", estimated_overhead );
#ifdef MY_SCALE
    printf( " MY_SCALE          : %d (%d extra bits)\n", MY_SCALE, (int)log2(MY_SCALE));
#endif
    printf( "---------------------------------------------------------------\n" );
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
  
  scene_change =
    ( enum FLAG ** )getarray( info.tPyrLev, sizeof( enum FLAG * ),
                              "scene_change" );
  GOPsz = info.GOPsz;
  for( i = 0; i < info.tPyrLev; i++ ) {
    scene_change[i] =
      ( enum FLAG * )getarray( GOPsz + 1, sizeof( enum FLAG ), "scene_change" );
    GOPsz /= 2;
  }
  
  curr = info.start;
  last = info.last;

  int GOPcounter = 0;
  if( pulled_bitstream == NO ){
    sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", seq_name, info.bitrate );
  } else {
    sprintf( bit_alloc_name, "%s_%d.bytes_per_GOP", output_name, info.bitrate );
  }
  
  if( CBR == NO ) {
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
    if( display_operations ){
      printf( " pulling GOP %d: frame %d - %d\n", GOPcounter,
              curr, curr + info.eff_GOPsz - 1 );
    }
    
    GOPheader_bytes = 0;
    if( info.tPyrLev >= 1 ) {

      /*************/
      /* GOPheader */
      /*************/

      GOPheader_bytes = read_GOPheader( scene_change, info );   //printf("%d\n", GOPheader_bytes);
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
      free( data );
      
      /******/
      /* MV */
      /******/
      
      GOPsz = 2;
      for( j = info.tPyrLev - 1; j >= 0; j-- ) {
        for( k = 0; k <= GOPsz; k++ ) {
          if( scene_change[j][k] == NO ) {
                
            mvBytes = read_number_core( fpbit ); // fpbit will move forward by 4 bytes 
            in_bytes += 4;
                
            if( j >= info.t_level ) { // temporal scalability
              write_number_core( mvBytes, outfp );
              out_bytes      += 4;
              total_mv_bytes += 4;
                  
              // printf("mvBytes %d\n", mvBytes);
              data =
                ( unsigned char * )getarray( mvBytes, sizeof( unsigned char ),
                                             "data" );
              if( fread( data, sizeof( unsigned char ), mvBytes, fpbit ) !=
                  ( unsigned )mvBytes ) {
                printf( "fread error2 %d\n", mvBytes );
                exit( 1 );
              }
              in_bytes += mvBytes;
              fwrite( data, sizeof( unsigned char ), mvBytes, outfp );
              out_bytes      += mvBytes;
              total_mv_bytes += mvBytes;
              free( data );
            } else {
              fseek( fpbit, mvBytes, SEEK_CUR );
              in_bytes += mvBytes;
            }
          }
        }
        GOPsz *= 2;
      }
    }   // if( info.tPyrLev >= 1 )

    /*************/
    /* Bitplanes */
    /*************/

    newfat[GOPcounter] = 0;
  
    eff_GOPsz = info.eff_GOPsz;
    if( info.denoise_flag == YES ) {
      eff_GOPsz *= ( YUV420 == 1 ) ? 2 : 4;
    }
    
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
        in_bytes += sizeof( unsigned char ) * subband_fat[nm];  //printf("%d\n", in_bytes);

        // write bit plane      
        fscanf( fpbytes, "%d\n", &tmp ); // Read "bytes_per_GOP" (bitalloc.c)
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
    }
        
    curr += info.GOPsz;
    GOPcounter++;
  }                           // while

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
    free( scene_change[i] );
  }
  free( scene_change );
  
  if( display_operations ){
    printf( "---------------------------------------------------------------\n" );
    printf( " total overhead    : %ld bytes\n", overhead );
    printf( " + motion data     : %ld bytes\n", total_mv_bytes );
    printf( " + subband data    : %ld bytes\n", total_subband_bytes );
    printf( " == total output   : %ld bytes (%ld bytes read)\n", 
            out_bytes, in_bytes );
    printf( "---------------------------------------------------------------\n" );
  }
  
  fclose( fpbytes );
  fclose( fpbit );
  fclose( outfp );

  free( subband_fat );
  free( newsubband_fat );

  return info.qp;
}