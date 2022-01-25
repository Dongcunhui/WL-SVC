#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#define EXTERN extern
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

void write_number_core( int outputbyte, FILE * fp );
int read_GOPheader( enum FLAG **dec_scene_change, videoinfo info );
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
			  char *pre_stream, enum FLAG *CBR, enum FLAG *list, 
              char *list_name)
{
  int i, j, k, argnum = 1;

  *CBR = NO;
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

  // determine effective GOP size in level j
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
  
  while( curr <= last ) {     
    remaining_frs = last - curr + 1; 
    
    for( j = info.tPyrLev - 1; j >= 0; j-- ) {

	  if ( info.bi_exist[j] )  // bi-directional motion field 
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