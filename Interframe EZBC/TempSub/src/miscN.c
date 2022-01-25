#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "structN.h"
#include "basic.h"
#include <iostream>

#define CLOCKS   1e6            /* for CPU time */

long clock(  );



/*
 *    variance()   
 *  calculate variance of a sub_block inside an frame.
 *  frame : a pointer pointing the upper-left corner of this block
 *  lenx  : width of sub_block.
 *  leny  : height of sub_block.
 *  hor   : width of the original frame
 */
float
variance( float *frame, int lenx, int leny, int hor )
{
  int x, y;                     //, pos
  double mean, mean2;           //data, 
  float *iptr;

  if( !( lenx * leny ) )
    return 0.;
  iptr = frame;

  mean = mean2 = 0.;
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      mean += *iptr;
      mean2 += ( *iptr ) * ( *iptr );
      iptr++;
    }
    iptr += hor - lenx;
  }
  mean /= lenx * leny;
  mean2 /= lenx * leny;

  mean2 -= mean * mean;

  return ( ( float )mean2 );

}

/*
 *    get_mean()   
 *  calculate mean of a sub_block inside a frame.
 *  frame : a pointer pointing the upper-left corner of this block
 *  lenx  : width of sub_block.
 *  leny  : height of sub_block.
 *  hor   : width of the original frame
 */

float
get_mean( float *frame, int lenx, int leny, int hor )
{
  int x, y;                     //, pos
  double mean;
  float *iptr;

  if( !( lenx * leny ) )
    return 0.;
  iptr = frame;

  mean = 0.;
  for( y = 0; y < leny; y++ ) {
    for( x = 0; x < lenx; x++ ) {
      mean += *iptr;
      iptr++;
    }
    iptr += hor - lenx;
  }
  mean /= lenx * leny;
  return ( ( float )mean );
}


double
variance_mean( float *frame, int lenx, int leny, float *mean )
{
  int x, y, pos;
  double data, sum, sum2;

  if( !( lenx * leny ) )
    return 0.;
  sum = sum2 = 0.;
  for( y = 0; y < leny; y++ ) {
    pos = y * lenx;
    for( x = 0; x < lenx; x++ ) {
      data = frame[pos + x];
      sum += data;
      sum2 += data * data;
    }
  }
  *mean = ( float )( sum / ( lenx * leny ) );
  sum2 /= lenx * leny;
  sum2 -= ( *mean ) * ( *mean );

  return ( sum2 );
}


void
varianceGOP( YUVimage_ptr frames, videoinfo info, Rate_ptr alloc )
{
  int i, yhor, yver, chor, cver;

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;

  for( i = 0; i < info.GOPsz; i++ ) {
    alloc->yvar[i] = variance( frames[i].Y, yhor, yver, yhor );
    alloc->uvar[i] = variance( frames[i].U, chor, cver, chor );
    alloc->vvar[i] = variance( frames[i].V, chor, cver, chor );
  }

  return;
}


long int
get_GOP_num( videoinfo info )
{
  long int num;

  num = (long int) ( (info.last - info.start + info.GOPsz) / info.GOPsz );

  return num;
}
/*****************************************************************************/
/*                              header2info()                                */
/*****************************************************************************/
void
header2info( videoinfo_ptr info, videoheader header )
{
  int hor, ver, xblk, yblk, i;

  info->framerate = header.framerate;
  info->format = YUV;

  info->ywidth = header.ywidth;
  info->yheight = header.yheight;
  info->cwidth = header.cwidth;
  info->cheight = header.cheight;
  info->start = header.start;
  info->last = header.last;

  hor  = info->ywidth;
  ver  = info->yheight;
  info->maxMBnum = 0;
  for (i = 0; i < MAX_TLEVELS; i++) {
    info->level[i]        = header.level[i] & 0x0f ;		   // bit: 0-3 for level
	info->AGP_exist[i]    = ( header.level[i] & 0xf0 )>>4;     // bit: 4-7 for AGP_exist

    info->subpel[i]       =  header.subpel[i] & 0x03;          // bit: 0-1 for subpel
	info->AGP_level[i]    = (header.subpel[i] & 0x0C )>>2;     // bit: 2-3 for AGP_level
	info->bi_mv[i]        = (header.subpel[i] & 0x10 )>>4;     // bit: 4   for bi_mv
	info->bi_exist[i]     = (header.subpel[i] & 0x20 )>>5;     // bit: 5   for bi_exist
	info->layer_mv[i]     = (header.subpel[i] & 0x40 )>>6;     // bit: 6   for layer_mv
	info->layer_exist[i]  = (header.subpel[i] & 0x80 )>>7;     // bit: 7   for layer_exist
    info->xblk[i] = header.xblk[i];
    info->yblk[i] = header.yblk[i];
    xblk = info->xblk[i];
    yblk = info->yblk[i];
    info->xnum[i] = ( !( hor % xblk ) ) ? hor / xblk : hor / xblk + 1;
    info->ynum[i] = ( !( ver % yblk ) ) ? ver / yblk : ver / yblk + 1;

    if(info->xnum[i] * info->ynum[i] > info->maxMBnum) {
      info->maxMBnum = info->xnum[i] * info->ynum[i];
    }
  }

  info->bitrate = header.bitrate;

  info->t_level = header.t_level;
  info->s_level = header.s_level;
  info->tPyrLev = header.tPyrLev;
  info->GOPsz = 0x1 << header.tPyrLev;
  info->bigGOP = ( 0x1 << ( header.tPyrLev + 1 ) ) - 1;
  info->denoise_flag = header.denoise_flag;

  info->SLTF_range = header.SLTF_range;
  info->SHTF_range = header.SHTF_range;

}

/*****************************************************************************/
/*                              info2header()                                */
/*****************************************************************************/
void
info2header( videoheader_ptr header, videoinfo info )
{
  int i;

  header->framerate = ( short int )info.framerate;

  header->ywidth = ( short int )info.ywidth;
  header->yheight = ( short int )info.yheight;
  header->cwidth = ( short int )info.cwidth;
  header->cheight = ( short int )info.cheight;
  header->start = ( short int )info.start;
  header->last = ( short int )info.last;

  for (i = 0; i < MAX_TLEVELS; i++) {
    header->level[i] = ( short int ) (  info.level[i] +             // bit: 0-3   for level
		                                (info.AGP_exist[i]<<4) );   // bit: 4-7   for AGP_exist
    header->subpel[i] = ( short int )( (info.subpel[i]<<0)+         // bit: 0-1   for subpel
									   (info.AGP_level[i]<<2) +     // bit: 2-3   for AGP_level
 									   (info.bi_mv[i]<<4)+          // bit: 4     for bi_mv
									   (info.bi_exist[i]<<5)+       // bit: 5     for bi_exist
									   (info.layer_mv[i]<<6)+       // bit: 6     for layer_mv
									   (info.layer_exist[i]<<7));   // bit: 7     for layer_exist
    header->xblk[i] = ( short int )info.xblk[i];
    header->yblk[i] = ( short int )info.yblk[i];
  }

  header->bitrate = info.bitrate;

  header->t_level = info.t_level ;   


  header->s_level = info.s_level;
  header->tPyrLev = ( unsigned char )info.tPyrLev;
  header->denoise_flag = info.denoise_flag;

  header->SLTF_range = info.SLTF_range;
  header->SHTF_range = info.SHTF_range;

}

/*****************************************************************************/
/*                                read_header()                              */
/*****************************************************************************/
void
read_header( char *name, videoinfo * info )
{
  FILE *fpio;
  videoheader header;

  if( !( fpio = fopen( name, "rb" ) ) ) {
    printf( "read_header: %s\n", name );
    exit( 1 );
  }
  fread( &header, sizeof( videoheader ), 1, fpio );
  fclose( fpio );
  header2info( info, header );// 头信息转换为info

}

/*****************************************************************************/
/*                               write_header()                              */
/*****************************************************************************/
void
write_header( char *name, videoinfo info )
{
  FILE *fpio;
  videoheader header;

  info2header( &header, info );

  if( !( fpio = fopen( name, "wb" ) ) ) {
    printf( "write_header: %s\n", name );
    exit( 1 );
  }
  fwrite( &header, sizeof( videoheader ), 1, fpio );
  fclose( fpio );
}

/*****************************************************************************/
/*                                get_mvBytes()                              */
/*****************************************************************************/

int
get_mvBytes( FILE * fp_mv )
{
  int bytes;
  char iline[80];

  fgets( iline, 78, fp_mv );
  sscanf( iline, "%d", &bytes );
  return bytes;
}

/*****************************************************************************/
/*                                 timecheck()                               */
/*****************************************************************************/
void
timecheck( long *pmark, int mode )
{
  double sec;
  long cmark;

  cmark = clock(  );            /* current clock */
  sec = ( cmark - *pmark ) / CLOCKS;
  *pmark = cmark;

  switch ( mode ) {
  default:
    printf( "error in timecheck()\n" );
    exit( 1 );
  case 0:
    printf( "ME   = %5.2f seconds.\n", sec );
    break;
  case 1:
    printf( "MCTF = %5.2f seconds.\n", sec );
    break;
  case 2:
    printf( "WAVE = %5.2f seconds.\n", sec );
    break;
  case 3:
    printf( "FSSQ = %5.2f seconds.\n", sec );
    break;
  }
}

/*****************************************************************************/
/*                                 showtime()                                */
/*****************************************************************************/
void
showtime( long pmark )
{
  int hour, min;
  double sec;
  long cmark;

  cmark = clock(  );            /* current clock */
  sec = ( cmark - pmark ) / CLOCKS;
  hour = ( int )( sec / 3600.0 );
  sec -= 3600.0 * hour;
  min = ( int )( sec / 60.0 );
  sec -= 60.0 * min;

  printf( "%d hour %d min %5.2f seconds.\n", hour, min, sec );
}


void
computegain( YUVimage lowband, YUVimage highband, videoinfo info )
{
  int yhor, yver, chor, cver;
  double var0, var1, gain;

  yhor = info.ywidth;
  yver = info.yheight;
  chor = info.cwidth;
  cver = info.cheight;

  var0 = variance( lowband.Y, yhor, yver, yhor );
  var1 = variance( highband.Y, yhor, yver, yhor );
  gain = ( var0 + var1 ) / 2.0;
  gain /= sqrt( var0 * var1 );
  fprintf( stdout, "lowband yvar = %f highband yvar = %f coding gain = %f\n",
           ( float )var0, ( float )var1, ( float )gain );

  if( chor ) {
    var0 = variance( lowband.U, chor, cver, chor );
    var1 = variance( highband.U, chor, cver, chor );
    gain = ( var0 + var1 ) / 2.0;
    gain /= sqrt( var0 * var1 );
    fprintf( stdout,
             "lowband uvar = %f highband uvar = %f coding gain = %f\n",
             ( float )var0, ( float )var1, ( float )gain );

    var0 = variance( lowband.V, chor, cver, chor );
    var1 = variance( highband.V, chor, cver, chor );
    gain = ( var0 + var1 ) / 2.0;
    gain /= sqrt( var0 * var1 );
    fprintf( stdout,
             "lowband vvar = %f highband vvar = %f coding gain = %f\n",
             ( float )var0, ( float )var1, ( float )gain );
  }

}


void
copyframe( YUVimage * source, YUVimage * dest, videoinfo info )
{
  int i;

  for( i = 0; i < info.ywidth * info.yheight; i++ ) {
    dest->Y[i] = source->Y[i];
  }
  if( info.cwidth && info.cheight ) {
    for( i = 0; i < info.cwidth * info.cheight; i++ ) {
      dest->U[i] = source->U[i];
      dest->V[i] = source->V[i];
    }
  }

}


void
wcopyframe( YUVimage * source, YUVimage * dest, float weight, videoinfo info )
{
  int i;

  for( i = 0; i < info.ywidth * info.yheight; i++ ) {
    dest->Y[i] = source->Y[i] * weight;
  }
  if( info.cwidth && info.cheight ) {
    for( i = 0; i < info.cwidth * info.cheight; i++ ) {
      dest->U[i] = source->U[i] * weight;
      dest->V[i] = source->V[i] * weight;
    }
  }

}


void 
child_mv_copy(vector_ptr source, vector_ptr dest, videoinfo info)
{
  int i;

  assert( dest->child == 0 );
  
  if( source->child ) {
    dest->child = 1;

    // create empty children
    dest->child0 = ( vector_ptr ) getarray( 1, sizeof( vector ), "dest->child0" );
    dest->child1 = ( vector_ptr ) getarray( 1, sizeof( vector ), "dest->child1" );
    dest->child2 = ( vector_ptr ) getarray( 1, sizeof( vector ), "dest->child2" );
    dest->child3 = ( vector_ptr ) getarray( 1, sizeof( vector ), "dest->child3" );
   
    dest->child0->child = 0;
    dest->child1->child = 0;
    dest->child2->child = 0;
    dest->child3->child = 0;
   
    dest->child0->parent = dest;
    dest->child1->parent = dest;
    dest->child2->parent = dest;
    dest->child3->parent = dest;
    
    // go through the tree
    child_mv_copy( source->child0, dest->child0, info );
    child_mv_copy( source->child1, dest->child1, info ); 
    child_mv_copy( source->child2, dest->child2, info );
    child_mv_copy( source->child3, dest->child3, info );
  } else {
    // fill nodes
    dest->child = source->child;
    dest->bi_mode = source->bi_mode;
    dest->lifting_mode = source->lifting_mode;
//    dest->Ymean = source->Ymean;
//    dest->Umean = source->Umean;
//    dest->Vmean = source->Vmean;
    dest->meandepth = source->meandepth;
    dest->mvx = source->mvx;
    dest->mvy = source->mvy;
    dest->mad = source->mad;
//    dest->bi_mvx = source->bi_mvx;
//    dest->bi_mvy = source->bi_mvy;
//    dest->bi_mad = source->bi_mad; 
    dest->is_predictor = source->is_predictor;
    dest->sad_cost = source->sad_cost; 
    dest->bit_cost = source->bit_cost;
    dest->total_cost = source->total_cost;
    dest->mse = source->mse;
    dest->var = source->var;
//    dest->dL = source->dL;                       
//    dest->dMAP = source->dMAP;                   
//    dest->dMV = source->dMV;                      
//    dest->dD = source->dD;                       
//    dest->slope = source->slope;                  
//    dest->mslope = source->mslope;
    dest->merge = source->merge;

	dest->med_idx = source->med_idx;	//Added by Yuan Liu on 01.23.2016
	dest->aff_mrg = source->aff_mrg;

	dest->two_comp_src = source->two_comp_src;
	dest->skip_sign = source->skip_sign;  //Added on 08.12.2018
/*
	for(i=0;i<=3;i++){
		dest->mrg_mvx[i] = source->mrg_mvx[i];
		dest->mrg_mvy[i] = source->mrg_mvy[i];

		dest->mrg_aff_mvx1[i] = source->mrg_aff_mvx1[i];
		dest->mrg_aff_mvy1[i] = source->mrg_aff_mvy1[i];
		dest->mrg_aff_mvx2[i] = source->mrg_aff_mvx2[i];
		dest->mrg_aff_mvy2[i] = source->mrg_aff_mvy2[i];
		dest->mrg_aff_mvx3[i] = source->mrg_aff_mvx3[i];
		dest->mrg_aff_mvy3[i] = source->mrg_aff_mvy3[i];
	}
*/
	for(i=0;i<3;i++){
		dest->aff1_pred_mvx[i] = source->aff1_pred_mvx[i];
		dest->aff1_pred_mvy[i] = source->aff1_pred_mvy[i];
		dest->aff2_pred_mvx[i] = source->aff2_pred_mvx[i];
		dest->aff2_pred_mvy[i] = source->aff2_pred_mvy[i];
		dest->aff3_pred_mvx[i] = source->aff3_pred_mvx[i];
		dest->aff3_pred_mvy[i] = source->aff3_pred_mvy[i];
	}
	////////	Added by Yuan Liu	///////////
	if( (source->bi_mode >= 9 && source->bi_mode <=11) || (source->bi_mode == 7 && source->aff_mrg == YES) ){

		dest->direct_idx = source->direct_idx;
		dest->aff_idx = source->aff_idx;
//		printf("\naff1x = %f, aff1y = %f \n aff2x = %f, aff2y = %f \n aff3x = %f, aff3y = %f\n",source->aff1_mvx,source->aff1_mvy,source->aff2_mvx,source->aff2_mvy,source->aff3_mvx,source->aff3_mvy);

		dest->aff1_mvx = source->aff1_mvx;
		dest->aff1_mvy = source->aff1_mvy;
		dest->aff2_mvx = source->aff2_mvx;
		dest->aff2_mvy = source->aff2_mvy;
		dest->aff3_mvx = source->aff3_mvx;
		dest->aff3_mvy = source->aff3_mvy;

		dest->merge_idx = source->merge_idx;
		dest->merge_dir = source->merge_dir;
		dest->trans_pred_idx = source->trans_pred_idx;   //Added on 01.16.2017

		dest->aff1_dmvx = source->aff1_dmvx;
		dest->aff1_dmvy = source->aff1_dmvy;
		dest->aff2_dmvx = source->aff2_dmvx;
		dest->aff2_dmvy = source->aff2_dmvy;
		dest->aff3_dmvx = source->aff3_dmvx;
		dest->aff3_dmvy = source->aff3_dmvy;

		dest->dmvx = source->dmvx;
		dest->dmvy = source->dmvy;

	}
	///////////////////////////////////////////

	dest->iblock_spatial_mode  = source->iblock_spatial_mode;   // for spatial prediction mode 

    for (i = 0; i < NUMBER_OF_BI_MODES; i++) {
      dest->mode_info[i] = source->mode_info[i];
    }
  }
}

void 
mv_copy( vector_ptr source, vector_ptr dest, videoinfo info )
{
  int i, size;
  
  size = info.maxMBnum;
  for( i = 0; i < size; i++ ) {
    child_mv_copy( &source[i], &dest[i], info );  
  }
}

