#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "structN.h"
#include "memoryN.h"

extern long int totalY[5], totalU[5], totalV[5], totalmap[5], totalMV[5];       /*initiated in initN.c */
extern float avgpsnr[3][16], avgvar[3][16];
extern long int *unconnectedL;

void read_frame( YUVimage_ptr frame, videoinfo info, char *inname, int index,
                 enum FORMAT format );

/****************************************************************************/
/*                                 snr_frame()                              */
/****************************************************************************/
void
snr_frame( float *ysnr, float *usnr, float *vsnr, YUVimage_ptr codeframe,
           YUVimage_ptr inframe, videoinfo info )
{
  int i, ypix, cpix;
  double sum, diff, peak;

  ypix = (info.ywidth * info.yheight); // >> (2 * info.s_level);
  cpix = (info.cwidth * info.cheight); // >> (2 * info.s_level);

  switch ( info.format ) {
  case YUV:
  case RAS:
    peak = 255.;
    break;
  case DPX:
    peak = 1023.;
    break;
  default:
    printf( "image format error format = %d(pstatN.c)\n", info.format );
    exit( 1 );
  }

  sum = 0.;
  for( i = 0; i < ypix; i++ ) {
    diff = inframe->Y[i] - codeframe->Y[i];
    sum += diff * diff;
  }
  *ysnr =
    ( ypix ) ? ( float )( 20. * log10( peak / sqrt( sum / ypix ) ) ) : 0;

  sum = 0.;
  for( i = 0; i < cpix; i++ ) {
    diff = inframe->U[i] - codeframe->U[i];
    sum += diff * diff;
  }
  *usnr =
    ( cpix ) ? ( float )( 20. * log10( peak / sqrt( sum / cpix ) ) ) : 0;

  sum = 0.;
  for( i = 0; i < cpix; i++ ) {
    diff = inframe->V[i] - codeframe->V[i];
    sum += diff * diff;
  }
  *vsnr =
    ( cpix ) ? ( float )( 20. * log10( peak / sqrt( sum / cpix ) ) ) : 0;
}


/****************************************************************************/
/*                                 snr_crop()                               */
/****************************************************************************/
void
snr_crop( float *ysnr, float *usnr, float *vsnr, YUVimage_ptr codeframe,
          YUVimage_ptr inframe, videoinfo info )
{
  int i, j, start_x, start_y, width, height, pos, ypix, cpix;
  double sum, diff, peak;

  start_x = 65;
  start_y = 0;
  width = 59;
  height = 240;
  ypix = width * height;
  cpix = ypix / 2;

  switch ( info.format ) {
  case YUV:
  case RAS:
    peak = 255.;
    break;
  case DPX:
    peak = 1023.;
    break;
  default:
    printf( "image format error format = %d(pstatN.c)\n", info.format );
    exit( 1 );
  }

  sum = 0.;
  for( i = 0; i < height; i++ ) {
    for( j = 0; j < width; j++ ) {
      pos = i * info.ywidth + start_x + j;
      diff = inframe->Y[pos] - codeframe->Y[pos];
      sum += diff * diff;
    }
  }

  *ysnr = ( float )( 20. * log10( peak / sqrt( sum / ypix ) ) );
  //*ysnr = sum;

  sum = 0.;
  for( i = 0; i < height / 2; i++ ) {
    for( j = 0; j < width / 2; j++ ) {
      pos = i * info.cwidth + start_x / 2 + j;
      diff = inframe->U[pos] - codeframe->U[pos];
      sum += diff * diff;
    }
  }
  *usnr = ( float )( 20. * log10( peak / sqrt( sum / cpix ) ) );
  //*usnr = sum;

  sum = 0.;
  for( i = 0; i < height / 2; i++ ) {
    for( j = 0; j < width / 2; j++ ) {
      pos = i * info.cwidth + start_x / 2 + j;
      diff = inframe->V[pos] - codeframe->V[pos];
      sum += diff * diff;
    }
  }
  *vsnr = ( float )( 20. * log10( peak / sqrt( sum / cpix ) ) );
  //*vsnr = sum;
}



/****************************************************************************/
/*                                calsnr()                                  */
/****************************************************************************/
float
calsnr( int start, int last, videoinfo info )
{
  int i, num; // j
  float ysnr, usnr, vsnr, mean1, mean2, mean3;
  YUVimage fr0, fr1;
  FILE *fpstat;


  //TR  info.coding_domain = LOG;     // calculate PSNR in LOG domain

  frame_alloc( &fr0, info );
  frame_alloc( &fr1, info );
  fpstat = fopen( info.statname, "at+" );
  fprintf( fpstat, "\n <psnr>\n" );
  fprintf( fpstat, " ysnr      usnr    vsnr\n" );
  printf( "\n" );
  printf( " calculate PSNR frame %d ~ frame %d\n", start, last );

  num   = 0;
  mean1 = 0.;
  mean2 = 0.;
  mean3 = 0.;
  for( i = start; i <= last; i = i + (1 << info.t_level) ) {

/*
    if ((strchr (info.decname, '%'))==0){ //check for frame or video mode 
      read_frame(&fr0, info, info.decname, (i - start) >> info.t_level, info.format); 
	}else
*/
    read_frame(&fr0, info, info.decname, i, info.format);  /* coded frame */   
    
    read_frame( &fr1, info, info.inname,  i, info.format );  /* original frame */

    snr_frame( &ysnr, &usnr, &vsnr, &fr0, &fr1, info );
    mean1 += ysnr;
    mean2 += usnr;
    mean3 += vsnr;  
    num++;
    fprintf( fpstat, "%.2f\t  %.2f\t  %.2f\n", ysnr, usnr, vsnr ); 
    printf( " frame %3d:  %.2f\t  %.2f\t  %.2f\n", i, ysnr, usnr, vsnr );
  }
 
  fprintf( fpstat, "==================================================\n" );
  printf( "==================================================\n" );
  // num = last - start + 1;
  fprintf( fpstat, "%6s(%03d) ysnr = %.2f usnr = %.2f vsnr = %.2f\n", "avg",
           num, mean1 / num, mean2 / num, mean3 / num );
  printf( " avg (%03d) ysnr = %.2f usnr = %.2f vsnr = %.2f\n\n", 
          num, mean1 / num, mean2 / num, mean3 / num );
  fclose( fpstat );

  free_frame( fr0 );
  free_frame( fr1 );

  mean1 = mean1 / num;

  return mean1;
}



/*
 * calsnr_seq()
 * calculate psnr of the whole sequence
 */
void
calsnr_seq( YUVimage_ptr fr0, YUVimage_ptr fr1, int start, int last,
            videoinfo info )
{
  int i, num;
  float ysnr, usnr, vsnr, mean1, mean2, mean3;
  FILE *fpstat;

  mean1 = 0.;
  mean2 = 0.;
  mean3 = 0.;
  fpstat = fopen( info.statname, "at+" );       /* this is the difference from calsnr */
  fprintf( fpstat, "\n <psnr>\n" );
  fprintf( fpstat, " ysnr      usnr    vsnr\n" );
  printf( "start %d last %d\n", start, last );

  for( i = start; i <= last; i++ ) {
    read_frame( fr1, info, info.inname, i, info.format );       /* original frame */
    read_frame( fr0, info, info.decname, i, info.format );      /* coded frame */
    snr_frame( &ysnr, &usnr, &vsnr, fr0, fr1, info );
    mean1 += ysnr;
    mean2 += usnr;
    mean3 += vsnr;

    fprintf( fpstat, "%.2f\t  %.2f\t  %.2f\n", ysnr, usnr, vsnr );
  }

  fprintf( fpstat, "=================================================\n" );

  num = last - start + 1;

  fprintf( fpstat, "%6s(%03d) ysnr = %.2f usnr = %.2f vsnr = %.2f\n", "avg",
           num, mean1 / num, mean2 / num, mean3 / num );
  fclose( fpstat );

}

/*
 * calsnr_frame()
 * calculate psnr of a frame
 */

void
calsnr_frame( YUVimage_ptr fr0, YUVimage_ptr fr1, int curr, videoinfo info )
{
  int num;
  float ysnr, usnr, vsnr;
  static float mean1 = 0., mean2 = 0., mean3 = 0.;
  FILE *fpstat;

  //mean1=0.; mean2=0.; mean3=0.;
  fpstat = fopen( info.statname, "at+" );       /* this is the difference from calsnr */
  if( curr == info.start ) {
    fprintf( fpstat, "\n <psnr>\n" );
    fprintf( fpstat, " ysnr      usnr    vsnr\n" );
  }

  snr_frame( &ysnr, &usnr, &vsnr, fr0, fr1, info );

  mean1 += ysnr;
  mean2 += usnr;
  mean3 += vsnr;

  fprintf( fpstat, "%.2f\t  %.2f\t  %.2f\n", ysnr, usnr, vsnr );

  if( curr == info.last ) {
    fprintf( fpstat, "=================================================\n" );
    num = info.last - info.start + 1;
    fprintf( fpstat, "%6s(%03d) ysnr = %.2f usnr = %.2f vsnr = %.2f\n", "avg",
             num, mean1 / num, mean2 / num, mean3 / num );
  }
  fclose( fpstat );
}


void
print_mvbits( videoinfo info, Rate FrsRate )
{
  int i, j, nlev, nfrs, count, itmp, sum;
  FILE *fp_mv;

  if( !( fp_mv = fopen( info.mvstatname, "at+" ) ) ) {
    printf( "can not open %s\n", info.mvstatname );
    exit( 1 );
  }
 
  nlev = info.tPyrLev;
  
  nfrs = 2;
  count = 1;
  sum = 0; 
  itmp = 0; 
  for( i = nlev - 1; i >= 0; i-- ) {
    for( j = 0; j <= nfrs; count++, j++ ) 
    {
//      itmp = (FrsRate.map[count] + FrsRate.mv[count]);
      sum += itmp;
      printf  ("Fr%.2d: mvbits sum = %5d\n", count, itmp);
      fprintf (fp_mv, "%d\n", itmp >> 3);
    }
    nfrs *= 2;
  }
  fclose( fp_mv );
  printf( "---> %.2d MV sets -- total sum = %7d\n\n", count, sum);
}
