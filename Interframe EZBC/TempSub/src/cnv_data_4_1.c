/*   file: cnv_data_4_1.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "basic.h"

void
cnv_data_4_1( float *y, unsigned char *x, int n )
{
  int i, tmp;

  for( i = 0; i < n; i++ ) {
    tmp = nint( y[i] + 128.0 );
    if( tmp < 0 )
      tmp = 0;
    if( tmp > 255 )
      tmp = 255;
    x[i] = ( unsigned char )tmp;
  }
  return;
}
