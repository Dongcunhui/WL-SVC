#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define EXTERN  extern
#include "basic.h"
#include "rasterfile.h"
#include "structN.h"
#include "coderN.h"
#include "analsyn.h"
#include "ioN.h"
#include "miscN.h"
#include "memoryN.h"
#include "chrom.h"
#include "mvcodingN.h"
#include "ezbc_dec_3d.h"
#include "bmeN.h"
#include "pstatN.h"


void ezbc3d_dec_GOP( YUVimage * pyrFrs, videoinfo info, long total_bytes_past,
                     long int GOP_counter, int curr );


void
intra_decode( int curr, int GOP_counter, long int *total_bytes_past,
              videoinfo info )
{
  int i;

  ezbc3d_dec_GOP( dec_pyrFrs, info, *total_bytes_past, GOP_counter, curr );

  for( i = 0; i < info.GOPsz; i++ ) {
    write_frame( dec_pyrFrs[i], info, info.decname, curr + i, info.format );
  }

  /*
  for( i = 0; i < info.GOPsz; i++ ) {
    free_frame( dec_pyrFrs[i] );
  }
  */
}
