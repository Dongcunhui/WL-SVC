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
#include "ezbc_enc_3d.h"


long int
intra_encode( int curr, int GOP_counter, videoinfo info )
{
  int i;
  long int output_GOP_bytes;

  for( i = 0; i < info.GOPsz; i++ ) {
    read_frame( &pyrFrs[i], info, info.inname, curr + i, info.format ); /*ioN.c */
  }
  output_GOP_bytes = ezbc3d_enc_GOP( curr, pyrFrs, info, GOP_counter );
  return output_GOP_bytes;

}
