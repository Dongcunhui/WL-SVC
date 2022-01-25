/* ========================================================================= */
/* Description: video uitils for i/o                                         */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <stdio.h>
#include "structN.h"
#include "general.h"
#include "video_utils.h"



/*****************************************************************************/
/*                         alloc_pyr_frames()                                */
/*****************************************************************************/
int
alloc_pyr_frames( YUVimage * pyrFrs, videoinfo & info )
{
  int i;
  int sz  = info.ywidth * info.yheight;
  int csz = info.cwidth * info.cheight;

  for( i = info.GOPsz - 1; i >= 0; i-- ) {
    NEW_VECTOR( pyrFrs[i].Y, sz, float,  "Y" );
    NEW_VECTOR( pyrFrs[i].U, csz, float, "U" );
    NEW_VECTOR( pyrFrs[i].V, csz, float, "V" );
  }
  return ( 0 );
}

void
free_frames( YUVimage * pyrFrs, int n )
{
  for( int i = n - 1; i >= 0; i-- ) {
    free( pyrFrs[i].Y );        //DELETE_VECTOR(pyrFrs[i].Y);
    free( pyrFrs[i].U );        //DELETE_VECTOR(pyrFrs[i].U);
    free( pyrFrs[i].V );        //DELETE_VECTOR(pyrFrs[i].V);
  }
}
