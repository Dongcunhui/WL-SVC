/* ========================================================================= */
/* Description: function implementation of wavelet decomposotion             */
/* Author: Shih-Ta Hsiang                                                    */
/* Version: v0.a                                                             */
/* Last Revised: Aug. 15, 2000                                               */
/* ========================================================================= */

#include <iostream>
#include <string.h>
#include <assert.h>
#include "subband.h"

#include "structN.h"


int
  SubbandLayer::lsb = - EXTRA_BIT;

int
Subband::

initialization( Image_Coord dimen, Image_Coord memory_org,
                Image_Coord memory_dim )
{

  dim = dimen;
  mem_org = memory_org;
  mem_dim = memory_dim;
  alloc_coeff(  );

  transpose_flag = 0;

  return 0;
}


Subband::Subband( Image_Coord dimen ):coeff( NULL )
{
  initialization( dimen );
}


Subband::Subband( Image_Coord dimen, Image_Coord org_ord, Image_Coord mem_dimen ):coeff
  ( NULL )
{
  initialization( dimen, org_ord, mem_dimen );
}

void
Subband::reset_coeff(  )
{
  memset( &( coeff[mem_org.x][mem_org.y] ), 0,
          sizeof( SUB_COEFF_TYPE ) * mem_dim.x * mem_dim.y );
}

void
Subband::alloc_coeff( void )
{

  free_coeff(  );

  NEW_VECTOR( coeff, mem_dim.x, SUB_COEFF_TYPE *, "matrix_row" );

  SUB_COEFF_TYPE **dp = coeff;

  coeff -= mem_org.x;

  int nn = mem_dim.x * mem_dim.y;

  NEW_VECTOR( *dp, nn, SUB_COEFF_TYPE, "matrix_row" );
  memset( *dp, 0, nn * sizeof( SUB_COEFF_TYPE ) );

  SUB_COEFF_TYPE *sp = *dp -= mem_org.y;

  for( int i = mem_dim.x - 1; i > 0; i--, *( ++dp ) = ( sp += mem_dim.y ) );

}

void
Subband::free_coeff(  )
{
  if( coeff ) {
    delete[] & ( coeff[mem_org.x][mem_org.y] );
    delete[] & ( coeff[mem_org.x] );
    coeff = NULL;
  }
}
int
Subband::get_band_level(  )
{
  assert( wp );
  if( band_idx < 4 )
    return wp->get_pyr_levels(  );
  return wp->get_pyr_levels(  ) - ( band_idx - 1 ) / 3;
}

//-----------------------------------------------------------------------------

//  get_child_band_idx()

//    ATTENTION to band 0 and bottom bands
//-----------------------------------------------------------------------------
int
Subband::get_child_band_idx( void )
{

  if( get_band_level(  ) > 1 )
    return band_idx + 3;
  else
    return NO_CHILD;
}

//-----------------------------------------------------------------------------

//  get_child_band_idx()

//    ATTENTION to band 0 and bottom bands
//-----------------------------------------------------------------------------
int
Subband::get_grand_child_band_idx( void )
{

  if( get_band_level(  ) > 2 )
    return band_idx + 6;
  else
    return NO_CHILD;
}

void
Subband::transpose( void )
{

  transpose_flag = !transpose_flag;
  Image_Coord old_dim = dim;
  Image_Coord old_mem_org = mem_org;
  Image_Coord old_mem_dim = mem_dim;
  SUB_COEFF_TYPE **old_coeff = coeff;

  dim = old_dim.t(  );
  mem_org = old_mem_org.t(  );
  mem_dim = old_mem_dim.t(  );
  coeff = NULL;
  alloc_coeff(  );

  for( int i = old_mem_dim.x + old_mem_org.x - 1; i >= old_mem_org.x; i-- )
    for( int j = old_mem_dim.y + old_mem_org.y - 1; j >= old_mem_org.y; j-- )
      coeff[j][i] = old_coeff[i][j];

  delete[] & ( old_coeff[old_mem_org.x][old_mem_org.y] );
  delete[] & ( old_coeff[old_mem_org.x] );
}

/*-----------------------------------------------------------------

The bitplane corresponds to the integral parts of the DWT coeffs.
When floating parts are needed (lsb < 0), the bitplane is shifted
to the left according to the digits of floating number needed.

------------------------------------------------------------------*/

SubbandLayer::
SubbandLayer( Image_BW & image_pyr, Image_Coord orig, Image_Coord dimen,
              int par_n )
  :
Subband( dimen, -1, dimen + 2 ),
org( orig ),
par_ind( par_n )
{

  int i, j, m, n;


  int max = 0;

  if( lsb < 0 ) {
    int scale = 1 << ( -lsb );
    for( i = 0, m = org.x; i < dim.x; i++, m++ )
      for( j = 0, n = org.y; j < dim.y; j++, n++ ) {
        coeff[i][j] = int ( ABS( image_pyr( m, n ) * scale ) );
        if( max < coeff[i][j] )
          max = coeff[i][j];
        if( image_pyr( m, n ) < 0 )
          coeff[i][j] |= sign;
      }
  } else {
    for( i = 0, m = org.x; i < dim.x; i++, m++ )
      for( j = 0, n = org.y; j < dim.y; j++, n++ ) {
        coeff[i][j] = int ( ABS( image_pyr( m, n ) ) );
        if( max < coeff[i][j] )
          max = coeff[i][j];
        if( image_pyr( m, n ) < 0 )
          coeff[i][j] |= sign;
      }
  }

  max_msb = log2( max );
  if( lsb < 0 ){
    max_msb += lsb;
    // printf("lsb %d, max %d\n ", lsb, max / (1<<(-lsb)));
  }

}

/*------------------------------------------------------------

  initialization for decoding
_____________________________________________________________*/


SubbandLayer::
SubbandLayer( Image_Coord orig, Image_Coord dimen, int par_n, int msb )
  :
Subband( dimen, -1, dimen + 2 ),
org( orig ),
par_ind( par_n ),
max_msb( msb )
{
//  int i, j, m, n;

  reset_coeff(  );

}

enum ORIENTATION
SubbandLayer::get_orientation(  )
{
  if( band_idx )
    return ORIENTATION( ( band_idx - 1 ) % 3 + 1 );
  else
    return LL;
}


SubSetLayer::SubSetLayer( Image_BW & image_pyr, int s_level, int *msb_layers )
:SubSet(  )
{
  int i;
  Image_Coord dim, org, small_dim;


  pyr_levels = image_pyr.pyramid_levels(  );
  pdim = image_pyr.pyramid_dim(  );
  if( pyr_levels < 0 ) {
    fprintf( stderr, "SubbandSet: Input Image object is empty\n" );
    exit( -1 );
  } else if( pyr_levels == 0 ) {
    image_pyr.transform( 0 );   // smoothing = 0
    pyr_levels = image_pyr.pyramid_levels(  );
  }

#ifdef ROLL_STRUCTURE_TWO
  if (s_level<2)
  {
	nband = 1 + 3 * pyr_levels+9;  // there are additional 9 subbands for structure two 
  }else
	nband = 1 + 3 * pyr_levels;

#else
  nband = 1 + 3 * pyr_levels; // 子带个数
#endif


  NEW_VECTOR( bands, nband, SUBBAND_TYPE *, "bands*" ); // 为每个子带申请位置

  dim = pdim;
  dim.x >>= pyr_levels;
  dim.y >>= pyr_levels;
  org.x = org.y = 0;

  //parent index

  int par_n = -1;
  bands[0] = ( !msb_layers ) ? new SUBBAND_TYPE( image_pyr, org, dim, par_n )
    : new SUBBAND_TYPE( org, dim, par_n, msb_layers[0] );
  assert( bands[0] );
  bands[0]->set_wp( this );
  bands[0]->set_band_idx( 0 );

  int cur = 1;
  for( i = 0; i < pyr_levels; i++ ) {
#ifdef  ROLL_STRUCTURE_TWO
	int m, n; 
	if (i!= pyr_levels-2+s_level ) // pyr_levels-2 in encoder, pyr_levels-2+s_level in decoder 
	{
		org.x = 0; org.y = dim.y;
		par_n = i ? cur-3 : 0;
		bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, dim, par_n)
								  : new SUBBAND_TYPE(org, dim, par_n, msb_layers[cur]);
		assert(bands[cur]);
		bands[cur]->set_wp(this);  
		bands[cur]->set_band_idx(cur);
		cur++;

		org.x = dim.x; org.y = 0;
		par_n = i ? cur-3 : 0;
		bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, dim, par_n)
								  : new SUBBAND_TYPE(org, dim, par_n, msb_layers[cur]);
		assert(bands[cur]);
		bands[cur]->set_wp(this);  
		bands[cur]->set_band_idx(cur);
		cur++;

		org.x = dim.x; org.y = dim.y;
		par_n = i ? cur-3 : 0;
		bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, dim, par_n)
								  : new SUBBAND_TYPE(org, dim, par_n, msb_layers[cur]);
		assert(bands[cur]);
		bands[cur]->set_wp(this);  
		bands[cur]->set_band_idx(cur);
		cur++;

		dim.x <<= 1; dim.y <<= 1;
	}else
	{
		// now there are four small bands in HL subband at level 5 
		small_dim.x = dim.x/2;   small_dim.y = dim.y/2; 
		par_n = 100;
		for (m=0; m<2; m++)
		for (n=0; n<2; n++)
		{
			org.x = m*small_dim.x;  org.y = dim.y+n*small_dim.y;
			bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, small_dim, par_n)
								      : new SUBBAND_TYPE(org, small_dim, par_n, msb_layers[cur]);
			assert(bands[cur]);
			bands[cur]->set_wp(this);  
			bands[cur]->set_band_idx(cur);
			cur++;
		}
 
		// now there are four small bands in LH subband at level 5 
		for (m=0; m<2; m++)
		for (n=0; n<2; n++)
		{
			org.x = dim.x+m*small_dim.x;  org.y = n*small_dim.y;
			bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, small_dim, par_n)
								      : new SUBBAND_TYPE(org, small_dim, par_n, msb_layers[cur]);
			assert(bands[cur]);
			bands[cur]->set_wp(this);  bands[cur]->set_band_idx(cur);
			cur++;
		}

		// now there are four small bands in HH subband at level 5
		for (m=0; m<2; m++)
		for (n=0; n<2; n++)
		{
			org.x = dim.x+m*small_dim.x;  org.y = dim.y+n*small_dim.y;
			bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, small_dim, par_n)
								      : new SUBBAND_TYPE(org, small_dim, par_n, msb_layers[cur]);
			assert(bands[cur]);
			bands[cur]->set_wp(this);  bands[cur]->set_band_idx(cur);
			cur++;
		}

		dim.x <<= 1; dim.y <<= 1;
	}

#else
    org.x = 0; org.y = dim.y;
    par_n = i ? cur-3 : 0;
    bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, dim, par_n)
      : new SUBBAND_TYPE(org, dim, par_n, msb_layers[cur]);
    assert(bands[cur]);
    bands[cur]->set_wp(this);  bands[cur]->set_band_idx(cur);
    cur++;

    org.x = dim.x; org.y = 0;
    par_n = i ? cur-3 : 0;
    bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, dim, par_n)
      : new SUBBAND_TYPE(org, dim, par_n, msb_layers[cur]);
    assert(bands[cur]);
    bands[cur]->set_wp(this);  bands[cur]->set_band_idx(cur);
    cur++;

    org.x = dim.x; org.y = dim.y;
    par_n = i ? cur-3 : 0;
    bands[cur] = (!msb_layers)? new SUBBAND_TYPE(image_pyr, org, dim, par_n)
      : new SUBBAND_TYPE(org, dim, par_n, msb_layers[cur]);
    assert(bands[cur]);
    bands[cur]->set_wp(this);  bands[cur]->set_band_idx(cur);
    cur++;

    dim.x <<= 1; dim.y <<= 1;
#endif 

  }

  for( i = 1, max_msb = bands[0]->get_max_msb(  ); i < nband; i++ ) {
    if( max_msb < bands[i]->get_max_msb(  ) )
      max_msb = bands[i]->get_max_msb(  );
  }
}

//-----------------------------------------------------------------------------

//    reset_max_msb

//-----------------------------------------------------------------------------
void
SubSetLayer::reset_max_msb(  )
{
  int i;

  for( i = 1, max_msb = bands[0]->get_max_msb(  ); i < nband; i++ ) {
    if( max_msb < bands[i]->get_max_msb(  ) )
      max_msb = bands[i]->get_max_msb(  );

  }

}


/*-----------------------------------------------------------------------------

This function tries to reconstructs the image pyramid from the
received bit planes. The reconstructed coeffs take the input bitplane
values after bitplane shifting.

-----------------------------------------------------------------------------*/

void
SubSetLayer::SubSetLayer2Image_BW( Image_BW & image_pyr )
{
  int shift, lsb, sig, i, j, k;
  Image_Coord dim, org;
  SUB_COEFF_TYPE **base_coeff;

  for( k = 0; k < nband; k++ ) {
    if( bands[k]->get_transpose_flag(  ) )
      bands[k]->transpose(  );
    org = bands[k]->get_org(  );
    dim = bands[k]->get_dim(  );
    lsb = bands[k]->get_lsb(  );
    base_coeff = bands[k]->get_coeff(  );
    shift = ( lsb < 0 ) ? (-lsb) : 0;
    for( i = 0; i < dim.x; i++ ) {
      for( j = 0; j < dim.y; j++ ) {
        sig = ( base_coeff[i][j] & SIGN_BIT ) ? (-1) : 1;
        image_pyr( org.x + i, org.y + j ) = ( float )
          ( sig * ( ( base_coeff[i][j] & MAG_MASK ) >> shift ) );
      }
    }
  }
}

int
SubSetLayer::cal_band_level( int band_idx )
{
  if( band_idx )
    return pyr_levels - ( band_idx - 1 ) / 3;
  else
    return pyr_levels;
}
