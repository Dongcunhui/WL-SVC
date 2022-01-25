
#include "structN.h"
#include "ezbc_codec_3d.h"


EzbcCodec3d::EzbcCodec3d(  )
{
}

EzbcCodec3d::EzbcCodec3d( videoinfo & info )
{
  int i;

  dim.y = info.ywidth;
  dim.x = info.yheight;
  cdim.y = info.cwidth;
  cdim.x = info.cheight;
  GOPsz = info.GOPsz;           //byte_budget = info.GOPbytes;

  // by Yongjun Wu, save the resolution information
  // s_level=0: full resolution, s_level=1: half resolution, ...
  s_level = info.s_level;   
  if (info.org_yheight == 1080 && info.org_ywidth==1920)
	  frame_res = HD; 
  else if (info.org_yheight == 480 && info.org_ywidth == 832)
	  frame_res = SD; 
  else if (info.org_yheight == 288 && info.org_ywidth == 352)
	  frame_res = CIF;
  else if (info.org_yheight == 576 && info.org_ywidth == 704)
	  frame_res = SD2;
  else
	  frame_res = INVALID_RES; 


  ncomps = ( cdim.x > 0 ) ? 3 : 1;

  NEW_VECTOR( imgY, GOPsz, Image_BW, "imgY" );
  NEW_VECTOR( pyrY, GOPsz, SUBBAND_TREE_TYPE *, "pyrY" );
  for( i = GOPsz - 1; i >= 0; pyrY[i--] = NULL );

  if( cdim.x > 0 ) {
    NEW_VECTOR( imgU, GOPsz, Image_BW, "imgU" );
    NEW_VECTOR( imgV, GOPsz, Image_BW, "imgV" );
    NEW_VECTOR( pyrU, GOPsz, SUBBAND_TREE_TYPE *, "pyrU" );
    NEW_VECTOR( pyrV, GOPsz, SUBBAND_TREE_TYPE *, "pyrV" );
    for( i = GOPsz - 1; i >= 0; i-- ) {
      pyrU[i] = pyrV[i] = NULL;
    }
  } else {
    imgU = imgV = NULL;
    pyrU = pyrV = NULL;
  }
}

EzbcCodec3d::~EzbcCodec3d( void )
{
  int i, j;

  if( imgY ) {
    DELETE_VECTOR( imgY );
  }

  if( imgU ) {
    DELETE_VECTOR( imgU );
  }

  if( imgV ) {
    DELETE_VECTOR( imgV );
  }

  if( pyrY ) {
    for( i = GOPsz - 1; i >= 0; i-- ) {
      if( pyrY[i] )
        delete pyrY[i];
    }
    DELETE_VECTOR( pyrY );
  }

  if( pyrU ) {
    for( i = GOPsz - 1; i >= 0; i-- ) {
      if( pyrU[i] ) {
        DELETE_OBJECT( pyrU[i] );
      }
    }
    DELETE_VECTOR( pyrU );
  }
  if( pyrV ) {
    for( i = GOPsz - 1; i >= 0; i-- ) {
      if( pyrV[i] ) {
        DELETE_OBJECT( pyrV[i] );
      }
    }
    DELETE_VECTOR( pyrV );
  }

  for( i = 0; i < GOPsz; i++ ) {
    for( j = 0; j < ncomps; j++ ) {
      DELETE_VECTOR( subband_ACcoder[j][i] );
    }
  }
  for( i = 0; i < ncomps; i++ ) {
    DELETE_VECTOR( subband_ACcoder[i] );
  }

}

void EzbcCodec3d::initialize(  )     //0616
{
  int i, j, m, n;
  SUBBAND_TREE_TYPE **GOPpyr[3] = { pyrY, pyrU, pyrV };
  int roll_level1, roll_level2;   // for frequency roll off 

  max_Slev = 0;

  for( i = 0; i < ncomps; i++ ) {
    num_Slev[i] = GOPpyr[i][0]->get_pyr_levels(  ) + 1;
    if( max_Slev < num_Slev[i] )
      max_Slev = num_Slev[i];

    nbands[i] = GOPpyr[i][0]->get_nband(  );
  }

  for( i = 0; i < ncomps; i++ ) {
    NEW_VECTOR( subband_ACcoder[i], GOPsz, int *, "subband_ACcoder[i]" );
  }
  for( i = 0; i < GOPsz; i++ ) {
    for( j = 0; j < ncomps; j++ ) {
      NEW_VECTOR( subband_ACcoder[j][i], nbands[j], int,
                  "subband_ACcoder[j][i]" );
    }
  }
  int YUVband[3] = { 0, 0, 0 };

  total_num_ACcoder = 0;

  assert(frame_res != INVALID_RES ); 
  if (frame_res == HD)  // for HD sequence there are two resolution scalability, added by Yuan Liu
  {
	  roll_level1 = 2-s_level; 
	  roll_level2 = 3-s_level; 
  }else if (frame_res == SD2)  // for SD sequence there are two resolution scalability
  {
	  roll_level1 = 2-s_level; 
	  roll_level2 = 3-s_level; 
  }else if (frame_res == SD)  // for SD sequence there are two resolution scalability
  {
	  roll_level1 = 2-s_level; 
	  roll_level2 = 2-s_level; 
  }else if (frame_res == CIF )
	  roll_level1 = roll_level2 = 2-s_level; // for CIF sequence there is only one resolution scalability
     

  for( i = 0; i < GOPsz; i++ ) { // Öð¸öÖ¡½øÐÐ
	  for( m = max_Slev; m >= 1; m-- ) {  // the highest resolution is m=1, the second resolution is m=2, ...
      for( j = 0; j < ncomps; j++ ) {

        if( num_Slev[j] >= m ) {
          if( num_Slev[j] == m ) {
            subband_ACcoder[j][i][YUVband[j]] = total_num_ACcoder;
            YUVband[j]++;
          } else {

#ifdef       FREQUENCY_ROLL_OFF  // by Yongjun Wu
			  int small_band; 

#ifdef       ROLL_STRUCTURE_ONE
			  if (m!=roll_level1  && m!=roll_level2)  
			  {
				  // the three subbands are coded by the same arithematic coder 
				  for(n=0; n<3; n++){	
					subband_ACcoder[j][i][YUVband[j]] = total_num_ACcoder;
					YUVband[j]++;
				  }
			  }else
			  {
				  // the three subbands are coded separately
				  for (small_band=0; small_band<3; small_band++)
				  {
					  subband_ACcoder[j][i][YUVband[j]] = total_num_ACcoder+small_band;
					  YUVband[j]++;
				  }
			  }
#endif 

#ifdef       ROLL_STRUCTURE_TWO
			  if (m!=roll_level1  && m!=roll_level2)
			  {
				  // the other levels have three subbands
				  for(n=0; n<3; n++){	   //spatial subbands at the same resolution level
					subband_ACcoder[j][i][YUVband[j]] = total_num_ACcoder;
					YUVband[j]++;
				  }
			  }else if (m==roll_level1) // the 12 small subbands are coded separately 
			  {
				  for (small_band=0; small_band<12; small_band++)
				  {
					  subband_ACcoder[j][i][YUVband[j]] = total_num_ACcoder+small_band;
					  YUVband[j]++;
				  }
			  } else if (m==roll_level2)  // the 3 small subbands are coded by different arithematic coders
			  {
				  for (small_band=0; small_band<3; small_band++)
				  {
					  subband_ACcoder[j][i][YUVband[j]] = total_num_ACcoder+small_band;
					  YUVband[j]++;
				  }
			  }

#endif 

#else
			  // the three subbands are coded together
			  for(n=0; n<3; n++){	   
				subband_ACcoder[j][i][YUVband[j]] = total_num_ACcoder;
				YUVband[j]++;
			  }
#endif
		  
		  }
        }
      }                         //j  

#ifdef 	FREQUENCY_ROLL_OFF
	  total_num_ACcoder++;

#ifdef ROLL_STRUCTURE_ONE
	   if (m==roll_level1  || m==roll_level2)
			total_num_ACcoder += 2;
#endif 

#ifdef ROLL_STRUCTURE_TWO
	   if (m==roll_level1)
			total_num_ACcoder += 11;
	   if ( m==roll_level2)
			total_num_ACcoder += 2;
#endif 


#else
	   total_num_ACcoder++; // increase identity of Arithematic Coder
#endif

    }
  }

}
