
#define EXTERN  extern
#include "general.h"
#include "structN.h"
#include "coderN.h"
#include "bme_tools.h"
#include "basic.h"
#include "analsyn.h"

// 5/3 filter normalization
const float alpha = (float)sqrt(23.0 / 32.0);
const float beta  = (float)sqrt( 3.0 /  2.0);
#define DRIFT_PENALTY_IBLOCK	2.5
#define PROPAGATION_PENALTY_IBLOCK	10.0


     // block neighbors: ( example for 4x4 block )
  	 //
	 // neighbor A ... H  for some block X
	 //                 cx-1          cx+blksize
 	 //                  |			  |
     //                  H   A A A A   B B B B   --- cy-1 row
     //                  G   X X X X   C 
     //                  G   X X X X   C
     //                  G   X X X X   C
     //                  G   X X X X   C
     //                  F   E E E E   D D D D   --- cy+blksize row
     //                  F
	 //                  F
	 //                  F
	 //


// SPATIAL_VERTICAL prediction for block X
char spatial_ver_pre(float *predict_blk, float *neighborA, float *neighborE, int xblk, int yblk, char blkthresh)
{
	int   i, x, y;
	char  A_sign, E_sign;

	if (xblk!=yblk)
	{
		assert(0);
		return 0;
	}

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder ( U V components)
		return 0;
	else
	{
		E_sign = A_sign = 1;  // assume neighbor A and E are available
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				E_sign = 0;
		}
	}


	// the processing is the same for 2x2, 4x4, and 8x8 blocks
	if (!A_sign && !E_sign)
		return 0;  // no vertical spatial prediction, since neighbor A and neighbor E are not available

	if (A_sign && !E_sign)
		// neighbor A available and neighbor E not
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = neighborA[x];

	if (!A_sign && E_sign)
		// neighbor A not available and neighbor E available
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = neighborE[x];

	if (A_sign && E_sign)
		// neighbor A available and neighbor E available
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)  // this is the best case: linear interpolation using neighbor A and neighbor E
			predict_blk[y*xblk+x] = (neighborA[x]*((yblk+1)-(y+1)) + neighborE[x]*(y+1))/(float)(yblk+1);

	return 1;  // finsih spatial prediction successfully
}


char spatial_hor_pre(float *predict_blk, float *neighborG, float *neighborC, int xblk, int yblk, char blkthresh)
{
	int  i, x, y;
	char G_sign, C_sign;


	if (xblk!=yblk)
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder (U V components)
		return 0;
	else
	{
		G_sign = C_sign = 1;  // assume neighbor G and C are available
		for (i=0; i<xblk; i++)
		{
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
		}
	}

	// the processing is the same for 1x1, 2x2, 4x4, and 8x8 blocks 
	if (!G_sign && !C_sign)
		return 0;  // no horizontal spatial prediction, since neighbor G and neighbor C are not available

	if (G_sign && !C_sign)
		// neighbor G available and neighbor C not
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = neighborG[y];

	if (!G_sign && C_sign)
		// neighbor G not available and neighbor C available
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = neighborC[y];

	if (G_sign && C_sign)
		// neighbor G available and neighbor C available
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)  // linear interpolation using neighbor G and neighbor C
			predict_blk[y*xblk+x] = (neighborG[x]*((xblk+1)-(x+1)) + neighborC[x]*(x+1))/(float)(xblk+1);

	return 1;  // finish spaitial horizontal prediction successfully 
}


char spatial_dc_pre(float *predict_blk, float *neighborA, float *neighborC, 
					float *neighborE,   float *neighborG, int xblk, int yblk, char blkthresh)
{
	int   i, x, y;
	float dc_value;
	char  A_sign, E_sign, G_sign, C_sign;

	if (xblk!=yblk)
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder  (U V components)
		return 0;
	else
	{
		E_sign = A_sign = 1;  // assume neighbor A and E are available
		G_sign = C_sign = 1;  // assume neighbor G and C are available
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor E is not available
				E_sign = 0;
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
		}
	}

	// the processing is the same for 1x1, 2x2, 4x4, 8x8 blocks
	// 0000
    if (!A_sign && !C_sign && !E_sign && !G_sign)
		return 0;  // no neighbors available

	// neighbor A available 1000: average neighbor A only
    if (A_sign && !C_sign && !E_sign && !G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		dc_value /= (float)xblk;
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor C available  0100: average neighbor C only
    if (!A_sign && C_sign && !E_sign && !G_sign)
	{
		dc_value = 0;
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		dc_value /= (float)yblk;
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor A and C available 1100: average neighbor A and C
    if (A_sign && C_sign && !E_sign && !G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		dc_value /= (float)(xblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor E available 0010: average neighbor E only
    if (!A_sign && !C_sign && E_sign && !G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		dc_value /= (float)xblk;
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor A and E available  1010: average neighbor A and E
    if (A_sign && !C_sign && E_sign && !G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		dc_value /= (float)(xblk+xblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor C and E available  0110: average neighbor C and E
    if (!A_sign && C_sign && E_sign && !G_sign)
	{
		dc_value = 0;
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		dc_value /= (float)(xblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor A, C and E available 1110: average neighbor A, C and E
    if (A_sign && C_sign && E_sign && !G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		dc_value /= (float)(2*xblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor G available 0001: average neighbor G only
    if (!A_sign && !C_sign && !E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}


	// neighbor A and G available 1001: average neighbor A and G
    if (A_sign && !C_sign && !E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(xblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor C and G available 0101: average neighbor C and G
    if (!A_sign && C_sign && !E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(yblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor A, C and G available 1101: average A, C and G
    if (A_sign && C_sign && !E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(xblk+yblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}


	// neighbor E and G available 0011: average neighbor E and G
    if (!A_sign && !C_sign && E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(xblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor A, E and G available  1011: average neighbor A, E and G
    if (A_sign && !C_sign && E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(xblk+xblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor C, E and G available  0111: average neighbor C, E and G
    if (!A_sign && C_sign && E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(xblk+yblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	// neighbor A, C, E and G available  1111: the best case --- average neighbor A, C, E and G
    if (A_sign && C_sign && E_sign && G_sign)
	{
		dc_value = 0;
		for (i=0; i<xblk; i++)
			dc_value += neighborA[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborC[i];
		for (i=0; i<xblk; i++)
			dc_value += neighborE[i];
		for (i=0; i<yblk; i++)
			dc_value += neighborG[i];
		dc_value /= (float)(xblk+xblk+yblk+yblk);
		for (x=0; x<xblk; x++)
	    for (y=0; y<yblk; y++)
			predict_blk[y*xblk+x] = dc_value;
	}

	return 1;
}


// SPATIAL_DOWN_LEFT prediction for block X
char spatial_downleft_pre(float *predict_blk, float *neighborA, float *neighborG, float *neighborB, 
						  float *neighborF,    float *neighborC, float *neighborE, int xblk, int yblk, char blkthresh)
{
	int i, x, y, diag_dist, down_dist, upper_dist;
	// the sign used to indicate whether the corresponding neighbor is available or not
	char A_sign, E_sign, C_sign, G_sign, B_sign1, F_sign1, B_sign2, F_sign2;

	// add some extension for the neighbors in this mode
	// that means when neighbor C not available, use neighbor B, 
	// when neighbor E not avialable use neighbor F  --- use as much information as possible

	// so we neighbor B and neighbor F are divided into two parts
	// first part is used to interpolate/predict the diagonal part
	// second part is the alternative for neighbor C and neighbor E

	if (xblk!=yblk)  // it's not square, not do spatial interpolation/prediction in this case
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder ( U V components )
		return 0;
	else
	{
		A_sign  = E_sign  = 1;  // assume neighbor A and E are available
		G_sign  = C_sign  = 1;  // assume neighbor G and C are available
		B_sign2 = F_sign2 = 1;  
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor E is not available
				E_sign = 0;
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
			if (neighborB[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor B is not available
				B_sign2 = 0;
			if (neighborF[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor F is not available
				F_sign2 = 0;
		}
	}
    B_sign1 = (neighborB[0]!=(float)HUGE_VAL)? 1 : 0;
	F_sign1 = (neighborF[0]!=(float)HUGE_VAL)? 1:  0;

	// spatial prediction in this mode is divided into 3 parts
	// left-upper part, diagonal,  right-down parts

	// diagonal interpolate/predict from neighbor B and neighbor F

	// neighbor B and neighbor F are not available
	if (!B_sign1 && !F_sign1)
		return 0;

	// neighbor B available,  neighbor F not available
	if (B_sign1 && !F_sign1)
	{
	   // do the prediction from upper-right to down-left
	   x = xblk-1; y = 0;  
	   do 
	   {
	 	  // predict according to the values in neighbor B
		  predict_blk[y*xblk+x] = neighborB[0];
          x--; y++;
	   }while(x>=0);
	}
	

	// neighbor B not available,  neighbor F available
	if (!B_sign1 && F_sign1)
	{
	   // do the prediction from upper-right to down-left
	   x = xblk-1; y = 0;  
	   do 
	   {
	 	  // predict according to the values in neighbor B
		  predict_blk[y*xblk+x] = neighborF[0];
          x--; y++;
	   }while(x>=0);
	}

	// neighbor B and neighbor F are available: --- the best case
	if (B_sign1 && F_sign1)
	{
		// do the prediction from upper-right to down-left
		x = xblk-1; y = 0;  
		upper_dist = 1;  down_dist = (xblk+1)-upper_dist;
		do 
		{
			// linear interpolation according to the values in neighbor B and neighbor F
			predict_blk[y*xblk+x] = (neighborB[0]*down_dist+neighborF[0]*upper_dist)/(float)(upper_dist+down_dist);
            x--; y++;
			upper_dist++;  down_dist--;	
		}while(x>=0);
	}


	if (xblk>1)  // when blksize is larger than 1, then the block has upper-left and down-right parts
	{
	// left-upper part is predicted from neighbor A and neighbor G
	// neighbor A and neighbor G are not available
	if (!A_sign && !G_sign)
		return 0;

	// neighbor A available,  neighbor G not available
	if (A_sign && !G_sign)
	{
		for (diag_dist=0; diag_dist<xblk-1; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			x = diag_dist; y = 0;  
			do 
			{
				// predict according to the values in neighbor A
				predict_blk[y*xblk+x] = neighborA[diag_dist+1];
                x--; y++;
			}while(x>=0);
		}
	}

	// neighbor A not available,  neighbor G available
	if (!A_sign && G_sign)
	{
		for (diag_dist=0; diag_dist<xblk-1; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			x = diag_dist; y = 0;  
			do 
			{
				// predict according to the values in neighbor G
				predict_blk[y*xblk+x] = neighborG[diag_dist+1];
                x--; y++;
			}while(x>=0);
		}
	}

	// neighbor A and neighbor G are available: --- the best case
	if (A_sign && G_sign)
	{
		for (diag_dist=0; diag_dist<xblk-1; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			x = diag_dist; y = 0;  
			upper_dist = 1;  down_dist = (diag_dist+2)-upper_dist;
			do 
			{
				// linear interpolation according to the values in neighbor A and neighbor G
				predict_blk[y*xblk+x] = (neighborA[diag_dist+1]*down_dist+neighborG[diag_dist+1]*upper_dist)/(float)(diag_dist+2);
                x--; y++;
				upper_dist++;  down_dist--;
			}while(x>=0);
		}
	}

	// right-down diagonal part is predicted from C , E  and the second parts of neighbor B and F
	
	// C, E, B2 and F2 are not available: 0000
	if (!C_sign && !E_sign && !B_sign2 && !F_sign2)
		return 0;

	// C available, E, F2 not available (do not care about B2): 1000 or 1010 (C E B2 F2)
	if (C_sign && !E_sign && !F_sign2)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			do 
			{
				// predict according to the values in neighbor C
				predict_blk[y*xblk+x] = neighborC[diag_dist-1];
                x--; y++;
			}while(y<yblk);
		}
	}

	// C not available, E available, B2 not available ( not care F): 0100 or 0101
	if (!C_sign && E_sign && !B_sign2)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			do 
			{
				// predict according to the values in neighbor E
				predict_blk[y*xblk+x] = neighborE[diag_dist-1];
                x--; y++;
			}while(y<yblk);
		}
	}

	// C, E are available: the best case (not care B2, F2) 1100, 1101, 1110, 1111
	if (C_sign && E_sign)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			upper_dist = 1;  down_dist = yblk-diag_dist;
			do 
			{
				// linear interpolation according to the values in neighbor C and E
				predict_blk[y*xblk+x] = (neighborC[diag_dist-1]*down_dist+neighborE[diag_dist-1]*upper_dist)/(float)(upper_dist+down_dist);
                x--; y++; upper_dist++; down_dist--;
			}while(y<yblk);
		}
	}

	// C, E, F2 not available, B2 available: 0010
	if (!C_sign && !E_sign && B_sign2 && !F_sign2)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			do 
			{
				// predict according to the values in the second part of neighbor B
				predict_blk[y*xblk+x] = neighborB[diag_dist];
                x--; y++;
			}while(y<yblk);
		}
	}
	

	// C not available, E, B2 available(not care about neighbor F2): 0110 or 0111
	if (!C_sign && E_sign && B_sign2)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			upper_dist = diag_dist+1;  down_dist = (xblk+1)-upper_dist;
			do 
			{
				// linear interpolation according to the values in neighbor C and E
				predict_blk[y*xblk+x] = (neighborB[diag_dist]*down_dist+neighborE[diag_dist-1]*upper_dist)/(float)(upper_dist+down_dist);
                x--; y++; upper_dist++; down_dist--;
			}while(y<yblk);
		}
	}

	// neighbor C, E, B2 not available, F2 available: 0001
	if (!C_sign && !E_sign && !B_sign2 && F_sign2)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			do 
			{
				// predict according to the values in the second part of neighbor F
				predict_blk[y*xblk+x] = neighborF[diag_dist];
                x--; y++;
			}while(y<yblk);
		}
	}

	// neighbor C available, E not available, F2 available (not care about neighbor B2): 1001 or 1011
	if (C_sign && !E_sign && F_sign2)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			upper_dist = 1;  down_dist = (xblk+1)-upper_dist;
			do 
			{
				// linear interpolation according to the values in neighbor C and the second part of neighbor F
				predict_blk[y*xblk+x] = (neighborF[diag_dist]*upper_dist+neighborC[diag_dist-1]*down_dist)/(float)(upper_dist+down_dist);
                x--; y++; upper_dist++; down_dist--;
			}while(y<yblk);
		}
	}

	// neighbor C not available, E not available, B2, F2 available: 0011
	if (!C_sign && !E_sign && B_sign2 && F_sign2)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-right to down-left
			y = diag_dist; x = xblk-1;  
			upper_dist = diag_dist+1;  down_dist = xblk;
			do 
			{
				// linear interpolation according to the values in the second parts of neighbor F and B
				predict_blk[y*xblk+x] = (neighborB[diag_dist]*down_dist+neighborF[diag_dist]*upper_dist)/(float)(upper_dist+down_dist);
                x--; y++; upper_dist++; down_dist--;
			}while(y<yblk);
		}
	}
  }

  return 1;  // successfully interpolate/predict block X
}


// NEW VERSION: SPATIAL_DOWNRIGHT prediction for block X
// it's very similar to SPATIAL_DOWNLEFT prediction for block X
char spatial_downright_pre(float *predict_blk, float *neighborA, float *neighborC, 
						   float *neighborG,  float *neighborE, float neighborH, float *neighborD,
						   float *neighborD1, int xblk, int yblk, char blkthresh)
{
	int  i, x, y, diag_dist, down_dist, upper_dist;
	char A_sign, E_sign, G_sign, C_sign, H_sign, D_sign, D_sign1, D_sign2;

	if (xblk!=yblk)  // it's not square, not do spatial prediction in this case
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder
		return 0;
	else
	{
		A_sign  = E_sign  = 1;  // assume neighbor A and E are available
		G_sign  = C_sign  = 1;  // assume neighbor G and C are available
		D_sign1 = D_sign2 = 1;  
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor E is not available
				E_sign = 0;
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
			if (neighborD[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor B is not available
				D_sign2 = 0;
			if (neighborD1[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor F is not available
				D_sign1 = 0;
		}
	}
    H_sign = (neighborH   !=(float)HUGE_VAL)? 1 : 0;
	D_sign = (neighborD[0]!=(float)HUGE_VAL)? 1:  0;

	// spatial prediction in this mode is divided into three parts
	// right-upper diagonal part, diagnal line  and  left-down diagonal part

	// diagonal is predict from neighbor H and neighbor D

	// neighbor H and neighbor D are not available
	if (!H_sign && !D_sign)
		return 0;

	// neighbor H available,  neighbor D not available
	if (H_sign && !D_sign)
	{
	   // do the prediction from upper-left to down-right
	   x = 0; y = 0;  
	   do 
	   {
	 	  // predict according to the values in neighbor H
		  predict_blk[y*xblk+x] = neighborH;
          x++; y++;
	   }while(x<xblk);
	}
	

	// neighbor H not available,  neighbor D available
	if (!H_sign && D_sign)
	{
	   // do the prediction from upper-left to down-right
	   x = 0; y = 0;  
	   do 
	   {
	 	  // predict according to the values in the first part of neighbor D (neighborD[0])
		  predict_blk[y*xblk+x] = neighborD[0];
          x++; y++;
	   }while(x<xblk);
	}

	// neighbor H and neighbor D are available
	if (H_sign && D_sign)
	{
		// do the prediction from upper-left to down-right
		x = 0; y = 0;  
		upper_dist = 1;  down_dist = xblk;
		do 
		{
			// linear interpolation according to the values in neighbor H and 
			// the first part of neighbor D (neighborD[0]
			predict_blk[y*xblk+x] = (neighborH*down_dist+neighborD[0]*upper_dist)/(float)(upper_dist+down_dist);
            x++; y++;
			upper_dist++;  down_dist--;	
		}while(x<xblk);
	}


	if (xblk>1)  // when blksize is larger than 1 then it will have the following right-upper and down-right parts
	{
	// right-upper diagonal part is predicted from neighbor A, C and the second part of neighbor D

	// neighbor A and neighbor C are not available: 000
	if (!A_sign && !C_sign && !D_sign2)
		return 0;

	// neighbor A available,  C, D2 not available: 100
	if (A_sign && !C_sign && !D_sign2)
	{
		for (diag_dist=1; diag_dist<xblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; y = 0;  
			do 
			{
				// predict according to the values in neighbor A
				predict_blk[y*xblk+x] = neighborA[diag_dist-1];
                x++; y++;
			}while(x<xblk);
		}
	}
	

	// neighbor A not available,  C available ( not care about D2 ): 010 or 011
	if (!A_sign && C_sign)
	{
		for (diag_dist=1; diag_dist<xblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; y = 0;  
			do 
			{
				// predict according to the values in neighbor C
				predict_blk[y*xblk+x] = neighborC[yblk-diag_dist];
                x++; y++;
			}while(x<xblk);
		}
	}

	// neighbor A and neighbor C are available (not care about D2) : 110 or 111
	if (A_sign && C_sign)   //---BUG HERE!!! CORRECTED!!!
	{
		for (diag_dist=1; diag_dist<xblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; y = 0;  
			upper_dist = 1; down_dist = yblk-diag_dist;
			do 
			{
				// interpolate according to the values in neighbor A and neighbor C
				predict_blk[y*xblk+x] = (neighborA[diag_dist-1]*down_dist+neighborC[yblk-diag_dist]*upper_dist)/(float)(upper_dist+down_dist);
                x++; y++; upper_dist++; down_dist--;
			}while(x<xblk);
		}
	}


	// neighbor A not available,  C not available, D2 available: 001
	if (!A_sign && !C_sign && D_sign2)
	{
		for (diag_dist=1; diag_dist<xblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; y = 0;  
			do 
			{
				// predict according to the values in the second part of neighbor D
				// from neighborD[1] to neighborD[7]
				predict_blk[y*xblk+x] = neighborD[diag_dist];
                x++; y++;
			}while(x<xblk);
		}
	}


	// neighbor A, D2 are available, C not available : 101
	if (A_sign && !C_sign && D_sign2)
	{
		for (diag_dist=1; diag_dist<xblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; y = 0;  
			upper_dist = 1; down_dist = yblk;
			do 
			{
				// interpolate according to the values in neighbor A and neighbor C
				predict_blk[y*xblk+x] = (neighborA[diag_dist-1]*down_dist+neighborD[diag_dist]*upper_dist)/(float)(upper_dist+down_dist);
                x++; y++; upper_dist++; down_dist--;
			}while(x<xblk);
		}
	}


	// right-down diagonal part is predicted from neighbor  G, neighbor E and neighborD1
	
	// neighbor G, E, D1 are not available: 000
	if (!G_sign && !E_sign && !D_sign1)
		return 0;

	// neighbor G available, neighbor E, D1 not available: 100
	if (G_sign && !E_sign  && !D_sign1)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; x = 0;  
			do 
			{
				// predict according to the values in neighbor G
				predict_blk[y*xblk+x] = neighborG[diag_dist-1];
                x++; y++;
			}while(y<yblk);
		}
	}

	// neighbor G not available, E available( not care about neighbor D1): 010 or 011
	if (!G_sign && E_sign)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; x = 0;  
			do 
			{
				// predict according to the values in neighbor E
				predict_blk[y*xblk+x] = neighborE[xblk-diag_dist];
                x++; y++;
			}while(y<yblk);
		}
	}

	// neighbor G and neighbor E are available ( not care about neighbor D1): 110 or 111
	if (G_sign && E_sign)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; x = 0;  
			upper_dist = 1;  down_dist = yblk-diag_dist;
			do 
			{
				// linear interpolation according to the values in neighbor C and E
				predict_blk[y*xblk+x] = (neighborG[diag_dist-1]*down_dist+neighborE[xblk-diag_dist]*upper_dist)/(float)(upper_dist+down_dist);
                x++; y++; upper_dist++; down_dist--;
			}while(y<yblk);
		}
	}

	// neighbor G, E not available, D1 available: 100
	if (!G_sign && !E_sign  && D_sign1)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; x = 0;  
			do 
			{
				// predict according to the values in neighbor D1
				predict_blk[y*xblk+x] = neighborD1[diag_dist];
                x++; y++;
			}while(y<yblk);
		}
	}


	// neighbor G available, E not available, D1 available: 101 
	if (G_sign && !E_sign && D_sign1)
	{
		for (diag_dist=1; diag_dist<yblk; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; x = 0;  
			upper_dist = 1;  down_dist = xblk;
			do 
			{
				// linear interpolation according to the values in neighbor G 
				// and neighborD1
				predict_blk[y*xblk+x] = (neighborG[diag_dist-1]*down_dist+neighborD1[diag_dist]*upper_dist)/(float)(upper_dist+down_dist);
                x++; y++; upper_dist++; down_dist--;
			}while(y<yblk);
		}
	}
	}  // if (xblk>1)

	return 1;
}


// NEW VERSION: SPATIAL_VERTICAL_RIGHT prediction for block X
char spatial_verright_pre(float *predict_blk, float *neighborA, float *neighborC, 
						  float *neighborG, float *neighborE, float neighborH, float *neighborD, 
						  int xblk, int yblk, char blkthresh)
{
	int  i, x, y, diag_dist, down_dist, upper_dist;
	char A_sign, E_sign, G_sign, C_sign, H_sign, D_sign, D_sign2, index;

	// it's not square, not do spatial prediction in this case --- simplify the spatial interpolation/prediction
	if (xblk!=yblk)  
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder
		return 0;
	else
	{
		A_sign  = E_sign  = 1;  // assume neighbor A and E are available
		G_sign  = C_sign  = 1;  // assume neighbor G and C are available
		D_sign2 = 1;  
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor E is not available
				E_sign = 0;
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
			if (neighborD[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor B is not available
				D_sign2 = 0;
		}
	}
    H_sign = (neighborH!=(float)HUGE_VAL)? 1 : 0;
	D_sign = (neighborD[0]!=(float)HUGE_VAL)? 1:0;

	// spatial prediction in this mode is divided into 5 parts for 4x4 and 8x8 blocks
    // however, spatial prediction in this mode only need neighbor H and E, neighbor A and D
	// for 2x2 block, so 2x2 block is a special case for this mode

	if (xblk>1)
	{
	// neighbor H and neighbor E are not available
	if (!H_sign && !E_sign)
		return 0;

	// neighbor H available,  neighbor E not available
	if (H_sign && !E_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; 
			y = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborH;
			y = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborH;
		}
	}
	

	// neighbor H not available,  neighbor E available
	if (!H_sign && E_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; 
			y = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborE[xblk/2];
			y = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborE[xblk/2];
		}
	}

	// neighbor H and neighbor E are available: the best case
	if (H_sign && E_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			x = diag_dist; 	y = diag_dist*2;  
			upper_dist = 2*diag_dist+1;  down_dist = (xblk+1)-upper_dist;
			predict_blk[y*xblk+x] = (neighborH*down_dist+neighborE[xblk/2]*upper_dist)/(float)(xblk+1);
			y = diag_dist*2+1;  upper_dist++; down_dist--;
			predict_blk[y*xblk+x] = (neighborH*down_dist+neighborE[xblk/2]*upper_dist)/(float)(xblk+1);
		}
	}


	// neighbor A and neighbor D are not available
	if (!A_sign && !D_sign)
		return 0;

	// neighbor A available,  neighbor D not available
	if (A_sign && !D_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-right to upper-left
			x = xblk-1-diag_dist; 
			y = yblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborA[xblk/2-1];
			y = yblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborA[xblk/2-1];
		}
	}
	

	// neighbor A not available,  neighbor D available
	if (!A_sign && D_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-right to upper-left
			x = xblk-1-diag_dist; 
			y = yblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborD[0];
			y = yblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborD[0];
		}
	}

	// neighbor A and neighbor D are available
	if (A_sign && D_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-right to upper-left
			x = xblk-1-diag_dist; 
			y = yblk-1-diag_dist*2;  
			down_dist = 2*diag_dist+1; upper_dist = (xblk+1)-down_dist;
			predict_blk[y*xblk+x] = (neighborA[xblk/2-1]*down_dist+neighborD[0]*upper_dist)/(float)(xblk+1);
			y = yblk-1-(diag_dist*2+1);
			down_dist++;  upper_dist--;
			predict_blk[y*xblk+x] = (neighborA[xblk/2-1]*down_dist+neighborD[0]*upper_dist)/(float)(xblk+1);
		}
	}

	
	if (xblk>2)  // only xblk>2, there will be the following cases: you can check the block diagram!
	{
		// neighbor A and E are not available
		if (!A_sign && !E_sign)
			return 0;

		// neighbor A available and neighbor E not available  
		if (A_sign && !E_sign)
		{
			for (index=0; index<(xblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborA[index];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborA[index];
				}
		}

		// neighbor A not available and  neighbor E available
		if (!A_sign && E_sign)
		{
			for (index=0; index<(xblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborE[xblk/2+1+index];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborE[xblk/2+1+index];
				}
		}

		// neighbor A and neighbor E are available: the best case
		if (A_sign && E_sign)
		{
			for (index=0; index<(xblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					y = diag_dist*2;  
					upper_dist = 2*diag_dist+1; down_dist = (xblk+1)-upper_dist; 
					predict_blk[y*xblk+x] = (neighborA[index]*down_dist+neighborE[xblk/2+1+index]*upper_dist)/(float)(xblk+1);
					y = diag_dist*2+1;
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborA[index]*down_dist+neighborE[xblk/2+1+index]*upper_dist)/(float)(xblk+1);
				}
		}
	

		// neighbor A, C and D2 not available: 000
		if (!A_sign && !C_sign && !D_sign2)
			return 0;

		// neighbor A available, C and D2 not available: 100
		if (A_sign && !C_sign  && !D_sign2)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					if (x>=xblk)  break;
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborA[index];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborA[index];
				}
		}

		// neighbor A not available, C available (not care about D2): 010 or 011
		if (!A_sign && C_sign)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					if (x>=xblk)  break;
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborC[yblk-2-2*(index-xblk/2)];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborC[yblk-2-2*(index-xblk/2)];
			}
		}


        // neighbor A and C available (not care about D2): 110 or 111
		if (A_sign && C_sign)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					if (x>=xblk)  break;
					upper_dist = 2*diag_dist+1;  down_dist = (yblk-1-2*(index-xblk/2))-upper_dist;
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = (neighborA[index]*down_dist+neighborC[yblk-2-2*(index-xblk/2)]*upper_dist)/(float)(upper_dist+down_dist);
					y = diag_dist*2+1;
					upper_dist++;     down_dist--;
					predict_blk[y*xblk+x] = (neighborA[index]*down_dist+neighborC[yblk-2-2*(index-xblk/2)]*upper_dist)/(float)(upper_dist+down_dist);
				}
		}

		// neighbor A, C not available, and D2 available: 001
		if (!A_sign && !C_sign  && D_sign2)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					if (x>=xblk)  break;
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborD[index-xblk/2+1];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborD[index-xblk/2+1];
				}
		}


        // neighbor A available, C not available, D2 available: 101
		if (A_sign && !C_sign && D_sign2)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					x = index+1+diag_dist; 
					if (x>=xblk)  break;
					upper_dist = 2*diag_dist+1;  down_dist = (yblk+1)-upper_dist;
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = (neighborA[index]*down_dist+neighborD[index-xblk/2+1]*upper_dist)/(float)(upper_dist+down_dist);
					y = diag_dist*2+1;
					upper_dist++;     down_dist--;
					predict_blk[y*xblk+x] = (neighborA[index]*down_dist+neighborD[index-xblk/2+1]*upper_dist)/(float)(upper_dist+down_dist);
				}
		}


		if (!E_sign && !G_sign)
			return 0;

		if (E_sign && !G_sign)
		{
			for (index=1; index<=xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from down_right to upper_left
					x = index-1-diag_dist; 
					if (x<0)  break;
					y = yblk-1-diag_dist*2;  
					predict_blk[y*xblk+x] = neighborE[index];
					y = yblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborE[index];
				}
		}

		if (!E_sign && G_sign)
		{
			for (index=1; index<=xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from down_right to upper_left
					x = index-1-diag_dist; 
					if (x<0)  break;
					y = yblk-1-diag_dist*2;  
					predict_blk[y*xblk+x] = neighborG[yblk-1-2*index];
					y = yblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborG[yblk-1-2*index];
				}
		}


		if (E_sign && G_sign)
		{
			for (index=1; index<=xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from down_right to upper_left
					x = index-1-diag_dist; 
					if (x<0)  break;
					y = yblk-1-diag_dist*2;
					down_dist = 2*diag_dist+1;  upper_dist = (2*index+1)-down_dist;
					predict_blk[y*xblk+x] = (neighborE[index]*upper_dist+neighborG[yblk-1-2*index]*down_dist)/(float)(upper_dist+down_dist);
					y = yblk-1-(diag_dist*2+1);
					down_dist++; upper_dist--;
					predict_blk[y*xblk+x] = (neighborE[index]*upper_dist+neighborG[yblk-1-2*index]*down_dist)/(float)(upper_dist+down_dist);
				}
		}

	}  // if (xblk>2)  only when xblk>2, it needs these parts
	} // if (xblk>1)
	else  // xblk==1
	{
		// neighbor H and neighbor E are not available
		if (!H_sign && !E_sign)
			return 0;

		// neighbor H available,  neighbor E not available
		if (H_sign && !E_sign)
			predict_blk[0] = neighborH;

		// neighbor H not available,  neighbor E available
		if (!H_sign && E_sign)
			predict_blk[0] = neighborE[0];

		// neighbor H and neighbor E are available: the best case
		if (H_sign && E_sign)
			predict_blk[0] = (neighborH+neighborE[0])/(float)2.0;
	}

	return 1;
}


// NEW VERSION: SPATIAL_HORIZONTAL_DOWN prediction for block X
// it's quite similar to SPATIAL_VERTICAL_RIGHT prediction for block X
char spatial_hordown_pre(float *predict_blk, float *neighborA, float *neighborC, 
						 float *neighborG, float *neighborE, float neighborH, float neighborD, 
						 float *neighborD1,
						 int xblk, int yblk, char blkthresh)
{
	int  i, x, y, diag_dist, down_dist, upper_dist;
	char A_sign, E_sign, G_sign, C_sign, H_sign, D_sign, D_sign1, index;

	// it's not square, not do spatial prediction in this case  --- simplify spatial interpolation/prediction
	if (xblk!=yblk)  
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder
		return 0;
	else
	{
		A_sign  = E_sign  = 1;  // assume neighbor A and E are available
		G_sign  = C_sign  = 1;  // assume neighbor G and C are available
		D_sign1 = 1;  
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor E is not available
				E_sign = 0;
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
			if (neighborD1[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor B is not available
				D_sign1 = 0;
		}
	}
    H_sign = (neighborH!=(float)HUGE_VAL)? 1 : 0;
	D_sign = (neighborD!=(float)HUGE_VAL)? 1:  0;


	// spatial prediction in this mode is divided into 5 parts for 4x4 and 8x8 blocks
    // however, spatial prediction in this mode only need neighbor H and E, neighbor A and D
	// for 2x2 block, so 2x2 block is a special case for this mode

	if (xblk>1)
	{
	// neighbor H and neighbor C are not available
	if (!H_sign && !C_sign)
		return 0;

	// neighbor H available,  neighbor C not available
	if (H_sign && !C_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; 
			x = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborH;
			x = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborH;
		}
	}

	// neighbor H not available,  neighbor C available
	if (!H_sign && C_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; 
			x = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborC[xblk/2];
			x = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborC[xblk/2];
		}
	}

	// neighbor H and neighbor C are available
	if (H_sign && C_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-left to down-right
			y = diag_dist; 	x = diag_dist*2;  
			upper_dist = 2*diag_dist+1;  down_dist = (xblk+1)-upper_dist;
			predict_blk[y*xblk+x] = (neighborH*down_dist+neighborC[xblk/2]*upper_dist)/(float)(xblk+1);
			x = diag_dist*2+1;  upper_dist++; down_dist--;
			predict_blk[y*xblk+x] = (neighborH*down_dist+neighborC[xblk/2]*upper_dist)/(float)(xblk+1);
		}
	}

    // the second part

	// neighbor G and neighbor D are not available
	if (!G_sign && !D_sign)
		return 0;

	// neighbor G available,  neighbor D not available
	if (G_sign && !D_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-right to upper-left
			y = yblk-1-diag_dist; 
			x = xblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborG[yblk/2-1];
			x = xblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborG[yblk/2-1];
		}
	}
	

	// neighbor G not available,  neighbor D available
	if (!G_sign && D_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-right to upper-left
			y = yblk-1-diag_dist; 
			x = xblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborD;
			x = yblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborD;
		}
	}

	// neighbor G and neighbor D are available
	if (G_sign && D_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-right to upper-left
			y = yblk-1-diag_dist; 
			x = xblk-1-diag_dist*2;  
			down_dist = 2*diag_dist+1; upper_dist = (xblk+1)-down_dist;
			predict_blk[y*xblk+x] = (neighborG[xblk/2-1]*down_dist+neighborD*upper_dist)/(float)(xblk+1);
			x = xblk-1-(diag_dist*2+1);
			down_dist++;  upper_dist--;
			predict_blk[y*xblk+x] = (neighborG[xblk/2-1]*down_dist+neighborD*upper_dist)/(float)(xblk+1);
		}
	}

	
	if (xblk>2)
	{
		// neighbor G and C are not available
		if (!G_sign && !C_sign)
			return 0;

		// neighbor G available and neighbor C not available  
		if (G_sign && !C_sign)
		{
			for (index=0; index<(yblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					y = index+1+diag_dist; 
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborG[index];
					x = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborG[index];
				}
		}

		// neighbor G not available and  neighbor C available
		if (!G_sign && C_sign)
		{
			for (index=0; index<(xblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					y = index+1+diag_dist; 
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborC[xblk/2+1+index];
					x = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborC[xblk/2+1+index];
				}
		}

		// neighbor G and neighbor C are available
		if (G_sign && C_sign)
		{
			for (index=0; index<(xblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					y = index+1+diag_dist; 
					x = diag_dist*2;  
					upper_dist = 2*diag_dist+1; down_dist = (xblk+1)-upper_dist; 
					predict_blk[y*xblk+x] = (neighborG[index]*down_dist+neighborC[xblk/2+1+index]*upper_dist)/(float)(xblk+1);
					x = diag_dist*2+1;
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborG[index]*down_dist+neighborC[xblk/2+1+index]*upper_dist)/(float)(xblk+1);
				}
		}
	

		if (!A_sign && !C_sign)
			return 0;

		// neighbor A  available and neighbor C not available
		if (A_sign && !C_sign)
		{
			for (index=xblk/2-1; index>0; index--)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from down-right to upper-left
					y = index-1-diag_dist;  
					if (y<0)  break;
					x = xblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborA[(xblk/2-index)*2-1];
					x = xblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborA[(xblk/2-index)*2-1];
			}
		}

		// neighbor A not available and neighbor C available
		if (!A_sign && C_sign)
		{
			for (index=xblk/2-1; index>0; index--)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from down-right to upper-left
					y = index-1-diag_dist;  
					if (y<0)  break;
					x = xblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborC[index];
					x = xblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborC[index];
			}
		}

        // both neighbor A and neighbor C are available
		if (A_sign && C_sign)
		{
			for (index=xblk/2-1; index>0; index--)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from down-right to upper-left
					y = index-1-diag_dist;  
					if (y<0)  break;
					x = xblk-1-2*diag_dist;  
					// 2*index+1 is the total distance, down_dist is the distance to the down-right neighbor 
					down_dist = 1+2*diag_dist; upper_dist = (2*index+1)-down_dist;
					predict_blk[y*xblk+x] = (neighborC[index]*upper_dist+neighborA[(xblk/2-index)*2-1]*down_dist)/(float)(upper_dist+down_dist);
					x = xblk-1-(diag_dist*2+1);
					down_dist++;  upper_dist--;  
					// BUG here! Corrected!
					predict_blk[y*xblk+x] = (neighborC[index]*upper_dist+neighborA[(xblk/2-index)*2-1]*down_dist)/(float)(upper_dist+down_dist);
					//neighborA[(xblk/2-index)*2-1];
			}
		}

		// neighbor E, G and D1 not available: 000
		if (!E_sign && !G_sign && !D_sign1)
			return 0;

		// neighbor E available, G not available (not care about D1):100 or 101
		if (E_sign && !G_sign)
		{
			for (index=yblk/2; index<yblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					// y position 
					y = index+1+diag_dist; 
					if (y>=yblk)  break;
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborE[(xblk-1-index)*2];
					x = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborE[(xblk-1-index)*2];
				}
		}

		// neighbor E not available, G available, D1 not available: 010
		if (!E_sign && G_sign && !D_sign1)
		{
			for (index=yblk/2; index<yblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					// y position 
					y = index+1+diag_dist; 
					if (y>=yblk)  break;
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborG[index];
					x = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborG[index];
				}
		}

        // neighbor E and neighbor G are available (not care about D1): 110 or 111
		if (E_sign && G_sign)
		{
			for (index=yblk/2; index<yblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					// y position 
					y = index+1+diag_dist; 
					// 2*(yblk-1-index)+1 is the total distance 
					upper_dist = 1+2*diag_dist;  down_dist = (2*(yblk-1-index)+1)-upper_dist;
					if (y>=yblk)  break;
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = (neighborG[index]*down_dist+neighborE[(xblk-1-index)*2]*upper_dist)/(float)(upper_dist+down_dist);
					x = diag_dist*2+1;
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborG[index]*down_dist+neighborE[(xblk-1-index)*2]*upper_dist)/(float)(upper_dist+down_dist);
				}
		}

		// neighbor E, G not available, D1 available: 001
		if (!E_sign && !G_sign && D_sign1)
		{
			for (index=yblk/2; index<yblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					// y position 
					y = index+1+diag_dist; 
					if (y>=yblk)  break;
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborD1[index-yblk/2+1];
					x = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborD1[index-yblk/2+1];
				}
		}

        // neighbor E not available, G and D1 are available: 011
		if (!E_sign && G_sign && D_sign1)
		{
			for (index=yblk/2; index<yblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do the prediction from upper-left to down-right
					// y position 
					y = index+1+diag_dist; 
					// xblk+1 is the total distance 
					upper_dist = 1+2*diag_dist;  down_dist = (xblk+1)-upper_dist;
					if (y>=yblk)  break;
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = (neighborG[index]*down_dist+neighborD1[index-yblk/2+1]*upper_dist)/(float)(upper_dist+down_dist);
					x = diag_dist*2+1;
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborG[index]*down_dist+neighborD1[index-yblk/2+1]*upper_dist)/(float)(upper_dist+down_dist);
				}
		}
	}  // if (xblk>2)  only when xblk>2, it needs these parts
	} // if (xblk>1)
	else  // xblk==1
	{
		// neighbor H and neighbor C are not available
		if (!H_sign && !C_sign)
			return 0;

		// neighbor H available,  neighbor C not available
		if (H_sign && !C_sign)
			predict_blk[0] = neighborH;

		// neighbor H not available,  neighbor C available
		if (!H_sign && C_sign)
			predict_blk[0] = neighborC[0];

		// neighbor H and neighbor C are available
		if (H_sign && C_sign)
			predict_blk[0] = (neighborH+neighborC[0])/(float)2.0;
	}

	return 1;
}


// NEW VERSION: SPATIAL_VERTICAL_LEFT prediction for block X
// it's quite similar to SPATIAL_VERTICAL_RIGHT prediction for block X
char spatial_verleft_pre(float *predict_blk, float *neighborA, float *neighborC, 
						 float *neighborG, float *neighborE, float *neighborB,  float neighborF, 
						 int xblk, int yblk, char blkthresh)
{
	int  i, x, y, diag_dist, down_dist, upper_dist;
	char A_sign, E_sign, G_sign, C_sign, B_sign, B_sign1, F_sign, index;

	// it's not square, not do spatial prediction in this case  --- simplify spatial interpolation/prediction
	if (xblk!=yblk)  
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder
		return 0;
	else
	{
		A_sign  = E_sign  = 1;  // assume neighbor A and E are available
		G_sign  = C_sign  = 1;  // assume neighbor G and C are available
		B_sign1 = 1;  
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor E is not available
				E_sign = 0;
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
			if (neighborB[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor B is not available
				B_sign1 = 0;
		}
	}
    B_sign = (neighborB[0]!=(float)HUGE_VAL)? 1 : 0;
	F_sign = (neighborF!=(float)HUGE_VAL)?    1:  0;

	// spatial prediction in this mode is divided into 5 parts for 4x4 and 8x8 blocks
    // however, spatial prediction in this mode only need neighbor H and E, neighbor A and D
	// for 2x2 block, so 2x2 block is a special case for this mode

	if (xblk>1)
	{
	// neighbor A and neighbor F are not available
	if (!A_sign && !F_sign)
		return 0;

	// neighbor A available,  neighbor F not available
	if (A_sign && !F_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-left to upper-right
			x = diag_dist; 
			y = yblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborA[xblk/2];
			y = yblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborA[xblk/2];
		}
	}
	

	// neighbor A not available,  neighbor F available
	if (!A_sign && F_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from down-left to upper-right
			x = diag_dist; 
			y = yblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborF;
			y = yblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborF;
		}
	}

	// neighbor A and neighbor F are available
	if (A_sign && F_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do linear spatial interpolation from down-left to upper-right
			x = diag_dist; 
			y = yblk-1-diag_dist*2;  
			down_dist = 2*diag_dist+1; upper_dist = (yblk+1)-down_dist;
			predict_blk[y*xblk+x] = (neighborA[xblk/2]*down_dist+neighborF*upper_dist)/(float)(yblk+1);
			y = yblk-1-(diag_dist*2+1);
			down_dist++;  upper_dist--;
			predict_blk[y*xblk+x] = (neighborA[xblk/2]*down_dist+neighborF*upper_dist)/(float)(yblk+1);
		}
	}


	// neighbor B and neighbor E are not available
	if (!B_sign && !E_sign)
		return 0;

	// neighbor B available,  neighbor E not available
	if (B_sign && !E_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-right to down-left (B->E)
			x = xblk-1-diag_dist; 
			y = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborB[0];
			y = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborB[0];
		}
	}
	

	// neighbor B not available,  neighbor E available
	if (!B_sign && E_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do spatial prediction from upper-right to down-left (B->E)
			x = xblk-1-diag_dist; 
			y = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborE[xblk/2-1];
			y = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborE[xblk/2-1];
		}
	}

	// neighbor B and neighbor E are available
	if (B_sign && E_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do spatial interpolation from upper-right to down-left (B->E)
			x = xblk-1-diag_dist; 
			y = diag_dist*2;  
			upper_dist = diag_dist*2+1;  down_dist = (yblk+1)-upper_dist;
			predict_blk[y*xblk+x] = (neighborE[xblk/2-1]*upper_dist+neighborB[0]*down_dist)/(float)(yblk+1);
			y = diag_dist*2+1;
			upper_dist++; down_dist--;
			predict_blk[y*xblk+x] = (neighborE[xblk/2-1]*upper_dist+neighborB[0]*down_dist)/(float)(yblk+1);
		}
	}

	
	if (xblk>2)
	{
		// neighbor A and E are not available
		if (!A_sign && !E_sign)
			return 0;

		// neighbor A available and neighbor E not available  
		if (A_sign && !E_sign)
		{
			for (index=0; index<(yblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-left to down-right (A->E)
					x = xblk/2+index-diag_dist; 
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborA[xblk/2+1+index];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborA[xblk/2+1+index];;
				}
		}

		// neighbor A not available and  neighbor E available
		if (!A_sign && E_sign)
		{
			for (index=0; index<(yblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-left to down-right (A->E)
					x = xblk/2+index-diag_dist; 
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborE[index];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborE[index];
				}
		}

		// neighbor A available and  neighbor E available
		if (A_sign && E_sign)
		{
			for (index=0; index<(yblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial interpolation from upper-left to down-right (A->E)
					x = xblk/2+index-diag_dist; 
					y = diag_dist*2;
					upper_dist = 2*diag_dist+1;  down_dist = (yblk+1)-upper_dist;
					predict_blk[y*xblk+x] = (neighborA[xblk/2+1+index]*down_dist+neighborE[index]*upper_dist)/(float)(yblk+1);
					y = diag_dist*2+1;
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborA[xblk/2+1+index]*down_dist+neighborE[index]*upper_dist)/(float)(yblk+1);
				}
		}
	

		// C, E and B1 not available: 000
		if (!C_sign && !E_sign && !B_sign1)
			return 0;

		// neighbor C available, E not available (not care about B1): 100 or 101
		if (C_sign && !E_sign)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from down-right to upper-left (E->C)
					x = index+1+diag_dist;  					
					if (x>=xblk)  break;
					y = yblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborC[(index-xblk/2)*2+1];
					y = yblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborC[(index-xblk/2)*2+1];
			}
		}

		// neighbor C not available,  E available, B1 not available: 010
		if (!C_sign && E_sign && !B_sign1)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from down-right to upper-left (E->C)
					x = index+1+diag_dist;  					
					if (x>=xblk)  break;
					y = yblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborE[index];
					y = yblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborE[index];
			}
		}

        // both neighbor C, E are available (not care about B1): 110 or 111
		if (C_sign && E_sign)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial interpolation from down-right to upper-left (E->C)
					x = index+1+diag_dist;  					
					if (x>=xblk)  break;
					y = yblk-1-2*diag_dist; 
					// 2*(xblk-index)-1 is the total distance
					down_dist = 2*diag_dist+1;  upper_dist = (2*(xblk-index)-1)-down_dist;
					predict_blk[y*xblk+x] = (neighborC[(index-xblk/2)*2+1]*down_dist+neighborE[index]*upper_dist)/(float)(upper_dist+down_dist);
					y = yblk-1-(diag_dist*2+1);
					down_dist++;   upper_dist--;
					predict_blk[y*xblk+x] = (neighborC[(index-xblk/2)*2+1]*down_dist+neighborE[index]*upper_dist)/(float)(upper_dist+down_dist);
			}
		}

		// neighbor C not available,  E not available, B1 available: 001
		if (!C_sign && !E_sign && B_sign1)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from down-right to upper-left (E->B1)
					x = index+1+diag_dist;  					
					if (x>=xblk)  break;
					y = yblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborB[index-xblk/2+1];
					y = yblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborB[index-xblk/2+1];
			}
		}


        // C not available, E and B1 are available: 011
		if (!C_sign && E_sign && B_sign1)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial interpolation from down-right to upper-left (E->C)
					x = index+1+diag_dist;  					
					if (x>=xblk)  break;
					y = yblk-1-2*diag_dist; 
					// yblk+1 is the total distance
					down_dist = 2*diag_dist+1;  upper_dist = (yblk+1)-down_dist;
					predict_blk[y*xblk+x] = (neighborB[index-xblk/2+1]*down_dist+neighborE[index]*upper_dist)/(float)(upper_dist+down_dist);
					y = yblk-1-(diag_dist*2+1);
					down_dist++;   upper_dist--;
					predict_blk[y*xblk+x] = (neighborB[index-xblk/2+1]*down_dist+neighborE[index]*upper_dist)/(float)(upper_dist+down_dist);
			}
		}


		if (!A_sign && !G_sign)
			return 0;

		// neighbor A available and neighbor G not available
		if (A_sign && !G_sign)
		{
			for (index=0; index<xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-left to down-right (A->G)
					// y position 
					x = index-diag_dist; 
					if (x<0)  break;
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborA[index+1];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborA[index+1];
				}
		}

		// neighbor A not available and neighbor G available
		if (!A_sign && G_sign)
		{
			for (index=0; index<xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-left to down-right (A->G)
					// y position 
					x = index-diag_dist; 
					if (x<0)  break;
					y = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborG[(index+1)*2];
					y = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborG[(index+1)*2];
				}
		}

        // neighbor A and neighbor G are available -- spatial interpolation
		if (A_sign && G_sign)
		{
			for (index=0; index<xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial interpolation from upper-left to down-right
					// y position 
					x = index-diag_dist; 
					if (x<0)  break;
					y = diag_dist*2; 
					// 2*index+3 is the total distance between neighbor A and neighbor G
					upper_dist = diag_dist*2+1;  down_dist = (2*index+3)-upper_dist;
					predict_blk[y*xblk+x] = (neighborA[index+1]*down_dist+neighborG[(index+1)*2]*upper_dist)/(float)(upper_dist+down_dist);
					y = diag_dist*2+1;
					upper_dist++;   down_dist--;
					predict_blk[y*xblk+x] = (neighborA[index+1]*down_dist+neighborG[(index+1)*2]*upper_dist)/(float)(upper_dist+down_dist);
				}
		}
	}  // if (xblk>2)  only when xblk>2, it needs these parts
	}  //if (xblk>1)
	else  // xblk==1
	{
		// neighbor B and neighbor E are not available
		if (!B_sign && !E_sign)
			return 0;

		// neighbor B available,  neighbor E not available
		if (B_sign && !E_sign)
			predict_blk[0] = neighborB[0];

		// neighbor B not available,  neighbor E available
		if (!B_sign && E_sign)
			predict_blk[0] = neighborE[0];

		// neighbor B and neighbor E are available
		if (B_sign && E_sign)
			predict_blk[0] = (neighborE[0]+neighborB[0])/(float)2.0;

	}

	return 1;
}


// NEW VERSION
char spatial_horup_pre(float *predict_blk, float *neighborA, float *neighborC, float *neighborG, 
					   float *neighborE, float neighborB, float *neighborF, int xblk, int yblk, char blkthresh)
{
	int  i, x, y, diag_dist, down_dist, upper_dist;
	char A_sign, E_sign, G_sign, C_sign, B_sign, F_sign, F_sign1, index;

	// it's not square, not do spatial prediction in this case  --- simplify spatial interpolation/prediction
	if (xblk!=yblk)  
		return 0;

	// check whether the block size is in the allowed value set 
	// xblk==1 is in the case of resolution scalability
	if (!(xblk==1 || xblk==2 || xblk==4 || xblk==8 || xblk==16))
		return 0; 

	if (xblk==0)  // this is for resolution scalability in decoder
		return 0;
	else
	{
		A_sign  = E_sign  = 1;  // assume neighbor A and E are available
		G_sign  = C_sign  = 1;  // assume neighbor G and C are available
		F_sign1 = 1;  
		for (i=0; i<xblk; i++)
		{
			if (neighborA[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor A is not available
				A_sign = 0;
			if (neighborE[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor E is not available
				E_sign = 0;
			if (neighborG[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor G is not available
				G_sign = 0;
			if (neighborC[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor C is not available
				C_sign = 0;
			if (neighborF[i]==(float)HUGE_VAL)  
				// if there is one element not available, then neighbor B is not available
				F_sign1 = 0;
		}
	}
    B_sign = (neighborB!=(float)HUGE_VAL)? 1 : 0;
	F_sign = (neighborF[0]!=(float)HUGE_VAL)? 1:  0;

	// spatial interpolation/prediction in this mode is divided into 5 parts for 4x4 and 8x8 blocks
    // however, spatial prediction in this mode only need neighbor H and E, neighbor A and D
	// for 2x2 block, so 2x2 block is a special case for this mode

	if (xblk>1)
	{
	// neighbor B and neighbor G are not available
	if (!B_sign && !G_sign)
		return 0;

	// neighbor B available, neighbor G not available
	if (B_sign && !G_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-right to down-left (B->G)
			y = diag_dist; 
			x = xblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborB;
			x = xblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborB;
		}
	}
	

	// neighbor B not available, neighbor G available
	if (!B_sign && G_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do the prediction from upper-right to down-left (B->G)
			y = diag_dist; 
			x = xblk-1-diag_dist*2;  
			predict_blk[y*xblk+x] = neighborG[yblk/2];
			x = xblk-1-(diag_dist*2+1);
			predict_blk[y*xblk+x] = neighborG[yblk/2];
		}
	}

	// neighbor B available, neighbor G available: the best case
	if (B_sign && G_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do spatial interpolation from upper-right to down-left (B->G)
			y = diag_dist; 
			x = yblk-1-diag_dist*2;  
			upper_dist = 2*diag_dist+1;  down_dist = (xblk+1)-upper_dist;
			predict_blk[y*xblk+x] = (neighborB*down_dist+neighborG[yblk/2]*upper_dist)/(float)(xblk+1);
			x = yblk-1-(diag_dist*2+1);
			upper_dist++; down_dist--;
			predict_blk[y*xblk+x] = (neighborB*down_dist+neighborG[yblk/2]*upper_dist)/(float)(xblk+1);
		}
	}


	// neighbor C and neighbor F are not available
	if (!C_sign && !F_sign)
		return 0;

	// neighbor C available,  neighbor F not available
	if (C_sign && !F_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do spatial prediction from down-left to upper-right (F->C)
			y = xblk-1-diag_dist; 
			x = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborC[yblk/2-1];
			x = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborC[yblk/2-1];
		}
	}
	

	// neighbor C not available,  neighbor F available
	if (!C_sign && F_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do spatial prediction from down-left to upper-right (F->C)
			y = xblk-1-diag_dist; 
			x = diag_dist*2;  
			predict_blk[y*xblk+x] = neighborF[0];
			x = diag_dist*2+1;
			predict_blk[y*xblk+x] = neighborF[0];
		}
	}

	// neighbor C available,  neighbor F available
	if (C_sign && F_sign)
	{
		for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
		{
			// do spatial prediction from down-left to upper-right (F->C)
			y = xblk-1-diag_dist; 
			x = diag_dist*2;  
			down_dist = 2*diag_dist+1;  upper_dist = (xblk+1)-down_dist;
			predict_blk[y*xblk+x] = (neighborC[yblk/2-1]*down_dist+neighborF[0]*upper_dist)/(float)(xblk+1);
			x = diag_dist*2+1;
			down_dist++; upper_dist--;
			predict_blk[y*xblk+x] = (neighborC[yblk/2-1]*down_dist+neighborF[0]*upper_dist)/(float)(xblk+1);
		}
	}

	
	if (xblk>2)
	{
		// neighbor G and C are not available
		if (!G_sign && !C_sign)
			return 0;

		// neighbor G available and neighbor C not available  
		if (G_sign && !C_sign)
		{
			for (index=0; index<(yblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from down-left to upper-right (G->C)
					y = yblk/2+index-diag_dist; 
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborG[yblk/2+1+index];
					x = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborG[yblk/2+1+index];;
				}
		}

		// neighbor G not available and neighbor C available  
		if (!G_sign && C_sign)
		{
			for (index=0; index<(yblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from down-left to upper-right (G->C)
					y = yblk/2+index-diag_dist; 
					x = diag_dist*2;  
					predict_blk[y*xblk+x] = neighborC[index];
					x = diag_dist*2+1;
					predict_blk[y*xblk+x] = neighborC[index];
				}
		}

		// neighbor G available and neighbor C not available  
		if (G_sign && C_sign)
		{
			for (index=0; index<(yblk/2-1); index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial interpolation from down-left to upper-right (G->C)
					y = yblk/2+index-diag_dist; 
					x = diag_dist*2;  
					down_dist = diag_dist*2+1;  upper_dist = (xblk+1)-down_dist;
					predict_blk[y*xblk+x] = (neighborC[index]*down_dist+neighborG[yblk/2+1+index]*upper_dist)/(float)(xblk+1);
					x = diag_dist*2+1;
					down_dist++; upper_dist--;
					predict_blk[y*xblk+x] = (neighborC[index]*down_dist+neighborG[yblk/2+1+index]*upper_dist)/(float)(xblk+1);
				}
		}
	

		if (!C_sign && !E_sign && !F_sign1)
			return 0;

		// neighbor C available and E and F1 not available: 100
		if (C_sign && !E_sign && !F_sign1)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-right to down-left (C->E)
					y = index+1+diag_dist;  					
					if (y>=yblk)  break;
					x = xblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborC[index];
					x = xblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborC[index];
			}
		}

		// neighbor C not available, E available( not care about F1): 010 or 011
		if (!C_sign && E_sign)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-right to down-left (C->E)
					y = index+1+diag_dist;  					
					if (y>=yblk)  break;
					x = xblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborE[(index-xblk/2)*2+1];
					x = xblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborE[(index-xblk/2)*2+1];
			}
		}

		// neighbor C available and E available( not care about F1): 110 or 111
		if (C_sign && E_sign)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-right to down-left (C->E)
					y = index+1+diag_dist;  					
					if (y>=yblk)  break;
					x = xblk-1-2*diag_dist; 
					// (xblk-1-index)*2+1 is the total distance
					upper_dist = diag_dist*2+1;  down_dist = (xblk-1-index)*2+1-upper_dist;
					predict_blk[y*xblk+x] = (neighborC[index]*down_dist+neighborE[(index-xblk/2)*2+1]*upper_dist)/(float)(upper_dist+down_dist);
					x = xblk-1-(diag_dist*2+1);
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborC[index]*down_dist+neighborE[(index-xblk/2)*2+1]*upper_dist)/(float)(upper_dist+down_dist);
			}
		}


		// neighbor C and E not available and F1 available: 001
		if (!C_sign && !E_sign && F_sign1)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-right to down-left (C->E)
					y = index+1+diag_dist;  					
					if (y>=yblk)  break;
					x = xblk-1-2*diag_dist;  
					predict_blk[y*xblk+x] = neighborF[index-xblk/2+1];
					x = xblk-1-(diag_dist*2+1);
					predict_blk[y*xblk+x] = neighborF[index-xblk/2+1];
			}
		}

		// neighbor C available and E not available, F1 available: 101
		if (C_sign && !E_sign && F_sign1)
		{
			for (index=xblk/2; index<xblk-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-right to down-left (C->F1)
					y = index+1+diag_dist;  					
					if (y>=yblk)  break;
					x = xblk-1-2*diag_dist; 
					// (xblk+1) is the total distance
					upper_dist = diag_dist*2+1;  down_dist = (xblk+1)-upper_dist;
					predict_blk[y*xblk+x] = (neighborC[index]*down_dist+neighborF[index-xblk/2+1]*upper_dist)/(float)(upper_dist+down_dist);
					x = xblk-1-(diag_dist*2+1);
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborC[index]*down_dist+neighborF[index-xblk/2+1]*upper_dist)/(float)(upper_dist+down_dist);
			}
		}


		if (!A_sign && !G_sign)
			return 0;

		// neighbor A available and neighbor G not available
		if (A_sign && !G_sign)
		{
			for (index=0; index<xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-right to down-left (A->G)
					// y position 
					y = diag_dist; 
					x = (index+1)*2-1-2*diag_dist;
					if (x<0)  break;
					predict_blk[y*xblk+x] = neighborA[(index+1)*2];
					x = (index+1)*2-1-(2*diag_dist+1);
					predict_blk[y*xblk+x] = neighborA[(index+1)*2];
				}
		}

		// neighbor A not available and neighbor G available
		if (!A_sign && G_sign)
		{
			for (index=0; index<xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial prediction from upper-right to down-left (A->G)
					// y position 
					y = diag_dist; 
					x = (index+1)*2-1-2*diag_dist;
					if (x<0)  break;
					predict_blk[y*xblk+x] = neighborG[index+1];
					x = (index+1)*2-1-(2*diag_dist+1);
					predict_blk[y*xblk+x] = neighborG[index+1];
				}
		}

		// neighbor A available and neighbor G available
		if (A_sign && G_sign)
		{
			for (index=0; index<xblk/2-1; index++)
				for (diag_dist=0; diag_dist<yblk/2; diag_dist++)
				{
					// do spatial interpolation from upper-right to down-left (A->G)
					// y position 
					y = diag_dist; 
					x = (index+1)*2-1-2*diag_dist;
					if (x<0)  break;
					// (index+1)*2+1 is the total distance between neighbor A and neighbor G
					upper_dist = 2*diag_dist+1;  down_dist = (index+1)*2+1-upper_dist;
					predict_blk[y*xblk+x] = (neighborA[(index+1)*2]*down_dist+neighborG[index+1]*upper_dist)/(float)(upper_dist+down_dist);
					x = (index+1)*2-1-(2*diag_dist+1);
					upper_dist++;  down_dist--;
					predict_blk[y*xblk+x] = (neighborA[(index+1)*2]*down_dist+neighborG[index+1]*upper_dist)/(float)(upper_dist+down_dist);
				}
		}
	}  // if (xblk>2)  only when xblk>2, it needs these parts
	} //  if (xblk>1)
	else  // xblk==1
	{
		// neighbor B and neighbor G are not available
		if (!B_sign && !G_sign)
			return 0;

		// neighbor B available, neighbor G not available
		if (B_sign && !G_sign)
			predict_blk[0] = neighborB;

		// neighbor B not available, neighbor G available
		if (!B_sign && G_sign)
			predict_blk[0] = neighborG[0];

		// neighbor B available, neighbor G available: the best case
		if (B_sign && G_sign)
			predict_blk[0] = (neighborB+neighborG[0])/(float)2.0;
	}


	return 1;
}

// cx, cy, xblk, yblk, yhor, yver, chor and cver are in full resolution 
void GetSpaitialNeighborPixels8x8(videoinfo info, float *fadditional_penalty, vector_ptr fmv, 
								  float *fr_cur, int cx, int cy, int xblk, int yblk, int yhor, int yver, 
								  int t_level, float *fr_curU, float *fr_curV)
{
	int   i, chor, cver, uvx, uvy, uvblk, x_dest, y_dest; 
	float  mvx_d, mvy_d;
	enum  BiMode block_mode1, block_mode2;
	enum  FLAG is_pred;
	int   s_level, blksize, uvblksize, start_x, start_y, s_yhor, s_yver; 
	int   s_chor, s_cver, start_uvx, start_uvy; 
	int   propagate_iblk1=0, propagate_iblk2=0; 
	
	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
    if( s_level > 0 ) {
		info.ywidth  <<= s_level;
		info.yheight <<= s_level;
		info.cwidth  <<= s_level;	
		info.cheight <<= s_level;
	}
	chor  = info.cwidth; 
	cver  = info.cheight;
	uvx   = cx/2;        
	uvy  =  cy/2;            
	uvblk = xblk/2; 

	// the starting position, block size and dimension size at specific resolution 
	blksize = xblk>>s_level;
	uvblksize = blksize/2;
	start_x = cx>>s_level;
	start_y = cy>>s_level;
	s_yhor   = yhor>>s_level; 
	s_yver   = yver>>s_level; 
	s_chor   = chor>>s_level; 
	s_cver   = cver>>s_level; 
	start_uvx = uvx>>s_level; 
	start_uvy = uvy>>s_level; 

	if (cy-1>=0) // neighbor A is valid
	{
	    x_dest = cx;
		y_dest = cy-1;
		// get the block mode for neighbor A (first part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
					  &block_mode1, &propagate_iblk1, 0);
		x_dest = cx+4;  // since the basic block size is 4x4
		y_dest = cy-1; 
		// get the block mode for neighbor A (second part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
					 &block_mode2, &propagate_iblk2, 0);
		assert( block_mode1!= UNDEFINED && block_mode2!= UNDEFINED );  
		if (block_mode1==DIRECTIONAL_IBLOCK  || block_mode2==DIRECTIONAL_IBLOCK)
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK; 
			if (propagate_iblk1  || propagate_iblk2 ) 
				*fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}
		for (i=0; i<blksize; i++) 
		{
			neighborA[i] = (float)fr_cur[(start_y-1)*s_yhor+(start_x+i)];
			if (i%2==0  && fr_curU!=NULL && fr_curV!=NULL)
			{
				neighborUA[i/2] = (float)fr_curU[(start_uvy-1)*s_chor+(start_uvx+i/2)];
				neighborVA[i/2] = (float)fr_curV[(start_uvy-1)*s_chor+(start_uvx+i/2)];
			}
		}
	}

	 if (cy-1>=0  && cx+xblk+xblk-1<=yhor-1)   // neighbor B is in the frame
	 {
		x_dest = cx+xblk;
		y_dest = cy-1; 
		// get the block mode for neighbor B (first part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
					  &block_mode1, &propagate_iblk1, 0 );
		x_dest = cx+xblk+4; 
		y_dest = cy-1; 
		// get the block mode for neighbor B (second part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
					  &block_mode2, &propagate_iblk2, 0);
		if ( block_mode1 != UNDEFINED  && block_mode2 != UNDEFINED  && 
			(block_mode1==DIRECTIONAL_IBLOCK  || block_mode2==DIRECTIONAL_IBLOCK) )
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK;
			if ( propagate_iblk1  || propagate_iblk2  ) 
				*fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}
		// neighbor B has been visited before current block 
		if ( block_mode1 != UNDEFINED  && block_mode2 != UNDEFINED) 
			for (i=0; i<blksize; i++)
			{
				neighborB[i] = (float)fr_cur[(start_y-1)*s_yhor+(start_x+blksize+i)];
				if (i%2==0 && fr_curU!=NULL && fr_curV!=NULL)
				{
					neighborUB[i/2] = (float)fr_curU[(start_uvy-1)*s_chor+(start_uvx+uvblksize+i/2)];
					neighborVB[i/2] = (float)fr_curV[(start_uvy-1)*s_chor+(start_uvx+uvblksize+i/2)];
				}
			}
	 }


	if (cx-1>=0)  // neighbor G
	{
		x_dest = cx-1; 
		y_dest = cy; 
		// get the block mode for neighbor B (first part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
				      &block_mode1, &propagate_iblk1, 0);
		x_dest = cx-1;
		y_dest = cy+4; 
		// get the block mode for neighbor B (second part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
					  &block_mode2, &propagate_iblk2, 0);
		assert( block_mode1!= UNDEFINED && block_mode2!= UNDEFINED );  
		if (block_mode1==DIRECTIONAL_IBLOCK  || block_mode2==DIRECTIONAL_IBLOCK)
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK; 
			if ( propagate_iblk1  || propagate_iblk2  ) 
				*fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}
		for (i=0; i<blksize; i++)  // yblk instead of blksize
		{
			neighborG[i] = fr_cur[(start_y+i)*s_yhor+(start_x-1)];
			if (i%2==0 && fr_curU!=NULL && fr_curV!=NULL)
			{
				neighborUG[i/2] = (float)fr_curU[(start_uvy+i/2)*s_chor+(start_uvx-1)];
				neighborVG[i/2] = (float)fr_curV[(start_uvy+i/2)*s_chor+(start_uvx-1)];
			}
		}
	}


	if (cy-1>=0 && cx-1>=0)  // neighbor H 
	{
		neighborH  = fr_cur[(start_y-1)*s_yhor+(start_x-1)];
		x_dest = cx-1; 
		y_dest = cy-1; 
		// get the block mode for neighbor H
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
				      &block_mode1, &propagate_iblk1, 0);
		if (block_mode1==DIRECTIONAL_IBLOCK)
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK; 
			if ( propagate_iblk1) 
				*fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}
		if (fr_curU!=NULL && fr_curV!=NULL)
		{
			neighborUH = (float)fr_curU[(start_uvy-1)*s_chor+(start_uvx-1)];
			neighborVH = (float)fr_curV[(start_uvy-1)*s_chor+(start_uvx-1)];
		}
	}
	return; 
}


void GetSpaitialNeighborPixels4x4(videoinfo info, float *fadditional_penalty, vector_ptr fmv, 
								  float *fr_cur, int cx, int cy, int xblk, int yblk, int yhor, int yver, 
								  int t_level, float *fr_curU, float *fr_curV)
{
	int   i, chor, cver, uvx, uvy, uvblk, x_dest, y_dest; 
	float  mvx_d, mvy_d;
	enum  BiMode block_mode1;
	enum  FLAG is_pred;
	int   s_level, spatial_scale, blksize, uvblksize, start_x, start_y, s_yhor, s_yver; 
	int   s_chor, s_cver, start_uvx, start_uvy; 
	int   propagate_iblk1=0; 
	
	// the dimensions and starting positions for Y U V at full resolution 
	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
    if( s_level > 0 ) {
		info.ywidth  <<= s_level;
		info.yheight <<= s_level;
		info.cwidth  <<= s_level;	
		info.cheight <<= s_level;
		spatial_scale = 1<<s_level; 
	}
	chor  = info.cwidth; 
	cver  = info.cheight;
	uvx   = cx/2;        
	uvy   = cy/2;            
	uvblk = xblk/2; 

	// the starting position, block size and dimension size at specific resolution 
	blksize = xblk>>s_level;
	uvblksize = blksize/2;
	if (uvblksize==0)   // for the subsampling effect in spatial scalability 
		                // ( more the one leve resolution reduction )
	{
		if ( uvx%spatial_scale==0  && uvy%spatial_scale==0 )
			uvblksize = 1; 
	}
	start_x = cx>>s_level;
	start_y = cy>>s_level;
	s_yhor   = yhor>>s_level; 
	s_yver   = yver>>s_level; 
	s_chor   = chor>>s_level; 
	s_cver   = cver>>s_level; 
	start_uvx = uvx>>s_level; 
	start_uvy = uvy>>s_level; 

	if (cy-1>=0) // neighbor A is valid
	{
	    x_dest = cx;
		y_dest = cy-1;
		// get the block mode for neighbor A (first part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
			          &block_mode1, &propagate_iblk1, 0);
		assert( block_mode1!= UNDEFINED );  
		if (block_mode1==DIRECTIONAL_IBLOCK )
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK; 
			if (propagate_iblk1) *fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}
		for (i=0; i<blksize; i++) 
		{
			neighborA[i] = (float)fr_cur[(start_y-1)*s_yhor+(start_x+i)];
			if (i%2==0  && fr_curU!=NULL && fr_curV!=NULL  && uvblksize!=0 )
			{
				neighborUA[i/2] = (float)fr_curU[(start_uvy-1)*s_chor+(start_uvx+i/2)];
				neighborVA[i/2] = (float)fr_curV[(start_uvy-1)*s_chor+(start_uvx+i/2)];
			}
		}
	}

	 if (cy-1>=0  && cx+xblk+xblk-1<=yhor-1)   // neighbor B is in the frame
	 {
		x_dest = cx+xblk;
		y_dest = cy-1; 
		// get the block mode for neighbor B (first part)
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
			          &block_mode1, &propagate_iblk1, 0);
		if ( block_mode1 != UNDEFINED  && block_mode1==DIRECTIONAL_IBLOCK  )
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK;
			if (propagate_iblk1) *fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}

		// neighbor B has been visited before current block 
		if ( block_mode1 != UNDEFINED ) 
			for (i=0; i<blksize; i++)
			{
				neighborB[i] = (float)fr_cur[(start_y-1)*s_yhor+(start_x+blksize+i)];
				if (i%2==0 && fr_curU!=NULL && fr_curV!=NULL && uvblksize!=0)
				{
					neighborUB[i/2] = (float)fr_curU[(start_uvy-1)*s_chor+(start_uvx+uvblksize+i/2)];
					neighborVB[i/2] = (float)fr_curV[(start_uvy-1)*s_chor+(start_uvx+uvblksize+i/2)];
				}
			}
	 }


	if (cx-1>=0)  // neighbor G
	{
		x_dest = cx-1; 
		y_dest = cy; 
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
			&block_mode1, &propagate_iblk1, 0);
		assert( block_mode1!= UNDEFINED );  
		if (block_mode1==DIRECTIONAL_IBLOCK )
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK; 
			if (propagate_iblk1) *fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}
		for (i=0; i<blksize; i++)  // yblk instead of blksize
		{
			neighborG[i] = fr_cur[(start_y+i)*s_yhor+(start_x-1)];
			if (i%2==0 && fr_curU!=NULL && fr_curV!=NULL && uvblksize!=0)
			{
				neighborUG[i/2] = (float)fr_curU[(start_uvy+i/2)*s_chor+(start_uvx-1)];
				neighborVG[i/2] = (float)fr_curV[(start_uvy+i/2)*s_chor+(start_uvx-1)];
			}
		}
	}


	if (cy-1>=0 && cx-1>=0)  // neighbor H 
	{
		neighborH  = fr_cur[(start_y-1)*s_yhor+(start_x-1)];
		x_dest = cx-1; 
		y_dest = cy-1; 
		// get the block mode for neighbor H
		get_predictor(&mvx_d, &mvy_d, &is_pred, fmv, x_dest, y_dest, info, t_level, 
			&block_mode1, &propagate_iblk1, 0);
		if (block_mode1==DIRECTIONAL_IBLOCK)
		{
			*fadditional_penalty += DRIFT_PENALTY_IBLOCK; 
			if (propagate_iblk1) *fadditional_penalty += PROPAGATION_PENALTY_IBLOCK; 
		}
		if (fr_curU!=NULL && fr_curV!=NULL && uvblksize!=0)
		{
			neighborUH = (float)fr_curU[(start_uvy-1)*s_chor+(start_uvx-1)];
			neighborVH = (float)fr_curV[(start_uvy-1)*s_chor+(start_uvx-1)];
		}
	}

	return; 
}


void Spatial_prediction_Decision(videoinfo info, float *fadditional_penalty, vector_ptr fmv, 
								 float *fr_cur, int cx, int cy, int xblk, int yblk, int yhor, int yver, 
								 int t_level, float *fr_curU, float *fr_curV, vector *tmv1)
{
 
    int     spatial_mode;
	int     x, y, blkthresh;
	char    pre_sign; 
	float   pre_sum, pre_mse, diff; 

	tmv1->total_cost = tmv1->sad_cost = (float)HUGE_VAL; 
	blkthresh = xblk;
	for (spatial_mode=0; spatial_mode<PRE_MODE_NUM; spatial_mode++)
	{
		switch  (spatial_mode)
		{
			case SPATIAL_VERTICAL: // interpolate/predict vertically
				// we only need neighbor A and neighbor E for spatial vertical prediction
				// there are 4 cases for their existence
				pre_sign  = spatial_ver_pre(predict_blkY, neighborA, neighborE, xblk, xblk, blkthresh);
				break;

			case SPATIAL_HORIZONTAL: // interpolate/predict horizontally
				// we only need neighbor G and neighbor C for spatial horizontal prediction
				// there are 4 cases for their existence
				pre_sign = spatial_hor_pre(predict_blkY, neighborG, neighborC, xblk, xblk, blkthresh);
				break;

			case SPATIAL_DC:
				// we need neighbor A, C, G and E for spatial DC prediction
				// there are 16 cases for their existence
				pre_sign = spatial_dc_pre(predict_blkY, neighborA, neighborC, neighborE, 
										  neighborG, xblk, xblk, blkthresh);
				break;

			case SPATIAL_DOWN_LEFT:
				// we need neighbor A, G, B, F, C, E for spatial down-left prediction
				// the block is divided into 3 parts
				// each part has 4 cases, so totally there are 12 cases for their existence
				pre_sign = spatial_downleft_pre(predict_blkY, neighborA, neighborG,
					                            neighborB, neighborF,								
											    neighborC, neighborE, xblk, xblk, blkthresh);
				break;

			case SPATIAL_DOWN_RIGHT:
				// we need neighbor A, C, G, E, and H, D for spatial down-right prediction
                // the block is divided into 3 parts
				// each part has 4 cases, so totally there are 12 cases for their existence
				pre_sign = spatial_downright_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
				 								neighborH, neighborD, neighborD1, xblk, xblk, blkthresh);
				break;

			case SPATIAL_VERTICAL_RIGHT:
				// we need neighbor A, C, G, E, and H, D for spatial vertical-right prediction
				// for 2x2 block, there are only 2 parts: 8 cases
				// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
				pre_sign = spatial_verright_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
												neighborH, neighborD, xblk, xblk, blkthresh);
				break;

			case SPATIAL_HORIZONTAL_DOWN:
				// we need neighbor A, C, G, E, and H, D for spatial horizontal prediction
				// for 2x2 block, there are only 2 parts: 8 cases
				// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
				pre_sign = spatial_hordown_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
											neighborH, neighborD[0], neighborD1, xblk, xblk, blkthresh);
				break;
						
			case SPATIAL_VERTICAL_LEFT:
				// we need neighbor A, C, G, E, B, F for spatial vertical-left prediction
				// for 2x2 block there are only 2 parts: 8 cases
				// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
				pre_sign = spatial_verleft_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
					 						neighborB, neighborF[0], xblk, xblk, blkthresh);
				break;
			case SPATIAL_HORIZANTAL_UP:
				// we need neighbor A, C, G, E, B, F for spatial vertical-left prediction
				// for 2x2 block there are only 2 parts: 8 cases
				// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
				pre_sign = spatial_horup_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
											 neighborB[0], neighborF, xblk, xblk, blkthresh);
				break;
		} // switch  (spatial_mode)  

		if (pre_sign)  
		{
			pre_sum = pre_mse = 0.0; 
			for (y=0; y<xblk; y++)
			for (x=0; x<xblk; x++)
			{
				diff     = fr_cur[(cy+y)*yhor+cx+x] - predict_blkY[y*xblk+x];
				pre_sum  = (diff<0)?  pre_sum-diff: pre_sum+diff;
				pre_mse += diff*diff;
			}
		}else
			pre_sum = (float)HUGE_VAL;  // for invalid spatial mode 

		if ( pre_sum < tmv1->sad_cost)
		{
			tmv1->bi_mode = DIRECTIONAL_IBLOCK;
			tmv1->iblock_spatial_mode = (spatialMODE)spatial_mode; 
			tmv1->sad_cost   = pre_sum; 
			tmv1->bit_cost   = (*fadditional_penalty+DIRECTIONAL_IBLOCK_BIAS)*info.lambda[t_level]*3; 
			tmv1->total_cost = tmv1->sad_cost+tmv1->bit_cost;
		}
    } //for modes
}


void directional_iblock_detection(videoinfo info, vector_ptr fmv1_array, float *fr_cur,  
								  int cx, int cy, int xblock, int yblock, int hor, int ver, int t_level, 
								  vector *tmv1,  float *fadditional_penalty)
{
	int x, y, xblk, yblk, n_index; 

	xblk = ( cx + xblock <= hor ) ? xblock : hor - cx;
	yblk = ( cy + yblock <= ver ) ? yblock : ver - cy;

	if (xblk!=yblk || ( xblk!=8 && xblk!=4) )
		return; // we skip IBLOCK detection for non-squred blocks 

	// initialize all the neighbors to be HUGE_VAL
	neighborH = (float)HUGE_VAL; 
	for ( n_index=0; n_index<IBLOCK_MAX_SIZE; n_index++)
	{
		neighborA[n_index] = neighborB[n_index] = neighborC[n_index] = (float)HUGE_VAL; 
		neighborD[n_index]=  neighborD1[n_index]= neighborE[n_index] = (float)HUGE_VAL; 
		neighborF[n_index] = neighborG[n_index] = (float)HUGE_VAL; 
	}
	for (x=0; x<IBLOCK_MAX_SIZE; x++)
	for (y=0; y<IBLOCK_MAX_SIZE; y++)
		predict_blkY[y*IBLOCK_MAX_SIZE+x] = (float)HUGE_VAL; 

	if (xblk==8)
	{
		GetSpaitialNeighborPixels8x8(info, fadditional_penalty, fmv1_array, fr_cur, 
			                         cx, cy, xblk, yblk, hor, ver, t_level, NULL, NULL); 
		Spatial_prediction_Decision(info, fadditional_penalty, fmv1_array, fr_cur, 
									cx, cy, xblk, yblk, hor, ver, t_level, NULL, NULL, tmv1);
	}
	else if (xblk==4)
	{
		GetSpaitialNeighborPixels4x4(info, fadditional_penalty, fmv1_array, fr_cur, 
			                         cx, cy, xblk, yblk, hor, ver, t_level, NULL, NULL); 
		Spatial_prediction_Decision(info, fadditional_penalty, fmv1_array, fr_cur, 
									cx, cy, xblk, yblk, hor, ver, t_level, NULL, NULL, tmv1);
	}
	else
	{
		assert(0); 
		printf("Error in directional_iblock_detection()!\n"); 
	}
}

// the set of analysis functions is always working in full resolution (original resolution )
// in decoder (cx, cy) is the position at a specific resolution 
// hor and ver are the dimension at a specific resolution
void get_neighbor_predict_values(vector_ptr fmv, int cx, int cy, int xblock, int yblock, int hor, int ver, 
								 YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
					             YUVimage_ptr  fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
								 vector_ptr    fmv1, vector_ptr fmv2, vector_ptr fmv3,
								 vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
								 int level, int remaining_frs, videoinfo info,
								 ImageMEinfo *imagemeinfo, Varblkarrayinfo *varblkarray, int uvblksize)
{
	int   cpos, pos, row, col, i, j; 
	float self_weightY, self_weightU, self_weightV;
	float hor_weightY,  hor_weightU,  hor_weightV;
	float ver_weightY,  ver_weightU,  ver_weightV;
	float neighbor_valueY, neighbor_valueU, neighbor_valueV;
	int   leftx, topy, xblk, yblk, disx, disy; 
	int   s_level; 

	// the spatial resolution reduction in decoder 
	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));

	row = cy;  // y position of current pixel 
	for(i = 0; i < yblock; i++) {
		col = cx;  // x position of current pixel 
		for(j=0 ; j < xblock ; j++){
			// the position in Y coordinate with full resolution 
			cpos = ( row<<s_level )* ( hor<<s_level) + ( col<<s_level);  
			// information for current block in Y coordinate with full resolution 
			leftx     = imagemeinfo[cpos].leftx;   
			topy      = imagemeinfo[cpos].topy;
			xblk      = imagemeinfo[cpos].blksize;   // block size in full resolution 
			yblk      = imagemeinfo[cpos].blksize;
			// the distance from this pixel to the corner of the block in Y coordinate with full resolution
			disx = j<<s_level;   
			disy = i<<s_level;
			
			// weight coefficients for self motion vector and neighbor motion vectors 
			self_weightY = self_weightU = self_weightV = imagemeinfo[cpos].self_weight; 
			hor_weightY  = hor_weightU  = hor_weightV  = imagemeinfo[cpos].h_weight; 
			ver_weightY  = ver_weightU  = ver_weightV  = imagemeinfo[cpos].v_weight;
			assert(self_weightY+hor_weightY+ver_weightY==(float)1.0);
			neighbor_valueY = neighbor_valueU = neighbor_valueV = (float)0.0;

			// current pixel position ( col, row ) in Y coordinate with full resolution 
			// contributions from vertical neighbor motion vector
			if (disy < (yblk >> 1))  {     //  top neighbor is effective
				// the position of top neighbor at full resolution
				pos = (topy-3)*( hor<<s_level ) + ( col<<s_level );  
				if ( (topy > 0) && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK )
				{
					get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->Y, fr3->Y, &neighbor_valueY, 
											  &self_weightY, ver_weightY, col, row, info, hor, ver, 
											  0, level);
					if (i%2==0 && j%2==0  && uvblksize)
					{
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->U, fr3->U, &neighbor_valueU, 
												&self_weightU, ver_weightU, col/2, row/2, info, hor/2, ver/2, 
												1, level);
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->V, fr3->V, &neighbor_valueV, 
												&self_weightV, ver_weightV, col/2, row/2, info, hor/2, ver/2, 
												1, level);
					}
				}
			} else if (disy >= (yblk >> 1))  { // bottom neighbor is effective 
				pos = (topy+yblk+2)*( hor<<s_level ) + ( col<<s_level ); 
				if  ( (topy+yblk < ( ver<<s_level ) )  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK ) 
				{
					get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->Y, fr3->Y, &neighbor_valueY, 
											  &self_weightY, ver_weightY, col, row, info, hor, ver, 
											  0, level);
					if (i%2==0 && j%2==0  && uvblksize)
					{
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->U, fr3->U, &neighbor_valueU, 
												&self_weightU, ver_weightU, col/2, row/2, info, hor/2, ver/2, 
												1, level);
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->V, fr3->V, &neighbor_valueV, 
												&self_weightV, ver_weightV, col/2, row/2, info, hor/2, ver/2, 
												1, level);
					}
				}
			}

			// contribution from horizontal neighbor motion vectors
			if (disx < (xblk >> 1))  { // left neighbor is effective 
				pos = ( row<<s_level ) * ( hor<<s_level ) + leftx-3; 
				if ( (leftx > 0)  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK )
				{
					get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->Y, fr3->Y, &neighbor_valueY, 
											  &self_weightY, hor_weightY, col, row, info, hor, ver, 
											  0, level);
					if (i%2==0 && j%2==0  && uvblksize)
					{
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->U, fr3->U, &neighbor_valueU, 
												&self_weightU, hor_weightU, col/2, row/2, info, hor/2, ver/2, 
												1, level);
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->V, fr3->V, &neighbor_valueV, 
												&self_weightV, hor_weightV, col/2, row/2, info, hor/2, ver/2, 
												1, level);
					}
				}
			}else if (disx >= (xblk >> 1))  { // right neighbor is effective 
				pos = ( row<<s_level ) * ( hor<<s_level ) + leftx+xblk+2; 
				if ( (leftx+xblk < ( hor<<s_level ) )  && imagemeinfo[pos].bi_mode != DIRECTIONAL_IBLOCK )
				{
					get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->Y, fr3->Y, &neighbor_valueY, 
											  &self_weightY, hor_weightY, col, row, info, hor, ver, 
											  0, level);
					if (i%2==0 && j%2==0  && uvblksize)
					{
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->U, fr3->U, &neighbor_valueU, 
												&self_weightU, hor_weightU, col/2, row/2, info, hor/2, ver/2, 
												1, level);
						get_value_for_neighbor_mv(cpos, pos, imagemeinfo, fr1->V, fr3->V, &neighbor_valueV, 
												&self_weightV, hor_weightV, col/2, row/2, info, hor/2, ver/2, 
												1, level);
					}
				}
			}

			// neighbor_valueY, neighbor_valueU and neighbor_valueV have already been weighted by alpha
			// and neighbor weighting coefficients ( the set of neighbor values are negative! )
			neighbor_predictY[i*xblock+j] = neighbor_valueY;  
			self_weight_matrix[i*xblock+j]= self_weightY; 
			if ( i%2==0 && j%2==0  && uvblksize)
			{
				assert( self_weightU==self_weightV );
				self_weight_matrixUV[i/2*xblock/2+j/2] = self_weightV; 
				neighbor_predictU[i/2*xblock/2+j/2]    = neighbor_valueU; 
				neighbor_predictV[i/2*xblock/2+j/2]    = neighbor_valueV; 
			}
			col++;
	}
	row++;
 }
	

}



// the set of analysis functions is always working in full resolution (original resolution )
void Spatial_prediction_analysis_with_OBMC(vector_ptr fmv, int cx, int cy, int xblk, int yblk, int hor, int ver, 
								 YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
					             YUVimage_ptr  fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
								 vector_ptr    fmv1, vector_ptr fmv2, vector_ptr fmv3,
								 vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
								 int t_level, int remaining_frs, videoinfo info,
								 ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
{
    enum spatialMODE     spatial_mode;
	int     x, y, blkthresh;
	char    pre_sign, pre_signU, pre_signV; 
	float   diffY, diffU, diffV; 
	int     chor, cver;

	chor = hor/2;  cver = ver/2; 
	spatial_mode = fmv->iblock_spatial_mode;  // encoder side 
	blkthresh = xblk;  
	assert( spatial_mode>=0 && spatial_mode<=8 ); 
	switch  (spatial_mode)
	{
	case SPATIAL_VERTICAL: // interpolate/predict vertically
		// we only need neighbor A and neighbor E for spatial vertical prediction
		pre_sign  = spatial_ver_pre(predict_blkY, neighborA, neighborE, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_ver_pre(predict_blkU, neighborUA, neighborUE, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_ver_pre(predict_blkV, neighborVA, neighborVE, xblk/2, xblk/2, blkthresh);
		break;

	case SPATIAL_HORIZONTAL: // interpolate/predict horizontally
		// we only need neighbor G and neighbor C for spatial horizontal prediction
		pre_sign = spatial_hor_pre(predict_blkY, neighborG, neighborC, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_hor_pre(predict_blkU, neighborUG, neighborUC, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_hor_pre(predict_blkV, neighborVG, neighborVC, xblk/2, xblk/2, blkthresh);
		break;

	case SPATIAL_DC:
		// we need neighbor A, C, G and E for spatial DC prediction
		pre_sign = spatial_dc_pre(predict_blkY, neighborA, neighborC, neighborE, neighborG, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_dc_pre(predict_blkU, neighborUA, neighborUC, neighborUE, neighborUG, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_dc_pre(predict_blkV, neighborVA, neighborVC, neighborVE, neighborVG, xblk/2, xblk/2, blkthresh);
		break;

	case SPATIAL_DOWN_LEFT:
		// we need neighbor A, G, B, F, C, E for spatial down-left prediction
		// the block is divided into 3 parts
		// each part has 4 cases, so totally there are 12 cases for their existence
		pre_sign = spatial_downleft_pre(predict_blkY, neighborA, neighborG,
			                            neighborB, neighborF, neighborC, neighborE, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_downleft_pre(predict_blkU, neighborUA, neighborUG, 
	                              neighborUB, neighborUF, neighborUC, neighborUE, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_downleft_pre(predict_blkV, neighborVA, neighborVG, 
						          neighborVB, neighborVF, neighborVC, neighborVE, xblk/2, xblk/2, blkthresh);
		break;

	case SPATIAL_DOWN_RIGHT:
		// we need neighbor A, C, G, E, and H, D for spatial down-right prediction
        // the block is divided into 3 parts
		// each part has 4 cases, so totally there are 12 cases for their existence
		pre_sign  = spatial_downright_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
			 								neighborH, neighborD, neighborD1, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_downright_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUH, neighborUD, neighborUD1, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_downright_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
										neighborVH, neighborVD, neighborVD1, xblk/2, xblk/2, blkthresh);
		break;

	case SPATIAL_VERTICAL_RIGHT:
		// we need neighbor A, C, G, E, and H, D for spatial vertical-right prediction
		// for 2x2 block, there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_verright_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
										neighborH, neighborD, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_verright_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUH, neighborUD, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_verright_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
									neighborVH, neighborVD, xblk/2, xblk/2, blkthresh);
		break;

	case SPATIAL_HORIZONTAL_DOWN:
		// we need neighbor A, C, G, E, and H, D for spatial horizontal prediction
		// for 2x2 block, there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_hordown_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
										neighborH, neighborD[0], neighborD1, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_hordown_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUH, neighborUD[0], neighborUD1, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_hordown_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
										neighborVH, neighborVD[0], neighborVD1, xblk/2, xblk/2, blkthresh);
		break;
						
	case SPATIAL_VERTICAL_LEFT:
		// we need neighbor A, C, G, E, B, F for spatial vertical-left prediction
		// for 2x2 block there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_verleft_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
				 					neighborB, neighborF[0], xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_verleft_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUB, neighborUF[0], xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_verleft_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
										neighborVB, neighborVF[0], xblk/2, xblk/2, blkthresh);
		break;

	case SPATIAL_HORIZANTAL_UP:
		// we need neighbor A, C, G, E, B, F for spatial vertical-left prediction
		// for 2x2 block there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_horup_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
									 neighborB[0], neighborF, xblk, xblk, blkthresh);
		blkthresh = xblk/2;  // this is for format 4:2:0 
		pre_signU = spatial_horup_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
									  neighborUB[0], neighborUF, xblk/2, xblk/2, blkthresh);
		pre_signV = spatial_horup_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
									  neighborVB[0], neighborVF, xblk/2, xblk/2, blkthresh);
		break;
	} // switch  (spatial_mode)  

	if ( pre_sign )  
	{
		diffY = diffU = diffV = 0; 


		// for IBLOCK in OBMC framework 
		get_neighbor_predict_values(  fmv, cx, cy,  xblk, yblk, hor, ver, 
									  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									  mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, info,
									  frameMEinfo, varblkarray, xblk/2);

		for (y=0; y<xblk; y++)
		for (x=0; x<xblk; x++)
		{
			// alpha is the scaling factor for 5/3 filter normalization
			// neighbor_valueY, neighbor_valueU and neighbor_valueV have already been weighted by alpha
			// and neighbor weighting coefficients 
			H1->Y[(cy+y)*hor+(cx+x)] = alpha* (   fr2->Y[(cy+y)*hor+(cx+x)]
										        - predict_blkY[y*xblk+x]*self_weight_matrix[y*xblk+x] )
												+ neighbor_predictY[y*xblk+x];
			diffY += (float)fabs(H1->Y[(cy+y)*hor+(cx+x)]);
			assert( pre_signU && pre_signV ); 
			if ( x%2==0 && y%2==0 )
			{   // Note: we use sub-sampled version of self_weight_matrix, i.e. [y*xblk+x]
				H1->U[(cy+y)/2*chor+(cx+x)/2] = alpha * ( fr2->U[(cy+y)/2*chor+(cx+x)/2]
					                 -predict_blkU[y/2*xblk/2+x/2]*self_weight_matrixUV[y/2*xblk/2+x/2] )
									 +neighbor_predictU[y/2*xblk/2+x/2] ;
				diffU += (float)fabs(H1->U[(cy+y)/2*chor+(cx+x)/2]);
				H1->V[(cy+y)/2*chor+(cx+x)/2] = alpha * ( fr2->V[(cy+y)/2*chor+(cx+x)/2]
									-predict_blkV[y/2*xblk/2+x/2]*self_weight_matrixUV[y/2*xblk/2+x/2] )
									+neighbor_predictV[y/2*xblk/2+x/2] ;
				diffV += (float)fabs(H1->V[(cy+y)/2*chor+(cx+x)/2]);
			}
		}

#ifdef DEGUB_DIRECTIONAL_IBLOCK_OBMC
		FILE *fiblk = fopen("iblock_info_encoder.txt", "at");
		fprintf(fiblk, "x=%03d\t y=%03d\t blk=%02d\t %.2f\t %.2f\t %.2f\n", cx, cy, xblk, diffY, diffU, diffV); 
		fclose(fiblk); 
#endif 

		
	}else
		assert(0);   // in this function pre_sign must be valid
	
}


void diretional_iblock_spatial_interpolation_with_OBMC(vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
											 int hor, int ver, 
											 YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
											 YUVimage_ptr  fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
											 vector_ptr    fmv1, vector_ptr fmv2, vector_ptr fmv3,
											 vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
											 int t_level, int remaining_frs, videoinfo info, 
											 vector_ptr fmv_root,
											 ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
{
	int x, y, n_index; 
	float fadditional_penalty=0;

	// initialize all the neighbors to be HUGE_VAL
	neighborH = neighborUH = neighborVH = (float)HUGE_VAL; 
	for ( n_index=0; n_index<IBLOCK_MAX_SIZE; n_index++)
	{
		neighborA[n_index] = neighborB[n_index] = neighborC[n_index] = (float)HUGE_VAL; 
		neighborD[n_index]=  neighborD1[n_index]= neighborE[n_index] = (float)HUGE_VAL; 
		neighborF[n_index] = neighborG[n_index] = (float)HUGE_VAL; 
		if ( n_index%2==0)
		{
			neighborUA[n_index/2] = neighborUB[n_index/2] = neighborUC[n_index/2] = (float)HUGE_VAL; 
			neighborUD[n_index/2]=  neighborUD1[n_index/2]= neighborUE[n_index/2] = (float)HUGE_VAL; 
			neighborUF[n_index/2] = neighborUG[n_index/2] = (float)HUGE_VAL; 
			neighborVA[n_index/2] = neighborVB[n_index/2] = neighborVC[n_index/2] = (float)HUGE_VAL; 
			neighborVD[n_index/2]=  neighborVD1[n_index/2]= neighborVE[n_index/2] = (float)HUGE_VAL; 
			neighborVF[n_index/2] = neighborVG[n_index/2] = (float)HUGE_VAL; 
		}
	}
	for (x=0; x<IBLOCK_MAX_SIZE; x++)
	for (y=0; y<IBLOCK_MAX_SIZE; y++)
	{
		predict_blkY[y*IBLOCK_MAX_SIZE+x]      = (float)HUGE_VAL; 
		neighbor_predictY[y*IBLOCK_MAX_SIZE+x] = (float)HUGE_VAL; 
		self_weight_matrix[y*IBLOCK_MAX_SIZE+x]= (float)HUGE_VAL; 
		if ( x%2==0 && y%2==0)
		{
			predict_blkU[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
			predict_blkV[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
			neighbor_predictU[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
		    neighbor_predictV[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
			self_weight_matrixUV[y/2*IBLOCK_MAX_SIZE/2+x/2]= (float)HUGE_VAL; 
		}
	}

	if (xblk==8)
	{
		// fmv2 is the root of the motion vector quad-tree, fr2 is the current high temporal frame
		// all the other information is only for debug use
		GetSpaitialNeighborPixels8x8(info, &fadditional_penalty, fmv_root, fr2->Y,  
									 cx,  cy,  xblk,  yblk,  hor,  ver, t_level, fr2->U, fr2->V);

		Spatial_prediction_analysis_with_OBMC(  fmv, cx, cy,  xblk, yblk, hor, ver, 
									  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									  mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, info,
									  frameMEinfo, varblkarray);
	}
	else if (xblk==4)
	{
		GetSpaitialNeighborPixels4x4(info, &fadditional_penalty, fmv_root, fr2->Y, 
			                         cx, cy, xblk, yblk, hor, ver, t_level, fr2->U, fr2->V);
		Spatial_prediction_analysis_with_OBMC(  fmv, cx, cy,  xblk, yblk, hor, ver, 
									  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									  mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, info,
									  frameMEinfo, varblkarray);
	}
	else
	{
		assert(0); 
		printf("Error in diretional_iblock_spatial_interpolation()!\n"); 
	}
}


// this function is always working in full resolution ( original video resolution )
// fmv2 is the root of the non-NULL motion vector quad-tree (fmv)  递归调用，划分成最小的块
void  rec_directional_iblock_analysis_with_OBMC(vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
									  int hor, int ver, 
									YUVimage_ptr  L1, YUVimage_ptr  H1, YUVimage_ptr  H0, 
					                YUVimage_ptr  fr1, YUVimage_ptr fr2, YUVimage_ptr fr3,
									vector_ptr    fmv1, vector_ptr fmv2, vector_ptr fmv3,
									vector_ptr mv_ref1, vector_ptr mv_ref2, vector_ptr mv_ref3,
									int t_level, int remaining_frs, videoinfo info,
									vector_ptr fmv_root,
									ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
{
  int xblock, yblock;

  if( fmv->child ) {
    rec_directional_iblock_analysis_with_OBMC(  fmv->child0, cx,               cy,            xblk/2, yblk/2, hor, ver, 
									  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									  mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, info, fmv_root,
									  frameMEinfo, varblkarray);
    rec_directional_iblock_analysis_with_OBMC(  fmv->child1, cx + xblk / 2,     cy,           xblk/2, yblk/2, hor, ver, 
									  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									  mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, info, fmv_root,
									  frameMEinfo, varblkarray);
    rec_directional_iblock_analysis_with_OBMC(  fmv->child2, cx           ,     cy+ yblk / 2, xblk/2, yblk/2, hor, ver, 
									  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									  mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, info, fmv_root,
									  frameMEinfo, varblkarray);
    rec_directional_iblock_analysis_with_OBMC(  fmv->child3, cx + xblk / 2,     cy+ yblk / 2, xblk/2, yblk/2, hor, ver, 
									  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									  mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, info, fmv_root,
									  frameMEinfo, varblkarray);
  } else { // 最小的块
    /* consider the small block around the boundaries */
    xblock = ( cx + xblk <= hor ) ? xblk : hor - cx;
    yblock = ( cy + yblk <= ver ) ? yblk : ver - cy;

    if( xblock <= 0 || yblock <= 0 )		return;

	if (fmv->bi_mode != DIRECTIONAL_IBLOCK)	return; 

	// double check the validity of this directional IBLOCK 
	assert(fmv->lifting_mode == SPATIAL_PREDICTED ); 
	assert(xblock==yblock && (xblock==8 || xblock==4)); 
	
	// do the iblock spatial interpolation for this block 为这个块进行空间插值
	diretional_iblock_spatial_interpolation_with_OBMC(  fmv, cx,  cy,  xblock, yblock, hor, ver, 
											  L1,   H1,   H0, fr1,  fr2,  fr3,  fmv1,  fmv2,  fmv3, 
									          mv_ref1,  mv_ref2,  mv_ref3, t_level,  remaining_frs, 
											  info, fmv_root, frameMEinfo, varblkarray);
  }
}


void Spatial_prediction_synthesis_with_OBMC( vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
								   int hor, int ver, 
								   YUVimage_ptr fr0, YUVimage_ptr fr1, YUVimage H0,
								   YUVimage L1, YUVimage H1, YUVimage frp,
                                   vector_ptr fmv0, vector_ptr fmv1, vector_ptr fmv2,
                                   vector_ptr mv_ref0, vector_ptr mv_ref1, 
                                   vector_ptr mv_ref2, int level, videoinfo info, 
								   vector_ptr fmv_root,
								   ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
{
    enum spatialMODE     spatial_mode;
	int     x, y, blkthresh;
	char    pre_sign, pre_signU, pre_signV; 
	float   diffY, diffU, diffV; 
	int     chor, cver, s_level, uvx, uvy, uvblk, spatial_scale;
	int     blksize, uvblksize; 

	uvx   = cx/2;        
	uvy   = cy/2;            
	uvblk = xblk/2; 
	s_level = MY_MAX (0, info.s_level - (info.denoise_flag == YES));
	spatial_scale =1<<s_level; 

	// the starting position, block size and dimension size at specific resolution 
	blksize = xblk>>s_level;
	uvblksize = blksize/2;
	if (uvblksize==0)   // for the subsampling effect in spatial scalability 
	{
		if ( uvx%spatial_scale==0  && uvy%spatial_scale==0 )
			uvblksize = 1; 
	}
	// the starting position, block size and dimension size at specific resolution 
	cx   >>= s_level;
	cy   >>= s_level;
	hor  >>= s_level; 
	ver  >>= s_level; 
	chor = hor/2;  cver = ver/2; 

	spatial_mode = fmv->iblock_spatial_mode;  // encoder side 
	blkthresh = blksize;  
	assert( spatial_mode>=0 && spatial_mode<=8 && blksize<=8 ); 
	switch  (spatial_mode)
	{
	case SPATIAL_VERTICAL: // interpolate/predict vertically
		// we only need neighbor A and neighbor E for spatial vertical prediction
		pre_sign  = spatial_ver_pre(predict_blkY, neighborA, neighborE, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_ver_pre(predict_blkU, neighborUA, neighborUE, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_ver_pre(predict_blkV, neighborVA, neighborVE, uvblksize, uvblksize, blkthresh);
		break;

	case SPATIAL_HORIZONTAL: // interpolate/predict horizontally
		// we only need neighbor G and neighbor C for spatial horizontal prediction
		pre_sign = spatial_hor_pre(predict_blkY, neighborG, neighborC, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_hor_pre(predict_blkU, neighborUG, neighborUC, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_hor_pre(predict_blkV, neighborVG, neighborVC, uvblksize, uvblksize, blkthresh);
		break;

	case SPATIAL_DC:
		// we need neighbor A, C, G and E for spatial DC prediction
		pre_sign = spatial_dc_pre(predict_blkY, neighborA, neighborC, neighborE, neighborG, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_dc_pre(predict_blkU, neighborUA, neighborUC, neighborUE, neighborUG, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_dc_pre(predict_blkV, neighborVA, neighborVC, neighborVE, neighborVG, uvblksize, uvblksize, blkthresh);
		break;

	case SPATIAL_DOWN_LEFT:
		// we need neighbor A, G, B, F, C, E for spatial down-left prediction
		// the block is divided into 3 parts
		// each part has 4 cases, so totally there are 12 cases for their existence
		pre_sign = spatial_downleft_pre(predict_blkY, neighborA, neighborG,
			                            neighborB, neighborF, neighborC, neighborE, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_downleft_pre(predict_blkU, neighborUA, neighborUG, 
	                              neighborUB, neighborUF, neighborUC, neighborUE, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_downleft_pre(predict_blkV, neighborVA, neighborVG, 
						          neighborVB, neighborVF, neighborVC, neighborVE, uvblksize, uvblksize, blkthresh);
		break;

	case SPATIAL_DOWN_RIGHT:
		// we need neighbor A, C, G, E, and H, D for spatial down-right prediction
        // the block is divided into 3 parts
		// each part has 4 cases, so totally there are 12 cases for their existence
		pre_sign  = spatial_downright_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
			 								neighborH, neighborD, neighborD1, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_downright_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUH, neighborUD, neighborUD1, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_downright_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
										neighborVH, neighborVD, neighborVD1, uvblksize, uvblksize, blkthresh);
		break;

	case SPATIAL_VERTICAL_RIGHT:
		// we need neighbor A, C, G, E, and H, D for spatial vertical-right prediction
		// for 2x2 block, there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_verright_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
										neighborH, neighborD, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_verright_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUH, neighborUD, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_verright_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
									neighborVH, neighborVD, uvblksize, uvblksize, blkthresh);
		break;

	case SPATIAL_HORIZONTAL_DOWN:
		// we need neighbor A, C, G, E, and H, D for spatial horizontal prediction
		// for 2x2 block, there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_hordown_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
										neighborH, neighborD[0], neighborD1, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_hordown_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUH, neighborUD[0], neighborUD1, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_hordown_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
										neighborVH, neighborVD[0], neighborVD1, uvblksize, uvblksize, blkthresh);
		break;
						
	case SPATIAL_VERTICAL_LEFT:
		// we need neighbor A, C, G, E, B, F for spatial vertical-left prediction
		// for 2x2 block there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_verleft_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
				 					neighborB, neighborF[0], blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_verleft_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
										neighborUB, neighborUF[0], uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_verleft_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
										neighborVB, neighborVF[0], uvblksize, uvblksize, blkthresh);
		break;

	case SPATIAL_HORIZANTAL_UP:
		// we need neighbor A, C, G, E, B, F for spatial vertical-left prediction
		// for 2x2 block there are only 2 parts: 8 cases
		// for 4x4 and 8x8 blocks, there are 5 parts: 20 cases
		pre_sign  = spatial_horup_pre(predict_blkY, neighborA, neighborC, neighborG, neighborE,
									 neighborB[0], neighborF, blksize, blksize, blkthresh);
		blkthresh = uvblksize;  // this is for format 4:2:0 
		pre_signU = spatial_horup_pre(predict_blkU, neighborUA, neighborUC, neighborUG, neighborUE,
									  neighborUB[0], neighborUF, uvblksize, uvblksize, blkthresh);
		pre_signV = spatial_horup_pre(predict_blkV, neighborVA, neighborVC, neighborVG, neighborVE,
									  neighborVB[0], neighborVF, uvblksize, uvblksize, blkthresh);
		break;
	} // switch  (spatial_mode)  

	if ( pre_sign )  
	{
		diffY = diffU = diffV = 0; 


		// for IBLOCK in OBMC framework 
		get_neighbor_predict_values(  fmv, cx, cy,  blksize, blksize, hor, ver, 
									  &L1,   &H1,   &H0, &frp,  fr0,  fr1,  fmv0,  fmv1,  fmv2, 
									  mv_ref0,  mv_ref1,  mv_ref2, level,  0, info,
									  frameMEinfo, varblkarray, uvblksize);

		for (y=0; y<blksize; y++)
		for (x=0; x<blksize; x++)
		{
			// alpha is the scaling factor for 5/3 filter normalization
			fr0->Y[(cy+y)*hor+(cx+x)] = 1/alpha* ( H0.Y[(cy+y)*hor+(cx+x)]- neighbor_predictY[y*blksize+x] )
				                        +predict_blkY[y*blksize+x]*self_weight_matrix[y*blksize+x] ;
			diffY += (float)fabs(H0.Y[(cy+y)*hor+(cx+x)]);
			if ( x%2==0 && y%2==0 && pre_signU && pre_signV  && uvblksize )
			{
				fr0->U[(cy+y)/2*chor+(cx+x)/2] = 1/alpha * ( H0.U[(cy+y)/2*chor+(cx+x)/2]-neighbor_predictU[y/2*uvblksize+x/2])
					                            +predict_blkU[y/2*uvblksize+x/2]* self_weight_matrixUV[y/2*uvblksize+x/2];
				diffU += (float)fabs(H0.U[(cy+y)/2*chor+(cx+x)/2]);
				fr0->V[(cy+y)/2*chor+(cx+x)/2] = 1/alpha * ( H0.V[(cy+y)/2*chor+(cx+x)/2]-neighbor_predictV[y/2*uvblksize+x/2])
					                            +predict_blkV[y/2*uvblksize+x/2]* self_weight_matrixUV[y/2*uvblksize+x/2];
				diffV += (float)fabs(H1.V[(cy+y)/2*chor+(cx+x)/2]);
			}
		}

//		printf("x=%03d\t y=%03d\t blk=%02d\t %.2f\t %.2f\t %.2f\n", cx, cy, xblk, diffY, diffU, diffV); 
	}else
		assert(0);   // in this function pre_sign must be valid
	
}


void diretional_iblock_spatial_synthesis_with_OBMC( vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
									   int hor, int ver, 
									   YUVimage_ptr fr0, YUVimage_ptr fr1, YUVimage H0,
									   YUVimage L1, YUVimage H1, YUVimage frp,
                                       vector_ptr fmv0, vector_ptr fmv1, vector_ptr fmv2,
                                       vector_ptr mv_ref0, vector_ptr mv_ref1, 
                                       vector_ptr mv_ref2, int level, videoinfo info, 
									   vector_ptr fmv_root,
									   ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
{
	int x, y, n_index; 
	float fadditional_penalty=0;

	// initialize all the neighbors to be HUGE_VAL
	neighborH = neighborUH = neighborVH = (float)HUGE_VAL; 
	for ( n_index=0; n_index<IBLOCK_MAX_SIZE; n_index++)
	{
		neighborA[n_index] = neighborB[n_index] = neighborC[n_index] = (float)HUGE_VAL; 
		neighborD[n_index]=  neighborD1[n_index]= neighborE[n_index] = (float)HUGE_VAL; 
		neighborF[n_index] = neighborG[n_index] = (float)HUGE_VAL; 
		if ( n_index%2==0)
		{
			neighborUA[n_index/2] = neighborUB[n_index/2] = neighborUC[n_index/2] = (float)HUGE_VAL; 
			neighborUD[n_index/2]=  neighborUD1[n_index/2]= neighborUE[n_index/2] = (float)HUGE_VAL; 
			neighborUF[n_index/2] = neighborUG[n_index/2] = (float)HUGE_VAL; 
			neighborVA[n_index/2] = neighborVB[n_index/2] = neighborVC[n_index/2] = (float)HUGE_VAL; 
			neighborVD[n_index/2]=  neighborVD1[n_index/2]= neighborVE[n_index/2] = (float)HUGE_VAL; 
			neighborVF[n_index/2] = neighborVG[n_index/2] = (float)HUGE_VAL; 
		}
	}

	for (x=0; x<IBLOCK_MAX_SIZE; x++)
	for (y=0; y<IBLOCK_MAX_SIZE; y++)
	{
		predict_blkY[y*IBLOCK_MAX_SIZE+x]      = (float)HUGE_VAL; 
		neighbor_predictY[y*IBLOCK_MAX_SIZE+x] = (float)HUGE_VAL; 
		self_weight_matrix[y*IBLOCK_MAX_SIZE+x]= (float)HUGE_VAL; 
		if ( x%2==0 && y%2==0)
		{
			predict_blkU[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
			predict_blkV[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
			neighbor_predictU[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
		    neighbor_predictV[y/2*IBLOCK_MAX_SIZE/2+x/2] = (float)HUGE_VAL; 
			self_weight_matrixUV[y/2*IBLOCK_MAX_SIZE/2+x/2]= (float)HUGE_VAL; 
		}
	}

	if (xblk==8)
	{
		GetSpaitialNeighborPixels8x8(info, &fadditional_penalty, fmv_root, fr0->Y,  
									 cx,  cy,  xblk,  yblk,  hor,  ver, level, fr0->U, fr0->V);
		Spatial_prediction_synthesis_with_OBMC(   fmv,  cx,  cy,  xblk,  yblk, 
									    hor,  ver, fr0,  fr1,  H0, L1,  H1,  frp, fmv0,  fmv1,  fmv2,
                                        mv_ref0,  mv_ref1,  mv_ref2,  level,  info, fmv_root,
										frameMEinfo, varblkarray);
		
	}
	else if (xblk==4)
	{
		GetSpaitialNeighborPixels4x4(info, &fadditional_penalty, fmv_root, fr0->Y, 
			                         cx, cy, xblk, yblk, hor, ver, level, fr0->U, fr0->V);
		Spatial_prediction_synthesis_with_OBMC(   fmv,  cx,  cy,  xblk,  yblk, 
									    hor,  ver, fr0,  fr1,  H0, L1,  H1,  frp, fmv0,  fmv1,  fmv2,
                                        mv_ref0,  mv_ref1,  mv_ref2,  level,  info, fmv_root,
										frameMEinfo, varblkarray);
	}
	else
	{
		assert(0); 
		printf("Error in diretional_iblock_spatial_interpolation()!\n"); 
	}
}


void rec_directional_iblock_synthesis_with_OBMC( vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
									   int hor, int ver, 
									   YUVimage_ptr fr0, YUVimage_ptr fr1, YUVimage H0,
									   YUVimage L1, YUVimage H1, YUVimage frp,
                                       vector_ptr fmv0, vector_ptr fmv1, vector_ptr fmv2,
                                       vector_ptr mv_ref0, vector_ptr mv_ref1, 
                                       vector_ptr mv_ref2, int level, videoinfo info, 
									   vector_ptr fmv_root,
									   ImageMEinfo *frameMEinfo, Varblkarrayinfo *varblkarray)
{
  int xblock, yblock;

  if( fmv->child ) {
	rec_directional_iblock_synthesis_with_OBMC(  fmv->child0,  cx,         cy      ,  xblk/2,  yblk/2, hor,  ver, 							    
		                               fr0,  fr1,  H0,  L1,  H1,  frp,  
									   fmv0,  fmv1,  fmv2,  mv_ref0,  mv_ref1, mv_ref2,  
									   level,  info,  fmv_root,
									   frameMEinfo, varblkarray);
	rec_directional_iblock_synthesis_with_OBMC(  fmv->child1,  cx+xblk/2,  cy      ,  xblk/2,  yblk/2, hor,  ver, 							    
		                               fr0,  fr1,  H0,  L1,  H1,  frp,  
									   fmv0,  fmv1,  fmv2,  mv_ref0,  mv_ref1, mv_ref2,  
									   level,  info,  fmv_root,
									   frameMEinfo, varblkarray);
	rec_directional_iblock_synthesis_with_OBMC(  fmv->child2,  cx       ,  cy+yblk/2,  xblk/2,  yblk/2, hor,  ver, 							    
		                               fr0,  fr1,  H0,  L1,  H1,  frp,  
									   fmv0,  fmv1,  fmv2,  mv_ref0,  mv_ref1, mv_ref2,  
									   level,  info,  fmv_root,
									   frameMEinfo, varblkarray);
	rec_directional_iblock_synthesis_with_OBMC(  fmv->child3,  cx+xblk/2 ,  cy+yblk/2,  xblk/2,  yblk/2, hor,  ver, 							    
		                               fr0,  fr1,  H0,  L1,  H1,  frp,  
									   fmv0,  fmv1,  fmv2,  mv_ref0,  mv_ref1, mv_ref2,  
									   level,  info,  fmv_root,
									   frameMEinfo, varblkarray);
  } else {
    /* consider the small block around the boundaries */
    xblock = ( cx + xblk <= hor ) ? xblk : hor - cx;
    yblock = ( cy + yblk <= ver ) ? yblk : ver - cy;

    if( xblock <= 0 || yblock <= 0 )		return;

	if (fmv->bi_mode != DIRECTIONAL_IBLOCK)	return; 

	// double check the validity of this directional IBLOCK 
	assert(fmv->lifting_mode == SPATIAL_PREDICTED ); 
	assert(xblock==yblock && (xblock==8 || xblock==4)); 

	// do the iblock spatial synthesis for this block 
	diretional_iblock_spatial_synthesis_with_OBMC(  fmv,  cx,  cy,  xblk,  yblk, 
									    hor,  ver, fr0,  fr1,  H0, L1,  H1,  frp, fmv0,  fmv1,  fmv2,
                                        mv_ref0,  mv_ref1,  mv_ref2,  level,  info, fmv_root,
										frameMEinfo, varblkarray);
  }


}

