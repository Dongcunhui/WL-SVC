
#include "stdio.h"
#include "stdlib.h"
#include "structN.h"
#include <assert.h>

#define SELECT_MERGE_BLOCK 32
#define MV_MERGE_THRESH    7
#define MERGENCE_DEBUG

inline float median(float a, float b, float c)
{
  return ((a>b) ? ((a>c) ? ((b>c) ? b : c) : a) : 
          ((b>c) ? ((a>c) ? a : c) : b));
}


void
replace_enhancement_mv_subtree( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count, float sample_mvx, float sample_mvy )
{
	int cx, cy; 

	if (!fmv1->child)
	{
		if ((fmv1->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv1->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
		fmv1->mvx = sample_mvx;  fmv1->mvy = sample_mvy; 
	}else
	{
	    cx = x;
		cy = y;
		replace_enhancement_mv_subtree( fmv1_array, fmv2_array, fmv1->child0, fmv2? fmv2->child0:NULL,
										pmvx, pmvy, num_symbol, subpel,
										cx,  cy,  xblk/2,  yblk/2,  hor,  ver,
										info,  decode_parallelmv,  t_level,  bidir_exist, 
										blk_thresh,  count,  sample_mvx,  sample_mvy );
	    cx = x + xblk / 2;
		cy = y;
		replace_enhancement_mv_subtree( fmv1_array, fmv2_array, fmv1->child1, fmv2? fmv2->child1:NULL,
										pmvx, pmvy, num_symbol, subpel,
										cx,  cy,  xblk/2,  yblk/2,  hor,  ver,
										info,  decode_parallelmv,  t_level,  bidir_exist, 
										blk_thresh,  count,  sample_mvx,  sample_mvy );
	    cx = x;
		cy = y + yblk / 2;
		replace_enhancement_mv_subtree( fmv1_array, fmv2_array, fmv1->child2, fmv2? fmv2->child2:NULL,
										pmvx, pmvy, num_symbol, subpel,
										cx,  cy,  xblk/2,  yblk/2,  hor,  ver,
										info,  decode_parallelmv,  t_level,  bidir_exist, 
										blk_thresh,  count,  sample_mvx,  sample_mvy );
		cx = x + xblk / 2;
		cy = y + yblk / 2;
		replace_enhancement_mv_subtree( fmv1_array, fmv2_array, fmv1->child3, fmv2? fmv2->child3:NULL,
										pmvx, pmvy, num_symbol, subpel,
										cx,  cy,  xblk/2,  yblk/2,  hor,  ver,
										info,  decode_parallelmv,  t_level,  bidir_exist, 
										blk_thresh,  count,  sample_mvx,  sample_mvy );
	}
}


void
replace_enhancement_mv_further( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count )
{
  int cx, cy; 

  if( fmv1->child && xblk>blk_thresh ) {     
    cx = x;
    cy = y;
    replace_enhancement_mv_further(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    replace_enhancement_mv_further(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    replace_enhancement_mv_further(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    replace_enhancement_mv_further(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);
  } else {
	  if( x >= hor || y >= ver )     return;

	  if (!fmv1->child)   		     return;

	  if (!fmv1->mv_exist)           return; 

	  // irregular motion area
	  if ( fmv1->child0->child || fmv1->child1->child ||
		   fmv1->child2->child || fmv1->child3->child)
		   return;

	  if ( fmv1->child && fmv1->mv_exist && fmv1->merge_sign )
		  replace_enhancement_mv_subtree(fmv1_array, fmv2_array, fmv1, fmv2 ,
                    pmvx, pmvy, num_symbol, subpel, x, y, xblk, yblk,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count, 
					fmv1->sample_mvx, fmv1->sample_mvy);

  }
}


void
replace_enhancement_mv( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count )
{
  int cx, cy; 

  if( fmv1->child && xblk>blk_thresh ) {     
    cx = x;
    cy = y;
    replace_enhancement_mv(fmv1_array, fmv2_array, fmv1->child0, fmv2 ? fmv2->child0 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y;
    replace_enhancement_mv(fmv1_array, fmv2_array, fmv1->child1, fmv2 ? fmv2->child1 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x;
    cy = y + yblk / 2;
    replace_enhancement_mv(fmv1_array, fmv2_array, fmv1->child2, fmv2 ? fmv2->child2 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);

    cx = x + xblk / 2;
    cy = y + yblk / 2;
    replace_enhancement_mv(fmv1_array, fmv2_array, fmv1->child3, fmv2 ? fmv2->child3 : NULL,
                    pmvx, pmvy, num_symbol, subpel, cx, cy, xblk / 2, yblk / 2,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count);
  } else {
	  if( x >= hor || y >= ver )     return;

	  if (!fmv1->child)   		 return;

	  if (!fmv1->mv_exist)       return; 

	  if ( fmv1->child && fmv1->mv_exist && fmv1->merge_sign )
		  replace_enhancement_mv_subtree(fmv1_array, fmv2_array, fmv1, fmv2 ,
                    pmvx, pmvy, num_symbol, subpel, x, y, xblk, yblk,
                    hor, ver, info, decode_parallelmv, t_level, bidir_exist, blk_thresh, count, 
					fmv1->sample_mvx, fmv1->sample_mvy);

	  if (fmv1->child && fmv1->mv_exist && !fmv1->merge_sign)     
	  {
		  replace_enhancement_mv_further(  fmv1_array,  fmv2_array,  fmv1,  fmv2,
                 pmvx, pmvy,  num_symbol,  subpel, x,  y,  xblk,  yblk,  hor,  ver,
                  info,  decode_parallelmv,  t_level,  bidir_exist, blk_thresh/2,  count );

	  }

  }

}


// 就是进行了一些设置fmv->mv_exist、 fmv->sample_mvx、 fmv->sample_mvy
// recursively visit the child, do the subsampling  递归访问孩子，做下采样
// we subsmaple for blocks with size ( blk_thresh x blk_thresh )
void scalable_quad_tree_child(vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
							  int hor, int ver, int blk_thresh)
{
  // blkszie>blk_thresh and have children 有孩子且块大于阈值
  if(fmv->child  && yblk>blk_thresh && xblk>blk_thresh){  // there are children and blksize>blk_thresh
    scalable_quad_tree_child(fmv->child0, cx,        cy,        xblk/2, yblk/2, hor, ver, blk_thresh);
    scalable_quad_tree_child(fmv->child1, cx+xblk/2, cy,        xblk/2, yblk/2, hor, ver, blk_thresh);
    scalable_quad_tree_child(fmv->child2, cx,        cy+yblk/2, xblk/2, yblk/2, hor, ver, blk_thresh);
    scalable_quad_tree_child(fmv->child3, cx+xblk/2, cy+yblk/2, xblk/2, yblk/2, hor, ver, blk_thresh);
  }
  else{
		if (cx >= hor || cy >= ver) return; // boundary check 边界返回 

		// blksize>blk_thresh and have no children 块大但是没有孩子，返回
		if  (!(fmv->child) && yblk>blk_thresh  && xblk>blk_thresh ) return; 

		// blksize==blk_thresh and no children 相同大小但是没有孩子
		if  (!(fmv->child) && yblk==blk_thresh && xblk==blk_thresh ){
			fmv->mv_exist = 0;
			// there is motion vector for this block on this side 
			if ( fmv->lifting_mode != IGNORED && fmv->bi_mode != DIRECTIONAL_IBLOCK)
			{
				fmv->mv_exist = 1;
				fmv->sample_mvx  = fmv->mvx;
				fmv->sample_mvy  = fmv->mvy;
			}
			return;
		}
			          
		// blksize==blk_thresh and have children 剩下的最后一种情况 块大小相同但是有孩子
		fmv->mv_exist = 0;
		// child0:
		if ( !(cx >= hor || cy >= ver) ) // boundary check for the child0
		{
			// there is motion vector on this side for child0
			if ( fmv->child0->mv_exist )
			{
				fmv->mv_exist = 1;
				fmv->sample_mvx   = fmv->child0->sample_mvx;
				fmv->sample_mvy   = fmv->child0->sample_mvy;
				return;
			}
			
		}

		// child1:
		if ( !(cx+xblk/2>=hor || cy >=ver) ) // boundary check for the child1
		{
			// there is motion vector on this side for child1
			if ( fmv->child1->mv_exist )
			{
				fmv->mv_exist = 1;
				fmv->sample_mvx   = fmv->child1->sample_mvx;
				fmv->sample_mvy   = fmv->child1->sample_mvy;
				return;
			}
		}

		// child2:
		if ( !(cx >= hor || cy+yblk/2>ver) ) // boundary check for the child2
		{
			// there is motion vector on this side for child0
			if ( fmv->child2->mv_exist )
			{
				fmv->mv_exist = 1;
				fmv->sample_mvx   = fmv->child2->sample_mvx;
				fmv->sample_mvy   = fmv->child2->sample_mvy;
				return;
			}
		}

		// child3:
		if ( !(cx+xblk/2>=hor || cy+yblk/2>=ver) ) // boundary check for the child3
		{
			// there is motion vector on this side for child0
			if ( fmv->child3->mv_exist )
			{
				fmv->mv_exist = 1;
				fmv->sample_mvx   = fmv->child3->sample_mvx;
				fmv->sample_mvy   = fmv->child3->sample_mvy;
				return;
			}
		}
	}
}


void whether_merge_this_sub_tree(vector_ptr fmv, int cx, int cy, int xblk, int yblk, int hor, int ver, 
								 int *merge_sign, float sample_mvx, float sample_mvy, 
								 int *sub_blocks, float scale_factor, 
								 int tlevel, int tindex, int blk_thresh)
{
	if (fmv->child)
	{
		whether_merge_this_sub_tree(fmv->child0, cx,        cy,        xblk/2, yblk/2, hor, ver, merge_sign, 
			                        sample_mvx, sample_mvy, sub_blocks, scale_factor, 
									tlevel, tindex, blk_thresh);
		whether_merge_this_sub_tree(fmv->child1, cx+xblk/2, cy,        xblk/2, yblk/2, hor, ver, merge_sign, 
			                        sample_mvx, sample_mvy, sub_blocks, scale_factor,
									tlevel, tindex, blk_thresh);
		whether_merge_this_sub_tree(fmv->child2, cx       , cy+yblk/2, xblk/2, yblk/2, hor, ver, merge_sign, 
			                        sample_mvx, sample_mvy, sub_blocks, scale_factor,
									tlevel, tindex, blk_thresh);
		whether_merge_this_sub_tree(fmv->child3, cx+xblk/2, cy+yblk/2, xblk/2, yblk/2, hor, ver, merge_sign, 
			                        sample_mvx, sample_mvy, sub_blocks, scale_factor,
									tlevel, tindex, blk_thresh);
	}
	else
	{
		if (cx>=hor  || cy>=ver) 	return;
	
		(*sub_blocks)++;

		if ((fmv->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 
		
#ifdef  MERGENCE_DEBUG
		FILE *fsub_tree; 
		char  subtree_file[80];
		sprintf(subtree_file, "subtree%d_i%d_j%d.txt", blk_thresh, tlevel, tindex);
		fsub_tree = fopen(subtree_file, "at"); 
		fprintf(fsub_tree, "x=%d\t y=%d blk=%d\t mvx=%.2f\t mvy=%.2f\n",
			    cx, cy, xblk, fmv->mvx, fmv->mvy);
		fclose(fsub_tree); 
#endif

		if (fabs(fmv->mvx - sample_mvx)>MV_MERGE_THRESH*scale_factor  || 
			fabs(fmv->mvy - sample_mvy)>MV_MERGE_THRESH*scale_factor )
			*merge_sign = 0;
		return;
	}
}

// selectively merge quad-treee
void select_merge_quad_tree_child_further(vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
										  int hor, int ver, int *merge_block_num2, int *checked_block2, 
										  int blk_thresh, int tlevel, int tindex)
{
	int merge_sign, sub_blocks;
	
    // if a child is NULL, then there is nothing to do 
	// there is children and blksize>blk_thresh
	if(fmv->child  && yblk>blk_thresh && xblk>blk_thresh){  
		select_merge_quad_tree_child_further(fmv->child0, cx,        cy,        xblk/2, yblk/2, hor, ver, 
									merge_block_num2, checked_block2, blk_thresh, tlevel, tindex);
		select_merge_quad_tree_child_further(fmv->child1, cx+xblk/2, cy,        xblk/2, yblk/2, hor, ver, 
									merge_block_num2, checked_block2, blk_thresh, tlevel, tindex);
		select_merge_quad_tree_child_further(fmv->child2, cx,        cy+yblk/2, xblk/2, yblk/2, hor, ver, 
									merge_block_num2, checked_block2, blk_thresh, tlevel, tindex);
		select_merge_quad_tree_child_further(fmv->child3, cx+xblk/2, cy+yblk/2, xblk/2, yblk/2, hor, ver, 
									merge_block_num2, checked_block2, blk_thresh, tlevel, tindex);
	}
	else{
		if (cx >= hor || cy >= ver) 	return;

		// blksize>=blk_threh but no children
		if  (!(fmv->child))  return;          // motion vectors in base layer 

#ifdef  MERGENCE_DEBUG
		FILE *fsub_tree; 
		char  subtree_file[80];
		sprintf(subtree_file, "subtree%d_i%d_j%d.txt", blk_thresh, tlevel, tindex);
		fsub_tree = fopen(subtree_file, "at"); 
		fprintf(fsub_tree, "checking subtree x=%d\t y=%d blk=%d\n",  cx, cy, xblk);
		fclose(fsub_tree); 
#endif 

		// irregular motion area
		if ( fmv->child0->child || fmv->child1->child ||
			 fmv->child2->child || fmv->child3->child)
			 return;

		(*checked_block2)++;
		merge_sign = 1;
		if (fmv->mv_exist)  // check whether the block is allowed to do subsampling for motion vectors
		{
			sub_blocks  = 0;
			whether_merge_this_sub_tree(fmv, cx, cy, xblk, yblk, hor, ver, &merge_sign, 
										fmv->sample_mvx, fmv->sample_mvy, &sub_blocks, 1.5, 
										tlevel, tindex, blk_thresh);
		}else
		{
#ifdef  MERGENCE_DEBUG
			sprintf(subtree_file, "subtree%d_i%d_j%d.txt", blk_thresh, tlevel, tindex);
			fsub_tree = fopen(subtree_file, "at"); 
			fprintf(fsub_tree, "NO motion vector in this subtree\n\n\n");
			fclose(fsub_tree); 
#endif 
			fmv->merge_sign = 0;
			return;
		}

#ifdef  MERGENCE_DEBUG
		sprintf(subtree_file, "subtree%d_i%d_j%d.txt", blk_thresh, tlevel, tindex);
		fsub_tree = fopen(subtree_file, "at"); 
		if (merge_sign)
			fprintf(fsub_tree, "subtree MERGED\n");
		else
			fprintf(fsub_tree, "subtree NOT MERGED\n");
		fprintf(fsub_tree, "\n\n"); 
		fclose(fsub_tree);
#endif		

		fmv->merge_sign = merge_sign;  // the final decision whether merge this sub-tree
		if (merge_sign)  
		   (*merge_block_num2)++;

	}
}


// selectively merge quad treee
void select_merge_quad_tree_child(vector_ptr fmv, int cx, int cy, int xblk, int yblk, int hor, int ver, 
								  int *merge_block_num, int *merge_block_num2,int *checked_block1, 
								  int *checked_block2, int blk_thresh, int t_level, int tindex)
{
	int merge_sign, sub_blocks;
	
	// there is children and blksize>blk_thresh
	if(fmv->child  && yblk>blk_thresh && xblk>blk_thresh){  
		select_merge_quad_tree_child(fmv->child0, cx,        cy,        xblk/2, yblk/2, hor, ver, 
									merge_block_num, merge_block_num2,  checked_block1, checked_block2, 
									blk_thresh, t_level, tindex);
		select_merge_quad_tree_child(fmv->child1, cx+xblk/2, cy,        xblk/2, yblk/2, hor, ver, 
									merge_block_num, merge_block_num2, checked_block1, checked_block2,  
									blk_thresh, t_level, tindex);
		select_merge_quad_tree_child(fmv->child2, cx,        cy+yblk/2, xblk/2, yblk/2, hor, ver, 
									merge_block_num, merge_block_num2, checked_block1, checked_block2, 
									blk_thresh, t_level, tindex);
		select_merge_quad_tree_child(fmv->child3, cx+xblk/2, cy+yblk/2, xblk/2, yblk/2, hor, ver, 
									merge_block_num, merge_block_num2, checked_block1, checked_block2, 
									blk_thresh, t_level, tindex);
	}
	else{
		if (cx >= hor || cy >= ver) 	return;

		// blksize>=blk_threh but no children
		if  (!(fmv->child))  return;          // motion vectors in base layer 

#ifdef  MERGENCE_DEBUG
		FILE *fsub_tree; 
		char  subtree_file[80];
		sprintf(subtree_file, "subtree%d_i%d_j%d.txt", blk_thresh, t_level, tindex);
		fsub_tree = fopen(subtree_file, "at"); 
		fprintf(fsub_tree, "checking subtree x=%d\t y=%d blk=%d\n",  cx, cy, xblk);
		fclose(fsub_tree); 
#endif 

		merge_sign = 1;
		if (fmv->mv_exist)  // check whether the block is allowed to do subsampling 
		{
			sub_blocks  = 0;
			whether_merge_this_sub_tree(fmv, cx, cy, xblk, yblk, hor, ver, &merge_sign, 
										fmv->sample_mvx, fmv->sample_mvy, &sub_blocks, 1.0,
										t_level, tindex, blk_thresh);
		}else
		{
#ifdef  MERGENCE_DEBUG
			sprintf(subtree_file, "subtree%d_i%d_j%d.txt", blk_thresh, t_level, tindex);
			fsub_tree = fopen(subtree_file, "at"); 
			fprintf(fsub_tree, "NO motion vector in this subtree\n\n\n");
			fclose(fsub_tree); 
#endif 
			fmv->merge_sign = 0;
			return; 
		}

#ifdef  MERGENCE_DEBUG
		sprintf(subtree_file, "subtree%d_i%d_j%d.txt", blk_thresh, t_level, tindex);
		fsub_tree = fopen(subtree_file, "at"); 
		if (merge_sign)
			fprintf(fsub_tree, "subtree MERGED\n");
		else
			fprintf(fsub_tree, "subtree NOT MERGED\n");
		fprintf(fsub_tree, "\n\n"); 
		fclose(fsub_tree);
#endif		

		fmv->merge_sign = merge_sign;  // the final decision whether merge this sub-tree
		if (merge_sign)  
		{  // for statistical use
		   *merge_block_num += sub_blocks;
		   (*checked_block1)++;  // merged blk_thresh x blk_thresh blocks 
		}else
		{
			select_merge_quad_tree_child_further(fmv, cx, cy, xblk, yblk, hor, ver, 
								  merge_block_num2, checked_block2, blk_thresh/2, t_level, tindex);

		}

	}
}


// motion vector subsample from 4x4->8x8->16x16->32x32 
// there is no redundunt motion vector for coding 编码不需要冗余的运动向量 就是设置了一些变量
void layer_structure_mv_subsample(videoinfo info, int t_level, int tindex, vector_ptr fmv, int encoder_side)
{

    int i, x, y, X, Y, xnum, ynum, xblk, yblk, hor, ver, blk_thresh;
    int merge_block_num, merge_block_num2, checked_block1, checked_block2, total_blk_thresh;

	if (!info.layer_mv[t_level])  return; // no layer structure for this temporal level
		
	xnum = info.xnum[t_level];
	ynum = info.ynum[t_level];
	xblk = info.xblk[t_level];
	yblk = info.yblk[t_level];
	hor  = info.ywidth;
	ver  = info.yheight;

	// recursively subsample motion vectors from bottom to top, i.e. 4x4->8x8->16x16->32x32 递归的
	for (i=0; i<4; i++)
	{
		blk_thresh = 4*(1<<i); // 块大小 从小到大
		for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {   
			for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
				scalable_quad_tree_child(fmv+Y*xnum+X, x, y, xblk, yblk,  hor, ver, blk_thresh); 
			}
		}
	}

	if (encoder_side)
	{
		// do the selective merge 
		// for 32x32 blocks we check the difference among the children motion vectors
		// if the 32x32 block can not be merged, check the four 16x16 block
		merge_block_num2 = merge_block_num = checked_block1 = checked_block2 = 0;  
		blk_thresh = SELECT_MERGE_BLOCK; 
		for(y=0, Y=0 ; Y< ynum ; y+= yblk, Y++){
			for(x=0, X=0 ; X< xnum ; x+= xblk, X++){
				select_merge_quad_tree_child(fmv+Y*xnum+X, x, y, xblk, yblk, hor, ver, 
											 &merge_block_num, &merge_block_num2, &checked_block1, 
											 &checked_block2,  blk_thresh, t_level, tindex); 
			}
		}

		total_blk_thresh = (info.yheight/blk_thresh * info.ywidth/blk_thresh);
		printf("total merged small blocks: %03d \n",  merge_block_num);
		printf("total merged %dx%d blocks: %03d, percentage: %f\n",  blk_thresh, blk_thresh, checked_block1,
				checked_block1/(float)total_blk_thresh);
		printf("total merged %dx%d blocks: %03d, percentage: %f \n", blk_thresh/2, blk_thresh/2, 
			merge_block_num2, checked_block2? merge_block_num2/(float)checked_block2:0 );
	}	

}

void recursive_layer_structure_mv_trim(vector_ptr fmv, int cx, int cy, int xblk, int yblk, 
							  int hor, int ver)
{
	int xblk2, yblk2; 

  // blkszie>blk_thresh and have children
  if(fmv->child){  // there are children and blksize>blk_thresh
    recursive_layer_structure_mv_trim(fmv->child0, cx,        cy,        xblk/2, yblk/2, hor, ver);
    recursive_layer_structure_mv_trim(fmv->child1, cx+xblk/2, cy,        xblk/2, yblk/2, hor, ver);
    recursive_layer_structure_mv_trim(fmv->child2, cx,        cy+yblk/2, xblk/2, yblk/2, hor, ver);
    recursive_layer_structure_mv_trim(fmv->child3, cx+xblk/2, cy+yblk/2, xblk/2, yblk/2, hor, ver);
  }
  else{
		if (cx >= hor || cy >= ver) return; // boundary check 

		if ((fmv->lifting_mode == IGNORED))   return;
#ifdef  DIRECTIONAL_IBLOCK_EMPLOYED
		if (fmv->bi_mode==DIRECTIONAL_IBLOCK) return; 
#endif 

        // trim the replaced motion vector, and make it in the frame boundary
		xblk2 = ( cx + xblk <= hor ) ? xblk : hor - cx;
		yblk2 = ( cy + yblk <= ver ) ? yblk : ver - cy;
		if (cx + xblk2 - fmv->mvx > hor )
			fmv->mvx = (float) ( hor - (cx+xblk2)); 
		if ( cx-fmv->mvx<0)
			fmv->mvx = (float) cx; 
		if ( cy+yblk2-fmv->mvy > ver )
			fmv->mvy = (float) (ver- (cy+yblk2));
		if ( cy-fmv->mvy<0)
			fmv->mvy = (float) cy; 

  }
}


void layer_structure_mv_trim(videoinfo info, int t_level, int tindex, vector_ptr fmv, int count)
{

    int x, y, X, Y, xnum, ynum, xblk, yblk, hor, ver;
		
	xnum = info.xnum[t_level];
	ynum = info.ynum[t_level];
	xblk = info.xblk[t_level];
	yblk = info.yblk[t_level];
	hor  = info.ywidth;
	ver  = info.yheight;

	// recursively trim the replaced motion vectors
	for( y = 0, Y = 0; Y < ynum; y += yblk, Y++ ) {   
		for( x = 0, X = 0; X < xnum; x += xblk, X++ ) {
			recursive_layer_structure_mv_trim(fmv+Y*xnum+X, x, y, xblk, yblk,  hor, ver); 
		}
	}
}


void  clear_frame_motion_field(FRAME_MOTION_FIELD *frame_motion_field, videoinfo info)
{
    int i, j; 
	for (i=0; i<info.ywidth; i++)
		for (j=0; j<info.yheight; j++)
		{
			frame_motion_field[j*info.ywidth+i].available = 0;
			frame_motion_field[j*info.ywidth+i].dmvx      = (float)HUGE_VAL;
			frame_motion_field[j*info.ywidth+i].dmvy      = (float)HUGE_VAL;
			frame_motion_field[j*info.ywidth+i].mvx       = (float)HUGE_VAL;
			frame_motion_field[j*info.ywidth+i].mvy       = (float)HUGE_VAL;
		}
}

void  clear_frame_motion_field_simp(SIMP_FRAME_MOTION_FIELD *frame_motion_field, videoinfo info)
{
    int i, j; 
	for (i=0; i<info.ywidth; i++)
		for (j=0; j<info.yheight; j++)
		{
			frame_motion_field[j*info.ywidth+i].available = 0;
			frame_motion_field[j*info.ywidth+i].mvx       = (float)HUGE_VAL;
			frame_motion_field[j*info.ywidth+i].mvy       = (float)HUGE_VAL;
		}
}


void update_frame_motion_field(FRAME_MOTION_FIELD *frame_motion_field, int x, int y, 
							   int xblk, int yblk, videoinfo info, vector_ptr fmv, float dmvx, float dmvy,
							   float mvx, float mvy)
{
	int i, j, xblk2, yblk2; 

	int subpel = 2;

	float accu = 0.25;
	int addx,addy;
	float int_mvx,int_mvy;

	xblk2 = ( x + xblk <= info.ywidth ) ?  xblk : info.ywidth - x;
	yblk2 = ( y + yblk <= info.yheight ) ? yblk : info.yheight - y;

	for (i=0; i<xblk2; i++)
		for (j=0; j<yblk2; j++)
		{
			if( (fmv->bi_mode <= 6 && fmv->bi_mode >= 0) || fmv->bi_mode == 8 || (fmv->bi_mode == 7 && fmv->aff_mrg == NO) ){
				frame_motion_field[(y+j)*info.ywidth+(x+i)].available = 1;
				frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = dmvx;
				frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = dmvy;
				frame_motion_field[(y+j)*info.ywidth+(x+i)].mvx       = mvx;
				frame_motion_field[(y+j)*info.ywidth+(x+i)].mvy       = mvy;
			}
			else if( (fmv->bi_mode <= 11 && fmv->bi_mode >= 9) || (fmv->bi_mode == 7 && fmv->aff_mrg == YES) ){

				if( (i+1)/xblk2 == 1 )
					addx = 1;
				else
					addx = 0;

				if( (j+1)/yblk2 == 1 )
					addy = 1;
				else
					addy = 0;

				frame_motion_field[(y+j)*info.ywidth+(x+i)].available = 1;
				frame_motion_field[(y+j)*info.ywidth+(x+i)].mvx       = (fmv->aff2_mvx - fmv->aff1_mvx)*(i + addx)/xblk2 + (fmv->aff3_mvx - fmv->aff1_mvx)*(j + addy)/yblk2 + fmv->aff1_mvx;
				frame_motion_field[(y+j)*info.ywidth+(x+i)].mvy       = (fmv->aff2_mvy - fmv->aff1_mvy)*(i + addx)/xblk2 + (fmv->aff3_mvy - fmv->aff1_mvy)*(j + addy)/yblk2 + fmv->aff1_mvy;
				
				if(fmv->bi_mode == 7 && fmv->aff_mrg == YES){
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = 0;
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = 0;
				}else if(fmv->direct_idx == DIRECT ){
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = 0;
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = 0;
				}
				else if(fmv->direct_idx == INDIRECT && fmv->merge_idx == MERGE && fmv->merge_dir == UP){
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = (0 - 0)*((float)i)/((float)xblk2) + (fmv->aff3_dmvx - 0)*((float)j)/((float)yblk2) + 0;
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = (0 - 0)*((float)i)/((float)xblk2) + (fmv->aff3_dmvy - 0)*((float)j)/((float)yblk2) + 0;
				}
				else if(fmv->direct_idx == INDIRECT && fmv->merge_idx == MERGE && fmv->merge_dir == LEFT){
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = (fmv->aff2_dmvx - 0)*((float)i)/((float)xblk2) + (0 - 0)*((float)j)/((float)yblk2) + 0;
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = (fmv->aff2_dmvy - 0)*((float)i)/((float)xblk2) + (0 - 0)*((float)j)/((float)yblk2) + 0;
				}
				else if(fmv->direct_idx == INDIRECT && fmv->merge_idx == INTER){
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = (fmv->aff2_dmvx - fmv->aff1_dmvx)*((float)i)/((float)xblk2) + (fmv->aff3_dmvx - fmv->aff1_dmvx)*((float)j)/((float)yblk2) + fmv->aff1_dmvx;
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = (fmv->aff2_dmvy - fmv->aff1_dmvy)*((float)i)/((float)xblk2) + (fmv->aff3_dmvy - fmv->aff1_dmvy)*((float)j)/((float)yblk2) + fmv->aff1_dmvy;
				}
				else if(fmv->direct_idx == INDIRECT && fmv->merge_idx == MERGE && fmv->merge_dir == PAL_L ){
					if(fmv->aff_idx >= 0){
						frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = (fmv->aff2_dmvx - fmv->aff1_dmvx)*((float)i)/((float)xblk2) + (fmv->aff3_dmvx - fmv->aff1_dmvx)*((float)j)/((float)yblk2) + fmv->aff1_dmvx;
						frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = (fmv->aff2_dmvy - fmv->aff1_dmvy)*((float)i)/((float)xblk2) + (fmv->aff3_dmvy - fmv->aff1_dmvy)*((float)j)/((float)yblk2) + fmv->aff1_dmvy;
					}else{
						frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx      = 0;
						frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy      = 0;
					}
				}
				else if(fmv->direct_idx == INDIRECT && fmv->merge_idx == MERGE && fmv->merge_dir == TRAN_P ){
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvx	= (fmv->aff2_dmvx - fmv->aff1_dmvx)*((float)i)/((float)xblk2) + (fmv->aff3_dmvx - fmv->aff1_dmvx)*((float)j)/((float)yblk2) + fmv->aff1_dmvx;
					frame_motion_field[(y+j)*info.ywidth+(x+i)].dmvy	= (fmv->aff2_dmvy - fmv->aff1_dmvy)*((float)i)/((float)xblk2) + (fmv->aff3_dmvy - fmv->aff1_dmvy)*((float)j)/((float)yblk2) + fmv->aff1_dmvy;
				}else
					assert(0);
			}else
				assert(0);
		}
}

///////////// Added on 02.20.2016 ////////////////////
void copy_frame_motion_field( FRAME_MOTION_FIELD *src_frame_motion_field, SIMP_FRAME_MOTION_FIELD *dest_frame_motion_field, videoinfo info ){
	int i,j;

	for (i=0; i<info.ywidth; i++)
		for (j=0; j<info.yheight; j++)
		{
			dest_frame_motion_field[j*info.ywidth+i].available = src_frame_motion_field[j*info.ywidth+i].available;
			dest_frame_motion_field[j*info.ywidth+i].mvx       = src_frame_motion_field[j*info.ywidth+i].mvx;
			dest_frame_motion_field[j*info.ywidth+i].mvy       = src_frame_motion_field[j*info.ywidth+i].mvy;
		}
}
//////////////////////////////////////////////////////

///////////// Added on 02.20.2016 ////////////////////
void copy_frame_motion_field2( SIMP_FRAME_MOTION_FIELD *src_frame_motion_field, SIMP_FRAME_MOTION_FIELD *dest_frame_motion_field, videoinfo info ){
	int i,j;

	for (i=0; i<info.ywidth; i++)
		for (j=0; j<info.yheight; j++)
		{
			dest_frame_motion_field[j*info.ywidth+i].available = src_frame_motion_field[j*info.ywidth+i].available;
			dest_frame_motion_field[j*info.ywidth+i].mvx       = src_frame_motion_field[j*info.ywidth+i].mvx;
			dest_frame_motion_field[j*info.ywidth+i].mvy       = src_frame_motion_field[j*info.ywidth+i].mvy;
		}
}
//////////////////////////////////////////////////////

// get median predictor for block E, using blocks A, B, C, and D
// block scheme:  CBD
//                AN
//				  E
void get_median_predictor_motion_field(float *pmvx, float *pmvy,
									   FRAME_MOTION_FIELD *frame_motion_field, SIMP_FRAME_MOTION_FIELD *prev_frame_motion_field,
									  int x_pos, int y_pos, int xblock, int yblock,
			                          videoinfo info, int t_level, int blk_thresh)
{  //Modified by Yuan Liu on 01.12.2016 regarding the subpel accuracy part of median prediction
   //Modified by Yuan Liu on 01.23.2016, median prediction is replaced by neighbor prediction competition.
  int hor, ver, x_dest, y_dest,i,k;
  int is_pred_a, is_pred_b, is_pred_c, is_pred_d, is_pred_e;
  float mvx_a, mvy_a;
  float mvx_b, mvy_b;
  float mvx_c, mvy_c;
  float mvx_d, mvy_d;
  float mvx_e, mvy_e;

  float mvx_tmp, mvy_tmp;
  int is_pred_tmp;

  int blocks_avail,count; 

  int it_mvx,it_mvy, subpel = 2;

  float med_px[4], med_py[4];
  float int_mvx,int_mvy,accu = 0.25;

  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }
  
  // initialization
  hor = info.ywidth;
  ver = info.yheight;
  
  xblock = ( x_pos + xblock <= hor ) ? xblock : hor - x_pos;
  yblock = ( y_pos + yblock <= ver ) ? yblock : ver - y_pos;

  assert(x_pos >= 0 && x_pos < hor && y_pos >= 0 && y_pos < ver);

  //Modified by Yuan Liu on 01.12.2016 regarding the subpel accuracy part of median prediction

  // try to get predictors for blocks A, B, and C
  // A
  x_dest = x_pos - 1;
  y_dest = y_pos + yblock - 1;
  if (x_dest < 0 || y_dest >= ver) {
    is_pred_a = 0;
  } else {
	is_pred_a = frame_motion_field[y_dest*hor+x_dest].available;
	if ( is_pred_a )
	{
		mvx_a = frame_motion_field[y_dest*hor+x_dest].mvx;
		mvy_a = frame_motion_field[y_dest*hor+x_dest].mvy;

		mvx_a = mvx_a * (1 << subpel);
		mvy_a = mvy_a * (1 << subpel);

		it_mvx = (int)(mvx_a);
		it_mvy = (int)(mvy_a);

		mvx_a = (float)(it_mvx);
		mvy_a = (float)(it_mvy);

		mvx_a = mvx_a / (1 << subpel);
		mvy_a = mvy_a / (1 << subpel);

	}
  }

  // B
  y_dest = y_pos - 1;
  x_dest = x_pos + xblock - 1;
  if (y_dest < 0 || x_dest >= hor) {
    is_pred_b = 0;
  } else {
	is_pred_b = frame_motion_field[y_dest*hor+x_dest].available;
	if ( is_pred_b )
	{
		mvx_b = frame_motion_field[y_dest*hor+x_dest].mvx;
		mvy_b = frame_motion_field[y_dest*hor+x_dest].mvy;

		mvx_b = mvx_b * (1 << subpel);
		mvy_b = mvy_b * (1 << subpel);

		it_mvx = (int)(mvx_b);
		it_mvy = (int)(mvy_b);

		mvx_b = (float)(it_mvx);
		mvy_b = (float)(it_mvy);

		mvx_b = mvx_b / (1 << subpel);
		mvy_b = mvy_b / (1 << subpel);

	}
  }

  // C
  x_dest = x_pos - 1;
  y_dest = y_pos - 1;
  if (x_dest < 0 || y_dest < 0) {
    is_pred_c = 0;
  } else {
	is_pred_c = frame_motion_field[y_dest*hor+x_dest].available;
	if ( is_pred_c )
	{
		mvx_c = frame_motion_field[y_dest*hor+x_dest].mvx;
		mvy_c = frame_motion_field[y_dest*hor+x_dest].mvy;

		mvx_c = mvx_c * (1 << subpel);
		mvy_c = mvy_c * (1 << subpel);

		it_mvx = (int)(mvx_c);
		it_mvy = (int)(mvy_c);

		mvx_c = (float)(it_mvx);
		mvy_c = (float)(it_mvy);

		mvx_c = mvx_c / (1 << subpel);
		mvy_c = mvy_c / (1 << subpel);
	}
  }

  // blocks A, B, and C available?
  if (is_pred_a == 1 && is_pred_b == 1 && is_pred_c == 1) {
	pmvx[0] = mvx_a;
	pmvy[0] = mvy_a;
	pmvx[1] = mvx_b;
	pmvy[1] = mvy_b;
	pmvx[2] = mvx_c;
	pmvy[2] = mvy_c;
  } else {
// try to get predictor for block D
// D
	x_dest = x_pos + xblock;
    y_dest = y_pos - 1;
    if (x_dest >= hor || y_dest < 0) {
      is_pred_d = 0;
    } else {
		is_pred_d = frame_motion_field[y_dest*hor+x_dest].available;
		if ( is_pred_d )
		{
			mvx_d = frame_motion_field[y_dest*hor+x_dest].mvx;
			mvy_d = frame_motion_field[y_dest*hor+x_dest].mvy;

			mvx_d = mvx_d * (1 << subpel);
			mvy_d = mvy_d * (1 << subpel);

			it_mvx = (int)(mvx_d);
			it_mvy = (int)(mvy_d);

			mvx_d = (float)(it_mvx);
			mvy_d = (float)(it_mvy);

			mvx_d = mvx_d / (1 << subpel);
			mvy_d = mvy_d / (1 << subpel);

		}
    }

////////////////////	Added by Yuan Liu	////////////////////
	if(is_pred_a == 1){
		pmvx[0] = mvx_a;
		pmvy[0] = mvy_a;
	}
	if(is_pred_b == 1){
		pmvx[1] = mvx_b;
		pmvy[1] = mvy_b;
	}
	if(is_pred_c == 1){
		pmvx[2] = mvx_c;
		pmvy[2] = mvy_c;
	}
	if(is_pred_d == 1){
		pmvx[3] = mvx_d;
		pmvy[3] = mvy_d;
	}
////////////////////////////////////////////////////////////////
  }

  // check if predictor points outside the reference frame
  // only use integer part of the prediction motion vector for checking
  // which means sub-pixel distance such as 0.25, 0.50 and 0.75 beyond the 
  // boundary is allowed. 
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
      pmvx[i] = (float)HUGE_VAL;
      pmvy[i] = (float)HUGE_VAL;
    }
  }

  //////////////  Added on 02.02.2016  /////////////////////////
  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 3);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }
  assert(pmvx[3] == (float)HUGE_VAL && pmvy[3] == (float)HUGE_VAL);
//////////////////////////////////////////////////////////////
  //Set median prediction candidate if available
  if(pmvx[0] != (float)HUGE_VAL && pmvy[0] != (float)HUGE_VAL && pmvx[1] != (float)HUGE_VAL && pmvy[1] != (float)HUGE_VAL
	  && pmvx[2] != (float)HUGE_VAL && pmvy[2] != (float)HUGE_VAL){
	pmvx[3] = median( pmvx[0],pmvx[1],pmvx[2]);
	pmvy[3] = median( pmvy[0],pmvy[1],pmvy[2]);
  }

  //////////////////////////////////////////////////////////////
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
      pmvx[i] = (float)HUGE_VAL;
      pmvy[i] = (float)HUGE_VAL;
    }
  }

  //////////////  Added on 02.21.2016  /////////////////////////
  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 4);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }

//  if(pmvx[3]==(float)HUGE_VAL && pmvy[3]==(float)HUGE_VAL && pmvx[2]!=(float)HUGE_VAL && pmvy[2]!=(float)HUGE_VAL){
//	  assert( pmvx[0]!=(float)HUGE_VAL && pmvy[0]!=(float)HUGE_VAL && pmvx[1]!=(float)HUGE_VAL && pmvy[1]!=(float)HUGE_VAL );
  for(count = 0; count <= 3; count++){
    if(pmvx[count] == (float)HUGE_VAL && pmvy[count] == (float)HUGE_VAL){
	  x_dest = x_pos;
	  y_dest = y_pos;
	  if (x_dest < 0 || y_dest < 0) {
		assert(0);
	  } else {
		  is_pred_tmp = prev_frame_motion_field[y_dest*hor+x_dest].available;
		  if(is_pred_tmp == 1){
			mvx_tmp = prev_frame_motion_field[y_dest*hor+x_dest].mvx;
			mvy_tmp = prev_frame_motion_field[y_dest*hor+x_dest].mvy;

			mvx_tmp = mvx_tmp * (1 << subpel);
			mvy_tmp = mvy_tmp * (1 << subpel);

			it_mvx = (int)(mvx_tmp);
			it_mvy = (int)(mvy_tmp);

			mvx_tmp = (float)(it_mvx);
			mvy_tmp = (float)(it_mvy);

			mvx_tmp = mvx_tmp / (1 << subpel);
			mvy_tmp = mvy_tmp / (1 << subpel);

			if(mvx_tmp != (float)HUGE_VAL && mvy_tmp != (float)HUGE_VAL){
				pmvx[count] = mvx_tmp;
				pmvy[count] = mvy_tmp;
			}
		  }
	  }
	  break;
	}
  }

////////////////////////////////////////////////////////////////////
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
      pmvx[i] = (float)HUGE_VAL;
      pmvy[i] = (float)HUGE_VAL;
    }
  }

  //////////////  Added on 02.27.2016  /////////////////////////
  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 4);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }
//////////////////////////////////////////////////////////////
  for(count = 0; count <= 3; count++){
    if(pmvx[count] == (float)HUGE_VAL && pmvy[count] == (float)HUGE_VAL){
		x_dest = x_pos - 1;
		y_dest = y_pos + yblock;
		if (x_dest < 0 || y_dest >= ver) {
		  is_pred_e = 0;
		} else {
			is_pred_e = frame_motion_field[y_dest*hor+x_dest].available;
			if ( is_pred_e == 1 )
			{
				mvx_e = frame_motion_field[y_dest*hor+x_dest].mvx;
				mvy_e = frame_motion_field[y_dest*hor+x_dest].mvy;

				mvx_e = mvx_e * (1 << subpel);
				mvy_e = mvy_e * (1 << subpel);

				it_mvx = (int)(mvx_e);
				it_mvy = (int)(mvy_e);

				mvx_e = (float)(it_mvx);
				mvy_e = (float)(it_mvy);

				mvx_e = mvx_e / (1 << subpel);
				mvy_e = mvy_e / (1 << subpel);

				if(mvx_e != (float)HUGE_VAL && mvy_e != (float)HUGE_VAL){
					pmvx[count] = mvx_e;
					pmvy[count] = mvy_e;
				}
			}
		}
		break;
	}
  }

////////////////////////////////////////////////////////////////////
  for(i = 0; i <= 3; i++){
    if (x_pos - (int)(pmvx[i]) < 0 || x_pos - (int)(pmvx[i]) + xblock > hor  || 
        y_pos - (int)(pmvy[i]) < 0 || y_pos - (int)(pmvy[i]) + yblock > ver) {
      pmvx[i] = (float)HUGE_VAL;
      pmvy[i] = (float)HUGE_VAL;
    }
  }

  //////////////  Added on 02.27.2016  /////////////////////////
  //check for repeated candidates
  for(i = 3;i >= 1;i --){
	  for(k = i-1;k >= 0;k --){
		if(pmvx[i] == pmvx[k] && pmvy[i] == pmvy[k] && pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
			pmvx[i] = (float)HUGE_VAL;
			pmvy[i] = (float)HUGE_VAL;
		}
	  }
  }
  //check for HUGE_VALs
  k = 0;
  for(i = 0;i <= 3;i++){
	  if(pmvx[i] != (float)HUGE_VAL && pmvy[i] != (float)HUGE_VAL){
		med_px[k] = pmvx[i];
		med_py[k] = pmvy[i];
		k++;
	  }
  }
  assert(k <= 4);
  for(i = 0;i <= 3;i++){
	pmvx[i] = med_px[i];
	pmvy[i] = med_py[i];
  }
  for(i=0;i<=3;i++){
	med_px[i] = (float)HUGE_VAL;
	med_py[i] = (float)HUGE_VAL;
  }

}


void ec_get_contexts_motion_field(int *ctx_x, int *ctx_y, FRAME_MOTION_FIELD *frame_motion_field,
								  int x, int y, videoinfo info, int t_level, int blk_thresh)
{
  float dmvx_hor, dmvy_hor, dmvx_ver, dmvy_ver, dmvx_dia, dmvy_dia;
  int ispred_hor, ispred_ver, ispred_dia;
  float e_x, e_y;

  if (x-1>=0)
  {
	  ispred_hor = frame_motion_field[y*info.ywidth+(x-1)].available;
	  if ( ispred_hor )
	  {
		  dmvx_hor = frame_motion_field[y*info.ywidth+(x-1)].dmvx;
		  dmvy_hor = frame_motion_field[y*info.ywidth+(x-1)].dmvy;
	  }
  }else
	  ispred_hor = 0; 

  if (y-1>=0)
  {
	  ispred_ver = frame_motion_field[(y-1)*info.ywidth+x].available;
	  if ( ispred_ver )
	  {
		  dmvx_ver = frame_motion_field[(y-1)*info.ywidth+x].dmvx;
		  dmvy_ver = frame_motion_field[(y-1)*info.ywidth+x].dmvy;
	  }
  }else
	ispred_ver = 0; 

  if ( x-1>=0 && y-1>=0)
  {
	  ispred_dia = frame_motion_field[(y-1)*info.ywidth+(x-1)].available;
	  if ( ispred_dia )
	  {
		  dmvx_dia = frame_motion_field[(y-1)*info.ywidth+(x-1)].dmvx;
		  dmvy_dia = frame_motion_field[(y-1)*info.ywidth+(x-1)].dmvy;
	  }
  }else
	ispred_dia = 0; 

  if(ispred_hor == 0) dmvx_hor = dmvy_hor = 0.;
  if(ispred_ver == 0) dmvx_ver = dmvy_ver = 0.;
  if(ispred_dia == 0) dmvx_dia = dmvy_dia = 0.;

  e_x = float(fabs(dmvx_hor) + fabs(dmvx_ver));
  e_y = float(fabs(dmvy_hor) + fabs(dmvy_ver));

  if (e_x < 3) {
    *ctx_x = 0;
  } else if (e_x > 15) {
    *ctx_x = 1;
  } else {
    *ctx_x = 2;
  }

  if (e_y < 3) {
    *ctx_y = 0;
  } else if (e_y > 15) {
    *ctx_y = 1;
  } else {
    *ctx_y = 2;
  }

  //  mvStat_setCTX(dmvx_hor, dmvy_hor, 
  //                dmvx_ver, dmvy_ver, 
  //                dmvx_dia, dmvy_dia);
}
