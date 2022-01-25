void update_frame_motion_field(FRAME_MOTION_FIELD *frame_motion_field, int x, int y, 
							   int xblk, int yblk, videoinfo info, vector_ptr fmv, float dmvx, float dmvy,
							   float mvx, float mvy);

void  clear_frame_motion_field(FRAME_MOTION_FIELD *frame_motion_field, videoinfo info);

void  clear_frame_motion_field_simp(SIMP_FRAME_MOTION_FIELD *frame_motion_field, videoinfo info);

void  copy_frame_motion_field(FRAME_MOTION_FIELD *src_frame_motion_field, SIMP_FRAME_MOTION_FIELD *dest_frame_motion_field, videoinfo info);

void  copy_frame_motion_field2(SIMP_FRAME_MOTION_FIELD *src_frame_motion_field, SIMP_FRAME_MOTION_FIELD *dest_frame_motion_field, videoinfo info);

void layer_structure_mv_subsample(videoinfo info, int t_level, int tindex, 
								  vector_ptr fmv, int encoder_side);

void
replace_enhancement_mv( vector_ptr fmv1_array, vector_ptr fmv2_array, vector_ptr fmv1, vector_ptr fmv2,
                 float *pmvx, float *pmvy, int num_symbol, int subpel,
                 int x, int y, int xblk, int yblk, int hor, int ver,
                 videoinfo info, int decode_parallelmv, int t_level, int bidir_exist, 
				 int blk_thresh, int count );

void layer_structure_mv_trim(videoinfo info, int t_level, int tindex, vector_ptr fmv, 
							 int count);


void get_median_predictor_motion_field(float *pmvx, float *pmvy,
									   FRAME_MOTION_FIELD *frame_motion_field, SIMP_FRAME_MOTION_FIELD *prev_frame_motion_field,
									  int x_pos, int y_pos, int xblock, int yblock,
			                          videoinfo info, int t_level, int blk_thresh);

void ec_get_contexts_motion_field(int *ctx_x, int *ctx_y, FRAME_MOTION_FIELD *frame_motion_field,
								  int x, int y, videoinfo info, int t_level, int blk_thresh);


void test( videoinfo info, unsigned char *qp, int estimated_overhead, unsigned long int *gop_mv, int curr_last, int type );

