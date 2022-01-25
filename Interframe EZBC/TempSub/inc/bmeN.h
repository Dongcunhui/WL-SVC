#ifndef BME_H
#define BME_H

#define NOISE_VAR 3.0

#define  FIR88     88
#define   BILINEAR  22
#ifdef use_BILINEAR_type
#define   TYPE       BILINEAR   //    BILINEAR
#else
#define   TYPE       FIR88   //    BILINEAR
#endif

#define   AFF_TYPE  BILINEAR


void block_matching( vector_ptr fmv1, vector_ptr fmv2, vector_ptr fmv3, vector_ptr fmv4, 
		     YUVimage *fr_cur, YUVimage *fr_ref1, YUVimage *fr_ref2,
                     enum FLAG *sc1, enum FLAG *sc2,
                     videoinfo info, int level, int dist, int subpel, float *upframe1, float *upframe2, float *upsamp_x );
void position( int *px, int *py, float pfx, float pfy, float mvx, float mvy,
               int hor, int ver );
int inbound( float x, float y, int hor, int ver );
void blockmv2pixelmv( vector_ptr fmv, float **ymvx, float **ymvy,
                      float **cmvx, float **cmvy, enum LiftingMode block_mode,
                      videoinfo info, int t_level );

void get_cvector( float *cmvx, float *cmvy, float *ymvx, float *ymvy,
                  int yhor, int yver, int chor, int cver, videoinfo info, int t_level );

void full_search_fast( float *mvx1, float *mvy1, 
                       float *frame_cur, float *frame_ref1, float *frame_ref2, 
                       float *upframe1, float *upframe2,
                       float mvx2, float mvy2, 
                       int cx, int cy, int xblock, int yblock,
                       int maxx, int maxy, int hor, int ver, int half, 
                       float lambda, float pmvx, float pmvy, 
                       int ctx_x, int ctx_y,
                       float *sad_cost, float *bit_cost, float *total_cost,
                       int do_parallel);  // add parallel mode. mwi 

/* int overlap_size( int x ); */
float MCP_Error( float *frame1, int m, float *frame0, int n, int oxblk,
                 int oyblk, int hor, float *first_pred, int k, int xblk_max );
float MCP_Error2(float *frame1, int m, float *frame0, int n, int oxblk,
                 int oyblk, int hor, float *first_pred, float min_cost,
                  float *frame2, int n2); // add parallel mode. mwi 
float Subpel_MCP_Error( float *frame1, int m, float *frame0, float px,
                        float hx, float py, float hy, int xblk, int yblk,
                        int hor, int ver, 
                        float *first_pred, int k, int xblk_max,
                        float *frame2, float ppx, float ppy); // add parallel mode. mwi 
float Subpel_MCP_Error2( float *frame1, int m, float *upframe, int px, int hx,
                         int py, int hy, int xblk, int yblk, int hor, int ver,
                         int uphor, int upver, int step, 
			 float *first_pred, int k, int xblk_max  );
float Subpel_MCP_Error3( float *frame1, int m, float *upframe, int px, int py, 
                         int xblk, int yblk, int hor, int ver,
                         int uphor, int upver, int step, float *first_pred,
                         float min_cost,
                         float *upframe2, int ppx, int ppy);
float Subpel_MCP_Error4( float *frame_cur,  
                         int cx, int cy, int xblock, int yblock, 
                         float *frame_ref1, float *upframe_ref1, 
                         float mvx1, float mvy1,
                         float *frame_ref2, float *upframe_ref2, 
                         float mvx2, float mvy2,
                         int hor, int ver, int subpel );
void generate_child( vector_ptr fmv, float *mvx, float *mvy, float *mad );

void find_MSE( vector_ptr fmv1, vector_ptr fmv2,
               float *frame_cur, float *frame_ref1, float *frame_ref2,
               int cx, int cy, int xblock, int yblock, int hor, int ver, 
               int t_level );

void initialize_spiral_search(int search_range);
void clean_up_spiral_search(int search_range);

float getFrameSATD(float *fr, int hor, int ver);

void find_diff( vector_ptr fmv, float *frame1, float *frame2, float *D,
          int cx, int cy, int xblock, int yblock, int hor, int ver, 
          int t_level, int type);

float local_diff_search( float *mvx, float *mvy, float *D1, float *D2, int *Dx,
          int xblk, int yblk, int hor, int ver, int t_level );

float get_sad( float mvx, float mvy, float *D1, float *D2, int *Dx,
          int xblk, int yblk, int t_level, int type );

int two_comp_est( float *mvx, float *mvy, vector_ptr fmv, float *frame_cur, float *frame_ref1, float *frame_ref2, float *SAD,
				int cx, int cy, int xblk, int yblk, int hor, int ver, int t_level );

#endif // BME_H 
