#define USED   1
#define UNUSED 0
#define ANAL   0
#define SYN    1
#define PRUNE  1        /* overlapped block matching */
#define Global 0
#define Local  1
#define UN 0
#define CT 1            /* central */
#define HO 2            /* horizontal */
#define VE 3            /* vertical */
#define DI 4            /* diagonal */
#define Bi 5            /* Bi-direction or intra mode */
#define ER 6            /* error */
#define orign 0

long int denoise_mctf_anal_ezbc( int curr, int GOP_counter,videoinfo info,
                                 enum FLAG first_GOP, enum FLAG Level_change,
                                 int remaining_frs );
long int mctf_anal_ezbc( int curr, int GOP_counter, videoinfo info,
                         enum FLAG first_GOP, enum FLAG Level_change,
                         int remaining_frs, long int *sum_mv, int simul_enc, float *upframe1, float *upframe2, float *upsamp_x);

void save_enc_status(videoinfo *info);

void resume_enc_status(videoinfo *info);

void denoise_mctf_syn_ezbc( int curr, int GOP_counter, long int *total_bytes_past,
                            videoinfo info, enum FLAG first_GOP,
                            enum FLAG Level_change, int remaining_frs );
void mctf_syn_ezbc( int curr, int GOP_counter, long int *total_bytes_past,
                    videoinfo info, enum FLAG first_GOP,
                    enum FLAG Level_change, int remaining_frs, int simul_dec, int theo_dec );

long int intra_encode( int curr, int GOP_counter, videoinfo info );
void intra_decode( int curr, int GOP_counter, long int *total_bytes_past,
                   videoinfo info );

long int three_D_anal( int curr, int GOP_counter, videoinfo info );
void three_D_syn( int curr, int GOP_counter, long int *total_bytes_past,
                  videoinfo info );

long int mc_anal( int curr, int GOP_counter, videoinfo info, enum FLAG first_GOP,
              enum FLAG Level_change, int remaining_frs, int simul_enc, float *upframe1, float *upframe2, float *upsamp_x );
void mc_syn( int curr, long int *total_bytes_past, videoinfo info,
             enum FLAG first_GOP, enum FLAG Level_change, int remaining_frs, int simul_dec, int theo_dec );

void mc_analysis( float *L1, float *H1, float *H0, float *fr1, float *fr2, float *fr3, 
                  float *mvx1, float *mvy1, float *mvx2, float *mvy2, float *mvx3, float *mvy3,
                  float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
                  float *mvx3_int, float *mvy3_int, vector_ptr mv_ref1, vector_ptr mv_ref2,
                  vector_ptr mv_ref3, int hor, int ver, int level, int remaining_frs,
                  videoinfo info );

void   get_value_for_neighbor_mv(int cpos, int pos, ImageMEinfo *imagemeinfo, float *fr1, float *fr3, 
								 float *neighbor_value, float *self_weight, float neighbor_weight,
								 int cx, int cy, videoinfo info, int hor, int ver, 
								 int  UVcom, int t_level);


// by Yongjun Wu MCTF analysis with OBMC
void
mc_analysis_with_OBMC( float *L1, float *H1, float *H0, float *fr1, float *fr2, float *fr3,
             float *mvx1, float *mvy1, float *mvx2, float *mvy2, float *mvx3, float *mvy3,
             float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
             float *mvx3_int, float *mvy3_int, vector_ptr mv_ref1, vector_ptr mv_ref2,
             vector_ptr mv_ref3, int hor, int ver, int level, int remaining_frs,
             videoinfo info, ImageMEinfo *imagemeinfo, Varblkarrayinfo *varblkarray, int UVcom,
			 enum FLAG left_scene, enum FLAG right_scene, int type);


void mc_synthesis( float *fr0, float *fr1, float *H0, float *L1, float *H1, float *frp,
                   float *mvx0, float *mvy0, float *mvx1, float *mvy1, 
                   float *mvx2, float *mvy2, float *mvx0_int, float *mvy0_int, 
                   float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
                   vector_ptr mv_ref0, vector_ptr mv_ref1, vector_ptr mv_ref2,
                   int hor, int ver, int level, videoinfo info );

void
mc_synthesis_with_OBMC( float *fr0, float *fr1, float *H0, float *L1, float *H1, float *frp,
              float *mvx0, float *mvy0, float *mvx1, float *mvy1, 
              float *mvx2, float *mvy2, float *mvx0_int, float *mvy0_int, 
              float *mvx1_int, float *mvy1_int, float *mvx2_int, float *mvy2_int,
              vector_ptr mv_ref0, vector_ptr mv_ref1, vector_ptr mv_ref2,
              int hor, int ver, int level, videoinfo info,
			  ImageMEinfo *imagemeinfo, Varblkarrayinfo *varblkarray, int UVcom, int type );


void spatial_anal( float *full, int hor, int ver, float *ll, float *lh,
                   float *hl, float *hh, int filterType );
void spatial_syn( float *full, int hor, int ver, float *ll, float *lh,
                  float *hl, float *hh, int filterType );

void spatial_anal_frame( YUVimage_ptr fr, videoinfo info, YUVimage_ptr low,
                         YUVimage_ptr high0, YUVimage_ptr high1,
                         YUVimage_ptr high2, int filterType );
void spatial_syn_frame( YUVimage_ptr fr, videoinfo info, YUVimage_ptr low,
                        YUVimage_ptr high0, YUVimage_ptr high1,
                        YUVimage_ptr high2, int filterType );

void
synscheme3( int curr, videoinfo info, enum FLAG first_GOP,
            enum FLAG Level_change, int remaining_frs, int simul_dec, int theo_dec );
