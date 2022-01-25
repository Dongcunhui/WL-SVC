void snr_frame( float *ysnr, float *usnr, float *vsnr, YUVimage_ptr codeframe,
                YUVimage_ptr inframe, videoinfo info );
float calsnr( int start, int last, videoinfo info );
void calsnr1( YUVimage_ptr fr0, YUVimage_ptr fr1, int start, int last,
              videoinfo info );
/* void print_stat( videoinfo info ); */
void print_mvbits( videoinfo info, Rate FrsRate );
