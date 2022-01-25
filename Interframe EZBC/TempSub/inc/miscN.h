#ifndef _MISCN_H_
#define _MISCN_H_

void info2header( videoheader_ptr header, videoinfo info );     /* miscN.c */
void header2info( videoinfo_ptr info, videoheader header );
void read_header( char *name, videoinfo * info );
void write_header( char *name, videoinfo info );
long int get_GOP_num( videoinfo info );

int get_mvBytes( FILE * fp_mv );

float get_mean( float *frame, int lenx, int leny, int hor );

int subband_location( int *sx, int *sy, int *lenx, int *leny, int curr,
                      int band, int hor, int ver );
float variance( float *frame, int lenx, int leny, int hor );

double variance_mean( float *frame, int lenx, int leny, float *mean );
void varianceGOP( YUVimage_ptr frames, videoinfo info, Rate_ptr alloc );
void wcopyframe( YUVimage * source, YUVimage * dest, float weight,
                 videoinfo info );
void copyframe( YUVimage * source, YUVimage * dest, videoinfo info );
void computegain( YUVimage lowband, YUVimage highband, videoinfo info );

void mv_copy( vector_ptr source, vector_ptr dest, videoinfo info );

#endif // _MISCN_H_
