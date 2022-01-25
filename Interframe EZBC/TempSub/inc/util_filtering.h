#ifndef UTIL_FILTERING_H
#define UTIL_FILTERING_H

#define ANAL 0
#define SYN 1

extern const float static_prequant_weight_high[];
extern const float static_prequant_weight_low[];

extern const float prequant_weight_high[];
extern const float prequant_weight_low[];
extern const float copycomp_weight_high[];
extern const float copycomp_weight_low[];


int next_roundtrip_index( int curr_i, float *array, int length, int *direction );
void extend_line( float *orig_line, int orig_len, int margin,
             int algorithm, int flag, float *new_line );

void line_convolve( float *input, float *extension, int length, float *fco,
                    int flength, int analsyn_flag, float *out );
void interpolate_filter( );
void temporal_filter( );

float FIRinterpolate( float fx, float fy, float *frame, int length, int hor,
                      int ver );

float interpolate( float fx, float fy, float *frame, int hor, int ver,
                   int type );

void Interpolate_frame( float *frame, int hor, int ver, float **upframe,
                        int *uph, int *upv, int subpel );

void Interpolate_frame2( float *frame, int hor, int ver, float **upframe, float *upsamp_x,
                         int subpel );

#endif // UTIL_FILTERING_H
