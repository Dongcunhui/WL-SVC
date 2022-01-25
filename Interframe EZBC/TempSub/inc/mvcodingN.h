#ifndef MVCODINGN_H
#define MVCODINGN_H

#include "structN.h"
#include <stdio.h>

extern VLCtable blockmodeVLC[2];

int read_GOPheader( enum FLAG **dec_scene_change, videoinfo info );
int write_GOPheader( enum FLAG **scene_change, videoinfo info );
int mv_encoding( videoinfo info, Rate FrsRate, vector_ptr * yfmv, 
                 int GOP_counter, int simul_enc, int curr );
long int mv_decoding( videoinfo info, vector_ptr * yfmv, int GOP_counter, long int starting_pos, int simul_dec, int theo_dec );
int get_mvBytes( FILE * fp_mv );
int write_GOPheader( enum FLAG **scene_change, videoinfo info );

void decode_MV( long int *total_bytes_past, videoinfo info, int GOP_counter, long int starting_pos, int simul_dec, int theo_dec );

int get_mode_coding_cost(enum BiMode bi_mode, vector_ptr fmv2, int bi_sign, int t_level);

void clean_mrg_mv(vector_ptr fmv);

#endif
