#ifndef MEMORY_H
#define MEMORY_H

void enc_destructor( videoinfo info );
void dec_destructor( videoinfo info );

void alloc_vector( vector_ptr fmv, videoinfo info );
void free_vector( vector_ptr fmv, videoinfo info );
void free_mvs( vector_ptr * yfmv, videoinfo info );

void frame_alloc( YUVimage_ptr frame, videoinfo info );
void free_frame( YUVimage frame );

void alloc_transmap( transmap ** map, int hor, int ver );
void free_transmap( transmap * map );

#endif // MEMORY_H
