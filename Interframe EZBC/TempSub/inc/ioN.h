void read_frame( YUVimage_ptr frame, videoinfo info, char *inname, int index,
                 enum FORMAT format );
void write_frame( YUVimage frame, videoinfo info, char *inname, int index,
                  enum FORMAT format );

void print_mv( vector_ptr fmv, int cx, int cy, int xblk, int yblk, int hor,
               int ver );
void write_GOP( int curr, YUVimage ** pyrTemp, YUVimage * pyrFrs,
                videoinfo info, enum FLAG first_GOP, int remaining_frs );

void write_block_mode_motion_vector(char *direction, int GOP_counter, int count, int  frame_type, 
									vector_ptr fmv,  videoinfo info, int t_level);

void write_frame_into_file(YUVimage frame, videoinfo info,  int index, int GOP_counter);

void read_jp2k_frame( YUVimage_ptr frame, videoinfo info, char *inname, int index,
                 enum FORMAT format );
void write_jp2k_frame( YUVimage frame, videoinfo info, char *inname, int index,
                  enum FORMAT format );

