//EXTERN rdsum rd[3];
/*EXTERN mvnode_ptr  mvtop; APR23*/
//EXTERN int yfclass_flag;

/* New */
extern long int totalY[5], totalU[5], totalV[5], totalmap[5], totalMV[5];
/* New */

long int ezbc3d_enc( int curr, YUVimage * pyrFrs, videoinfo info,
                     long int GOP_counter, float &gop_mse );
void release_3dezbc_enc(  );
void print_time( double sc );

void init_enc( videoinfo * info );
