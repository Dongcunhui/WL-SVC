void f420_444( videoinfo info, YUVimage * frame );
void f444_420( videoinfo info, YUVimage * frame );

void f444_422( videoinfo info, YUVimage * frame );
void f422_444( videoinfo info, YUVimage * frame );
void YUV2RGB( float Y, float U, float V, float *R, float *G, float *B );
void RGB2YUV( float R, float G, float B, float *Y, float *U, float *V );
