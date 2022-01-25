void analysis( float *in, int sx, int sy, int lenx, int leny, int hor,
               int ver, int filterType );
void synthesis( float *in, int sx, int sy, int lenx, int leny, int hor,
                int ver, int filterType );
int filter_coeff( float **lpf, float **hpf, int FILTER_TYPE, int analsyn );
