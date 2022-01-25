#ifndef BINARCODE_H
#define BINARCODE_H

#include "structN.h"

void binar_enc_init(int contexts, const int numSepBins, 
                    const int *multiContexts, const int **freq_init = NULL);
void binar_enc_end(int context = 0);
void binar_dec_init(int contexts, const int numSepBins, 
                    const int *multiContexts, const int **freq_init = NULL);
void binar_dec_end();

void binar_encode_word(int word, int context);
int binar_decode_word(int context);
void binar_update_model(int word, int context);
float binar_get_expected_length(int word, int context);

#endif // BINARCODE_H
