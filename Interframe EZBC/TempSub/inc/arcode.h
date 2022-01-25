#ifndef ARCODE_H
#define ARCODE_H

#include "structN.h"

#define Top_value        65535  /* largest code value */
#define Max_frequency    16383  /* largest frequency count */
#define First_qtr        Top_value/4+1
#define Half             2*First_qtr
#define Third_qtr        3*First_qtr

EXTERN int low, high, value, **freq, **cmf, **int2sym, **sym2int;
EXTERN int bits_to_follow;

EXTERN int aff_low, aff_high, aff_value, **aff_freq, **aff_cmf, **aff_int2sym, **aff_sym2int;
EXTERN int aff_bits_to_follow;

void encode_symbol(int symbol, int context = 0);
void ardecode_init(int symbols, int contexts = 1, 
		   const int **freq_init = NULL);
void ardecode_end();
void arencode_init(int symbols, int contexts = 1, 
		   const int **freq_init = NULL);
void arencode_end(int context = 0);
int decode_symbol(int context = 0);

void encode_word(int word, int context = 0);
int decode_word(int context = 0);
void update_model(int word, int context = 0);
float get_expected_length(int word, int context = 0);

#endif // ARCODE_H
