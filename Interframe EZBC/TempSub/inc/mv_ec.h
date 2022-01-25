#ifndef MV_EC_H
#define MV_EC_H

#include <stdio.h>
#include <stdlib.h>
#include "structN.h"

#define Code_value_bits  16

enum coderType {AR_NARY = 0, GOLOMB, AR_BINARY};

extern int buffer, bits_to_go, garbage_bits, outbyte, inbyte;
extern FILE *fptmp;
extern FILE *fpbit;

void ec_enc_init(videoinfo info, coderType coder, 
		 int symbols, int contexts = 1, int write = 1);
void ec_enc_end1();
void ec_enc_end2();
long int ec_dec_preinit(FILE *fp);
void ec_dec_init(coderType coder, int symbols, int contexts = 1);
void ec_dec_end1();
long int ec_dec_end2();
void ec_encode_word(int word, int context = 0);
int ec_decode_word(int context = 0);
void ec_update_model(int word, int context = 0);
float ec_get_expected_length(int word, int context = 0);

void ec_freeze();
void ec_unfreeze(int keep_current_states);

void ec_get_contexts(int *ctx_x, int *ctx_y, vector_ptr fmv, int x, int y,
                     videoinfo info, int t_level, int blk_thresh);

void done_splitted_bytes(unsigned char   store_splitted_bytes[65535], 
						 unsigned int    splitted_byte_num, videoinfo info,
						 int t_level);

void get_splitted_bytes(unsigned char    *store_splitted_bytes, 
						   unsigned int  *splitted_byte_num, videoinfo info,
						   int t_level, FILE *fpbit, long int starting_pos);


/****************************************************************************/
/****************************************************************************/

void encode_init(videoinfo info);
void encode_end(int with, videoinfo info);
void decode_init();
void decode_end();

/****************************************************************************/
/*                             output_bit()                                 */
/****************************************************************************/
#define output_bit(bit)                                                      \
 /**  int bit;  **/                                                          \
{                                                                            \
  buffer >>=1;                                                               \
  if(bit) buffer |= 0x80;     /* put bit in top of buffer */                 \
    bits_to_go--;                                                            \
    if(!bits_to_go){            /* output buffer if it is full */            \
      putc(buffer, fptmp);                                                   \
      bits_to_go=8;                                                          \
      outbyte++;                                                             \
    }                                                                        \
}


/****************************************************************************/
/*                             input_bit()     ¶ÁÒ»¸öbit                             */
/****************************************************************************/
#define input_bit(bit)                                                       \
/*** int bit ***/                                                            \
{                                                                            \
  if(!bits_to_go){   /* Read the next byte if no bits are left in buffer */  \
    if(inbyte>0) {   /* prevent reading for continous one bit files */       \
      buffer=getc(fpbit);                                                    \
      inbyte--;                                                              \
    }                                                                        \
    else{                                                                    \
      buffer = EOF;                                                          \
      garbage_bits++;                        /* return arbitrary bits */     \
      if(garbage_bits>Code_value_bits-2){    /*  after EOF, but check */     \
		printf("Bad input file\n");    /*  for too many such */              \
      }                                                                      \
    }                                                                        \
    bits_to_go = 8;                                                          \
  }                                                                          \
  bit = buffer&1;     /* return the next bit from */                         \
  buffer>>=1;         /*  the bottom of the byte */                          \
  bits_to_go--;                                                              \
}

#endif // MV_EC_H
