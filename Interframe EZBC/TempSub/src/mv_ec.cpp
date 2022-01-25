#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "mv_ec.h"
#include "golomb.h"
#include "binarcode.h"
#define EXTERN
#include "arcode.h"
#include "bme_tools.h"
#include "mv_statistics.h"

int buffer, bits_to_go, garbage_bits, outbyte, inbyte;
FILE *fptmp;
extern FILE *fpbit;

int numSymbols = -1;
coderType theCoder = coderType(-1);
videoinfo theInfo;
int theWrite = -1;
long int mvBytes;

/**************************************************************************
 *     Initialization for context adaptive binary arithmetic coder        *
 **************************************************************************/

const int numSeparateBins = 3; // number of bins with own context models
                               // additional models for
                               //   a) remaining bins
                               //   b) sign bit
const int multiContexts[numSeparateBins + 2] = {1, 0, 0, 0, 0};

const int initFreqDef[2] = {1, 1};
const int *initFreq[12] = {initFreqDef, initFreqDef, initFreqDef, 
                          initFreqDef, initFreqDef, initFreqDef,
                          initFreqDef, initFreqDef, initFreqDef,
                          initFreqDef, initFreqDef, initFreqDef};

/****************************************************************************/

void write_number();
long int read_number(int with);

// save 2 bytes for the splitted motion vector bits
void write_splitted_number(int outputbyte, FILE *fp)
{
  int tmpbyte;
  tmpbyte = (outputbyte & 0xff00) >>8;    putc(tmpbyte, fp);
  tmpbyte =  outputbyte & 0x00ff;         putc(tmpbyte, fp);
  return;
}





void done_splitted_bytes(unsigned char   store_splitted_bytes[65535], 
						   unsigned int  splitted_byte_num, videoinfo info,
						   int t_level)
{
	int i;
	FILE *fpbit_split;

	// open the bitfile and copy the number of data + bit data 
	if (!(fpbit_split = fopen(info.bitname, "a+b")))
	{
		printf("append_bit: %s\n", info.bitname);
		exit(1);
	}
  
	if(splitted_byte_num<0 || splitted_byte_num>=0xffff){
		printf("splitted motion vector bits are not valid---%d!\n",  splitted_byte_num); 
		getchar();
		exit(1);
	}
	
	write_splitted_number(splitted_byte_num, fpbit_split);

	for(i=0 ; i<(int)splitted_byte_num ; i++) 
		putc(store_splitted_bytes[i], fpbit_split);

	fclose(fpbit_split); 
}

// 读取两个字节
unsigned int read_splitted_number(FILE *fp)
{
  int tmpbyte;
  unsigned int inputbyte;

  tmpbyte = getc( fp );// 获得一个字符
  inputbyte = tmpbyte & 0xff;

  inputbyte <<= 8;
  inputbyte += getc( fp );

  return inputbyte;
}


void get_splitted_bytes(unsigned char    *store_splitted_bytes, 
						   unsigned int  *splitted_byte_num, videoinfo info,
						   int t_level, FILE *fpbit, long int starting_pos)
{
	int i;

	fseek(fpbit, starting_pos, SEEK_SET); 
    *splitted_byte_num = read_splitted_number(fpbit);// 数量
	for(i=0 ; i<(int)(*splitted_byte_num) ; i++) // 
		store_splitted_bytes[i] = getc(fpbit);
}


/****************************************************************************/
/*                               append_bit()                               */
/****************************************************************************/
void append_bit(char *bitname, char *tmpname, int with)
{
  int i;

  /* open the bitfile and copy the number of data + bit data */
  if (!(fpbit = fopen(bitname, "a+b")))
  {
    printf("append_bit: %s\n", bitname);
    exit(1);
  }
  if (!(fptmp = fopen(tmpname, "rb")))
  {
    printf("append_bit: %s\n", tmpname);
    exit(1);
  }

  /* write the number of byte and copy data */
  if (with)
    write_number();     

  for (i = 0; i < outbyte; i++)
    putc(getc(fptmp), fpbit);

  fclose(fpbit);
  fclose(fptmp);

  if (with)
  {
    outbyte += 4;
  }
}

/****************************************************************************/
/*                               arencode_init()                            */
/****************************************************************************/
void ec_enc_init(videoinfo info, coderType coder, int symbols, int contexts,
		 int write)
{
  assert(theWrite < 0 && theCoder < 0 && numSymbols < 0);
  assert(write >= 0);

  numSymbols = symbols;
  theCoder = coder;
  theInfo  = info;
  theWrite = write;

  switch(theCoder) {
    case AR_NARY:
      arencode_init(symbols, contexts, NULL);
      break;
    case AR_BINARY:
      binar_enc_init(contexts, numSeparateBins, multiContexts, initFreq);
      break;
    case GOLOMB:
      golomb_init();
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }

  /* open the bitfile */
  if (theWrite > 0)
  {
    if (!(fptmp = fopen(info.tmpname, "wb")))
    {
      printf("ec_enc_init: %s\n", info.tmpname);
      exit(1);
    }
  }

  /* initialize variables */
  buffer = 0;
  bits_to_go = 8;
  outbyte = 0;
}

/****************************************************************************/
/*                               arencode_end()                             */
/****************************************************************************/
void ec_enc_end1()
{
  assert(theWrite >= 0 && theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      arencode_end();
      break;
    case AR_BINARY:
      binar_enc_end();
      break;
    case GOLOMB:
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
}

void ec_enc_end2()
{
  if (theWrite > 0)
  {
    if (bits_to_go < 8)
    {
      putc(buffer >> bits_to_go, fptmp);
      outbyte++;
//	  mv_res_bits += (8 - bits_to_go);
    }
    fclose(fptmp);

    /* write to the bit file */
    append_bit(theInfo.bitname, theInfo.tmpname, 1);
//	mv_res_bits += 32;
  }

  theCoder = coderType(-1);
  theWrite = -1;
  numSymbols = -1;
}

/****************************************************************************/
/*                              ardecode_init()                             */
/****************************************************************************/
long int ec_dec_preinit(FILE *fp)
{
  assert(fp != NULL);

  fpbit = fp;

  mvBytes = read_number(1);

  /* initialize variables */
  // buffer = 0;
  bits_to_go = 0;
  garbage_bits = 0;

  return mvBytes; 
}


void ec_dec_init(coderType coder, int symbols, int contexts)
{
  assert(theWrite < 0 && theCoder < 0 && numSymbols < 0);

  theCoder = coder;
  numSymbols = symbols;

  switch(theCoder) {
    case AR_NARY:
      ardecode_init(symbols, contexts, NULL);
      break;
    case AR_BINARY:
      binar_dec_init(contexts, numSeparateBins, multiContexts, initFreq);
      break;
    case GOLOMB:
      golomb_init();
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
}

/****************************************************************************/
/*                              ardecode_end()                              */
/****************************************************************************/
void ec_dec_end1(){
  assert(theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      ardecode_end();
      break;
    case AR_BINARY:
      binar_dec_end();
      break;
    case GOLOMB:
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
}
long int ec_dec_end2()
{
  /* read if there is something left */
  while (inbyte > 0)
  {
    getc(fpbit);
    inbyte--;
  }

  theCoder = coderType(-1);
  theWrite = -1;
  numSymbols = -1;

  return mvBytes;
}

void ec_encode_word(int word, int context)
{
  assert(theWrite >= 0 && theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      encode_word(word + numSymbols / 2, context);
      break;
    case AR_BINARY:
      binar_encode_word(word, context);
      break;
    case GOLOMB:
      mv_res_bits += golomb_encode_word((word > 0) ? (word * 2 - 1) : (-word * 2), 1);
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
}

int ec_decode_word(int context)
{
  int word;

  assert(theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      return decode_word(context) - numSymbols / 2;
      break;
    case AR_BINARY:
      return binar_decode_word(context);
      break;
    case GOLOMB:
      word = golomb_decode_word();
      return (word % 2) ? (word / 2 + 1) : (-word / 2);
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
  
  return -1;
}

/****************************************************************************/
/*                            update_model()                                */
/****************************************************************************/
void ec_update_model(int word, int context)
{
  assert(theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      update_model(word + numSymbols / 2, context);
      break;
    case AR_BINARY:
      binar_update_model(word, context);
      break;
    case GOLOMB:
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
}

float ec_get_expected_length(int word, int context)
{
  assert(theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      return get_expected_length(word + numSymbols / 2, context);
      break;
    case AR_BINARY:
      return binar_get_expected_length(word, context);
      break;
    case GOLOMB:
      return (float)golomb_encode_word((word > 0) ? (word * 2 - 1) : (-word * 2), 0);
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }

  return -1.;
}

void ec_freeze()
{
  assert(theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      printf("freezing is not yet implemented for AR_NARY\n");
      exit(1);
      break;
    case AR_BINARY:
      printf("freezing is not yet implemented for AR_BINARY\n");
      exit(1);
      break;
    case GOLOMB:
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
}
  
void ec_unfreeze(int keep_current_states)
{
  assert(theCoder >= 0 && numSymbols >= 0);

  switch(theCoder) {
    case AR_NARY:
      printf("unfreezing is not yet implemented for AR_NARY\n");
      exit(1);
      break;
    case AR_BINARY:
      printf("unfreezing is not yet implemented for AR_BINARY\n");
      exit(1);
      break;
    case GOLOMB:
      break;
    default:
      printf("illegal codec in ec_enc_init.\n");
      exit(1);
      break;
  }
}

void ec_get_contexts(int *ctx_x, int *ctx_y, vector_ptr fmv, int x, int y,
                     videoinfo info, int t_level, int blk_thresh)
{
  float dmvx_hor, dmvy_hor, dmvx_ver, dmvy_ver, dmvx_dia, dmvy_dia;
  enum FLAG ispred_hor, ispred_ver, ispred_dia;
  float e_x, e_y;

  get_dmv(&dmvx_hor, &dmvy_hor, &ispred_hor, fmv, x - 1, y, info, t_level, blk_thresh);
  get_dmv(&dmvx_ver, &dmvy_ver, &ispred_ver, fmv, x, y - 1, info, t_level, blk_thresh);
  get_dmv(&dmvx_dia, &dmvy_dia, &ispred_dia, fmv, x - 1, y - 1, info, t_level, blk_thresh);

  if(ispred_hor == NO) dmvx_hor = dmvy_hor = 0.;
  if(ispred_ver == NO) dmvx_ver = dmvy_ver = 0.;
  if(ispred_dia == NO) dmvx_dia = dmvy_dia = 0.;

  e_x = float(fabs(dmvx_hor) + fabs(dmvx_ver));
  e_y = float(fabs(dmvy_hor) + fabs(dmvy_ver));

  if (e_x < 3) {
    *ctx_x = 0;
  } else if (e_x > 15) {
    *ctx_x = 1;
  } else {
    *ctx_x = 2;
  }

  if (e_y < 3) {
    *ctx_y = 0;
  } else if (e_y > 15) {
    *ctx_y = 1;
  } else {
    *ctx_y = 2;
  }

  //  mvStat_setCTX(dmvx_hor, dmvy_hor, 
  //                dmvx_ver, dmvy_ver, 
  //                dmvx_dia, dmvy_dia);
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/

void encode_init(videoinfo info)
{
  /* open the bitfile */

  if (!(fptmp = fopen(info.tmpname, "wb")))
  {
    printf("can not open: %s\n", info.tmpname);
    exit(1);
  }

  /* initialize variables */
  buffer = 0;
  bits_to_go = 8;
  outbyte = 0;

}

void encode_end(int with, videoinfo info)
{
  if (bits_to_go < 8)
  {
    putc(buffer >> bits_to_go, fptmp);
    outbyte++;
  }
  fclose(fptmp);

  /* write to the bit file */
  append_bit(info.bitname, info.tmpname, with);
}

void decode_init()
{
//  int i, bit;

  garbage_bits = 0;
  bits_to_go = 0;
}

void decode_end()
{
  /* read if there is something left */
  while (inbyte > 0)
  {
    getc(fpbit);
    inbyte--;
  }
}
