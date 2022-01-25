#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EXTERN extern
#include "arcode.h"
#include <assert.h>
#include "mv_ec.h"
#include "mv_statistics.h"

//#define SEND_EOF_SYMBOL

int num_contexts;
int num_symbol;

/*H----------------------------*/
double log2(double f)
{
  double x;

  x = log10(f);
  x /= log10(2.0);
  return (x);
}

/*H----------------------------*/

/****************************************************************************/
/*                          bit_plus_follow()                               */
/****************************************************************************/
void bit_plus_follow(int bit)
{
  output_bit(bit);      /* output the bit */
  mv_bits ++;
//  printf("ONE BIT!\n");
  while (bits_to_follow > 0)
  {
    output_bit(!bit);   /* output bits_to_follow opposite bits */
	mv_bits ++;
//	printf("ONE BIT!\n");
    bits_to_follow--;   /* set bits_to_follow to zero */
  }
}



/****************************************************************************/
/*                           encode_symbol()                                */
/****************************************************************************/
void encode_symbol(int symbol, int context)
{
  int range;    /* size of the current code region */

  //  printf("encode_symbol\n");

  if (symbol < 0 || symbol > num_symbol + 1)
  {
    printf("error in encode_symbol() %d\n", symbol);
    exit(1);
  }

  if (context < 0 || context >= num_contexts)
  {
    printf("illegal context number in encode_symbol() (%d)\n", context);
    exit(1);
  }

  range = high - low + 1;
  high = low + (range * cmf[context][symbol - 1]) / cmf[context][0] - 1;  /* narrow the code */
  low = low + (range * cmf[context][symbol]) / cmf[context][0];   /* region */

  for (;;)
  {     /* loop to output bits */
    if (high < Half)
    {
	  mv_bits = 0;
      bit_plus_follow(0);       /* output 0 if in low half */
	  mv_res_bits += mv_bits;
    }
    else if (low >= Half)
    {   /* output 1 if in high half */
	  mv_bits = 0;
      bit_plus_follow(1);
	  mv_res_bits += mv_bits;
      low -= Half;      /* subtract offset to top */
      high -= Half;
    }
    else if (low >= First_qtr && high < Third_qtr)
    {   /* output an opposite bit */
      bits_to_follow++; /* later if in middle half */
      low -= First_qtr; /* subtract offset to middle */
      high -= First_qtr;
    }
    else
      break;    /* otherwise exit loop */

    low = 2 * low;      /* scale up code range */
    high = 2 * high + 1;
  }
}



void encode_word(int word, int context)
{
  if (word < 0 || word > num_symbol - 1)
  {
    printf("illegal word in encode_word() (%d)\n", word);
    exit(1);
  }
  if (context < 0 || context > num_contexts - 1)
  {
    printf("illegal context number in encode_word() (%d)\n", context);
    exit(1);
  }    

  encode_symbol(int2sym[context][word], context);
}


/****************************************************************************/
/*                               arencode_end()                             */
/****************************************************************************/
void arencode_end(int context)
{
  int i, j;

#ifdef SEND_EOF_SYMBOL
  encode_symbol(num_symbol + 1, context);    /* EOF symbol */
#endif

  bits_to_follow++;     /* output two bits for underflow */
  (low < First_qtr) ? bit_plus_follow(0) : bit_plus_follow(1);

  /* release the memory */
  for (i = 0; i < num_contexts; i++)
  {
    mvStat_writeFreqNumber(i, cmf[i][0], num_symbol);
    for (j = 0; j < num_symbol; j++) {
      mvStat_writeFreq(freq[i][int2sym[i][j]]);
    }

    free(freq[i]);
    free(cmf[i]);
    free(int2sym[i]);
    free(sym2int[i]);
  }
  free(freq);
  free(cmf);
  free(int2sym);
  free(sym2int);
}

/****************************************************************************/
/*                            start_model()                                 */
/****************************************************************************/
void start_model(int context = 0, const int *freq_init = NULL)
{
  int i;

  if (context < 0 || context >= num_contexts)
  {
    printf("illegal context number in start_model() (%d)\n", context);
    exit(1);
  }

  for (i = 0; i < num_symbol; i++)
  {     /* setup the tables that */
    int2sym[context][i] = i + 1; /* translate between symbol */
    sym2int[context][i + 1] = i; /*  indexes and characters */
  }

  if (freq_init)
  {
#ifdef SEND_EOF_SYMBOL
    freq[context][num_symbol + 1] = 1;
#else
    freq[context][num_symbol + 1] = 0;
#endif
    cmf[context][num_symbol + 1] = 0;
    for (i = num_symbol; i >= 0; i--)
    {
      freq[context][i] = freq_init[sym2int[context][i]];
      cmf[context][i] = cmf[context][i + 1] + freq[context][i + 1];
    }

    freq[context][0] = 0;  
  }
  else
  {

#ifdef SEND_EOF_SYMBOL
    for (i = 0; i <= num_symbol + 1; i++)
    {
      freq[context][i] = 1;
      cmf[context][i] = num_symbol + 1 - i;
    }
#else
    freq[context][0] = 0; // EOF
    for (i = 0; i <= num_symbol; i++)
    {
      freq[context][i] = 1;
      cmf[context][i] = num_symbol - i;
    }
#endif // SEND_EOF_SYMBOL

    freq[context][0] = 0;
  }
}

/****************************************************************************/
/*                               arencode_init()                            */
/****************************************************************************/
void arencode_init(int symbols, int contexts, const int **freq_init)
{
  int i;

  if (symbols < 1)
  {
    printf("illegal number of symbols in arencode_init() (%d)\n", symbols);
    exit(1);
  }

  if (contexts < 1)
  {
    printf("illegal number of contexts in arencode_init() (%d)\n", contexts);
    exit(1);
  }

  num_symbol = symbols;  //num_symbol = 2
  num_contexts = contexts;  //num_contexts = 5

  /* memory allocation */
  freq = (int **) calloc(num_contexts, sizeof(int *));
  cmf = (int **) calloc(num_contexts, sizeof(int *));
  int2sym = (int **) calloc(num_contexts, sizeof(int *));
  sym2int = (int **) calloc(num_contexts, sizeof(int *));
  for (i = 0; i < num_contexts; i++)
  {
    freq[i] = (int *) calloc(num_symbol + 2, sizeof(int));
    cmf[i] = (int *) calloc(num_symbol + 2, sizeof(int));
    int2sym[i] = (int *) calloc(num_symbol, sizeof(int));
    sym2int[i] = (int *) calloc(num_symbol + 2, sizeof(int));
  }

  /* initialize variables */
  bits_to_follow = 0;

  /* start model */
  for (i = 0; i < num_contexts; i++)
  {
    start_model(i, freq_init ? freq_init[i] : NULL);
  }
  low = 0;
  high = Top_value;
}

/****************************************************************************/
/*                              ardecode_end()                              */
/****************************************************************************/
void ardecode_end()
{
  int i, j;

  for (i = 0; i < num_contexts; i++)
  {
    mvStat_writeFreqNumber(i, cmf[i][0], num_symbol);
    for (j = 0; j < num_symbol; j++) {
      mvStat_writeFreq(freq[i][int2sym[i][j]]);
    }

    free(freq[i]);
    free(cmf[i]);
    free(int2sym[i]);
    free(sym2int[i]);
  }
  free(freq);
  free(cmf);
  free(int2sym);
  free(sym2int);
}

/****************************************************************************/
/*                              ardecode_init()                             */
/****************************************************************************/
void ardecode_init(int symbols, int contexts, const int **freq_init)
{
  int i, bit;

  if (symbols < 1)
  {
    printf("illegal number of symbols in ardecode_init() (%d)\n", symbols);
    exit(1);
  }

  if (contexts < 1)
  {
    printf("number of contexts number in ardecode_init() (%d)\n", contexts);
    exit(1);
  }

  num_symbol = symbols;
  num_contexts = contexts;

  /* memory allocation */
  freq = (int **) calloc(num_contexts, sizeof(int *));
  cmf = (int **) calloc(num_contexts, sizeof(int *));
  int2sym = (int **) calloc(num_contexts, sizeof(int *));
  sym2int = (int **) calloc(num_contexts, sizeof(int *));
  for (i = 0; i < num_contexts; i++)
  {
    freq[i] = (int *) calloc(num_symbol + 2, sizeof(int));
    cmf[i] = (int *) calloc(num_symbol + 2, sizeof(int));
    int2sym[i] = (int *) calloc(num_symbol, sizeof(int));
    sym2int[i] = (int *) calloc(num_symbol + 2, sizeof(int));

    start_model(i, freq_init ? freq_init[i] : NULL);
  }

  bits_to_follow = 0;

  low = 0;
  high = Top_value;
  value = 0;    /* input bits ot fill the */

  for (i = 1; i <= Code_value_bits; i++)
  {     /* code value */
    input_bit(bit);
    value = 2 * value + bit;
  }
}

/****************************************************************************/
/*                             decode_symbol()                              */
/****************************************************************************/
int decode_symbol(int context)
{
  int range;    /* size of current code region */
  int cum;      /* cumulative frequency calculated */
  int symbol;   /* symbol decoded */
  int bit;

  if (context < 0 || context >= num_contexts)
  {
    printf("illegal context number in decode_symbol() (%d)\n", context);
    exit(1);
  }

  range = high - low + 1;
  cum = (((value - low) + 1) * cmf[context][0] - 1) / range;     /*find cum freq for value */

  for (symbol = 1; cmf[context][symbol] > cum; symbol++);        /* then find symbol */

  high = low + (range * cmf[context][symbol - 1]) / cmf[context][0] - 1;
  low = low + (range * cmf[context][symbol]) / cmf[context][0];
  for (;;)
  {     /* loop to get rid of bits */
    if (high < Half);   /*  do nothing */
    else if (low >= Half)
    {   /*  expand low half */
      value -= Half;
      low -= Half;      /* subtract offset to top */
      high -= Half;
    }
    else if (low >= First_qtr && high < Third_qtr)
    {   /* expand middle half */
      value -= First_qtr;
      low -= First_qtr; /* subtract offset to middle */
      high -= First_qtr;
    }
    else
      break;    /* otherwise exit loop */

    low = 2 * low;
    high = 2 * high + 1;        /* scale up code range */
    input_bit(bit);
    value = 2 * value + bit;    /* move in next input bit */
  }

  return symbol;
}

int decode_word(int context)
{
  int symbol = decode_symbol(context);

  if (symbol < 1 || symbol > num_symbol)
  {
    return -1;
  }
  else
  {
    return sym2int[context][symbol];
  }
}


/****************************************************************************/
/*                            update_model()                                */
/****************************************************************************/
void update_model(int word, int context)
{
  int i;        /* new index for symbol */
  int cum;
  int ch_i, ch_symbol;
  int symbol;

  if (word < 0 || word > num_symbol - 1)
  {
    printf("illegal word in update_model() (%d)\n", word);
    exit(1);
  }

  if (context < 0 || context >= num_contexts)
  {
    printf("illegal context number in update_model() (%d)\n", context);
    exit(1);
  }

  symbol = int2sym[context][word];

  if (cmf[context][0] >= Max_frequency)
  {     /* see if frequency counts */
    mvStat_writeFreqNumber(context, cmf[context][0], num_symbol);
    for (i = 0; i < num_symbol; i++) {
      mvStat_writeFreq(freq[context][int2sym[context][i]]);
    }
    cum = 0;
    for (i = num_symbol + 1; i >= 0; i--)
    {   /* if so, halve all the counts */
      freq[context][i] = (freq[context][i] + 1) / 2;      /*   (keeping them non-zero)   */
      cmf[context][i] = cum;
      cum += freq[context][i];
    }
  }

  for (i = symbol; freq[context][i] == freq[context][i - 1]; i--);        /* find symbol's new index */

  if (i < symbol)
  {
    ch_i = sym2int[context][i];
    ch_symbol = sym2int[context][symbol];        /* update the translation */
    sym2int[context][i] = ch_symbol;     /*  tables if the symbol has */
    sym2int[context][symbol] = ch_i;     /*  moved */
    int2sym[context][ch_i] = symbol;
    int2sym[context][ch_symbol] = i;
  }
  freq[context][i]++;    /* increment the frequency */
  while (i > 0)
  {
    i--;        /*  count for the symbol and */
    cmf[context][i]++;   /* update the cumulative freq. */
  }
}

float get_expected_length(int word, int context)
{
  if (word < 0 || word > num_symbol - 1)
  {
    printf("illegal word in get_expected_length() (%d)\n", word);
    exit(1);
  }

  if (context < 0 || context >= num_contexts)
  {
    printf("illegal context number in get_expected_length() (%d)\n", context);
    exit(1);
  }
  
  return -float(log2(double(freq[context][int2sym[context][word]]) / double(cmf[context][0])));
}


/****************************************************************************/
/*                               est_number()                               */
/****************************************************************************/
int est_number(int inum)
{
  int numbyte;

  if (inum < 0)
  {
    printf("error in est_number() %d\n", inum);
    exit(1);
  }

  if (inum < 128)
    numbyte = 1;
  /*else if(inum<32786) numbyte=2; */
  else if (inum <= (int)0xffff)
    numbyte = 2;
  else if (inum <= (int)0xffffffff)
    numbyte = 4;
  else
  {
    printf("error in est_number() too large\n");
    exit(1);
  }

  return numbyte;
}



/*****************************************************************************/
/*                                  entropy()                                */
/*****************************************************************************/
float entropy(int codeMAX, int total, int *pmf)
{
  int i, sum;
  float H, tprob;

  if (!total)
    return 0.;

  H = 0.;
  sum = 0;
  for (i = 0; i < codeMAX; i++)
  {
    tprob = (float) pmf[i] / total;
    if (tprob > 0. && tprob <= 1.)
      H += -tprob * (float) log2(tprob);
    else if (tprob < 0. || tprob > 1.)
    {
      printf("error in entropy() 1\n");
      exit(1);
    }
    sum += pmf[i];
  }
  if (sum != total)
    printf("error in entropy() sum = %d, total = %d\n", sum, total);

  return H;
}
