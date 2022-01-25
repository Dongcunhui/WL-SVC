#include <assert.h>
#include "mv_ec.h"
#define EXTERN extern
#include "arcode.h"

int **contextMap = NULL;
int numSeparateBins = -1;

int generateContextMap(int contexts, const int *multiContexts)
{
  int i, j;
  int contextCounter;

  assert(contextMap == NULL);

  contextCounter = 0;
  contextMap = (int**)malloc((numSeparateBins + 2) * sizeof(int*));
  for (i = 0; i < numSeparateBins + 2; i++)
  {
    contextMap[i] = (int*)malloc(contexts * sizeof(int));

    if (multiContexts[i])
    {
      for (j = 0; j < contexts; j++)
      {
        contextMap[i][j] = contextCounter;
        contextCounter++;
      }
    } else {
      for (j = 0; j < contexts; j++)
      {
        contextMap[i][j] = contextCounter;
      }
      contextCounter++;
    }
  }

  return contextCounter;
}

void freeContextMap()
{
  int i;

  assert(contextMap != NULL);

  for (i = 0; i < numSeparateBins + 2; i++)
  {
    free(contextMap[i]);
  }
  free(contextMap);

  contextMap = NULL;
}

/****************************************************************************/
/*                               arencode_init()                            */
/****************************************************************************/
void binar_enc_init(int contexts, const int numSepBins, 
                    const int *multiContexts, const int **freq_init)
{
  int totalContexts;

  numSeparateBins = numSepBins;

  totalContexts = generateContextMap(contexts, multiContexts);

  arencode_init(2, totalContexts, freq_init);
}

/****************************************************************************/
/*                               arencode_end()                             */
/****************************************************************************/
void binar_enc_end(int context)
{
  assert(contextMap != NULL);

  arencode_end();

  freeContextMap();
}

/****************************************************************************/
/*                              ardecode_init()                             */
/****************************************************************************/
void binar_dec_init(int contexts, const int numSepBins, 
                    const int *multiContexts, const int **freq_init)
{
  int totalContexts;

  numSeparateBins = numSepBins;

  totalContexts = generateContextMap(contexts, multiContexts);

  ardecode_init(2, totalContexts, freq_init);
}

/****************************************************************************/
/*                              ardecode_end()                              */
/****************************************************************************/
void binar_dec_end()
{
  assert(contextMap != NULL);

  ardecode_end();

  freeContextMap();
}

void binar_encode_word(int word, int context)
{
  int i;
  int abs_word, bin;
  
  assert(contextMap != NULL);

  abs_word = (word < 0) ? (-word) : word;
  
  bin = 0;
  for (i = 0; i < abs_word; i++) {
    encode_word(0, contextMap[bin][context]);
    if (bin < numSeparateBins) bin++;
  }
  encode_word(1, contextMap[bin][context]);
  
  if (word != 0) {
    encode_word(word < 0, contextMap[numSeparateBins + 1][context]);
  }
}

int binar_decode_word(int context)
{
  int bin, abs_word;

  assert(contextMap != NULL);

  bin = 0;
  abs_word = 0;
  for (;;) {
    if (decode_word(contextMap[bin][context])) break;
    abs_word++;
    if (bin < numSeparateBins) bin++;
  }

  if (abs_word == 0)
  {
    return 0;
  } else {
    return decode_word(contextMap[numSeparateBins + 1][context]) ? 
      (-abs_word) : abs_word;
  }
}

/****************************************************************************/
/*                            update_model()                                */
/****************************************************************************/
void binar_update_model(int word, int context)
{
  int i;
  int abs_word, bin;
  
  assert(contextMap != NULL);

  abs_word = (word < 0) ? (-word) : word;
  
  bin = 0;
  for (i = 0; i < abs_word; i++) {
    update_model(0, contextMap[bin][context]);
    if (bin < numSeparateBins) bin++;
  }
  update_model(1, contextMap[bin][context]);
  
  if (word != 0) {
    update_model(word < 0, contextMap[numSeparateBins + 1][context]);
  }
}

float binar_get_expected_length(int word, int context)
{
  float length;
  int i;
  int abs_word, bin;
  
  assert(contextMap != NULL);

  length = 0.;

  abs_word = (word < 0) ? (-word) : word;
  
  bin = 0;
  for (i = 0; i < abs_word; i++) {
    length += get_expected_length(0, contextMap[bin][context]);
    if (bin < numSeparateBins) bin++;
  }
  length += get_expected_length(1, contextMap[bin][context]);
  
  if (word != 0) {
    length += get_expected_length(word < 0, 
                                  contextMap[numSeparateBins + 1][context]);
  }

  return length;
}


