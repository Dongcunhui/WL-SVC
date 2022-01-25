#include "golomb.h"
#include "mv_ec.h"
#include <assert.h>

int theGrad0;
int theMaxLevels;

void golomb_init(int grad0, int maxLevels)
{
  theGrad0 = grad0;
  theMaxLevels = maxLevels;
}

int golomb_encode_word(int word, int write)
{
  int level, res, numbits;
  int i;

  assert(word >= 0);
         //  if (word < 0) {
         //    printf("illegal word in golomb_encode_word (%d)\n", word);
         //  }

  res = 1 << theGrad0;
  level = 1;
  numbits = 1 + theGrad0;

  while (word >= res && level < theMaxLevels) {
    word -= res;
    res <<= 1;
    level++;
    numbits += 2;
  }

  if (level >= theMaxLevels) {
    if (word >= res) {
      word = res - 1; // crop if too large
    }
    res = 0;
    numbits -= 1;
  }

  if (write) {
    // data bits
    res |= word;

    // write bits
    for (i = numbits - 1; i >= 0; i--) {
      output_bit((res >> i) & 1);
    }
  }

  return numbits;
}

int golomb_decode_word()
{
  int level;
  int bit, databits;
  int i;

  level = 0;
  while (level + 1 < theMaxLevels) {
    input_bit(bit);
    if (bit) break;
    level++;
  }

  databits = 0;
  for (i = 0; i < theGrad0 + level; i++) {
    input_bit(bit);
    databits = (databits << 1) | bit;
  }

  return (((1 << level) - 1 << theGrad0) + databits);
}
