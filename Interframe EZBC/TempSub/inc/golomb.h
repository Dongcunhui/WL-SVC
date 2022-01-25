#ifndef GOLOMB_H
#define GOLOMB_H

#define MAX_LEVELS 100
#define GRAD0 0

void golomb_init(int grad0 = GRAD0, int maxLevels = MAX_LEVELS);
int golomb_encode_word(int word, int write);
int golomb_decode_word();

#endif // GOLOMB_H
