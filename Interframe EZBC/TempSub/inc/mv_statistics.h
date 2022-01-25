#ifndef MV_STATISTICS_H
#define MV_STATISTICS_H

#include "structN.h"

void mvStat_open(char *filename);
void mvStat_setFrame(int GOP, int level, int frame);
void mvStat_setPos(int x, int y);
void mvStat_setDMV(float dmvx, float dmvy);
void mvStat_setCTX(float dmvx_hor, float dmvy_hor,
                   float dmvx_ver, float dmvy_ver,
                   float dmvx_dia, float dmvy_dia);
void mvStat_writeDMVCTX();
void mvStat_writeFreqNumber(int ctx, int cum, int number);
void mvStat_writeFreq(int freq);
void mvStat_close();

#ifdef BLOCKMODE_STATISTICS
void record_blockmode_statistics(int xblk, int yblk, enum BiMode,
                                 int bi, int hor, int ver);
void dump_blockmode_statistics(char *filename);
#endif

#endif
