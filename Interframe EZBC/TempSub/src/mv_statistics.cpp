#include "mv_statistics.h"
#include "basic.h"
#include "general.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

const int lazy = 1;

FILE *mvStat_fp = NULL; 
int mvStat_GOP = -1;
int mvStat_level = -1;
int mvStat_frame = -1;
int mvStat_xPos = -1;
int mvStat_yPos = -1;
int mvStat_numFreq = -1;

float mvStat_dmvx, mvStat_dmvy;
float mvStat_dmvx_hor, mvStat_dmvy_hor;
float mvStat_dmvx_ver, mvStat_dmvy_ver;
float mvStat_dmvx_dia, mvStat_dmvy_dia;

void mvStat_open(char *filename)
{
  if (lazy) return;

  assert(mvStat_fp == NULL);

  mvStat_fp = fopen(filename, "w");

  assert(mvStat_fp != NULL);

  fprintf(mvStat_fp, "dmv: GOP level frame xpos ypos dmvx dmvy "
          "dmvx_hor dmvy_hor dmvx_ver dmvy_ver dmvx_dia dmvy_dia\n");
  fprintf(mvStat_fp, "frq: GOP level frame ctx elements freq0 ...\n");
}

void mvStat_setFrame(int GOP, int level, int frame)
{
  if (mvStat_fp == NULL) return;

  assert(GOP >= 0 && level >= 0 && frame >= 0);

  mvStat_GOP = GOP;
  mvStat_level = level;
  mvStat_frame = frame;
}

void mvStat_setPos(int x, int y)
{
  if (mvStat_fp == NULL) return;

  assert(x >= 0 && y >= 0);

  mvStat_xPos = x;
  mvStat_yPos = y;
}

void mvStat_setDMV(float dmvx, float dmvy)
{
  mvStat_dmvx = dmvx;
  mvStat_dmvy = dmvy;
}

void mvStat_setCTX(float dmvx_hor, float dmvy_hor,
                   float dmvx_ver, float dmvy_ver,
                   float dmvx_dia, float dmvy_dia)
{
  mvStat_dmvx_hor = dmvx_hor;
  mvStat_dmvy_hor = dmvy_hor;
  mvStat_dmvx_ver = dmvx_ver;
  mvStat_dmvy_ver = dmvy_ver;
  mvStat_dmvx_dia = dmvx_dia;
  mvStat_dmvy_dia = dmvy_dia;
}

void mvStat_writeDMVCTX()
{
  if (mvStat_fp == NULL) return;
  
  fprintf(mvStat_fp, "dmv: %02d %02d %02d %03d %03d "
          "%.3f %.3f  %.3f %.3f  %.3f %.3f  %.3f %.3f\n",
          mvStat_GOP, mvStat_level, mvStat_frame, 
          mvStat_xPos, mvStat_yPos, 
          mvStat_dmvx, mvStat_dmvy, mvStat_dmvx_hor, mvStat_dmvy_hor, 
          mvStat_dmvx_ver, mvStat_dmvy_ver, mvStat_dmvx_dia, mvStat_dmvy_dia);
}

void mvStat_writeFreqNumber(int ctx, int cum, int number)
{
  if (mvStat_fp == NULL) return;

  fprintf(mvStat_fp, "frq: %02d %02d %02d %02d %4d ",
          mvStat_GOP, mvStat_level, mvStat_frame, ctx, cum);

  mvStat_numFreq = number;
}

void mvStat_writeFreq(int freq)
{
  if (mvStat_fp == NULL) return;

  assert(mvStat_numFreq > 0);

  fprintf(mvStat_fp, "%4d ", freq);

  mvStat_numFreq--;

  if (mvStat_numFreq == 0) {
    fprintf(mvStat_fp, "\n");
  }
}

void mvStat_close()
{
  if (mvStat_fp == NULL) return;

  fclose(mvStat_fp);

  mvStat_fp = NULL;
}


#ifdef BLOCKMODE_STATISTICS

// block size histogram
// 0: 1x1,   1: 2x2,   2: 4x4,     3: 8x8,     4: 16x16, 
// 5: 32x32, 6: 64x64, 7: 128x128, 8: 256x256
#define NUMBER_OF_BLOCKSIZES 9
long int block_count[NUMBER_OF_BLOCKSIZES] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

// block mode histogram
long int mode_count[NUMBER_OF_BI_MODES] = {0, 0, 0, 0, 0, 0, 0, 0};

// total area of block modes
double mode_area[NUMBER_OF_BI_MODES] = {0., 0., 0., 0., 0., 0., 0., 0.};


void record_blockmode_statistics(int xblk, int yblk, enum BiMode mode,
                                 int bi, int hor, int ver)
{
  block_count[int(rint(log2(double(MY_MAX(xblk, yblk)))))]++;
  if (bi) {
    mode_count[mode]++;
    mode_area[mode] += double(xblk * yblk) / double(hor * ver);
  }
}

void dump_blockmode_statistics(char *filename)
{
  FILE *outfile;
  char tmpname[1024];
  int i;

  sprintf(tmpname, "%s-blocksizes", filename);
  outfile = fopen(tmpname, "w");
  assert(outfile);
  for (i = 0; i < NUMBER_OF_BLOCKSIZES; i++) {
    fprintf(outfile, "%ld\n", block_count[i]);
  }
  fclose(outfile);

  sprintf(tmpname, "%s-blockmodes", filename);
  outfile = fopen(tmpname, "w");
  assert(outfile);
  for (i = 0; i < NUMBER_OF_BI_MODES; i++) {
    fprintf(outfile, "%ld %lf\n", mode_count[i], mode_area[i]);
  }
  fclose(outfile);
}
  
#endif
