#include "stdio.h"
#include "stdlib.h"
#include <assert.h>


void
write_number_core( int outputbyte, FILE * fp )
{
  int tmpbyte;
  tmpbyte = ( outputbyte & 0xff000000 ) >> 24;
  putc( tmpbyte, fp );
  tmpbyte = ( outputbyte & 0x00ff0000 ) >> 16;
  putc( tmpbyte, fp );
  tmpbyte = ( outputbyte & 0x0000ff00 ) >> 8;
  putc( tmpbyte, fp );
  tmpbyte = outputbyte & 0x000000ff;
  putc( tmpbyte, fp );
  return;
}

int
read_number_core( FILE * fp )
{
  int i, tmpbyte, inputbyte;

  tmpbyte = getc( fp );
  inputbyte = tmpbyte & 0x7f;

//  printf("inputbyte = %d\n",inputbyte);

  for( i = 0; i < 3; i++ ) {
    inputbyte <<= 8;
    inputbyte += getc( fp );
//	printf("inputbyte = %d\n",inputbyte);
  }

  return ( inputbyte );
}


unsigned int
read_number_sub( FILE * fp )
{
  int tmpbyte; 
  unsigned int inputbyte;

  tmpbyte = getc( fp );
  inputbyte = tmpbyte & 0xff;

  inputbyte <<= 8;
  inputbyte += getc( fp );

  return ( inputbyte );
}


void
write_number_sub(unsigned int outputbyte, FILE * fp )
{
  int tmpbyte;

  tmpbyte = ( outputbyte & 0xff00 ) >> 8;
  putc( tmpbyte, fp );
  
  tmpbyte = outputbyte & 0x00ff;
  
  putc( tmpbyte, fp );
  return;
}

int
write_substream_length( int length, FILE * fp )
{
  int i, tmpbyte;

  if( length <= 0x7fff ) {
    for( i = 1; i >= 0; i-- ) {
      tmpbyte = length >> 8 * i;
      putc( tmpbyte, fp );
    }

    return 2;
  } else if( length <= 0x7fffffff ) {
    length = length | 0x80000000;
    for( i = 3; i >= 0; i-- ) {
      tmpbyte = length >> 8 * i;
      putc( tmpbyte, fp );
    }
    return 4;
  } else {
    printf( "substream %d is too long (dataio.c)\n", length );
    exit( 0 );
  }
}

int
read_substream_length( int *length, FILE * fp )
{
  int i, tmpbyte;

  tmpbyte = getc( fp );

  if( tmpbyte & 0x80 ) {
    for( i = 2; i >= 0; i-- ) {
      tmpbyte <<= 8;
      tmpbyte += getc( fp );
    }
    *length = tmpbyte & 0x7fffffff;;
    return 4;
  } else {
    tmpbyte <<= 8;
    tmpbyte += getc( fp );
    *length = tmpbyte;
    return 2;
  }
}
