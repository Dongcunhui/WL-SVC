
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//         G E N E R A L   P U R P O S E   F U N C T I O N S

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//           > > > >    C++ version  1.08 -  02/01/96   < < < <

// Amir Said - amir@densis.fee.unicamp.br
// University of Campinas (UNICAMP)
// Campinas, SP 13081, Brazil

// William A. Pearlman - pearlman@ecse.rpi.edu
// Rensselaer Polytechnic Institute
// Troy, NY 12180, USA

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// Copyright (c) 1995, 1996 Amir Said & William A. Pearlman

// This program is Copyright (c) by Amir Said & William A. Pearlman.
// It may not be redistributed without the consent of the copyright
// holders. In no circumstances may the copyright notice be removed.
// The program may not be sold for profit nor may they be incorporated
// in commercial programs without the written permission of the copyright
// holders. This program is provided as is, without any express or
// implied warranty, without even the warranty of fitness for a
// particular purpose.


// - - Inclusion - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "general.h"
#include <time.h>
#include <string.h>

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1e6
#endif


// - - Static variable - - - - - - - - - - - - - - - - - - - - - - - - -

//static Char_Line line;



// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//  Error-handling functions

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

void
Error( char *s )
{
  fprintf( stderr, "\n\n\a -> Error: " );
  fputs( s, stderr );
  fputs( "\n Execution terminated!\n", stderr );
  exit( 1 );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Warning( char *s )
{
  fprintf( stderr, "\n\n\a -> Warning: " );
  fputs( s, stderr );
  //  Pause(  );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Test_Pointer( void *p, char *orig )
{
  if( p == NULL ) {
    fputs( "\n\n\a -> Error: insufficient memory.", stderr );
    if( orig != NULL )
      fprintf( stderr, "\nOrigin = %s\n", orig );
    fputs( "\n Execution terminated!\n", stderr );
    exit( 1 );
  }
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//  Functions of the class < Chronometer >

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

void
Chronometer::start( char *s )
{
  if( s != NULL )
    puts( s );
  if( stat )
    Warning( "chronometer already on!" );
  else {
    mark = clock(  );
    stat = 1;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Chronometer::stop( void )
{
  if( stat ) {
    elp += clock(  ) - mark;
    stat = 0;
  } else
    Warning( "chronometer already off!" );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

float
Chronometer::read( void )
{
  return float ( stat ? elp + ( clock(  ) - mark ) : elp ) / CLOCKS_PER_SEC;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
Chronometer::display( char *s )
{
  float sc =
    float ( stat ? elp + ( clock(  ) - mark ) : elp ) / CLOCKS_PER_SEC;
  int hr = int ( sc / 3600.0 );
  sc -= ( float )3600.0 *hr;
  int mn = int ( sc / 60.0 );
  sc -= ( float )60.0 *mn;
  if( s != NULL )
    printf( " %s ", s );
  if( hr ) {
    printf( "%d hour", hr );
    if( hr > 1 )
      printf( "s, " );
    else
      printf( ", " );
  }
  if( ( hr ) || ( mn ) ) {
    printf( "%d minute", mn );
    if( mn > 1 )
      printf( "s, and " );
    else
      printf( ", and " );
  }
  printf( "%5.2f seconds.\n", sc );
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// end of file < General.C >
