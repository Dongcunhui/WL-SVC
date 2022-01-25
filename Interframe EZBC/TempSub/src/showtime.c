#include <stdio.h>
#include <stdlib.h>

void
print_time( double sc )
{
  int hr, mn;
  hr = ( int )( sc / 3600.0 );
  sc -= 3600.0 * hr;
  mn = ( int )( sc / 60.0 );
  sc -= 60.0 * mn;

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
