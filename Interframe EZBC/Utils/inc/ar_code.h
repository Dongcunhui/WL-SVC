
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//            A R I T H M E T I C   C O D E   C L A S S E S

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//         > > > >    C++ version 3.02  -  02/01/96    < < < <

// Amir Said - amir@densis.fee.unicamp.br
// University of Campinas (UNICAMP)
// 13081 Campinas, SP, Brazil

// William A. Pearlman - pearlman@ecse.rpi.edu
// Rensselaer Polytechnic Institute
// Troy, NY 12180, USA


// C++ implementation of the arithmetic-coding algorithm by I.H. Witten,
// R.M. Neal, and J.G. Cleary, published in ``Arithmetic coding for data
// compression,'' {\em Commun. ACM}, vol.~30, pp.~520--540, June 1987.


// - - Inclusion - - - - - - - - - - - - - -------- - - - - - - - - - - - - - -
#ifndef AR_CODE
#define AR_CODE
#include <iostream>
#ifndef General_H
#include "general.h"
#endif

#include "structN.h"


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//  Class definitions

// = = = = = = = = = = = = = = = = = = = = = = = = = = = == = = = = = = = = = =

const double ZERO_ENT = -1.0;

double ent( double p0 );
double xent( double p0, double pp0 );

class Histo
{
protected:
  static const double OneLog2;
  static const int MaxFreq;
};

class HistoBi:public Histo
{
protected:
  int c0;
  int c1;
public:
    HistoBi( void )
  {
    c0 = c1 = 0;
  }
  HistoBi( HistoBi & hist )
  {
    reset( hist );
  }
  void reset( int i )
  {
    c0 = c1 = i;
  }
  void reset( HistoBi & hist )
  {
    c0 = hist.c0;
    c1 = hist.c1;
  }
  int get_c0( void )
  {
    return c0;
  }
  double get_p0( void )
  {
    if( c0 && c1 )
      return double ( c0 ) / double ( c0 + c1 );
    else
    return ZERO_ENT;
  }
  double get_p1( void )
  {
    if( c0 && c1 )
      return double ( c1 ) / double ( c0 + c1 );
    else
    return ZERO_ENT;
  }
  int get_c1( void )
  {
    return c1;
  }
  int get_total( void )
  {
    return c0 + c1;
  }
  double entropy0( void )
  {
    if( c0 && c1 )
      return -log( double ( c0 ) / double ( c0 + c1 ) )*OneLog2;
    else
      return ZERO_ENT;
  }
  double entropy1( void )
  {
    if( c0 && c1 )
      return -log( double ( c1 ) / double ( c0 + c1 ) )*OneLog2;
    else
      return ZERO_ENT;
  }
  double bit_cost( int bit )
  {
    if( bit )
      return entropy1(  );
    else
      return entropy0(  );
  }
  double cost0( void )
  {
    if( c0 && c1 )
      return c0 * entropy0(  );
    else
      return ZERO_ENT;
  }
  double cost1( void )
  {
    if( c0 && c1 )
      return c1 * entropy1(  );
    else
      return ZERO_ENT;
  }
  double total_cost( void )
  {
    if( c0 && c1 )
      return cost1(  ) + cost0(  );
    else
      return ZERO_ENT;
  }
  double entropy_avg( void )
  {
    if( c0 && c1 )
      return ( cost0(  ) + cost1(  ) ) / ( double )( c0 + c1 );
    else
      return ZERO_ENT;
  }
  void update( int c )
  {
    if( c )
      c1++;
    else
      c0++;
  }
  void reset( void )
  {
    c0 = c1 = 1;
  }
  void add_to_c0( int n )
  {
    c0 += n;
  }
  void add_to_c1( int n )
  {
    c1 += n;
  }
  int get_min( void )
  {
    return ( c0 < c1 ) ? c0 : c1;
  }
};


class HistoBiModel:public HistoBi
{
protected:
  friend class Encoder;
  friend class Decoder;
  int max_f;
  static float ThreshRecScale;  //max allowed XH after recursives_cale()
  static float ThreshUpdate;    //max allowed XH after update()
  static int MinMaxFreq;
  static int InitMaxFreq;
  static int UpdateScheme;
  static int MAX_MIN;
  void update( int c );
  void update_test_total( void );
  void update_test_total_and_min( void );
public:
    HistoBiModel( int max_freq = InitMaxFreq ):HistoBi(  ), max_f( max_freq )
  {
  }
  HistoBiModel( HistoBi & hist ):HistoBi( hist )
  {
  }
  void set_symbols( int s )
  {
    max_f = InitMaxFreq;
    if( s == 2 )
      reset( 1 );
    else {
      fprintf( stderr, "set_symbols: model is binary.\n" );
      exit( 1 );
    }
  }
  void set_symbols( int s, int max_freq )
  {
    set_symbols( s );
    max_f = max_freq;
  }
  void scale( void )
  {
    c0 = ( c0 + 1 ) >> 1;
    c1 = ( c1 + 1 ) >> 1;
  }
  void recursive_scale( void );
  void taub_scale( int min_thre = 2, int total_thre = 32, int trunc_flag =
                   0 );
  void set_max_freq( int max_freq )
  {
    while( get_total(  ) >= max_freq )
      scale(  );
    max_f = max_freq;
  }
  void reset_maxf( int max_freq )
  {
    max_f = max_freq;
  }
  int get_maxf( void )
  {
    return max_f;
  }
  void set_model( HistoBiModel & model )
  {
    reset( model );
    max_f = model.max_f;
  }
  void set_update_scheme( int scheme )
  {
    UpdateScheme = scheme;
  }
};


class Adaptive_Model
{
  // . private data .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  int max_f, nsb, *cum_f, *freq, *indx_to_sb, *sb_to_indx;

  // . friend classes  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  friend class Encoder;

  friend class Decoder;

  // . private functions  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  void update( int index );

  int select_symbol( long value, long *l, long *h );

  void new_interval( int symb, long *l, long *h );

  // . constructors and destructor .  .  .  .  .  .  .  .  .  .  .  .  .

    public:Adaptive_Model( void )
  {
    nsb = 0;
    cum_f = NULL;
  }

  Adaptive_Model( int ns )
  {
    nsb = 0;
    cum_f = NULL;
    set_symbols( ns );
  }

  ~Adaptive_Model( void )
  {
    nsb = 0;
    delete cum_f;
  }

  // . public functions   .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  void reset( void );

  void set_symbols( int ns );

  void set_symbols( int s, int max_freq );

  float get_cost( int symbol );

  void scale( void );

  void set_max_freq( int max_freq );

  void set_model( Adaptive_Model & model );
};                              // end definition of class < Adaptive_Model >

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
class Encoder
{
  // . private data .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
  FILE *out_file;

  int bit_buffer, bit_index, closed, temp;

  long byte_counter, symbol_counter, low, high, bits_to_follow;

  // . private functions  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  void bit_plus_follow( int b );

  void update_interval( void );

  void reset( char *file_name );

  void reset_append( char *file_name );

  void stop( void );

  // . constructor  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

    public:Encoder( void )
  {
    byte_upto_bitplane = NULL;
    closed = 1;
  }                             //062602
   ~Encoder( void )
  {                             //printf("run destructor of Encoder (ar_code.h)\n"); 
    if( byte_upto_bitplane )
      delete[]byte_upto_bitplane;
  }

  // . public data   .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  long int *byte_upto_bitplane; // allocated in 062602

  char temp_name[250];

  static int object_count;      //062602


  // . public functions   .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  long bytes_used( void )
  {                             //062602
    if( bit_index == 8 )
      return byte_counter;
    else
      return byte_counter + 1;

  }

  long symbols_encoded( void )
  {
    return symbol_counter;
  }

  //int getCount(void) {return object_count;}

  void new_open_file( void );   //062602

  void open_file( void );

  void new_open_file( char *file_name );

  void open_file( char *file_name );

  void close_file( char *file_name );

  void close_file( void );

  void code_symbol( int s, Adaptive_Model & );

  void code_bits( int bits, int word );

  int code_bit( int bit );
  //   void code_bit_arith(int bit, HistoBiModel &);
  void code_symbol( int bit, HistoBiModel & );
};                              // end definition of class < Encoder >

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
class Decoder
{
  // . private data .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .


  int bit_buffer, bit_index, extra_bits, ext_stop, closed;

  long low, high, value, symbol_counter, byte_counter, mss_symbols;

  // . private functions  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  void input_byte( void );

  void update_interval( void );

  // . constructor  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

    public:int byte_budget;     // not include header
  FILE *in_file;


    Decoder( void )
  {
    closed = 1;
  }

  // . public functions   .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

  int end_of_file( void )
  {
    return ( symbol_counter >= mss_symbols );
  }

  long bytes_used( void )
  {
    return byte_counter;
  }

  long symbols_decoded( void )
  {
    return symbol_counter;
  }

  long message_symbols( void )
  {
    return mss_symbols;
  }

  void new_open_file( char *file_name, long shift_bytes );
  void new_open_file( char *file_name, long shift_bytes, long budget ); //062602

  void open_file( char *file_name );

  void close_file( void );

  int decode_symbol( Adaptive_Model & );

  int decode_bits( int bits );

  int decode_bit( void );
  //    int decode_bit_arith(HistoBiModel &hist);
  int decode_symbol( HistoBiModel & hist );
};                              // end definition of class < Decoder >

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

//  Inline functions

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


inline void
HistoBiModel::update_test_total( void )
{
  if( get_total(  ) == max_f )
    scale(  );
}

inline void
HistoBiModel::update_test_total_and_min( void )
{
  if( get_min(  ) == MAX_MIN || get_total(  ) == max_f )
    scale(  );
}

//#if 0

inline void
HistoBiModel::update( int c )
{
  //update_test_total_and_min(); //suggested by Taubman
  update_test_total(  );
  if( c == 0 )
    c0++;
  else
    c1++;
}

//#endif

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

inline int
Encoder::code_bit( int bit )
{
  long lm1 = low - 1, range = high - lm1;
  if( bit )
    low += range >> 1;
  else
    high = lm1 + ( range >> 1 );

  update_interval(  );
  return bit;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

inline void
Encoder::code_symbol( int bit, HistoBiModel & hist )
{
  long lm1 = low - 1;
  long range = ( high - lm1 ) * hist.get_c0(  ) / hist.get_total(  );

  if( bit )
    low += range;
  else
    high = lm1 + range;

  hist.update( bit );
  update_interval(  );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline int
Decoder::decode_bit( void )
{

  long lm1 = low - 1, range = high - lm1;
  int bit = int ( ( ( ( value - lm1 ) << 1 ) - 1 ) / range );

  if( bit ) {
    low += range >> 1;
  } else
    high = lm1 + ( range >> 1 );

  update_interval(  );
  return bit;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline int
Decoder::decode_symbol( HistoBiModel & hist )
{

  long lm1 = low - 1, range = high - lm1;

  int cum = int ( ( ( value - lm1 ) * hist.get_total(  ) - 1 ) / range );

  int bit = ( cum >= hist.get_c0(  ) );
  range = range * hist.get_c0(  ) / hist.get_total(  );

  if( bit )
    low += range;
  else
    high = lm1 + range;

  hist.update( bit );
  update_interval(  );
  return bit;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#endif
// end of file < Ar_Code.H >
