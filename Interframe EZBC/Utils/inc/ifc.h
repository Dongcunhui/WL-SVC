#ifndef _IFC_H
#define _IFC_H

#include <limits.h>

/*****************************************************************************/
/*                                  std_int                                  */
/*****************************************************************************/

#if (INT_MAX == 2147483647)
typedef int std_int;
#elif (LONG_MAX == 2147483647)
typedef long int std_int
#else
# error "Platform does not support 32-bit integers"
#endif

#define MAX_STD_INT ((std_int) 0x7FFFFFFF)
#define MIN_STD_INT ((std_int)-0x80000000)

/*****************************************************************************/
/*                                 std_short                                 */
/*****************************************************************************/

#if (SHRT_MAX == 32767)
typedef short int std_short;
#else
# error "Platform does not support 16-bit integers"
#endif

/*****************************************************************************/
/*                                std_ushort                                 */
/*****************************************************************************/

#if (SHRT_MAX == 32767)
typedef unsigned short int std_ushort;
#else
# error "Platform does not support 16-bit integers"
#endif

/*****************************************************************************/
/*                                  std_byte                                 */
/*****************************************************************************/

typedef unsigned char std_byte;


/* ========================================================================= */
/* -------------------- IMPLEMENTATION PRECISION --------------------------- */
/* ========================================================================= */

#ifdef IMPLEMENT_32
#  define IMPLEMENTATION_PRECISION 32
#endif /* IMPLEMENT_32 */

#ifdef IMPLEMENT_16
#  define IMPLEMENTATION_PRECISION 16
#endif /* IMPLEMENT_16 */

#ifndef IMPLEMENTATION_PRECISION
#  define IMPLEMENTATION_PRECISION 32   /* Must be one of 16 or 32 currently. */
#endif /* IMPLEMENTATION_PRECISION */

#if (IMPLEMENTATION_PRECISION == 16)
typedef short int ifc_int;
#elif (IMPLEMENTATION_PRECISION == 32)
#  if (INT_MAX == 2147483647)
typedef int ifc_int;
#  elif (LONG_MAX == 2147483647)
typedef long int ifc_int;
#  else
#    error "Platform does not support 32 bit integers!"
#  endif
#else
#  error "IMPLEMENTATION_PRECISION must be one of 16 or 32!"
#endif

#define MIN_IFC_INT ((ifc_int)((-1)<<(IMPLEMENTATION_PRECISION-1)))
#define MAX_IFC_INT (~MIN_IFC_INT)

#endif
