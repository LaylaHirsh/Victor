/**  
@Description  system independent data types

*/

#ifndef _TYPES_STD_H_
#define _TYPES_STD_H_

#define OLD_IOSTREAMS		// all systems seem to have "oldstyle" IO streams

#include <stddef.h>		// for "ptrdiff_t" and others


typedef signed char    vg_int8;	// signed integer with exactly 8 bits
typedef unsigned char  vg_uint8; // unsigned integer with exactly 8 bits
typedef signed short   vg_int16; // signed integer with exactly 16 bits
typedef unsigned short vg_uint16; // unsigned integer with exactly 16 bits
typedef signed int     vg_int32; // signed integer with exactly 32 bits
typedef unsigned int   vg_uint32; // unsigned integer with exactly 32 bits
typedef float          vg_ieee32; // IEEE floating point with exactly 32 bits
typedef double         vg_ieee64; // IEEE floating point with exactly 64 bits

extern const char* const vglVersion;

// simple structure to provide some information about classes.
// (for internal use only)
template <class T> 
class vgClassInfo
{
public:
  static const char* name;	// name of the class
  static const size_t size;	// size in bytes
  static const bool numeric;	// true if +-*/ defined
  static const bool ordered;	// true if comparison defined
  static const T min;		// representable minimum
  static const T max;		// representable maximum
};

#endif /* _TYPES_STD_H_ */
