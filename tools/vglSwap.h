/**  
@Description functions for value exchange and byte swapping

*/

#ifndef _TOOLS_SWAP_H_
#define _TOOLS_SWAP_H_

#include <std.h>

namespace Biopool {


/* exchange values */

template <class T> inline void rotate(T& a, T& b, T&c)
{
  T tmp = a;
  a = b;
  b = c;
  c = tmp;
}

template <class T> inline void swap(T& a, T& b)
{
  T tmp = a;
  a = b;
  b = tmp;
}

/* little/big endian conversion */

inline  void swap(vg_uint16& a)
{
  a = (a << 8) | (a >> 8);
}

inline  void swap(vg_int16& a)
{
  a = (vg_int16)((vg_uint16)a << 8) | ((vg_uint16)a >> 8);
}

inline  void swap(vg_uint32& a)
{
  vg_uint16* p = (vg_uint16*)(&a);
  *p = (*p << 8) | (*p >> 8);
  p++;
  *p = (*p << 8) | (*p >> 8);
  a = (a >> 16) | (a << 16);
}

inline  void swap(vg_int32& a)
{
  vg_uint16* p = (vg_uint16*)(&a);
  *p = (*p << 8) | (*p >> 8);
  p++;
  *p = (*p << 8) | (*p >> 8);
  a = (vg_int32)(((vg_uint32)a >> 16) | ((vg_uint32)a << 16));
}

template <class T> inline  void swap(T& a)
{
  switch (sizeof(T)) {
  case 2:
    swap(*(vg_uint16*)&a);
    break;
  case 4:
    swap(*(vg_uint32*)&a);
    break;
  default:
    break;
  }
}

} // end-namespace

#endif /* _TOOLS_SWAP_H_ */
