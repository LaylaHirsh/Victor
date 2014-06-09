/**  
@Description  system independent byte order definitions

*/

#ifndef _TYPES_ENDIAN_H_
#define _TYPES_ENDIAN_H_

#if defined _OS_Linux_
#include <endian.h>
#endif
#if defined _OS_IRIX_ || defined _OS_IRIX64_
#include <sys/endian.h>
#endif
#if defined _OS_OSF1_
#include <machine/endian.h>
#endif
#if defined _OS_AIX_
#include <sys/machine.h>
#endif
#if defined __APPLE__
#include <sys/param.h>
//#else
//#include <i386/endian.h>
#endif

#if !defined BIG_ENDIAN
#define BIG_ENDIAN 4321
#endif

#if !defined LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#endif

#if !defined BYTE_ORDER
#if defined __i386 || defined WIN32 || defined __x86_64
#define BYTE_ORDER LITTLE_ENDIAN
#endif
#if defined __sparc__ || defined __mips__
#define BYTE_ORDER BIG_ENDIAN
#endif
#endif

#if !defined BYTE_ORDER || !defined LITTLE_ENDIAN || !defined BIG_ENDIAN
#error FATAL compile error: byte order not defined
#endif

#endif /* _TYPES_ENDIAN_H_ */
