/**  
@Description  compiler dependent definitions

*/

#ifndef _TYPES_COMPILERS_H_
#define _TYPES_COMPILERS_H_

// SGI CC seems to have no own "marker"
#if (defined _OS_IRIX_ || defined _OS_IRIX64_) && !defined __GNUC__ && !defined __SGICC__
#define __SGICC__
#endif

// HP CC seems to have no own "marker"
#if defined _OS_HPUX_ && !defined __GNUC__ && !defined __HPCC__
#define __HPCC__
#endif

// Microsoft Visual C++ compiler
#if defined _MSC_VER && !defined __MSCC__
#define __MSCC__
#endif

// HP CC doesn't support "signed" and "volatile" declarations
#ifdef __HPCC__
#define signed
#define volatile
#endif

// only GCC and SGI CC (in a higher version!) seems to know "bool"
#if !defined __GNUC__ && !defined __SGICC__
// SGI CC in a newer version knows "_BOOL" as marker for bool type
#if !defined __SGICC__ || !defined _BOOL
#define HAS_NO_BOOL
#endif
#endif

#if defined HAS_NO_BOOL
#define bool vg_boolean_t
#define true vg_true_value
#define false vg_false_value

typedef int bool;
static const bool true = 1;
static const bool false = 0;

#undef HAS_NO_BOOL
#endif

#endif /* _TYPES_COMPILERS_H_ */
