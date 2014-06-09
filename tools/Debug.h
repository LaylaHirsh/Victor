// -*- C++ -*-
// $Id: Debug.h,v 1.6 2008-05-09 13:58:49 biocomp Exp $
// debug.h
// General error-handling and debugging facilities in C++
// 
// (c) 1998 by Matthias Heiler
//
// This header file provides some macros or inlined functions:
//
//      ASSERT(assertion, exception)
//          checks if assertion is true and throws exeption otherwise.
//          Additionally, a message containing LOCATION is displayed.
//
//      PRECOND(assertion, exception)
//      POSTCOND(assertion, exception)
//      INVARIANT(assertion, exception)
//          alias name for ASSERT. Used to produce more expressive
//          code and to make it easier to switch off pre- or
//          postconditions or invariants selectively.
//
//      ERROR(message, exception)
//          prints message and throws exception.
//
//      ERROR_IF(condition, message, exception)
//          throws exception if condition is true.
//
//      DEBUG_MSG(message)
//          prints a debugging message.
//
//      DEBUG_BREAK(condition, message)
//          if condition holds, prints message and waits for a keystroke. 
//
//      DUMP(x)
//          prints the variable x.
//
//      DEBUG_TRACE
//          creates object tracing and counting instances of a class.
//
//      PRINT_NAME
//          prints name of current function.
//
//      PRINT_LOCATION
//          prints location of current function.
//
// Some conditionals influence the behaviour:
//
//      NEXCEPTION
//          if defined, all exception throws are ignored and messages 
//          to cerr are wrote instead. ASSERT and ERROR will terminate
//          the program.
//
//      NDEBUG
//          if defined, assertions are switched off and DEBUG_MSGs 
//          and DUMPs are removed.
//
//      TERMINATION
//          defines how to terminate the program. Defaults to "exit(1)".
//          Maybe you prefer "abort()"  to abort with a core dump.
//      
//      VERBOSE=1
//      VERBOSE=2
//      VERBOSE=3
//          defines levels of verbosity.
//          VERBOSE = 1 enables DEBUG_MSG and DEBUG_BREAK
//          VERBOSE = 2 enables above and DUMP 
//          VERBOSE = 3 enables every possible output
//
//      OSTREAM
//          defines where to write the debug messages to (standard: "cerr").
//         
//  Note: you can avoid any runtime penalty by defining NEXCEPTIONS and
//        NDEBUG together.
//
//  Note: the header file supports internationalization with the GNU
//        "gettext()"-utility. Please to include <libintl.h> *before*
//        debug.h if you want to use this feature.
//


#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <exception>
using namespace std;

extern "C" {
#include <stdlib.h>
}

#include <iostream>
#include <fstream>
using namespace std;

// The user may specify other output streams than cerr.
#ifndef OSTREAM
#define OSTREAM cerr
#endif

#if VERBOSE >= 3
#define V3OUT OSTREAM
#else
#define V3OUT ofstream("/dev/null")
#endif

#if VERBOSE >= 2
#define V2OUT OSTREAM
#else
#define V2OUT ofstream("/dev/null")
#endif

#if VERBOSE >= 1
#define V1OUT OSTREAM
#else
#define V1OUT ofstream("/dev/null")
#endif

// If internalization is used, please include <libintl.h> *before* <debug.h>
#ifndef gettext
#define gettext(aString) aString
#endif

// The user may redefine TERMINATION to his needs (e.g. "abort()").
#ifndef TERMINATION
#define TERMINATION exit(1)
#endif

#ifndef NEXCEPTIONS
#define TERMINATION throw(exception)
#endif

#define LOCATION SourceLocator(__FILE__, __LINE__, __FUNCTION__)

/** SourceLocator: class to provide information about
    C++-Source-Location.

    SourceLocator is used for error reports to "OSTREAM". */
class SourceLocator {
public:
  SourceLocator(const char* f, long l, const char* fu=0) 
    : file(f), func(fu), line(l) { }
  friend ostream& operator<<(ostream& os, const SourceLocator& loc) 
  { if (loc.func != 0) { return os << loc.file << ":" << loc.line << ": " 
				   << loc.func << ": "; } 
  else { return os << loc.file << ":" << loc.line << ": "; } };
private:  
  const char* file;
  const char* func;
  long line;
};

/** Class helping to trace number of instances of a class. */
class ObjectTrace {
public:
  ObjectTrace() : ct(++count)
    { V3OUT << "DEBUG_TRACE: Object " << ct << " constructed." << endl; }
  // commented to remove useless compiler warning - ST 23.2.99:
  // ObjectTrace(const char* n) :  ct(++count)
  //  { V3OUT << "DEBUG_TRACE: Object " << ct << " constructed." << endl; }
  //  ObjectTrace(const ObjectTrace& orig) : ct(++count)
  //  { V3OUT << "DEBUG_TRACE: Object " << ct << " constructed." << endl; }
  ~ObjectTrace() 
    { V3OUT << "DEBUG_TRACE: Object " << ct << " destroyed." << endl; }
private:
  static unsigned count;
  unsigned        ct;
};

#ifndef NDEBUG

#define PRINT_LOCATION V3OUT << "PRINT_LOCATION: " << LOCATION << endl  
#define PRINT_NAME V3OUT << "PRINT_NAME: " << __PRETTY_FUNCTION__ << endl

#define ASSERT(assertion, exception) {                              \
   V3OUT << "ASSERTION: " << #assertion << endl;                    \
   if (!(assertion)) {                                              \
      OSTREAM << (LOCATION)                                         \
           << gettext("Assert: assertion violated.") << endl        \
           << #assertion << endl                                    \
           << gettext("Aborting.") << endl << flush;                \
      TERMINATION;                                                  \
   } }

#if VERBOSE >= 1
#define DEBUG_BREAK(aCondition, aString)                            \
   if (aCondition) {                                                \
      V1OUT << "DEBUG_BRK:  " << aString << endl << flush;          \
      getchar(); } 
#else
#define DEBUG_BREAK(a, b)
#endif

#define DEBUG_TRACE ObjectTrace _debug_trace_

#define DEBUG_MSG(aString) V1OUT << "DEBUG_MSG:  " << (aString) << endl
#define DUMP(x) V2OUT << "DEBUG_DUMP: " << #x << " == " << (x) << endl
#define PRECOND(a, b)   ASSERT(a, b)
#define POSTCOND(a, b)  ASSERT(a, b)
#define INVARIANT(a, b) ASSERT(a, b)

#else /* NDEBUG */

#define ASSERT(a, b)
#define PRECOND(a, b) 
#define POSTCOND(a, b) 
#define INVARIANT(a, b) 
#define DEBUG_MSG(a)
#define DEBUG_BREAK(a, b)
#define DEBUG_TRACE 
#define DUMP(a)
#define PRINT_LOCATION
#define PRINT_NAME

#endif /* NDEBUG */
#define exception string("")
#define ERROR(message, exception) {                                 \
  OSTREAM << (LOCATION)                                             \
       << gettext("Error: an error occured.") << endl               \
       << gettext(message) << endl                                  \
       << gettext("Aborting.") << endl << flush;                    \
  TERMINATION; }

#define ERROR_IF(condition, message, exception) {                   \
  if (condition) { ERROR(message, exception); } }

#endif /* __DEBUG_H__ */


