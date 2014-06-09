/**
* @Class:              IoTools
* @Author:             Silvio Tosatto
* @Project Name:       Victor
* @Description:        Library of auxiliary I/O functions.
*/

#ifndef _IO_TOOLS_H_

#define _IO_TOOLS_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector3.h>
#include <matrix3.h>
#include <string>
#include <String2Number.h>

/// Writes a line separator.
inline void 
fillLine(ostream& out)
{
  out << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

/// Checks if subsequent line is empty.
int checkForBlankLine(istream& is, bool returnbuffer = false);

/// Checks if keyword is present at current position in stream.
int checkForKeyword(istream& is, string key);

/// Reads a number (int, double, etc) from same line to result 
/// (and returns 1) or returns 0 if none present.
template<class T> int readOnSameLine(istream& is, T& result);

/// Skip whitespace.
void eatWhite(istream& is);
// Skips over comments (starting with '#').
void eatComment(istream& is);

/// Peek if number follows, if so return it (1) else not (0).
template<class T> int readNumber(istream& is, T& result);

/// Read a whole line into a string.
string readLine(istream& is);

/// Skips to new line.
void skipToNewLine(istream& is);

template<class T> ostream& operator <<(ostream &os, const vgVector3<T> &vec);
template<class T> istream& operator >>(istream &is, vgVector3<T> &vec);

template<class T> ostream& operator <<(ostream &os, const vgMatrix3<T> &vec);
template<class T> istream& operator >>(istream &is, vgMatrix3<T> &vec);

#endif
