// stringtools.h

#ifndef __stringtools_H__
#define __stringtools_H__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


void strip(std::string &str);

std::string strip_and_return_string(std::string);

int strip_and_return_int(std::string);

unsigned int strip_and_return_unsigned_int(std::string);

float strip_and_return_float(std::string);

std::vector<std::string> split(std::string, char);

int seq_length(const std::string);

// template<class T>
// std::string to_string(const T&);

template<class T> std::string to_string(const T &x)
{
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

std::string StringToLower(std::string);
std::string StringToUpper(std::string);

#endif
