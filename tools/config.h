/**
 *  @Author: Marcus Pr�mmer
*  @Class: config
*  @Description: This class makes it possible to handle all parameters over 
*               a simple plain text config file.
*               You can easily ask for parameters in the config file,
*               change the value while the program is running, create new
*               parameters, delete someone or set a new value of a parameter
*               and save it.
*	
*/ 
#ifndef _config_h_
#define _config_h_

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <typeinfo>
#include <String2Number.h>
#include <vector>
#include <vector3.h>

struct parameter {
  string name;
  string value;
  string file;

public:
  const parameter& operator= (const parameter& src) { name=src.name; value=src.value; file=src.file;  return *this; };

};

/** @brief class implements methods to convert strings into another type of data.And another data type to a string
* 
* @Description  
 **/
class config {
 
public:
config(string fileName);
config();
template<class Tset>    int  setParameter(string paramName,Tset value);
template<class Tget>    Tget getParameter(string paramName,Tget tmp);
template<class Tchange> int  changeParameter(string,Tchange);
template<class Tnew>    int  newParameter(string,string,Tnew);
int                          delParameter(string);
vector<parameter>            scan_configFile(string fileName);
int                          existParameter(string paramName);

 protected:

 private:
  int          convertString2Type(string value,int tmp);
  unsigned int convertString2Type(string value,unsigned int tmp);
  float        convertString2Type(string value,float tmp);
  long         convertString2Type(string value,long tmp);
  double       convertString2Type(string value,double tmp);
  string       convertString2Type(string value,string tmp);

  string convertTypes(int value);
  string convertTypes(unsigned int value);
  string convertTypes(double value);
  string convertTypes(float value);
  string convertTypes(long value);
  string convertTypes(string value);
};
/** @example ConfigTest.cc
   *  A simple program to test class Config's features.
 */
#endif

