/**
 * @Class:              SEQ Constructor
 * @Author:             Silvio Tosatto
 * @Description:
*    This class builds a spacer by concatenating the same aminoacid type for 
*    n times.
*/
#ifndef _SIDECHAIN_CONSTRUCTOR_H_
#define _SIDECHAIN_CONSTRUCTOR_H_

// Includes:
#include <Spacer.h>
#include <AminoAcid.h>
#include <Debug.h>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief This class builds a spacer by concatenating the same aminoacid type for 
*    n times.
 * 
* @Description  
* @This 
 * */
class SeqConstructor{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  SeqConstructor(istream& _refInput = cin) :
    refInput(_refInput), loaded(false), refAmino() { }
  virtual ~SeqConstructor() { PRINT_NAME; }  

// MODIFIERS:
  virtual Spacer& makeSpacer(string type, unsigned int n); 
  // make spacer by concatenating code n-times 

protected:
  // HELPERS:
  virtual void loadReference();
  void searchReference(AminoAcid& aa, string type);
  void buildAminoAcid(AminoAcid& aa, AminoAcid* prev);

private:
  istream& refInput;             // reference stream
  bool loaded;                   // is reference loaded yet?
  vector<AminoAcid> refAmino;    // reference types
};
/** @example BuildTest.cc
   *  A simple program to test class SeqConstructor's features.
 */
} // namespace
#endif //_SEQ_CONSTRUCTOR_H_
