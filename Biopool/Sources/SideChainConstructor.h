/**
 * @Class:              SideChainConstructor
 * @Author:             Silvio Tosatto
*/
#ifndef _SIDECHAIN_CONSTRUCTOR_H_
#define _SIDECHAIN_CONSTRUCTOR_H_

// Includes:
#include <Spacer.h>
#include <SideChain.h>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief base reference stream
* @This 
 * */
class SideChainConstructor {
public: 

// CONSTRUCTORS/DESTRUCTOR:
  SideChainConstructor(istream& _refInput = cin) :
    refInput(_refInput), loaded(false), refSpacer() { }
  virtual ~SideChainConstructor() { PRINT_NAME; }  

// MODIFIERS:
  virtual SideChain& makeSideChain(string code); 

protected:
  // HELPERS:
  virtual void loadReference();

private:
  istream& refInput;             // reference stream
  int loaded;                    
  Spacer refSpacer;  // reference types
};

} // namespace
#endif //_SIDECHAIN_CONSTRUCTOR_H_
