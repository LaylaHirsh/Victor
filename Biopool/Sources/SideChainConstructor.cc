/**
 * @Class:             SideChainConstructor
 * @Author:            Silvio Tosatto
 * @Description:
 * */

// Includes:
#include <SideChainConstructor.h>
#include <PdbLoader.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

void SideChainConstructor::loadReference(){
  PdbLoader pl(refInput);
  refSpacer.load(pl);
  loaded = true;
}

SideChain& SideChainConstructor::makeSideChain(string type){
  if (!loaded)
    loadReference();
  for (unsigned int i = 0; i < refSpacer.sizeAmino(); i++)
    if (refSpacer.getAmino(i).getSideChain().getType() == type)      {
	SideChain& sc = refSpacer.getAmino(i).getSideChain();
	return sc;
      }
  ERROR("No reference structure found.", exception);
  SideChain* sc = new SideChain;
  return *sc;
}

