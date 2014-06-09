/**
 * @Class:              SEQ Constructor
 * @Author:             Silvio Tosatto
 * @Project Name:       Victor
 * @Description:
*    This class builds a spacer by concatenating the same aminoacid type for 
*    n times.
*/

// Includes:
#include <SeqConstructor.h>
#include <RelLoader.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqConstructor::loadReference()
//
//  Author:        Silvio Tosatto 
//
//  Date:          11/99
//
//  Description:
//    Helper function to load the reference file. (same as in SeqLoader)
//
// ----------------------------------------------------------------------------
void SeqConstructor::loadReference(){
  RelLoader relL(refInput);
  relL.connectSegments(false);

  while (!loaded)    {
      AminoAcid* tAA = new AminoAcid;
      tAA->load(relL);
      if (tAA->size())
	refAmino.push_back(*tAA);
      else
	loaded = true;
    }
}


Spacer& SeqConstructor::makeSpacer(string type, unsigned int n){
  if (!loaded)
    loadReference();

  Spacer* sp = new Spacer;
  AminoAcid aa;
  searchReference(aa, type);

  for (unsigned int i = 0; i < n; i++)    {
      AminoAcid* tmpAA = new AminoAcid(aa);
      AminoAcid* tmpPrev = i ? &(sp->getAmino(i-1)) : NULL;
      buildAminoAcid(*tmpAA, tmpPrev);
      sp->insertComponent(tmpAA);
    }
  
  return *sp;
}


void SeqConstructor::buildAminoAcid(AminoAcid& aa, AminoAcid* prev){
  IntCoordConverter icc;

  aa.setBondsFromPdbCode(false);

  if (prev != NULL)
    icc.connectStructure(aa, *prev);  // connect aa's structure with pred

  // apply standard torsion angles
  aa.setPhi(-63);
  aa.setPsi(-42);
  aa.setOmega(180);
  aa.sync();
}


void SeqConstructor::searchReference(AminoAcid& aa, string type){
  for (unsigned int i = 0; i < refAmino.size(); i++)
    if (refAmino[i].getType() == type)      {
	aa = refAmino[i];
	return;
      }

  ERROR("No reference structure found.", exception);
}
