/**@Class:             SeqLoader
 * @Base Class(es):    Loader
 * @Author:            Silvio Tosatto
 * @Description:
*    Loads components (Atoms, Groups, etc.) in SEQ format.
*    SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. In order to construct the protein, a 
*    reference file with sample aminoacids has to be loaded first.
*/
// Includes:
#include <SeqLoader.h>
#include <RelLoader.h>
#include <IoTools.h>
#include <vector3.h>
#include <IntCoordConverter.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqLoader::loadAminoAcid
//
//  Author:        Silvio Tosatto 
//
//  Date:          11/99
//
//  Description:
//    Loads an aminoacid in SEQ format. 
//
// ----------------------------------------------------------------------------
void SeqLoader::loadAminoAcid(AminoAcid& aa, AminoAcid* prev){
  PRINT_NAME;
  if (!loaded)
    loadReference();
  IntCoordConverter icc;
  string type;
  double tPhi, tPsi, tOmega;
  input >> type >> tPhi >> tPsi >> tOmega;  // read an entry
  aa.setType(type);
  aa.getSideChain().setType(type);
  aa.getSideChain().setBackboneRef(&aa);
  setStructure(aa, type);  // set structure data as per reference file

  aa.setBondsFromPdbCode(false);
  
  if (prev != NULL)
    icc.connectStructure(aa, *prev);  // connect aa's structure with pred
  
  // apply torsion angles read from file
  if (tPhi < 990)     
    aa.setPhi(tPhi);
  if (tPsi < 990)
    aa.setPsi(tPsi);
  if (tOmega < 990)
    aa.setOmega(tOmega);
  unsigned int i = 0;
  double tmp;
  while (readOnSameLine(input, tmp))    {
      if ((tmp < 990) && (i < aa.getSideChain().getMaxChi()))
    	aa.setChi(i, tmp);
      i++;
    }
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqLoader::loadSpacer
//
//  Author:        Silvio Tosatto 
//
//  Date:          11/99
//
//  Description:
//    Loads a spacer in SEQ format. 
//
// ----------------------------------------------------------------------------
void SeqLoader::loadSpacer(Spacer& sp){
  PRINT_NAME;
  if (!loaded)   // load reference file if it hasn't been loaded yet
    loadReference();
  string type = readLine(input);
  sp.setType(type);

  while (input) { // load each aminoacid    
      AminoAcid* aa = new AminoAcid;
      AminoAcid* tmp = sp.size() ? &(sp.getAmino(sp.size()-1)) : NULL;
      loadAminoAcid(*aa, tmp );
      if (aa->size())
	sp.insertComponent(aa);
      eatComment(input);
    }

}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqLoader::loadLigand
//
//
//  Date:          06/2000
//
//  Description:
//    Loads a Ligand in SEQ format. 
//
// ----------------------------------------------------------------------------
void SeqLoader::loadLigand(Ligand& l){
  PRINT_NAME;
  ERROR("Not implemented yet",exception);
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqLoader::loadReference
//
//  Author:        Silvio Tosatto 
//
//  Date:          11/99
//
//  Description:
//    Private helper function to load the reference file.
//
// ----------------------------------------------------------------------------
void SeqLoader::loadReference(){
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


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqLoader::setStructure
//
//  Author:        Silvio Tosatto 
//
//  Date:          11/99
//
//  Description:
//    Private helper function to set aa to the correct structure as per 
//    reference file.
//
// ----------------------------------------------------------------------------
void SeqLoader::setStructure(AminoAcid& aa, string type){
  for (unsigned int i = 0; i < refAmino.size(); i++)
    if (refAmino[i].getType() == type)      {
	aa = refAmino[i];
	return;
      }
  ERROR("No reference structure found.", exception);
}

