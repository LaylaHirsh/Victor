/**
* @Class:             RelLoader
* @Base Class(es):    Loader
* @Derived Class(es): -
* @Containing:        -
* @Author:            Silvio Tosatto
* @Project Name:       Victor
* @Description:
*    Loads components (Atoms, Groups, etc.) in relative format.
*    Relative format is similiar in structure to XYZ format.
*    The only difference is that the coordinates here are relative 
*    to the previous atom rather absolute.
*/

// Includes:
#include <RelLoader.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        XyzLoader::loadGroup
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Loads a group in relative format. 
//
// ----------------------------------------------------------------------------
void RelLoader::loadGroup(Group& gr){  
  PRINT_NAME;
  if (!input)
    return;

  gr.setType(readLine(input));   // read title line

  int n = -1;
  unsigned int i = 0;
  string typeName;
  vgVector3<double> _trans;
  while (readNumber(input, n))  { // read all entries
      Atom at;
      input >> typeName >> _trans;
      at.setType(typeName);
      gr.addAtom(at);
      unsigned int size = gr.size()-1;
      while (readOnSameLine(input,i))     // read all connections
	if ((i-1) <= size) 
	  gr[size].bindIn(gr[i-1]);
      gr[size].setTrans(_trans);
    }
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        XyzLoader::loadSideChain
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Loads a sidechain in relative format. 
//
// ----------------------------------------------------------------------------
void RelLoader::loadSideChain(SideChain& sc, AminoAcid* aaRef){
  PRINT_NAME;
  if (!input)
    return;

  loadGroup(sc);

  if (aaRef != NULL)
    sc.setBackboneRef(aaRef);
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        XyzLoader::loadAminoAcid
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Loads an aminoacid in relative format. 
//
// ----------------------------------------------------------------------------
void RelLoader::loadAminoAcid(AminoAcid& aa){
  PRINT_NAME;
  if (!input)
    return;

  loadGroup(aa);

  if (checkForKeyword(input, "sidechain"))    {
    loadSideChain(aa.getSideChain(), &aa);
    if (connect)
      if (aa.getSideChain().size())
	for (unsigned int i = 0; i < aa.size(); i++)
	  if (aa[i].getType().c_str()[0] == 'C')   {
	      aa[i].bindOut(aa.getSideChain()[0]);
	      break;
	    }
    }
  else
    DEBUG_MSG("RelLoader::loadAminoAcid: No sidechain found.");
  
  aa.adjustLeadingN();
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        XyzLoader::loadSpacer
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Loads a spacer in relative format. 
//
// ----------------------------------------------------------------------------
void RelLoader::loadSpacer(Spacer& sp){
  PRINT_NAME;
  if (!input)
    return;

  sp.setType(readLine(input));   // read title line
  while (checkForKeyword(input, "aminoacid"))    {
      AminoAcid* aa = new AminoAcid();
      loadAminoAcid(*aa);
      // connect aa to previous chain segment
      if (connect)
	aa->bindIn((*aa)[N], sp.getAmino(sp.size()-1),
		   sp.getAmino(sp.size()-1)[C]);
      sp.insertComponent(aa);
    }
  if (sp.sizeAmino())
    sp.getAmino(0).adjustLeadingN();
}

