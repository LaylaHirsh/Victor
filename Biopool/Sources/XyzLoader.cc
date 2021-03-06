/** Class:             XyzLoader
*@Author:            Silvio Tosatto
 * @Description:
*    Loads components (Atoms, Groups, etc.) in XYZ (carthesian) format.
*    Internal format is made of type, coords & bonds of each atom,
*    one per line. Keywords "aminoacid" and "sidechain" delimit these
*    structures.
*/

// Includes:
#include <XyzLoader.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:
/**
 * @Description Loads the information from the input stream into an atom group considering the xyz format. 
 * @param  Group reference(Group&)
 * @return changes are made internally
 */  
void XyzLoader::loadGroup(Group& gr){  
  PRINT_NAME;
  if (!input)
    return;

  gr.setType(readLine(input));   // read title line

  int n = -1;
  unsigned int i = 0;
  string typeName;
  vgVector3<double> _coord;
  while (readNumber(input, n)) {  // read all entries
    
      Atom at;
      input >> typeName >> _coord;
      at.setType(typeName);
      gr.addAtom(at);
      unsigned int size = gr.size()-1;
      while (readOnSameLine(input,i))     // read all connections
	if ((i-1) <= size) 
	  gr[size].bindIn(gr[i-1]);
      gr[size].setCoords(_coord);
    }
}


/**
 * @Description  Loads the side chain from the file, and returns the reference to the corresponding amino acid 
 * @param  SideChain reference(SideChain&), AminoAcid pointer(AminoAcid&)
 * @return void
 */  
 void XyzLoader::loadSideChain(SideChain& sc, AminoAcid* aaRef){
  PRINT_NAME;
  if (!input)
    return;

  loadGroup(sc);

  if (aaRef != NULL)
    sc.setBackboneRef(aaRef);
}

/**
 * @Description  Loads an aminoacid in xyz format. 
 *  
 * @param  AminoAcid reference
 * @return void
 */  
void XyzLoader::loadAminoAcid(AminoAcid& aa){
  PRINT_NAME;
  if (!input)
    return;

  loadGroup(aa);

  if (checkForKeyword(input, "sidechain"))    {
    loadSideChain(aa.getSideChain(), &aa);
    if (connect)
      if (aa.getSideChain().size())
	for (unsigned int i = 0; i < aa.size(); i++)
	  if (aa[i].getType().c_str()[0] == 'C') { // connect backbone
	      if (!aa[i].isOutBond(aa.getSideChain()[0]))
		  aa[i].bindOut(aa.getSideChain()[0]);
	      break;
	    }
    }
  else
    DEBUG_MSG("XyzLoader::loadAminoAcid: No sidechain found.");

  aa.adjustLeadingN();
}


/**
 * @Description Loads a spacer in xyz format. 
 *  
 * @param  Spacer reference
 * @return void
 */   
void XyzLoader::loadSpacer(Spacer& sp){
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


/**
 * @Description Loads a Ligand in xyz format. 
 *  
 * @param  Ligand reference
 * @return void
 */    
 
void XyzLoader::loadLigand(Ligand& l){
  PRINT_NAME;
  ERROR("Not implemented yet",exception);
}






