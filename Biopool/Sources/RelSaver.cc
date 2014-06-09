/**
* @Class:             RelSaver
* @Base Class(es):    Saver
* @Derived Class(es): -
* @Containing:        -
* @Author:            Silvio Tosatto
* @Project Name:      Victor
* @Description:
*    Loads components (Atoms, Groups, etc.) in relative format.
*    Relative format is similiar in structure to XYZ format.
*    The only difference is that the coordinates here are relative 
*    to the previous atom rather absolute.
*/


// Includes:
#include <RelSaver.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        RelSaver::saveGroup
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Saves a group in relative format. 
//
// ----------------------------------------------------------------------------
void RelSaver::saveGroup(Group& gr){  
  PRINT_NAME;
  gr.sync();
  output << gr.getType() << "\n";
  pSaveAtomVector(gr.giveAtoms()); 
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        RelSaver::saveSideChain
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Saves a sidechain in relative format. 
//
// ----------------------------------------------------------------------------
void RelSaver::saveSideChain(SideChain& sc){
  PRINT_NAME;
  sc.sync();
  output << sc.getType() << "\n";
  pSaveAtomVector(sc.giveAtoms()); 
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        RelSaver::saveAminoAcid
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Saves an aminoacid in relative format. 
//
// ----------------------------------------------------------------------------
void RelSaver::saveAminoAcid(AminoAcid& aa){
  PRINT_NAME;
  aa.sync();
  output << aa.getType() << "\n";
  pSaveAtomVector(aa.giveAtoms()); 

  // write relative position of atom following C, if it exists:
  // can be implemented more efficiently..
  if (aa.isMember(C) && !aa.isMember(OXT))
    for (unsigned int i = 0; i < aa[C].sizeOutBonds(); i++)
      if (aa[C].getOutBond(i).getCode() == N)
	{
  	  output << "  " << setw(4) << aa[C].getOutBond(i).getNumber() 
		 << "    OXT  " << aa[C].getOutBond(i).getTrans() 
		 << "      " << setw(3) << aa[C].getNumber() << "\n";
	  break;
	}
  output << "  sidechain\n  ";
  saveSideChain(aa.getSideChain());
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        RelSaver::saveSpacer
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Saves a spacer in relative format. 
//
// ----------------------------------------------------------------------------
void RelSaver::saveSpacer(Spacer& sp){
  PRINT_NAME;
  output << sp.getType() << "\n";
  for (unsigned int i = 0; i < sp.size(); i++)
    {
      output << "aminoacid\n";
      sp[i].save(*this);
    }
}


// HELPER:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        RelSaver::pSaveAtomVector
//
//  Author:        Silvio Tosatto 
//
//  Date:          10/99
//
//  Description:
//    Private helper method to save a vector of atoms.
//    Attention: ID realignment is currently commented out.
//
// ----------------------------------------------------------------------------
void RelSaver::pSaveAtomVector(vector<Atom>& va){
  unsigned old_prec = output.precision();
  ios::fmtflags old_flags = output.flags();
  output.setf(ios::fixed, ios::floatfield);

  for (unsigned int k = 0; k < va.size(); k++)  // write all entries
    {
      string atName = va[k].getType();
      if (!isdigit(atName[0]))
	atName = ' ' + atName;
      while (atName.size() < 4)
	atName += ' ';

      output << "  " << setw(4) << va[k].getNumber() << "   " << atName << "  "
	     << va[k].getTrans() << "   ";
      for (unsigned int i = 0; i < va[k].sizeInBonds(); i++)
	  output << "   " << setw(3) << va[k].getInBond(i).getNumber();
      for (unsigned int i = 0; i < va[k].sizeOutBonds(); i++)

	  output << "   " << setw(3) << va[k].getOutBond(i).getNumber();
      output << "\n";
    }

  output.precision(old_prec);
  output.flags(old_flags);
}
