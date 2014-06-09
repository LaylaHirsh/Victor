/**
*@Class:             XyzSaver
*@Base Class(es):    Saver
*@Author:            Silvio Tosatto
*@Project Name:      Victor
*@Description:
*    Saves components (Atoms, Groups, etc.) in XYZ (carthesian) format.
*    Internal format is made of type, coords & bonds of each atom,
*    one per line. Keywords "aminoacid" and "sidechain" delimit these
*    structures.
*/

// Includes:
#include <XyzSaver.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

/**
 * @Description  Obtains the atom from the given group and saves it as a vector format. 
 * @param   Group reference(Group&)
 * @return  changes are made internally (void)
 */

void XyzSaver::saveGroup(Group& gr){  
  PRINT_NAME;
  gr.sync();
  output << gr.getType() << "\n";
  pSaveAtomVector(gr.giveAtoms()); 
}

/**
 * @Description  Obtains the atom from the given side chain and saves it as a vector format. 
 * @param   Side chain reference(SideChain&)
 * @return  changes are made internally (void)
 */
 
void XyzSaver::saveSideChain(SideChain& sc){
  PRINT_NAME;
  sc.sync();
  if (delimit)
    output << sc.getType() << "\n";
  pSaveAtomVector(sc.giveAtoms()); 
}

/**
 * @Description  Obtains the atom and the said chain from the given amino acid and saves them as a vector format. 
 * @param  Amino acid reference(AminoAcid&)
 * @return  changes are made internally (void)
 */
void XyzSaver::saveAminoAcid(AminoAcid& aa){
  PRINT_NAME;
  aa.sync();
  if (delimit)
    output << aa.getType() << "\n";
  pSaveAtomVector(aa.giveAtoms()); 
  if (delimit)
    output << "  sidechain\n  ";
  saveSideChain(aa.getSideChain());
}

/**
 * @Description  Obtains the atoms and side chains  from each of the amino acids in the spacer/protein and saves them  as a vector format. 
 * @param   Spacer reference(Spacer&)
 * @return  changes are made internally (void)
 */
      
void XyzSaver::saveSpacer(Spacer& sp){
  PRINT_NAME;
  if (delimit){
    output << sp.getType() << "\n";
  }
  else    {  // write number of atoms
      unsigned int size = 0;
      for (unsigned int i = 0; i < sp.sizeAmino(); i++)
	size += sp.getAmino(i).size();
      output << size << "\n\n";
    }
  for (unsigned int i = 0; i < sp.size(); i++)    {
      if (delimit)
	output << "aminoacid\n";
      sp[i].save(*this);
    }
}
void XyzSaver::saveLigand(Ligand& l){
  PRINT_NAME;
  ERROR("Not implemented yet",exception);
}

// HELPER:

/**
 * @Description  Creates a list of the atom types contained in the vector and saves it.
 * @param   reference to the vector that contains atoms(vector<Atom>&)
 * @return  changes are made internally (void)
 */
 
void XyzSaver::pSaveAtomVector(vector<Atom>& va){  // warning: don't copy atom vector as it would lose the original bonds
  unsigned old_prec = output.precision();
  ios::fmtflags old_flags = output.flags();
  output.setf(ios::fixed, ios::floatfield);
 
  for (unsigned int k = 0; k < va.size(); k++) { // write all entries
      string atName = va[k].getType();
      if (!isdigit(atName[0]))
	atName = ' ' + atName;
      while (atName.size() < 4)
	atName += ' ';

      if (delimit)
	output << "  " << setw(4) << va[k].getNumber() << "   ";
      output << atName << "  " << va[k].getCoords() << "   ";
      for (unsigned int i = 0; i < va[k].sizeInBonds(); i++)
 	  output << "   " << setw(3) << va[k].getInBond(i).getNumber();
      for (unsigned int i = 0; i < va[k].sizeOutBonds(); i++)
 	  output << "   " << setw(3) << va[k].getOutBond(i).getNumber();
      output << "\n";
    }
   output.precision(old_prec);
  output.flags(old_flags);
}