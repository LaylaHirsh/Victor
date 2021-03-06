/**
* @Class:              PdbSaver
* @Base Class(es):     Saver
* @Author:             Silvio Tosatto
* @Project Name:       Victor
*/

#ifndef _PDB_SAVER_H_
#define _PDB_SAVER_H_

// Includes:
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <LigandSet.h>
#include <Ligand.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>
#include <Protein.h>
// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief Saves components (Atoms, Groups, etc.) in standard PDB format
 * 
* @Description 
* @This 
 * */
class PdbSaver : public Saver
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  PdbSaver(ostream& _output = cout) 
    : output(_output), writeSeq(true), writeSecStr(true), writeTer(true), 
      atomOffset(0), aminoOffset(0), ligandOffset(0), chain(' ') { }
  // this class uses the implicit copy operator.
  virtual ~PdbSaver() { PRINT_NAME; }  

// PREDICATES:
    void endFile() { output << "END\n"; }

// MODIFIERS:
  void setWriteSecondaryStructure() { writeSecStr = true; }
  void setDoNotWriteSecondaryStructure() { writeSecStr = false; }

  void setWriteSeqRes() { writeSeq = true; }
  void setDoNotWriteSeqRes() { writeSeq = false; }

  void setWriteAtomOnly() { writeSecStr = false; writeSeq = false; 
                            writeTer = false; }
  void setWriteAll()  { writeSecStr = true; writeSeq = true; 
                        writeTer = true; }
  void setChain(char _ch) { chain = _ch; }
  virtual void saveGroup(Group& gr);
  virtual void saveSideChain(SideChain& sc);
  virtual void saveAminoAcid(AminoAcid& aa);
  virtual void saveSpacer(Spacer& sp);
  virtual void saveLigand(Ligand& l);
  virtual void saveLigandSet(LigandSet& l);
  virtual void saveProtein(Protein& prot);

protected:

private:

  // HELPERS:
  void writeSeqRes(Spacer& sp); // writes SEQRES entry
  void writeSecondary(Spacer& sp); 
     // writes secondary entries (SHEET, HELIX, etc.)
  // ATTRIBUTES 
  ostream& output;   // output stream
  bool writeSeq, writeSecStr, writeTer;  
  unsigned int atomOffset, ligandOffset;
  int aminoOffset;
  char chain;      // chain ID
  // offsets that determine at which atom, aminoacid and ligand number to start
};

} // namespace
#endif //_PDB_SAVER_H_


