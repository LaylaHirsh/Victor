/**
* @Class:              Saver
* @Derived Class(es):  PdbSaver, XyzSaver, IntSaver, RelSaver, SeqSaver
* @Author:             Silvio Tosatto
* @Project Name:       Victor
*/
#ifndef _SAVER_H_
#define _SAVER_H_

// Includes:
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

class Group;
class SideChain;
class AminoAcid;
class Spacer;
class Ligand;
class LigandSet;
class Protein;
class Nucleotide;
/** @brief Base class for saving components (Atoms, Groups, etc.).
 * 
* @Description  
* @This 
 * */
class Saver{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  Saver() {};
  // this class uses the implicit copy operator.
  virtual ~Saver() {};  

// MODIFIERS:
  virtual void saveGroup(Group& group) {};
  virtual void saveSideChain(SideChain& node) {};
  virtual void saveAminoAcid(AminoAcid& node) {};
  virtual void saveSpacer(Spacer& node) {};
  virtual void saveLigand(Ligand& node) {};
  virtual void saveLigandSet(LigandSet& node) {};
  virtual void saveProtein(Protein& node) {};
  virtual void saveNucleotide(Nucleotide& node) {};
protected:

private:

};

} // namespace
#endif //_SAVER_H_










