/** 
* @Class:              Loader
* @Derived Class(es):  PdbLoader, XyzLoader, IntLoader, RelLoader, SeqLoader
* @Author:             Silvio Tosatto
* @Project Name:       Victor
*/
#ifndef _LOADER_H_
#define _LOADER_H_

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
/** @brief Base class for loading components (Atoms, Groups, etc.).
 * 
* @Description  
* @This 
 * */
class Loader
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  Loader()  { PRINT_NAME; }
  // this class uses the implicit copy operator.
  virtual ~Loader() { PRINT_NAME; };  

// MODIFIERS:
  virtual void loadGroup(Group& group) { PRINT_NAME; };
  virtual void loadSideChain(SideChain& node) { PRINT_NAME; };
  virtual void loadAminoAcid(AminoAcid& node) { PRINT_NAME; };
  virtual void loadSpacer(Spacer& node) { PRINT_NAME; };
  virtual void loadLigand(Ligand& node) { PRINT_NAME; };
  virtual void loadLigandSet(LigandSet& node) { PRINT_NAME; };
  virtual void loadProtein(Protein& node) { PRINT_NAME; };
  virtual void loadNucleotide(Nucleotide& node) { PRINT_NAME; };
protected:

private:

};


} // namespace
#endif //_LOADER_H_









