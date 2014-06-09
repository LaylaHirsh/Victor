/**
 * @Class:              SeqLoader
 * @Base Class(es):     Loader
 * @Author:             Silvio Tosatto
 * @Description:
*    Loads components (Atoms, Groups, etc.) in SEQ format.
*    SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. In order to construct the protein, a 
*    reference file with sample aminoacids has to be loaded first.
*/

#ifndef _SEQ_LOADER_H_
#define _SEQ_LOADER_H_

// Includes:
#include <SideChain.h>
#include <AminoAcid.h>
#include <Loader.h>
#include <Debug.h>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief Loads components (Atoms, Groups, etc.) in SEQ format. 
 * 
* @Description SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. In order to construct the protein, a 
*    reference file with sample aminoacids has to be loaded first. .
* @This 
 * */
class SeqLoader : public Loader{
public: 
  
  // CONSTRUCTORS/DESTRUCTOR:
  SeqLoader(istream& _input = cin, istream& _refInput = cin) :
    input(_input), refInput(_refInput), loaded(false), refAmino() { }
  // this class uses the implicit copy operator.
  ~SeqLoader() { PRINT_NAME; }

  // MODIFIERS:
  virtual void loadAminoAcid(AminoAcid& node, AminoAcid* prev = NULL);
  virtual void loadSpacer(Spacer& node);
  virtual void loadLigand(Ligand& node);
  
protected:
  // HELPERS:
  virtual void loadReference();
  virtual void setStructure(AminoAcid& aa, string type);
  
  // ATTRIBUTES:
  istream& input;             // input stream
  istream& refInput;          // reference stream
  int loaded;                 // is reference loaded yet?
  vector<AminoAcid> refAmino; // reference amino acids

private:

};

} // namespace
#endif //_SEQ_LOADER_H_
