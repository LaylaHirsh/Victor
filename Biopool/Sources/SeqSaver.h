/**
 * @Class:              SeqSaver
 * @Base Class(es):     Saver
 * @Author:             Silvio Tosatto
*@Description:
*    Loads components (Atoms, Groups, etc.) in SEQ format.
*    SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. 
*    Note: saveGroup() is not implemented, as it has no valid use.
*/

#ifndef _SEQ_SAVER_H_
#define _SEQ_SAVER_H_

// Includes:
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>
#include <Ligand.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief Loads components (Atoms, Groups, etc.) in SEQ format.
 * 
* @Description SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. 
* @This 
 * */
class SeqSaver : public Saver{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  SeqSaver(ostream& _output = cout) 
    : output(_output), writeChi(true) { }
  // this class uses the implicit copy operator.
  virtual ~SeqSaver() { PRINT_NAME; }  

// MODIFIERS:
  void setWriteChi(bool _w) { writeChi = _w; }
  virtual void saveSideChain(SideChain& sc, bool header = 1);
  virtual void saveAminoAcid(AminoAcid& aa);
  virtual void saveSpacer(Spacer& sp);
  virtual void saveLigand(Ligand& l);

protected:

private:
  ostream& output;   // output stream
  bool writeChi;     // switch: write chi angles?
};

} // namespace
#endif //_SEQ_SAVER_H_
