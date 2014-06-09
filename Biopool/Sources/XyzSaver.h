/**
*@Class:              XyzSaver
*@Base Class(es):     Saver
*@Author:             Silvio Tosatto
*@Project Name:       Victor
*
*/

#ifndef _XYZ_SAVER_H_
#define _XYZ_SAVER_H_

// Includes:
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <Ligand.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief Saves components (Atoms, Groups, etc.) in carthesian format 
 * 
* @Description 
*    Internal format is made of type, coords & bonds of each atom,
*    one per line. Keywords "aminoacid" and "sidechain" delimit these
*    structures.
 * */
class XyzSaver : public Saver{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  XyzSaver(ostream& _output = cout, int _offset = 1) 
    : output(_output), delimit(true), offset(_offset) { }
  // this class uses the implicit copy operator.
  virtual ~XyzSaver() { PRINT_NAME; }  

// MODIFIERS:
  virtual void saveGroup(Group& gr);
  virtual void saveSideChain(SideChain& sc);
  virtual void saveAminoAcid(AminoAcid& aa);
  virtual void saveSpacer(Spacer& sp);
  virtual void saveLigand(Ligand& l);
  void setDelimit(bool _d) { delimit = _d; }

protected:
  virtual void pSaveAtomVector(vector<Atom>& va);

private:
  ostream& output;   // output stream
  bool delimit;      // write delimiters ("aminoacid", "sidechain", etc.)  
  int offset;        // ID offset for saving (optional) 
};
}
#endif //_XYZ_SAVER_H_
