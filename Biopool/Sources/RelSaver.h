/**
* @Class:              RelSaver
* @Base Class(es):     Saver
* @Author:             Silvio Tosatto
* @Project Name:       Victor
*/

#ifndef _REL_SAVER_H_
#define _REL_SAVER_H_

// Includes:
#include <vector>
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief Loads components (Atoms, Groups, etc.) in relative format 
 * 
* @Description Relative format is similiar in structure to XYZ format.
*    The only difference is that the coordinates here are relative 
*    to the previous atom rather absolute.
* @This 
 * */
class RelSaver : public Saver
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  RelSaver(ostream& _output = cout, int _offset = 1) 
    : output(_output), offset(_offset) { }
  // this class uses the implicit copy operator.
  virtual ~RelSaver() { PRINT_NAME; }  

// MODIFIERS:
  virtual void saveGroup(Group& gr);
  virtual void saveSideChain(SideChain& sc);
  virtual void saveAminoAcid(AminoAcid& aa);
  virtual void saveSpacer(Spacer& sp);

protected:
  virtual void pSaveAtomVector(vector<Atom>& va);

private:
  ostream& output;   // output stream
  int offset;        // ID offset for saving (optional) -- currently disabled
};

} // namespace
#endif //_REL_SAVER_H_
