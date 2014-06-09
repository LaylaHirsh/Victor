/**   
 * @Class:              IntSaver
 * @Author:             Silvio Tosatto
 * @Project Name:       Victor
*    Attention: This class is *NOT* finished yet!
*/

#ifndef _INT_SAVER_H_
#define _INT_SAVER_H_

// Includes:
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief Loads components (Atoms, Groups, etc.) in internal format.
 * 
* @Description Internal format is defined by listing type, bond length partner &
*    bond length, bond angle partner & bond angle, torsion angle partner
*    & torsion angle plus a chirality (0 if it is a 'true' torsion angle,
*    +1 or -1 if the 'torsion angle' is a second bond angle), for each
*    atom, one per line.
*    NB: Only chirality 0 is currently supported.
* 
 **/
class IntSaver : public Saver
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  IntSaver(ostream& _output = cout) : output(_output) { };
  // this class uses the implicit copy operator.
  virtual ~IntSaver() { PRINT_NAME; };  

// MODIFIERS:
  virtual void saveGroup(Group& group);
  virtual void saveSideChain(SideChain& node);
  virtual void saveAminoAcid(AminoAcid& node);
  virtual void saveSpacer(Spacer& sp);
  virtual void saveLigand(Ligand& l);

protected:
  void pSaveAtomVector(vector<Atom>& va, unsigned int offset = 0);

private:
  ostream& output;   // output stream
  // unsigned int count;
};

} // namespace
#endif //_INT_SAVER_H_
