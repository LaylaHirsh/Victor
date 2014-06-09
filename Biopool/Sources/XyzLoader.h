/**
 * @Class:              XyzLoader
* @Base Class(es):     Loader
@Author:             Silvio Tosatto
@Description:
*    Loads components (Atoms, Groups, etc.) in XYZ (carthesian) format.
*    Internal format is made of type, coords & bonds of each atom,
*    one per line. Keywords "aminoacid" and "sidechain" delimit these
*    structures.
*
*/

#ifndef _XYZ_LOADER_H_
#define _XYZ_LOADER_H_

// Includes:
#include <Group.h>
#include <Ligand.h>
#include <SideChain.h>
#include <Loader.h>
#include <AminoAcid.h>
#include <Spacer.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief class allows you to load different components in XYZ (carthesian) format
 * 
* @Description Includes methods that allow to Loads components (Atoms, Groups, etc.) in XYZ (carthesian) format.
*    Internal format is made of type, coords & bonds of each atom,
*    one per line. Keywords "aminoacid" and "sidechain" delimit these
*    structures.
* @This 
 * */
class XyzLoader : public Loader{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  XyzLoader(istream& _input = cin) : input(_input), connect(true) { }
  // this class uses the implicit copy operator.
  virtual ~XyzLoader() { PRINT_NAME; }  

// MODIFIERS:
  void connectSegments(bool c) { connect = c; } 
  // determines if segments (aminoacids, sidechains) have to be connected
  virtual void loadGroup(Group& gr);
  virtual void loadSideChain(SideChain& sc, AminoAcid* aaRef = NULL);
  virtual void loadAminoAcid(AminoAcid& aa);
  virtual void loadSpacer(Spacer& sp);
  virtual void loadLigand(Ligand& l);

protected:

private:
  istream& input;   // input stream
  bool connect; // are segments to be connected to each other?
};
} // namespace
#endif //_XYZ_LOADER_H_
