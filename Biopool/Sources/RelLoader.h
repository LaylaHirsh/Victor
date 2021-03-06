/**
* @Class:              RelLoader
* @Base Class(es):     Loader
* @Author:             Silvio Tosatto
* @Project Name:       Victor
*/

#ifndef _REL_LOADER_H_
#define _REL_LOADER_H_

// Includes:
#include <Group.h>
#include <SideChain.h>
#include <Loader.h>
#include <AminoAcid.h>
#include <Spacer.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief Loads components (Atoms, Groups, etc.) in relative format.
 * 
* @Description Relative format is similiar in structure to XYZ format.
*    The only difference is that the coordinates here are relative 
*    to the previous atom rather absolute.
* @This 
 * */
class RelLoader : public Loader
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  RelLoader(istream& _input = cin) : input(_input), connect(true) { }
  // this class uses the implicit copy operator.
  virtual ~RelLoader() { PRINT_NAME; }  

// MODIFIERS:
  void connectSegments(bool c) { connect = c; }
  // determines if segments (aminoacids, sidechains) have to be connected
  virtual void loadGroup(Group& gr);
  virtual void loadSideChain(SideChain& sc, AminoAcid* aaRef = NULL);
  virtual void loadAminoAcid(AminoAcid& aa);
  virtual void loadSpacer(Spacer& sp);

protected:

private:
  istream& input;   // input stream
  bool connect; // are segments to be connected to each other?
};

} // namespace
#endif //_REL_LOADER_H_
