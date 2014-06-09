/**
 * @Class:             XyzSaver
 * @Author:            Silvio Tosatto
 * @Description:
*    Loads components (Atoms, Groups, etc.) in SEQ format.
*    SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. 
*    Note: saveGroup() is not implemented, as it has no valid use.
*/
// Includes:
#include <SeqSaver.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveSideChain
//
//  Author:        Silvio Tosatto 
//
//  Date:          09/99
//
//  Description:
//    Saves a sidechain in SEQ format. 
//
// ----------------------------------------------------------------------------
void SeqSaver::saveSideChain(SideChain& sc, bool header){
  PRINT_NAME;
  unsigned old_prec = output.precision();
  ios::fmtflags old_flags = output.flags();
  output.setf(ios::fixed, ios::floatfield);
  if (header)
    output << sc.getType() << "   ";

  if (writeChi)
    for (unsigned int i = 0; i < sc.getMaxChi(); i++) // write torsion angles
      output << "   " << setw(8) << setprecision(3) << sc.getChi(i);
  output << "\n";
  output.precision(old_prec);
  output.flags(old_flags);
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveAminoAcid
//
//  Author:        Silvio Tosatto 
//
//  Date:          09/99
//
//  Description:
//    Saves an aminoacid in SEQ format. 
//
// ----------------------------------------------------------------------------
void SeqSaver::saveAminoAcid(AminoAcid& aa){
  PRINT_NAME;
  unsigned old_prec = output.precision();
  ios::fmtflags old_flags = output.flags();
  output.setf(ios::fixed, ios::floatfield);
  output << aa.getType() << "   " // write torsion angles
	 << setw(8) << setprecision(3) << aa.getPhi() 
	 << "   " << setw(8) << setprecision(3) << aa.getPsi() 
	 << "   " << setw(8) << setprecision(3) << aa.getOmega();
  saveSideChain(aa.getSideChain(), false);
  output.precision(old_prec);
  output.flags(old_flags);
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveSpacer
//
//  Author:        Silvio Tosatto 
//
//  Date:          09/99
//
//  Description:
//    Saves a spacer in SEQ format. 
//
// ----------------------------------------------------------------------------

void SeqSaver::saveSpacer(Spacer& sp){
  PRINT_NAME;
  output << sp.getType() << "\n";
  for (unsigned int i = 0; i < sp.size(); i++)
    sp[i].save(*this);
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveLigand
//
//  Author:        Marcus Pruemmer
//
//  Date:          06/2000
//
//  Description:
//    Saves a Ligand in SEQ format. 
//
// ----------------------------------------------------------------------------
void SeqSaver::saveLigand(Ligand& l){
  PRINT_NAME;
  ERROR("Not implemented yet",exception);
}
