/**   
* @Class:             IntSaver
* @Base Class(es):    Saver
*  @Author:            Silvio Tosatto
* @Project Name:      Victor
* @Description:
 *   Loads components (Atoms, Groups, etc.) in internal format.
*    Internal format is defined by listing type, bond length partner &
*    bond length, bond angle partner & bond angle, torsion angle partner
*    & torsion angle plus a chirality (0 if it is a 'true' torsion angle,
*    +1 or -1 if the 'torsion angle' is a second bond angle), for each
*    atom, one per line.
*    NB: Only chirality 0 is currently supported.
*
*    Attention: This class is *NOT* finished yet!
*/

// Includes:
#include <IntSaver.h>
#include <IoTools.h>
#include <vector3.h>
#include <IntCoordConverter.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

 /**
 * @Description Saves a group in INT format. 
 * @param group reference
 * @return  void 
 */   
 void IntSaver::saveGroup(Group& gr)
{  
  gr.sync();
  pSaveAtomVector(gr.giveAtoms(),0);
}

 /**
 * @Description   Saves a Ligand in INT format. 
 * @param Ligand reference
 * @return  void
 */     
 void IntSaver::saveLigand(Ligand& l)
{  
  ERROR("Not implemented yet",exception);
}

 /**
 * @Description   Saves a sidechain in INT format. 
 * @param SideChain reference
 * @return  void
 */     
 void IntSaver::saveSideChain(SideChain& sc)
{
  sc.sync();
  pSaveAtomVector(sc.giveAtoms(), 0);
}

 /**
 * @Description   Saves an aminoacid in INT format. 
 * @param AminoAcid reference
 * @return  void
 */     
 void IntSaver::saveAminoAcid(AminoAcid& aa)
{
  aa.sync();
  pSaveAtomVector(aa.giveAtoms(), 0); 
      saveSideChain(aa.getSideChain());
}

 /**
 * @Description    Saves a spacer in INT format. 
 * @param Spacer reference
 * @return  void
 */   

    void IntSaver::saveSpacer(Spacer& sp)
{
  sp.sync();
  Spacer tmpSp = sp;   // make a local copy for realignment
  unsigned int count = 0;
  for(unsigned int i = 0; i < tmpSp.sizeAmino(); i++)
    {
      for(unsigned int j = 0; j < tmpSp.getAmino(i).sizeBackbone(); j++)
	tmpSp.getAmino(i)[j].setNumber(++count);	  
      for(unsigned int j = 0; j < (tmpSp.getAmino(i).size() 
				   - tmpSp.getAmino(i).sizeBackbone()); j++)
	  tmpSp.getAmino(i).getSideChain()[j].setNumber(++count);	  
    }
  output << setw(4) << count << "  " << tmpSp.getType() << "\n";
  for (unsigned int i = 0; i < tmpSp.size(); i++)
    tmpSp[i].save(*this);
}


// HELPER:

 /**
 * @Description    Private helper method to save a vector of atoms.
*    Attention: ID realignment is currently commented out.
 * @param vector<Atom> reference, unsigned int
 * @return  void
 */   
 void IntSaver::pSaveAtomVector(vector<Atom>& va, unsigned int offset)
{  // warning: don't copy atom vector as it would lose the original bonds
  unsigned old_prec = output.precision();
  ios::fmtflags old_flags = output.flags();
  output.setf(ios::fixed, ios::floatfield);
  IntCoordConverter icc;
  
  for (unsigned int k = 0; k < va.size(); k++)  // write all entries
    {
      string atName = va[k].getType();
      if (!isdigit(atName[0]))
	atName = ' ' + atName;
      while (atName.size() < 4)
	atName += ' ';
      
      output << setw(4) << va[k].getNumber() << "   " << atName << "  "
	     << setw(3) << static_cast<int>(va[k].getCode()); 
      if (!va[k].sizeInBonds())
	{
	  output << "\n";
	  continue;
	}
      output << "  " << setw(4) << va[k].getInBond(0).getNumber() 
	     << "  " << setw(8) << setprecision(5) 
	     << icc.getBondLength(va[k].getInBond(0), va[k]);
      if (va[k].getInBond(0).sizeInBonds())
	{
	  output << "  " << setw(4) 
		 << va[k].getInBond(0).getInBond(0).getNumber() 
		 << "  " << setw(8) << setprecision(4) 
		 << RAD2DEG * 
	    icc.getBondAngle(va[k].getInBond(0).getInBond(0), 
			     va[k].getInBond(0), va[k]);
	  
	  if ( (va[k].getType() == "CB" ) || (va[k].getType() == "1HA") ||
	       (va[k].getType() == "2HA") || (va[k].getType() == "HA" ) )
	    {
	      output << "  " << setw(4) 
		     << va[k].getInBond(0).getOutBond(0).getNumber() << "  ";
	      unsigned int c_Index=0;
	      for(unsigned int w = 0; w < va[k].getInBond(0).sizeOutBonds(); w++)
		if(va[k].getInBond(0).getOutBond(w).getType() == "C") 
		  c_Index = w;

	      if((va[k].getType() == "CB") || (va[k].getType() == "2HA") )
		  output << setw(8) << setprecision(3) 
			 << RAD2DEG * icc.getBondAngle( 
				va[k].getInBond(0).getOutBond(c_Index),  
				va[k].getInBond(0), va[k]) << "    1"; 
	      else
		output << setw(8) << setprecision(3) 
		       << RAD2DEG * icc.getBondAngle( 
				  va[k].getInBond(0).getOutBond(c_Index),  
				  va[k].getInBond(0), va[k]) << "   -1"; 
	    }
	  else
	    if (va[k].getInBond(0).getInBond(0).sizeInBonds())
	      output << "  " << setw(4) 
		     << va[k].getInBond(0).getInBond(0).getInBond(0).getNumber() 
		     << "  " << setw(8) << setprecision(3) 
		     << RAD2DEG * icc.getTorsionAngle(
				  va[k].getInBond(0).getInBond(0).getInBond(0),
				  va[k].getInBond(0).getInBond(0), 
				  va[k].getInBond(0), va[k]) << "    0";
	}
      else
	if(isHAtom(va[k].getCode()))  // special case for H atoms
	  {
	    output << "  " << setw(4) 
		   << va[k].getInBond(0).getOutBond(0).getNumber() 
		   << "  " << setw(8) << setprecision(4) << RAD2DEG * 
  	               icc.getBondAngle(va[k].getInBond(0).getOutBond(0), 
			       va[k].getInBond(0), va[k]);
	    
	    if (va[k].getInBond(0).getOutBond(0).sizeOutBonds())
	      output << "  " << setw(4) << 
		va[k].getInBond(0).getOutBond(0).getOutBond(0).getNumber() 
		     << "  " << setw(8) << setprecision(4) 
		     << RAD2DEG * icc.getTorsionAngle(
			    va[k].getInBond(0).getOutBond(0).getOutBond(0),
			    va[k].getInBond(0).getOutBond(0), 
			    va[k].getInBond(0), va[k]) 
		     << "    0";
	  }
      output << "\n";
    }      
  
  output.precision(old_prec);
  output.flags(old_flags);
}







