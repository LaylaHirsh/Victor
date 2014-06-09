/** 
  * 
 * @Class:              LoopModel 
  *      
*/

#ifndef _LOOPMODEL_H_
#define _LOOPMODEL_H_

// Includes:
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector3.h>
#include <matrix3.h>
#include <RamachandranData.h>
#include <LoopTableEntry.h>
#include <Debug.h>
#include <LoopTable.h>
#include <VectorTransformation.h>
#include <Spacer.h>
#include <SeqConstructor.h>
#include <set>
#include <ranking_helper.h>
#include <ranking_helper2.h>
#include <SolvationPotential.h>
#include <RapdfPotential.h>
#include <PhiPsi.h>
#include <ctype.h>

namespace Biopool {
  
// Global constants, typedefs, etc. (to avoid):
 /** @brief  This class implements methods that allow to create a model and also evaluate it.
 * 
 * @Description  Includes methods to calculate the RMS, propensity, etc.
 * */
class LoopModel{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  LoopModel();
  LoopModel(const LoopModel& orig);
  virtual ~LoopModel();

// PREDICATES:
  std::vector<string>& getTableFileName();

  std::vector<int> vdwValues(Spacer& sp, unsigned int index1, 
			unsigned int index2, std::vector<Spacer>& solVec);

  void rankRawScore(Spacer& sp, unsigned int index1, unsigned int index2,
	       std::vector<Spacer>& solVec, unsigned int maxWrite = 9999);

  void doScatterPlot(Spacer& sp, unsigned int index1, unsigned int index2,
	       std::vector<Spacer>& solVec, bool withOxygen = true);

  void refineModel(Spacer& sp, unsigned int index1, unsigned int index2,
	       std::vector<Spacer>& solVec);

  void optimizeModel(Spacer& sp, unsigned int index1, unsigned int index2,
	       std::vector<Spacer>& solVec, bool verbose = true);

  std::vector<double> rankRms2(Spacer& sp, unsigned int index1, unsigned int index2,
			  std::vector<Spacer>& solVec);

  double calcOrigEnergy(const Spacer& sp, unsigned int index1, 
			unsigned int index2);
  double calculateEndRms(const Spacer& sp, unsigned int index1, 
			 unsigned int index2, const Spacer& sp2);
  double calculateRms(const Spacer& sp, unsigned int index1, 
		      unsigned int index2, const Spacer& sp2, 
		      bool output = true, bool withOxygen = false);
  double calculateRms2(const Spacer& sp, unsigned int index1, 
		       unsigned int index2, const Spacer& sp2, bool output, 
		       bool withOxygen = false);
  void setStructure(Spacer& sp, Spacer& sp2, unsigned int index1, 
		    unsigned int index2);

  std::vector<int> consistencyValues(Spacer& sp, unsigned int index1, 
				unsigned int index2, std::vector<Spacer>& solVec);

  double calculatePropensities(Spacer& sp); 
  double calculatePropensities(const Spacer& sp, unsigned int index1, 
			       unsigned int index2, Spacer& sp2); 

  double getENDRMS_WEIGHT(unsigned int len = 0);
  void saveENDRMS_WEIGHT(ostream& output);

  static string getSCWRLConservedSequence(const Spacer& sp, 
				unsigned int index1, unsigned int index2);

  void defineLoopAnchors(Spacer& sp, unsigned int& index1, 
				unsigned int& index2, bool isDeletion);

// MODIFIERS:
  void setInterpolated() { pInter = true; }
  void setNotInterpolated() { pInter = false; }
  void setLimitedVerbose() { pVerbose = 1; } 
  void setVerbose(unsigned int _verb = 10) { pVerbose = _verb; } 
  void setQuiet() { pVerbose = 0; pPlot = false; }
  void setScatterPlot(ostream* _sc) { pPlot = true; pScatter = _sc; }

  void releaseTables(unsigned int index = 0); 
  void setTableFileName(string basename, string ending = ".lt");
  void setTableFileName(std::vector<string>& _tFN);
  std::vector<Spacer> createLoopModel(const AminoAcid& start,
		   const vgVector3<double>& startN, const AminoAcid& end, 
		   const vgVector3<double>& endNAtom, unsigned int indexS, 
		   unsigned int indexE, unsigned int numLoops, 
		   unsigned int numLoops2, std::vector<string> typeVec);
  void clusterLoops(std::vector<Spacer>& vsp);
  
  void copy(const LoopModel& orig);

  int amino_amino_collision(AminoAcid& aa1, AminoAcid& aa2, int pos1, 
			    int pos2, double distance);
  double compactness(Spacer& loop, unsigned int index1, unsigned int index2, 
		     Spacer& proteine);

  double calculateEnergy(const Spacer& sp, unsigned int index1, 
		      unsigned int index2, Spacer& sp2);
  double calculateSecondaryPreference(const Spacer& sp, unsigned int index1, 
		      unsigned int index2, Spacer& sp2);
  double calculatePacking(const Spacer& sp, unsigned int index1, 
		      unsigned int index2, Spacer& sp2);
  double calculateSolvation(const Spacer& sp, unsigned int index1,
		      unsigned int index2, Spacer& sp2);
  double calculateHydrogen(const Spacer& sp, unsigned int index1, 
		      unsigned int index2, Spacer& sp2);
  double calculateCompactness(const Spacer& sp, unsigned int index1, 
		      unsigned int index2, Spacer& sp2);
  double calculateFlankingConformation(const Spacer& sp, unsigned int index1, 
		      unsigned int index2, Spacer& sp2);
  double calculateConsistency(Spacer& sp, unsigned int index1, 
			      unsigned int index2, Spacer& sp2);

  void setENDRMS_WEIGHT(double val, unsigned int len = 0);
  void setAllENDRMS_WEIGHT(double val, unsigned int max = 13);
  void loadENDRMS_WEIGHT(istream& input);

// OPERATORS:
  LoopModel& operator=(const LoopModel& orig);

  static unsigned int OPT_MAX1;
  static unsigned int OPT_MAX2;
  static unsigned int OPT_NUM;
  static double VDW_LIMIT;
  static double SIM_LIMIT;
  static double ENERGY_LIMIT;
  static double ENERGY_WEIGTH;
  static double SECPREF_WEIGTH;
  static double SECPREF_TOL;
  static double PACKING_WEIGTH;

//ATTRIBUTES
  static unsigned int MAX_SPAN;
    static unsigned int MAX_ITER_SOL;
    //static string TABLE_PARAM_FILE; This causes a runtime error on 64 bit OS, use the path variable to avoid it if you are using a 64 bit SO

    static Biopool::SolvationPotential solv;
    static Biopool::RapdfPotential rapdf;
    static Biopool::PhiPsi tor;

protected:
// HELPERS: 
  void convertCoords(AminoAcid start, vgVector3<float> startN, AminoAcid end,
		     vgVector3<float> endN, LoopTableEntry& startEntry, 
		     LoopTableEntry& endEntry, VectorTransformation& vt, 
		     unsigned int nAmino);
  bool ringClosureBase(const LoopTableEntry&, const LoopTableEntry&, 
		       unsigned int, unsigned int, double, 
		       VectorTransformation vt, unsigned int num, 
		       unsigned int depth, 
		       std::vector<vgVector3<float> >& partialSolution);
  bool ringClosure(const LoopTableEntry&, const LoopTableEntry&, unsigned int,
		   unsigned int, double, VectorTransformation vt,
		   std::vector<vgVector3<float> >& partialSolution, 
		   unsigned int currentSelection = 1);
  std::vector<Spacer> calculateLoop(const vgVector3<float>& startN, 
			       unsigned int length, std::vector<string> typeVec);

private:
// HELPERS: 
  void pAddRot(AminoAcid& start, vgVector3<float>& startN, 
	       AminoAcid& end, vgVector3<float>& endN, 
	       const vgMatrix3<float>& rotMat);
  void pAddTrans(AminoAcid& start, vgVector3<float>& startN, 
		 AminoAcid& end, vgVector3<float>& endN, 
		 const vgVector3<float>& trans);
  void loadTables(unsigned int nAmino);
  void pAminoAcidSetup(AminoAcid* aa, string type, double bfac = 90.0);
  void pSetSideChain(AminoAcid& aa);
  double pCalculateLoopRms(Spacer& sp1, Spacer& sp2);

  int loop_loop_vdw(Spacer& sp, unsigned int index1);
  int loop_spacer_vdw(Spacer& loop, unsigned int index1, unsigned int index2, 
		      Spacer& proteine);

  double sCalcEn(AminoAcid& aa, AminoAcid& aa2);
  double sICalcEn(AminoAcid& aa, AminoAcid& aa2);
  double sCalcEnInt(AminoAcid& aa, AminoAcid& aa2);
  double sICalcEnInt(AminoAcid& aa, AminoAcid& aa2);

// ATTRIBUTES:
  bool pInter, pPlot;
  unsigned int pVerbose;
  ostream* pScatter;
  std::vector<LoopTable*> table;
  std::vector<string> tableFileName; 
  std::vector<vgVector3<float> > solution;
  static unsigned int MAX_CHAIN_LENGTH;
  static double BOND_ANGLE_N_TO_CB;
  static double BOND_ANGLE_CB_TO_C;
  static double BOND_LENGTH_CA_TO_CB;

  std::vector<double> ENDRMS_WEIGTH;

};

// ---------------------------------------------------------------------------
//                                  LoopModel
// -----------------x-------------------x-------------------x-----------------


// PREDICATES:
/**
 * @Description  Obtains the names of the files that contains the tables
 * @param   none
 * @return  reference to the table file name(vector<string>&)
 */
inline std::vector<string>& LoopModel::getTableFileName(){
  return tableFileName;
}

// MODIFIERS:
/**
 * @Description  sets the names of the file that contains the table, based on the chain possible lengths, creates a file name for each one.
 * @param   the base name for the file (string) , the ending part for the name(string)
 *  * @return   changes the object internally (void)  
 */
inline void LoopModel::setTableFileName(string basename, string ending){

  tableFileName.clear();
  for (unsigned int i = 0; i < MAX_CHAIN_LENGTH; i++)
    tableFileName.push_back(basename + (uitos(i)).c_str() + ending);
}
/**
 * @Description  sets the names of the files, based on the given names as parameter.
 * @param   a vector containing the file names to load
 * @return   changes the object internally (void)  
 */
inline void LoopModel::setTableFileName(std::vector<string>& _tFN){
  PRECOND( _tFN.size() == MAX_CHAIN_LENGTH, exception);
  tableFileName.clear();
  for (unsigned int i = 0; i < _tFN.size(); i++)
    tableFileName.push_back(_tFN[i]);
}
/**
 * @Description  Adds a rotation, considering a portion of amino acids
 * @param  reference of inicial amino acid and reference for the corresponding coords(AminoAcid&, vgVector3<float>&), 
 *         final reference amino acid and the corresponding coords(AminoAcid&, vgVector3<float>&),
 *         reference for the rotation matrix(const vgMatrix3<float>& )
 * @return   changes the object internally (void)  
 */
inline void LoopModel::pAddRot(AminoAcid& start, vgVector3<float>& startN, 
AminoAcid& end, vgVector3<float>& endN, const vgMatrix3<float>& rotMat){
  for (unsigned int i = 0; i < start.sizeBackbone(); i++)
    start[i].setCoords( convert(rotMat) * start[i].getCoords() );
  startN  =  rotMat * startN;
  for (unsigned int i = 0; i < end.sizeBackbone(); i++)
    end[i].setCoords( convert(rotMat) * end[i].getCoords() );
  endN  =  rotMat * endN;
}
/**
 * @Description  Adds a translation, considering a portion of amino acids
 * @param  reference of inicial amino acid and reference for the corresponding coords(AminoAcid&, vgVector3<float>&), 
 *         final reference amino acid and the corresponding coords(AminoAcid&, vgVector3<float>&),
 *         reference for the translation matrix(const vgMatrix3<float>& )
 * @return   changes the object internally (void)  
 */
inline void LoopModel::pAddTrans(AminoAcid& start, vgVector3<float>& startN, 
AminoAcid& end, vgVector3<float>& endN, const vgVector3<float>& trans){
  for (unsigned int i = 0; i < start.sizeBackbone(); i++)
    start[i].setCoords( start[i].getCoords() - convert(trans));
  for (unsigned int i = 0; i < end.sizeBackbone(); i++)
    end[i].setCoords( end[i].getCoords() - convert(trans));
  startN -= trans;
  endN -= trans;
}
/**
 * @Description  Calculates the root mean square for the loop
 * @param   the two spacers' references, both should have the same length 
 * @return  corresponding value ( double)
 */
inline double LoopModel::pCalculateLoopRms(Spacer& sp1, Spacer& sp2){
  if (sp1.sizeAmino() != sp2.sizeAmino())
    ERROR("Arguments do not match.", exception);
  
  double res = 0.0;
  for (unsigned int i = 0; i < sp1.sizeAmino(); i++)
    res += sqr(sp1.getAmino(i)[N].distance(sp2.getAmino(i)[N]))
           + sqr(sp1.getAmino(i)[CA].distance(sp2.getAmino(i)[CA]))
           + sqr(sp1.getAmino(i)[C].distance(sp2.getAmino(i)[C]));

  return sqrt(res / (3 * sp1.sizeAmino()));
}
/**
 * @Description  Obtains the conserved sequence in the SCWRL format 
 * @param   reference to the spacer that contains all the sequence(Spacer&), the starting and ending index for the conserved part(unsigned int, unsigned int)
 * @return  conserved sequence ( string )
 */
inline string LoopModel::getSCWRLConservedSequence(const Spacer& sp, unsigned int index1, 
unsigned int index2){
  string tmp = "";

  for (unsigned int i = 0; i <= index1; i++)
    tmp += static_cast<char>(tolower(threeLetter2OneLetter(
                       sp.getAmino(i).getType())));

  for (unsigned int i = index1+1; i <= index2; i++)
    tmp += threeLetter2OneLetter(sp.getAmino(i).getType());

  for (unsigned int i = index2+1; i < sp.sizeAmino(); i++)
    tmp += static_cast<char>(tolower(threeLetter2OneLetter(
		       sp.getAmino(i).getType())));

  return tmp;
}

} // namespace

#endif //_LOOPMODEL_H_



