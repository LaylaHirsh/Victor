/** 
  * 
 * @Class:           LoopModel   
  * 
 * @Description: This class implements methods that allow to create a model and also evaluate it. 
 *  Includes methods to calculate the RMS, propensity, etc.
  *      
*/

// Includes:
#include <LoopModel.h>
#include <LoopTableEntry.h>
#include <IntCoordConverter.h>
#include <String2Number.h>
#include <queue>
 
#include <Debug.h>
#include <limits.h>
#include <float.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


unsigned int LoopModel::MAX_CHAIN_LENGTH = 64;
unsigned int LoopModel::OPT_MAX1 = 5;
unsigned int LoopModel::OPT_MAX2 = 5;
unsigned int LoopModel::OPT_NUM = 5;
double LoopModel::VDW_LIMIT = 100.0;
double LoopModel::SIM_LIMIT = 0.1;
double LoopModel::ENERGY_LIMIT = 8.0;
double LoopModel::ENERGY_WEIGTH = 1.0;
double LoopModel::SECPREF_WEIGTH = 0.05;
double LoopModel::SECPREF_TOL = 30.0;
double LoopModel::PACKING_WEIGTH = 1.0;

double LoopModel::BOND_ANGLE_N_TO_CB = 110.5;
double LoopModel::BOND_ANGLE_CB_TO_C = 110.1;
double LoopModel::BOND_LENGTH_CA_TO_CB = 1.54;

unsigned int LoopModel::MAX_SPAN = 15;
unsigned int LoopModel::MAX_ITER_SOL = 500;
//NOTE:
//This causes a runtime error on 64 bit OS, use the path variable to avoid it if you are using a 64 bit SO
//string LoopModel::TABLE_PARAM_FILE = "data/aa";

const double WEIGHT_SEC = 10.0;
const double L_HELIX = 0.414;
const double L_STRAND = 0.42;
const double L_COIL = 0.13;
SolvationPotential LoopModel::solv;
RapdfPotential LoopModel::rapdf;
PhiPsi LoopModel::tor;


/**
 * @Description Calculates distance between two points
 * @param  two point's coords (vgVector3<double>,vgVector3<double>)
 * @return  the corresponding value(static double)
 */
static double sDistance( vgVector3<double> xv,  vgVector3<double> yv){
  return  sqrt(sqr(xv.x-yv.x) + sqr(xv.y-yv.y) + sqr(xv.z-yv.z));
}


// CONSTRUCTORS/DESTRUCTOR:
/**
 * @Description  Basic constructor
 */
LoopModel::LoopModel() : pInter(false), pPlot(false), pVerbose(0), pScatter(&cout), table(), solution(){ 
    char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string path = "data/aa";//comment this if you are using a no 64bits OC
    string tableFileR = getenv("VICTOR_ROOT");
  if (tableFileR.length() < 3)
    ERROR("Environment variable VICTOR_ROOT was not found.", exception);

  //string tableFile = tableFileR + TABLE_PARAM_FILE;
  string tableFile = tableFileR + path;
  setTableFileName(tableFile);
 
  ENDRMS_WEIGTH.push_back(125);
}


/**
 * @Description  Basic destructor
 */
LoopModel::~LoopModel(){
  PRINT_NAME;
  while (table.size() > 0)    {
      delete table[table.size()-1];
      table.pop_back();
    }

  solution.clear();
}
/**
 * @Description   releases memory occupied by loop tables bigger than the index 
 * @param   maximum size to consider(unsigned int)
 * @return   changes the object internally (void)  
 */
void LoopModel::releaseTables(unsigned int index){
  while (table.size() > index)    {
      delete table[table.size()-1];
      table.pop_back();
    }
  solution.clear();
}


// PREDICATES:
/**
 * @Description Sets the vdw_forces
 * @param   reference for the spacer(Spacer& sp),the wto index to consider
 *          (unsigned int,unsigned int) vector that contains the solutions(std::vector<Spacer>&)
 * @return   a vector containing the corresponding values (std::vector<int>)
 */
std::vector<int>LoopModel::vdwValues(Spacer& sp, unsigned int index1, unsigned int index2,
std::vector<Spacer>& solVec){
  std::vector<int> ret_vector; 
  int count = 0;          

  for (unsigned int loop = 0; loop < solVec.size(); loop++)    {
      count = loop_loop_vdw(solVec[loop], index1 + 1);
      count += loop_spacer_vdw(solVec[loop], index1 + 1, index2 + 1, sp);
      
      ret_vector.push_back(count); 
    }
 return ret_vector; 
}

 
/**
 * @Description  Calculates the consistency values for each generated solution. 
 * @param   Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solutions/loops( vector<Spacer>&)
 * @return   the consistency values(std::vector<int>).
 */
std::vector<int>LoopModel::consistencyValues(Spacer& sp, unsigned int index1, 
unsigned int index2, vector<Spacer>& solVec) {
  std::vector<int> ret_vector; 
  for (unsigned int loop = 0; loop < solVec.size(); loop++) 
      ret_vector.push_back(calculateConsistency(sp, index1, index2, 
						solVec[loop]) 
	    + 100 * loop_loop_vdw(solVec[loop], index1 + 1)
	    + 100 * loop_spacer_vdw(solVec[loop], index1 + 1, index2 + 1, sp)
			   );
  return ret_vector;
}

/**
 * @Description  Calculates the consistency between two proteins 
 * @param   Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   Spacer reference that contains the protein (Spacer&)
 * @return   the consistency value(double)
 */
double LoopModel::calculateConsistency(Spacer& sp, unsigned int index1, 
unsigned int index2, Spacer& sp2) {
  double tmp, count = 0;
  IntCoordConverter icc;

  tmp = calculateRms(sp, index1, index2, sp2, false);
  if (tmp > 1.0)  {
      if (pVerbose > 2)
	  cout << "ENDRMS check failed. ( " 
	       << tmp << " ) \t\t\t "
	       << calculateRms2(sp, index1, index2, sp2, false) << "\n";

      if (tmp < 1.5)
	  count += 2000;
      else
	  count += 3000;
  }
  
  tmp = calculatePropensities(sp, index1, index2, sp2);

  if (tmp > 0.5)  {
      if (pVerbose > 2)
	  cout << "PROPENSITY check failed. ( " 
	       << tmp << " ) \t\t\t"
	       << calculateRms2(sp, index1, index2, sp2, false) <<"\n";
      count += 4000;
  }
      
  // *** consistency check for the loop:

  tmp = sDistance(sp2.getAmino( sp2.sizeAmino()-1)[N].getCoords(), 
		  sp.getAmino(index2)[CA].getCoords());
  if ( (tmp < 0.9) || (tmp > 1.9) )  {
      if (pVerbose > 2)
	  cout << "loop closure length check failed. ( " << tmp << " ) \t "
	       << calculateRms2(sp, index1, index2, sp2, false) << "\n";
      count += 8000;
  }
  
  tmp = RAD2DEG * icc.getBondAngle(sp2.getAmino(sp2.sizeAmino()-2)[C], 
	   sp2.getAmino(sp2.sizeAmino()-1)[N], sp.getAmino(index2)[CA]);
  if (fabs(tmp - 121.5) > 20.0)  {
      if (pVerbose > 2)
	  cout << "loop closure angle check failed. ( " << tmp << " ) \t "
	       << calculateRms2(sp, index1, index2, sp2, false) << "\n";
      count += 8000;
  }
     
  tmp = RAD2DEG * icc.getTorsionAngle(sp2.getAmino(
         sp2.sizeAmino()-2)[CA], sp2.getAmino(sp2.sizeAmino()-2)[C], 
	 sp2.getAmino(sp2.sizeAmino()-1)[N], sp.getAmino(index2)[CA]);
  if ( (fabs(tmp - 180.0 ) > 40.0) && (fabs(tmp + 180.0 ) > 40.0) )  {
      if (pVerbose > 2)
	  cout << "loop closure torsion check failed. ( " << tmp << " ) \t "
	       << calculateRms2(sp, index1, index2, sp2, false) << "\n";
      count += 8000;
  }

  // proline filter:
  
  if (sp2.getAmino(0).getType() == "PRO")  {
      double tPhi = RAD2DEG * icc.getTorsionAngle(
	  const_cast<Atom&>(sp.getAmino(index1)[C]), sp2.getAmino(0)[N], 
	  sp2.getAmino(0)[CA], sp2.getAmino(0)[C]);
      double tPsi = sp2.getAmino(0).getPsi();

      if ( (tPhi < -140.0) || (tPhi > -10.0) 
	   || ( (tPsi < -90.0) && (tPsi > -140.0) ) )      {
	  count += 16000;
	  if (pVerbose > 2)
	      cout << "PROLINE torsion angle check failed. ( " 
		   << tPhi << " / " << tPsi << " ) \t " 
		   << calculateRms2(sp, index1, index2, sp2, false) <<"\n";
      }
  }
  
  for (unsigned int i = 1; i < sp2.sizeAmino()-1; i++)
  if (sp2[i].getType() == "PRO")   {
      double tPhi = sp2.getAmino(i).getPhi();
      double tPsi = sp2.getAmino(i).getPsi(); 
      
      if ( (tPhi < -140.0) || (tPhi > -10.0)
	   ||  ( (tPsi < -90.0) && (tPsi > -140.0) ) )      {
	  count += 16000;
	  if (pVerbose > 2)
	      cout << "PROLINE torsion angle check failed. ( " 
		   << tPhi << " / " << tPsi << " ) \t " 
		   << calculateRms2(sp, index1, index2, sp2, false) <<"\n";
	  break;
      }
  }  
  
  unsigned int tNum = sp2.sizeAmino()-1;
  if (sp2.getAmino(tNum).getType() == "PRO")  {
      double tPhi = sp2.getAmino(tNum).getPhi();

      double tPsi =  RAD2DEG * icc.getTorsionAngle( sp2.getAmino(tNum)[N], 
	 sp2.getAmino(tNum)[CA], sp2.getAmino(tNum)[C], 
	 const_cast<Atom&>(sp.getAmino(index2+1)[N]));

      if ( (tPhi < -140.0) || (tPhi > -10.0) 
	   || ( (tPsi < -90.0) && (tPsi > -140.0) ) )  {
	  count += 16000;
	  if (pVerbose > 2)
	      cout << "PROLINE torsion angle check failed. ( " 
		   << tPhi << " / " << tPsi << " ) \t " 
		   << calculateRms2(sp, index1, index2, sp2, false) <<"\n";
      }
  }
  
 return count;
}
 
/**
 * @Description  Calculates a different random value inside a range, based in a original value
 * @param   original value(int), limit(int)
 * @return   corresponding value(double)
 */
double sRandom(int c, int m){
  double MAX_DISP = 0.1 + 0.5 * (m/2 - fabs( c - (m/2)));

  return MAX_DISP * rand()/RAND_MAX - (MAX_DISP/2);

}

/**
 * @Description  Refines the possible solutions,refine end rms:
 * @param   Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solutions/loops( vector<Spacer>&)
 * @return   changes the object internally (void)  
 */
void LoopModel::refineModel(Spacer& sp, unsigned int index1, unsigned int index2,
std::vector<Spacer>& solVec){
  
  unsigned int offset = index2 - index1 - 1;
  unsigned int max = index2-index1+1;
  
  for (unsigned int svOffset = 0; svOffset < solVec.size(); svOffset++)    {
      vgVector3<double> tmp = ( (solVec[svOffset].getAmino(offset)[N].getCoords() 
				 - sp.getAmino(index2)[N].getCoords()) );

      solVec[svOffset].getAmino(offset)[N].setCoords(
	     solVec[svOffset].getAmino(offset)[N].getCoords() - tmp/ 2 );

      for (unsigned int i = 1; i < max/2; i++)
	{
	  solVec[svOffset].getAmino(offset - i)[N].setCoords(
		solVec[svOffset].getAmino(offset - i)[N].getCoords() 
		- tmp * ((max-i-1)/max));
	  solVec[svOffset].getAmino(offset - i)[CA].setCoords(
		 solVec[svOffset].getAmino(offset - i)[CA].getCoords() 
		 - tmp * ((max-i-1)/max));
	  solVec[svOffset].getAmino(offset - i)[C].setCoords(
		solVec[svOffset].getAmino(offset - i)[C].getCoords() 
		- tmp * ((max-i-1)/max));
	}
    }
}
/**
 * @Description  Optimizes the models, can be a verbose process
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solutions/loops( vector<Spacer>&),verbose flag (bool)
 * @return   changes the object internally (void)  
 */
void LoopModel::optimizeModel(Spacer& sp, unsigned int index1, unsigned int index2,
std::vector<Spacer>& solVec, bool verbose){
  if (solVec.size() == 0)
    return;

  if (verbose)
    {
      cout << "-----------------------------------------\n";
      cout << "Optimizing:\n";
    }

  for (unsigned svOffset = 0; svOffset < (solVec.size() > OPT_NUM ? OPT_NUM : 
					  solVec.size() ); svOffset++)    {
      if (verbose)	{
	  cout << "-----------------------------------------\n";
	  cout << ">>>>> optimizing solution # " << svOffset << "\n";
	}

      double minEn = calculateEnergy(sp, index1, index2, solVec[0]);
      double minScore = getENDRMS_WEIGHT() * calculateRms(sp, index1, index2, 
		         solVec[svOffset], false) + ENERGY_WEIGTH * minEn;

      double bestScore = minScore;
      int max = index2-index1+1;

      if (verbose){
	  cout << "en= " << minEn << "\t";
	  calculateRms(sp, index1, index2, solVec[svOffset]);

	  cout << "-----------------------------------------\n";
	  cout << "Attempting local optimization:\n";
	}

      Spacer bestSp = solVec[svOffset];

      for (unsigned int i = 0; i < OPT_MAX1; i++)
	for (unsigned int j = 0; j < OPT_MAX2; j++) {
	    Spacer tmpSp = solVec[svOffset];
	    minScore = getENDRMS_WEIGHT() * calculateRms(sp, index1, index2, 
		         solVec[svOffset], false) 
	     + ENERGY_WEIGTH *  calculateEnergy(sp, index1, index2, 
						solVec[svOffset]);

	    // attempt local modifications
	    int curr = rand() % (max-1);
	
	    vgVector3<double> disp(0.0, 0.0, 0.0);
	
	    disp[0] += sRandom(curr, max-1);
	    disp[1] += sRandom(curr, max-1);
	    disp[2] += sRandom(curr, max-1);
	
	    tmpSp.getAmino(curr)[N].setCoords(
		tmpSp.getAmino(curr)[N].getCoords() + disp);
	    tmpSp.getAmino(curr)[CA].setCoords(
		tmpSp.getAmino(curr)[CA].getCoords() + disp);
	    tmpSp.getAmino(curr)[C].setCoords(
		tmpSp.getAmino(curr)[C].getCoords() + disp);

	    if (curr > 1) {
		tmpSp.getAmino(curr-1)[N].setCoords(
		  tmpSp.getAmino(curr-1)[N].getCoords() + disp/2);
		tmpSp.getAmino(curr-1)[CA].setCoords(
		  tmpSp.getAmino(curr-1)[CA].getCoords() + disp/2);
		tmpSp.getAmino(curr-1)[C].setCoords(
		  tmpSp.getAmino(curr-1)[C].getCoords() + disp/2);
	      }

	    if (curr < max-2) {
		tmpSp.getAmino(curr+1)[N].setCoords(
		  tmpSp.getAmino(curr+1)[N].getCoords() + disp/2);
		tmpSp.getAmino(curr+1)[CA].setCoords(
		  tmpSp.getAmino(curr+1)[CA].getCoords() + disp/2);
		tmpSp.getAmino(curr+1)[C].setCoords(
		  tmpSp.getAmino(curr+1)[C].getCoords() + disp/2);
	      }

	    double actEn = calculateEnergy(sp, index1, index2, tmpSp); 
	    double actScore = getENDRMS_WEIGHT() * calculateRms(sp, index1, index2, 
		         tmpSp, false) 
	      + ENERGY_WEIGTH * actEn;
	    if (actScore < minScore) {
		if (verbose)  {
		    cout << i << " " << j << " new minimum= " << actEn << "\t";
		    calculateRms(sp, index1, index2, tmpSp);
		  }

		minScore = actScore;

		if (actScore < bestScore) {
		    bestScore = actScore;
		    bestSp = tmpSp;
		  }
	      }
	  }
      solVec[svOffset] = bestSp;
      if (verbose){
	  cout << "-----------------------------------------\n";
	  cout << "selected = " << calculateEnergy(sp, index1, index2, bestSp) 
	       << "\t";
	  calculateRms(sp, index1, index2, bestSp);
	}
    }
  if (verbose) {
      cout << "-----------------------------------------\n";
      unsigned int size1 = sp.sizeAmino();
      double en = 0.0;
            
      for (unsigned int i = 0; i < index1; i++)
	for (unsigned int j = index1; j < index2; j++)
	  en += rapdf.calculateEnergy(const_cast<AminoAcid&>(sp.getAmino(i)), 
				      sp.getAmino(j));
      
      for (unsigned int i = index2+1; i < size1; i++)
	for (unsigned int j = index1; j < index2; j++)
	  en += rapdf.calculateEnergy(const_cast<AminoAcid&>(sp.getAmino(i)), 
				      sp.getAmino(j));
      
      for (unsigned int i = 0; i < size1; i++)
	for (unsigned int j = index1; j < index2; j++)
	  en += rapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
      
      cout << "Original energy:\t" << en << "\n";
    }
}

/**
 * @Description  Calculates the energy 
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution(const Spacer&)
 * @return   corresponding energy value(double)
 */
double LoopModel::calculateEnergy(const Spacer& sp, unsigned int index1, 
unsigned int index2, Spacer& sp2){
  long double en = 0.0;
  unsigned int size1 = sp.sizeAmino();
  unsigned int size2 = sp2.sizeAmino();

  for (unsigned int i = 0; i < index1; i++)
    for (unsigned int j = 0; j < size2; j++)
      en += rapdf.calculateEnergy(const_cast<AminoAcid&>(sp.getAmino(i)), 
				  sp2.getAmino(j));
  
  for (unsigned int i = index2+1; i < size1; i++)
    for (unsigned int j = 0; j < size2; j++)
      en += rapdf.calculateEnergy(const_cast<AminoAcid&>(sp.getAmino(i)), 
				  sp2.getAmino(j));
  
  for (unsigned int i = 0; i < size2; i++)
    for (unsigned int j = i; j < size2; j++)
      en += rapdf.calculateEnergy(sp2.getAmino(i), sp2.getAmino(j));
  
  return en;

}
/**
 * @Description  Calculates the secondary Preference
 * @param    Spacer reference that contains the protein (const Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution ( Spacer&)
 * @return   corresponding value(double)
 */
double LoopModel::calculateSecondaryPreference(const Spacer& sp, unsigned int index1, unsigned int index2, Spacer& sp2){


  double res = sqr(sp2.getAmino(1).getPhi() 
		   - const_cast<Spacer&>(sp).getAmino(index1).getPhi()) 
    + sqr(sp2.getAmino(0).getPsi() 
	           - const_cast<Spacer&>(sp).getAmino(index1).getPsi())
    + sqr(sp2.getAmino(sp2.sizeAmino()-1).getPhi() 
		   - const_cast<Spacer&>(sp).getAmino(index2+1).getPhi()) 
    + sqr(sp2.getAmino(sp2.sizeAmino()-2).getPsi() 
	           - const_cast<Spacer&>(sp).getAmino(index2+1).getPsi());

  if (sqrt(res/2)  > SECPREF_TOL)
    return sqrt(res/2) - SECPREF_TOL;

  return 0;
}

/**
 * @Description  calculates the estimation of the accessible surface area of a protein molecule
 * @param   two amino acids to consider
 * @return   corresponding value(doube)
 */
double pCalculateOoi(const AminoAcid& aa1, const AminoAcid& aa2){
  return (const_cast<Atom&>(aa1[CA]).distance(
	     const_cast<Atom&>(aa2[CA])) > 14.0 ? 0.0 : 
	   const_cast<Atom&>(aa1[CA]).distance(
	      const_cast<Atom&>(aa2[CA])) > 8.0 ? -0.5 : -1.0); 
}

/**
 * @Description  Calculates the packing
 * @param    Spacer reference that contains the protein (const Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&)
 * @return   
 */
double LoopModel::calculatePacking(const Spacer& sp, unsigned int index1, unsigned int index2, Spacer& sp2){

  double res = 0.0;
  for (unsigned int i = 0; i < index2 - index1-1; i++){
      for (unsigned int j = 0; j < index1-1; j++)
	res += pCalculateOoi(sp.getAmino(j), sp2.getAmino(i));
      for (unsigned int j = index2+1; j < sp.sizeAmino(); j++)
	res += pCalculateOoi(sp.getAmino(j), sp2.getAmino(i));
    }

  return res / (index2 - index1 - 1);
}

/**
 * @Description calculates the solvation value  
 * @param    Spacer reference that contains the protein (const Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&)
 * @return   corresponding value(double)
 */
double LoopModel::calculateSolvation(const Spacer& sp, unsigned int index1, unsigned int index2, Spacer& sp2){
  long double res = 0.0;

  for (unsigned int i = 0; i < sp2.sizeAmino(); i++)
    res += solv.calculateSolvation(const_cast<AminoAcid&>(sp2.getAmino(i)), 
				   const_cast<Spacer&>(sp));

  return res;
}

/**
 * @Description  Calculates the hidrogen values
 * @param   the two amino acids to consider
 * @return   corresponding value(double)
 */
double sCalcHydrogen(AminoAcid& aa1, AminoAcid& aa2){
  if (!aa2.isMember(O))
    return 0.0;
  
  double dist = aa1[N].distance(aa2[O]);
  double dist2 = aa1[N].distance(aa2[C]); 
  double dist3 = aa1[CA].distance(aa2[O]); 
	  
  if ((dist <= 4.0) && (dist >= 2.0) && (dist < dist2)
      && (dist < dist3))
      return -3.0;
  else
    return 0.0;
}

/**
 * @Description  calculates the hidrogen value for the protein
 * @param    Spacer reference that contains the protein (const Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&)
 * @return   corresponding value
 */
double LoopModel::calculateHydrogen(const Spacer& sp, unsigned int index1, unsigned int index2, Spacer& sp2){
  long double res = 0.0;

  for (unsigned int i = 0; i < sp2.sizeAmino(); i++)
    for (unsigned int j = 0; j < sp2.sizeAmino(); j++)
      if (abs((long)(i - j)) >= 2)
	res += sCalcHydrogen(const_cast<AminoAcid&>(sp2.getAmino(i)), 
		      const_cast<AminoAcid&>(sp2.getAmino(j)));

  for (unsigned int i = 0; i <= index1-2; i++)
    for (unsigned int j = 0; j < sp2.sizeAmino(); j++)
	res += sCalcHydrogen(const_cast<AminoAcid&>(sp.getAmino(i)), 
		      const_cast<AminoAcid&>(sp2.getAmino(j)));

  for (unsigned int j = 1; j < sp2.sizeAmino(); j++)
    res += sCalcHydrogen(const_cast<AminoAcid&>(sp.getAmino(index1-1)), 
		  const_cast<AminoAcid&>(sp2.getAmino(j)));

  for (unsigned int j = 0; j < sp2.sizeAmino()-1; j++)
    res += sCalcHydrogen(const_cast<AminoAcid&>(sp2.getAmino(j)), 
		  const_cast<AminoAcid&>(sp.getAmino(index2+1)));

  for (unsigned int i = index2+2; i < sp.sizeAmino(); i++)
    for (unsigned int j = 0; j < sp2.sizeAmino(); j++)
	res += sCalcHydrogen(const_cast<AminoAcid&>(sp2.getAmino(j)), 
		      const_cast<AminoAcid&>(sp.getAmino(i)));

  return res;
}

/**
 * @Description  Calculates the compactness besides the protein and the possible solution
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&)
 * @return   corresponding value(double)
 */
double LoopModel::calculateCompactness(const Spacer& sp, unsigned int index1, unsigned int index2, Spacer& sp2){
  long double res = 10000.0;

  for (unsigned int i = 0; i < index2 - index1-1; i++)    {
	unsigned int limit = (index1 > 3) ? 3 : index1-1;
      for (unsigned int j = 0; j < limit; j++)
	  if (const_cast<AminoAcid&>(sp.getAmino(j))[CA].distance(
	      sp2.getAmino(i)[CA]) < res)
	      res = const_cast<AminoAcid&>(sp.getAmino(j))[CA].distance(
		  sp2.getAmino(i)[CA]);
      for (unsigned int j = index2+3; j < sp.sizeAmino(); j++)
	  if (const_cast<AminoAcid&>(sp.getAmino(j))[CA].distance(
	      sp2.getAmino(i)[CA]) < res)
	      res = const_cast<AminoAcid&>(sp.getAmino(j))[CA].distance(
		  sp2.getAmino(i)[CA]);
    }

  return res;
}


/**
 * @Description Calculates the Flanking conformation  
 * @param    Spacer reference that contains the protein (const Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&)
 * @return   corresponding value(double)
 */
double LoopModel::calculateFlankingConformation(const Spacer& sp, 
unsigned int index1, unsigned int index2, Spacer& sp2){
    double res = 1.0;

    // get flanks:

    unsigned int flankN = COIL;
    unsigned int limit = (index1 > 5) ? 5 : index1-1;
    for (unsigned int i = 1; i < limit; i++)
	if ((flankN != HELIX) && (flankN != STRAND))
	    flankN = const_cast<AminoAcid&>(sp.getAmino(index1-i)).getState();

    unsigned int flankC = COIL;
    limit = (sp.size()-index2 > 5) ? 5 : sp.size()-index2;
    for (unsigned int i = 1; i < 5; i++)
	if ((flankC != HELIX) && (flankC != STRAND))
	    flankC = const_cast<AminoAcid&>(sp.getAmino(index2+i)).getState();


    unsigned int L1;
    unsigned int Ln_1 = COIL;
    unsigned int Ln = COIL;

    // get Ramachandran area:

    IntCoordConverter icc;

    double tPhi = RAD2DEG * icc.getTorsionAngle(
		     const_cast<Atom&>(sp.getAmino(index1)[C]), 
		     sp2.getAmino(0)[N], 
		     sp2.getAmino(0)[CA], sp2.getAmino(0)[C]);

    double tPsi = sp2.getAmino(0).getPsi(); 

    if (((tPhi >= -140) && (tPhi <= -40)) && ((tPsi >= -80) && (tPsi <= 40)))
	L1 = HELIX;
    else if ((((tPhi >= -180) && (tPhi <= -40)) 
	     && ( ((tPsi >= 60) && (tPsi <= 180)))) || (tPsi <= -140))
	L1 = STRAND;
    else 
	L1 = COIL;

    if (sp2.sizeAmino() > 3)    {
	double tPhi = sp2.getAmino(sp2.sizeAmino()-3).getPhi(); 
	double tPsi = sp2.getAmino(sp2.sizeAmino()-3).getPsi(); 

	if (((tPhi >= -140) && (tPhi <= -40)) && ((tPsi >= -80) 
						  && (tPsi <= 40)))
	    Ln_1 = HELIX;
	else if (((tPhi >= -180) && (tPhi <= -40))
	     && ( (((tPsi >= 60) && (tPsi <= 180))) || (tPsi <= -140)))
	    Ln_1 = STRAND;
	else 
	    Ln_1 = COIL;
    }

    if (sp2.sizeAmino() > 2)    {
	double tPhi = sp2.getAmino(sp2.sizeAmino()-2).getPhi(); 
	double tPsi = sp2.getAmino(sp2.sizeAmino()-2).getPsi(); 

	if (((tPhi >= -140) && (tPhi <= -40)) && ((tPsi >= -80) 
						  && (tPsi <= 40)))
	    Ln = HELIX;
	else if (((tPhi >= -180) && (tPhi <= -40)) 
	     && (( ((tPsi >= 60) && (tPsi <= 180))) || (tPsi <= -140)))
	    Ln = STRAND;
	else 
	    Ln = COIL;
    }

    // now do the prop's:

    if (flankN == HELIX)    {
	if (L1 == HELIX)
	    res *= 0.38 / L_HELIX;
	else if (L1 == STRAND)
	    res *= 0.294 / L_STRAND;
	else
	    res *= 0.138 / L_COIL;
    }
    else if (flankN == STRAND)    {
	if (L1 == HELIX)
	    res *= 0.484 / L_HELIX;
	else if (L1 == STRAND)
	    res *= 0.362 / L_STRAND;
	else
	    res *= 0.155 / L_COIL;
    }

    if (sp2.sizeAmino() > 3)    {
	if (flankC == HELIX)	{
	    if (Ln_1 == HELIX)
		res *= 0.388 / L_HELIX;
	    else if (Ln_1 == STRAND)
		res *= 0.523 / L_STRAND;
	    else
		res *= 0.089 / L_COIL;
	    if (Ln == HELIX)
		res *= 0.268 / L_HELIX;
	    else if (Ln == STRAND)
		res *= 0.633 / L_STRAND;
	    else
		res *= 0.098 / L_COIL;
	}
	else if (flankC == STRAND)	{
	    if (Ln_1 == HELIX)
		res *= 0.388 / L_HELIX;
	    else if (Ln_1 == STRAND)
		res *= 0.523 / L_STRAND;
	    else
		res *= 0.089 / L_COIL;
	    if (Ln == HELIX)
		res *= 0.268 / L_HELIX;
	    else if (Ln == STRAND)
		res *= 0.633 / L_STRAND;
	    else
		res *= 0.098 / L_COIL;
	}
    }

    return res;
}


/**
 * @Description  Ranks the possible solutions
  * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solutions/loops( vector<Spacer>&), the maximum quantity of solutions to consider(unsigned int))
 * @return   changes the object internally (void)  
 */
void LoopModel::rankRawScore(Spacer& sp, unsigned int index1, unsigned int index2,
std::vector<Spacer>& solVec, unsigned int maxWrite){
  std::vector<Spacer> tmpSolVec;
  priority_queue<solutionQueueElem> solutionQueue;

  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
  cout << "A total of " << setw(5) << solVec.size() 
       << " solutions were generated.\n";
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";


  std::vector<double> tmpScore;
  for (unsigned int i = 0; i < solVec.size(); i++)    {
      solutionQueueElem sqe;
      sqe.dev = calculateConsistency(sp, index1, index2, solVec[i])
                 + 100 * loop_loop_vdw(solVec[i], index1 + 1)
	         + 100 * loop_spacer_vdw(solVec[i], index1 + 1, 
					  index2 + 1, sp);

      sqe.index1 = i;
      solutionQueue.push(sqe);
    }
  unsigned int tmpIndex = 0;
  unsigned int numNotFiltered = 0;

  while (solutionQueue.size() > 0)    {
	if ( ( (tmpIndex > 50) && (solutionQueue.top().dev > 1000) )
	     || (tmpIndex > 1000) )
	  break;

	tmpIndex++;
	tmpSolVec.push_back(solVec[solutionQueue.top().index1]);
	tmpScore.push_back(solutionQueue.top().dev);
	if (solutionQueue.top().dev < 1000)
	  numNotFiltered++;
	
	solutionQueue.pop();
    }

  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
  cout << "A total of " << setw(5) << numNotFiltered
       << " solutions survived the filtering procedure.\n";
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";

  solVec.clear();
  solVec = tmpSolVec;

  tmpSolVec.clear();
  while (solutionQueue.size() > 0)
      solutionQueue.pop();

  for (unsigned int i = 0; i < solVec.size(); i++)    {
      solutionQueueElem sqe;
      sqe.dev = getENDRMS_WEIGHT() * calculateRms(sp, index1, index2, 
						  solVec[i], false)
	  + ENERGY_WEIGTH * calculateEnergy(sp, index1, index2, solVec[i])
	  + tmpScore[i];

      sqe.index1 = i;
      solutionQueue.push(sqe);
    }
  
  unsigned int oldPrec = cout.precision();
  ios::fmtflags oldFlags = cout.flags();
  cout.setf(ios::fixed, ios::floatfield);

  tmpIndex = 0;
  if (pVerbose > 1)
    cout << "Raw-score: sum \t endrms \t energy \t  prop \n"
	 << "               \t " << setw(6) <<"-"<< setprecision(2) 
 	<<"-"<< getENDRMS_WEIGHT() <<"-"<< "\t\t " 
	 << setw(6) << setprecision(2) << ENERGY_WEIGTH << "\n";

  unsigned int counter = 0;
  while (solutionQueue.size() > 0)    {
      if (pVerbose > 1)      {    
	  double tmpEndRms = calculateEndRms(sp, index1, index2, 
			      solVec[solutionQueue.top().index1]);
	  double tmpEnergy = calculateEnergy(sp, index1, 
 			      index2, solVec[solutionQueue.top().index1]);
	  double tmpProp = calculatePropensities(sp, index1, index2,
		 	      solVec[solutionQueue.top().index1]);
	  cout << tmpIndex 
	       << " \t " << setw(7) << setprecision(3) 
	       << solutionQueue.top().dev 
	       << " \t " << setw(6) << setprecision(3) << tmpEndRms
	       << " \t " << setw(6) << setprecision(3) << tmpEnergy
	       << " \t " << setw(6) << setprecision(3) << tmpProp
	       << "\n";
      }

     tmpIndex++;
      
      tmpSolVec.push_back(solVec[solutionQueue.top().index1]);
      solutionQueue.pop();

      counter++;
      if (counter >= maxWrite)
	break;
    }

  cout.precision(oldPrec);
  cout.flags(oldFlags);

  solVec.clear();
  solVec = tmpSolVec;

  tmpSolVec.clear();
}

/**
 * @Description  Creates the Scatter Plot
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solutions/loops( vector<Spacer>&),verbose for oxygen (bool)
 * @return   changes the object internally (void)  
 */
void LoopModel::doScatterPlot(Spacer& sp, unsigned int index1, unsigned int index2,
vector<Spacer>& solVec, bool withOxygen){
  if (!pPlot)
      return;

  cout << "WARNING: Scatter plot is incompatible with old-style analysis! "
       << " Updated for SVM classification - ST; 31/12/02" << endl;

  *pScatter << "LOOP \t" << index2-index1 << "\t" 
	   << solVec.size() << "\n";

  for (unsigned int i = 0; i < solVec.size(); i++)  {
      double tmpEndRms = calculateRms(sp, index1, index2, 
			     solVec[i], false);
      double tmpEnergy = calculateEnergy(sp, index1, 
 			     index2, solVec[i]);
      double tmpCons = calculateConsistency(sp, index1, index2, solVec[i]);

	
      double tmpSolv = calculateSolvation(sp, index1, index2, 
					  solVec[i]);
      double tmpSecPref = calculateSecondaryPreference(sp, index1, index2, 
					  solVec[i]);
      double tmpPack = calculatePacking(sp, index1, index2, 
					  solVec[i]);
      double tmpHydro = calculateHydrogen(sp, index1, index2,
					  solVec[i]);
      double tmpComp = calculateCompactness(sp, index1, index2,
					  solVec[i]);
      double tmpFlank = calculateFlankingConformation(sp, index1, index2,
					  solVec[i]);
      double tmpProp = calculatePropensities(sp, index1, index2,
					  solVec[i]);

      *pScatter << setw(5) << calculateRms2(sp, index1, index2, 
				     solVec[i], false, withOxygen)
	       << " \t 1:" << (getENDRMS_WEIGHT() * tmpEndRms
					+ ENERGY_WEIGTH * tmpEnergy)
	       << " \t 2:" << tmpEndRms
	       << " \t 3:" << tmpEnergy
	       << " \t 4:" << tmpSolv
	       << " \t 5:" << tmpSecPref
	       << " \t 6:" << tmpPack
	       << " \t 7:" << tmpHydro
	       << " \t 8:" << tmpComp
	       << " \t 9:" << tmpFlank
	       << " \t 10:" << tmpProp
	       << " \t 11:" << tmpCons
	       << "\n";
  }
  
}

/**
 * @Description  returns the end-RMS values for the solutions in solVec as an array
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solutions/loops( vector<Spacer>&) 
 * @return  corresponding values of final RMS(vector<double>)) 
 */
std::vector<double> LoopModel::rankRms2(Spacer& sp, unsigned int index1, unsigned int index2,
		    std::vector<Spacer>& solVec){
  std::vector<double> tmpSolVec;

  for (unsigned int i = 0; i < solVec.size(); i++)    {
      tmpSolVec.push_back(calculateRms(sp, index1, index2, solVec[i], false));
    }
  return tmpSolVec;
}

/**
 * @Description  calculates the RMS as a global RMS
 * @param    Spacer reference that contains the protein (const Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution ( const Spacer&), output flag (bool), oxygen flag (bool)
 * @return   corresponding value(double)
 */
double LoopModel::calculateRms2(const Spacer& sp, unsigned int index1, 
unsigned int index2, const Spacer& sp2, bool output, bool withOxygen){
  PRECOND( index2 - index1 == sp2.size(), exception);
 
  double rms2 = 0.0;
  unsigned int rmsDiv = 0;

  for (unsigned int i = 0; i < index2 - index1-1; i++)    {
      double distN = const_cast<Atom&>(sp.getAmino(index1+i+1)[N]).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[N]));
      double distCA = const_cast<Atom&>(sp.getAmino(index1+i+1)[CA]).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[CA]));
      double distC = const_cast<Atom&>(sp.getAmino(index1+i+1)[C]).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[C]));
      rms2 += sqr(distN) + sqr(distCA) + sqr(distC);
      rmsDiv += 3;
      
      if (withOxygen)	{
	  double distO = const_cast<Atom&>(sp.getAmino(index1+i+1)[O]
					   ).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[O]));
	  rms2 += sqr(distO);
	  rmsDiv += 1;
	}
    }
  
  unsigned int index = index2 - index1 - 1;
  double distN = const_cast<Atom&>(sp.getAmino(index1+index+1)[N]).distance(
				   const_cast<Atom&>(sp2.getAmino(index)[N]));
  double distCA = const_cast<Atom&>(sp.getAmino(index1+index+1)[CA]).distance(
				   const_cast<Atom&>(sp2.getAmino(index)[CA]));
  double distC = const_cast<Atom&>(sp.getAmino(index1+index+1)[C]).distance(
				   const_cast<Atom&>(sp2.getAmino(index)[C]));
  double rms = rms2 + sqr(distN) + sqr(distCA) + sqr(distC);
  rmsDiv += 3;

  if (output)
    if (pVerbose)      {
	cout << "global RMS= " << setw(6) << setprecision(3)
	     << sqrt(rms2/(rmsDiv-3)) << "   (" << setw(6) << setprecision(3) 
	     << sqrt(rms/rmsDiv) << ")\n";
      }
	
  return sqrt(rms2/(rmsDiv-3));
}



/**
 * @Description   calculates the end RMS  
 * @param    Spacer reference that contains the protein (const Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution ( const Spacer&), output flag (bool), oxygen flag (bool)
 * @return   corresponding value(double)
 */
double LoopModel::calculateEndRms(const Spacer& sp, unsigned int index1, 
unsigned int index2, const Spacer& sp2){

  IntCoordConverter icc;
  unsigned int offset = index2 - index1 - 1;
  double rmsE = 
    sqr(const_cast<Atom&>(sp.getAmino(index2)[N]).distance(
			  const_cast<Atom&>(sp2.getAmino(offset)[N]))) 
    + sqr( const_cast<Atom&>(sp.getAmino(index2)[CA]).distance(
				const_cast<Atom&>(sp2.getAmino(offset)[CA])))
    + sqr( const_cast<Atom&>(sp.getAmino(index2)[C]).distance(
				const_cast<Atom&>(sp2.getAmino(offset)[C])));

  return sqrt(rmsE/3);
}

/**
 * @Description  calculates the  RMS  and if set can print the values, considering or not the oxygen
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&) , flag for the outout(bool), flag for the oxygen(bool)
 * @return   corresponding value(double)
 */
double LoopModel::calculateRms(const Spacer& sp, unsigned int index1, 
unsigned int index2, const Spacer& sp2, bool output, bool withOxygen){
  PRECOND( index2 - index1 == sp2.size(), exception);
 
  double rms2 = 0.0;
  unsigned int rmsDiv = 0;

  for (unsigned int i = 0; i < index2 - index1-1; i++)    {
      double distN = const_cast<Atom&>(sp.getAmino(index1+i+1)[N]).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[N]));
      double distCA = const_cast<Atom&>(sp.getAmino(index1+i+1)[CA]).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[CA]));
      double distC = const_cast<Atom&>(sp.getAmino(index1+i+1)[C]).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[C]));
      rms2 += sqr(distN) + sqr(distCA) + sqr(distC);
      rmsDiv += 3;

      if (withOxygen)
	{
	  double distO = const_cast<Atom&>(sp.getAmino(index1+i+1)[O]
					   ).distance(
				       const_cast<Atom&>(sp2.getAmino(i)[O]));
	  rms2 += sqr(distO);
	  rmsDiv += 1;
	}
    }

  unsigned int index = index2 - index1 - 1;
  double distN = const_cast<Atom&>(sp.getAmino(index1+index+1)[N]).distance(
				   const_cast<Atom&>(sp2.getAmino(index)[N]));
  double distCA = const_cast<Atom&>(sp.getAmino(index1+index+1)[CA]).distance(
				   const_cast<Atom&>(sp2.getAmino(index)[CA]));
  double distC = const_cast<Atom&>(sp.getAmino(index1+index+1)[C]).distance(
				   const_cast<Atom&>(sp2.getAmino(index)[C]));
  double rms = rms2 + sqr(distN) + sqr(distCA) + sqr(distC);
  rmsDiv += 3;

  if (output)
    if (pVerbose > 0)      {
	cout << "global RMS= " << setw(6) << setprecision(3)
	     << sqrt(rms2/(rmsDiv-3)) << "   (" << setw(6) << setprecision(3) 
	     << sqrt(rms/rmsDiv) << ")\t";
      }

  IntCoordConverter icc;
  unsigned int offset = index2 - index1 - 1;
  double rmsE = calculateEndRms(sp, index1, index2, sp2);
   
  if (output)
    if (pVerbose > 0)      {
	cout << "end-RMS= " << setw(6) << setprecision(3) << rmsE << "\t";
	cout << "  " << setw(6) << setprecision(3) 
	     << icc.getBondLength(const_cast<Atom&>(sp2.getAmino(offset)[N]), 
				const_cast<Atom&>(sp2.getAmino(offset)[CA]));
	cout << "  " << setw(6) << setprecision(3) << RAD2DEG * 
	  icc.getBondAngle(const_cast<Atom&>(sp2.getAmino(offset-1)[C]), 
			 const_cast<Atom&>(sp2.getAmino(offset)[N]),
			 const_cast<Atom&>(sp2.getAmino(offset)[CA]) );
	cout << "  " << setw(6) << setprecision(3) << RAD2DEG * 
	  icc.getTorsionAngle(const_cast<Atom&>(sp2.getAmino(offset-1)[CA]),  
			    const_cast<Atom&>(sp2.getAmino(offset-1)[C]), 
			    const_cast<Atom&>(sp2.getAmino(offset)[N]),
			    const_cast<Atom&>(sp2.getAmino(offset)[CA]) );
	cout << "\n";
      }

  return rmsE;
}

/**
 * @Description  Sets the structure
 * @param    Spacer reference that contains the protein (Spacer&), reference that contains the possible solution (  Spacer&) ,
 * Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int) 
 * @return   changes the object internally (void)  
 */
void LoopModel::setStructure(Spacer& sp, Spacer& sp2, unsigned int index1, 
			unsigned int index2){
  PRECOND(sp2.size() == index2 - index1, exception);

  for (unsigned int i = 0; i < sp2.size(); i++)    {
      sp.getAmino(index1+1+i).setType(sp2.getAmino(i).getType());

      sp.getAmino(index1+1+i)[N].setCoords(sp2.getAmino(i)[N].getCoords());
      sp.getAmino(index1+1+i)[CA].setCoords(sp2.getAmino(i)[CA].getCoords());
      sp.getAmino(index1+1+i)[C].setCoords(sp2.getAmino(i)[C].getCoords());

      pSetSideChain(sp.getAmino(index1+i));
    }

  for (unsigned int i = 0; i < sp2.size(); i++)    {
      if (sp.getAmino(index1+1+i).isMember(O))
	sp.getAmino(index1+1+i).removeAtom(
		    sp.getAmino(index1+1+i)[O]);

      sp.getAmino(index1+1+i).addMissingO();

      if (sp.getAmino(index1+1+i).isMember(O))
	sp.getAmino(index1+1+i)[O].setNumber(sp.getAmino(
	      index1+1+i)[C].getNumber()+1);
    }

}


/**
 * @Description  Calculates the propensities for a protein
 * @param   reference to the protein(Spacer&)
 * @return   corresponding value(double)
 */
double LoopModel::calculatePropensities(Spacer& sp2) {
  return tor.calculateEnergy(sp2);
}


/**
 * @Description  Calculates the propensities for the loop
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&)
 * @return   corresponding value(double)
 */
double LoopModel::calculatePropensities(const Spacer& sp, 
unsigned int index1, unsigned int index2, Spacer& sp2) {
  double partRes = 0.0; 
  IntCoordConverter icc;

  // propensity of sp2(0)
  partRes += tor.calculateEnergy(const_cast<AminoAcid&>(sp.getAmino(0)),
				 RAD2DEG * icc.getTorsionAngle(const_cast<Atom&>(
					       sp.getAmino(index1)[C]), sp2.getAmino(0)[N], 
					       sp2.getAmino(0)[CA], sp2.getAmino(0)[C]), 
				 sp2.getAmino(0).getPsi(true));

  return partRes + tor.calculateEnergy(sp2);
}
 

  


/**
 * @Description  This routine calculates the compactness of an given solution. It simply
*  determines for each atom in the loop the three least distances to the 
*  framework. If the given solution is compact the calculated distances
*  (and thus the return value which is simply the sum of all distances)
*  will be small.
 *  @param    Spacer reference that contains the loop (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the protein (  Spacer&). index1 and index2 are expected to start counting from 0! Max loop lenght=25
 * @return   corresponding value(double)
 */
double LoopModel::compactness(Spacer& loop, unsigned int index1, 
unsigned int index2, Spacer& proteine) {

  double ret_value = 0.0;          // used to store the return value
  const int NUMBER_OF_ATOMS = 100; // the maximal number of atoms in the loop
  const int BEST_VALUES = 3;    
  // the n best distance values are added to the result
  double best_distances[NUMBER_OF_ATOMS][BEST_VALUES];       
  // here we store the best distances for each atom in the loop
  std::vector<Atom*> loopAtoms;    
  // contains all the atoms in the loop (just the backbone)
  std::vector<Atom*> frameAtoms;    
  // contains all the atoms of an aminoacid in the framework
  int countLoopAtoms = 0;                                    
  // the number of the atoms in the loop (just the backbone)
  int countFrameAtoms = 0;                                   
  // the number of the atoms in an framework aminoacid
  double distance = 0;                                       
  // distance between two atoms
  double swap = 0;                                           
  // used to update the best_distance array
  double swap2 = 0;                // dito


  //initialize the array best_distances
  for (int i = 0; i < NUMBER_OF_ATOMS; i += 1)     {
      for (int j = 0; j < BEST_VALUES; j += 1)	{
	  best_distances[i][j] = 100000;
	}
    }

  // get the number of atoms in the loop and fill the vector loopAtoms
  PRECOND(index2 - index1 == loop.sizeAmino(), exception);
  for (unsigned int i = 0; i < loop.sizeAmino(); i += 1)     {
      AminoAcid& aa = loop.getAmino(i); 
      
      //put all the atoms of the aminoacid into the vector loopAtoms
      if (aa.isMember(N)) {
	loopAtoms.push_back(&(aa[N]));
	countLoopAtoms += 1;
      }
      if (aa.isMember(C)) {
	loopAtoms.push_back(&(aa[C]));
      	countLoopAtoms += 1;
      }
      if (aa.isMember(CA)) {
	loopAtoms.push_back(&(aa[CA]));
	countLoopAtoms += 1;
      }
      if (aa.isMember(O)) {
	loopAtoms.push_back(&(aa[O]));
	countLoopAtoms += 1;
      }
    }

  // calculate the compactness by iterating over all the atoms in the framework
  for (unsigned int i = 0; i < proteine.size(); i++)    {
      if (i + 2 > index1 && i < index2 + 1) {                
	// we are currently in the loop section of the proteine
	continue;                                            
	// or at the edge of the proteine to the loop
      }

      // now we put all the Atoms of the current aminoacid in the 
      // framework into an vector
      AminoAcid& aa = proteine.getAmino(i);
       
      if (aa.isMember(N)) {
	frameAtoms.push_back(&(aa[N]));
	countFrameAtoms += 1;
      }
      if (aa.isMember(C)) {
	frameAtoms.push_back(&(aa[C]));
      	countFrameAtoms += 1;
      }
      if (aa.isMember(CA)) {
	frameAtoms.push_back(&(aa[CA]));
	countFrameAtoms += 1;
      }
      if (aa.isMember(O)) {
	frameAtoms.push_back(&(aa[O]));
	countFrameAtoms += 1;
      }
      
      //iterate over all atoms in the loop and the current frame amino acid
      for (int j = 0; j < countLoopAtoms; j += 1)	{
	  for (int k = 0; k < countFrameAtoms; k += 1) 	    {
              // get the coordinates of the two relevant atoms
	      vgVector3<double> fv = frameAtoms[k]->getCoords();
	      vgVector3<double> lv = loopAtoms[j]->getCoords();
	      distance = sDistance(fv, lv);
	      
	      // update the best_values array...
	      if (best_distances[j][0] > distance) {
		  // distance gets into the new top position
		  swap = best_distances[j][0];
		  swap2 = best_distances[j][1];
		  best_distances[j][0] = distance;
		  best_distances[j][1] = swap;
		  best_distances[j][2] = swap2;
		  continue;
		}
	      if (best_distances[j][1] > distance){
		  // distance gets into the middle position
		  swap = best_distances[j][1];
		  best_distances[j][1] = distance;
		  best_distances[j][2] = swap;
		  continue;
		}
	      if (best_distances[j][2] > distance)	{
		  // distance gets into the last position
		  best_distances[j][2] = distance;
		  continue;
		}
	    }
	}
      countFrameAtoms = 0;
      frameAtoms.clear();
    }
  // now calculate the return value
  for (int j = 0; j < countLoopAtoms; j += 1)   {
      for (int k = 0; k < 3; k += 1) {
	  ret_value += best_distances[j][k];
	}
    }
  return ret_value;
}


/**
 * @Description  prints the matrix
 * @param   structure that constains the matrix info(vector<double>)
 * @return   prints the results(void)
 */
static void sPrintMatrix(std::vector<double> span){
  cout << ">";
  for (unsigned int i = 0; i < span.size(); i++)
    cout << setw(5) << i << "\t";
  cout << "\n>";
  for (unsigned int i = 0; i < span.size(); i++)
    cout << setw(5) << setprecision(3) << span[i] << "\t"; 
  cout << "\n";
}

/**
 * @Description  calculates the minimum distance considering the given data
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *  spanning value(unsigned int), distance(double), verbose flag(bool)
 * @return   corresponding value (static double)
 */
static double sDeletionAnchors(Spacer& sp, unsigned int& index1, 
unsigned int& index2, unsigned int span, double dist, bool pVerbose){
  unsigned int span2 = 2 + static_cast<unsigned int>( dist / 2.0 );

  std::vector<double> spanning;
  unsigned int minI = (index1 > span2) ? index1 
    - span2 + 2 : 1;
  unsigned int maxI = (index1 + span2 < sp.sizeAmino()) 
    ? index1 - 1 : sp.sizeAmino() - span2 - 1;

  if (pVerbose)
    cout << "minI = " << sp.getPdbNumberFromIndex(minI) 
         << "\tmaxI = " << sp.getPdbNumberFromIndex(maxI) << endl;

  for (unsigned int i = minI; i <= maxI; i++)    {
      double tmp = sp.getAmino(i)[CA].distance(sp.getAmino(i + span2)[CA]);

      if ((sp.getAmino(i).getState() == HELIX) 
	  || (sp.getAmino(i).getState() == STRAND)){
	  tmp += WEIGHT_SEC;

	  for (unsigned j = 0; j < spanning.size(); j++)
	    spanning[j] += WEIGHT_SEC;
	}

      if ((sp.getAmino(i + span2).getState() == HELIX) 
	  || (sp.getAmino(i + span2).getState() == STRAND)){
	  for (unsigned j = i; j <= maxI; j++)
	    tmp += WEIGHT_SEC;
	}

      spanning.push_back(tmp);
      
    }

  if (pVerbose)
    sPrintMatrix(spanning);
  
  unsigned int optI = 0;
  double minDist = 999.9;
  
  for (unsigned int i = 0; i < spanning.size(); i++)
      if (spanning[i] < minDist){
	  optI = minI + i;
	  minDist = spanning[i];
	}
  
  // set 'correct' anchors:
  index1 = optI;
  index2 = optI+span2;  
  return minDist;
}

/**
 * @Description  calculates the maximum distance considering the given data
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int),
 *  spanning value(unsigned int), distance(double), verbose flag(bool)
 * @return   corresponding value(static double)
 */
static double sInsertionAnchors(Spacer& sp, unsigned int& index1, 
unsigned int& index2, unsigned int span, double dist, bool pVerbose){
  unsigned int span2 = static_cast<unsigned int>( (span + 1) / 2 );

  std::vector<double> spanning;
  unsigned int minI = (index1 > span2) ? index1 
    - span + 2 : 1;
  unsigned int maxI = (index1 + span2 < sp.sizeAmino()) 
    ? index1 - 1 : sp.sizeAmino() - span2 - span - 1;

  cout << "minI = " << sp.getPdbNumberFromIndex(minI) 
       << "\tmaxI = " << sp.getPdbNumberFromIndex(maxI) << endl;

  for (unsigned int i = minI; i <= maxI; i++) {
      if (i + span + span2 + 2 > sp.sizeAmino())
	break;

      double tmp = sp.getAmino(i)[CA].distance(sp.getAmino(i 
						     + span + span2)[CA]);

      if ((sp.getAmino(i).getState() == HELIX) 
	  || (sp.getAmino(i).getState() == STRAND)){
	  tmp -= WEIGHT_SEC;

	  for (unsigned j = 0; j < spanning.size(); j++)
	    spanning[j] -= WEIGHT_SEC;
	}

      if (i + span2 >= sp.sizeAmino())
	ERROR("Index out of scope.", exception);

      if ((sp.getAmino(i + span2).getState() == HELIX) 
	  || (sp.getAmino(i + span2).getState() == STRAND)){
	  for (unsigned j = i; j <= maxI; j++)
	    tmp -= WEIGHT_SEC;
	}

      spanning.push_back(tmp);
    }
  
  sPrintMatrix(spanning);
  
  unsigned int optI = 0;
  double maxDist = -999.9;
  
  for (unsigned int i = 0; i < spanning.size(); i++)
      if (spanning[i] > maxDist){
	  optI = minI + i;
	  maxDist = spanning[i];
	}

  // set 'correct' anchors:
  index1 = optI;
  index2 = optI + span + span2;

  return maxDist;
}

/**
 * @Description  calculates the maximum/minimum distance considering the given data depending in the flag
 * @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int), 
 * flag for defining as deletion(bool)
 * @return   changes the object internally (void)  
 */
void LoopModel::defineLoopAnchors(Spacer& sp, unsigned int& index1, unsigned int& index2, bool isDeletion){
  unsigned int span = index2 - index1; 
  double dist = 999.9;
  if (isDeletion)
    dist = sp.getAmino(index1)[CA].distance(sp.getAmino(index1+1)[CA]);
  else
    dist = sp.getAmino(index1)[CA].distance(sp.getAmino(index2)[CA]);

  double maxDist = 0.0;

  if (pVerbose)
    cout << "Index1= " << index1 << "  Index2= " << index2
      << " Span= " << span << "\t Dist= " << dist << "\n";
  
  if (isDeletion)
    maxDist = sDeletionAnchors(sp, index1, index2, span, dist, pVerbose);
  else // isInsertion
    maxDist = sInsertionAnchors(sp, index1, index2, span, dist, pVerbose);

  cout << "---------------------------------------------------\n";
  cout << "Optimized: \t dist = " << setw(5) << setprecision(3) 
       << dist << " (Pen. = " << setw(5) << setprecision(3) << maxDist
       << ")\t loop length = " << index2-index1
       << "\n Anchors: \t N-term. = " <<  sp.getPdbNumberFromIndex(index1) 
       << "\t\t C-term. = " << sp.getPdbNumberFromIndex(index2) << endl;
}


// MODIFIERS:

/**
 * @Description  Copies the original loop in to a new object
 * @param   reference to the original loop model(const LoopModel&)
 * @return   changes the object internally (void)  
 */
void LoopModel::copy(const LoopModel& orig){
  PRINT_NAME; 

  pInter = orig.pInter; 
  pVerbose = orig.pVerbose; 
  pPlot = orig.pPlot; 

  //  deep copy of: table = orig.table;
  table.clear();
  for (unsigned int i = 0; i < orig.table.size(); i++)    {
      LoopTable* tmp = new LoopTable;
      tmp = orig.table[i];
      table.push_back(tmp);
    }

  tableFileName.clear();
  tableFileName = orig.tableFileName;
  solution = orig.solution;

  //  deep copy of: ENDRMS_WEIGTH = orig.ENDRMS_WEIGTH;
  ENDRMS_WEIGTH.clear();
  for (unsigned int i = 0; i < orig.ENDRMS_WEIGTH.size(); i++)
    ENDRMS_WEIGTH.push_back(orig.ENDRMS_WEIGTH[i]);

}
 
/**
 * @Description  Creates the loop model 
 * @param  the start amino acid reference and the corresponding N coord, the start amino acid reference and the corresponding N coord, 
 * indexes for the anchor regions(unsigned int, unsigned int) , numbers of loops(unsigned int ,unsigned int)
 * @return  the posibles solutions (vector<Spacer>  )
 */
std::vector<Spacer>  LoopModel::createLoopModel(const AminoAcid& start, const vgVector3<double>& 
startN, const AminoAcid& end, const vgVector3<double>& endN, unsigned int 
indexS, unsigned int indexE, unsigned int numLoops, unsigned int numLoops2,
std::vector<string> typeVec){
  PRECOND((indexS > 0) && (indexE > 0), exception);
  LoopTableEntry startEntry, endEntry;

  VectorTransformation vt;
  // prepare input data:
  convertCoords(start, convert(startN), end, convert(endN), startEntry, 
		endEntry, vt, indexS);
  
  if (pVerbose)
    cout << "loading...\n";
  loadTables(indexE-indexS);
  solution.clear();

  if (pVerbose)
    cout << "searching...\n";
  ringClosureBase(startEntry, endEntry, indexS, indexE-indexS, 0.0, 
		  vt, numLoops, numLoops2, solution);

  if (pVerbose)
    cout << "returned...\n";
  // prepare output data:
  return calculateLoop(convert(startN), indexE-indexS, typeVec);

}


/**
 * @Description  Creates temporary clusters for the loops
 * @param   the posibles loops (vector<Spacer>  )
 * @return   changes the object internally (void)  
 */
void LoopModel::clusterLoops(std::vector<Spacer>& solVec){
  std::vector<Spacer> tmpSolVec;

  tmpSolVec.push_back(solVec[0]);

  for (unsigned int i = 1; i < solVec.size(); i++)    {
      bool contained = false;

      for (unsigned int j = 0; j < tmpSolVec.size(); j++)
	if (pCalculateLoopRms(solVec[i], tmpSolVec[j]) < SIM_LIMIT) {
	    contained = true;
	    break;
	  }

      if (!contained)
	tmpSolVec.push_back(solVec[i]);
    }
  
  if (pVerbose > 1)  {
   cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
   cout << "A total of " << setw(5) << (solVec.size() - tmpSolVec.size()) 
	<< " solutions were clustered at " << setw(4) << setprecision(2) 
	<< SIM_LIMIT << " A distance.\n";
   cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
  }

  solVec.clear();
  solVec = tmpSolVec;

}
 
/**
 * @Description  Calculates the van-der-Waals forces between all the atoms in the two amino
*  acids. If two atoms are two close to each other (=> collision) a large 
*  van-der-Waals term is returned. If two atoms are bonded, no van-der-Waals term is calculated.
 * @param    pos1 is the position of amino acid 1 in the protein and pos2 is the
 *  positition of amino acid 2 in the protein(int). Start counting from 1, the threshold distance below which we have a collision(double).
 * @return   number of collisions(int)
 */
int LoopModel::amino_amino_collision(AminoAcid& aa1, AminoAcid& aa2, int pos1, 
int pos2, double distThreshold) {
  std::vector<Atom*> v1;  // this vector contains all the atoms of amino acid 1
  std::vector<Atom*> v2;  // vector contains all the atoms of amino acid 2
  double distance;   // contains the distance between two atoms
  int ret_value = 0; // adds up the returned van-der-Waals forces
  const int COLLISON_WEIGHT = 1;

  // get the relevant atoms of amino acid 1 and amino acid 2 into their vectors
  if (aa1.isMember(N))
    v1.push_back(&(aa1[N]));
  if (aa1.isMember(C))
    v1.push_back(&(aa1[C]));
  if (aa1.isMember(CA))
    v1.push_back(&(aa1[CA]));
  if (aa1.isMember(O))
    v1.push_back(&(aa1[O]));
   if (aa1.getSideChain().isMember(CB))
     v1.push_back(&(aa1.getSideChain()[CB]));

  if (aa2.isMember(N))
    v2.push_back(&(aa2[N]));
  if (aa2.isMember(C))
    v2.push_back(&(aa2[C]));
  if (aa2.isMember(CA))
    v2.push_back(&(aa2[CA]));
  if (aa2.isMember(O))
    v2.push_back(&(aa2[O]));
  if (aa2.getSideChain().isMember(CB))
    v2.push_back(&(aa2.getSideChain()[CB]));

  for (unsigned int i = 0; i < v1.size(); i += 1)     {
      vgVector3<double> xv = v1[i]->getCoords();
      for (unsigned int j = 0; j < v2.size(); j += 1) 	{
	  // we have to check whether C is bound to N (in this case we don't 
	  // check collisions)
	  if(pos1 + 1 == pos2 && v1[i]->getCode() == C && v2[j]->getCode() == N)
	    continue;
	  if(pos2 + 1 == pos1 && v1[i]->getCode() == N && v2[j]->getCode() == C)
	    continue;
	  
	  vgVector3<double> yv = v2[j]->getCoords();
	  distance = sDistance(xv, yv);
	  
	  if( distance < distThreshold) 
	    ret_value += COLLISON_WEIGHT;
	}
    }
  return ret_value;
}

// HELPERS: 
 
/**
 * @Description  converts the input coords to something processable.
 *    Converts the aminoacids start and end into the looptable entries 
 *    startEntry and endEntry, which can be passed on to ringClosure.
 * @param  starting amino acid and N coords (AminoAcid,vector3<float>, ending amino acid and N coords (AminoAcid,vector3<float>,
 * reference to the starting and ending Loop Table Entry(LoopTableEntry&), reference to the transformation vector( VectorTransformation&), unsigned int
 * @return   changes the object internally (void)  
 */
void LoopModel::convertCoords(AminoAcid start, vgVector3<float> startN, 
AminoAcid end, vgVector3<float> endN, LoopTableEntry& startEntry, 
LoopTableEntry& endEntry, VectorTransformation& vt, unsigned int nAmino){
  if (startN.length() != 0.0) {// special case if there is no 'start' offset
     // adjust positions:

      vgMatrix3<float> rotMat;
      
      // 1. set N of start to (0,0,0):

      vgVector3<float> newOrigin = convert(start[N].getCoords());
      pAddTrans(start, startN, end, endN, newOrigin);
      vt.addTrans(newOrigin);
      
      // 2. apply rotation to bring C'prev-N-CA into the x,y plane:

      vgVector3<float> zNormal(0,0,1);
      vgVector3<float> startNormal = convert(start[N].getTrans()).cross(
		  convert(start[CA].getCoords() - start[N].getCoords())
		  ).normalize();
      alignVectors(zNormal, startNormal, rotMat);
      vt.addAlignVectors(startNormal, zNormal);
      pAddRot(start, startN, end, endN, rotMat);
      
      // 3. apply rotation to bring CA to (0, bl, 0):
      
      vgVector3<float> refAxis(0,1.46,0);
      alignVectors(refAxis, convert(start[CA].getTrans()), rotMat);
      vt.addAlignVectors(convert(start[CA].getTrans()), refAxis);
      pAddRot(start, startN, end, endN, rotMat);
 
    }

  // set data...

  startEntry.endPoint = convert(start[C].getCoords());
  startEntry.endDirection = startN - convert(start[C].getCoords());
  startEntry.endNormal = (convert(start[CA].getCoords())
			  - convert(start[C].getCoords())).cross(
			   startEntry.endDirection).normalize();

  endEntry.endPoint = convert(end[C].getCoords());
  endEntry.endDirection = endN - convert(end[C].getCoords());
  endEntry.endNormal = (convert(end[CA].getCoords()-end[C].getCoords())).cross(
			   endEntry.endDirection).normalize();
}

 
/**
 * @Description  tries to find a ring closure between to locations.Tries to find a ring closure between source and destination
*    num is the starting residue number offset; offsetPhi is used 
*    (internally) to represent the already computed phi angle offset. 
 * 
 * @param  reference of the source and destination loop table entry (const LoopTableEntry&,LoopTableEntry&) destination, 
 * chain number (unsigned int),number of amino acids( unsigned int ),phi offset( double ),Transformation vector( VectorTransformation),
 * partial solution(std::vector<vgVector3<float> >&),current selection(unsigned int).
 * @return   flag that confirms the  existance of a ring closure (bool)
 */
bool LoopModel::ringClosure(const LoopTableEntry& source, 
const LoopTableEntry& destination, unsigned int nStart, 
unsigned int nAminoAcids, double offsetPhi, VectorTransformation vt,
std::vector<vgVector3<float> >& partialSolution, unsigned int currentSelection){
  PRECOND(((nAminoAcids > 0) && (nAminoAcids <= MAX_CHAIN_LENGTH)), exception);
  
  LoopTableEntry tmpSrc = source;
  LoopTableEntry tmpEnd = destination;

  vgVector3<float> refNull(0,0,0);

  if (tmpSrc.endDirection != refNull)
    tmpEnd = tmpSrc.setToOrigin(tmpEnd, nStart, vt);

  // rotate the destination back into the X,Y plane relative to the origin:
  vg_ieee64 phi = -(tmpEnd.rotateIntoXYPlane(vt));

  if (nAminoAcids == 1)    { 
      // calculate result and store it
      vgVector3<float> posCA(0,0,0), posC(0,0,0), posN(0,0,0);

      posC = vt.transform(tmpEnd.endPoint);
      posN = vt.transform(tmpEnd.endPoint + tmpEnd.endDirection);
      posCA = vt.transform(tmpEnd.endPoint + 
		 LoopTableEntry::BOND_LENGTH_CALPHA_TO_CPRIME * 
		 ( vgMatrix3<float>::createRotationMatrix(tmpEnd.endNormal,
		    -(DEG2RAD * LoopTableEntry::BOND_ANGLE_AT_CPRIME_TO_N)) 
		    * tmpEnd.endDirection).normalize() );

      partialSolution.push_back(posCA);
      partialSolution.push_back(posC);
      partialSolution.push_back(posN);

      return 1;
    }

  // divide ...
  int midAminoAcids = (nAminoAcids / 2) + (nAminoAcids % 2); // eg. 5 = 3 + 2

  INVARIANT( (nAminoAcids < table.size()) && (table[nAminoAcids] != NULL),
	     exception);
  std::vector<LoopTableEntry> tmpMidVec = table[nAminoAcids]->getNClosest(tmpEnd, 
					currentSelection+1, nAminoAcids);
  LoopTableEntry newMiddle;

  if (tmpMidVec.size() > currentSelection)
    newMiddle = tmpMidVec[currentSelection];
  else
     newMiddle = tmpMidVec[0];
   
  newMiddle.endPoint = newMiddle.midPoint;
  newMiddle.endDirection = newMiddle.midDirection;
  newMiddle.endNormal = newMiddle.midNormal;

  // define the origin:
  LoopTableEntry origin;
  origin.endDirection = refNull;

  // ... & conquer:
  return ( ringClosure(origin, newMiddle, 1, midAminoAcids, offsetPhi+phi, vt,
		       partialSolution)
	   * ringClosure(newMiddle, tmpEnd, midAminoAcids+1, 
			 (nAminoAcids / 2), 0, vt, partialSolution) );
}

/**
 * @Description  tries to find a ring closure between to locations.Tries to find a ring closure between source and destination
*    num is the starting residue number offset; offsetPhi is used 
*    (internally) to represent the already computed phi angle offset. 
 * 
 * @param  reference of the source and destination loop table entry (const LoopTableEntry&,LoopTableEntry&) destination, 
 * chain number (unsigned int),number of amino acids( unsigned int ),phi offset( double ),Transformation vector( VectorTransformation),
 * partial solution(std::vector<vgVector3<float> >&),current selection(unsigned int).
 * @return   flag that confirms the  existance of a ring closure base(bool)
 */
bool LoopModel::ringClosureBase(const LoopTableEntry& source, 
const LoopTableEntry& destination, unsigned int nStart, 
unsigned int nAminoAcids, double offsetPhi, VectorTransformation vt, 
unsigned int num, unsigned int depth, std::vector<vgVector3<float> >& partialSolution){
  PRECOND(((nAminoAcids > 1) && (nAminoAcids <= MAX_CHAIN_LENGTH)), exception);
  
  if (pVerbose)    {
      cout << "------------------------------------\n";
      cout << "generating " << num << " models...\n";
      cout << "--> rc: nStart = " << nStart << "   nAmino = " 
	   << nAminoAcids << "\n";
    }

  LoopTableEntry tmpSrc = source;
  LoopTableEntry tmpEnd = destination;

  vgVector3<float> refNull(0,0,0);

  if (tmpSrc.endDirection != refNull)
    tmpEnd = tmpSrc.setToOrigin(tmpEnd, nStart, vt);

  // rotate the destination back into the X,Y plane relative to the origin:
  vg_ieee64 phi = -(tmpEnd.rotateIntoXYPlane(vt));

  if (nAminoAcids == 1)    { 
      // calculate result and store it
      vgVector3<float> posCA(0,0,0), posC(0,0,0), posN(0,0,0);

      posC = vt.transform(tmpEnd.endPoint);
      posN = vt.transform(tmpEnd.endPoint + tmpEnd.endDirection);
      posCA = vt.transform( tmpEnd.endPoint + 
		 LoopTableEntry::BOND_LENGTH_CALPHA_TO_CPRIME * 
		 ( vgMatrix3<float>::createRotationMatrix(tmpEnd.endNormal,
		    -(DEG2RAD * LoopTableEntry::BOND_ANGLE_AT_CPRIME_TO_N)) 
		    * tmpEnd.endDirection).normalize() );

      partialSolution.push_back(posCA);
      partialSolution.push_back(posC);
      partialSolution.push_back(posN);

      return 1;
    }

  // divide ...
  int midAminoAcids = (nAminoAcids / 2) + (nAminoAcids % 2); // eg. 5 = 3 + 2

  INVARIANT( (nAminoAcids < table.size()) && (table[nAminoAcids] != NULL),
	     exception);

  std::vector <LoopTableEntry> middle; 
  middle = table[nAminoAcids]->getNClosest(tmpEnd, num, nAminoAcids);

  // define the origin:
  LoopTableEntry origin;
  origin.endDirection = refNull;

  for (unsigned int i = 0; i < middle.size(); i++)    {
      middle[i].endPoint = middle[i].midPoint;
      middle[i].endDirection = middle[i].midDirection;
      middle[i].endNormal = middle[i].midNormal;
      
      // ... & conquer:

      for (unsigned int j = 0; j < depth; j++)	{
	  ringClosure(origin, middle[i], 1, midAminoAcids, offsetPhi+phi, 
		      vt, partialSolution, j);
	  ringClosure(middle[i], tmpEnd, midAminoAcids+1, (nAminoAcids / 2), 
		      0, vt, partialSolution, j);
	}
    }
	return 0;
}



/**
 * @Description  Calculates the Loop
 * @param   N coords (const vgVector3<float>& ), length (unsigned int) , vector containing the amino acid types (std::vector<string> )
 * @return  possible solutions(vector<Spacer>) 
 */
std::vector<Spacer> LoopModel::calculateLoop(const vgVector3<float>& startN, unsigned int length,
std::vector<string> typeVec){
  // build solution spacer:
  
  std::vector<Spacer> result;
  unsigned int num = solution.size() / (length*3);

  if (length != typeVec.size())
    ERROR("typeVec does not match number of aminoacids.", exception);

  for (unsigned int i = 0; i < num; i++)    {
      // determine rough ''quality'' of loop to store in B-factors:
      double bfac = 50.0;

      unsigned int offset = i * (length*3);

      Spacer sp;
      
      sp.setType("Loop Model");
      
      AminoAcid* prev = new AminoAcid;
      pAminoAcidSetup(prev, typeVec[0], bfac);

      (*prev)[N].setCoords(convert(startN));
      (*prev)[CA].setCoords(convert(solution[offset + 0]));
      (*prev)[C].setCoords(convert(solution[offset + 1]));

      if (typeVec[0] != "GLY")
	(*prev).patchBetaPosition();

      sp.insertComponent(prev);
      
      for (unsigned int j = 1; j < length; j++)	{
	  // determine rough ''quality'' of loop to store in B-factors:
	  bfac = (50 + 5 * length) - 5 * fabs(j - length);

	  AminoAcid* aa = new AminoAcid;
	  pAminoAcidSetup(aa, typeVec[j], bfac);
	  
	  (*aa)[N].bindIn((*prev)[C]);

	  (*aa)[N].setCoords(convert(solution[offset + j*3-1]));
	  (*aa)[CA].setCoords(convert(solution[offset + j*3]));
	  (*aa)[C].setCoords(convert(solution[offset + j*3+1]));
	 
  	  (*prev).bindOut((*prev)[C],(*aa),(*aa)[N]);
	  (*prev).addMissingO();

	  if (typeVec[j] != "GLY")
	    (*prev).patchBetaPosition();

	  sp.insertComponent(aa);
	  prev = aa;
	}
      
      result.push_back(sp);
    }

  return result;
}
  

/**
 * @Description  Load the really used tables .
 * @param  number of amino acids (unsigned int )
 * @return   changes the object internally (void)  
 */
void LoopModel::loadTables(unsigned int nAmino){
  PRECOND(tableFileName.size() >= nAmino, exception);

  if (nAmino <= 1) 
    return;  // no loading needed 

  unsigned int tmpSize = (table.size() > 0) ? table.size()-1 : 0;
  for(unsigned int loop = tmpSize; loop <= nAmino; loop++)
    table.push_back(NULL);

  std::vector<unsigned int> index;
  unsigned int offset = 1;

  while (nAmino / offset >= 1)    {
      if (nAmino % offset >= 1)
	index.push_back(nAmino / offset + 1);
      index.push_back(nAmino / offset);
      offset *= 2;
    }

  for (unsigned int i = 0; i < index.size(); i++)
    if ((table[index[i]] == NULL) && (index[i] > 1)){
	  if (pVerbose)
	    cout << "--> " << index[i] << "\n";

	  LoopTable* lt = new LoopTable;
	  lt->read(tableFileName[index[i]]); 
	  
	  VectorTransformation vt;
	  
	  for (unsigned int j = 0; j < lt->size(); j++)
	    (*lt)[j].rotateIntoXYPlane(vt);
	  
	  table[index[i]] = lt;        
	}
}

/**
 * @Description  Sets the side chain for an amino acid
 * @param   the reference of the amion acid(AminoAcid& )
 * @return   changes the object internally (void)  
 */
void LoopModel::pSetSideChain(AminoAcid& aa){
  if (!aa.getSideChain().isMember(CB))
    return;

  aa.removeSideChain();
  aa.patchBetaPosition();

  aa.getSideChain()[CB].setNumber(aa[C].getNumber()+2);
}


/**
 * @Description  This routine calculates the van-der-Waals forces of the loops contained
* in the vector. The van-der-Waals forces are calculated only between
* the amino acids in the loop. Not between the loop and the surrounding
* proteine. This is done in: LoopModel::loop_spacer_vdw()
 * @param reference of the protein(Spacer&),index (unsigned int)
 * @return   the corresponding value
 */
int LoopModel::loop_loop_vdw(Spacer& sp, unsigned int index1){
  int ret_value = 0;           // used to store the resulting vdw-forces
  const double CHECK_THRESHOLD = 6.0; 
  // all amino acids with a distance lower than that are checked

  for (unsigned int i = 0; i < sp.size(); i++)
    for (unsigned int j = i+1; j < sp.size(); j++){
	double dist = sDistance(sp.getAmino(i)[CA].getCoords(),
				sp.getAmino(j)[CA].getCoords());
	if (dist <= CHECK_THRESHOLD){
	    ret_value += amino_amino_collision(sp.getAmino(i), 
			  sp.getAmino(j), i + 1 + index1, j + 1 + index1, 2.0);
	  }
      }
  return ret_value;
}

/**
 * @Description  This routine calculates the van-der-Waals forces between the loop and the 
* surrounding proteine. Index1 and index2 demarcate the area of the loop in the proteine and thus
* these amino acids in the proteine are not used in the van-der-Waals
* calculation. It is assumed, that index1 and index2 start counting from
* 1. The van-der Waals forces between the amino acids in the loop itself is calculated in: LoopModel::loop_loop_vdw().
 *   @param    Spacer reference that contains the protein (Spacer&),Index1 and index2 demarcate the area of the loop(unsigned int,unsigned int)
 *   reference that contains the possible solution (  Spacer&)
 * @return   corresponding value (int)
 */
int LoopModel::loop_spacer_vdw(Spacer& loop, unsigned int index1, 
unsigned int index2, Spacer& proteine) 
{
  int ret_value = 0;            // used to store the resulting vdw-forces
  const double CHECK_THRESHOLD = 6.0;
  // all amino acids with a distance lower than that are checked

  for (unsigned int i = 0; i < proteine.size(); i++){
      if (i + 2 > index1 && i < index2 + 1)
	          // we are currently in the loop section of the proteine
	continue; // or at the edge of the proteine to the loop

      for (unsigned int j = 0; j < loop.size(); j++){
	  double dist = sDistance(proteine.getAmino(i)[CA].getCoords(),
				  loop.getAmino(j)[CA].getCoords());
	  if (dist <= CHECK_THRESHOLD)
	    ret_value += amino_amino_collision(proteine.getAmino(i), 
			    loop.getAmino(j), i + 1, index1 + j + 1, 2.0);
	}
    }
  return ret_value;
}
      

/**
 * @Description  Sets the amino acid type,bfactor,etc 
 * @param  amino acid pointer( AminoAcid*),amino acid type(string ),bfactor (double)
 * @return   changes the object internally (void)  
 */
void LoopModel::pAminoAcidSetup(AminoAcid* aa, string type, double bfac){
  Atom aN;
  aN.setCode(N);
  aN.setBFac(bfac);
  aa->addAtom(aN);
  Atom aCA;
  aCA.setCode(CA);
  aCA.setBFac(bfac);
  aa->addAtom(aCA);
  Atom aC;
  aC.setCode(C);
  aC.setBFac(bfac);
  aa->addAtom(aC);
  aa->setType(type);

  (*aa)[CA].bindIn((*aa)[N]);
  (*aa)[C].bindIn((*aa)[CA]);
}

/**
 * @Description  returns the ENDRMS_WEIGTH
 * @param   length(unsigned int)
 * @return   corresponding value(double)
 */
double LoopModel::getENDRMS_WEIGHT(unsigned int len){
  if (len >= ENDRMS_WEIGTH.size())
    ERROR("Length out of scope for ENDRMS_WEIGTH array.", exception);

  return ENDRMS_WEIGTH[len];
}

/**
 * @Description  sets the ENDRMS_WEIGTH
 * @param   value to set(double), length (unsigned int)
 * @return   changes the object internally (void)  
 */
void LoopModel::setENDRMS_WEIGHT(double val, unsigned int len){
  if (ENDRMS_WEIGTH.size() == 0)
    ERROR("ENDRMS_WEIGTH array is invalid, ie. empty?!", exception);

  if (len >= ENDRMS_WEIGTH.size())
    for (unsigned int i = ENDRMS_WEIGTH.size(); i <= len; i++)
      ENDRMS_WEIGTH.push_back(ENDRMS_WEIGTH[0]);

  ENDRMS_WEIGTH[len] = val;
}

/**
 * @Description  loads the ENDRMS_WEIGTH value from a istream
 * @param   file reference(istream&)
 * @return   changes the object internally (void)  
 */
void LoopModel::loadENDRMS_WEIGHT(istream& input){
  ENDRMS_WEIGTH.clear();

  while (input){
      double tmp;
      input >> tmp;

      if (!input) // Warning: This 'if' serves to remove a bug
	break;    // which causes the last line in a file to be read twice.

      ENDRMS_WEIGTH.push_back(tmp);
    }
}

/**
 * @Description  sets all the ENDRMS_WEIGTH
 * @param   vavlue to set, queantity of ENDRMS_WEIGTH to set (unsigned int)
 * @return   changes the object internally (void)  
 */
void LoopModel::setAllENDRMS_WEIGHT(double val, unsigned int max){
  ENDRMS_WEIGTH.clear();

  for (unsigned int i = 0; i <= max; i++)
    ENDRMS_WEIGTH.push_back(val);
    
}

/**
 * @Description  saves the ENDRMS_WEIGTH value in a ostream
 * @param   file reference(ostream&)
 * @return   changes the object internally (void)  
 */
void LoopModel::saveENDRMS_WEIGHT(ostream& output){
  for (unsigned int i = 0; i < ENDRMS_WEIGTH.size(); i++)
    output << ENDRMS_WEIGTH[i] << "\t";
  output << "\n";
}


/**
 * @Description calculates the energy between two amino acids  
 * @param   amino acids references( AminoAcid&, AminoAcid&)
 * @return   corresponding value(double)
 */
double LoopModel::sCalcEn(AminoAcid& aa, AminoAcid& aa2){
  double en = 0.0;

  string aaType = aa.getType();
  string aaType2 = aa2.getType();

  for (unsigned int j = 0; j < aa.size(); j++)
	for (unsigned int k = 0; k < aa2.size(); k++)
	  if ( (aa2[k].getCode() == N) || (aa2[k].getCode() == CA) 
	       || (aa2[k].getCode() == C) || (aa2[k].getCode() == O) 
	       || (aa2[k].getCode() == CB) )
	    en += rapdf.calculateEnergy(aa[j], aa2[k], aaType, aaType2);
  
  return en;
}

/**
 * @Description calculates the energy between two amino acids  
 * @param   amino acids references( AminoAcid&, AminoAcid&)
 * @return   corresponding value(double)
 */
double LoopModel::sICalcEn(AminoAcid& aa, AminoAcid& aa2){
  double en = 0.0;

  string aaType = aa.getType();
  string aaType2 = aa2.getType();

  for (unsigned int j = 0; j < aa.size(); j++)
	for (unsigned int k = 0; k < aa2.size(); k++)
	  if ( (aa2[k].getCode() == N) || (aa2[k].getCode() == CA) 
	       || (aa2[k].getCode() == C) || (aa2[k].getCode() == O) 
	       || (aa2[k].getCode() == CB) )
	    en += rapdf.calculateEnergy(aa[j], aa2[k], aaType, aaType2);
  
  return en;
}

/**
 * @Description calculates the Internal energy between two amino acids  
 * @param   amino acids references( AminoAcid&, AminoAcid&)
 * @return   corresponding value(double)
 */
double LoopModel::sICalcEnInt(AminoAcid& aa, AminoAcid& aa2){
  double en = 0.0;

  string aaType = aa.getType();
  string aaType2 = aa2.getType();

  for (unsigned int j = 0; j < aa.size(); j++)
    if ( (aa[j].getCode() == N) || (aa[j].getCode() == CA) 
	 || (aa[j].getCode() == C) || (aa[j].getCode() == O) 
	 || (aa[j].getCode() == CB) )
	for (unsigned int k = 0; k < aa2.size(); k++)
	  if ( (aa2[k].getCode() == N) || (aa2[k].getCode() == CA) 
	       || (aa2[k].getCode() == C) || (aa2[k].getCode() == O) 
	       || (aa2[k].getCode() == CB) )
	    en += rapdf.calculateEnergy(aa[j], aa2[k], aaType, aaType2);
  
  return en;
}


/**
 * @Description calculates the Internal energy between two amino acids  
 * @param   amino acids references( AminoAcid&, AminoAcid&)
 * @return   corresponding value(double)
 */
double LoopModel::sCalcEnInt(AminoAcid& aa, AminoAcid& aa2){
  double en = 0.0;

  string aaType = aa.getType();
  string aaType2 = aa2.getType();

  for (unsigned int j = 0; j < aa.size(); j++)
    if ( (aa[j].getCode() == N) || (aa[j].getCode() == CA) 
	 || (aa[j].getCode() == C) || (aa[j].getCode() == O) 
	 || (aa[j].getCode() == CB) )
	for (unsigned int k = 0; k < aa2.size(); k++)
	  if ( (aa2[k].getCode() == N) || (aa2[k].getCode() == CA) 
	       || (aa2[k].getCode() == C) || (aa2[k].getCode() == O) 
	       || (aa2[k].getCode() == CB) )
	    en += rapdf.calculateEnergy(aa[j], aa2[k], aaType, aaType2);
  
  return en;
}



/**
 * @Description calculates the energy  
 * @param   protein reference(spacer&), indexes for the Loop section(unsigned int, unsigned int))
 * @return   corresponding value(double)
 */
double LoopModel::calcOrigEnergy(const Spacer& sp, unsigned int index1, 
			unsigned int index2){
  long double en = 0.0;

  if (pInter)    {
      for (unsigned int i = 0; i < index1; i++)
	for (unsigned int j = index1+1; j <= index2; j++)
	  en += sICalcEn(const_cast<AminoAcid&>(sp.getAmino(i)),
				       const_cast<AminoAcid&>(sp.getAmino(j)));
      
      for (unsigned int i = index2+1; i < sp.sizeAmino(); i++)
	for (unsigned int j = index1+1; j <= index2; j++)
	  en += sICalcEn(const_cast<AminoAcid&>(sp.getAmino(i)),
				       const_cast<AminoAcid&>(sp.getAmino(j)));

      for (unsigned int i = index1+1; i < index2; i++)
	for (unsigned int j = i; j <= index2; j++)
	  en += sICalcEnInt(const_cast<AminoAcid&>(sp.getAmino(i)),
				       const_cast<AminoAcid&>(sp.getAmino(j)));
    }
  else    {
      for (unsigned int i = 0; i < index1; i++)
	for (unsigned int j = index1+1; j <= index2; j++)
	  en += sCalcEn(const_cast<AminoAcid&>(sp.getAmino(i)),
				       const_cast<AminoAcid&>(sp.getAmino(j)));
      
      for (unsigned int i = index2+1; i < sp.sizeAmino(); i++)
	for (unsigned int j = index1+1; j <= index2; j++)
	  en += sCalcEn(const_cast<AminoAcid&>(sp.getAmino(i)),
				       const_cast<AminoAcid&>(sp.getAmino(j)));

      for (unsigned int i = index1+1; i < index2; i++)
	for (unsigned int j = i; j <= index2; j++)
	  en += sCalcEnInt(const_cast<AminoAcid&>(sp.getAmino(i)),
				       const_cast<AminoAcid&>(sp.getAmino(j)));
    }

  return en;
}

