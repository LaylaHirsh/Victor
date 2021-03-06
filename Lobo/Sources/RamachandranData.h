/** 
 * @Class:              RamachandranData 
 * @Project Name:        
  * 
 * @Description:This class implements the container for ramachandran plot-like 
*    phi/psi angle combinations for the LoopTable class(es).
  *      
*/

#ifndef _RAMACHANDRANDATA_H_
#define _RAMACHANDRANDATA_H_

// Includes:
#include <vector>
#include <Debug.h>

namespace Biopool {
  
// Global constants, typedefs, etc. (to avoid):
/**
* @Description 
* @param  
* @return  
*/
inline double sqr(double x){
  return x*x;
}

    /**
 * @brief This struct implements the container for ramachandran plot-like 
*    phi/psi angle combinations for the LoopTable class(es).
     * */
struct RamachandranData
{
public:
// CONSTRUCTORS/DESTRUCTOR:
  RamachandranData();
  RamachandranData(const RamachandranData& orig);
  virtual ~RamachandranData();

// PREDICATES:
  double getRandomPhi(bool noAdvance = false);
  double getRandomPsi(bool noAdvance = false);
  static double getAngleTol() 
  { return (PHI_ANGLE_TOL + PSI_ANGLE_TOL) / 2.0; }

// MODIFIERS:
  void copy(const RamachandranData& orig);
  void load(istream& input);
  void save(ostream& output);
  void cluster(double cutoff);

  static void setAngleTol(float t) { PHI_ANGLE_TOL = t; PSI_ANGLE_TOL = t; }
  static void setPhiAngleTol(float t) { PHI_ANGLE_TOL = t; }
  static void setPsiAngleTol(float t) { PSI_ANGLE_TOL = t; }

// OPERATORS:
  RamachandranData& operator=(const RamachandranData& orig);
  
private:

// ATTRIBUTES:
  static float PHI_ANGLE_TOL;
  static float PSI_ANGLE_TOL;

  unsigned long nextRama;	
  vector<double> ramaPhi; 
  vector<double> ramaPsi; 

// HELPERS:

  double pGetRand();

};
 
 
// ---------------------------------------------------------------------------
//                               RamachandranData
// -----------------x-------------------x-------------------x-----------------
/**
* @Description returns a random value
* @param  none
* @return  the corresponding value( double)
*/
inline double RamachandranData::pGetRand(){
    double tmp = 0.0;
    for (unsigned int i = 0; i < 12; i++)
	tmp += static_cast<double>(rand())/RAND_MAX;

    return tmp - 6;
}

} // namespace

#endif //_RAMACHANDRANDATA_H_
