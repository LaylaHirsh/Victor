 
#ifndef _STRUCTURALALIGNMENT_H_
#define _STRUCTURALALIGNMENT_H_

// Includes:
#include <string>
#include <vector>
#include <Debug.h>
#include <AlignmentBase.h>
#include <Spacer.h>
#include <vector3.h>
#include <matrix3.h>

using namespace Biopool;
using namespace std;

namespace Biopool 
{
  /** @brief   Class for structural alignments.
 * 
* @Description  
* @This 
 **/ 
struct EData 
{
    EData(unsigned int _o = 0, double _d = 0.0) : other(_o), dist(_d) {}
    unsigned int other;
    double dist;
};

// -----------------x-------------------x-------------------x-----------------
/// @b NB: This code is being discontinued.
///
/// @author Silvio Tosatto
/// @date 03/2002
// -----------------x-------------------x-------------------x-----------------
struct FData 
{
    FData(unsigned int _mi = 0, unsigned int _ma = 0, double _d = 0.0) : 
	min(_mi), max(_ma), dist(_d) { }
    unsigned int min;
    unsigned int max;
    double dist;
};

// -----------------x-------------------x-------------------x-----------------
/// Class for structural alignments.
/// @b NB: This code is being discontinued.
///
/// @author Silvio Tosatto
/// @date 03/2002
// -----------------x-------------------x-------------------x-----------------
class StructuralAlignment : public AlignmentBase
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  StructuralAlignment();
  StructuralAlignment(const StructuralAlignment& orig); 
  virtual ~StructuralAlignment();

// PREDICATES:
    Spacer& getTarget() { return spTarget; }
    Spacer& getTemplate() { return spTemplate; }

// MODIFIERS:

  void setTarget(Spacer& sp) { spTarget = sp; }
  void setTemplate(Spacer& sp) { spTemplate = sp; }

  void loadCE(istream& input, Spacer& spNew);
  void buildEquivalenceNetwork(); 
  void buildFragmentNetwork(double maxDist = 4.0);
       // maximum CA distance for "equivalence"
  void writeData(double maxDist = 4.0);
  virtual void copy(const StructuralAlignment& orig);

// OPERATORS:

protected:

private:

// HELPERS: 

  void pExecStructAli(vgMatrix3<double> rot, vgVector3<double> trans);  

// ATTRIBUTES:
  Spacer  spTarget;
  Spacer  spTemplate;

public:
  vector<EData> equivData;
  vector<FData> fragData;

};

} // namespace

#endif //_STRUCTURALALIGNMENT_H_
