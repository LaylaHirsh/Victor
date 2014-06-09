/**
* @Class:             Nucleotide
* @Author:            Damiano Piovesan
* @Project Name:      Victor
*/

// Includes:
#include <Nucleotide.h>
#include <NucleotideCode.h>
#include <Debug.h>
#include <limits.h>
#include <float.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Biopool;

/**
 * @Description Constructor
* @param none
*/
Nucleotide::Nucleotide() : Group(1,1), type(XX), icc() { }


/**
 * @Description Constructor
* @param Nucleotide
*/
Nucleotide::Nucleotide(const Nucleotide& orig)
{
  PRINT_NAME;
  this->copy(orig);
}
/**
 * @Description DESTRUCTOR
 * @param none
 */
Nucleotide::~Nucleotide()
{ PRINT_NAME; } 
/**
 * @Description Returns the atom corresponding to N, 
 * aminoacids have only a single possible, hard coded, open in-bond
 * @param unsigned int 
 * @return const Atom& 
 */


// MODIFIERS:

/**
* @Description  Copies an aa
 * @param const Nucleotide& (copy from the orig)
 * @return  void
 */
void Nucleotide::copy(const Nucleotide& orig)
{
  PRINT_NAME; 
  Group::copy(orig);

  type = orig.type;
  
  // set absolute position to orig's:
  if (orig[0].sizeInBonds())
    {
      setTrans(const_cast<Nucleotide&>(orig)[0].getInBond(0).getCoords());
    }
}

/**
* @Description  Clone the aa
 * @param none
 * @return  Component* 
 */
Component* Nucleotide::clone(){
  Nucleotide* tmp = new Nucleotide;
  tmp->copy(*this);
  return tmp;
}

// OPERATORS:
/**
* @Description  Operator =, assign the aa
 * @param Nucleotide reference
 * @return  Nucleotide
 */
Nucleotide& Nucleotide::operator=(const Nucleotide& orig){
  PRINT_NAME;
  if (&orig != this)
    copy(orig);
  return *this;
}


