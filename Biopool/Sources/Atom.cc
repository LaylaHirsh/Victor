/**
* @Class:             Atom
* @Base Class(es):    SimpleBond
* @Project Name:      Victor

*/

// Includes:
#include <Atom.h>
#include <Group.h>
#include <matrix3.h>
#include <math.h>
#include <IoTools.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:
/**
 * @description Constructor 
 */
Atom::Atom(unsigned int mI, unsigned int mO) : SimpleBond(mI,mO), 
  superior(NULL), type(X), coords(0,0,0), Bfac(0.0), trans(0,0,0), rot(1), 
  modified(false){ 
  PRINT_NAME;
}
/**
 * @description Constructor based in another object
 */
Atom::Atom(const Atom& orig){
  PRINT_NAME;
  this->copy(orig);
}
/**
 * @description Destructor 
 */
Atom::~Atom(){
  PRINT_NAME;
}

// PREDICATES:

/**
 * @description finds the distance between two atoms
 * @param reference to an atom(Atom&)
 * @return distance value(double)
 */
double Atom::distance(Atom& other){
  sync();
  other.sync();
  
  double tmp = sqrt( 
		    (coords[0]-other.coords[0]) * (coords[0] -other.coords[0])
		    + (coords[1]-other.coords[1]) * (coords[1]-other.coords[1])
		    + (coords[2]-other.coords[2]) * (coords[2]-other.coords[2])
		    );
  return tmp;
}


// MODIFIERS:
/**
 * @description Defines the atom coords 
 * @param values for the x,y,z coords(double,double,double)
 * @return changes are made internally(void)
 */
void Atom::setCoords(double _x, double _y, double _z){
  vgVector3<double> tmp(_x, _y, _z);
  setCoords(tmp);
}
/**
 * @description Defines the atom coords 
 * @param a vector containing the values for the x,y,z coords(vgVector3<double>)
 * @return changes are made internally(void)
 */
void Atom::setCoords(vgVector3<double> c){
  sync();
  vgVector3<double> oldTrans = trans;
  coords = c;
  
  if (isNotFirstAtomInStructure())
    trans = coords - getInBond(0).coords;
  else
    if (hasSuperior())
      trans = coords - getSuperior().getTrans();
    else 
      trans = coords;

  // adjust coords of following atoms:
  for (unsigned int i = 0; i < sizeOutBonds(); i++)
    {
      getOutBond(i).setTrans( getOutBond(i).getTrans() + (oldTrans - trans));
    }

  setModified();
}
/**
 * @description Copies an atom into another
 * @param reference to the atom(atom&)
 * @return changes are made internally(void)
 */
void Atom::copy(const Atom& orig){
  PRINT_NAME; 
  SimpleBond::copy(orig);

  setNumber(orig.getNumber());

  superior = orig.superior;
  type = orig.type;
  coords = orig.coords;
  Bfac = orig.Bfac;

  trans = orig.trans;
  rot = orig.rot;
  modified = orig.modified;
}
/**
 * @description synchronize coords with structure
 * @param none
 * @return changes are made internally(void)
 */
void Atom::sync(){ // 
  if (inSync())
    return;

  vgMatrix3<double> tmpMatrix(1);    
  vgVector3<double> supTrans(0,0,0);

  if (isNotFirstAtomInStructure())
    getInBond(0).sync();
  else 
    if (hasSuperior())   {    // use superior's trans & rot
	supTrans = getSuperior().getTrans();
	rot = getSuperior().getRot() * rot;
	const_cast<Group&>(getSuperior()).setRot(tmpMatrix);
      }
  
  trans = rot * trans;
  coords = trans + supTrans;         // set the relative position

  if (isNotFirstAtomInStructure())
    coords += getInBond(0).getCoords();  // make absolute position

  propagateRotation();
  modified = false;
}
/**
 * @description sets the modified flag for the atom, and for all the out bonds
 * @param none
 * @return changes are made internally(void)
 */
void Atom::setModified(){
  if (modified)
    return;
  modified = true;

  for (unsigned int i = 0; i < sizeOutBonds(); i++)
    getOutBond(i).setModified();

}
/**
 * @description sets the unmodified flag for the atom 
 * @param none
 * @return changes are made internally(void)
 */
void Atom::setUnModified(){
  modified = false;
}


// OPERATORS:
/**
 * @description copies one atom into another 
 * @param reference to the atom to copy(const Atom&)
 * @return reference to the new atom(Atom&)
 */
Atom& Atom::operator=(const Atom& orig){
  PRINT_NAME;

  if (&orig != this)
    copy(orig);
  return *this;
}


// HELPERS:
/**
 * @description sets the rotation matrix data
 * @param none
 * @return changes are made internally(void)
 */
void Atom::propagateRotation(){
  vgMatrix3<double> tmpMatrix(1);    
  if (rot == tmpMatrix)  
    return;
  
  for (unsigned int i = 0; i < sizeOutBonds(); i++)  
    if ( !hasSuperior() || 
	 ( !( (getSuperior().getType() == "PRO") && (getCode() == CD)
	      && (getOutBond(i).getCode() == N) ) 
	   && !( (getSuperior().getType() == "PHE") && (getCode() == CE2)
		 && (getOutBond(i).getCode() == CZ) )
	     && !( (getSuperior().getType() == "TYR") && (getCode() == CE2)
		   && (getOutBond(i).getCode() == CZ) )
	   && !( (getSuperior().getType() == "TRP") && (getCode() == NE1)
		 && (getOutBond(i).getCode() == CE2) )
	   && !( (getSuperior().getType() == "TRP") && (getCode() == CZ3)
		 && (getOutBond(i).getCode() == CH2) )
	   && !( (getSuperior().getType() == "HIS"  && (getCode() == CE1)
		  && (getOutBond(i).getCode() == NE2)) )
	   )
	 )
      {
	getOutBond(i).addRot(rot);
      }
  // Cyclic structures (PRO, PHE, TYR, TRP & HIS) are special cases 
  // when propagating the rotation, because the the re-conjunction of 
  // two branches the rotation would be added twice.
  rot = tmpMatrix;
}

 /**
 * @description verifies if its not the first atom in the structure
 * @param none
 * @return true if its not the first atom in structure
 */
inline bool Atom::isNotFirstAtomInStructure(){
  if (sizeInBonds() && !((getSuperior().getType() == "PRO") 
	    && (getCode() == N) && (getInBond(0).getCode() == CD)))
    // NB: Proline's N is a special case because unbound PRO 
    // (ie. not part of a Spacer) has CD as the 1st inBond of N,
    // which would crash the program
    return true;
  else 
    return false;
}
