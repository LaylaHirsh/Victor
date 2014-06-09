/**
* @Class:              Atom
* @Base Class(es):     SimpleBond
* @Containing:         Identity, AtomCode
* @Project Name:       Victor
*/

#ifndef _ATOM_H_
#define _ATOM_H_

// Includes:
#include <AtomCode.h>
#include <vector3.h>
#include <matrix3.h>
#include <SimpleBond.h>

namespace Biopool {

    
// Global constants, typedefs, etc. (to avoid):

class Group;
/** @brief class implements a simple atom type.
* 
* @Description Includes methods that allow to get and set type, bind, unbind,coordinates , code,  etc.NB Angles are in degrees.
 **/
class Atom : public SimpleBond{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  Atom(unsigned int mI = 1, unsigned int mO = 4);
  Atom(const Atom& orig); 
  virtual ~Atom();

// PREDICATES:
  AtomCode getCode() const;
  unsigned long getNumber() const;

  vgVector3<double> getCoords();

  double getBFac() { return Bfac; }

  double distance(Atom& other);

  bool hasSuperior();
  Group& getSuperior();
  const Group& getSuperior() const;

  vgVector3<double> getTrans() const;
  vgMatrix3<double> getRot() const;
  virtual bool inSync(); // coords in-sync with changes?

// MODIFIERS:
  void clear();
  void copy(const Atom& orig);
  void bindStructure(Atom& before, bool connect); 
  // sets this relative to before if connect == true, otherwise like bindIn
  virtual void bindIn(SimpleBond& c);
  virtual void bindOut(SimpleBond& c);
  virtual void unbindIn(SimpleBond& c);
  virtual void unbindOut(SimpleBond& c);
  void setType(string _name);
  void setCode(AtomCode ac);
  void setNumber(unsigned long _number);

  void setCoords(double _x, double _y, double _z);
  void setCoords(vgVector3<double> c);

  void setBFac(double _b) { Bfac = _b; }

  void setTrans(vgVector3<double> t);
  void addTrans(vgVector3<double> t);
  void setRot(vgMatrix3<double> r);
  void addRot(vgMatrix3<double> r);
  virtual void sync(); // synchronize coords with structure

  virtual const Atom& getInBond(unsigned n) const;
  virtual Atom& getInBond(unsigned n);
  virtual const Atom& getOutBond(unsigned n) const;
  virtual Atom& getOutBond(unsigned n);

  void setSuperior(Group* gr);
  void setModified();
  void setUnModified();   //used in Qmean: this flag cause problem there.

// OPERATORS:
  Atom& operator=(const Atom& orig);

protected:

private:

// HELPERS: 
  bool isNotFirstAtomInStructure();
  void propagateRotation(); // pass on rotation matrix to following atoms


// ATTRIBUTES:
  Group* superior;           // structure to which atom belongs,
  AtomCode type;	         // Atom type
  vgVector3<double> coords;      // xyz-Coords

  double Bfac;                   // B-factor

  vgVector3<double> trans;       // relative translation
  vgMatrix3<double> rot;         // relative rotation
  bool modified;                 // --""--  modified?  
};

/** @example AtomTest.cc
 *  A simple program to test class Atom's features.
 */

// ---------------------------------------------------------------------------
//                                    Atom
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:
/**
 * @Description returns the corresponding atom type
 */
inline AtomCode Atom::getCode() const{ 
  return type; 
}
/**
 * @Description returns the corresponding atom id

 */
inline unsigned long Atom::getNumber() const { 
  return id; 
}
/**
 * @Description returns the corresponding coordinates of the atom in a vector
 */
inline vgVector3<double> Atom::getCoords(){
  if (!inSync())
    sync();
  return coords;
}

/**
 * @Description returns a reference to the group of the atom
 */
inline Group& Atom::getSuperior(){
  PRECOND( superior != NULL, exception);
  return *superior;
}
/**
 * @Description returns a reference to the group of the atom
 */
inline const Group& Atom::getSuperior() const{
  PRECOND( superior != NULL, exception);
  return *superior;
}
 /**
 * @Description velifies if the atom is part of a group
 */
inline bool Atom::hasSuperior(){
  return (superior != NULL);
}
/**
 * @Description returns the translation vector of the atom
 */
inline vgVector3<double> Atom::getTrans() const{
  return trans;
}
/**
 * @Description return the rotation matrix
 */
inline vgMatrix3<double> Atom::getRot() const{
  return rot;
}
/**
 * @Description verifies if its syncronized
 */
inline bool Atom::inSync(){
  return (!modified); 
}

// MODIFIERS:
 /**
 * @Description connects the atom to the given structure
 * @param reference to the previous atom(Atom&), flag to set the translation vector(bool)
 * @return changes are made internally(void)
 */
inline void  Atom::bindStructure(Atom& before, bool connect){
  if (connect)
    this->setTrans( this->getCoords() - before.getCoords() );
  this->bindIn(before);
}
/**
 * @Description bind in the atom 
 * @param reference to a simple bond(SimpleBond&)
 * @return changes are made internally(void)
 */
inline void Atom::bindIn(SimpleBond& c){
  SimpleBond::bindIn(c);
  setModified();  
}
/**
 * @Description bind out the atom 
 * @param reference to a simple bond(SimpleBond&)
 * @return changes are made internally(void)
 */
inline void Atom::bindOut(SimpleBond& c){
  SimpleBond::bindOut(c);
  setModified();  
}
/**
 * @Description unbind in the atom 
 * @param reference to a simple bond(SimpleBond&)
 * @return changes are made internally(void)
 */
inline void Atom::unbindIn(SimpleBond& c){
  SimpleBond::unbindIn(c);
  setModified();  
}
/**
 * @Description unbind out the atom 
 * @param reference to a simple bond(SimpleBond&)
 * @return changes are made internally(void)
 */
inline void Atom::unbindOut(SimpleBond& c){
  SimpleBond::unbindOut(c);
  setModified();  
}
/**
 * @Description sets the atom type and id
 * @param atom name(string)
 * @return changes are made internally(void)
 */
inline void Atom::setType(string _name) { 
  id.setName(_name); 
  type = AtomTranslator(_name); 
}
/**
 * @Description sets the atom type and id
 * @param atom type(AtomCode)
 * @return changes are made internally(void)
 */
inline void Atom::setCode(AtomCode ac) { 
  type = ac;
  id.setName(AtomTranslator(ac)); 
}
/**
 * @Description sets the atom   id
 * @param atom number(unsigned long)
 * @return changes are made internally(void)
 */
inline void Atom::setNumber(unsigned long _number) { 
  id.setNumber(_number); 
}
/**
 * @Description sets the atom translation vector
 * @param a vector containing the values for the translation vector(vgVector3<double>)
 * @return changes are made internally(void)
 */
inline void Atom::setTrans(vgVector3<double> t){
  trans = t;
  setModified();
}
/**
 * @Description add some values to the translation vector
 * @param values to be add(vgVector3<double> )
 * @return changes are made internally(void)
 */
inline void Atom::addTrans(vgVector3<double> t){
  trans += t;
  setModified();
}
/**
 * @Description defines the rotation matrix
 * @param values to set(vgMatrix3<double>)
 * @return changes are made internally(void)
 */
inline void Atom::setRot(vgMatrix3<double> r){
  rot = r;
  setModified();
}
/**
 * @Description add some values to the rotation matrix
 * @param values to be add(vgMatrix3<double> )
 * @return changes are made internally(void)
 */
inline void Atom::addRot(vgMatrix3<double> r){
  rot = r * rot;
  setModified();
}
/**
 * @Description obtains the inBond 
 * @param index (unsigned)
 * @return reference to the atom(const Atom&)
 */
inline const Atom& Atom::getInBond(unsigned n) const{
  return dynamic_cast<const Atom&>(SimpleBond::getInBond(n));
}
/**
 * @Description obtains the inBond 
 * @param index (unsigned)
 * @return reference to the atom( Atom&)
 */
inline Atom& Atom::getInBond(unsigned n){
  return dynamic_cast<Atom&>(SimpleBond::getInBond(n));
}
/**
 * @Description obtains the outBond 
 * @param index (unsigned)
 * @return reference to the atom(const Atom&)
 */
inline const Atom& Atom::getOutBond(unsigned n) const{
  return dynamic_cast<const Atom&>(SimpleBond::getOutBond(n));
}
/**
 * @Description obtains the outBond 
 * @param index (unsigned)
 * @return reference to the atom( Atom&)
 */
inline Atom& Atom::getOutBond(unsigned n){
  return dynamic_cast<Atom&>(SimpleBond::getOutBond(n));
}

/**
 * @Description sets the group for the atom
 * @param pointer to the group(group*)
 * @return changes are made internally(void)
 */
inline void Atom::setSuperior(Group* gr){
  this->superior = gr;
}

// OPERATORS:


} // namespace
#endif //_ATOM_H_
