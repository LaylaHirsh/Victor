/**
 
*                  Attention: copy() strips orig from its bonds and 
*                  attaches them to the new bond. 
*/

#ifndef _BOND_H_
#define _BOND_H_
 
// Includes:
#include <algorithm>
#include <vector>
#include <SimpleBond.h>
#include <Atom.h>

namespace Biopool {
    class Atom;
/** @brief Defines chemical and abstract bonds between objects which are compositions of atoms.
 * 
* @Description Eg. 'bond' between two amino acids.
 * */
// Global constants, typedefs, etc. (to avoid):



class Bond : public SimpleBond
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  Bond(unsigned int mI = 1, unsigned int mO = 1);
  Bond(const Bond& orig); 
  virtual ~Bond();

  // PREDICATES:
  virtual const Atom& getInBondRef(unsigned int n) const;
  virtual Atom& getInBondRef(unsigned int n);
  virtual const Atom& getOutBondRef(unsigned int n) const;
  virtual Atom& getOutBondRef(unsigned int n);

  virtual const Atom& getOpenInBondRef(unsigned int n = 0) const;
  virtual Atom& getOpenInBondRef(unsigned int n = 0);
  virtual const Atom& getOpenOutBondRef(unsigned int n = 0) const;
  virtual Atom& getOpenOutBondRef(unsigned int n = 0);

  unsigned int sizeOpenInBonds() const;
  unsigned int sizeOpenOutBonds() const;

  // MODIFIERS:
  void copy(const Bond& orig);

  virtual void bindIn(Atom& _this, Bond& c, Atom& _other);
  virtual void bindOut(Atom& _this, Bond& c, Atom& _other);

  virtual void unbindIn(Bond& c);
  virtual void unbindOut(Bond& c);

  // OPERATORS:
  
protected:
  
private:
  
  // HELPERS:
  virtual void pUnbindIn(Bond& c, bool unbind = false);
  virtual void pUnbindOut(Bond& c, bool unbind = false);

  // ATTRIBUTES:
  vector<Atom*> inRef, outRef;  
  // reference to the atom containg the in- and out-bonds

};


// ---------------------------------------------------------------------------
//                                    Bond
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:

/**
 * @Description  Returns reference to the i-th in-bond atom.
 * @param unsigned int 
 * @return  Atom reference
 */
 inline Atom&
Bond::getInBondRef(unsigned int n) 
{
  PRECOND(n < inRef.size(), exception);
  return *inRef[n];
}

inline const Atom&
Bond::getInBondRef(unsigned int n) const
{
  PRECOND(n < inRef.size(), exception);
  return *inRef[n];
}

/**
 * @Description  Returns reference to the i-th in-bond atom.
 * @param unsigned int 
 * @return  Atom reference
 */ 
 inline Atom&
Bond::getOutBondRef(unsigned int n) 
{
  PRECOND(n < outRef.size(), exception);
  return *outRef[n];
}

inline const Atom&
Bond::getOutBondRef(unsigned int n) const
{
  PRECOND(n < outRef.size(), exception);
  return *outRef[n];
}



/**
 * @Description  Returns the number of open in-bonds, ie. how many in-bonds can still be added.
 * @param void
 * @return unsigned int 
 */ 
 inline unsigned int 
Bond::sizeOpenInBonds() const
{
  return (getMaxInBonds() - sizeInBonds());
}


/**
 * @Description  Returns the number of open out-bonds, ie. how many out-bonds can still  be added.
 * @param void
 * @return unsigned int 
 */  
 inline unsigned int Bond::sizeOpenOutBonds() const
{
  return (getMaxOutBonds() - sizeOutBonds());
}


// MODIFIERS:

// OPERATORS:

} // namespace
#endif //_BOND_H_


