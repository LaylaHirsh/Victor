/**
 * @Class:           SimpleBond
 * @Derived classes: Bond, Atom
 * @Author:          Silvio Tosatto
 * @Description:     
*                  Attention: copy() strips orig from its SimpleBonds and 
*                  attaches them to the new SimpleBond. 
* 
*/
#ifndef _SIMPLEBOND_H_
#define _SIMPLEBOND_H_

// Includes:
#include <algorithm>
#include <vector>
#include <string>
//#include <Debug.h>
#include <Identity.h>
using namespace std;

// Global constants, typedefs, etc. (to avoid):
/** @brief Defines chemical and abstract bonds between objects.
*                  eg.: covalent bonds.
 * */


class SimpleBond{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  SimpleBond(unsigned int mI = 1, unsigned int mO = 4);
  SimpleBond(const SimpleBond& orig); 
  virtual ~SimpleBond();

// PREDICATES:
  virtual string getType() const;

  bool isBond(const SimpleBond& c) const; 
  bool isInBond(const SimpleBond& c) const;
  bool isOutBond(const SimpleBond& c) const;
  bool isIndirectBond(const SimpleBond& c) const; 
  // eg. A to C if A bond B & B bond C
  bool isIndirectInBond(const SimpleBond& c) const;
  bool isIndirectOutBond(const SimpleBond& c) const;
  bool isTorsionBond(const SimpleBond& c) const; 
  // eg. A to D if A indirect bond C & C bond D

  virtual const SimpleBond& getInBond(unsigned int n) const;
  virtual const SimpleBond& getOutBond(unsigned int n) const;
  virtual SimpleBond& getInBond(unsigned int n);
  virtual SimpleBond& getOutBond(unsigned int n);

  unsigned int sizeInBonds() const;
  unsigned int sizeOutBonds() const;

  unsigned int getMaxInBonds() const;
  unsigned int getMaxOutBonds() const;

// MODIFIERS:
  virtual void copy(const SimpleBond& orig);

  virtual void setType(string _name);

  virtual void bindIn(SimpleBond& c);
  virtual void bindOut(SimpleBond& c);

  virtual void unbindIn(SimpleBond& c);
  virtual void unbindOut(SimpleBond& c);

  void setMaxInBonds(unsigned int m);
  void setMaxOutBonds(unsigned int m);

// OPERATORS:
  bool operator==(const SimpleBond& other) const;
  SimpleBond& operator=(const SimpleBond& orig);
  
protected:
  
// HELPERS:
  virtual void pUnbindIn(SimpleBond& c);
  virtual void pUnbindOut(SimpleBond& c);

// ATTRIBUTES:
  vector<SimpleBond*> inBonds, outBonds; // current in and out bonds
  unsigned int  maxIn, maxOut;           // maximum number of in and out bonds
  Identity id;                           // Object Id and name

private:
  
};


// ---------------------------------------------------------------------------
//                                    SimpleBond
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:


// -*- C++ -*-----------------------------------------------------------------
//
//  Description:
//    Returns type (eg. C-Atom, GLY amino acid, etc.).  
//
// ----------------------------------------------------------------------------
inline string SimpleBond::getType() const{ 
  return this->id; 
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::isBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Is c bonded to this ?
//
// ----------------------------------------------------------------------------
inline bool SimpleBond::isBond(const SimpleBond& c) const{
  return (isInBond(c) || isOutBond(c));
} 


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::isInBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Is c in-bonded to this ?
//
// ----------------------------------------------------------------------------
inline bool SimpleBond::isInBond(const SimpleBond& c) const{
  return find(inBonds.begin(), inBonds.end(), &c) != inBonds.end();
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::isOutBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Is c out-bonded to this ?
//
// ----------------------------------------------------------------------------
inline bool
SimpleBond::isOutBond(const SimpleBond& c) const
{
  return find(outBonds.begin(), outBonds.end(), &c) != outBonds.end();
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::sizeInBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    How many in-bonds are there ?
//
// ----------------------------------------------------------------------------
inline unsigned int  SimpleBond::sizeInBonds() const{
  return inBonds.size();
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::sizeOutBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    How many out-bonds are there ?
//
// ----------------------------------------------------------------------------
inline unsigned int SimpleBond::sizeOutBonds() const{
  return outBonds.size();
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::getInBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Returns i-th in-bond.
//
// ----------------------------------------------------------------------------
inline SimpleBond& SimpleBond::getInBond(unsigned int n) {
  PRECOND(n < sizeInBonds(), exception);
  return *inBonds[n];
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::getOutBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Returns i-th out-bond.
//
// ----------------------------------------------------------------------------
inline SimpleBond& SimpleBond::getOutBond(unsigned int n) {
  PRECOND(n < sizeOutBonds(), exception);
  return *outBonds[n];
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::getInBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Returns i-th in-bond.
//
// ----------------------------------------------------------------------------
inline const SimpleBond& SimpleBond::getInBond(unsigned int n) const {
  PRECOND(n < sizeInBonds(), exception);
  return *inBonds[n];
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::getOutBond
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Returns i-th out-bond.
//
// ----------------------------------------------------------------------------
inline const SimpleBond& SimpleBond::getOutBond(unsigned int n) const {
  PRECOND(n < sizeOutBonds(), exception);
  return *outBonds[n];
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::getMaxInBonds
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Get maximum of in-bonds.
//
// ----------------------------------------------------------------------------
inline unsigned int SimpleBond::getMaxInBonds() const{
  return maxIn;
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::getMaxOutBonds
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Get maximum of out-bonds.
//
// ----------------------------------------------------------------------------
inline unsigned int SimpleBond::getMaxOutBonds() const{
  return maxOut;
}

// MODIFIERS:


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::setType
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Sets type (eg. C atom, GLY amino acid, etc.) of this.
//
// ----------------------------------------------------------------------------
inline void SimpleBond::setType(string _name) { 
  id.setName(_name); 
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::setMaxInBonds
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Sets maximum of in-bonds.
//
// ----------------------------------------------------------------------------
inline void SimpleBond::setMaxInBonds(unsigned int m){
  PRECOND((m >= inBonds.size()), exception);
  maxIn = m;
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SimpleBond::setMaxOutBonds
//
//  Author:        Silvio Tosatto 
//
//  Date:          08/99
//
//  Description:
//    Sets maximum of out-bonds.
//
// ----------------------------------------------------------------------------
inline void SimpleBond::setMaxOutBonds(unsigned int m){
  PRECOND((m >= outBonds.size()), exception);
  maxOut = m;
}

// OPERATORS:

inline SimpleBond& SimpleBond::operator=(const SimpleBond& orig){
  PRINT_NAME;
  if (&orig != this)
    copy(orig);
  return *this;
}

inline bool SimpleBond::operator==(const SimpleBond& other) const{
  return id == (other.id);
}

#endif //_SimpleBond_H_

