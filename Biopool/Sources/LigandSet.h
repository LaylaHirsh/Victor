/** @Class            LigandSet
* @Author:            Silvio Tosatto, previous work: Marcus Pruemmer
* @Project Name:      Victor
* @ATTENTION:         This class is *NOT* finished yet.
*                     There are many methods which have only a head
*                     declaration and no code inside.
*                     2012-Francesco Lovo: insert, delete, getLigand, load and save
**                    methods are completed
*/
#ifndef _LigandSet_H_
#define _LigandSet_H_


// Includes:
#include <Bond.h>
#include <Polymer.h>
#include <Ligand.h>
#include <Visitor.h>
#include <AminoAcidCode.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief class implements methods to manage the ligandSet
 * 
* @Description Includes methods that allow to set, get and modify ligand set content.
* @This 
 * */
class LigandSet : public Polymer
{
public: 
 
  
  // CONSTRUCTORS/DESTRUCTOR:
  LigandSet();
  LigandSet(const LigandSet& orig);
  virtual ~LigandSet(); 
  
  // PREDICATES:
  int getStartOffset() { return startOffset; }
  virtual string getClassName() const
    { return "LigandSet"; }
  Ligand& getLigand(unsigned int n);
  const Ligand& getLigand(unsigned int n) const;
  bool isGap(int index);
  const unsigned int  sizeLigand() const;
  unsigned int maxPdbNumber() { 
    return startOffset + gaps.size() + sizeLigand(); }

  // MODIFIERS:
  void addGap(int index);
  void insertComponent(Component* g);
  void removeLigand(Ligand* g);
  void deleteLigand(Ligand* g);
  void setStartOffset(int _offset) {startOffset = _offset;} 
  
  void copy(const LigandSet& orig);

  inline void setTrans(vgVector3<double> t);
  inline void addTrans(vgVector3<double> t);
  inline void setRot(vgMatrix3<double> r);
  inline void addRot(vgMatrix3<double> r); 

  void load(Loader& l);  // data loader
  void save(Saver&  s);  // data saver

  virtual LigandSet* clone();
  
  // OPERATORS:
  LigandSet& operator=(const LigandSet& orig);
  Ligand&  operator[](unsigned int n);
  const Ligand&  operator[](unsigned int n) const;
  
protected:

  // HELPERS:
  void resetBoundaries();
  // ATTRIBUTES
  int startOffset;
  vector<int> gaps;       //to keep trace of gaps in the pdb file
  
private:
/** @example LigandSetTest.cc
   *  A simple program to test class LigantSet's features.
 */
};

// ---------------------------------------------------------------------------
//                                    LigandSet
// -----------------x-------------------x-------------------x-----------------

// MODIFIERS

inline void
LigandSet::setTrans(vgVector3<double> t)
{
  if ( sizeLigand() )
    getLigand(0).setTrans(t);
}

inline void
LigandSet::addTrans(vgVector3<double> t)
{
  if ( sizeLigand() )
    getLigand(0).addTrans(t);
}

inline void 
LigandSet::setRot(vgMatrix3<double> r)
{
  if (sizeLigand())
    getLigand(0).setRot(r);
}

inline void 
LigandSet::addRot(vgMatrix3<double> r)
{
  if (sizeLigand())
    getLigand(0).addRot(r);
}

inline void
LigandSet::load(Loader& l)
{
  l.loadLigandSet(*this);
  resetBoundaries();
}

inline void 
LigandSet::save(Saver& s)
{
  s.saveLigandSet(*this);
}

inline const unsigned int
LigandSet::sizeLigand() const
{
  return components.size();
}

inline Ligand&
LigandSet::getLigand(unsigned int n)
{
  // ERROR("Not implemented yet!",exception);
  if ( n > components.size() - 1 )
    ERROR("Index out of range",exception);
  return *(dynamic_cast<Ligand*> (components[n]));
}

inline const Ligand&
LigandSet::getLigand(unsigned int n) const
{
  // ERROR("Not implemented yet!",exception);
  if ( n > components.size() - 1 )
    ERROR("Index out of range",exception);
  return *( dynamic_cast<const Ligand*> (components[n]) );
}

inline void
LigandSet::removeLigand(Ligand* l)
{
   Polymer::removeComponent(l);
  setModified();
}

inline void
LigandSet::deleteLigand(Ligand* l)
{
   Polymer::deleteComponent(l);
  setModified();
}


// OPERATORS:

inline Ligand&
LigandSet::operator[](unsigned int n)
{
  return getLigand(n);
}

inline const Ligand&
LigandSet::operator[](unsigned int n) const
{
  return getLigand(n);
}

// HELPERS:

} // namespace
#endif //_LigandSet_H_







