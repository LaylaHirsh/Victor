/**
* @Class:             Ligand
* @Author:            Silvio Tosatto, previous work: Marcus Pruemmer
* @Project Name:      Victor
* @ATTENTION:         This class is *NOT* finished yet.
*                     There are many methods which have only a head
*                     declaration and no code inside.
*                     2012-Francesco Lovo: load and save methods are completed
*/

#ifndef _Ligand_H_
#define _Ligand_H_


// Includes:
#include <Bond.h>
#include <Polymer.h>
#include <Group.h>
#include <Visitor.h>
#include <AminoAcidCode.h>
#include <NucleotideCode.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief class implements methods to verify the ligand qualities
 * 
* @Description  
* @This 
 * */
class Ligand : public Group
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  Ligand();
  Ligand(const Ligand& orig);
  virtual ~Ligand(); 
  virtual string getClassName() const
    { return "Ligand"; }

  // PREDICATES:
  bool isMetalCompound();  //I don't understans which was the purpose of thoose
  bool isCommonMetal();    //two functions (their lists of IDs appears uncorrect)
 
  bool isSimpleMetalIon(); //so I introduced thoose to distinguish between 
  bool isWater();          //water, metal ions and other cofactors
  bool isCofactor();
  virtual void save(Saver& s);  // data saver
 
  
  // MODIFIERS:
  void copy(const Ligand& orig);

  virtual void load(Loader& l);  // data loader

  // OPERATORS:
  Ligand& operator=(const Ligand& orig);
  virtual Atom& operator[](unsigned int n);
  virtual const Atom& operator[](unsigned int n) const;
  

protected:
  // HELPERS:

  // ATTRIBUTES:
  
private:
};

// ---------------------------------------------------------------------------
//                                    Ligand
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:


inline void 
Ligand::save(Saver& s)
{
  s.saveLigand(*this);
}

// MODIFIERS:

inline void 
Ligand::load(Loader& l)
{
  l.loadLigand(*this);
  resetBoundaries();
}


// OPERATORS:

inline Atom& 
Ligand::operator[](unsigned int n)
{
    return Group::operator[](n);
}

inline const Atom& 
Ligand::operator[](unsigned int n) const
{
    return Group::operator[](n);
}


} // namespace
#endif //_Ligand_H_

