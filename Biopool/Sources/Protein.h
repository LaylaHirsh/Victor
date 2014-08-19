/**
 * @Class            Protein
 * @Base Class(es):    Polymer
 * @Containing:        LigandSet, Spacer
* @Project Name:      Victor
 *  Warning:           Copy and delete methods needs more work
 */

#ifndef PROTEIN_H
#define	PROTEIN_H

// Includes:
#include <Component.h>
#include <Polymer.h>
#include <Spacer.h>
#include <LigandSet.h>

namespace Biopool{
    /** @brief This class is a container of Polymers
 * 
* @Description Each one storing a Spacer and (eventually) a LigandSet for each
 *                     valid chain in the PDB input file.
 * */
class Protein : public Polymer{
public: 
  
  // CONSTRUCTORS/DESTRUCTOR:
  Protein();                                                    
  Protein(const Protein& orig);                                 
  virtual ~Protein();                                          
  
  // PREDICATES:
  virtual string getClassName() const                           
    { return "Protein"; }
  Polymer& getPolymer(unsigned int n);
  Spacer* getSpacer(char c);     //return a pointer to the Spacer with the correct chainID
  Spacer* getSpacer(unsigned int n);            
  
  LigandSet* getLigandSet(char c);//return a pointer to the LigandSet with the correct chainID(NULL if not exist)
  LigandSet* getLigandSet(unsigned int n);     
  
  const unsigned int sizeProtein() const;                                  
  unsigned int getChainNum (char c);                                     
  char getChainLetter (unsigned int i);                                  
  vector <char> getAllChains();                                          
  
  void save(Saver&  s);  // data saver                  
  
  // MODIFIERS:
  void addChain(char c);                                       
  void insertComponent(Component* c); 
  void setChainSelection();                                    
  void removeComponent(Component* c);
  void deleteComponent(Component* c);
  
  void copy(const Protein& orig);                       
  void load(Loader& l);  // data loader                 
  
  virtual Protein* clone();                            
  
  // OPERATORS:
  Protein& operator=(const Protein& orig);            
  Polymer& operator[](unsigned int n);                 
  //const Polymer& operator[](unsigned int n) const;     
  
protected:
  // HELPERS:
  void printComponents();
  // ATTRIBUTES 
private:
  vector<char> chains;
    //true: user selected a particular chain or all chains
                        //false: first chain selected by default
};  
// ---------------------------------------------------------------------------
//                                    Protein
// -----------------x-------------------x-------------------x-----------------

//PREDICATES
inline const unsigned int Protein::sizeProtein() const{
  return components.size();
}

inline vector<char>Protein::getAllChains()  {
    return chains;
}
inline void   Protein::save(Saver& s) {
  s.saveProtein(*this);
}

// MODIFIERS
inline void Protein::addChain(char c){
    chains.push_back(c);
}

inline void  Protein::load(Loader& l){
  l.loadProtein(*this);
}

// OPERATORS:
inline Polymer&
Protein::operator[](unsigned int n)                     
{
    return getPolymer(n); 
}



// HELPERS:
/** @example ProteinTest.cc
   *  A simple program to test class Protein's features.
 */
}//namespace
#endif	/* PROTEIN2_H */

