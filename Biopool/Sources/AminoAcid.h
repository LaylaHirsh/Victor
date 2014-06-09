/**

* @Class:              AminoAcid
* @Base Class(es):     Group
* @Containing:         SideChain
* @Author:             Silvio Tosatto
* @Project Name:       Victor
* 
*/

#ifndef _AMINOACID_H_
#define _AMINOACID_H_



// Includes:
#include <Group.h>
#include <SideChain.h>
#include <AminoAcidCode.h>
#include <IntCoordConverter.h>
#include <Visitor.h>


// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief class implements a simple amino acid.
 * 
* @Description Includes methods that allow to manages an amino acid: get and set angles, connections, bonds, side chain info, etc .
* @This NB Angles are in degrees.
 * */
class AminoAcid : public Group{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  AminoAcid();
  AminoAcid(const AminoAcid& orig);
  virtual ~AminoAcid(); 
  
  // PREDICATES:
  virtual string getClassName() const
    { return "AminoAcid"; }
  virtual unsigned int getCode() const;
  char getType1L();
  unsigned int size() const;
  unsigned int sizeBackbone() const;
 
  double getPhi(bool override = false);
  double getPsi(bool override = false);
  double getOmega(bool override = false);
  double getChi(unsigned int n);
  vector<double> getChi();
  unsigned int getMaxChi();

  StateCode getState();  // returns AA state (H, S, T, C)
  
  SideChain& getSideChain();
  const SideChain& getSideChain() const;

  bool isMember(const AtomCode& ac) const;

  void save(Saver& s);   // data saver

  // MODIFIERS:
  virtual void connectIn(AminoAcid* c, unsigned int offset = 1);
  virtual void connectOut(AminoAcid* c, unsigned int offset = 1);
  virtual AminoAcid* unconnectIn();
  virtual AminoAcid* unconnectOut();

  void copy(const AminoAcid& orig);
  void setType1L(char _name);
  void setType(string _name);

  void setPhi(double a);
  void setPsi(double a);
  void setOmega(double a);
  void setChi(unsigned int n, double a);
  void setChi(vector<double> cv);

  void setState(StateCode sc);  // sets AA state (H, S, T, C)

  void setDefault();
  void setBonds(double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng); 
  
  void constructSideChain(SideChain& sc, vector<double> chi);
  void constructSideChain(SideChain& sc, double chi1, double chi2, 
			  double chi3, double chi4, double chi5);
  void setSideChain(SideChain& sc);
  void setStateFromTorsionAngles();  // sets states from its phi/psi angles
  void adjustLeadingN();             // adjusts the translation of the N atom
  void addTerminalOXT();             // adds an OXT atom to the C terminus
  void addMissingO();                // adds an O atom where missing
  void removeHAtomsfromLeadingNH3(); // removes H atoms in excess from NH3+
  void removeSideChain();

  void patchBetaPosition(unsigned int n = 0);
  // adds a CB atom to the sidechain, if necessary
  void patchAminoAcidCode(); // determines the 3-letter code from the sidechain
  bool setBondsFromPdbCode(bool connect, AminoAcid* prev = NULL, 
  bool permissive = false); // returns false if failed to connect
  virtual void sync(); // synchronize coords with structure
  virtual void setModified();

  void load(Loader& l);  // data loader

  const AminoAcid& getInBond(unsigned int n) const;
  AminoAcid& getInBond(unsigned int n);
  const AminoAcid& getOutBond(unsigned int n) const;
  AminoAcid& getOutBond(unsigned int n);

  virtual const Atom& getOpenInBondRef(unsigned int n = 0) const;
  virtual Atom& getOpenInBondRef(unsigned int n = 0);
  virtual const Atom& getOpenOutBondRef(unsigned int n = 0) const;
  virtual Atom& getOpenOutBondRef(unsigned int n = 0);

  virtual Component* clone();

  virtual void acceptCalculator(EnergyVisitor* v);
  virtual void acceptOptimizer(OptimizationVisitor* v);
   
  // OPERATORS:
  bool operator==(const AminoAcid& other) const;
  bool operator!=(const AminoAcid& other) const;
  AminoAcid& operator=(const AminoAcid& orig);
  Atom& operator[](unsigned int n);
  const Atom& operator[](unsigned int n) const;
  Atom& operator[](const AtomCode& ac);
  const Atom& operator[](const AtomCode& ac) const;

protected:

  // HELPERS:
  virtual void resetBoundaries();

private:

  // ATTRIBUTES:
  static double BOND_ANGLE_N_TO_CB;
  static double BOND_ANGLE_CB_TO_C;
  static double BOND_LENGTH_CA_TO_CB;

  AminoAcidCode type;	         // AminoAcid type
  StateCode state;               //  -- '' -- state (H, S, T, C)
  double phi, psi, omega;        // torsion angles
  SideChain sideChain;
  // vector atoms (inherited from Group) contains the backbone
  IntCoordConverter icc; // should be static, but compiler won't accept it?
};

// ---------------------------------------------------------------------------
//                                    AminoAcid
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:
/*
 * @description returns the three letter code of the amino acid
 */
inline unsigned int AminoAcid::getCode() const{
  return aminoAcidThreeLetterTranslator(id);
}
/*
 * @description returns the one letter code of the amino acid
 */
inline char AminoAcid::getType1L() { 
  return threeLetter2OneLetter(id); 
}
/*
 * @description returns the size of the amino acid, sidechain and group size
 */
inline unsigned int AminoAcid::size() const {
  return (Group::size() + sideChain.size());
}
/*
 * @description returns the backbone size of the amino acid 
 * @param
 * @return
 */
inline unsigned int AminoAcid::sizeBackbone() const{
  return Group::size();
}
/*
 * @description returns the value of phi angle 
 * @param override flag, recalculates the phi angle and updates it internally(bool)
 * @return angle value(double)
 */
inline double AminoAcid::getPhi(bool override){
  if ((phi > 990) || (override))
    if (sizeInBonds())
      phi = RAD2DEG * icc.getTorsionAngle(getInBond(0)[C], (*this)[N], 
				      (*this)[CA], (*this)[C]);
  return phi;
}
/*
 * @description returns the value of psi angle 
 * @param override flag, recalculates the psi angle and updates it internally(bool)
 * @return angle value(double)
 */
inline double AminoAcid::getPsi(bool override){
  if ((psi > 990) || (override))
    if (sizeOutBonds())
      psi = RAD2DEG * icc.getTorsionAngle((*this)[N], (*this)[CA], 
				      (*this)[C], getOutBond(0)[N]);
  return psi;
}
/*
 * @description returns the value of omega angle 
 * @param override flag, recalculates the omega angle and updates it internally(bool)
 * @return angle value(double)
 */
inline double AminoAcid::getOmega(bool override){
  if ((omega > 990) || (override))
    if (sizeOutBonds())
      omega = RAD2DEG * icc.getTorsionAngle((*this)[CA], (*this)[C], 
				      getOutBond(0)[N], getOutBond(0)[CA]);
  return omega;
}
/*
 * @description returns the value of chi angle 
 * @param returns the angle corresponding to the n position of the side chain
 * @return angle value(double)
 */
inline double AminoAcid::getChi(unsigned int n){
  return sideChain.getChi(n);
}
/*
 * @description returns the values of chi angle 
 * @param none
 * @return angle values(vector<double>)
 */
inline vector<double> AminoAcid::getChi(){
  return sideChain.getChi();
}
/*
 * @description returns the maximum value of chi angle 
 * @param none
 * @return angle value( unsigned int)
 */
inline unsigned int AminoAcid::getMaxChi(){
  return sideChain.getMaxChi();
}

/*
 * @description return the amino acid state
 * @param none
 * @return the code of the state(StateCode)
 */
inline StateCode AminoAcid::getState(){
  return state;
}
/*
 * @description returns the amino acid side chain
 * @param none
 * @return reference to the side chain( const SideChain& )
 */
inline const SideChain& AminoAcid::getSideChain() const{
  return sideChain;
}
/*
 * @description returns the amino acid side chain
 * @param none
 * @return reference to the side chain(  SideChain& )
 */
inline SideChain& AminoAcid::getSideChain(){
  return sideChain;
}
/*
 * @description verifies if the atom is present in the amino acid
 * @param reference to the atom code(AtomCode$)
 * @return flag to verify the presence(bool)
 */
inline bool AminoAcid::isMember(const AtomCode& ac) const{
  return (pGetAtom(ac) != NULL);
}
/*
 * @description Saves the amino acid
 * @param reference to the saver there the amino acid will be saved(saver&)
 * @return changes are made internally(void)
 */
inline void AminoAcid::save(Saver& s){
  s.saveAminoAcid(*this);
}


// MODIFIERS:
/*
 * @description defines the amino acid type
 * @param amino acid one letter code(char)
 * @return changes are made internally(void)
 */
inline void AminoAcid::setType1L(char _name) { 
  type = aminoAcidOneLetterTranslator(_name);
  id.setName(oneLetter2ThreeLetter(_name));
}
/*
 * @description defines the amino acid type
 * @param amino acid three letter code(string)
 * @return changes are made internally(void)
 */
inline void AminoAcid::setType(string _name){ 
  type = aminoAcidThreeLetterTranslator(_name);
  id.setName(_name);
}
/*
 * @description defines the chi value for a specific chi position
 * @param chi index(unsigned int), chi value(double)
 * @return  changes are made internally(void)
 */
inline void AminoAcid::setChi(unsigned int n, double a){
  sideChain.setChi(n, a);
}
/*
 * @description defines the chi values of the side chain
 * @param  chi values(vector<double>)
 * @return  changes are made internally(void)
 */
inline void AminoAcid::setChi(vector<double> cv){
  sideChain.setChi(cv);
}
/*
 * @description sets the state of the amino acid
 * @param state (StateCode)
 * @return changes are made internally(void)
 */
inline void AminoAcid::setState(StateCode sc){
  state = sc;
}
/*
 * @description removes the side chain 
 * @param none
 * @return changes are made internally(void)
 */
inline void AminoAcid::removeSideChain(){
  SideChain sc;
  sideChain = sc;
  resetBoundaries();
}
/*
 * @description Loads the amino amino acid
 * @param reference to the loader(loader&);
 * @return changes are made internally(void)
 */
inline void AminoAcid::load(Loader& l){
  l.loadAminoAcid(*this);
  resetBoundaries();
}
/*
 * @description returns the In bond n
 * @param value for n (unsigned int)
 * @return reference to the amino acid(const AminoAcid&)
 */
inline const AminoAcid& AminoAcid::getInBond(unsigned int n) const{
  return dynamic_cast<const AminoAcid&>(Bond::getInBond(n));
}
/*
 * @description returns the In bond n
 * @param value for n (unsigned int)
 * @return reference to the amino acid( AminoAcid&)
 */
inline AminoAcid& AminoAcid::getInBond(unsigned int n){
  return dynamic_cast<AminoAcid&>(Bond::getInBond(n));
}
/*
 * @description returns the out bond n
 * @param value for n (unsigned int)
 * @return reference to the amino acid(const AminoAcid&)
 */
inline const AminoAcid& AminoAcid::getOutBond(unsigned int n) const{
  return dynamic_cast<const AminoAcid&>(Bond::getOutBond(n));
}
/*
 * @description returns the out bond n
 * @param value for n (unsigned int)
 * @return reference to the amino acid( AminoAcid&)
 */
inline AminoAcid& AminoAcid::getOutBond(unsigned int n){
  return dynamic_cast<AminoAcid&>(Bond::getOutBond(n));
}
/*
 * @description fixes the amino acid, determining the the 3-letter code from the sidechain
 * @param none
 * @return changes are made internally(void)
 */
inline void AminoAcid::patchAminoAcidCode() {  
  sideChain.patchAminoAcidCode();
  setType(sideChain.getType());
}
/*
 * @description defines a flag for the group and for the side chain
 * @param none
 * @return changes are made internally(void)
 */
inline void AminoAcid::setModified(){
  Group::setModified();
  sideChain.setModified();
}

/*
 * @description defines the energy visitor as an accepted calculator of the amino acid
 * @param pointer to the energy visitor
 * @return changes are made internally(void)
 */
inline void AminoAcid::acceptCalculator(EnergyVisitor* v){
  v->PrepareAminoAcid(*this);
}
/*
 * @description defines the energy visitor as an optimization visitor of the amino acid
 * @param pointer to the optimization visitor
 * @return changes are made internally(void)
 */
inline void AminoAcid::acceptOptimizer(OptimizationVisitor* v){
  v->PrepareAminoAcid(*this);
}


// OPERATORS:
/*
 * @description allows to compare an amino acid to another
 * @param reference to the amino acid(const AminoAcid&)
 * @return flag that verifies that the amino acids are equal(bool)
 */
inline bool AminoAcid::operator==(const AminoAcid& other) const{
  return (dynamic_cast<const Identity*>(this)) == 
    dynamic_cast<const Identity*>(&other);
}
/*
 * @description allows to compare an amino acid to another
 * @param reference to the amino acid(const AminoAcid&)
 * @return flag that verifies that the amino acids are not equal(bool)
 */
inline bool AminoAcid::operator!=(const AminoAcid& other) const{
  return (dynamic_cast<const Identity*>(this)) != 
    dynamic_cast<const Identity*>(&other);
}

/*
 * @description allows to see the amino acid as a vector of atoms, and obtain the n position of it
 * @param index of the atom in the amino acid(unsigned int)
 * @return reference to the atom(Atom&)
 */
inline Atom& AminoAcid::operator[](unsigned int n){
  PRECOND(n < this->size(), exception);
  return ((n < Group::size()) ? Group::getAtom(n) : sideChain[n-Group::size()]);
}
/*
 * @description allows to see the amino acid as a vector of atoms, and obtain the n position of it
 * @param index of the atom in the amino acid(unsigned int)
 * @return reference to the atom(cont Atom&)
 */
inline const Atom& AminoAcid::operator[](unsigned int n) const{
  PRECOND(n < this->size(), exception);
  return ((n < Group::size()) ? Group::getAtom(n) : sideChain[n-Group::size()]);
}
/*
 * @description allows to see the amino acid as a vector of atoms, and obtain the n position of it
 * @param code of the atom in the amino acid(const AtomCode& )
 * @return reference to the atom(Atom&)
 */
inline Atom& AminoAcid::operator[](const AtomCode& ac){
  Atom* a = pGetAtom(ac);
  INVARIANT(a != NULL, exception);
  return *a;
}
/*
 * @description allows to see the amino acid as a vector of atoms, and obtain the n position of it
 * @param code of the atom in the amino acid(const AtomCode& )
 * @return reference to the atom(const Atom&)
 */
inline const Atom& AminoAcid::operator[](const AtomCode& ac) const{
  Atom* a = pGetAtom(ac);
  INVARIANT(a != NULL, exception);
  return *a;
}

// HELPERS:
/** @example AminoAcidTest.cc
   *  A simple program to test class AminoAcid's features.
 */
} // namespace
#endif //_AMINOACID_H_

