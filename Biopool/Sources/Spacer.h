/**
* @Class:              Spacer
* @Base Class(es):     Identity, Bond
* @Containing:         Group, AminoAcid
* @Author:             Silvio Tosatto
* @Project Name:       Victor
*/
#ifndef _SPACER_H_
#define _SPACER_H_


// Includes:
#include <Bond.h>
#include <Polymer.h>
#include <AminoAcid.h>
#include <Visitor.h>
#include <AminoAcidCode.h>
#include <map>
#include <set>


// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief This class implements a "Spacer" for a protein chain. Include methods to obtain values from the atoms and its pdb information.
 * 
* @Description The current implementation allows for "1 to 1" spacers,
*    ie. spacers composed of a single aminoacid chain. Class that develop the virtual methods from the parent class
 *   Polymer.
 * */

class Spacer : public Polymer{
public: 
  //AK11
  friend class AtomEnergyCalculator;
  friend class SideChainPlacement;
  
  // CONSTRUCTORS/DESTRUCTOR:
  Spacer();
  Spacer(const Spacer& orig); 
  virtual ~Spacer(); 
  
  // PREDICATES:

  int getStartOffset() { return startOffset; }
  int getAtomStartOffset() { return startAtomOffset; }
  int getIndexFromPdbNumber(int index);  // in: pdb_num, out: array_num
  int getPdbNumberFromIndex(int index);  // in: array_num, out: pdb_num
  int maxPdbNumber() { return startOffset + gaps.size() + sizeAmino(); }
  bool isGap(int index);
  void printGaps();
  unsigned int sizeGaps()  { return gaps.size(); }

  virtual string getClassName() const    { return "Spacer"; }
  void save(Saver& s);   // data saver
  AminoAcid& getAmino(unsigned int n);
  const AminoAcid& getAmino(unsigned int n) const;
  const unsigned int  sizeAmino() const;
  
  Spacer& getSpacer(unsigned int n);  
  const Spacer& getSpacer(unsigned int n) const; 
  const unsigned int sizeSpacer() const;

  vector<pair<unsigned int, unsigned int> > getHelixData();
  vector<pair<unsigned int, unsigned int> > getStrandData();

  // MODIFIERS:
  
  ///inserts a new aminoacid after given position setting given torsion angles and atomic distances
  void insertAminoAfter(string code, unsigned int n = 9999, double phi = -62, 
			double psi = -41, double omega = 180, 
			double NToCaLen = 1.458, double CaToCLen = 1.52, double CToOLen = 1.231, 
			double atCaToCAng = 111.6, double CaToOAng = 120.80,
			double CtoNLen = 1.33, double atCToNAng = 116.4, 
			double OToNAng = 123.2, double atNToCaAng = 121.9);
  void insertAminoAfterWithGaps(string code, unsigned int n = 9999, double phi = -62, 
			double psi = -41, double omega = 180,
			int beginHole=0, int endHole=0,
			vgVector3 <double> ca1=vgVector3<double>(), vgVector3 <double> ca2=vgVector3<double>(),
			string target="", Spacer* refSpacer=0, Spacer* pOriginalSpacer=0,
			double NToCaLen = 1.458, double CaToCLen = 1.52, double CToOLen = 1.231, 
			double atCaToCAng = 111.6, double CaToOAng = 120.80,
			double CtoNLen = 1.33, double atCToNAng = 116.4, 
			double OToNAng = 123.2, double atNToCaAng = 121.9 );
  ///inserts a new aminoacid before given position setting given torsion angles and atomic distances
  void insertAminoBefore(string code, unsigned int p = 0, double phi = -62, 
			 double psi = -41, double omega = 180,
			 double NToCaLen = 1.458, double CaToCLen = 1.52, double CToOLen = 1.231, 
			 double atCaToCAng = 111.6, double CaToOAng = 120.80,
			 double CtoNLen = 1.33, double atCToNAng = 116.4, 
			 double OToNAng = 123.2, double atNToCaAng = 121.9);
  void setStartOffset(int _offset);
  void setAtomStartOffset(int _offset) { startAtomOffset = _offset; }

  void addGap(int index);
  void removeGap(int index);
  void removeAllGaps();
        
  void insertSubSpacerAfter(Spacer* sp, unsigned int pos);
 
  void insertFirstSpacer(Spacer* sp);
 

  void insertComponent(Component* c);
  void removeComponent(Component* c);
  //-------------
  void removeComponentFromIndex(unsigned int i);
  //--------------
  void deleteComponent(Component* c);

  void mergeSpacer( Spacer* s);
  Spacer* splitSpacer( unsigned int index, unsigned int count);

  void copy(const Spacer& orig);

  void load(Loader& l);  // data loader

  void setTrans(vgVector3<double> t);
  void addTrans(vgVector3<double> t);
  void setRot(vgMatrix3<double> r);
  void addRot(vgMatrix3<double> r);
  void sync(); 

  void setStateFromSecondary(string sec);  
  void setStateFromTorsionAngles();  
  void setDSSP(bool verbose);
  vector< set<char> > getDSSP(){return ss;};
  
 
  const Spacer& getInBond(unsigned int n) const;
  Spacer& getInBond(unsigned int n);
  const Spacer& getOutBond(unsigned int n) const;
  Spacer& getOutBond(unsigned int n);

  void  makeFlat();
  void  groupLikeStateCode();

  virtual const Atom& getOpenInBondRef(unsigned int n = 0) const;
  virtual Atom& getOpenInBondRef(unsigned int n = 0);
  virtual const Atom& getOpenOutBondRef(unsigned int n = 0) const;
  virtual Atom& getOpenOutBondRef(unsigned int n = 0);

  virtual Component* clone();

  virtual void acceptCalculator(EnergyVisitor* v);
  virtual void acceptOptimizer(OptimizationVisitor* v);
  
  // OPERATORS:
  Spacer& operator=(const Spacer& orig);

  // TESTERS AND DEVELOPERS
  void printSubSpacerList();
  void printComponents();
  vector< pair<int,int> > getHoles();

  static double BOND_ANGLE_AT_CPRIME_TO_N;
  static double BOND_LENGTH_CPRIME_TO_N;
  static double BOND_ANGLE_CA_TO_O;
  static double BOND_LENGTH_C_TO_O;
  static double BOND_LENGTH_CALPHA_TO_CPRIME;
  static double BOND_LENGTH_N_TO_CALPHA;
  static double BOND_ANGLE_AT_CALPHA_TO_CPRIME;
  static double BOND_ANGLE_O_TO_N;
  static double BOND_ANGLE_AT_N_TO_CALPHA;
  
  // 
  static bool NMRGetMinimumCADistanceVector( const string &strFileName, vector<double> *pNMR );
    
protected:

  // HELPERS:
  void resetBoundaries();
  void modifySubSpacerList(Spacer* ,int);
  void updateSubSpacerList();
  void getBackboneHbonds();  // backbone H bonds (for SS)

  // ATTRIBUTES

  int startOffset;       //number (as reported in the pdb file) of first amminoacid loaded 
  int startAtomOffset;   //number of first atom of the first amminoacid
  vector<int> gaps;

  pair<unsigned int, unsigned int> getSubSpacerListEntry(unsigned int); 
  pair<unsigned int, unsigned int> getSubSpacerListEntry(unsigned int) const;
  
  bool **backboneHbonds;
  vector<set<char > > ss;

private:
  vector<pair<unsigned int, unsigned int> > subSpacerList;  
};

// ---------------------------------------------------------------------------
//                                    Spacer
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:

/**
 * @Description Saves the information from the saver into the Spacer
 * @param reference to the saver(Saver&)
 * @return changes are made internally(void)
 */
inline void Spacer::save(Saver& s){
  s.saveSpacer(*this);
}
/**
 * @Description Allows to know if there is a spacer in the protein or not. 
 * @param none
 * @return quantity of spacers, possible values are 0 or 1(const unsigned int).
 */
inline const unsigned int Spacer::sizeSpacer() const{
  return subSpacerList.size();
}

// MODIFIERS:
/**
 * @Description Removes the corresponding component
 * @param pointer to the component to remove(Component*), usually the object that is calling the method.
 * @return  changes are made internally(void)
 */
inline void Spacer::removeComponent(Component* c){
  Polymer::removeComponent(c);
  setModified();
}

/**
 * @Description Removes the corresponding component in the corresponding index
 * @param index to the component to remove(unsigned int), usually the index of the object that is calling the method.
 * @return  changes are made internally(void)
 */
inline void Spacer::removeComponentFromIndex(unsigned int i){
  Polymer::removeComponentFromIndex(i);
  setModified();
}

/**
 * @Description Similar to removeComponent but also free the component's memory space
 * @param pointer to the Component to delete(Component*)
 * @return changes are made internally(void)
 */inline void Spacer::deleteComponent(Component* c){
  Polymer::deleteComponent(c);
  setModified();
}

/**
 * @Description Loads the Spacer from a Loader
 * @param reference of the Loader(Loader&)
 * @return changes are made internally(void)
 */
inline void Spacer::load(Loader& l){
  l.loadSpacer(*this);
  resetBoundaries();
}

/**
 * @Description Sets the translation vector
 * @param the translation vector to set(vgVector3<double>)
 * @return changes are made internally(void)
 */
inline void Spacer::setTrans(vgVector3<double> t){
  if (sizeAmino())
    getAmino(0).setTrans(t);
}

/**
 * @Description Adds a the translation vector
 * @param the translation vector to add(vgVector3<double>)
 * @return changes are made internally(void)
 */
inline void Spacer::addTrans(vgVector3<double> t){
  if (sizeAmino())
    getAmino(0).addTrans(t);
}

/**
 * @Description Sets the rotation matrix
 * @param the rotation matrix to set(vgMatrix3<double>)
 * @return changes are made internally(void)
 */
inline void Spacer::setRot(vgMatrix3<double> r){
  if (sizeAmino())
    getAmino(0).setRot(r);
}

/**
 * @Description Adds a rotation matrix
 * @param the rotation matrix to add(vgMatrix3<double>)
 * @return changes are made internally(void)
 */
inline void Spacer::addRot(vgMatrix3<double> r){
  if (sizeAmino())
    getAmino(0).addRot(r);
}


/**
 * @Description Returns the reference to the spacer that contains the inBond for the n Atom
 * @param index of the atom(unsigned int)
 * @return reference to the spacer that contains the inBond information (const Spacer& )
 */
inline const Spacer& Spacer::getInBond(unsigned int n) const{
  return dynamic_cast<const Spacer&>(Bond::getInBond(n));
}

/**
 * @Description Returns the reference to the spacer that contains the inBond for the n Atom
 * @param index of the atom(unsigned int)
 * @return reference to the spacer that contains the inBond information (const Spacer& )
 */
inline Spacer& Spacer::getInBond(unsigned int n){
  return dynamic_cast<Spacer&>(Bond::getInBond(n));
}

/**
 * @Description Returns the reference to the spacer that contains the outBond for the n Atom
 * @param index of the atom(unsigned int)
 * @return reference to the spacer that contains the outBond information (const Spacer& )
 */
inline const Spacer& Spacer::getOutBond(unsigned int n) const{
  return dynamic_cast<const Spacer&>(Bond::getOutBond(n));
}

/**
 * @Description Returns the reference to the spacer that contains the outBond for the n Atom
 * @param index of the atom(unsigned int)
 * @return reference to the spacer that contains the outBond information (const Spacer& )
 */
inline Spacer& Spacer::getOutBond(unsigned int n){
  return dynamic_cast<Spacer&>(Bond::getOutBond(n));
}

/**
 * @Description Sets the Spacer as part of the Energy Visitor object
 * @param pointer to the Energy visitor object(EnergyVisitor* )
 * @return changes are made internally(void)
 */
inline void Spacer::acceptCalculator(EnergyVisitor* v){
  v->PrepareSpacer(*this);
}

/**
 * @Description Sets the Spacer as part of the Optimization Visitor object
 * @param pointer to the Optimization Visitor(OptimizationVisitor*)
 * @return changes are made internally(void)
 */
inline void Spacer::acceptOptimizer(OptimizationVisitor* v){
  v->PrepareSpacer(*this);
}





// OPERATORS:

} // namespace
#endif //_SPACER_H_
