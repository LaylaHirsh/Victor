/** @mainpage Basic information
 *
 * @section intro_sec Introduction
 * The victor library .......
 * 
 * @section BiopoolFeatures Introduction of Biopool
 * The Biopool library (Biopolimer Object Oriented Library) was created in 2000 and developed over the years.
 * Its the central part for the Homer project.In it, all the basic classes for representing the protein's 
 * structure are developed, and include all the needed methods to manage it. The final purpose is to representate 
 * one aminoacidic string in a efficient way. Including the ableness of reading the lineal sequence of
 * aminoacids or to process one PDB structure. 
 *  
 * @section struct_sec Class structure
 * In this library we can identify 3 class groups:
 * 
 * @subsection step1 Main Group: 
 *  @image html BiopoolSimpleBond.jpg
  * @subsection step2 Saver Group: 
 *  @image html BiopoolSaver.jpg
  * @subsection step3 Loader Group: 
 *  @image html BiopoolLoader.jpg
 * 
 * @section Classstruct_sec Basic Concept
 *
 * @subsection step11 Scheme for Protein
 *    @image html SchemeProteinclass.jpg
 * Explanation: The protein has a vector of components, which are actually polymers. The first has the index 0 the second index 1, and so on.
 * Each polymer has also a vector of components but this time the first one(index(0)) is one spacer and the second one(index 1) is the ligand set.
 * Finally each Spacer is formed by AminoAcids, and a Ligand set is formed by Ligands. AminoAcids and Ligands are formed by Atoms.
 * @section Example Basic Sample
 * 
 * @subsection step12Protein Test:
 *       From now on the Biopool folder will be the folder were all the source codes will be.
 *       To run the example just write the next line on the terminal.
 *@li    ./ProteinTest -i 4GKV.pdb -o output.pdb (will create a file containing only the ATOMs from the pdb original file, assuming the first chain)
 *@li    ./ProteinTest -i 4GKV.pdb -o output.pdb --all (will create a file containing only the ATOMs from the pdb original file, for all chains, limited by TER tag)
 *@li    ./ProteinTest -i 4GKV.pdb -o output.pdb -c B (will create a file containing only the ATOMs from the pdb original file, for B chain)
 *@li    ./ProteinTest -i 2MKP.pdb -o output.pdb -m 1 (will create a file containing only the ATOMs from the first model)
 * You can find the executables in the project bin folder
 
 * @subsection step13 Options: 
 **@li-i <filename> 	 Input PDB file
 *@li -o <filename> 	 Output to file
 *@li -c <id>       	 Chain identifier to read(default is first chain)
 *@li --all         	 or type '-all' to select all chains analisys)
*@li -m <number>   	 Model number to read (NMR only, default is first model)
*@li --hetatm      	 type '-hetatm' if you want hetAtm
*@li --metalOnly   	 setOnlyMetalHetAtoms; require --hetatm option.
*@li --water       	 to enable water selection; require --hetatm option(and metalOnly not).
*Default options: verbose=true; noSecondary=true; noWater= true;
 *  etc...
 */
 
 /**  @namespace Biopool
 *   @brief This library contains methods that will allow you to representate one aminoacidic string in a efficient way. 
 *    
  *  @Description Including the ableness of reading the lineal sequence of aminoacids or to process one PDB structure.
 */


/**
* @Class:             AminoAcid
* @Author:            Silvio Tosatto
* @Project Name:      Victor
*/

// Includes:
#include <AminoAcid.h>
#include <Debug.h>
#include <limits.h>
#include <float.h>
#include <AtomCode.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Biopool;
double AminoAcid::BOND_ANGLE_N_TO_CB = 110.5;
double AminoAcid::BOND_ANGLE_CB_TO_C = 110.1;
double AminoAcid::BOND_LENGTH_CA_TO_CB = 1.53;
/**
 * @Description Constructor
* @param none
*/
AminoAcid::AminoAcid() : Group(1,1), type(XXX), state(COIL), phi(999), 
  psi(999), omega(999), sideChain(), icc()
{ sideChain.setBackboneRef(this); }


/**
 * @Description Constructor
* @param AminoAcid
*/
AminoAcid::AminoAcid(const AminoAcid& orig){
  PRINT_NAME;
  this->copy(orig);
}
/**
 * @Description Destructor
 * @param none
 */
AminoAcid::~AminoAcid()
{ PRINT_NAME; } 
/**
 * @Description Returns the atom corresponding to the open in bond  
 * aminoacids have only a single possible, hard coded, open in-bond
 * @param index (unsigned int) 
 * @return reference to the open in bond(const Atom&) 
 */
const Atom& AminoAcid::getOpenInBondRef(unsigned int n) const{
  PRECOND( (n == 0) && (sizeOpenInBonds() > 0), exception);
  return (*this)[N];
}
/**
 * @Description Returns the atom corresponding to the open in bond  
 * aminoacids have only a single possible, hard coded, open in-bond
 * @param index (unsigned int) 
 * @return reference to the open in bond( Atom&) 
 */
Atom& AminoAcid::getOpenInBondRef(unsigned int n){
  PRECOND( (n == 0) && (sizeOpenInBonds() > 0), exception);
  return (*this)[N];
}
/**
 * @Description Returns the atom corresponding to the open out bond  
 * aminoacids have only a single possible, hard coded, open out-bond
 * @param index (unsigned int) 
 * @return reference to the open in bond( Atom&) 
 */
const Atom& AminoAcid::getOpenOutBondRef(unsigned int n) const{
  // NB: aminoacids have only a single possible, hard coded, open out-bond
  PRECOND( (n == 0) && (sizeOpenOutBonds() > 0), exception);
  return (*this)[C];
}
/**
 * @Description Returns the atom corresponding to the open out bond  
 * aminoacids have only a single possible, hard coded, open out-bond
 * @param index (unsigned int) 
 * @return reference to the open in bond( Atom&) 
 */
Atom& AminoAcid::getOpenOutBondRef(unsigned int n){
  // NB: aminoacids have only a single possible, hard coded, open out-bond
  PRECOND( (n == 0) && (sizeOpenOutBonds() > 0), exception);
  return (*this)[C];
}

// MODIFIERS:
/**
 * @Description Connects the aa to the reverse structure
 * @param pointer to the aminoacid (AminoAcid*) , offset (unsigned int - not used )
 * @return  changes are made internally(void) 
 */
void AminoAcid::connectIn(AminoAcid* a, unsigned int offset){ 
  if (!(a->isMember(OXT)))
    a->addTerminalOXT();
  IntCoordConverter icc;
  a->setOmega(180.0);
  icc.connectReverseStructure((*this), (*a));
  sync();
  resetBoundaries();
}
// MODIFIERS:
/**
 * @Description Connects the structure  to the aa 
 * @param pointer to the aminoacid (AminoAcid*) , offset (unsigned int - not used )
 * @return  changes are made internally(void) 
 * */
void AminoAcid::connectOut(AminoAcid* a, unsigned int offset){
  if (!isMember(OXT))
    addTerminalOXT();
  IntCoordConverter icc;
//    if (offset % 2 == 1)
    setOmega(180.0);
  icc.connectStructure((*a), (*this));
  resetBoundaries();
}
/**
 * @Description unbind to structures from the predecessor
 * @param none
 * @return  pointer to the amino acid unbind (AminoAcid*)
 */
AminoAcid* AminoAcid::unconnectIn(){
  if (!sizeInBonds())    {
      DEBUG_MSG("Cannot unconnect aminoacid without predecessor.");
      return NULL;
    }
  AminoAcid* tmp = &getInBond(0);

  // set OXT for tmp:
  Atom at;
  at.setCode(OXT);
  at.setTrans((*this)[N].getTrans());
  at.bindIn((*tmp)[C]);
  tmp->addAtom(at);

  // unbind the two structures and reset the positions:
  unbindIn(*tmp);
  (*this)[N].unbindIn((*tmp)[C]);
  setTrans((*tmp)[C].getCoords());
  tmp->setOmega(999);

  resetBoundaries();
  tmp->resetBoundaries();
  return tmp;
}
/**
 * @Description unbind to structures from the follower
 * @param none
 * @return  pointer to the amino acid unbind (AminoAcid*)
 */
AminoAcid* AminoAcid::unconnectOut(){
  if (!sizeOutBonds())    {
      DEBUG_MSG("Cannot unconnect aminoacid without follower.");
      return NULL;
    }
  AminoAcid* tmp = &getOutBond(0);

  // set OXT for this:
  Atom at;
  at.setCode(OXT);
  at.setTrans((*tmp)[N].getTrans());
  at.bindIn((*this)[C]);
  addAtom(at);
  // unbind the two structures and reset the positions:
  unbindOut(*tmp);
  (*this)[C].unbindOut((*tmp)[N]);
  tmp->setTrans((*this)[C].getCoords());
  setOmega(999);

  resetBoundaries();
  tmp->resetBoundaries();
  return tmp;
}
/**
* @Description Sets phi angle
 * @param angle value(double)
 * @return  changes are made internally(void)
 */
void AminoAcid::setPhi(double a){
  if (a >= 990)    {
      phi = a;
      return;
    }

  if (a < -180)
    a += 360;
  else if (a > 180)
    a -= 360;

  PRECOND((a >= -180) && (a <= 180), exception);
  if (sizeInBonds())
    icc.setTorsionAngle( getInBond(0)[C], (*this)[N], (*this)[CA], 
			 (*this)[C], DEG2RAD * a);
  phi = a;
  sync();
}
/**
* @Description Sets psi angle
 * @param angle value(double)
 * @return  changes are made internally(void)
 */
void AminoAcid::setPsi(double a){
  if (a >= 990) {
      psi = a;
      return;
    }

  if (a < -180)
    a += 360;
  else if (a > 180)
    a -= 360;

  PRECOND((a >= -180) && (a <= 180), exception);
  if (sizeOutBonds())
    icc.setTorsionAngle( (*this)[N], (*this)[CA], (*this)[C],
			 getOutBond(0)[N], DEG2RAD * a);
  psi = a;
  sync();
}
/**
* @Description Sets omega angle
 * @param angle value(double)
 * @return  changes are made internally(void)
 */
void AminoAcid::setOmega(double a){
  if (a >= 990){
      omega = a;
      return;
    }

  if (a < -180)
    a += 360;
  else if (a > 180)
    a -= 360;

  PRECOND((a >= -180) && (a <= 180), exception);
  if (sizeOutBonds())

    icc.setTorsionAngle( (*this)[CA], (*this)[C], getOutBond(0)[N], 
			 getOutBond(0)[CA], DEG2RAD * a);
  omega = a;
  sync();
}

/**
* @Description Sets side chain
 * @param reference to the side chain (SideChain&) 
 * @return  changes are made internally(void)
 */
void AminoAcid::setSideChain(SideChain& sc){
  if (sc.getType() == "GLY")
    return;

  if (sideChain.size() == 0)
    patchBetaPosition();

  if (sideChain.size() == 0)
    ERROR("Sidechain undefined.", exception);

  vgVector3<double> tmpV = sideChain[0].getTrans();
  vgMatrix3<double> res(1);

  sideChain = sc;
  sideChain.setBackboneRef(this);

  alignVectors( tmpV, sideChain[0].getTrans(), res);

  sideChain[0].addRot(res);
  sync();
}
/**
* @Description Construct side chain
 * @param reference to the Sidechain(SideChain&), values for the chi angles as much as needed, max 5 (double, double , double,
 *  , double ,double)
 * @return   changes are made internally(void)
 */
void AminoAcid::constructSideChain(SideChain& sc, double chi1, double chi2, 
			      double chi3, double chi4, double chi5){
  setSideChain(sc);
  unsigned int maxChi = getMaxChi();
  switch(maxChi)    {
    case 1:
      {
	//control temporarily deactivated to allow computation in some erroneus cases 
	if (chi2 != 999) {
	    ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
	  }
	setChi(0, chi1);
	break;
      }
    case 2:
      {
	if (chi3 != 999) {
	    ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
	  }
	setChi(0, chi1);
	setChi(1, chi2);
	break;
      }
    case 3:
      {
	if (chi4 != 999) {
	    ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
	  }
	setChi(0, chi1);
	setChi(1, chi2);
	setChi(2, chi3);
	break;
      }
    case 4:
      {
	if (chi5 != 999)  {
	    ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
	  }	
	setChi(0, chi1);
	setChi(1, chi2);
	setChi(2, chi3);
	setChi(3, chi4);
	break;
      }
    case 5:
      {
	setChi(0, chi1);
	setChi(1, chi2);
	setChi(2, chi3);
	setChi(3, chi4);
	setChi(4, chi5);
	break;
      }
    }
  sync();
}
/**
* @Description Construct side chain
 * @param reference to the Sidechain(SideChain&), vector containing values for the chi angles ( vector<double>)
 * @return   changes are made internally(void)
 */
void AminoAcid::constructSideChain(SideChain& sc, vector<double> chi){
  setSideChain(sc);
  
  if (chi.size() > getSideChain().getMaxChi())
    ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception); 
  
  setChi(chi);
  sync();
}

/**
* @Description Sets the state from torsion angles from definition was taken from McGuffin et al.,
  // Bioinformatics (17):63-72 (2001)
 * @param none
 * @return   changes are made internally(void)
 */
void AminoAcid::setStateFromTorsionAngles(){

  double phi = getPhi();
  double psi = getPsi();

  if ( ((phi >= -80) && (phi <= -40) && (psi >= -65) && (psi <= -5))
       || ((phi >= -110) && (phi <= -40) && (psi >= -74) && (psi <= -0)) )
    setState(HELIX);
  else if ( ((phi >= -180) && (phi <= -60) && (psi >= 60) && (psi <= 180))
       || ((phi >= -180) && (phi <= -60) && (psi >= -180) && (psi <= -140)) )
    setState(STRAND);

}
/**
* @Description adjusts the translation of the N atom, set trans for N relative to CA & C,
 *         adjust this' translation, the phi angle in leading structures is always undefined
 * @param none
 * @return   changes are made internally(void)
 */
void AminoAcid::adjustLeadingN(){ 
    if ((*this)[N].getTrans().length() != 0)
    return;
  const double BOND_ANGLE_N_TO_CA = DEG2RAD * 116.5;
  const double BOND_LENGTH_N_TO_CA = 1.45;
  vgVector3<double> normal = (*this)[CA].getTrans().cross(
		       (*this)[C].getTrans());
  vgMatrix3<double> res = vgMatrix3<double>::createRotationMatrix( 
			      normal, (DEG2RAD * 180.0) - BOND_ANGLE_N_TO_CA);
  (*this)[N].setTrans( BOND_LENGTH_N_TO_CA 
	       * (res * ((*this)[CA].getTrans())).normalize() );
  addTrans( -(*this)[N].getTrans() );
}
/**
* @Description adds an OXT atom to the C terminus
 * @param none
 * @return   changes are made internally(void)
 */
void AminoAcid::addTerminalOXT(){ // adds an OXT atom to the C terminus
  if (isMember(OXT))
    return;

  if ((*this)[C].sizeOutBonds() > 1)  {
      DEBUG_MSG("Cannot add terminal OXT on bound aminoacids.");
      return;
    }

  const double BOND_ANGLE_CA_TO_OXT = DEG2RAD * 117.0;
  const double BOND_ANGLE_O_TO_N = 122.0;
  const double BOND_LENGTH_C_TO_OXT = 1.35;
  Atom at;
  at.setType("OXT");
  at.setBFac(15.0);
  at.bindIn((*this)[C]);
  addAtom(at);
  IntCoordConverter icc;
  icc.zAtomToCartesian((*this)[C], BOND_LENGTH_C_TO_OXT, (*this)[CA], 
		       BOND_ANGLE_CA_TO_OXT, (*this)[O], 
		       BOND_ANGLE_O_TO_N, 1, (*this)[OXT]);
  sync();
}
/**
* @Description  adds an O atom, if missing
 * @param none
 * @return   changes are made internally(void)
 */
void AminoAcid::addMissingO(){
  if (isMember(O))
    return;

  if (sizeOutBonds() < 1)
    return;

  const double BOND_ANGLE_CA_TO_O = 120.80;
  const double BOND_ANGLE_N_TO_O = 123.0;
  const double BOND_LENGTH_C_TO_O = 1.231;

  Atom at;
  at.setType("O");
  at.bindIn((*this)[C]);
  addAtom(at);
  
  IntCoordConverter icc;
  icc.zAtomToCartesian((*this)[C], BOND_LENGTH_C_TO_O, (*this)[CA], 
		       BOND_ANGLE_CA_TO_O, getOutBond(0)[N], 
		       BOND_ANGLE_N_TO_O, 1, (*this)[O]);
  sync();
}
/**
* @Description  remove H Atoms from Leading NH3,
 * check for NH3+ and, if so, remove superfluous H (i.e. 2H, 3H) atoms
 * @param none
 * @return   changes are made internally(void)
 */
void AminoAcid::removeHAtomsfromLeadingNH3(){
  if (!isMember(H))
    return;
  
  bool first = false;
  for (unsigned int i = 0; i < (*this)[N].sizeOutBonds(); i++)    {
      if ((*this)[N].getOutBond(i).getCode() == H){
	if (!first)
	  first = true;
	else	  {	// remove H atom
	    Atom tmp = (*this)[N].getOutBond(i);
	    (*this)[N].unbindOut(tmp);
	    (*this).removeAtom(tmp);
	  }
      }
    }
}


/**
* @Description  adds a CB atom to the sidechain, if necessary
 * @param value of number for CB to set(unsigned int )
 * @return   changes are made internally(void)
 */
void AminoAcid::patchBetaPosition(unsigned int n){   
    if (type == GLY)
    return;
    
  Atom at;
  at.setType("CB");

  if (n > 0)
    at.setNumber(n);

  getSideChain().addAtom(at);
  
  (*this)[CA].bindOut(getSideChain()[CB]);

  IntCoordConverter icc;
  icc.zAtomToCartesian((*this)[CA], BOND_LENGTH_CA_TO_CB, 
		       (*this)[N], BOND_ANGLE_N_TO_CB, 
		       (*this)[C], BOND_ANGLE_CB_TO_C, 
		       1, (*this).getSideChain()[CB]);
  sync();
}

/**
* @Description  Copies an amino acid
 * @param reference to the  original aa to copy(const AminoAcid&)
 * @return   changes are made internally(void)
 */
void AminoAcid::copy(const AminoAcid& orig){
  PRINT_NAME; 
  Group::copy(orig);

  type = orig.type;
  phi = orig.phi;
  psi = orig.psi;
  omega = orig.omega;

  sideChain = orig.sideChain;

  if ( (orig.getSideChain().size() > 0) 
      &&  (orig[CA].isBond(orig.getSideChain()[0])) )
      sideChain.setBackboneRef(this);

  // fix Proline CD to N bond:
  if ((getType() == "PRO"))  {
      (*this)[N].setMaxInBonds(2);
      if (getSideChain().isMember(CD) && orig.getSideChain().isMember(CD) && 
	  ( orig[N].isBond(orig[CD])  ) )
  	getSideChain()[CD].bindOut((*this)[N]);
    }

  // set absolute position to orig's:
  if (orig[0].sizeInBonds()) {
      setTrans(const_cast<AminoAcid&>(orig)[0].getInBond(0).getCoords());
    }
}
/**
* @Description  Syncronize and sets boundaries
 * @param none
 * @return   changes are made internally(void)
 */
void AminoAcid::sync(){ 
  Group::sync();
  sideChain.sync();
  resetBoundaries();
}
/**
* @Description  Clone the amino acid
 * @param none
 * @return  pointer to the new(Component* )
 */
Component* AminoAcid::clone(){
  AminoAcid* tmp = new AminoAcid;
  tmp->copy(*this);
  return tmp;
}

// OPERATORS:
/**
* @Description  Operator =, assign the amino acid
 * @param AminoAcid reference to the original amino acid(const AminoAcid&)
 * @return  reference to the amino acid(AminoAcid&)
 */
AminoAcid& AminoAcid::operator=(const AminoAcid& orig){
  PRINT_NAME;
  if (&orig != this)
    copy(orig);
  return *this;
}

// HELPERS:
/**
* @Description  Reset the boundaries 
 * @param none
 * @return   changes are made internally(void)
 */
void AminoAcid::resetBoundaries(){
  Group::resetBoundaries();
  sideChain.resetBoundaries();

  for (unsigned int i = 0; i < 3; i++)    {
      if (sideChain.getLowerBound()[i] < lowerBound[i])
	  lowerBound[i] =  sideChain.getLowerBound()[i];
      if (sideChain.getUpperBound()[i] > upperBound[i])
	  upperBound[i] = sideChain.getUpperBound()[i];
      }
}
/**
* @Description  set bond length and angles to cristallographic values
 * @param none
 * @return   changes are made internally(void)
 */
void  AminoAcid::setDefault() {
  const double BOND_ANGLE_CA_TO_O = 120.80;
  const double BOND_LENGTH_C_TO_O = 1.231;
  const double BOND_LENGTH_CALPHA_TO_CPRIME = 1.52;
  const double BOND_LENGTH_N_TO_CALPHA = 1.458;
  const double BOND_ANGLE_AT_CALPHA_TO_CPRIME = 111.6;

  if ((*this)[N].isBond((*this)[CA]) && (*this)[CA].isBond((*this)[C]) && (*this)[C].isBond((*this)[O])) { 
      IntCoordConverter icc;

      icc.setBondLength((*this)[N], (*this)[CA], BOND_LENGTH_N_TO_CALPHA);
      icc.setBondLength((*this)[CA], (*this)[C], BOND_LENGTH_CALPHA_TO_CPRIME);
      icc.setBondLength((*this)[C], (*this)[O], BOND_LENGTH_C_TO_O);
      
      icc.setBondAngle((*this)[N], (*this)[CA], (*this)[C],DEG2RAD*BOND_ANGLE_AT_CALPHA_TO_CPRIME);
      icc.setBondAngle((*this)[CA], (*this)[C], (*this)[O],DEG2RAD*BOND_ANGLE_CA_TO_O);
    }
  else
    {cout << "not in bond\n";}
  sync();
}
/**
* @Description  Set the bonds based on lengths and angles given as params
 * @param length of NToCa , CaToC, CToO(double, double, double) , Angle of   atCaToC, CaToOAng(double, double)
 * @return   changes are made internally(void)
 */
void 
AminoAcid::setBonds(double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng){

  icc.setBondLength((*this)[N], (*this)[CA], NToCaLen);
  icc.setBondLength((*this)[CA], (*this)[C], CaToCLen);
  icc.setBondLength((*this)[C], (*this)[O], CToOLen);
  icc.setBondAngle((*this)[N], (*this)[CA], (*this)[C], atCaToCAng);
  icc.setBondAngle((*this)[CA], (*this)[C], (*this)[O], CaToOAng);
  sync();
}   
/**
* @Description  sets the bond structure for an aminoacid, copnsidering the pdb values
 * @param flag of connection (bool), amino acid  pointer(AminoAcid*) permissive flag( bool)
 * @return  flag to verify if the bond could be set (bool)
 */
bool 
AminoAcid::setBondsFromPdbCode(bool connect, AminoAcid* prev, bool permissive) {
  
  if ( (!isMember(C)) || (!isMember(CA)) || (!isMember(N))) {
      if (permissive)	{
	  return false;
	}
      else
	ERROR("Cannot set bonds for residues with missing atoms.", exception);
    }

  if (getType() == "X")
    patchAminoAcidCode();

  if (prev != NULL) {
      if (connect) 
	(*this)[N].setTrans( (*this)[N].getCoords() - (*prev)[C].getCoords() );
      bindIn((*this)[N], (*prev), (*prev)[C]);
    }
  (*this)[CA].bindStructure((*this)[N], connect);
  if (isMember(H)) //HN
    (*this)[H].bindStructure((*this)[N], connect);
 
   
  (*this)[C].bindStructure((*this)[CA], connect);
  if (isMember(HA))
    (*this)[HA].bindStructure((*this)[CA], connect);
     
  if (isMember(O))
    (*this)[O].bindStructure((*this)[C], connect);

  if (isMember(OXT))
    (*this)[OXT].bindStructure((*this)[C], connect);

  return sideChain.setBondsFromPdbCode(connect, permissive);
}
