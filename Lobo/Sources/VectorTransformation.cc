/** 
 * @Class:              VectorTransformation
 * @Description:This class allows to store transformation steps for transforming a 
*     vector v into v' according to the series of steps performed earlier.
  *      
*/

// Includes:
#include <VectorTransformation.h>
#include <IoTools.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:
/**
* @Description basic constructor
*/
VectorTransformation::VectorTransformation() : rot(0), trans(0){
  addNewElem();
}

/**
* @Description Constructor that copies a original vector transformation 
*/
VectorTransformation::VectorTransformation(const VectorTransformation& orig){
 this->copy(orig);
}
/**
* @Description basic destructor
*/
VectorTransformation::~VectorTransformation(){
  PRINT_NAME;
  rot.clear();
  trans.clear();
}


// PREDICATES:
/**
* @Description transform the original vector
* @param  original vector(vgVector3<float> )
* @return  corresponding value after the transformation(vgVector3<float> )
*/
vgVector3<float> VectorTransformation::transform(vgVector3<float> orig){
  PRECOND(rot.size() == trans.size(), exception);
  
  for (unsigned int i = rot.size(); i > 0; i--)    {
      orig = rot[i-1] * orig;
      orig += trans[i-1];
    }

  return orig;
}



// MODIFIERS:
/**
* @Description Adds rotation 
* @param   matrix to rotate(vgMatrix3<float>)
* @return  changes are made internally(void)
*/
void VectorTransformation::addRot(vgMatrix3<float> rm){
  PRECOND(rot.size() == trans.size(), exception);

  if (trans[trans.size()-1].length() != 0)
    addNewElem();

  rot[rot.size()-1] = rot[rot.size()-1] * rm;

}


/**
* @Description Adds translation 
* @param   matrix to rotate(vgVector3<float>)
* @return  changes are made internally(void)
*/
void VectorTransformation::addTrans(vgVector3<float> t){
  PRECOND(rot.size() == trans.size(), exception);

  vgMatrix3<float> idM(1);
  if (rot[rot.size()-1] != idM)
    addNewElem();

  trans[trans.size()-1] += t;

}

/**
* @Description Clears the object data
* @param  none
* @return  changes are made internally(void)
*/
void VectorTransformation::clear(){
  rot.clear();
  trans.clear();
  addNewElem();
}

/**
* @Description Copies the original vector transformation
* @param  original vector transformation
* @return   changes are made internally(void)
*/
void VectorTransformation::copy(const VectorTransformation& orig){
  PRINT_NAME; 
  rot.clear();
  trans.clear();

  for (unsigned int i = 0; i < orig.rot.size(); i++)
    rot.push_back(orig.rot[i]);
  for (unsigned int i = 0; i < orig.trans.size(); i++)
    trans.push_back(orig.trans[i]);
}


// OPERATORS:
/**
* @Description Assigns a transformation vector into another one.
* @param  reference to the original transformation vector(VectorTransformation& )
* @return  reference to the new transformation vector(VectorTransformation& )
 */
VectorTransformation& VectorTransformation::operator=(const VectorTransformation& orig){
  if (&orig != this)
    copy(orig);
  return *this;
}

