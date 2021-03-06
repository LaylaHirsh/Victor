/** 
  * 
 * @Class:           ranking_helper   
  * 
 * @Description:Just contains an int and double value which gets sorted (by a STL set)
*    in order to determine the ranking of low rms solutions with special 
*    filters
  *      
*/

//Includes
#include<ranking_helper.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:
/**
* @Description constructor  sets the values of the index that contains the rms ranking of the solution , 
* and the value from a filter (like the propensity or the collision)
* @param index(int) , value(double)
*/
ranking_helper::ranking_helper(int ind, double val) {
    index = ind;
    value = val;
}
/**
* @Description Constructor that copies the info from another object
* @param  reference to original object(const ranking_helper& )
*/
ranking_helper::ranking_helper(const ranking_helper& c) {
  this->copy(c);
}

// PREDICATES:

// MODIFIERS:

// OPERATORS:
/**
* @Description Operator that allows to verify if its lower than other
* @param  reference to the object(const ranking_helper &)
* @return  result of the verification
*/
bool ranking_helper::operator< (const ranking_helper &name) const {
  return value > name.get_value();
}
/**
* @Description Operator that allows to assign one object to other
* @param  reference to the object(const ranking_helper &)
* @return  result of the verification
*/
ranking_helper& ranking_helper::operator=(const ranking_helper& orig) {
  if(&orig != this) {
    copy(orig);
  }
  return *this;
}
/**
* @Description copies the info from another object
* @param  reference to original object(const ranking_helper& )
* @return changes are made internally(void)
*/
void ranking_helper::copy(const ranking_helper& c) { 
  index = c.index;
  value = c.value;  
}

// HELPERS:
