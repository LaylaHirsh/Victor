/** 
 * @Class:              ranking_helper
  *      
*/

#ifndef _RANK_HELPER_H_
#define _RANK_HELPER_H_
	
// Includes:

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
  /** @brief the class contains methods to manage the ranking
 * 
* @Description  Just contains an int and double value which gets sorted (by a STL set)
*    in order to determine the ranking of low rms solutions with special 
*    filters
 * */
class ranking_helper {
public:
  
  ranking_helper(int ind, double val);
  ranking_helper(const ranking_helper& c);
  int get_index() const;
  double get_value() const;
  bool operator< (const ranking_helper &name) const;
  ranking_helper& operator=(const ranking_helper& orig);
  void copy(const ranking_helper& c);

private:  
  int index;           // contains the rms ranking of the solution
  double value;        // contains a value from a filter (like the propensity or the collision)
};
/**
 * @Description returns the rms ranking of the solution
 * @param  none
 * @return  the corresponding value
*/
 inline int ranking_helper::get_index() const {
   return index;
 }
/**
 * @Description returns a value from a filter (like the propensity or the collision)
 * @param  none
 * @return  the corresponding value
*/
 inline double ranking_helper::get_value() const {
   return value;
 }
 

} // Namespace
#endif
