/** 
 * @Class:              ranking_helper2
 * 
 * @Description:Just contains an int and double value which gets sorted (by a STL set)
*    in order to determine the ranking of low rms solutions with special 
*    filters
  *      
*/

#ifndef _RANK_HELPER2_H_
#define _RANK_HELPER2_H_
	
// Includes:

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
  /** @brief the class contains methods to manage the ranking
 * 
* @Description  Just contains an int and double value which gets sorted (by a STL set)
*    in order to determine the ranking of low rms solutions with special 
*    filters
 * */
class ranking_helper2 {
public:
  
  ranking_helper2(int ind, double val);
  ranking_helper2(const ranking_helper2& c);
  int get_index() const;
  double get_value() const;
  bool operator< (const ranking_helper2 &name) const;
  ranking_helper2& operator=(const ranking_helper2& orig);
  void copy(const ranking_helper2& c);

private:  
  int index;                                     
  double value;                                  
};
/**
 * @Description returns the rms ranking of the solution
 * @param  none
 * @return  the corresponding value
*/
 inline int ranking_helper2::get_index() const {
   return index;
 }
/**
 * @Description returns a value from a filter (like the propensity or the collision)
 * @param  none
 * @return  the corresponding value
*/
 inline double ranking_helper2::get_value() const {
   return value;
 }
 

} // Namespace
#endif
