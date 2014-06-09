/** 
 * @Class:              Loop extractor     
  *      
*/

#ifndef _LOOPEXTRACTOR_H_
#define _LOOPEXTRACTOR_H_
	
// Includes:
#include<Spacer.h>
#include<set>
#include<ranking_helper.h>
#include<ranking_helper2.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
 /** @brief  Extracts all the loop regions (by numbers) from a given spacer.

 * */
class LoopExtractor{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  LoopExtractor();
  LoopExtractor(Spacer* s);

 // PREDICATES:
  void nextLoop( int& start, int& end);
  double givePercentProp(multiset<ranking_helper> s, int count);
  double givePercentVdw(multiset<ranking_helper2> s, int count);
  void writeFile(vector<double> prop_percent, vector<double> vdw_percent, 
		 multiset<ranking_helper2> &rmsh, ofstream &statOut);

// MODIFIERS:
  void setSpacer(Spacer* s);
  void calcPropPercent(vector<double> &result, multiset<ranking_helper2> &rmsh,
		       multiset<ranking_helper> &rhs);
  void calcVdwPercent(vector<double> &result, multiset<ranking_helper2> &rmsh,
		       multiset<ranking_helper2> &rhs);
  static const int BEST_COUNT;  // determines the number of solutions which 
                                // are evaluated by propensities and collisions

// OPERATORS:

protected:

private:

// HELPERS:
int getPositionProp(multiset<ranking_helper> s, int index);
int getPositionVdw(multiset<ranking_helper2> s, int index);

// ATTRIBUTES:
  Spacer* sp;    // contains the spacer which is examined
  int position;  // current working position (from where we start 
                 // searching for the next loop)
};

} // Namespace
#endif
