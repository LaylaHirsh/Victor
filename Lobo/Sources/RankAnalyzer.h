/** 
 * @Class:              RankAnalyzer
 */

#ifndef _RANKANALYZER_H_
#define _RANKANALYZER_H_

// Includes:
#include <string>
#include <vector>
using namespace std;

namespace Biopool {
  
// Global constants, typedefs, etc. (to avoid):
 /** @brief  This class implements functions to analyze the ranking output of lobo.
 *  
 * */
class RankAnalyzer{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  RankAnalyzer();
  RankAnalyzer(const RankAnalyzer& orig); 
  virtual ~RankAnalyzer();

// PREDICATES:
  void printHeader();
  void printCorrelation(char* intro, unsigned int index);
  void printCorrelation(char* intro, vector<double> data);
  void printTopXResults(unsigned int top, string topFile, 
			double maxScore = 1000.0);

// MODIFIERS:
  void load(istream& inFile, unsigned int select);
  void calcStatistics();

  void copy(const RankAnalyzer& orig);

// OPERATORS:

//  protected:

//  private:

// HELPERS: 

// ATTRIBUTES:
  static const unsigned int MAX_COL = 11;
  static const unsigned int MAX_LOOP = 1000;

  double avg[MAX_COL];
  double sd[MAX_COL];
  vector< vector<double> > zScore;   // Z-scores of the data ("col")
  vector< vector<double> > col;      // raw input data (i.e. "column x")
  vector< vector<double> > result;   // RMSD result, used to get Top X 
  vector< vector<double> > resScore; // score of the RMSD results
};
 
} // namespace
#endif //_RANKANALYZER_H_


