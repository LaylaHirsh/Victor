// -*- C++ -*-----------------------------------------------------------------
//  $Id: globalStatistic.h,v 1.2 2007-12-17 09:23:14 biocomp Exp $
//
//  Class:              globalStatistic
//
//  Base Class(es):     -
//
//  Derived Class(es):  -
//
//  Containing:         -
//
//  Author:             Achim Trabold
//
//  Project Name:       Victor
//
//  Date:               08/00
//
//  Reviewed By:        <Name>
//
//  Description:
//   Used to work with the data for the global statistic in the 
//   prop_calibration program.
//
// -----------------x-------------------x-------------------x-----------------

#ifndef _globalStatistic_h_
#define _globalStatistic_h_
	
// Includes:
#include "LoopExtractor.h"
#include <set>
#include "ranking_helper2.h"

// Global constants, typedefs, etc. (to avoid):
namespace Biopool {

class globalStatistic
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  globalStatistic();
  ~globalStatistic();

// PREDICATES:
  void outputArray(ofstream &outF);

  static const int LOOP_MIN = 3;        // lower limit of examined Loops
  static const int LOOP_MAX = 13;       // upper limit of examined loops

// MODIFIERS:
  void updateVPArray(vector<double> &vdw_percent, 
	 vector<double> &prop_percent, int loopNr);
  void updateRMSArray( int loopNr, multiset<ranking_helper2> rmsh);
  void updateBadValuesArray( int loopNr, 
	 multiset<ranking_helper> prop_ranking, 
	 multiset<ranking_helper2> vdw_ranking, 
	 multiset<ranking_helper2> rmsh);
  void updateBadConsistencyArray( int loopNr, vector<int> consistency, 
				  multiset<ranking_helper2> rmsh);

  void updateEndRMSArray(vector<double> &endRMS_percent, int loopNr);
  void updateCompactnessArray(vector<double> &compactness_percent, 
			      int loopNr);
  void updateEnergyArray(vector<double> &energy_percent, int loopNr);

  void VdwCutoffGenerator(multiset<ranking_helper2> &rmsh,
			  multiset<ranking_helper2> &rhs, int count, 
			  int loopNr);
  void PropCutoffGenerator(multiset<ranking_helper2> &rmsh,
			   multiset<ranking_helper> &rhs, int count, 
			   int loopNr);
  void CompactnessCutoffGenerator(multiset<ranking_helper2> &rmsh,
				  multiset<ranking_helper2> &rhs, int count, 
				  int loopNr);
  void EndRMSCutoffGenerator(multiset<ranking_helper2> &rmsh,
			     multiset<ranking_helper2> &rhs, int count, 
			     int loopNr);
  void EnergyCutoffGenerator(multiset<ranking_helper2> &rmsh,
			     multiset<ranking_helper2> &rhs, int count, 
			     int loopNr);
  void FilterCutoffGenerator(vector<double> values, int count, int loopNr);
  void updateRankedRmsArray(int loopNr, vector<double> rms);
  void updateSidechainArray(vector<int> collisionCount, 
			    multiset<ranking_helper2> rmsh, int loopNr);

  //  void updateEndRMSArray( int loopNr, multiset<ranking_helper2> rmsh);
  
  // OPERATORS:
  
protected:

private:

  // HELPERS:
  int calcRmsPercent(int index, multiset<ranking_helper2> rmsh);
  double calcDeviation(vector<double> *values);
  double calcMedian(vector<double> *values);
  double giveRms(int index, multiset<ranking_helper2> &rmsh);
  
  // ATTRIBUTES:
  
  static const int BAD_PROPENSITY;  // value from which we have
                                 // a value to put into the bad prop statistic
  static const int BAD_VDW;         
    // contains the value from which we have a value to put into the bad 
    // VDW statistic
  static const int BAD_CONSISTENCY;  
    // contains the start value from which we put
    // a consistency value into the badConsistency array
/*
  int badConsistency[LOOP_MAX - LOOP_MIN + 1][20];
    // contains the bad consistency statistic
  int statArrayBadProp[LOOP_MAX - LOOP_MIN + 1][20];
    // contains the bad propensity statistic
  vector<double> * vdwCutoff[LOOP_MAX - LOOP_MIN + 1][20]; 
    // used to calculate the vdw cutoff value 
  vector<double> * energyCutoff[LOOP_MAX - LOOP_MIN + 1][20];
    // used to calculate the energy cutoff value 
  vector<double> * propCutoff[LOOP_MAX - LOOP_MIN + 1][20]; 
    // used to calculate the prop cutoff value 
  vector<double> * compactnessCutoff[LOOP_MAX - LOOP_MIN + 1][20];  
    // used to calculate the compactness cutoff value 
  vector<double> * endRMSCutoff[LOOP_MAX - LOOP_MIN + 1][20]; 
    // used to calculate the endRMS cutoff value 
*/
  vector<double> * filterCutoff[LOOP_MAX - LOOP_MIN + 1][20];
    // used to calculate the filter-combination cutoff 
/*
  int statArrayBadVdw[LOOP_MAX - LOOP_MIN + 1][20];
    // contains the bad vdw statistic
  int statArrayVdw[LOOP_MAX - LOOP_MIN + 1][6][20];
    // contains the overall solution statistic   @ 6!!
  int statArrayEndRMS[LOOP_MAX - LOOP_MIN + 1][6][20];
    // contains the overall solution statistic
  int statArrayCompactness[LOOP_MAX - LOOP_MIN + 1][6][20];
    // contains the overall solution statistic
  int statArrayProp[LOOP_MAX - LOOP_MIN + 1][6][20];
    // contains the overall solution statistic   @ 6!!
*/
  double statArrayRMS[LOOP_MAX - LOOP_MIN + 1][6];
    // contains the statistic for the RMS @ 6!!
  int statArrayRMSCount[LOOP_MAX - LOOP_MIN + 1][6];
    // counts the entries in statArrayRMS @ 6!!
/*
  int statArrayEnergy[LOOP_MAX - LOOP_MIN + 1][6][20];  
    // contains the overall solution statistic   @ 6!!
*/
  double statArrayRankedRms[LOOP_MAX - LOOP_MIN + 1][6];
    // contains the rms values of the surviving solutions 
    // ranked by the filter results
  int statArrayRankedRmsCount[LOOP_MAX - LOOP_MIN + 1][6];
    // contains the number of entries in statArrayRankedRms

  long solutionSum;    // contains the sum of all surviving solutions
  int loopCount;       // contains the number of loops we have examined

  vector<double> *deviationRms[LOOP_MAX - LOOP_MIN + 1][6];
    // used to calculate the standard deviation for the RMS
  vector<double> *deviationRankedRms[LOOP_MAX - LOOP_MIN + 1][6];
    // used to calculate the standard deviation for the ranked RMS
/*
  int statArraySidechain[LOOP_MAX - LOOP_MIN + 1][20][20];
    // used to determine how many collisions the sidechain has 
*/

  /*  double statArrayEndRMS[LOOP_MAX - LOOP_MIN + 1]
  [6];                                            // contains the statistic for the EndRMS @ 6!!
  int statArrayEndRMSCount[LOOP_MAX - LOOP_MIN + 1]
  [6]; */                                            // counts the entries in statArrayEndRMS @ 6!!

};

// ---------------------------------------------------------------------------
//                                 globalStatistic
// -----------------x-------------------x-------------------x-----------------

}
#endif
