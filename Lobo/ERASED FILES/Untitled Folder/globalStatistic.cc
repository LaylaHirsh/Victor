// -*- C++ -*------------------------------------------------------------------
//  $Id: globalStatistic.cc,v 1.2 2007-10-26 14:59:24 biocomp Exp $
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
//   Used to work with the data for the global statistic in the prop_calibration
//   program.
//
// -----------------x-------------------x-------------------x-----------------

// Includes:
#include <globalStatistic.h>


// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

const int globalStatistic::BAD_PROPENSITY = 0;
const int globalStatistic::BAD_VDW = 3000;
const int globalStatistic::BAD_CONSISTENCY = 1000;

// CONSTRUCTORS/DESTRUCTOR:
globalStatistic::globalStatistic() {

  solutionSum = 0;
  loopCount = 0;
  LoopExtractor le;
/*
  // initialize the arrays
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    for (int j = 0; j < le.BEST_COUNT; j += 1) {   
      for (int k = 0; k < 20; k += 1) {
	statArrayVdw[i][j][k] = 0;
	statArrayProp[i][j][k] = 0;
	statArrayEndRMS[i][j][k] = 0;
	statArrayCompactness[i][j][k] = 0;
	statArrayEnergy[i][j][k] = 0;
      }
    }
  }
*/
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    for (int j = 0; j < le.BEST_COUNT; j += 1) {   
      statArrayRMS[i][j] = 0;
      statArrayRMSCount[i][j] = 0;
      statArrayRankedRms[i][j] = 0;
      statArrayRankedRmsCount[i][j] = 0;
      deviationRms[i][j] = new vector<double>;
      deviationRankedRms[i][j] = new vector<double>;
    }
  }
  for(int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    for (int j = 0; j < 20; j += 1) {
//        statArrayBadProp[i][j] = 0;
//        statArrayBadVdw[i][j] = 0;
//        badConsistency[i][j] = 0;
//        vdwCutoff[i][j] = new vector<double>;
//        propCutoff[i][j] = new vector<double>;
//        energyCutoff[i][j] = new vector<double>;
//        compactnessCutoff[i][j] = new vector<double>;
//        endRMSCutoff[i][j] = new vector<double>;
      filterCutoff[i][j] = new vector<double>;
    }
  }
/*
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    for (int j = 0; j < 20; j += 1) {
      for (int k = 0; k < 20; k += 1) {
	statArraySidechain[i][j][k] = 0;
      }
    }
  }
*/

}

globalStatistic::~globalStatistic() {
  LoopExtractor le;
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    for (int j = 0; j < le.BEST_COUNT; j += 1) {   
      delete deviationRms[i][j];
      delete deviationRankedRms[i][j];
    }
  }
/*
  for(int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    for (int j = 0; j < 20; j += 1) {
      delete vdwCutoff[i][j];
      delete energyCutoff[i][j];
      delete propCutoff[i][j];
      delete compactnessCutoff[i][j];
      delete endRMSCutoff[i][j];
    }
  }
*/
}

// PREDICATES:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::outputArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  The arrays in this class are written to the file the filehandle refers to.
//  
// -----------------x-------------------x-------------------x-----------------

void
globalStatistic::outputArray(ofstream &outF) {

  LoopExtractor le;
//    outF << "Propensities" << "\n";
//    outF << "------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < le.BEST_COUNT; j += 1) {                 
//        outF << "Loesung Nr. " << j + 1 << "\n";
//        for (int k = 0; k < 20; k += 1) {
//  	outF << statArrayProp[i][j][k] << " ";
//        }
//        outF << "\n";
//      }
//    }
  
//    outF << "\n" << "Kollisionen" << "\n";
//    outF << "-----------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < le.BEST_COUNT; j += 1) {             
//        outF << "Loesung Nr. " << j + 1 << "\n";
//        for (int k = 0; k < 20; k += 1) {
//  	outF << statArrayVdw[i][j][k] << " ";
//        }
//        outF << "\n";
//      }
//    }


//    outF << "\n" << "EndRMS" << "\n";
//    outF << "------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < le.BEST_COUNT; j += 1) {         
//        outF << "Loesung Nr. " << j + 1 << "\n";
//        for (int k = 0; k < 20; k += 1) {
//  	outF << statArrayEndRMS[i][j][k] << " ";
//        }
//        outF << "\n";
//      }
//    }


//    outF << "\n" << "Compactness" << "\n";
//    outF << "-----------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < le.BEST_COUNT; j += 1) {               
//        outF << "Loesung Nr. " << j + 1 << "\n";
//        for (int k = 0; k < 20; k += 1) {
//  	outF << statArrayCompactness[i][j][k] << " ";
//        }
//        outF << "\n";
//      }
//    }

//    outF << "\n" << "Energy" << "\n";
//    outF << "------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < le.BEST_COUNT; j += 1) {     
//        outF << "Loesung Nr. " << j + 1 << "\n";
//        for (int k = 0; k < 20; k += 1) {
//  	outF << statArrayEnergy[i][j][k] << " ";
//        }
//        outF << "\n";
//      }
//    }


  outF << "\nRMS" << "\n---\n\n";
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) 
  {
    outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
    for (int j = 0; j < le.BEST_COUNT; j += 1) 
    {
      outF << "Durchschnittliche RMS fuer Loesung Nr. " << j + 1 << ": ";
      if (statArrayRMSCount[i][j] != 0) 
	  outF << ( statArrayRMS[i][j] / static_cast<double>(
	      statArrayRMSCount[i][j]) ) << "\n";
      else
	  outF << "kein Eintrag\n";
    }
  }

  outF << "\n" << "RMS der Loesungen nach dem Filterranking" << "\n";
  outF << "----------------------------------------" << "\n" << "\n";
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) 
  {
    outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
    for (int j = 0; j < 6; j += 1) 
    {
      outF << "Durchschnittliche RMS fuer Loesung Nr. " << j + 1 << ": ";
      if (statArrayRankedRmsCount[i][j] != 0) 
	outF << ( statArrayRankedRms[i][j] 
	    / static_cast<double>(statArrayRankedRmsCount[i][j]) ) << "\n";
      else 
	outF << "kein Eintrag\n";
    }
  }

//    outF << "\n" << "Bad prop values" << "\n";
//    outF << "---------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << statArrayBadProp[i][j] << " ";
//      }
//      outF << "\n";
//    }

//    outF << "\n" << "Bad vdw values" << "\n";
//    outF << "--------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << statArrayBadVdw[i][j] << " ";
//      }
//      outF << "\n";
//    }

//    outF << "\n" << "Good consistency values" << "\n";
//    outF <<"-----------------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << badConsistency[i][j] << " ";
//      }
//      outF << "\n";
//    }


//    outF << "\n" << "Best RMS of the n best energy values" << "\n";
//    outF << "---------------------------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << "Median for the " << j + 1 << " best values" << calcMedian(vdwCutoff[i][j]) << "\n";
//      }
//      outF << "\n";
//    }

//    outF << "\n" << "Best RMS of the n best propensity values" << "\n";
//    outF << "----------------------------------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << "Median for the " << j + 1 << " best values" << calcMedian(propCutoff[i][j]) << "\n";
//      }
//      outF << "\n";
//    }

//    outF << "\n" << "Best RMS of the n best compactness values" << "\n";
//    outF << "-----------------------------------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << "Median for the " << j + 1 << " best values" << calcMedian(compactnessCutoff[i][j]) << "\n";
//      }
//      outF << "\n";
//    }

//    outF << "\n" << "Best RMS of the n best EndRMS values" << "\n";
//    outF << "------------------------------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << "Median for the " << j + 1 << " best values" << calcMedian(endRMSCutoff[i][j]) << "\n";
//      }
//      outF << "\n";
//    }

//    outF << "\n" << "Best RMS of the n best Energy values" << "\n";
//    outF << "------------------------------------" << "\n" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {                                
//        outF << "Median for the " << j + 1 << " best values" << calcMedian(energyCutoff[i][j]) << "\n";
//      }
//      outF << "\n";
//    }

  outF << "\nBest RMSof the n best Filter ranked values\n"
       << "------------------------------------\n\n";
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) 
  {
    outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
    for (int j = 0; j < 20; j += 1) {                                
      outF << "Median for the " << j + 1 << " best values" 
	   << calcMedian(filterCutoff[i][j]) << "\n";
    }
    outF << "\n";
  }

  outF << "\nNumber of examined loops"
       << "\n------------------------\n";
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    outF << "Loop der Laenge: " << i + LOOP_MIN << " Anzahl: " 
	 << statArrayRMSCount[i][0] << "\n"; 
  }

  outF << "\n" << "As an average we have " << 
    static_cast<double>(solutionSum) / static_cast<double>(loopCount) << 
    " solutions remaining" << "\n";

  outF << "\n" << "Standard Deviation for the RMS" << "\n";
  outF << "------------------------------" << "\n";
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
    for (int j = 0; j < le.BEST_COUNT; j += 1) {
      outF << "Standard Deviation for solution nr. " << j + 1 << ": ";
      outF << calcDeviation(deviationRms[i][j]) << "\n";
    }
  }
  
  outF << "\n" << "Standard Deviation for the RankedRMS" << "\n";
  outF << "------------------------------------" << "\n";
  for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
    outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
    for (int j = 0; j < le.BEST_COUNT; j += 1) {
      outF << "Standard Deviation for solution nr. " << j + 1 << ": ";
      outF << calcDeviation(filterCutoff[i][j]) << "\n";
      //      outF << calcDeviation(deviationRankedRms[i][j]) << "\n";
    }
  }

//    outF << "\n" << "Number of sidechain collisions" << "\n";
//    outF <<"------------------------------" << "\n";
//    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
//      outF << "Loop der Laenge: " << i + LOOP_MIN << "\n"; 
//      for (int j = 0; j < 20; j += 1) {
//        for (int k = 0; k < 20; k += 1) {
//  	outF << statArraySidechain[i][j][k] << " ";
//        }
//        outF << "\n";
//      }
//    }

  outF.close();
}


/*

// MODIFIERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateVPArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  for each entry in the vectors, the both arrays (Van der Waals and 
//  Propensities) are updated (that means they are incremented by 1) at the 
//  proper position.
//
//  This routine depends on the fact, that the number of solutions generated
//  is less or equal to 100. Otherwise the last assertion may trigger!
//
// -----------------x-------------------x-------------------x-----------------
void
globalStatistic::updateVPArray(vector<double> &vdw_percent, vector<double> &prop_percent, int loopNr) {

  LoopExtractor le;
  ASSERT(prop_percent.size() == vdw_percent.size(), exception);
  ASSERT(le.BEST_COUNT >= prop_percent.size(), exception);
  for (unsigned int i = 0; i < prop_percent.size(); i += 1) {
    int vdwIndex = static_cast<int>(vdw_percent[i] * 100);
    if (vdwIndex == 100)
      vdwIndex = 99;
    vdwIndex -= 1;                   // changed in order to fix the array bug!
    vdwIndex /= 5;
    ASSERT(vdwIndex <= 19, exception);
    statArrayVdw[loopNr][i][vdwIndex] += 1;
    
    int propIndex = static_cast<int>(prop_percent[i] * 100);
    propIndex -= 1;                  // if we have more than 100 solutions 
    propIndex /= 5;                   // the variable could get negative
    ASSERT(propIndex <= 19 && propINdex >= 0, exception);
    statArrayProp[loopNr][i][propIndex] += 1;
  }	
};



// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateEndRMSArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  for each entry in the vector, the "statArrayEndRMS" is updated!
//  That means the array is incremented at the proper position!
//  
//  This routine depends on the fact, that the number of solutions generated
//  is less or equal to 100. Otherwise the last assertion may trigger!
//
// -----------------x-------------------x-------------------x-----------------
void
globalStatistic::updateEndRMSArray(vector<double> &endRMS_percent, int loopNr) {

  LoopExtractor le;
  ASSERT(le.BEST_COUNT >= endRMS_percent.size(), exception);
  for (unsigned int i = 0; i < endRMS_percent.size(); i += 1) {
    int endRMSIndex = static_cast<int>(endRMS_percent[i] * 100);
    endRMSIndex -= 1; 
    endRMSIndex /= 5;
    if (endRMSIndex > 19) {
      ERROR("endRMSIndex > 19!!! in globalStatistic::udateEndRMSArray", exception);
    }
    statArrayEndRMS[loopNr][i][endRMSIndex] += 1;
  }	
};


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateCompactnessArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  for each entry in the vector, the "statArrayCompactness" is updated!
//  That means the array is incremented at the proper position!
//  
//  This routine depends on the fact, that the number of solutions generated
//  is less or equal to 100. Otherwise the last assertion may trigger!
//
// -----------------x-------------------x-------------------x-----------------
void
globalStatistic::updateCompactnessArray(vector<double> &compactness_percent, int loopNr) {

  LoopExtractor le;
  ASSERT(le.BEST_COUNT >= compactness_percent.size(), exception);
  for (unsigned int i = 0; i < compactness_percent.size(); i += 1) {
    int compactnessIndex = static_cast<int>(compactness_percent[i] * 100);
    compactnessIndex -= 1; 
    compactnessIndex /= 5;
    if (compactnessIndex > 19) {
      ERROR("compactnessIndex > 19!!! in globalStatistic::udateCompactnessArray", exception);
    }
    statArrayCompactness[loopNr][i][compactnessIndex] += 1;
  }	
};


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateEnergyArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  for each entry in the vector, the "statArrayEnergy" is updated!
//  That means the array is incremented at the proper position!
//  
//  This routine depends on the fact, that the number of solutions generated
//  is less or equal to 100. Otherwise the last assertion may trigger!
//
// -----------------x-------------------x-------------------x-----------------
void
globalStatistic::updateEnergyArray(vector<double> &energy_percent, int loopNr) {

  LoopExtractor le;
  ASSERT(le.BEST_COUNT >= energy_percent.size(), exception);
  for (unsigned int i = 0; i < energy_percent.size(); i += 1) {
    int energyIndex = static_cast<int>(energy_percent[i] * 100);
    energyIndex -= 1; 
    energyIndex /= 5;
    if (energyIndex > 19) {
      ERROR("energyIndex > 19!!! in globalStatistic::udateEnergyArray", exception);
    }
    statArrayEnergy[loopNr][i][energyIndex] += 1;
  }	
};



// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateRMSArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  
//  This function iterates over the rmsh array until BEST_COUNT is reached.
//  In each loop the values of the RMS of the best BEST_COUNT solutions
//  are stored. Also the array which contains the number of entries is updated.
//
//  If we don't have le.BEST_COUNT solutions left, we fill the rest of the 
//  arrays with the worst rms value.
//
// -----------------x-------------------x-------------------x-----------------
*/
void 
globalStatistic::updateRMSArray( int loopNr, multiset<ranking_helper2> rmsh)
{
  int zaehler = 0;   
  set<ranking_helper2>::iterator pos;  // used to iterate over the rmsh
  LoopExtractor le;                    // needed in order to obtain BEST_COUNT
  double worst_value = 0;              // used to fill the rest of the arrays
 
  for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) 
  {
      if (zaehler >= le.BEST_COUNT)
	  break;   // we calculate these values only for the best solutions

      (deviationRms[loopNr][zaehler])->push_back(pos->get_value());
      statArrayRMS[loopNr][zaehler] += pos->get_value();
      statArrayRMSCount[loopNr][zaehler] += 1;
      worst_value = pos->get_value();
      zaehler += 1;
  }

  // now we check whether we have updated all le.BEST_COUNT solutions
  // if not, we do it be adding the worst of all solutions to it
  for (int i = zaehler; i < le.BEST_COUNT; i += 1) 
  {
      statArrayRMS[loopNr][i] += worst_value;
      (deviationRms[loopNr][i])->push_back(worst_value);
      statArrayRMSCount[loopNr][i] += 1;
  }
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateRankedRmsArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  
//  This function iterates over the rmsh array until BEST_COUNT is reached.
//  In each loop the values of the RMS of the best BEST_COUNT solutions
//  are stored. Also the array which contains the number of entries is updated
//.
//  If we have less than 6 surviving solutions, we take the worst solution
//  and add it to the remaining arrays.
//
//  Besides this, we update the statistic for the average number of surviving
//  solutions.
//
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::updateRankedRmsArray(int loopNr, vector<double> rms)
{
  if (rms.size() == 0)
      ERROR("We don't have a solution!", exception);

  unsigned int count = rms.size();
  double worstValue = 0;

  if (rms.size() >= 7) 
    count = 6;

  // now find a starting value for the worst value in the rms vector
  worstValue = rms[0];

  for (unsigned int i = 0; i < count; i += 1) 
  {
    statArrayRankedRms[loopNr][i] += rms[i];
    (deviationRankedRms[loopNr][i])->push_back(rms[i]);
    statArrayRankedRmsCount[loopNr][i] += 1;

    if (worstValue < rms[i]) 
      worstValue = rms[i];
  }

  // fill rest of array in case less than 6 solutions left
  for (unsigned int i = count; i < 6; i += 1) 
  {
    statArrayRankedRms[loopNr][i] += worstValue;
    (deviationRankedRms[loopNr][i])->push_back(worstValue);
    statArrayRankedRmsCount[loopNr][i] += 1;
  }

  // update statistics for average number of surviving solutions
  loopCount += 1;
  solutionSum += rms.size();
}



// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateEndRMSArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  
//  This function iterates over the rmsh array until BEST_COUNT is reached.
//  In each loop the values of the RMS of the best BEST_COUNT solutions
//  are stored. Also the array which contains the number of entries is updated.
//
// -----------------x-------------------x-------------------x-----------------

/*void 
globalStatistic::updateEndRMSArray( int loopNr, multiset<ranking_helper2> rmsh){

  int zaehler = 0;                             // used in a loop in order to calculate the propensities 
                                               // for the lowest rms solutions
  set<ranking_helper2>::iterator pos;          // used to iterate over the rmsh
  LoopExtractor le;                            // needed in order to obtain BEST_COUNT
 
  ASSERT(loopNr < LOOP_MAX - LOOP_MIN + 1, exception);
  for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) {
    if (zaehler >= le.BEST_COUNT) {
      break;                                              // we calculate these values only for the best solutions
    }
   statArrayEndRMS[loopNr][zaehler] += pos->get_value();
   statArrayEndRMSCount[loopNr][zaehler] += 1;
   zaehler += 1;
  }
  };*/

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateBadValuesArray()
//
//  Author:        Achim Trabold
//
//  Date:          09/00
//
//  Description:
//  
//  The both arrays statArrayBadProp and statArrayBadVdw are filled. These
//  arrays contain a statistic for each separate loop-length in which the matching
//  position of the bad value (that means a value which is equal or worse to 
//  BAD_PROPENSITY / BAD_VDW) is incremented by 1. The matching position means
//  here the position of the corresponding solution in the RMS ranking.
//
// -----------------x-------------------x-------------------x-----------------
/*
void 
globalStatistic::updateBadValuesArray( int loopNr, multiset<ranking_helper> prop_ranking, 
				       multiset<ranking_helper2> vdw_ranking, 
				       multiset<ranking_helper2> rmsh) {

  set<ranking_helper>::iterator pos;          // used to iterate over the prop_ranking set
  set<ranking_helper2>::iterator pos2;        // used to iterate over the vdw_ranking set

  //first we do the propensities
  for (pos = prop_ranking.begin(); pos != prop_ranking.end(); ++pos) {
    if ( pos->get_value() <= BAD_PROPENSITY) {
      int percent = calcRmsPercent(pos->get_index() - 1, rmsh);
      //      cout << "Kontrolle fuer calcRMSPercent Prop! Index: " << pos->get_index() - 1 << 
      //	 " percent: " << percent << "\n";
      int index = percent / 5;
      ASSERT (index >= 0 && index <= 20, exception);
      statArrayBadProp[loopNr][index] += 1;
    }
  }
  
  // now we do the vdw
  for (pos2 = vdw_ranking.begin(); pos2 != vdw_ranking.end(); ++pos2) {
    if (pos2->get_value() >= BAD_VDW) {
      int percent = calcRmsPercent(pos2->get_index() - 1, rmsh);
      int index = percent / 5;
      //      cout << "Kontrolle fuer calcRMSPercent VDW! Index: " << pos2->get_index() - 1<<
      //	 " percent: " << percent << "\n";
      ASSERT (index >= 0 && index <= 20, exception);
      statArrayBadVdw[loopNr][index] += 1;
    }
  }
}



// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateBadConsistencyArray()
//
//  Author:        Achim Trabold
//
//  Date:          09/00
//
//  Description:
//  
//  The arry "badConsistency" is updated. That means each solution which has
//  a bad consistency value has a ranking in the delivered sorted set of the rms 
//  (the variable rmsh). The position of the solution with a bad consistency
//  value in this sorted set is marked in the "badConsistency" arry.
//  
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::updateBadConsistencyArray( int loopNr, 
					    vector<int> consistency, 
					    multiset<ranking_helper2> rmsh) {

  for (unsigned int i = 0; i < consistency.size(); i += 1) {
    if ( consistency[i] < BAD_CONSISTENCY) {
      int percent = calcRmsPercent(i, rmsh);
      cout << "Kontrolle fuer badConsistency Index: " << i << 
	" percent: " << percent << "\n";
      int index = percent / 5;
      if (index < 0 || index >= 20) {
	ERROR(" index is out of bound!!!", exception);
      }
      badConsistency[loopNr][index] += 1;
    }
  }
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::VdwCutoffGenerator()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  This routine checks the first count values in vdw_percent for the 
//  minimum. This value is added into the proper position of the
//  vdwCutoff array.  
//  
//  If count is bigger than the size of vdw_percent we only check to the latter.
//  We expect count to start from 1
//
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::VdwCutoffGenerator(multiset<ranking_helper2> &rmsh,
				    multiset<ranking_helper2> &rhs, int count, int loopNr){

  int zaehler = 1;                      // used to count up to count
  double minRms;                        // contains the minimum found RMS
  double helper;                        // used to find the minimum

  set<ranking_helper2>::iterator pos;   // used to iterate over rhs
  
  // initialize minRms
  pos = rhs.begin();
  minRms = giveRms((*pos).get_index() - 1, rmsh);

  for (pos = rhs.begin(); pos != rhs.end(); ++pos) {
    if (zaehler > count) {
      break;
    }
    zaehler += 1;
    
    helper = giveRms((*pos).get_index() - 1, rmsh);
    if (minRms > helper) {
      minRms = helper;
    }
  }
  (vdwCutoff[loopNr][count-1])->push_back(minRms);
}



// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::PropCutoffGenerator()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  This routine checks the first count values in prop_percent for the 
//  minimum. This value is added into the proper position of the
//  propCutoff array.  
//  
//  If count is bigger than the size of prop_percent the latter is used.
//
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::PropCutoffGenerator(multiset<ranking_helper2> &rmsh,
				     multiset<ranking_helper> &rhs, int count, int loopNr){
  
  int zaehler = 1;                      // used to count up to count
  double minRms;                        // contains the minimum found RMS
  double helper;                        // used to find the minimum

  set<ranking_helper>::iterator pos;    // used to iterate over rhs
  
  // initialize minRms
  pos = rhs.begin();
  minRms = giveRms((*pos).get_index() - 1, rmsh);

  for (pos = rhs.begin(); pos != rhs.end(); ++pos) {
    if (zaehler > count) {
      break;
    }
    zaehler += 1;

    helper = giveRms((*pos).get_index() - 1, rmsh);
    if (minRms > helper) {
      minRms = helper;
    }
  }
  (propCutoff[loopNr][count-1])->push_back(minRms);
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::CompactnessCutoffGenerator()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  This routine checks the first count values in compactness_percent for the 
//  minimum. This value is added into the proper position of the
//  compactnessCutoff array.  
//  
//  If count is bigger than the size of compactness_percent the latter is used.
//
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::CompactnessCutoffGenerator(multiset<ranking_helper2> &rmsh,
					    multiset<ranking_helper2> &rhs, int count, int loopNr){
  
  int zaehler = 1;                      // used to count up to count
  double minRms;                        // contains the minimum found RMS
  double helper;                        // used to find the minimum

  set<ranking_helper2>::iterator pos;   // used to iterate over rhs
  
  // initialize minRms
  pos = rhs.begin();
  minRms = giveRms((*pos).get_index() - 1, rmsh);

  for (pos = rhs.begin(); pos != rhs.end(); ++pos) {
    if (zaehler > count) {
      break;
    }
    zaehler += 1;

    helper = giveRms((*pos).get_index() - 1, rmsh);
    if (minRms > helper) {
      minRms = helper;
    }
  }
  (compactnessCutoff[loopNr][count-1])->push_back(minRms);
}




// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::EnergyCutoffGenerator()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  This routine checks the first count values in energy_percent for the 
//  minimum. This value is added into the proper position of the
//  energyCutoff array.  
//  
//  If count is bigger than the size of energy_percent the latter is used.
//
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::EnergyCutoffGenerator(multiset<ranking_helper2> &rmsh,
				       multiset<ranking_helper2> &rhs, int count, int loopNr){
 
  int zaehler = 1;                      // used to count up to count
  double minRms;                        // contains the minimum found RMS
  double helper;                        // used to find the minimum

  set<ranking_helper2>::iterator pos;   // used to iterate over rhs
  
  // initialize minRms
  pos = rhs.begin();
  minRms = giveRms((*pos).get_index() - 1, rmsh);

  for (pos = rhs.begin(); pos != rhs.end(); ++pos) {
    if (zaehler > count) {
      break;
    }
    zaehler += 1;

    helper = giveRms((*pos).get_index() - 1, rmsh);
    if (minRms > helper) {
      minRms = helper;
    }
  }
  (energyCutoff[loopNr][count-1])->push_back(minRms); 
}



// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::EndRMSCutoffGenerator()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description: 
//  This routine checks the first count values in endRMS_percent for the 
//  minimum. This value is added into the proper position of the
//  endRMSCutoff array.  
//  
//  If count is bigger than the size of endRMS_percent, the latter is used..
//
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::EndRMSCutoffGenerator(multiset<ranking_helper2> &rmsh,
				       multiset<ranking_helper2> &rhs, int count, int loopNr){
 
  int zaehler = 1;                      // used to count up to count
  double minRms;                        // contains the minimum found RMS
  double helper;                        // used to find the minimum

  set<ranking_helper2>::iterator pos;   // used to iterate over rhs
  
  // initialize minRms
  pos = rhs.begin();
  minRms = giveRms((*pos).get_index() - 1, rmsh);

  for (pos = rhs.begin(); pos != rhs.end(); ++pos) {
    if (zaehler > count) {
      break;
    }
    zaehler += 1;
    
    helper = giveRms((*pos).get_index() - 1, rmsh);
    if (minRms > helper) {
      minRms = helper;
    }
  }
  (endRMSCutoff[loopNr][count-1])->push_back(minRms);
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::FilterCutoffGenerator()
//
//  Author:        Achim Trabold
//
//  Date:          11/00
//
//  Description: 
//  We check the first count values in the vector values for the minimum.
//  This minimum is put into the proper array.
//
// -----------------x-------------------x-------------------x-----------------
*/
void 
globalStatistic::FilterCutoffGenerator(vector<double> values, int count, 
int loopNr)
{
  if ((count > 20) || (values.size() == 0))
	return;

  double min = values[0];

  for (int i = 1; i < count; i += 1) 
      if (values[i] < min)
	  min = values[i];

  (filterCutoff[loopNr][count-1])->push_back(min);
}

/*
// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::updateSidechainArray()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description: 
//  
//  Here we fill the array statArraySidechain with values. For each solution 
//  (up to 20 solutions) we mark how many collisions the solution had (we 
//  permit up to 20 collisions).
//
//  The vector collisionCount needs to be filled with the collisions for each
//  solution and rmsh must contain the sorted set of rms.
//
// -----------------x-------------------x-------------------x-----------------

void 
globalStatistic::updateSidechainArray(vector<int> collisonCount, multiset<ranking_helper2> rmsh, 
				      int loopNr) {

   unsigned int zaehler = 0;                   // used in a loop in order to count to a max. of 20
  set<ranking_helper2>::iterator pos;          // used to iterate over the rmsh

  for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) {

    if (zaehler > 19) {
      break;
    }
    
    if (collisonCount.size() <= (*pos).get_index()) {
      ERROR("Zaehler wurde zu gross!!", exception);
    }
    int helper = collisonCount[(*pos).get_index()];
    if (helper >= 20 || helper == -1) {
      helper = 19;
    }
    statArraySidechain[loopNr][zaehler][helper] += 1;
    zaehler += 1;
  }
}

// OPERATORS:

// HELPERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::calcRMSPercent()
//
//  Author:        Achim Trabold
//
//  Date:          08/00
//
//  Description:
//  
//  This routine gives the percent-value of the index in the rmsh set back.
//  Percent value means, at which position in the rmsh set the entry with 
//  index is located where percent value of 0 means it is at the first position.
//  The return value runs from 0 to 99.
//
// -----------------x-------------------x-------------------x-----------------

int 
globalStatistic::calcRmsPercent(int index, multiset<ranking_helper2> rmsh){

  int counter = 1;                             // contains the position of the searched index
  int size = 0;                                // contains the size of the rmsh set
  bool found = false;                          // indicates whether we have found the index
  set<ranking_helper2>::iterator pos;          // used to iterate over rmsh

  for(pos = rmsh.begin(); pos != rmsh.end(); ++pos) {
    if (pos->get_index() == index) {
      found = true;
    }
    if (found == false) {
      counter += 1;
    }
    size += 1;
  }
  ASSERT (found == false, exception);          // otherwise something has terribly gone wrong
  return ( (counter * 100 / size) - 1);
}


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::calcDeviation()
//
//  Author:        Achim Trabold
//
//  Date:          10/00
//
//  Description:
//  
//  Calculates the standard deviation of the values in the given Vector
//
// -----------------x-------------------x-------------------x-----------------
*/
double
globalStatistic::calcDeviation(vector<double> *values) {
  
  double median = 0;                         
  double deviation = 0;

  // first we have to calculate the median
  for (unsigned int i = 0; i < values->size(); i += 1) {
    median += (*values)[i];
  }
  median = median / static_cast<double>(values->size());
  
  // with the median we can calculate the standard deviation
  for (unsigned int i = 0; i < values->size(); i += 1) {
    double helper = (*values)[i] - median;
    helper = helper * helper;
    deviation += helper;
  }
  deviation /= static_cast<double>(values->size() - 1);
  deviation = sqrt(deviation);
  return deviation;
}
/*

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::calcMedian()
//
//  Author:        Achim Trabold
//
//  Date:          10/00
//
//  Description:
//  
//  Calculates the median of the values in the given Vector
//
// -----------------x-------------------x-------------------x-----------------

*/
double
globalStatistic::calcMedian(vector<double> *values) 
{
  if (values->size() == 0) 
      return 0;

  double median = 0;                         
  for (unsigned int i = 0; i < values->size(); i += 1)
      median += (*values)[i];

  return ( median / static_cast<double>(values->size()) );
}
/*


// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        globalStatistic::giveRms()
//
//  Author:        Achim Trabold
//
//  Date:          11/00
//
//  Description:
//  
//  Extracts the rms value with index in rmsh
//
// -----------------x-------------------x-------------------x-----------------

double globalStatistic::giveRms(int index, multiset<ranking_helper2> &rmsh) {

  set<ranking_helper2>:: iterator pos;

  for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) {
    if ( (*pos).get_index() == index) {
      return (*pos).get_value();
    }
  }
  // if we reach this point index wasn't in rmsh thus we have an error
  ERROR("index is not in rmsh", exception);
}

*/

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        <Name>
//
//  Class:         <ClassName>
//  
//  Author:        <Name> 
//
//  Project Name:  Victor
//
//  Date:          mm/yy
//
//  Reviewed By:   <Name>
//
//  Description:
//    <Description>
//
// -----------------x-------------------x-------------------x-----------------
