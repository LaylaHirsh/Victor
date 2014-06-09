/**
 * @Description This program creates filter ranked solutions for a single loop.
 */
#include <string>
#include <GetArg.h>
#include <LoopModel.h>
#include <LoopTableEntry.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IntCoordConverter.h>
#include <VectorTransformation.h>
#include <LoopExtractor.h>
#include <set>
#include <ranking_helper.h>
#include <ranking_helper2.h>
#include <String2Number.h>
#include <globalStatistic.h>

#include <Cluster.h> 
#include <ClusterEntries.h>
#include <evaluateLoopSidechainCapability.h>


using namespace Biopool;


int main(int nArgs, char* argv[]){ 
  cout << "Start\n";
  
  multiset<ranking_helper> prop_ranking;   // used to sort the propensities
  multiset<ranking_helper2> rmsh;          // used to sort the rms values
  multiset<ranking_helper2> vdw_ranking;   // used to sort the vdw terms
  multiset<ranking_helper2> endRMS_ranking;// used to sort the endRMS values
  multiset<ranking_helper2> cmpct_ranking; // used to sort the compactness values
  multiset<ranking_helper2> energy_ranking;// used to sort the energy values
  
  unsigned int index1 = 0;                 // contains the start index of the loop
  unsigned int index2 = 0;                 // contains the end index of the loop
  unsigned int num;                        // determines how many solutions we are getting in the first level of the loop modeling
  unsigned int num2;                       // determines how many solutions we are getting in the second level of the loop modeling

  LoopModel lm;                            // here we create the loop
  globalStatistic gs;                      // contains the overall statistic for the whole thing (e.g. propensities) 

 
  static const int criticalConsistency = 1000; // all solutions with a worse consistency are rejected

  vector<int> v1, v2;
  Spacer sp2;
  evaluateLoopSidechainCapability elc(&sp2, v1, v2, index1, index2);

  if (getArg( "h", nArgs, argv))    {
      cout << "LoopModelFilter\n"
	   << "This program creates filter ranked solutions for a single loop.\n"
	   << " Options: \n"
	   << "\t-i <Input pdb-file> \t\t pdb-file of the protein with the loop\n"
	   << "\t-a <n-start> \t\t\t Aminoacid where loop starts\n" 
	   << "\t-b <n-end> \t\t\t Aminoacid where loop ends\n"
	   << "\t-v <Weight VDW> \t\t Weight of the vdw filter (0 is standard)\n"
	   << "\t-e <Weight Energy> \t\t Weight of the energy filter(1 is standard)\n"
	   << "\t-r <Weight EndRMS> \t\t Weight of the end-rms filter(80 is standard)\n"
	   << "\t-c <Weight compactness> \t Weight of the compactness filter(0 is standard)\n"
	   << "\t-p <Weight propensity> \t\t Weight of the propensity filter(0 is standard) \n"
	   << "\t-n <ep variant> \t\t variant to use for the ep parameter (1,2,3,4) (1 is standard)\n"
	   << "\t-s <Flag sidechains> \t\t 0 = no sidechain filter, 1 = sidechain filter (0 is standard)\n"
	   << "\t[-l forbidden-file] \t\t file which contains the forbidden areas of the protein\n"
	   << "\t[--noFullModel] \t\t Do not write full model, write only the loop\n";
      return 1;
    };
  
  string outputDir, inputFile, forbiddenFile;
  getArg("i", inputFile, nArgs, argv, "!");
  getArg("l", forbiddenFile, nArgs, argv, "!");
  cout << "input file: " << inputFile << endl;

  int epVariant;
  getArg( "n", epVariant, nArgs, argv, 1);

  double vdwWeight, energyWeight, endrmsWeight, compactnessWeight, propensityWeight, sidechainWeight;
  getArg( "v", vdwWeight, nArgs, argv, 0);
  getArg( "e", energyWeight, nArgs, argv, 1);
  getArg( "r", endrmsWeight, nArgs, argv, 80);
  getArg( "c", compactnessWeight, nArgs, argv, 0);
  getArg( "p", propensityWeight, nArgs, argv, 0);
  getArg( "s", sidechainWeight, nArgs, argv, 0);
  getArg( "a", index1, nArgs, argv, 0);
  getArg( "b", index2, nArgs, argv, 0);

  bool noFullModel = getArg( "-noFullModel", nArgs, argv);
  
  if ((inputFile == "!") || (index1 == 0) || (index2 == 0))     {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  if (forbiddenFile != "!") { 
    // now we have to set v1 and v2
    ifstream inForbiddenFile(forbiddenFile.c_str());
    if (!inForbiddenFile)
      ERROR("File not found.", exception);

    int anzEintraege;
    inForbiddenFile >> anzEintraege;

    for (int i = 0; i < anzEintraege; i += 1) {
      int forbPosStart, forbPosEnd;
      inForbiddenFile >> forbPosStart;
      inForbiddenFile >> forbPosEnd;
      v1.push_back(forbPosStart);
      v2.push_back(forbPosEnd);
    }
  }

  index1--;  // internal AA count starts at 0, not 1...
  index2--;

  cout << "We load the props" << endl;
  lm.load_propensities();

  // open the current file in order to process it
  ifstream inFile(inputFile.c_str());
  cout << " oeffne file mit dem Namen: " << inputFile.c_str() << endl;
  if (!inFile){
    cout <<"Error when opening the file: " << inputFile << endl;
  }

  Spacer sp;
  PdbLoader pl(inFile);
  pl.setNoHAtoms();
  sp.load(pl);
  
  cout << endl;
  LoopExtractor le;
  le.setSpacer(&sp);

    
  double looplaenge = index2 - index1 - 2;              // this stuff is soleley used for the propensities
  if (looplaenge <= 0)                                  // just that we don't divide by something stupid later on
    looplaenge = 1;                                     
  
  //set num and num2 to some reasonable values
  num = LoopModel::MAX_ITER_SOL;
  num2 = 1;

  EnergyCalculatorImpl eci;
  
  // now we set the proper wEP value
  if ( epVariant == 1) {
    cout << "we have variant 1" << endl;
    if (index2 - index1 == 3) 
      LoopTableEntry::LAMBDA_EP = 0.35;
    if (index2 - index1 == 4)
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 5) 
      LoopTableEntry::LAMBDA_EP = 1.0;
    if (index2 - index1 == 6)
      LoopTableEntry::LAMBDA_EP = 1.2;
    if (index2 - index1 == 7) 
      LoopTableEntry::LAMBDA_EP = 2.0;
    if (index2 - index1 == 8)
      LoopTableEntry::LAMBDA_EP = 3.0;
    if (index2 - index1 == 9) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 10)
      LoopTableEntry::LAMBDA_EP = 8.0;
    if (index2 - index1 == 11) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 >= 12)
      LoopTableEntry::LAMBDA_EP = 16.0;
  }
  if (epVariant == 2) {
    cout << "we have variant 2" << endl;
    if (index2 - index1 == 3) 
      LoopTableEntry::LAMBDA_EP = 0.35;
    if (index2 - index1 == 4)
      LoopTableEntry::LAMBDA_EP = 0.7;
    if (index2 - index1 == 5) 
      LoopTableEntry::LAMBDA_EP = 1.0;
    if (index2 - index1 == 6)
      LoopTableEntry::LAMBDA_EP = 1.2;
    if (index2 - index1 == 7) 
      LoopTableEntry::LAMBDA_EP = 2.0;
    if (index2 - index1 == 8)
      LoopTableEntry::LAMBDA_EP = 3.0;
    if (index2 - index1 == 9) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 10)
      LoopTableEntry::LAMBDA_EP = 8.0;
    if (index2 - index1 == 11) 
      LoopTableEntry::LAMBDA_EP = 8.0;
    if (index2 - index1 >= 12)
      LoopTableEntry::LAMBDA_EP = 16.0;
  }
  if (epVariant == 3) {
    cout << "we have variant 3" << endl;
    if (index2 - index1 == 3) 
      LoopTableEntry::LAMBDA_EP = 0.5;
    if (index2 - index1 == 4)
      LoopTableEntry::LAMBDA_EP = 5.0;
    if (index2 - index1 == 5) 
      LoopTableEntry::LAMBDA_EP = 1.0;
    if (index2 - index1 == 6)
      LoopTableEntry::LAMBDA_EP = 1.6;
    if (index2 - index1 == 7) 
      LoopTableEntry::LAMBDA_EP = 3.0;
    if (index2 - index1 == 8)
      LoopTableEntry::LAMBDA_EP = 5.0;
    if (index2 - index1 == 9) 
      LoopTableEntry::LAMBDA_EP = 3.0;
    if (index2 - index1 == 10)
      LoopTableEntry::LAMBDA_EP = 6.0;
    if (index2 - index1 == 11) 
      LoopTableEntry::LAMBDA_EP = 5.0;
    if (index2 - index1 >= 12)
      LoopTableEntry::LAMBDA_EP = 16.0;
  }
  if (epVariant == 4) {
    cout << "we have variant 4" << endl;
    if (index2 - index1 == 3) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 4)
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 5) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 6)
      LoopTableEntry::LAMBDA_EP =  4.0;
    if (index2 - index1 == 7) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 8)
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 9) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 10)
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 == 11) 
      LoopTableEntry::LAMBDA_EP = 4.0;
    if (index2 - index1 >= 12)
      LoopTableEntry::LAMBDA_EP = 4.0;
  }

  // here we calculate the loop
      vector<string> typeVec;
      
      for (unsigned int i = index1+1; i < index2+1; i++)
	typeVec.push_back(sp.getAmino(i).getType());

  vector<Spacer> vsp;
  vsp = lm.createLoopModel(sp.getAmino(index1), 
			   sp.getAmino(index1+1)[N].getCoords(), sp.getAmino(index2), 
			   sp.getAmino(index2+1)[N].getCoords(), index1, index2, num, num2, typeVec);
  
  
  // now we get the consistency values
  vector<int> consistency = lm.consistencyValues(sp, index1, index2, vsp);
  
  vector<Spacer> tempVsp;
  // remove all entries in vsp which have bad consistency values...
  for (unsigned int i = 0; i < vsp.size(); i += 1) {
    if (consistency[i] < criticalConsistency) {
	tempVsp.push_back(vsp[i]);
    }
  }
  // now correct the vsp vector:
  vsp.clear();
  consistency.clear();
  for (unsigned int i = 0; i < tempVsp.size(); i += 1) {
    vsp.push_back(tempVsp[i]);
  }
  
  // check whether we still have some solutions left...
  if (vsp.size() == 0) {
    cout << "sorry, we don't have any solutions left" << endl;
    exit(0);
  }
  
  
  // now we get the consistency values (again)
  consistency.clear();
  consistency = lm.consistencyValues(sp, index1, index2, vsp);
  
  // here we do the endRMS 
  vector<double> endRMS;
  endRMS.clear();
  endRMS = lm.rankRms2(sp, index1, index2, vsp);
  for (unsigned int i = 0; i < endRMS.size(); i += 1) 
    cout << "EndRMS of Loop Nr. " << i << ": " << endRMS[i] << endl;
  
  // here we get the Energy of the loop and the spacer
  vector<double> energy;
  for (unsigned int i = 0; i < vsp.size(); i += 1) {
    energy.push_back(eci.energyBetweenSpacer(sp, vsp[i], index1, index2));
    cout << "Energy of Loop Nr. " << i << ": " << energy[i] << endl;
  }
  
  // here we get the VDW values
  vector<int> vdw;
  vdw = lm.vdwValues(sp, index1, index2, vsp);
  for (unsigned int i = 0; i < vdw.size(); i += 1)
    cout << "gesamte VDW-Kraefte fuer Lsg. Nr.: " << i << " betraegt: " << vdw[i] << endl;
  cout << endl;
  
  // here we do the compactness
  vector<double> compactness;
  for (unsigned int i = 0; i < vsp.size(); i += 1) {
    compactness.push_back(lm.compactness(vsp[i], index1, index2, sp));
  }
  for (unsigned int i = 0; i < compactness.size(); i += 1) {
    cout << "Compactness Loop Nr. " << i << ": " << compactness[i] << endl;
  }
  
  // here we calculate the propensities
  vector<double> propensities_m;
  vector<double> propensities_a;
  lm.do_propensities(vsp, propensities_m, propensities_a, looplaenge);
  
  
  // here we do the sidechain stuff
  vector<int> sidechainCount;
  
  if (sidechainWeight == 0) {
    for (unsigned int i = 0; i < LoopModel::MAX_ITER_SOL; i += 1) {
      sidechainCount.push_back(0);
    }
  } else {
    for (unsigned int i = 0; i < vsp.size(); i += 1) {
      sp2 = sp;
      lm.setStructure(sp2, vsp[i], index1, index2);
      elc.setSpacer(&sp2);
      elc.setLoop(v1, v2, index1, index2);
      if (elc.createClusters() == true) {
	sidechainCount.push_back(elc.calculateSidechains());
      } else {
	sidechainCount.push_back(-1);
      }
    }
  }

  for(unsigned int i = 0; i < sidechainCount.size(); i += 1) {
    cout << "Sidechain collisions for Loop nr. " << i << ": " << sidechainCount[i] << endl;
  }
  
  int loopNr = index2 - index1 - gs.LOOP_MIN;
  
  // here we call the ranking function (that means the solutions in the spacer become ranked
  // by using the results of the filters)
  vector<Spacer> rankedSpacers;
  rankedSpacers = lm.rankSpacer(propensities_a, compactness, vdw, energy, 
				endRMS, consistency, vdwWeight, energyWeight,
				endrmsWeight, compactnessWeight, propensityWeight, vsp, 
				sidechainCount, sidechainWeight);
  
  vector<double> rankedSpacersRms; 
  // now we calculate the rms values for rankedSpacers
  for (unsigned int i = 0; i < rankedSpacers.size(); i += 1) {
    rankedSpacersRms.push_back(lm.calculateRms2(sp, index1, index2, rankedSpacers[i], false));
  }
  
  for (unsigned int i = 0; i < rankedSpacersRms.size(); i += 1) {
    gs.FilterCutoffGenerator(rankedSpacersRms, i + 1, loopNr);
    cout << "RMS der auf den " << i << ". Platz gerankten Loesung " << rankedSpacersRms[i] << endl;
  }

  for (unsigned int i = 0; i < vsp.size(); i++)    {
      double rms;
      cout << setw(3) << i << "   ";
      rms = lm.calculateRms2(sp, index1, index2, vsp[i], true);
    }

  for (unsigned int i = 0; i < rankedSpacers.size(); i++)    {
      string outputFile = "test.pdb";
      string tmpStr = outputFile + "." + itos(i);

      ofstream outFile(tmpStr.c_str());
      if (!outFile)
	ERROR("File not found.", exception);
      PdbSaver ps(outFile);

      if (!noFullModel)
	{
	  lm.setStructure(sp, rankedSpacers[i], index1, index2);
	  sp.save(ps);
	}
      else
	rankedSpacers[i].save(ps);
    };
  
  inFile.close();
  return 0; 
};

