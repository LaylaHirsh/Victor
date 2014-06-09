/**
 * @Description This program runs through filelist and generates the output  for the optimization
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

void sShowHelp(){
  cout << "LoopIteration\n"
       << "This program runs through filelist and generates the output " 
       << "for the optimization \n"
       << " Options: \n"
       << "\t-o <Directory> \t\t\t Output Directory (without /)\n"
       << "\t-i <Input File \t\t\t List of PDB-files (default = " 
       << "filelist)\n"
       << "\t-v <Weight VDW> \t\t Weight of the vdw filter " 
       << "(default = 0)\n"
       << "\t-e <Weight Energy> \t\t Weight of the energy filter " 
       << "(default = 1)\n"
       << "\t-r <Weight EndRMS> \t\t Weight of the end-rms filter "
       << "--- DISCONTINUED!\n"
       << "\t-c <Weight compactness> \t Weight of the compactness " 
       << "filter (default = 0)\n"
       << "\t-p <Weight propensity> \t\t Weight of the propensity " 
       << "filter (default = 0) \n"
       << "\t-n <ep variant> \t\t variant to use for the ep parameter " 
       << "(1,2,3,4) (default = 1)\n"
       << "\t-s <Flag sidechains> \t\t Sidechain filter: 0 = off, 1 = " 
       << "of (default = 0)\n"
       << "\t --sol1 <number> \t\t Number of primary solutions "
       << "(default = " << LoopModel::MAX_ITER_SOL << ")\n"
       << "\t --sol2 <number> \t\t Number of secondary solutions "
       << "(default = 1)\n"
       << "\t[--refine] \t\t Refines the solutions (ie. reduced endRms)\n"
       << "\t[--optall] \t\t Local optimization of all solutions pre-ranking\n"
       << "\t[--withOxygen] \t\t Include Oxygen atoms in RMSD calculation\n"
       << "\t--end <filename> \t\t Load length-dependent endrms weights from "
       << "file\n"
       << "\t[--interpol] \t\t Use interpolated RAPDF energy\n"
       << "\t[--scatter <filename>] \t Write a scatter plot file.\n"
       << endl;
  exit(1);
};

double sEnergyBetweenSpacer(Spacer& sp, unsigned int index1, unsigned int 
index2, bool inter){
  LoopModel lm;
  if ((index1 >= sp.sizeAmino()) || (index2 >= sp.sizeAmino()))
    ERROR("Index out of range.", exception);

  double en = 0.0;

if (inter)  {
    for (unsigned int i = index1; i < index2+1; i++)    {
      for (unsigned int j = 0; j < index1; j++)
	en += lm.irapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
      for (unsigned int j = index2+1; j < sp.sizeAmino(); j++)
	en += lm.irapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
    }

    for (unsigned int i = index1; i < index2; i++)
      for (unsigned int j = i; j < index2+1; j++)
      en += LoopModel::irapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
  }
 else   {
    for (unsigned int i = index1; i < index2+1; i++)    {
      for (unsigned int j = 0; j < index1; j++)
	en += lm.rapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
      for (unsigned int j = index2+1; j < sp.sizeAmino(); j++)
	en += lm.rapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
    }

    for (unsigned int i = index1; i < index2; i++)
      for (unsigned int j = i; j < index2+1; j++)
      en += LoopModel::rapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
   } 
  return en;
};

int main(int nArgs, char* argv[]){ 
  cout << "Start\n";
  
  multiset<ranking_helper> prop_ranking;   // used to sort the propensities
  multiset<ranking_helper2> rmsh;          // used to sort the rms values
  multiset<ranking_helper2> vdw_ranking;   // used to sort the vdw terms
  multiset<ranking_helper2> endRMS_ranking;// used to sort the endRMS values
  multiset<ranking_helper2> cmpct_ranking; // used to sort the compactness 
                                           // values
  multiset<ranking_helper2> energy_ranking;// used to sort the energy values
  
  unsigned int index1 = 0; // contains the start index of the loop
  unsigned int index2 = 0; // contains the end index of the loop
  unsigned int num;        // determines how many solutions we are getting 
                           // in the first level of the loop modeling
  unsigned int num2;       // determines how many solutions we are getting 
                           // in the second level of the loop modeling

  int start = 0;           // contains the beginning of another loop in 
                           // the pdb-file
  int end = 0;             // contains the end of another loop in the pdb-file
  string filename;         // contains the name of the current pdb file
  int pdb_count = 0;       // counts the number of examined pdb-files  
  globalStatistic gs;      // contains the overall statistic for the whole 
                           // thing (e.g. propensities) 

  static const int criticalConsistency = 1000; // all solutions with a worse 
                                               // consistency are rejected
  int createTable = 0;     // determines whether we create an output table

  // here we have all the results which become optimized later on...
  const unsigned int MAX_ARRAY = 1800;

  const unsigned int max_iter_sol = 200;

  unsigned int loopLength[MAX_ARRAY];

  double rmsResults[MAX_ARRAY][max_iter_sol];
  double endRmsResults[MAX_ARRAY][max_iter_sol];
  double energyResults[MAX_ARRAY][max_iter_sol];
  double propensityResults[MAX_ARRAY][max_iter_sol];
  double compactnessResults[MAX_ARRAY][max_iter_sol];

  int numberOfSolutions[MAX_ARRAY];

  int dummyArray[max_iter_sol];
  for (int i = 0; i < max_iter_sol; i += 1) {
    dummyArray[i] = 0;
  }
  unsigned int LoopCounter = 0;

  for (unsigned int i = 0; i < MAX_ARRAY; i += 1) {
    for (unsigned int j = 0; j < max_iter_sol; j += 1) {
      rmsResults[i][j] = -1;
      endRmsResults[i][j] = -1;
      energyResults[i][j] = -1;
      propensityResults[i][j] = -1;
      compactnessResults[i][j] = -1;
    }
  }
      
  vector<int> v1, v2;
  Spacer sp2;
  evaluateLoopSidechainCapability elc(&sp2, v1, v2, index1, index2);

   if (getArg( "h", nArgs, argv))
     sShowHelp();

   string outputDir, inputFile, endrmsFile;
   getArg( "o", outputDir, nArgs, argv, "!");
   getArg("i", inputFile, nArgs, argv, "filelist");
   cout << "input file: " << inputFile << "\n";

   if ( outputDir == "!")     {
      cout << "Missing file specification. Aborting. (-h for help)" << "\n";
      return -1;
    }

   string scatterFile;
  getArg( "-scatter", scatterFile, nArgs, argv, "!");

  ofstream scatFile(scatterFile.c_str());
  if (!scatFile)
      ERROR("File not found.", exception);
  
  LoopModel lm(scatFile);
  if (scatterFile != "!")
      lm.setScatterPlot();
  
   getArg( "-end", endrmsFile, nArgs, argv, "!");
   
   cout << "---------------------------------------------\n";
   if (endrmsFile == "!")
     lm.setAllENDRMS_WEIGHT(125);
   else     {
       cout << "Loading endrms weights from file:  " << endrmsFile << "\n";
       ifstream inF3(endrmsFile.c_str());
       if (!inF3)
	 ERROR("The endrms file could not be opened.", exception);
       lm.loadENDRMS_WEIGHT(inF3);
     }

   cout << "Starting ENDRMS_WEIGHT:\n";
   lm.saveENDRMS_WEIGHT(cout);
   cout << "---------------------------------------------\n";
   
   int epVariant;
   getArg( "n", epVariant, nArgs, argv, 1);

   double vdwWeight, energyWeight, compactnessWeight, // endrmsWeight, 
     propensityWeight, sidechainWeight;
   getArg( "v", vdwWeight, nArgs, argv, 0);
   getArg( "e", energyWeight, nArgs, argv, 1);
   getArg( "c", compactnessWeight, nArgs, argv, 0);
   getArg( "p", propensityWeight, nArgs, argv, 0);
   getArg( "s", sidechainWeight, nArgs, argv, 0);
   getArg( "t", createTable, nArgs, argv, 0);

   getArg( "-sol1", num, nArgs, argv, LoopModel::MAX_ITER_SOL);
   getArg( "-sol2", num2, nArgs, argv, 1);

   bool refine = getArg( "-refine", nArgs, argv);
   if (refine)
     cout << "refine\n";

   bool optall = getArg( "-optall", nArgs, argv);
   if (optall)     {
       cout << "optall\n";
       LoopModel::OPT_NUM = 1000;
     }

  bool withOxygen = getArg( "-withOxygen", nArgs, argv);
  if (withOxygen)
    cout << "withOxygen\n";

  bool inter = getArg( "-interpol", nArgs, argv);
  if (inter)    {
      cout << "interpol\n";
      lm.setInterpolated();
    }

  // open the file in which all the names of the relevant pdb-files are stored
  ifstream pdb_Filelist(inputFile.c_str());
  if (!pdb_Filelist)
    ERROR("Input file not found.", exception);

  cout << "We load the props" << "\n";
  lm.load_propensities();

  // as long as we have files to process in the filelist
  while(pdb_Filelist){
   
    // open the current file in order to process it
    filename = "";
    pdb_Filelist >> filename;
    ifstream inFile(filename.c_str());
    cout << " oeffne file mit dem Namen: " << filename.c_str() << "\n";
    if (!inFile){
      cout <<"pdb file from filelist not found!! Name: " << filename << "\n";
      continue;
    }
    
    // open the file for the statistic output for each pdb file
    string::size_type pos = filename.rfind("pdb"); 
    if (pos == string::npos) {
      cout << "Ist kein *.pdb File, gehe zum naechsten File" << "\n";
      continue;
    }

    Spacer sp;
    PdbLoader pl(inFile);
    pl.setNoHAtoms();
    sp.load(pl);
   
    cout << "\n";
    LoopExtractor le;
    le.setSpacer(&sp);

    while (start != -1) {
      le.nextLoop (start, end);
      cout << "Start: " << start << " end: " << end << "\n";
      if (end - start > gs.LOOP_MAX - 1 || end - start < gs.LOOP_MIN - 1 ) 	{  //gs.LOOP_MAX - 1  gs.LOOP_MIN - 1
	  continue;
	}
      if (start == 0 || (unsigned int)end == sp.sizeAmino() - 1)
	continue;

      
      // now we have to set index2 and index1
      index2 = end;                     // here we had end+1
      index1 = start - 1;               // here we had start
      cout << "index1 (-s): " << index1 << " index2 (-e) " << index2 << "\n";

      double looplaenge = index2 - index1 - 2;              
      // this stuff is soleley used for the propensities
      if (looplaenge <= 0)                                  
	looplaenge = 1;                                     
      // just that we don't divide by something stupid later on

      // output the energy between the original loop and the framework

	double originalEnergy;
	originalEnergy = sEnergyBetweenSpacer(sp, index1, index2, inter);

	cout << "Energy des Original-Loops: " << originalEnergy << "\n";
      
      // now we set the proper wEP value
//        if ( epVariant == 1) {
	cout << "we have variant 1" << "\n";
	if (index2 - index1 == 3) 
	  LoopTableEntry::LAMBDA_EP = 0.25;
	if (index2 - index1 == 4)
	  LoopTableEntry::LAMBDA_EP = 0.25;
	if (index2 - index1 == 5) 
	  LoopTableEntry::LAMBDA_EP = 0.5;
	if (index2 - index1 == 6)
	  LoopTableEntry::LAMBDA_EP = 0.5;
	if (index2 - index1 == 7) 
	  LoopTableEntry::LAMBDA_EP = 0.5;
	if (index2 - index1 == 8)
	  LoopTableEntry::LAMBDA_EP = 0.5;
	if (index2 - index1 == 9) 
	  LoopTableEntry::LAMBDA_EP = 0.75;
	if (index2 - index1 == 10)
	  LoopTableEntry::LAMBDA_EP = 0.75;
	if (index2 - index1 == 11) 
	  LoopTableEntry::LAMBDA_EP = 1.0;
	if (index2 - index1 >= 12)
	  LoopTableEntry::LAMBDA_EP = 1.5;
      // here we calculate the loop
      
      lm.releaseTables(4);

      vector<string> typeVec;
      
      for (unsigned int i = index1+1; i < index2+1; i++)
	typeVec.push_back(sp.getAmino(i).getType());

      vector<Spacer> vsp;
      vsp = lm.createLoopModel(sp.getAmino(index1), 
			       sp.getAmino(index1+1)[N].getCoords(), 
			       sp.getAmino(index2), 
			       sp.getAmino(index2+1)[N].getCoords(), 
			       index1, index2, num, num2, typeVec);
      if (refine)	{
	  lm.refineModel(sp, index1, index2, vsp);
	}
	
      if (optall)    { // attempt local minimization of all solutions
		  lm.optimizeModel(sp, index1, index2, vsp, false);
	}
	
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
	continue;
      }

      // here we do the endRMS
      vector<double> endRMS = lm.rankRms2(sp, index1, index2, vsp);

      // now we get the consistency values (again)
      consistency.clear();
      consistency = lm.consistencyValues(sp, index1, index2, vsp);

       // here we do the endRMS (again)
      endRMS.clear();
      endRMS = lm.rankRms2(sp, index1, index2, vsp);
      for (unsigned int i = 0; i < endRMS.size(); i += 1) 
	cout << "EndRMS of Loop Nr. " << i << ": " << endRMS[i] << "\n";

      // here we get the Energy of the loop and the spacer
      vector<double> energy;
	for (unsigned int i = 0; i < vsp.size(); i += 1) {
	  energy.push_back( (lm.calculateEnergy(sp, index1, 
						   index2, vsp[i])) );
	  cout << "Energy of Loop Nr. " << i << ": " << energy[i] << "\n";
      }

      // here we get the VDW values
      vector<int> vdw;
      vdw = lm.vdwValues(sp, index1, index2, vsp);
      for (unsigned int i = 0; i < vdw.size(); i += 1)
	cout << "gesamte VDW-Kraefte fuer Lsg. Nr.: " << i << " betraegt: " 
	     << vdw[i] << "\n";
      cout << "\n";

      // here we do the compactness
      vector<double> compactness;
      for (unsigned int i = 0; i < vsp.size(); i += 1) {
	compactness.push_back(lm.calculateSolvation(vsp[i], index1, index2, sp));
      }
      for (unsigned int i = 0; i < compactness.size(); i += 1) {
	cout << "Compactness Loop Nr. " << i << ": " << compactness[i] << "\n";
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
	cout << "Sidechain collisions for Loop nr. " << i << ": " 
	     << sidechainCount[i] << "\n";
      }

      int loopNr = index2 - index1 - gs.LOOP_MIN;

      // here we call the ranking function (that means the solutions 
      // in the spacer become ranked by using the results of the filters)
      vector<Spacer> rankedSpacers;
      rankedSpacers = lm.rankSpacer(propensities_m, compactness, vdw, energy, 
			   endRMS, consistency, vdwWeight, energyWeight, 
			   lm.getENDRMS_WEIGHT(index2-index1), // endrmsWeight, 
			   compactnessWeight, propensityWeight, vsp, 
			   sidechainCount, sidechainWeight);

      vector<double> rankedSpacersRms; 
      // now we calculate the rms values for rankedSpacers
      for (unsigned int i = 0; i < rankedSpacers.size(); i += 1) {
	rankedSpacersRms.push_back(lm.calculateRms2(sp, index1, index2, 
						    rankedSpacers[i], false, withOxygen));
      }

      for (unsigned int i = 0; i < rankedSpacersRms.size(); i += 1) {
	gs.FilterCutoffGenerator(rankedSpacersRms, i + 1, loopNr);
	cout << "RMS der auf den " << i << ". Platz gerankten Loesung " 
	     << rankedSpacersRms[i] << "\n";
      }

      // here we sort the solutions by rms
      cout << "non-filtered RMS (ie. original ranking):\n";
      for (unsigned int i = 0; i < vsp.size(); i++)	{
	  double rms;
  	  cout << setw(3) << i << "   ";
	  rms = lm.calculateRms2(sp, index1, index2, vsp[i], true, withOxygen);
	  // now fill the array for the search
	  rmsResults[LoopCounter][i] = rms;
	  ranking_helper2 rh = ranking_helper2(i, rms);
	  rmsh.insert(rh);
	}

      // now we have to check, whether it is a "bug" value
      set<ranking_helper2>::iterator pos;
      pos = rmsh.begin();
      if ( (index2 - index1) * 2 < pos->get_value()) {
	prop_ranking.clear();
	vdw_ranking.clear();
	rmsh.clear();
	endRMS_ranking.clear();
	cmpct_ranking.clear();
	energy_ranking.clear();
	continue;
      }
      
      // update the statistic arrays
     
      ASSERT(loopNr >= 0 && LoopNr < gs.LOOP_MAX - gs.LOOP_MIN + 1, exception);
      
      gs.updateRMSArray(loopNr, rmsh);
      gs.updateRankedRmsArray(loopNr, rankedSpacersRms);


      // fill the arrays for the search
      for (unsigned int i = 0; i < vsp.size(); i += 1) {
	endRmsResults[LoopCounter][i] = endRMS[i];
	energyResults[LoopCounter][i]= energy[i];
	propensityResults[LoopCounter][i] = propensities_a[i];
	compactnessResults[LoopCounter][i] = compactness[i];
	loopLength[LoopCounter] = index2 - index1;
      }
      numberOfSolutions[LoopCounter] = vsp.size();
      LoopCounter += 1;

      // clean up the arrays
      prop_ranking.clear();
      vdw_ranking.clear();
      rmsh.clear();
      endRMS_ranking.clear();
      cmpct_ranking.clear();
      energy_ranking.clear();
    }
    start = 0;
    end = 0;
    inFile.close();

    pdb_count += 1;                                      
    if (pdb_count % 1 == 0) {
      string outN;
      if (pdb_count < 100)
	{
	  if (pdb_count < 10)
	    outN = outputDir + "/count00" + itos(pdb_count);
	  else
	    outN = outputDir + "/count0" + itos(pdb_count);
	}
      else
	outN = outputDir + "/count" + itos(pdb_count);
      cout << "Name des Statistic Output Files welches erzeugt wird: " 
	   << outN << "\n";
      ofstream outF(outN.c_str());
      if (!outF)
	ERROR("File for the statistic output couldn't be created.", exception);
      gs.outputArray(outF);
    }
  }

  // here we output all the values of the filter and the loops into the 
  // file ParametersToOptimize
  if (createTable != 0) {
    ofstream outF2("ParametersToOptimize");
    if (!outF2)
      ERROR("The ParametersToOptimize file couldn't be created.", exception);
    for (unsigned int i = 0; i < LoopCounter; i += 1) {
      for (int j = 0; j < LoopModel::MAX_ITER_SOL; j += 1) {
	if (rmsResults[i][j] == -1) {
	  break;
	}
	outF2 << rmsResults[i][j] << " "; 
	outF2 << endRmsResults[i][j] << " "; 
	outF2 << energyResults[i][j] << " "; 
	outF2 << propensityResults[i][j] << " "; 
	outF2 << compactnessResults[i][j] << "\n";
      }
      outF2 << "\n";
    } 
  }

  // here we have a complete search !!
  double accumulatedRms = 0;
  double bestEndRms = 0; // bestPropensity = 0, bestCompactness = 0;
  int selectedNumber;
  cout << "we have " << LoopCounter << " Loops" << "\n";

  for (unsigned int i = 0; i < LoopCounter; i += 1)     {
      selectedNumber = lm.rankSpacerForOptimization(propensityResults[i], 
		        compactnessResults[i], dummyArray, energyResults[i], 
			endRmsResults[i], dummyArray, 0, 1, 
                        lm.getENDRMS_WEIGHT() , 0, 0, dummyArray, 0);
      accumulatedRms += rmsResults[i][selectedNumber];
  }

  cout << "for the parameters e = 1, r = " << lm.getENDRMS_WEIGHT() 
       << " we have " << accumulatedRms << " as the accumulatedRms " << "\n";

 // ---> main optimization loop follows:
  for (unsigned lc = gs.LOOP_MIN; lc <= gs.LOOP_MAX; lc++)    {
      cout << "---------------------------------------------\n";
      cout << ">>>>> lc = " << lc << "\n";

      accumulatedRms = 0;
      for (unsigned int i = 0; i < LoopCounter; i += 1) 
	if (loopLength[i] == lc)
	  {
	    selectedNumber = lm.rankSpacerForOptimumEndrms(
				 compactnessResults[i], energyResults[i], 
				 endRmsResults[i], lm.getENDRMS_WEIGHT(lc));
	    accumulatedRms += rmsResults[i][selectedNumber];
	  }

      double bestAccumulation = 0.995 * accumulatedRms;
      cout << "Default optimum:  EndRms: " 
	   << lm.getENDRMS_WEIGHT(lc) << " accumulated value of: " 
	   << accumulatedRms << "\n";
      
      for (double eWeight = 0; eWeight <= 500; eWeight += 0.25) 	{
	  accumulatedRms = 0;
	  for (unsigned int i = 0; i < LoopCounter; i += 1) 
	    if (loopLength[i] == lc)
	      {
		selectedNumber = lm.rankSpacerForOptimumEndrms(
				 compactnessResults[i], energyResults[i], 
				 endRmsResults[i], eWeight);
		accumulatedRms += rmsResults[i][selectedNumber];
	      }
	  
	  if (accumulatedRms < bestAccumulation) 	    {
	      bestAccumulation = 0.995 * accumulatedRms;
	      bestEndRms = eWeight;
	      cout << "\n New optimum:  EndRms: " 
		   << eWeight << " accumulated value of: " 
		   << accumulatedRms << "\n";

	      lm.setENDRMS_WEIGHT(eWeight, lc);
	    }

	  cout << "e = " << eWeight << "\t";
	}
      
    }
 // <--- end of main optimization loop 

    cout << "---------------------------------------------\n";
    cout << "Optimum ENDRMS_WEIGHT:\n";
    lm.saveENDRMS_WEIGHT(cout);
    cout << "---------------------------------------------\n";
      
    ofstream outF3("best_endrms_weights");
    if (!outF3)
      ERROR("The best_endrms_weights file could not be created.", exception);
    lm.saveENDRMS_WEIGHT(outF3);
      
  return 0; 
};


