/**  
@Description */
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

void sShowHelp()
{
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
       << "\t[--refine] \t\t\t Refines the solutions (ie. reduced endRms)\n"
       << "\t[--optall] \t\t\t Local optimization of all solutions pre-ranking\n"
       << "\t[--withOxygen] \t\t\t Include Oxygen atoms in RMSD calculation\n"
       << "\t--end <filename> \t\t Load length-dependent endrms weights from "
       << "file\n"
       << "\t[--interpol] \t\t\t Use interpolated RAPDF energy\n"
       << "\t[--scatter <filename>] \t\t Write a scatter plot file.\n"
       << "\t[--lambdaBase <filename>] \t Set base offset for lambda EP.\n"
       << "\t[--lambdaScale <filename>]\t Set (linear) scale for lambda EP.\n"
       << "\t[--lambdaScale <filename>]\t Set exponent for lambda EP.\n"
       << endl;
  exit(1);
}


void sSetEP(unsigned int epVariant, unsigned int index1, unsigned int index2)
{
    if ( epVariant == 1) 
    {
	cout << "we have variant 1" << "\n";
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
    else if (epVariant == 2) 
    {
	cout << "we have variant 2" << "\n";
	if (index2 - index1 == 3) 
	  LoopTableEntry::LAMBDA_EP = 0.4;
	if (index2 - index1 == 4)
	  LoopTableEntry::LAMBDA_EP = 0.8;
	if (index2 - index1 == 5) 
	  LoopTableEntry::LAMBDA_EP = 0.8;
	if (index2 - index1 == 6)
	  LoopTableEntry::LAMBDA_EP = 1.2;
	if (index2 - index1 == 7) 
	  LoopTableEntry::LAMBDA_EP = 2.0;
	if (index2 - index1 == 8)
	  LoopTableEntry::LAMBDA_EP = 3.0;
	if (index2 - index1 == 9) 
	  LoopTableEntry::LAMBDA_EP = 2.0;
	if (index2 - index1 == 10)
	  LoopTableEntry::LAMBDA_EP = 4.0;
	if (index2 - index1 == 11) 
	  LoopTableEntry::LAMBDA_EP = 4.0;
	if (index2 - index1 >= 12)
	  LoopTableEntry::LAMBDA_EP = 12.0;
    }
    else if (epVariant == 3) 
    {
	cout << "we have variant 3" << "\n";
	if (index2 - index1 == 3) 
	  LoopTableEntry::LAMBDA_EP = 0.5;
	if (index2 - index1 == 4)
	  LoopTableEntry::LAMBDA_EP = 1.0;
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
	  LoopTableEntry::LAMBDA_EP = 6.0;
	if (index2 - index1 >= 12)
	  LoopTableEntry::LAMBDA_EP = 16.0;
    }
    else if (epVariant == 4) 
    {
	cout << "we have variant 4" << "\n";
	if (index2 - index1 == 3) 
	  LoopTableEntry::LAMBDA_EP = 0.35;
	if (index2 - index1 == 4)
	  LoopTableEntry::LAMBDA_EP = 0.7;
	if (index2 - index1 == 5) 
	  LoopTableEntry::LAMBDA_EP = 1.2;
	if (index2 - index1 == 6)
	  LoopTableEntry::LAMBDA_EP = 2.0;
	if (index2 - index1 == 7) 
	  LoopTableEntry::LAMBDA_EP = 4.0;
	if (index2 - index1 == 8)
	  LoopTableEntry::LAMBDA_EP = 7.0;
	if (index2 - index1 == 9) 
	  LoopTableEntry::LAMBDA_EP = 16.0;
	if (index2 - index1 == 10)
	  LoopTableEntry::LAMBDA_EP = 20.0;
	if (index2 - index1 == 11) 
	  LoopTableEntry::LAMBDA_EP = 20.0;
	if (index2 - index1 >= 12)
	  LoopTableEntry::LAMBDA_EP = 20.0;
    }
}


double sEnergyBetweenSpacer(Spacer& sp, unsigned int index1, unsigned int 
index2, bool inter)
{
  LoopModel lm;
  if ((index1 >= sp.sizeAmino()) || (index2 >= sp.sizeAmino()))
    ERROR("Index out of range.", exception);

  double en = 0.0;

if (inter)
  {
    for (unsigned int i = index1; i < index2+1; i++)
    {
      for (unsigned int j = 0; j < index1; j++)
	en += lm.rapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
      for (unsigned int j = index2+1; j < sp.sizeAmino(); j++)
	en += lm.rapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
    }

    for (unsigned int i = index1; i < index2; i++)
      for (unsigned int j = i; j < index2+1; j++)
      en += LoopModel::rapdf.calculateEnergy(sp.getAmino(i), sp.getAmino(j));
  }
 else
   {
     for (unsigned int i = index1; i < index2+1; i++)
       {
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
}

int main(int nArgs, char* argv[])
{ 
   if (getArg( "h", nArgs, argv))
     sShowHelp();

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

//    static const int criticalConsistency = 1000; // all solutions with a worse 
//                                                 // consistency are rejected
  int createTable = 0;     // determines whether we create an output table

  // here we have all the results which become optimized later on...
  const unsigned int MAX_ARRAY = 900;

  const unsigned int max_iter_sol = 500;

  double rmsResults[MAX_ARRAY][max_iter_sol];
  double endRmsResults[MAX_ARRAY][max_iter_sol];

  int numberOfSolutions[MAX_ARRAY];

  unsigned int LoopCounter = 0;

  for (unsigned int i = 0; i < MAX_ARRAY; i += 1) {
    for (unsigned int j = 0; j < max_iter_sol; j += 1) {
      rmsResults[i][j] = -1;
      endRmsResults[i][j] = -1;
    }
  }
      
   string outputDir, inputFile, endrmsFile;
   getArg( "o", outputDir, nArgs, argv, "!");
   getArg("i", inputFile, nArgs, argv, "filelist");
   cout << "input file: " << inputFile << "\n";

   if ( outputDir == "!") 
    {
      cout << "Missing file specification. Aborting. (-h for help)" << "\n";
      return -1;
    }

   string scatterFile;
  getArg( "-scatter", scatterFile, nArgs, argv, "!");

  ofstream scatFile(scatterFile.c_str());
  if (!scatFile)
      ERROR("File not found.", exception);

  string rankFile = scatterFile + ".rank";
  ofstream outRank(rankFile.c_str());
  if (!outRank)
      ERROR("File not found.", exception);
  
  LoopModel lm;

  if (scatterFile != "!")
    {
      cout << "Doing scatter plot...\n";
      ofstream scatFile(scatterFile.c_str());
      if (!scatFile)
	ERROR("Could not open file for writing. " + scatterFile, exception);
      lm.setScatterPlot(&scatFile);
    }
  
   getArg( "-end", endrmsFile, nArgs, argv, "!");
   
   cout << "---------------------------------------------\n";
   if (endrmsFile == "!")
     lm.setAllENDRMS_WEIGHT(125);
   else
     {
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
//     getArg( "r", endrmsWeight, nArgs, argv, 80);
   getArg( "c", compactnessWeight, nArgs, argv, 0);
   getArg( "p", propensityWeight, nArgs, argv, 0);
   getArg( "s", sidechainWeight, nArgs, argv, 0);
   getArg( "t", createTable, nArgs, argv, 0);

   getArg( "-sol1", num, nArgs, argv, LoopModel::MAX_ITER_SOL);
   getArg( "-sol2", num2, nArgs, argv, 1);


   double lambdaBase, lambdaScale, lambdaExp;
   getArg( "-lambdaBase", lambdaBase, nArgs, argv, 0);
   getArg( "-lambdaScale", lambdaScale, nArgs, argv, 0);
   getArg( "-lambdaExp", lambdaExp, nArgs, argv, 1);

   if ((lambdaBase > 0) || (lambdaScale > 0))
       cout << "Using lambda EP function: Base= " << lambdaBase
	    << "\t Scale= " << lambdaScale 
	    << "\t Exp= " << lambdaExp 
	    << "\n";

   bool refine = getArg( "-refine", nArgs, argv);
   if (refine)
     cout << "refine\n";

   bool optall = getArg( "-optall", nArgs, argv);
   if (optall)
     {
       cout << "optall\n";
       LoopModel::OPT_NUM = 1000;
     }

  bool withOxygen = getArg( "-withOxygen", nArgs, argv);
  if (withOxygen)
    cout << "withOxygen\n";

  bool inter = getArg( "-interpol", nArgs, argv);
  if (inter)
    {
      cout << "interpol\n";
      lm.setInterpolated();
    }

  // open the file in which all the names of the relevant pdb-files are stored
  ifstream pdb_Filelist(inputFile.c_str());
  if (!pdb_Filelist)
    ERROR("Input file not found.", exception);

  // process all files in filelist
  while(pdb_Filelist)
  {
   
    // open the current file
    filename = "";
    pdb_Filelist >> filename;
    filename = filename + ".pdb";

    ifstream inFile(filename.c_str());
    cout << " oeffne file mit dem Namen: " << filename.c_str() << "\n";
    if (!inFile)
    {
      cout <<"pdb file from filelist not found!! Name: " << filename << "\n";
      continue;
    }
    
    // open the file for the statistic output for each pdb file
    string::size_type pos = filename.rfind("pdb"); 
    if (pos == string::npos) 
    {
      cout << "Ist kein *.pdb File, gehe zum naechsten File" << "\n";
      continue;
    }

    Spacer sp;
    PdbLoader pl(inFile);
    pl.setNoHAtoms();
    sp.load(pl);
   
    cout << "\n";

    // write secondary structure signature:
    for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    {
      if (sp.getAmino(i).getState() == HELIX)
	cout << "H";
      else if (sp.getAmino(i).getState() == STRAND)
	cout << "E";
      else
	cout << "-";
      if ((i+1) % 60 == 0)
	cout << "\n";
    }
    cout << "\n";


    LoopExtractor le;
    le.setSpacer(&sp);

    while (start != -1) 
    {
      le.nextLoop (start, end);
      cout << "Start: " << start << " end: " << end << "\n";
      if (end - start > gs.LOOP_MAX - 1 || end - start < gs.LOOP_MIN - 1 ) 
	  continue;
      if (start <= 1 || (unsigned int)end == sp.sizeAmino() - 1)
	  continue;

      // check for reasonable B-factors in loop:
      double maxB = 0.0;

      for (unsigned int i = start-2; i < end+1; i++)
	for (unsigned int j = 0; j < sp.getAmino(i).sizeBackbone(); j++)
	  if (sp.getAmino(i)[j].getBFac() > maxB)
	    maxB = sp.getAmino(i)[j].getBFac();

      if (maxB > 25.0)
	{
	  cout << "Bad B-factor (" << maxB<< ") found.\n";
	  continue;
	}

      // now we have to set index2 and index1
      index2 = end;                     // here we had end+1
      index1 = start - 1;               // here we had start
      cout << "index1 (-s): " << index1 << " index2 (-e) " << index2 << "\n";

      double looplaenge = index2 - index1 - 2;              
      // this stuff is soleley used for the propensities
      if (looplaenge <= 0)                                  
	looplaenge = 1;                                     
      // just that we don't divide by something stupid later on


      // write loop signature
      for (unsigned int i = 0; i < sp.sizeAmino(); i++)
      {
	  if ((i >= index1) && (i <= index2))
	      cout << "#";
	  else
	      cout << ".";
	  if ((i+1) % 60 == 0)
	      cout << "\n";
      }
      cout << "\n";

      
      // now we set the proper wEP value
      if ((lambdaBase > 0) || (lambdaScale > 0))
	  LoopTableEntry::LAMBDA_EP = lambdaBase 
	      + lambdaScale * pow((index2-index1),lambdaExp);
      else
	  sSetEP(epVariant, index1, index2);

      // here we calculate the loop
      
      lm.releaseTables(4);

      vector<string> typeVec;
      
      for (unsigned int i = index1+1; i < index2+1; i++)
	typeVec.push_back(sp.getAmino(i).getType());

      vector<Spacer> vsp;
      vsp.reserve(LoopModel::MAX_ITER_SOL);
      vsp = lm.createLoopModel(sp.getAmino(index1), 
			       sp.getAmino(index1+1)[N].getCoords(), 
			       sp.getAmino(index2), 
			       sp.getAmino(index2+1)[N].getCoords(), 
			       index1, index2, num, num2, typeVec);

      if (refine)
	  lm.refineModel(sp, index1, index2, vsp);
	
      if (optall)      // attempt local minimization of all solutions
	  lm.optimizeModel(sp, index1, index2, vsp, false);
	
      if (vsp.size() == 0) 
	  continue;

      int loopNr = index2 - index1 - gs.LOOP_MIN;

      vector<Spacer> rankedSpacers;
      rankedSpacers.reserve(LoopModel::MAX_ITER_SOL);

      for (unsigned int i = 0; i < vsp.size(); i++)
	  rankedSpacers.push_back(vsp[i]); 
      lm.rankRawScore(sp, index1, index2, rankedSpacers);

      vector<double> rankedSpacersRms; 

      // now we calculate the rms values for rankedSpacers
      for (unsigned int i = 0; i < rankedSpacers.size(); i += 1) 
	  rankedSpacersRms.push_back(lm.calculateRms2(sp, index1, index2, 
						      rankedSpacers[i], false,
						      withOxygen));

       // here we sort the solutions by rms
      for (unsigned int i = 0; i < vsp.size(); i++)
	{
	  double rms = lm.calculateRms2(sp, index1, index2, vsp[i], true, 
					withOxygen);
	  // now fill the array for the search
	  rmsResults[LoopCounter][i] = rms;
	  ranking_helper2 rh = ranking_helper2(i, rms);
	  rmsh.insert(rh);
	}

     // now we have to check, whether it is a "bug" value
      set<ranking_helper2>::iterator pos;
      pos = rmsh.begin();

      if ( pos->get_value() > (index2 - index1) * 2 ) 
      {
	prop_ranking.clear();
	vdw_ranking.clear();
	rmsh.clear();
	endRMS_ranking.clear();
	cmpct_ranking.clear();
	energy_ranking.clear();
	continue;
      }
      
      for (unsigned int i = 0; i < rankedSpacersRms.size(); i += 1) 
      {
	  if (i < 20) 
	      gs.FilterCutoffGenerator(rankedSpacersRms, i + 1, loopNr);
	  cout << "RMS of # " << i << " ranked solution = " 
	       << rankedSpacersRms[i] << "\n";
      }

      lm.doScatterPlot(sp, index1, index2, rankedSpacers, withOxygen);
      // calculate "virtual" rank of native loop:
      double nat = lm.calcOrigEnergy(sp, index1, index2);
      
      outRank << "L= " << index2-index1 << "\t";
      
      for (unsigned int i = 0; i < rankedSpacers.size(); i++)
      {
	  double tmpEndRms = lm.calculateRms(sp, index1, index2, 
					     rankedSpacers[i], false);
	  double tmpEnergy = lm.calculateEnergy(sp, index1, 
						index2, rankedSpacers[i]);
	  
	  if (nat <= (lm.getENDRMS_WEIGHT() * tmpEndRms 
		      + energyWeight * tmpEnergy))
	  {
	      outRank << i;
	      break;
	  }
      }
      outRank << endl;
    
      // update the statistic arrays
      gs.updateRMSArray(loopNr, rmsh);
      gs.updateRankedRmsArray(loopNr, rankedSpacersRms);


      // fill the arrays for the search
//        for (unsigned int i = 0; i < vsp.size(); i += 1) 
//        {
//  	  endRmsResults[LoopCounter][i] = endRMS[i];
//  	  energyResults[LoopCounter][i]= energy[i];
//  	  propensityResults[LoopCounter][i] = propensities_a[i];
//  	  compactnessResults[LoopCounter][i] = compactness[i];
//  	  loopLength[LoopCounter] = index2 - index1;
//        }
      numberOfSolutions[LoopCounter] = vsp.size();
      LoopCounter += 1;

      // clean up the arrays
      prop_ranking.clear();
      vdw_ranking.clear();
      rmsh.clear();
      endRMS_ranking.clear();
      cmpct_ranking.clear();
      energy_ranking.clear();
    
      vsp.clear();
      rankedSpacers.clear();
    }
    start = 0;
    end = 0;
    inFile.close();

    pdb_count += 1;                                      
    if (pdb_count % 1 == 0) 
    {
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
      cout << "Name of Statistic Output Files: " << outN << "\n";
      ofstream outF(outN.c_str());
      if (!outF)
	ERROR("File for the statistic output couldn't be created.", exception);
      gs.outputArray(outF);
    }

  }

  exit(0);
}
