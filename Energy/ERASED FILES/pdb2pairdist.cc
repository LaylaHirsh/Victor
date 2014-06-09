/**
 @Description This program give information Distance-dependent pairwise potential  
       * of a given protein structural model.
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <AminoAcid.h>
#include <String2Number.h>
#include <vector>
#include <cmath>

using namespace Biopool;

void sShowHelp(){
  cout << "pdb2pairdist - Pdb Distance-dependent pairwise potential\n\n"
       << "This program give information Distance-dependent pairwise potential" 
       << " of a given protein structural model. \n"
       << "Contact the author for further details.\n"
       << "Developed by AJMM <albertoj@bio.unipd.it> \n\n"
       << "Input options: \n"
       << "\t-i <filename> \t\t Input single pdb file.\n"
       << "\t[--chain <id>]  \t ID of chain to load from PDB file.\n"
       << "\t[--nmr] \t\t Calculate average over NMR ensemble ( only for NMR model ).\n"
       << "\t\t\t\t DEFAULT: The first chain and/or the first model.\n\n"
       << "\t-I <filelist> \t\t Input filelist (a file containing the list of PDB file)." <<"\n"
       << "\t\t\t\t For major information upon the format of the input filelist and the option\n"
       << "\t\t\t\t type: pdb2pairdist --filelist \n"
       << "\t --maxd X\t\t\t max dist X (default 25A int).\n"
       << "\t --mind X\t\t\t min dist X (default 3A int).\n"
       << "\t --bin X\t\t\t bin size X (default 0.5A double, reminder of (maxd-mind)/bin should be 0).\n"
       << "\t --seqsep X\t\t\t min seq separation to consider contacts X (default >= 7 positions int).\n"
       << "\t --ca\t\t\t C alpha as interaction centers.\n"
       << "\t --cb\t\t\t C beta as interaction centers (Ca for GLYs).\n"
       << "\t --bb\t\t\t C beta as interaction centers (pseudo Cb for GLYs as in HSE).\n"
       << "\t --cae\t\t\t C alpha as interaction centers in hemispheres as in HSE.\n"
       << "Output option: \n"
       << "\t -r \t\t\t Output in format for data file.\n"
       << "\n";
}
void sSpecification() {
  cout <<"\n\nINPUT FILELIST:\n"
       <<"The file list must contain the list of pdb file to examinate or their path.\n"
       <<"The file must be organized in column: the first containing the pdb file path\n"
       <<"and the second the Chain ( optional ).\n"
       <<"If you want to consider the chain specified in the second column you have to\n"
       <<"use the option --complete. If not specified the second column will be ignored\n"
       <<"and the programm will consider the first chain.\n"
       <<"If you use a filelist of nmr structure use the option --nmr for consider all models.\n\n"
       <<"FILELIST FORMAT EXAMPLE:\n\n"
       <<"PATH_OF_PDB_STRUCTURE1\tCHAIN\n"
       <<"PATH_OF_PDB_STRUCTURE2\tCHAIN\n"
       <<"....\n\n";
}

int getbin (double dist, double mind, double binnums){
	return (int)((dist-mind)/((double)binnums));

}

int getbins( double maxd, double mind, double bins){	
	if(0.0 != bins)
		return (static_cast<int>((maxd-mind)/bins));
}

int main(int nArgs, char* argv[]){
	
	if (getArg( "h", nArgs, argv)){
		sShowHelp();
		return 1;
	}
	if (getArg( "-filelist", nArgs, argv)){
		sSpecification();
		return 1;
	}
  
	string inputFile, inputFilelist, chainID, inputFilea, line, lineSS, lineSA, lineNUM;
	double maxd = 25.0;
	double mind = 3.0;
	double bin = 0.5;
	int seqsep = 7;
	getArg( "i", inputFile, nArgs, argv, "!");
	getArg( "I", inputFilelist, nArgs, argv, "!");
	getArg( "-chain", chainID, nArgs, argv, " ");
	getArg( "-maxd", maxd, nArgs, argv, 25.0);
	getArg( "-mind", mind, nArgs, argv, 3.0);
	getArg( "-bin", bin, nArgs, argv, 0.5);
	getArg( "-seqsep", seqsep, nArgs, argv, 7);
	bool nmr = getArg( "-nmr", nArgs, argv);
	bool complete = getArg( "-complete", nArgs, argv);
	bool ca = getArg( "-ca", nArgs, argv);
	bool cb = getArg( "-cb", nArgs, argv);
	bool bb = getArg( "-bb", nArgs, argv);
	bool cae = getArg( "-cae", nArgs, argv);
	bool verbose = getArg( "r", nArgs, argv);

	//~ //Control for ChaiId option
	//~ if ((chainID != " ") && (complete))
	//~ ERROR("You are using the 'chainID' option and the 'complete' option at the same time.", exception);

	//Control for input file
	if ((inputFile == "!") && (inputFilelist == "!")){
		cout << "Missing file specification. Aborting. (-h for help)" << endl;
		return -1;
	}
  
	if ((inputFile != "!") && (inputFilelist != "!")){
		cout << "Please choose between filelist and file mode. Aborting. " << "(-h for help)" << endl;
		return -2;     
	}
	if (( !ca ) && ( !cb ) && ( !bb ) && ( !cae ))
		ERROR("Choose a valid option for output.", exception);

	ifstream inFile(inputFilelist.c_str());
  
	if ((!inFile) && (inputFilelist != "!"))
		ERROR("File not found.", exception);

	if(0 != fmod((maxd-mind),bin)){
		cout << "wrong bin size, reminder of (maxd-mind)/bin should be 0" << endl;
		return -1;
	}
	
	int totalanalized = 0;
	int totalmodel = 0;
	int binnums = getbins( maxd, mind, bin);

	while ((inFile) || (inputFile != "!")){
		if (inputFilelist != "!"){
			inFile >> inputFile;
			if (!inFile)
				break;
			if (complete){
				int check = checkForBlankLine(inFile, true);	     
				if ( check == 1 )
					chainID = " ";
				else
					inFile >> chainID;
				skipToNewLine(inFile);
			}
		}
	      
		if ( inputFile == "!" ){
			cout << "Missing file specification. Aborting. (-h for help)" << endl;
			return -1;
		}
	    
		
		ifstream inFile2(inputFile.c_str());
		if (!inFile2)
			ERROR("File not found.", exception);
		PdbLoader pl(inFile2);
		pl.setNoHAtoms();
		pl.setNoVerbose();
		pl.setChain(chainID.c_str()[0]);
		pl.setPermissive();

		inputFilea = inputFile;
		inputFilea.append(".a");// its supposed to contain the adataset file for only the desired chain, and it's supossed to be in the same folder as the struct file
		ifstream inFile3(inputFilea.c_str());
		if (!inFile3)
			ERROR("File not found.", exception);
		getline (inFile3,line);//0
		getline (inFile3,line);//1
		int num = atoi(line.c_str());
		getline (inFile3,line);//2
		getline (inFile3,lineSS);//3
		getline (inFile3,line);//4
		getline (inFile3,line);//5
		getline (inFile3,lineSA);//6
		getline (inFile3,line);//7
		getline (inFile3,line);//8
		getline (inFile3,line);//9
		getline (inFile3,line);//10
		getline (inFile3,line);//11
		getline (inFile3,line);//12
		getline (inFile3,line);//13
		getline (inFile3,lineNUM);//14
		istringstream ss(lineSS);
		istringstream sa(lineSA);
		istringstream nu(lineNUM);

		int** adata;
		char ssc;
		int sanum;
		adata = new int*[num];
		for(int c = 0; c < num; c++){
			adata[c] = new int [3];
			ss >> ssc;
			if(ssc == 'H')
				adata[c][1] = 0;
			else if(ssc == 'G')
				adata[c][1] = 0;
			else if(ssc == 'I')
				adata[c][1] = 0;
			else if(ssc == 'E')
				adata[c][1] = 1;
			else if(ssc == 'B')
				adata[c][1] = 1;
			else
				adata[c][1] = 2;
			sa >> sanum;
			if(sanum < 4)
				adata[c][2] = 0;
			else if(sanum < 25)
				adata[c][2] = 1;
			else if(sanum < 50)
				adata[c][2] = 2;
			else
				adata[c][2] = 3;
			nu >> adata[c][0];
		}
		
		unsigned int max;
	      
		if (!pl.isValid()){
			if (!verbose)
				cout << "Warning: Invalid PDB file found "<<inputFile<<".\n";
			inputFile = "!";
			continue;
		}
	      
		if (!nmr)
			max = 1;
		else {
			max = pl.getMaxModels();
			if ( max == 0){
				max = 1;
				if ( !verbose )
					cout <<"Warning: the file "<<inputFile<<" probably is not an nmr structure.\n";
			}
		}

		totalanalized += 1;

		for (unsigned int ii = 1; ii <= max; ii++){
			pl.setModel(ii);
			Spacer sp;
			sp.load(pl);
			if (!pl.isValid()){
				if ( !verbose )
					cout << "Warning: Invalid PDB file found:"<<inputFile<<".\n";
				if ( ii == max ){
					inputFile = "!";
					totalanalized -= 1;
				}
				continue;
			}
			totalmodel += 1;
			int n , ini, end;
			end = (sp.sizeAmino());
			if(cae) {
				ini = 1;
				n = end - 1;
			}
			else {
				ini = 0;
				n = end;
			}
		  
			cout.setf(ios::fixed, ios::floatfield);
			
			vector<long unsigned int> count1;
			vector<long unsigned int> count2;
			//count1.assign ((binnums+1),0);//binnums+1 so contacts at <mind are also counted-> penalizes bad models (first bin [0,3.5) C VdW radious ~1.7A)
			//count2.assign ((binnums+1),0);
			
			if (verbose){
				cout <<"PDB file: "<<inputFile<<"\tChain: "<<chainID<<"\tMaxD: "<<maxd<<"\tMinD: "<<mind<<"\tBin: "<<bin<<"\n"; 
				cout <<"Type\tNumber\t";
			}
		      
		      
			if ( ca ){
				if(verbose)
					cout<<"Cas\t"<<chainID<<"\t"<<maxd<<"\t"<<mind<<"\t"<<bin<<"\n";
				int iii = 0;
				for ( int i = ini; i < n; i++ ){
					/*AminoAcid aa = sp.getAmino(i);
					if (aa.isMember(CA)){
						count1.assign ((binnums+1),0);
						for (int j = 0; j < (i + seqsep); j++){
							if (sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) <= (maxd) &&  aa[CA].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, binnums))++;
							else if (sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						for (int j = (i + 1 + seqsep); j < end; j++){
							if (sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) <= (maxd) &&  aa[CA].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, binnums))++;
							else if (sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						while(sp.getPdbNumberFromIndex(i) != adata[iii][0] && iii <= n)
							iii++;
						cout <<sp.getAmino(i).getType()<<"\t"<<sp.getPdbNumberFromIndex(i)<<"\t"<<adata[iii][1]<< "\t"<<adata[iii][2]<< "\t";
						for (unsigned int j = 0; j < count1.size(); j++)
							cout << j << " " << count1.at(j)<<"\t";
						cout  << seqsep << " seqsep" << endl;
					}*/

					if (sp.getAmino(i).isMember(CA)){
						count1.assign ((binnums+1),0);
						for (int j = ini; j < (i + seqsep) && j < n; j++){
							if (sp.getAmino(j).isMember(CA) && sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]) <= (maxd) &&  sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
							else if (sp.getAmino(j).isMember(CA) && sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						for (int j = (i + 1 + seqsep); j < n; j++){
							if (sp.getAmino(j).isMember(CA) && sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]) <= (maxd) &&  sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
							else if (sp.getAmino(j).isMember(CA) && sp.getAmino(i)[CA].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						while(sp.getPdbNumberFromIndex(i) != adata[iii][0] && iii <= n)
							iii++;
					}
				}
			}
			else if ( cb ){
				if(verbose)
					cout<<"Cbs\t"<<chainID<<"\t"<<maxd<<"\t"<<mind<<"\t"<<bin<<"\n";
				int iii = 0;
				for (int  i = ini; i < n; i++ ){
					AminoAcid aa = sp.getAmino(i);
					if (aa.getType() != "GLY" && aa.getSideChain().isMember(CB)){
						count1.assign ((binnums+1),0);
						for (int j = 0; j < (i + seqsep); j++){
							if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
								count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
							else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
								count1.at(binnums)++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]), mind, bin))++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						for (int j = (i + 1 + seqsep); j < end; j++){
							if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
								count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
							else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
								count1.at(binnums)++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]), mind, bin))++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa.getSideChain()[CB].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						cout <<sp.getAmino(i).getType()<<"\t"<<sp.getPdbNumberFromIndex(i)<<"\t";
						for (int j = 0; j < count1.size(); j++)
							cout << count1.at(j)<<"\t";
						cout << endl;
					}
					else if (aa.getType() == "GLY" && aa.isMember(CA)){
						count1.assign ((binnums+1),0);
						for (int j = 0; j < (i + seqsep); j++){
							if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa[CA].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa[CA].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
								count1.at(getbin(aa[CA].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
							else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa[CA].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
								count1.at(binnums)++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) <= (maxd) && aa[CA].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						for (int j = (i + 1 + seqsep); j < end; j++){
							if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa[CA].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa[CA].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
								count1.at(getbin(aa[CA].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
							else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa[CA].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
								count1.at(binnums)++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) <= (maxd) && aa[CA].distance(sp.getAmino(j)[CA]) >= (mind))
								count1.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
							else if(sp.getAmino(j).getType() == "GLY" && sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) < (mind))
								count1.at(binnums)++;
						}
						while(sp.getPdbNumberFromIndex(i) != adata[iii][0] && iii <= n)
							iii++;
						cout <<sp.getAmino(i).getType()<<"\t"<<sp.getPdbNumberFromIndex(i)<<"\t"<<adata[iii][1]<< "\t"<<adata[iii][2]<< "\t";
						for (int j = 0; j < count1.size(); j++)
							cout << count1.at(j)<<"\t";
						cout << endl;
					}
				}
			}
			else if ( bb ){
				if(verbose)
					cout<<"Cbs\t"<<chainID<<"\t"<<maxd<<"\t"<<mind<<"\t"<<bin<<"\n";
				int iii = 0;
				for (int i = ini; i < n; i++ ){
					AminoAcid aa = sp.getAmino(i);
					if (aa.getType() != "GLY" && aa.getSideChain().isMember(CB)){
						count1.assign ((binnums+1),0);
						for (int j = 0; j < (i + seqsep); j++){
							if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
								count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
							else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
								count1.at(binnums)++;
							else if(sp.getAmino(j).getType() == "GLY"){ 
								AminoAcid aa2 = sp.getAmino(j);
								aa2.patchBetaPositionGly();
								if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) >= (mind))
									count1.at(getbin(aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]), mind, binnums))++;
								else if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) < (mind))
									count1.at(binnums)++;
							}
						}
						for (int j = (i + 1 + seqsep); j < end; j++){
							if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
								count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
							else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
								count1.at(binnums)++;
							else if(sp.getAmino(j).getType() == "GLY"){ 
								AminoAcid aa2 = sp.getAmino(j);
								aa2.patchBetaPositionGly();
								if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) >= (mind))
									count1.at(getbin(aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]), mind, binnums))++;
								else if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) < (mind))
									count1.at(binnums)++;
							}
						}
						while(sp.getPdbNumberFromIndex(i) != adata[iii][0] && iii <= n)
							iii++;
						cout <<sp.getAmino(i).getType()<<"\t"<<sp.getPdbNumberFromIndex(i)<<"\t"<<adata[iii][1]<< "\t"<<adata[iii][2]<< "\t";
						for (int j = 0; j < count1.size(); j++)
							cout << count1.at(j)<<"\t";
						cout << endl;
					}
					else if (aa.getType() == "GLY"){// && aa.isMember(CA)){
						aa.patchBetaPositionGly();
						if(aa.getSideChain().isMember(CB)){
							count1.assign ((binnums+1),0);
							for (int j = 0; j < (i + seqsep); j++){
								if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
									count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
								else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
									count1.at(binnums)++;
								else if(sp.getAmino(j).getType() == "GLY"){ 
									AminoAcid aa2 = sp.getAmino(j);
									aa2.patchBetaPositionGly();
									if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) >= (mind))
										count1.at(getbin(aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]), mind, binnums))++;
									else if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) < (mind))
										count1.at(binnums)++; 
								}
							}
							for (int j = (i + 1 + seqsep); j < end; j++){
								if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) >= (mind))
									count1.at(getbin(aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]), mind, binnums))++;
								else if (sp.getAmino(j).getType() != "GLY" && sp.getAmino(j).getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) < (mind))
									count1.at(binnums)++;
								else if(sp.getAmino(j).getType() == "GLY"){ 
									AminoAcid aa2 = sp.getAmino(j);
									aa2.patchBetaPositionGly();
									if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) <= (maxd) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) >= (mind))
										count1.at(getbin(aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]), mind, binnums))++;
									else if(aa2.getSideChain().isMember(CB) && aa.getSideChain()[CB].distance(aa2.getSideChain()[CB]) < (mind))
										count1.at(binnums)++; 
								}
							}
							while(sp.getPdbNumberFromIndex(i) != adata[iii][0] && iii <= n)
								iii++;
							cout <<sp.getAmino(i).getType()<<"\t"<<sp.getPdbNumberFromIndex(i)<<"\t"<<adata[iii][1]<< "\t"<<adata[iii][2]<< "\t";
							for (int j = 0; j < count1.size(); j++)
								cout << count1.at(j)<<"\t";
							cout << endl;
						}
					}
				}
			}
			else if ( cae ){
				if(verbose)
					cout<<"HSE1\tHSE2\t"<<chainID<<"\t"<<maxd<<"\t"<<mind<<"\t"<<bin<<"\n";
				int iii = 0;
				for (int i = ini; i < n; i++ ){
					AminoAcid aa = sp.getAmino(i);
					AminoAcid aam = sp.getAmino(i-1);//aa minus 1
					AminoAcid aap = sp.getAmino(i+1);//aa plus 1
					int nums[3];
					nums[0] = sp.getPdbNumberFromIndex(i);
					nums[1] = sp.getPdbNumberFromIndex(i-1);
					nums[2] = sp.getPdbNumberFromIndex(i+1);
					if (aa.isMember(CA) && aam.isMember(CA) && aap.isMember(CA) && nums[1] == (nums[0] - 1) && nums[2] == (nums[0] + 1)){
						count1.assign ((binnums+1),0);
						count2.assign ((binnums+1),0);
						double pbeta[3];
						double dmag[2];
						vgVector3 <double> aac = aa[CA].getCoords();
						vgVector3 <double> aapc = aap[CA].getCoords();
						vgVector3 <double> aamc = aam[CA].getCoords();
						pbeta[0] = -(aac[0]+aac[0]-aamc[0]-aapc[0]);
						pbeta[1] = -(aac[1]+aac[1]-aamc[1]-aapc[1]);
						pbeta[2] = -(aac[2]+aac[2]-aamc[2]-aapc[2]);
						dmag[0] = -((pbeta[0]*aac[0]) + (pbeta[1]*aac[1]) + (pbeta[2]*aac[2]));
						dmag[1] = sqrt((pbeta[0]*pbeta[0])+(pbeta[1]*pbeta[1])+(pbeta[2]*pbeta[2]));
						for (int j = 0; j < (i + seqsep); j++){
							if(sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) <= (maxd) && aa[CA].distance(sp.getAmino(j)[CA]) >= (mind)){
								vgVector3 <double> aajc = sp.getAmino(j)[CA].getCoords();
								double tot = ((((pbeta[0])*(aajc[0]))+((pbeta[1])*(aajc[1]))+((pbeta[2])*(aajc[2]))+dmag[0])/dmag[3]);
								if(tot < 0)
									count1.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
								else
									count2.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
							}
							else if(sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) < (mind)){
								vgVector3 <double> aajc = sp.getAmino(j)[CA].getCoords();
								double tot = ((((pbeta[0])*(aajc[0]))+((pbeta[1])*(aajc[1]))+((pbeta[2])*(aajc[2]))+dmag[0])/dmag[3]);
								if(tot < 0)
									count1.at(binnums)++;
								else
									count2.at(binnums)++;
							}
						}
						for (int j = (i + 1 + seqsep); j < end; j++){
							if(sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) <= (maxd) && aa[CA].distance(sp.getAmino(j)[CA]) >= (mind)){
								vgVector3 <double> aajc = sp.getAmino(j)[CA].getCoords();
								double tot = ((((pbeta[0])*(aajc[0]))+((pbeta[1])*(aajc[1]))+((pbeta[2])*(aajc[2]))+dmag[0])/dmag[3]);
								if(tot < 0)
									count1.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
								else
									count2.at(getbin(aa[CA].distance(sp.getAmino(j)[CA]), mind, bin))++;
							}
							else if(sp.getAmino(j).isMember(CA) && aa[CA].distance(sp.getAmino(j)[CA]) < (mind)){
								vgVector3 <double> aajc = sp.getAmino(j)[CA].getCoords();
								double tot = ((((pbeta[0])*(aajc[0]))+((pbeta[1])*(aajc[1]))+((pbeta[2])*(aajc[2]))+dmag[0])/dmag[3]);
								if(tot < 0)
									count1.at(binnums)++;
								else
									count2.at(binnums)++;
							}
						}
						while(sp.getPdbNumberFromIndex(i) != adata[iii][0] && iii <= n)
							iii++;
						cout <<sp.getAmino(i).getType()<<"\t"<<sp.getPdbNumberFromIndex(i)<<"\t"<<adata[iii][1]<< "\t"<<adata[iii][2]<< "\t";
						for (int j = 0; j < count1.size(); j++)
							cout << count1.at(j)<< "\t" << count2.at(j)<<"\t";
						cout << endl;
					}
				}
				
			}
		}
		// reset variable to trigger break condition:
		inputFile = "!";
	}
	if ( verbose ){
		cout <<"Total file analized:"<<"\t"<<totalanalized<<"\n";
	if ( nmr )
		cout <<"Total nmr model analized:\t"<<totalmodel<<"\n";
	}
}
