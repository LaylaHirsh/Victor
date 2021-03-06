// -*- C++ -*-------x-----------------------------------------------------------
//
//
// Description:     Calculate a frequency profile or PSSM using Henikoff
//                  weighting scheme. Some explanations can be found in:
//
//                  Guoli Wang, Roland L. Dunbrack jr.
//                  Scoring profile-to-profile sequence alignments.
//                  Institute for Cancer Research, Fox Chase Cancer Center,
//                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
//
// -----------------x-----------------------------------------------------------

#include <HenikoffProfile.h>
#include <ctime>
namespace Biopool
{

// CONSTRUCTORS:

HenikoffProfile::HenikoffProfile() : Profile()
{ }


HenikoffProfile::HenikoffProfile(const HenikoffProfile &orig) : Profile(orig)
{
	copy(orig);
}


HenikoffProfile::~HenikoffProfile()
{ }


// OPERATORS:

HenikoffProfile&
HenikoffProfile::operator = (const HenikoffProfile &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// MODIFIERS:

void
HenikoffProfile::copy(const HenikoffProfile &orig)
{
	Profile::copy(orig);

	aliWeight.clear();
	for (unsigned int i = 0; i < orig.aliWeight.size(); i++)
	{
		vector<double> tmp;
		for (unsigned int j = 0; j < orig.aliWeight[0].size(); j++)
			tmp.push_back(orig.aliWeight[i][j]);
		aliWeight.push_back(tmp);
	}
}


HenikoffProfile*
HenikoffProfile::newCopy()
{
	HenikoffProfile *tmp = new HenikoffProfile(*this);
	return tmp;
}


// HELPERS:

void                  //suggested: use cLen=25 to save computational time
HenikoffProfile::pCalculateWeight(Alignment &ali, unsigned int cLen)
{	/*struct tm* newtime;
	time_t t;
	time(&t);newtime=localtime(&t);*/
	//cout << "blocco A1 inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;
	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";


	// --------------------------------------------------
	// 1. Calculate the positions of the first aminoacids
	//    for every sequence
	// --------------------------------------------------

	vector<unsigned int> firstAminoPos;

	firstAminoPos.push_back(0); // master sequence

        cout<<"target"<<ali.getTarget();
        cout<<"\ntemplate"<<ali.getTemplate(0)<<"\n";
	for (unsigned int j = 0; j < (numSeq - 1); j++) // other sequences
		for (unsigned int i = 0; i < seqLen; i++){
                    
                        if (ali.getTemplatePos(i, j) != '-')
			{
				firstAminoPos.push_back(i);
				break;
			}
                }

	/*time(&t);newtime=localtime(&t);
	cout << "blocco A2 inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/
	// --------------------------------------------------
	// 2. Calculate the positions of the last aminoacids
	//    for every sequence
	// --------------------------------------------------

	vector<unsigned int> lastAminoPos;

	lastAminoPos.push_back(seqLen-1); // master sequence

	for (unsigned int j = 0; j < (numSeq - 1); j++) // other sequences
		for (unsigned int i = seqLen; i > 0; i--)
			if (ali.getTemplatePos(i-1, j) != '-')
			{
				lastAminoPos.push_back(i-1);
				break;
			}


	// --------------------------------------------------
	// 3. Calculate weights for master sequence
	// --------------------------------------------------
	/*time(&t);newtime=localtime(&t);
	cout << "blocco A3 inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/
	vector<double> wSeq;
	vector<unsigned int> seqSubset;

	for (unsigned int i = 0; i < seqLen; i++)
		if (ali.getTargetPos(i) == '-')
			wSeq.push_back(0.00);
		else
		{
			// Calculate the subset of sequences

			seqSubset.clear();

			seqSubset.push_back(0);
			for (unsigned int k = 0; k < (numSeq - 1); k++)
				if (ali.getTemplatePos(i, k) != '-')
					seqSubset.push_back(k+1);

			unsigned int subsetSize = seqSubset.size();


			// Calculate Cleft and Cright

			unsigned int Cleft = firstAminoPos[seqSubset[0]];
			unsigned int Cright = lastAminoPos[seqSubset[0]];

			for (unsigned int k = 1; k < subsetSize; k++)
			{
				if (firstAminoPos[seqSubset[k]] > Cleft)
					Cleft = firstAminoPos[seqSubset[k]];
				if (lastAminoPos[seqSubset[k]] < Cright)
					Cright = lastAminoPos[seqSubset[k]];
			}

			double sum = 0.00;

			for (unsigned int p = Cleft; p <= Cright; p++)
			{
				// Calculate the number of different aminoacids

				unsigned int Ndiff = 0;

				for (unsigned int index = 0; index < 20; index++)
				{
					if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == residue_indices[index]))
					{
						Ndiff++;
						continue;
					}
					for (unsigned int k = 1; k < subsetSize; k++)
						if (ali.getTemplatePos(p, seqSubset[k]-1) == residue_indices[index])
						{
							Ndiff++;
							break;
						}
				}

				// Calculate the number of the same aminoacid

				unsigned int n = 0;

				if (seqSubset[0] == 0)
					n++;
				for (unsigned int k = 1; k < subsetSize; k++)
					if (ali.getTemplatePos(p, seqSubset[k]-1) == ali.getTargetPos(p))
						n++;


				sum += (1 / (double)(Ndiff * n));
			}

			wSeq.push_back((1 / (double)(Cright - Cleft + 1)) * sum);
		}
	aliWeight.push_back(wSeq);
	/*time(&t);newtime=localtime(&t);
	cout << "blocco A4 inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/

	// --------------------------------------------------
	// 4. Calculate weights for other sequences
	// --------------------------------------------------
	//cout << "	step: "<<numSeq-1<<"\n";
    //4.1 calculate seqSubset, Cleft, Crigth ,Ndiff for each position, outside main loop
	//this is a time optimization of section 4. Francesco Lovo 2012
	vector <pair<int, int> > Cvalues;         //Cleft and Cright for each subSeq
	vector <vector<unsigned int> > allseqSubseq;
	vector <vector<unsigned int> > allNdiff;  //for each position i, we store Ndiff
											 // for each column j between Cleft and Cright
	bool henikoffReduction=false;
	for (unsigned int i = 0; i < seqLen; i++)
	{
		// Calculate the subset of sequences
		seqSubset.clear();

		if (ali.getTargetPos(i) != '-')
			seqSubset.push_back(0);

		for (unsigned int k = 0; k < (numSeq - 1); k++)
			if (ali.getTemplatePos(i, k) != '-')
				seqSubset.push_back(k+1);

		unsigned int subsetSize = seqSubset.size();

		// Calculate Cleft and Cright
		unsigned int Cleft = firstAminoPos[seqSubset[0]];
		unsigned int Cright = lastAminoPos[seqSubset[0]];

		for (unsigned int k = 1; k < subsetSize; k++)
		{
			if (firstAminoPos[seqSubset[k]] > Cleft)
				Cleft = firstAminoPos[seqSubset[k]];
			if (lastAminoPos[seqSubset[k]] < Cright)
				Cright = lastAminoPos[seqSubset[k]];
		}
		//this is violation to Henikoff formula, but used to save some computational time
		if (i-Cleft>cLen)  {(Cleft=i-cLen);henikoffReduction=true;}
		if (Cright-i>cLen) {(Cright=i+cLen);henikoffReduction=true;}

		//copy to global variables
		Cvalues.push_back(make_pair(Cleft,Cright));
		allseqSubseq.push_back(seqSubset);

		// Calculate the number of different aminoacids
		vector<unsigned int> pDiff;

		for (unsigned int p = Cleft; p <= Cright; p++)
		{
			unsigned int Ndiff = 0;

			for (unsigned int index = 0; index < 20; index++)
			{
				if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == residue_indices[index]))
				{
					Ndiff++;
					continue;
				}
				for (unsigned int k = 1; k < subsetSize; k++)
					if (ali.getTemplatePos(p, seqSubset[k]-1) == residue_indices[index])
					{
						Ndiff++;
						break;
					}
			}
			pDiff.push_back(Ndiff);
		}
		allNdiff.push_back(pDiff);
	}
	if (henikoffReduction) cout<<"henikoff analisys reduced!\n";      //just a worning

	for (unsigned int j = 0; j < (numSeq - 1); j++)
	{
		wSeq.clear();
		for (unsigned int i = 0; i < seqLen; i++)
			if (ali.getTemplatePos(i, j) == '-')
				wSeq.push_back(0.00);
			else
			{
				unsigned int Cleft = Cvalues[i].first;
				unsigned int Cright =Cvalues[i].second;
				seqSubset= allseqSubseq[i];
				unsigned int subsetSize = seqSubset.size();
				double sum = 0.00;

				for (unsigned int p = Cleft; p <= Cright; p++)
				{
					// Calculate the number of different aminoacids
					unsigned int Ndiff = allNdiff[i][p-Cleft];

					// Calculate the number of the same aminoacid
					unsigned int n = 0;
					if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == ali.getTemplatePos(p, j)))
						n++;
					for (unsigned int k = 1; k < subsetSize; k++)
						if (ali.getTemplatePos(p, seqSubset[k]-1) == ali.getTemplatePos(p, j))
							n++;

					sum += (1 / (double)(Ndiff * n));
				}
				wSeq.push_back((1 / (double)(Cright - Cleft + 1)) * sum);
			}
		aliWeight.push_back(wSeq);
	}
	//4.1 fine
	/*
	for (unsigned int j = 0; j < (numSeq - 1); j++)
	{
		wSeq.clear();*/

		/*if (j%50==0){
			time(&t);newtime=localtime(&t);
			cout << "	check point, step: " <<j<<" "<< newtime->tm_hour << "/" << newtime->tm_min << endl;}*/

	/*	for (unsigned int i = 0; i < seqLen; i++)
			if (ali.getTemplatePos(i, j) == '-')
				wSeq.push_back(0.00);
			else
			{
				// Calculate the subset of sequences

				seqSubset.clear();

				if (ali.getTargetPos(i) != '-')
					seqSubset.push_back(0);

				for (unsigned int k = 0; k < (numSeq - 1); k++)
					if (ali.getTemplatePos(i, k) != '-')
						seqSubset.push_back(k+1);

				unsigned int subsetSize = seqSubset.size();


				// Calculate Cleft and Cright

				unsigned int Cleft = firstAminoPos[seqSubset[0]];
				unsigned int Cright = lastAminoPos[seqSubset[0]];

				for (unsigned int k = 1; k < subsetSize; k++)
				{
					if (firstAminoPos[seqSubset[k]] > Cleft)
						Cleft = firstAminoPos[seqSubset[k]];
					if (lastAminoPos[seqSubset[k]] < Cright)
						Cright = lastAminoPos[seqSubset[k]];
				}
        */      /*this is violation to Henikoff formula, but to save some computational time
                 *  we limit Crigth-Cleft<=50*/
		/*		if ((Cright-Cleft)>CLEN){
					((i-CLEN/2<0) ? (Cleft=0) : (Cleft=i-CLEN/2));
					((i+CLEN/2>(seqLen-1)) ? (Cright=(seqLen-1)) : (Cright=i+CLEN/2));
				}
		*/		//if ((i==1)&&(j==0)) cout<<"Crigth - CLeft: "<<Cright-Cleft<<"\n";
                /*done*/
		/*		double sum = 0.00;

				for (unsigned int p = Cleft; p <= Cright; p++)
				{
					// Calculate the number of different aminoacids

					unsigned int Ndiff = 0;

					for (unsigned int index = 0; index < 20; index++)
					{
						if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == residue_indices[index]))
						{
							Ndiff++;
							continue;
						}
						for (unsigned int k = 1; k < subsetSize; k++)
							if (ali.getTemplatePos(p, seqSubset[k]-1) == residue_indices[index])
							{
								Ndiff++;
								break;
							}
					}

					// Calculate the number of the same aminoacid

					unsigned int n = 0;

					if ((seqSubset[0] == 0) && (ali.getTargetPos(p) == ali.getTemplatePos(p, j)))
						n++;
					for (unsigned int k = 1; k < subsetSize; k++)
						if (ali.getTemplatePos(p, seqSubset[k]-1) == ali.getTemplatePos(p, j))
							n++;


					sum += (1 / (double)(Ndiff * n));
				}
				wSeq.push_back((1 / (double)(Cright - Cleft + 1)) * sum);
			}

		aliWeight.push_back(wSeq);
	}*/
	/*time(&t);newtime=localtime(&t);
	cout << "blocco A4 finito " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/

/*	for (unsigned int i = 0; i < seqLen; i++)
	{
		for (unsigned int j = 0; j < numSeq; j++)
			cout << aliWeight[j][i] << "  ";
		cout << endl;
	}*/
}

void
HenikoffProfile::pCalculateRawFrequency(vector<double> &freq, double &freqGap,
	Alignment &ali, unsigned int i)
{
	if (aminoAcidOneLetterTranslator(ali.getTargetPos(i)) != XXX)
		freq[aminoAcidOneLetterTranslator(ali.getTargetPos(i))] += aliWeight[0][i];
	else
		freqGap++;

	for (unsigned int j = 0; j < (numSeq - 1); j++)
		if (aminoAcidOneLetterTranslator(ali.getTemplatePos(i, j)) != XXX)
			freq[aminoAcidOneLetterTranslator(ali.getTemplatePos(i, j))] += aliWeight[j+1][i];
		else
			freqGap++;
}


void
HenikoffProfile::pConstructData(Alignment &ali)
{	/*struct tm* newtime;
	time_t t;
	time(&t);newtime=localtime(&t);
	cout << "blocco A inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/
	if (!gap)
	{
		ali.purgeTargetInsertions();
		seqLen = ali.getTarget().size();
	}
	/*time(&t);newtime=localtime(&t);
	cout << "blocco B inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/
	pCalculateWeight(ali);
	/*time(&t);newtime=localtime(&t);
	cout << "blocco C inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/
	gapFreq.reserve(seqLen);
	for (unsigned int i = 0; i < seqLen; i++)
	{
		vector<double> tmp;
		for (AminoAcidCode i = ALA; i <= TYR; i++)
			tmp.push_back(0.00);
		profAliFrequency.push_back(tmp);
		gapFreq.push_back(0.00);
	}
	/*time(&t);newtime=localtime(&t);
	cout << "blocco D inizia " << newtime->tm_hour << "/" << newtime->tm_min << endl;*/
	profAliFrequency.reserve(seqLen);
	for (unsigned int i = 0; i < seqLen; i++)
	{
		pCalculateRawFrequency(profAliFrequency[i], gapFreq[i], ali, i);

		double frequencySum = 0.00;
		for (AminoAcidCode j = ALA; j <= TYR; j++)
			frequencySum += profAliFrequency[i][j];

		for (AminoAcidCode j = ALA; j <= TYR; j++)
			profAliFrequency[i][j] /= frequencySum;
	}


/*	for (unsigned int i = 0; i < seqLen; i++)
	{
		for (unsigned int j = 0; j < 20; j++)
			cout << profAliFrequency[i][j] << "  ";
		cout << endl;
	}*/
	setSeq(ali.getTarget());
}

} // namespace
