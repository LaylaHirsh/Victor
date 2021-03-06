// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Calculate scores for profile to profile alignment using
//                  sequence metric factor. Some explanations can be found in:
//
//                  William R. Atchley, Jieping Zhao, Andrew D. Fernandes, Tanja Druke
//                  Solving the protein sequence metric problem.
//                  Edited by Walter M. Fitch, University of California, Irvine, CA,
//                  and approved March 22, 2005 (received for review December 14, 2004).
//
// -----------------x-----------------------------------------------------------

#include <AtchleyCorrelation.h>

namespace Biopool
{

// CONSTRUCTORS:

AtchleyCorrelation::AtchleyCorrelation(Profile *pro1, Profile *pro2)
	: ScoringFunction(), pro1(pro1), pro2(pro2), offset(8.00)
{
	pLoadFactor();
}


AtchleyCorrelation::AtchleyCorrelation(Profile *pro1, Profile *pro2,
	double offset) : ScoringFunction(), pro1(pro1), pro2(pro2), offset(offset)
{
	pLoadFactor();
}


AtchleyCorrelation::AtchleyCorrelation(const AtchleyCorrelation &orig)
	: ScoringFunction(orig)
{
	copy(orig);
}


AtchleyCorrelation::~AtchleyCorrelation()
{ }


// OPERATORS:

AtchleyCorrelation&
AtchleyCorrelation::operator = (const AtchleyCorrelation &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
AtchleyCorrelation::scoringSeq(int i, int j)
{
	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";


	// Calculate vectors f1 and f2

	vector<double> f1, f2;
	double m1 = 0.00;
	double m2 = 0.00;

	for (unsigned int z = 0; z < 5 ; z++)
	{
		double freq1 = 0.00;
		double freq2 = 0.00;
		double s1 = 0.00;
		double s2 = 0.00;

		for (unsigned int k = 0; k < 20; k++)
		{
			freq1 = pro1->getAminoFrequency(residue_indices[k], (i-1));
			freq2 = pro2->getAminoFrequency(residue_indices[k], (j-1));

			s1 += (freq1 * factor[k][z]);
			s2 += (freq2 * factor[k][z]);
		}

		f1.push_back(s1);
		f2.push_back(s2);

		m1 += s1;
		m2 += s2;
	}

	m1 /= f1.size();
	m2 /= f2.size();


	// Calculate correlation between vectors f1 and f2

	double num = 0.00;
	double den = 0.00;

	for (unsigned int z = 0; z < 5; z++)
	{
		num += (((f1[z] - m1) * (f2[z] - m2)) / f1.size());
		den += ((((f1[z] - m1) * (f1[z] - m1)) * ((f2[z] - m2) * (f2[z] - m2))) / f1.size());
	}


	return offset * (num / sqrt(den));
}


// MODIFIERS:

void
AtchleyCorrelation::copy(const AtchleyCorrelation &orig)
{
	ScoringFunction::copy(orig);
	pro1 = orig.pro1->newCopy();
	pro2 = orig.pro2->newCopy();
	offset = orig.offset;

	double factor[20][5];
	for (unsigned int k = 0; k < 20; k++)
		for (unsigned int z = 0; z < 5; z++)
			factor[k][z] = orig.factor[k][z];
}


AtchleyCorrelation*
AtchleyCorrelation::newCopy()
{
	AtchleyCorrelation *tmp = new AtchleyCorrelation(*this);
	return tmp;
}


// HELPERS:

void
AtchleyCorrelation::pLoadFactor()
{
	string path = getenv("VICTOR_ROOT");
	if (path.length() < 3)
		cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;
	string atchleyFileName = path + "data/atchley.dat";

	ifstream atchleyFile(atchleyFileName.c_str());
	if (!atchleyFile)
		ERROR("Error opening Atchley metric file.", exception);


	double f;
	for (unsigned int k = 0; k < 20; k++)
		for (unsigned int z = 0; z < 5; z++)
		{
			atchleyFile >> f;
			factor[k][z] = f;
		}
}

} // namespace
