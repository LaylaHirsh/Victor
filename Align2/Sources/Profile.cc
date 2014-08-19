// -*- C++ -*-------x-----------------------------------------------------------
//
// Description:     Calculate a frequency profile or PSSM.
//
// -----------------x-----------------------------------------------------------

#include <Profile.h>
#include <ctime>
namespace Biopool
{

// CONSTRUCTORS:

Profile::Profile() : profAliFrequency(), gapFreq(), seq(""), seqLen(0),
	numSeq(0), gap(false)
{ }


Profile::Profile(const Profile &orig)
{
	copy(orig);
}


Profile::~Profile()
{ }


// OPERATORS:

Profile&
Profile::operator = (const Profile &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// MODIFIERS:

void
Profile::copy(const Profile &orig)
{
	profAliFrequency.clear();
	for (unsigned int i = 0; i < orig.profAliFrequency.size(); i++)
	{
		vector<double> tmp;
		for (unsigned int j = 0; j < orig.profAliFrequency[0].size(); j++)
			tmp.push_back(orig.profAliFrequency[i][j]);
		profAliFrequency.push_back(tmp);
	}

	gapFreq.clear();
	for (unsigned int i = 0; i < orig.gapFreq.size(); i++)
		gapFreq.push_back(orig.gapFreq[i]);

	seq = orig.seq;
	seqLen = orig.seqLen;
	numSeq = orig.numSeq;
	gap = orig.gap;
}


Profile*
Profile::newCopy()
{
	Profile *tmp = new Profile(*this);
	return tmp;
}


void
Profile::setProfile(Alignment &ali)
{   struct tm* newtime;
    time_t t;
	pResetData();
	seqLen = ali.getTarget().size();
	numSeq = ali.size() + 1;
	time(&t);newtime=localtime(&t);
	cout << "ready for pConstructData " << newtime->tm_hour << "/" << newtime->tm_min << endl;
	pConstructData(ali);
}


void
Profile::setProfile(Alignment &ali, istream &is)
{   
	pResetData();
        cout<<"ali:  "<<ali.getTarget();
	seqLen = ali.getTarget().size();
	numSeq = ali.size() + 1;

	if (!gap)
	{
		ali.purgeTargetInsertions();
		seqLen = ali.getTarget().size();
	}

	for (unsigned int i = 0; i < seqLen; i++)
	{
		double f;

		vector<double> tmp;
		for (unsigned int j = 0; j < 20; j++)
		{
			is >> f;
			tmp.push_back(f);
		}
		profAliFrequency.push_back(tmp);

		is >> f;
		gapFreq.push_back(f);
	}
       
	setSeq(ali.getTarget());
}


void
Profile::reverse()
{
	vector< vector<double> > tmpPAF;
	for (unsigned int i = profAliFrequency.size(); i > 0; i--)
	{
		vector<double> tmp;
		for (unsigned int j = 0; j < profAliFrequency[i-1].size(); j++)
			tmp.push_back(profAliFrequency[i-1][j]);
		tmpPAF.push_back(tmp);
	}
	profAliFrequency.clear();
	profAliFrequency = tmpPAF;

	vector<double> tmpGF;
	for (unsigned int i = gapFreq.size(); i > 0; i--)
		tmpGF.push_back(gapFreq[i-1]);
	gapFreq.clear();
	gapFreq = tmpGF;

	string tmpS = "";
	for (unsigned int i = seq.length(); i > 0; i--)
		tmpS.push_back(seq[i-1]);
	seq = tmpS;
}


// HELPERS:

void
Profile::pCalculateRawFrequency(vector<double> &freq, double &freqGap,
	Alignment &ali, unsigned int i)
{
	if (aminoAcidOneLetterTranslator(ali.getTarget()[i]) != XXX)
		freq[aminoAcidOneLetterTranslator(ali.getTarget()[i])]++;
	else
		freqGap++;

	for (unsigned int j = 0; j < (numSeq - 1); j++)
		if (aminoAcidOneLetterTranslator(ali.getTemplate(j)[i]) != XXX)
			freq[aminoAcidOneLetterTranslator(ali.getTemplate(j)[i])]++;
		else
			freqGap++;
}


void
Profile::pConstructData(Alignment &ali)
{
	if (!gap)
	{
		ali.purgeTargetInsertions();
		seqLen = ali.getTarget().size();
	}
	gapFreq.reserve(seqLen);
	for (unsigned int i = 0; i < seqLen; i++)
	{
		vector<double> tmp;
		for (AminoAcidCode i = ALA; i <= TYR; i++)
			tmp.push_back(0.00);
		profAliFrequency.push_back(tmp);
		gapFreq.push_back(0.00);
	}
	for (unsigned int i = 0; i < seqLen; i++)
	{
		pCalculateRawFrequency(profAliFrequency[i], gapFreq[i], ali, i);
		for (AminoAcidCode j = ALA; j <= TYR; j++)
			profAliFrequency[i][j] /= (numSeq - gapFreq[i]);
	}
	setSeq(ali.getTarget());
}

void
Profile::pResetData()
{
	for (unsigned int i = 0; i < profAliFrequency.size(); i++)
		profAliFrequency[i].clear();
	profAliFrequency.clear();
	gapFreq.clear();
	seq = "";
	seqLen = 0;
	numSeq = 0;
}

} // namespace
