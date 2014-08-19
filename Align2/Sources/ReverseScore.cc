

#include <ReverseScore.h>

namespace Biopool
{

// CONSTRUCTORS:

ReverseScore::ReverseScore(Align *a)
{
	ali = a->newCopy();
	inv = a->newCopy();
	inv->getScoringScheme()->reverse();
}


ReverseScore::ReverseScore(const ReverseScore &orig)
{
	copy(orig);
}


ReverseScore::~ReverseScore()
{ }


// OPERATORS:

ReverseScore&
ReverseScore::operator = (const ReverseScore &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// MODIFIERS:

void
ReverseScore::copy(const ReverseScore &orig)
{
	ali = orig.ali;
	inv = orig.inv;
}


double
ReverseScore::getZScore(double &forward, double &reverse, unsigned int n)
{
	ali->recalculateMatrix();
	inv->recalculateMatrix();

	vector<double> score = inv->generateMultiMatchScore(n);

//	cout << ">>>>>>> " << ali->getScore() << endl;
//	for (unsigned int i = 0; i < score.size(); i++)
//		cout << score[i] << endl;

	forward = ali->getScore();
	reverse = average(score);
	return ((forward - reverse) / (standardDeviation(score) != 0 ? standardDeviation(score) : 1));
}

} // namespace
