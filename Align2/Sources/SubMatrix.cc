// --*- C++ -*------x-----------------------------------------------------------
//
 
//
// Description:     Implement a standard substitution matrix.
//
// -----------------x-----------------------------------------------------------

#include <SubMatrix.h>

namespace Biopool
{

// CONSTRUCTORS:

SubMatrix::SubMatrix() : Substitution()
{ }


SubMatrix::SubMatrix(istream &is) : Substitution()
{
	is >> *this;
}


SubMatrix::SubMatrix(const SubMatrix &orig)
{
	copy(orig);
}


SubMatrix::~SubMatrix()
{ }


// OPERATORS:

SubMatrix&
SubMatrix::operator = (const SubMatrix &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


ostream&
operator << (ostream &os, const SubMatrix &object)
{
	os << object.residues << endl;
	Substitution::pWriteDoubleVector(os, object.residuescores);
	return os;
}


istream&
operator >> (istream &is, SubMatrix &object)
{
	is >> object.residues;
	Substitution::pReadDoubleVector(is, object.residuescores);
	object.buildscore(object.residues, object.residuescores);
	return is;
}


// MODIFIERS:

void
SubMatrix::copy(const SubMatrix &orig)
{
	Substitution::copy(orig);

	residuescores.clear();
	for (unsigned int n = 0; n < orig.residuescores.size(); n++)
		residuescores.push_back(orig.residuescores[n]);

	residues = orig.residues;
}


SubMatrix*
SubMatrix::newCopy()
{
	SubMatrix *tmp = new SubMatrix(*this);
	return tmp;
}

} // namespace
