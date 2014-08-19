// --*- C++ -*------x---------------------------------------------------------
//
// Description:     This class implemenents a standard substitution matrix.
// 
// -----------------x-------------------x-------------------x-----------------

#include <Blosum.h>

// CONSTRUCTORS:

Blosum::Blosum() : Substitution()
{ }


Blosum::Blosum(const Blosum& orig) : Substitution(orig),
    residuescores(orig.residuescores), residues(orig.residues)
{
  copy(orig);
}


Blosum::Blosum(istream& is) : Substitution()
{
  is >> *this;
}


Blosum::~Blosum()
{ }


// OPERATORS:

Blosum& Blosum::operator = (const Blosum& orig)
{
  if (&orig != this)
    copy(orig);
  return *this;
}


ostream& 
operator << (ostream& os, const Blosum& object)
{
  os << object.residues << endl;
  Substitution::pWriteDoubleVector(os, object.residuescores);
  return os;
}


istream& 
operator >> (istream& is, Blosum& object)
{
  is >> object.residues; 
  Substitution::pReadDoubleVector(is, object.residuescores);
  object.buildscore(object.residues, object.residuescores); 
  return is;
}


// MODIFIERS:

void 
Blosum::copy(const Blosum& orig)
{
  Substitution::copy(orig); // copy score matrix
  residuescores = orig.residuescores;
  residues = orig.residues;
}




