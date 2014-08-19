 

#ifndef __SubMatrix_H__
#define __SubMatrix_H__

#include <Substitution.h>

namespace Biopool
{
/** @brief     Implement a standard substitution matrix.
 * 
* @Description  
* @This 
 **/
class SubMatrix : public Substitution
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	SubMatrix();

	/// istream constructor.
	SubMatrix(istream &is);

	/// Copy constructor.
	SubMatrix(const SubMatrix &orig);

	/// Destructor.
	virtual ~SubMatrix();


// OPERATORS:

	/// Assignment operator.
	SubMatrix& operator = (const SubMatrix &orig);

	/// Output operator.
	friend ostream& operator << (ostream &os, const SubMatrix &object);

	/// Input operator.
	friend istream& operator >> (istream &is, SubMatrix &object);


// PREDICATES:

	/// Implementation of abstract class method.
	virtual string getResidues() const;

	/// Return the dimension of the matrix.
	virtual unsigned int size() const;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const SubMatrix &orig);

	/// Construct a new "deep copy" of this object.
	virtual SubMatrix* newCopy();


protected:


private:

// ATTRIBUTES:

	vector< vector<int> > residuescores;    ///< Similarity scores.
	string residues;                        ///< Alphabet of allowed characters.

};

// -----------------------------------------------------------------------------
//                                 SubMatrix
// -----------------------------------------------------------------------------

// PREDICATES:

inline string
SubMatrix::getResidues() const
{
	return residues;
}


inline unsigned int
SubMatrix::size() const
{
	return residuescores.size();
}

} // namespace

#endif
