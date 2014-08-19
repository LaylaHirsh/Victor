 

#ifndef __PssmInput_H__
#define __PssmInput_H__

#include <Substitution.h>
#include <iostream>
#include <string>

namespace Biopool
{
/** @brief  Implement I/O objects for handling BLAST PSSM (Position
*                  Specific Score Matrix).
 * 
* @Description  
* @This 
 **/
class PssmInput
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	PssmInput();

	/// istream constructor.
	PssmInput(istream &is);

	/// Copy constructor.
	PssmInput(const PssmInput &orig);

	/// Destructor.
	virtual ~PssmInput();


// OPERATORS:

	/// Assignment operator.
	PssmInput& operator = (const PssmInput &orig);

	/// Output operator.
	friend ostream& operator << (ostream &os, const PssmInput &object);

	/// Input operator.
	friend istream& operator >> (istream &is, PssmInput &object);


// PREDICATES:

	/// Return the score of the aminoacid j in position i.
	double score(int i, int j);

	/// Return the size of the object referred as the dimension of the PSSM.
	virtual unsigned int size() const;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const PssmInput &orig);

	/// Construct a new "deep copy" of this object.
	virtual PssmInput* newCopy();


// HELPERS:

	/// Helper function used to write a vector<vector> construct.
	template<class T> static void pWriteDoubleVector(ostream &os,
		vector< vector<T> > data, vector<string> data1, vector<string> data2);

	/// Helper function used to read a vector<vector> construct.
	template<class T> static void pReadDoubleVector(istream &is,
		vector< vector<T> > &data, vector<string> &data1, vector<string> &data2);


protected:


private:

// ATTRIBUTES:

	vector< vector<double> > residuescores;    ///< PSSM scores.
	vector<string> allPosition;                ///< Companion variable for I/O class.
	vector<string> allAa;                      ///< Companion variable for I/O class.

};

// -----------------------------------------------------------------------------
//                                  PssmInput
// -----------------------------------------------------------------------------

// PREDICATES:

inline double
PssmInput::score(int i, int j)
{
	return residuescores[i][j];
}


inline unsigned int
PssmInput::size() const
{
	return residuescores.size();
}

} // namespace

#endif
