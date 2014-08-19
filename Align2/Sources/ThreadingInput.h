 

#ifndef __ThreadingInput_H__
#define __ThreadingInput_H__

#include <Substitution.h>
#include <iostream>
#include <string>

namespace Biopool
{
/** @brief    Implement I/O objects for handling threading files.
 * 
* @Description  
* @This 
 **/
class ThreadingInput
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ThreadingInput();

	/// istream constructor.
	ThreadingInput(istream &is);

	/// Copy constructor.
	ThreadingInput(const ThreadingInput &orig);

	/// Destructor.
	virtual ~ThreadingInput();


// OPERATORS:

	/// Assignment operator.
	ThreadingInput& operator = (const ThreadingInput &orig);

	/// Output operator.
	friend ostream& operator << (ostream &os, const ThreadingInput &object);

	/// Input operator.
	friend istream& operator >> (istream &is, ThreadingInput &object);


// PREDICATES:

	/// Return the threading score.
	double score(int i, int j);

	/// Return the size of the object referred as the dimension of the matrix.
	virtual unsigned int size() const;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ThreadingInput &orig);

	/// Construct a new "deep copy" of this object.
	virtual ThreadingInput* newCopy();


// HELPERS:

	/// Helper function used to write a vector<vector> construct.
	template<class T> static void pWriteDoubleVector(ostream &os,
		vector< vector<T> > data);

	/// Helper function used to read a vector<vector> construct.
	template<class T> static void pReadDoubleVector(istream &is,
		vector< vector<T> > &data);


protected:


private:

// ATTRIBUTES:

	vector< vector<double> > residuescores;    ///< Similarity scores.

};

// -----------------------------------------------------------------------------
//                                ThreadingInput
// -----------------------------------------------------------------------------

// PREDICATES:

inline double
ThreadingInput::score(int i, int j)
{
	return residuescores[j][i];
}


inline unsigned int
ThreadingInput::size() const
{
	return residuescores.size();
}

} // namespace

#endif
