 

#ifndef __Substitution_H__
#define __Substitution_H__

#include <Debug.h>
#include <iostream>
#include <string>
#include <vector>

namespace Biopool
{
/** @brief    Base class for deriving substitution matrices.
 * 
* @Description  
* @This 
 **/
class Substitution
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Substitution();

	/// Copy constructor.
	Substitution(const Substitution &orig);

	/// Destructor.
	virtual ~Substitution();


// OPERATORS:

	/// Assignment operator.
	Substitution& operator = (const Substitution &orig);

	/// Output operator.
	friend ostream& operator << (ostream &os, const Substitution &object);

	/// Input operator.
	friend istream& operator >> (istream &is, Substitution &object);


// PREDICATES:

	/// Dummy implementation.
	virtual string getResidues() const = 0;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Substitution &orig);

	/// Construct a new "deep copy" of this object.
	virtual Substitution* newCopy() = 0;

	/// Build scoring matrix from raw data.
	virtual void buildscore(const string &residues,
		const vector< vector<int> > &residuescores);


// HELPERS:

	/// Helper function used to write a vector<vector> construct.
	
        /*template<class T> static void pWriteDoubleVector(ostream &os,
		vector<vector<T> > data);
        There were problems compiling the previous fuction in a 64 bit architecture,
        but use of template was not necessary here, so we used the next one:*/
        static void pWriteDoubleVector(ostream &os,vector<vector<int> > data);
	/// Helper function used to read a vector<vector> construct.
	template<class T> static void pReadDoubleVector(istream &is,
		vector<vector<T> > &data);


// ATTRIBUTES:

	vector< vector<int> > score;    ///< Substitution score.


protected:


private:

};

} // namespace

#endif
