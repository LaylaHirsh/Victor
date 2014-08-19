 

#ifndef __SecSequenceData_H__
#define __SecSequenceData_H__

#include <AlignmentData.h>

namespace Biopool
{
/** @brief   Print alignment of two sequences considering also secondary
*                  structure.
 * 
* @Description  
* @This 
 **/
class SecSequenceData : public AlignmentData
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	SecSequenceData(int n, const string &seq1, const string &seq2,
		const string &sec1, const string &sec2, const string &name1 = "Seq1",
			const string &name2 = "Seq2");

	/// Copy constructor.
	SecSequenceData(const SecSequenceData &orig);

	/// Destructor.
	virtual ~SecSequenceData();


// OPERATORS:

	/// Assignment operator.
	SecSequenceData& operator = (const SecSequenceData &orig);


// PREDICATES:

	/// Return the sequence at position n of the vector.
	virtual string getSequence(int n);

	/// Calculate single match positions.
	virtual void calculateMatch(int i, int tbi, int j, int tbj);

	/// Reverse the strings of the vector.
	virtual void getMatch();

	/// Control if the strings of the vector are similar and print them.
	virtual void outputMatch(ostream &os, bool fasta = false);

	/// Generate and return an ensemble of suboptimal alignments.
	virtual Alignment& generateMatch(double score = 0.00);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const SecSequenceData &orig);

	/// Construct a new "deep copy" of this object.
	virtual SecSequenceData* newCopy();

	/// Set the sequence at position n of the vector.
	virtual void setSequence(string s, int n);


protected:


private:

// ATTRIBUTES:

	string seq1;    ///< Target sequence.
	string seq2;    ///< Template sequence.
	string sec1;    ///< Target secondary structure.
	string sec2;    ///< Template secondary structure.

};

} // namespace

#endif
