 

#ifndef __ScoringP2Ssec_H__
#define __ScoringP2Ssec_H__

#include <ScoringScheme.h>
#include <Profile.h>
#include <string>

namespace Biopool
{
/** @brief   Calculate score for profile to sequence alignment with
*                  secondary structure.
 * 
* @Description  
* @This 
 **/

class ScoringP2Ssec : public ScoringScheme
{
public:

// CONSTRUCTORS:
  /// Default constructor.
  ScoringP2Ssec(Blosum* sub, Blosum* sec, AlignmentData *a, 
		Profile* p);
  /// Copy constructor.
  ScoringP2Ssec(const ScoringP2Ssec& orig);
  /// Destructor.
  virtual ~ScoringP2Ssec();
  
// OPERATORS:
  /// Assignment operator.
  ScoringP2Ssec& operator = (const ScoringP2Ssec& orig);
  
// PREDICATES:
  /// Calculate scores to create Matrix values.
  double scoring(int i, int j);
  
// MODIFIERS:
  /// Copies orig object to this object ("deep copy").
  virtual void copy(const ScoringP2Ssec& orig);
  /// Constructs a new ("deep") copy of this object.
  virtual ScoringP2Ssec* newCopy();
  /// Reverse second component (sequence or profile)
  virtual void reverse();

protected:
  
private:

// ATTRIBUTES:
  Blosum* secstr;	///< Secondary structure substitution matrix.
  string seq1;          ///< First sequence.
  string seq2; 	        ///< Second sequence.
  Profile* profile;	///< (First) profile.
  string sec1;          ///< First  secondary structure.
  string sec2;	        ///< Second secondary structure.
};

}
#endif
