

#ifndef __Blosum_H__
#define __Blosum_H__

#include <iostream>
#include <string>
#include <Substitution.h>

/** @brief  This class implemenents a standard substitution matrix
 * 
* @Description  
* @This 
 **/
class Blosum : public Substitution 
{
public:
    
// CONSTRUCTORS:
  /// Default constructor.
  Blosum();
  /// Copy constructor.
  Blosum(const Blosum& orig);
  /// Constructor from istream.
  Blosum(istream& is);
  /// Destructor.
  virtual ~Blosum();
  
// OPERATORS:
  /// Assignment operator.
  Blosum& operator = (const Blosum& orig);
  /// Output operator.
  friend ostream& operator << (ostream& os, const Blosum& object);
  /// Input operator.
  friend istream& operator >> (istream& is, Blosum& object);

// PREDICATES:
  /// Implementation of abstract class method.
  virtual string getResidues() const { return residues; }
  /// Return dimension of matrix.
  virtual unsigned int size() const { return residuescores.size(); }

protected:
  
// MODIFIERS:
  /// Copies orig object to this object ("deep copy").
  virtual void copy(const Blosum& orig);

private:

// ATTRIBUTES:
  vector<vector<int> > residuescores; ///< Similarity scores.
  string residues;                    ///< Alphabet of allowed characters.

};

#endif /* __Blosum_H__ */

