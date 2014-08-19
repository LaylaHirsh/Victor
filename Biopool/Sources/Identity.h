/**
* @Class:           Identity
* @Project name:    Victor
* @Description:     Implements object id.
 */

#ifndef __IDENTITY_H__
#define __IDENTITY_H__

// Includes:

#include <iostream>
#include <string>
#include <Debug.h>

class Identity {
public:

  Identity(const string& name="", long number=0);
  Identity(const Identity& orig);

  /* OPERATORS */

  /** Assigment. */
  Identity& operator = (const Identity& orig);

  /** Comparison. */
  bool operator == (const Identity& other) const;
  bool operator != (const Identity& other) const { 
    return !((*this) == other); };

  /** Cast to int. */
  operator long() const;

  /** Cast to string. */
  //  operator string() const;

  /** Cast to string&. */
  operator const string&() const;

  /** Set a new name. */
  void setName(string _name);

  /** Set a new number. */
  void setNumber(long _number);



  /* FRIENDS */

  friend ostream& operator << (ostream& os, const Identity& id);

  friend istream& operator >> (istream& is, Identity& id);

private:
  long   number;
  string name;

  static long counter;
};

// ---------------------------------------------------------------------------
//                                  Identity
// -----------------x-------------------x-------------------x-----------------

inline
Identity::Identity(const string& s, long n) : name(s)
{
  if (n != 0)
    {
      number = n;
    }
  else
    {
      number = ++counter;
    }
}

inline
Identity::Identity(const Identity& orig) 
  : name(orig.name)
{
  number = ++counter;
}

inline
Identity& 
Identity::operator = (const Identity& orig)
{
  name   = orig.name;

  return *this;
}

inline
bool 
Identity::operator == (const Identity& other) const
{
  return number == other.number;
}

inline
Identity::operator long() const
{
  return number;
}



inline
Identity::operator const string&() const
{
  return name;
}

inline
void 
Identity::setName(string _name)
{
  name = _name;
}

inline
void 
Identity::setNumber(long _number)
{
  number = _number;
  if (_number >= counter)
    counter = _number;
}

// ---------------------------------------------------------------------------
//                                   Friends
// -----------------x-------------------x-------------------x-----------------

inline
ostream& 
operator << (ostream& os, const Identity& id)
{
  return os << "(Id " << id.name << " " << id.number << " ) ";
}

inline
istream& 
operator >> (istream& is, Identity& id)
{
  string tag;

  is >> tag;
  
  if (tag != "(Id")
    {
      DEBUG_MSG("Identity::operator >> : bad input.");

      is.clear(ios::badbit);
      return is;
    }

  return is >> id.name >> id.number >> tag;
}


#endif /* __IDENTITY_H__ */

