/**
* @Class:           -
* @Base class:      -
* @Derived classes: - 
* @Author:          Damiano Piovesan
* @Project name:    -
* @Description:     Translator: PDB names to internal one-word-code and 
*                  vice versa. Provides some simple predicates dealing with 
*                  one-word-code. 
* @ Copyright:       This file contains information from the Bioinformatics 
*                  Template Library (BTL).
*
*                  Copyright (C) 1997,1998 Birkbeck College, Malet Street, 
*                  London WC1E 7HX, U.K. (classlib@mail.cryst.bbk.ac.uk)
*/

#ifndef __Nucleotide_CODE_H__
#define __Nucleotide_CODE_H__

// Includes

#include <string>
#include <Debug.h>

/**
 */
enum NucleotideCode 
{ 
  ADENINE=0,
  THYMINE,
  CYTOSINE,
  GUANINE,
  URACIL,
  XX,  //  Corresponds to unknown Nucleotide type 
  Nucleotide_CODE_SIZE // number of Nucleotide types
};


/**
   true, if Nucleotide type name is known 
   magic code XXX corresponds to an unknown Nucleotide code
*/
inline
bool
isKnownNucleotide(NucleotideCode code)
{
  return !(code == XX);
}


/** Translate string into Nucleotide code enum. */
inline
NucleotideCode
nucleotideOneLetterTranslator(char name)
{
  DUMP(name);
  if (name == ' ')
    {
      return XX;
    }
  if (name == 'X')
    {
      return XX;
    }
  if (name == 'A')
    {
      return ADENINE;
    }
  else if (name == 'T')
    {
      return THYMINE;
    }
  else if (name == 'C')
    {
      return CYTOSINE;
    }
  else if (name == 'G')
    {
      return GUANINE;
    }
  else if (name == 'U')
    {
      return URACIL;
    }

  return XX; // Nucleotide type name is unknown
}


/** Translate string into Nucleotide code enum. */
inline
NucleotideCode
nucleotideThreeLetterTranslator(const string& name)
{
  DUMP(name);
  if (name == "")
    {
      return XX;
    }
  if (name == "X")
    {
      return XX;
    }
  if (name == "A")
    {
      return ADENINE;
    }
  else if (name == "T")
    {
      return THYMINE;
    }
  else if (name == "C")
    {
      return CYTOSINE;
    }
  else if (name == "G")
    {
      return GUANINE;
    }
  else if (name == "U")
    {
      return URACIL;
    }
   else if (name == "DA")
    {
      return ADENINE;
    }
   else if (name == "DT")
    {
      return THYMINE;
    }
  else if (name == "DC")
    {
      return CYTOSINE;
    }
  else if (name == "DG")
    {
      return GUANINE;
    }
  else if (name == "DU")
    {
      return URACIL;
    }
  
  return XX; // Nucleotide type name is unknown
}

/**
 * @Description Returns the corresponding three letter code as a string.
 * @param aminoacidCode 
 * @return string
 */
inline string nucleotideThreeLetterTranslator(NucleotideCode code) 
{
  switch (code)
    {
    case XX:
      return "X";
    case ADENINE:
      return "A";
    case THYMINE:
      return "T";
    case CYTOSINE:
      return "C";
    case GUANINE:
      return "G";
    case URACIL:
      return "U";
    case Nucleotide_CODE_SIZE: 
      ERROR("NucleotideTranslator(NucleotideCode code): unknown code", exception);
    }

  ERROR("NucleotideTranslator(NucleotideCode code): unknown code", eXXXception);

  return "X";
}


/**
 * @Description Returns the corresponding three letter code as a string.
 * @param one letter code (char)
 * @return string
 */
inline string nucleotideOneLetter2ThreeLetter(char oneLetter)
{
  NucleotideCode code;
  code = nucleotideOneLetterTranslator(oneLetter);
  return nucleotideThreeLetterTranslator(code);
}
/**
 * @Description Returns the corresponding one letter code.
 * @param three letter code (string)
 * @return char
 */
inline char nucleotideThreeLetter2OneLetter(const string& threeLetter)
{
  NucleotideCode code;
  code = nucleotideThreeLetterTranslator(threeLetter);
  return nucleotideOneLetterTranslator(code);
}
/**
 * @Description Returns one aminoacid code.
 * @param aminoacid code 
 * @return aminoacid code
 */
inline NucleotideCode& operator++(NucleotideCode& ac, int)
{
  return ac = ( (ac == XX) ? ADENINE : NucleotideCode(ac+1) );
}
/**
 * @Description Returns one aminoacid code.
 * @param aminoacid code 
 * @return aminoacid code
 */
inline NucleotideCode& operator--(NucleotideCode& ac, int)
{
  return ac = ( (ac == ADENINE) ? XX : NucleotideCode(ac-1) );
}

#endif /* __Nucleotide_CODE_H__ */

