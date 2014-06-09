/**
 * 
* @Class:             Monomer
* @Base Class(es):    Component
* @Derived Class(es): -
* @Containing:        -
* @ Author:            Silvio Tosatto
* @Project Name:      Victor
* @Description:       
*/ 
// Includes:
#include <Monomer.h>
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:
Monomer::Monomer(unsigned int mI, unsigned int mO) : Component(mI, mO) 
{ PRINT_NAME; }

Monomer::Monomer(const Monomer& orig)
{
  PRINT_NAME;
  this->copy(orig);  
}

Monomer::~Monomer()
{ 
  PRINT_NAME;
  if ( hasSuperior() )
    {
      getSuperior().removeComponent(this);
      setSuperior(NULL);
    }
} 


// PREDICATES:

// MODIFIERS:

void 
Monomer::copy(const Monomer& orig)
{
  PRINT_NAME; 
  

}

Component* 
Monomer::clone()
{
  return new Monomer;
}

// OPERATORS:

// HELPERS:

 