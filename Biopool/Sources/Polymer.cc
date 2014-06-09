/**
* @Class:             Polymer
* @Base Class(es):    -
* @Derived Class(es): -
* @Containing:        -
* @Author:            Silvio Tosatto
* @Project Name:      Victor
*/
// Includes:
#include <Polymer.h>
#include <Debug.h>
#include <AminoAcid.h> // for testing only

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:
Polymer::Polymer(unsigned int mI, unsigned int mO) : Component(mI, mO)
{ PRINT_NAME; }

Polymer::Polymer(const Polymer& orig)
{
  PRINT_NAME;
  this->copy(orig);  
}
Polymer::~Polymer()
{ 
  PRINT_NAME; 
  if ( hasSuperior() )
    getSuperior().removeComponent(this);
  setSuperior(NULL);
  while (size() > 0)
    deleteComponent(components[0]);
} 

// PREDICATES:

// MODIFIERS:
/**
 * @Description 
 * @param 
 */
void Polymer::copy(const Polymer& orig)
{
  PRINT_NAME; 
  Component::copy(orig);
}
/**
 * @Description 
 * @param 
 */
Component* Polymer::clone()
{
  return new Polymer;
}
/**
 * @Description 
 * @param 
 */
void Polymer::insertComponent(Component* c)
{
  PRECOND(c != NULL, exception);
  c->setSuperior(this);
  components.push_back(c);
} 
/**
 * @Description 
 * @param 
 */
void Polymer::removeComponent(Component* c)
{
  if ( c == NULL )
    return;
  for (unsigned int i = 0; i < size(); i++)
    if (components[i] == c)
      {
	components[i]->setSuperior(NULL);
	components.erase(components.begin()+i);
	return;
      }
  
  DEBUG_MSG("Polymer::removeComponent: component not found.");
} 

/**
 * @Description 
 * @param 
 */
void Polymer::removeComponentFromIndex(unsigned int i)
{
  if ( i > size())
    DEBUG_MSG("Index out of bound");
  components[i]->setSuperior(NULL);
  components.erase(components.begin()+i);
}
/**
 * @Description 
 * @param 
 */
void  Polymer::deleteComponent(Component* c)
{
  if ( c == NULL ) 
    return;
  for (unsigned int i = 0; i < size(); i++)
    if (components[i] == c)
      {
	
	components[i]->setSuperior(NULL);
	components.erase(components.begin()+i);
	delete c;
	return;
      }
  DEBUG_MSG("Polymer::removeComponent: component not found.");
} 


// OPERATORS:

// HELPERS:



