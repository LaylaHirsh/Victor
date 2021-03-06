/**
 
* @Description:       Base class for composite structures. 
*                     Implementing the Composite pattern.
*/

// Includes:
#include <Component.h>
#include <Debug.h>
#include <limits.h>
#include <float.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:
Component::Component(unsigned int mI, unsigned int mO) : Bond(mI, mO), 
  superior(NULL), components(), lowerBound(DBL_MAX-1,DBL_MAX-1,DBL_MAX-1), 
  upperBound(-(DBL_MAX-1),-(DBL_MAX-1),-(DBL_MAX-1)), modified(false)
{ }

Component::Component(const Component& orig) 
{
  PRINT_NAME;
  this->copy(orig);  
}

Component::~Component()
{
  PRINT_NAME; 
} 


// PREDICATES:

vgVector3<double> 
Component::getLowerBound(double dist)
{
  sync();
  vgVector3<double> tmp = lowerBound;
  for (unsigned int i = 0; i < 3; i++)
    tmp[i] -= dist;

  return lowerBound;
}

vgVector3<double> 
Component::getUpperBound(double dist)
{
  sync();
  vgVector3<double> tmp = upperBound;
  for (unsigned int i = 0; i < 3; i++)
    tmp[i] += dist;

  return upperBound;
}


// MODIFIERS:

void Component::connectIn(Component* c, unsigned int offset)
{
  ERROR("connectIn() undefined for this class.", exception);
}

void Component::connectOut(Component* c, unsigned int offset)
{
  ERROR("connectOut() undefined for this class.", exception);
}

Component* Component::unconnectIn()
{
  ERROR("unconnectIn() undefined for this class.", exception);
}

Component* Component::unconnectOut()
{
  ERROR("unconnectOut() undefined for this class.", exception);
}

void Component::copy(const Component& orig)
{
  // Attention: new elements have to be copied also in Monomer::copy()
  // since it does *NOT* invoke this function (segmentation fault if it does).

  PRINT_NAME; 

  Bond::copy(orig);

  // copy bounding box:
  lowerBound = orig.lowerBound;
  upperBound = orig.upperBound;

  // deep copy of components:

  components.clear();  // !!!!! <----- 15/12/00 Bug: does not release memory: 
                       //                            use delete

	for (unsigned int i = 0; i < orig.size(); i++)  // form new ones
    {
		//Component* c = orig.components[i]->clone();
		// ---
		Component* yyy = orig.components[i];//->clone();
		if ( 0==yyy )
			cerr <<"errore" <<endl;
		Component *c=yyy->clone();
		// ---
		c->superior = this;
		components.push_back(c);
    }

  superior = orig.superior;             // this is added as component 
  if (superior != NULL)
    superior->insertComponent(this);      // to its superior

  modified = orig.modified;
}


// OPERATORS:

Component& Component::operator=(const Component& orig)
{
  PRINT_NAME;
  if (&orig != this)
    copy(orig);
  return *this;
}

// HELPERS:

 