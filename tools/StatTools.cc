/**
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "StatTools.h"

/**
 * @Description get Random Number
 * @param double low, double high 
* @return double
*/
double getRandomNumber( double low, double high )
{
	double interval=max<double>( low, high )-min<double>(low, high );
	double delta=min<double>( low, high );
	double randomNumber=((double)rand())/(((double)RAND_MAX)+((double)1));
	return randomNumber*interval+delta;
}

/**
 * @Description get Gaussian Random Number
 * @param double low, double high double seed
* @return double
*/
double getGaussianRandomNumber( double low, double high, double seed )
{
	double sum=0.0;
	for ( long i=0; i<static_cast<long>(seed); ++i )
	{
		double rand=getRandomNumber(0,1);
		sum += rand;
	}
	double interval=max<double>( low, high )-min<double>(low, high );
	return interval*(sum-seed/2)/(seed/2);
}



