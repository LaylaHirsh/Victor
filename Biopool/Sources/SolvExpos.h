/**
 * @Description A module to determine the solvent exposure/accessibility of residues in a
*                              protein fragment.
*/
#ifndef __SolvExpos_H__
#define __SolvExpos_H__

#include <Spacer.h>

namespace Biopool{

	enum SolvExpos {CORE, EXPOSED};

	SolvExpos getSolvExpos(Spacer &chain, const unsigned int tgt,
		const unsigned int start, const unsigned int end);

	vector<SolvExpos>* getSolvExposVec(Spacer &chain,
		const unsigned int tgtS, const unsigned int tgtE,
		const unsigned int envS, const unsigned int envE);

	double getSolvAccess(Spacer &chain, unsigned int tgt,
		unsigned int start, unsigned int end);

	vector<double> getSolvAccessVec(Spacer &chain,
		unsigned int tgtS, unsigned int tgtE,
		unsigned int envS, unsigned int envE);

} // namespace

#endif
