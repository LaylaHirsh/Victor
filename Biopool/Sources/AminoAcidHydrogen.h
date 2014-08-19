 

#ifndef _AMINOACIDHYDROGEN_H_
#define	_AMINOACIDHYDROGEN_H_

#include <AminoAcid.h>
#include <AminoAcidCode.h>
#include <string.h>
#include <map>
#include <list>



namespace Biopool {
/** @brief class implements hydrogens.
 * 
* @Description Includes methods that allow to hydrogens to aminoacids.
 * */
class AminoAcidHydrogen {
public: 

    static void loadParam(string inputFile);
    static void setHydrogen(AminoAcid* aa, bool verbose);
    
private:
   
    static map<AminoAcidCode,vector<vector<string> > > paramH;
    
};
    

} // namespace


#endif	/* _AMINOACIDHYDROGEN_H_ */

