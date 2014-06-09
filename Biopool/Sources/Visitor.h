/**
*  @Class:              Visitor
* @Author:             Silvio Tosatto
*  @Project Name:       Victor
**/ 
#ifndef _VISITOR_H_
#define _VISITOR_H_

// Includes:
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
 
  class Group;
  class SideChain;
  class AminoAcid;
  class Spacer;
   /** @brief Base class implementing the visitor pattern  
 * 
* @Description   This class uses the implicit copy operator. This is mostly used for implementing optimizers not contained in  the Biopool  module (e.g. EnergyCalculator, SideChainPlacement).
 **/
  class EnergyVisitor  {
  public: 
    
    // CONSTRUCTORS/DESTRUCTOR:
    EnergyVisitor() {};
    virtual ~EnergyVisitor() {};  
    
    // MODIFIERS:
    virtual void PrepareGroup(Group& group) = 0;
    virtual void PrepareSideChain(SideChain& node) = 0;
    virtual void PrepareAminoAcid(AminoAcid& node) = 0;
    virtual void PrepareSpacer(Spacer& node) = 0;
    
  protected:
    
  private:

  };

/** @brief base class  Optimizacion Patter. 
 * 
* @Description   This class uses the implicit copy operator.This is mostly used for implementing optimizers not contained in  the Biopool  module (e.g. EnergyCalculator, SideChainPlacement).
 **/
  class OptimizationVisitor  {
  public: 
    
    // CONSTRUCTORS/DESTRUCTOR:
    OptimizationVisitor() {};
    // 
  virtual ~OptimizationVisitor() {};  
    
    // MODIFIERS:
    virtual void PrepareGroup(Group& group) = 0;
    virtual void PrepareSideChain(SideChain& node) = 0;
    virtual void PrepareAminoAcid(AminoAcid& node) = 0;
    virtual void PrepareSpacer(Spacer& node) = 0;
    
  protected:
    
  private:
    
  };
  
} // namespace
#endif //_SAVER_H_
