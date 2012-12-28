/// Brief description: solver mainloop for CellWave simulator, Li Rinzel model
///               
///

#ifndef SOLVER_SLEPCHENKO_2_BUFFER_H
#define SOLVER_SLEPCHENKO_2_BUFFER_H

#include <assert.h>
#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

//#include "CellWave.h"
#include "GenericReaction.h"
#include "ReactionFactory.h"
#include "GenericSolver.h"
#include "Info.h"
#include "getDiffusionDT.h"
#include "ReactionSlepchenko2Buffer.h"

namespace CellWave {
  //class ReactionSlepchenko2Buffer;

  class SolverSlepchenko2Buffer : public GenericSolver {
  public:
    static std::string solverType() //{ return "Slepchenko2Buffer"; }
      {
	return ReactionSlepchenko2Buffer::rxnType();
      };

    enum { cc=0, pc=1, hc=2,
	   b1c=3, b2c=4};  // hardwired components. must match RxnSlepchenko...

    SolverSlepchenko2Buffer( Info &data_ ) 
      : GenericSolver( data_ )   
      {
	this->setSolverName(  ReactionSlepchenko2Buffer::rxnType()  );	
      }
    ~SolverSlepchenko2Buffer()
      {
	//default
      };

    bool readParameterFile( CellWave::ParameterReader &param);
    void setup();
    void setup( CellWave::ParameterReader &param);    
    void initialData();
    void solve(); 
    //void finish();

    virtual real getTime();               // in base class
    virtual real getTotalWallTime();      // in base class
    virtual int  getNumberOfTimeSteps(); // in base class
    
    //.. solver data  
    CellWave::ReactionSlepchenko2Buffer *pChem;

    // -- IN BASE CLASS
    //CompositeGrid           cg;
    //Interpolant            interpolant;
    //CompositeGridOperators operators;
    //Ogshow show;
    //realCompositeGridFunction q;
    //GenericReaction *pChem;
    //real tcomp;   /// computational time
    //int  ktime;   /// timestep number
    

  };
}; //end namespace CellWave
#endif
