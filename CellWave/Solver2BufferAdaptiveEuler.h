/// Brief description: solver mainloop for CellWave simulator, Li Rinzel model
///               
///

#ifndef SOLVER_LI_RINZEL_H
#define SOLVER_LI_RINZEL_H

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

namespace CellWave {
  class ReactionLiRinzelWagner;

  class SolverLiRinzel : public GenericSolver {
  public:
    enum { cc=0, pc=1, hc=2,
	   bcfast=3, bcslow=4};  // hardwired components

    SolverLiRinzel( Info &data_ ) 
      : GenericSolver( data_ )   
      {
	//default
      }
    ~SolverLiRinzel()
      {
	//default
      };

    virtual bool readParameterFile( CellWave::ParameterReader &param);
    virtual void setup();
    virtual void setup( CellWave::ParameterReader &param);
    virtual void initialData();
    virtual void solve(); 

    virtual real getTime();               // in base class
    virtual real getTotalWallTime();      // in base class
    virtual int  getNumberOfTimeSteps(); // in base class
    
    //.. solver data  
    CellWave::ReactionLiRinzelWagner *pChem;

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
