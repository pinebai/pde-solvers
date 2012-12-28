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
#include "ReactionLiRinzelWagner.h"

namespace CellWave {
  //class ReactionLiRinzelWagner;

  class SolverLiRinzel : public GenericSolver {
  public:
    static std::string solverType()  //{ return "LiRinzel"; }
      { 
	return ReactionLiRinzelWagner::rxnType();
      };

    enum { cc=0, pc=1, hc=2,
	   bcfast=3, bcslow=4};  // hardwired components

    SolverLiRinzel( Info &data_ ) 
      : GenericSolver( data_ )   
      {
	this->setSolverName(  ReactionLiRinzelWagner::rxnType()  );
      }
    ~SolverLiRinzel()
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
