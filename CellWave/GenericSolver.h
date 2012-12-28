/// Brief description: solver mainloop for CellWave simulator
/// 
///

#ifndef GENERIC_SOLVER_H
#define GENERIC_SOLVER_H "GenericSolver.h"

//#include "GenericReaction.h"
//#include "ReactionFactory.h"
//#include "Info.h"
//#include "DPrintf.h"
#include "Overture.h"
#include "CompositeGridOperators.h"
#include "Ogshow.h"
#include "CellWave.h"
#include "NucleusGridFunction.h"

#include "Probes.h"
#include "FluxBC.h"

namespace CellWave {

  class FluxBCData {
  public:
    enum fluxBCType {
      noFluxBC          = 0,
      interpolateFluxBC = 1
    }; 

    FluxBCData() :
      physicalBoundaryID(1), fluxBoundaryID(2), useFluxBC(true)
      {
	//..do nothing
      }

    ~FluxBCData() { };

    int getPhysicalBoundaryID() {return physicalBoundaryID; };
    int getFluxBoundaryID() { return fluxBoundaryID; };

  private:
    FluxBCData( FluxBCData &X);             //do not use copy constructor
    FluxBCData &operator=( FluxBCData &X ); //do not use = on this class
  public:
    int physicalBoundaryID, fluxBoundaryID;
    FluxBC bc;
    bool useFluxBC; //..debug flag, allow =false and no FluxBC
    
  };

  class GenericSolver {
  private:
    std::string solverName;
  protected:
    void        setSolverName( const std::string &name ) { solverName = name; };
  public:
    //GenericSolver( ) 
    //{ };
    GenericSolver( Info &data_ ) 
      : data( data_ ), 
        solverName("Generic"),
      numberOfIP3Species(1),  // ip3 is one species, no ip3 buffers etc
      useDebugShow(false) // debug flags
      { };

    ~GenericSolver() 
      { };

    //virtual void setup( GenericReaction *pChem_) = 0; //NOT USED
    //{
    ////do nothing
    //  }

    virtual void setup( CellWave::ParameterReader &param)
      {
      	//do nothing
      }

    virtual void initialData()
    {
      //do nothing
    };

    virtual void solve()
      {
      	//do nothing
      };

    virtual void finish() 
      {
	genericFinish();
      }

    virtual real getTime()              { return tcomp; };
    virtual real getTotalWallTime()     { return totalWallTime; };
    virtual int  getNumberOfTimeSteps() { return ktime; };
    virtual bool
      readParameterFile( CellWave::ParameterReader &param);

    //.. base class services
    std::string getSolverName() { return solverName; };
    void genericSetup( GenericReaction *pChem_, 
		       CellWave::ParameterReader &param);
    void genericFinish(); // close files etc
    void genericInitialData();
    void collectProbeData( int ktime, double tcomp, realCompositeGridFunction &q);
    void updateModel();
    //void 
    void outputStepStart(int ktime_, double tcomp_);

    void printReactionParameters( int iPrintChannel );
    
    realArray &getNucleusMaskArray( int igrid ) { return( nucleusGF[igrid] ); };

    void turnOffFluxBoundaryConditions() { fluxBCData.useFluxBC=false;};
    void turnOnFluxBoundaryConditions() { fluxBCData.useFluxBC=true;};

    void setupFluxBoundaryConditions(  CompositeGrid          &cg,
				       Interpolant            &interpolant,
				       CompositeGridOperators &operators);
    void setFluxBCCoefficient( double newFlux ); //should be stored in FluxBC?

    void applyFluxBoundaryConditions( realCompositeGridFunction &u,     //to
				      int iSolutionComponent);
    void applyNoFluxBoundaryConditions(realCompositeGridFunction &u,
				       int iSolutionComponent);

    //.. solver data
    Info                   &data;
    CompositeGrid           cg;
    Interpolant            interpolant;
    CompositeGridOperators operators;
    Ogshow show;
    realCompositeGridFunction q;
    realCompositeGridFunction f, fp;

    const int numberOfIP3Species;      //used in solver loop, so as not to update IP3 2x 
    //realCompositeGridFunction ip3Next; // not used: ip3 at next time level

    GenericReaction *pChem;

    //..debug flags
    bool useDebugShow;

    //..probes
    Probes probes;
    FluxBCData fluxBCData;

    //.... nucleus data
    NucleusGridFunction         nucleusClass;
    doubleCompositeGridFunction nucleusGF;

    //.. parameters
    real tcomp;   /// computational time
    int  ktime;   /// timestep number

    //..Info for the nucleus
    

    //..log info

    real wallTime;
    int iWallTimeSteps;
    real totalWallTime;  
    
  };
}; //end namespace CellWave
#endif
