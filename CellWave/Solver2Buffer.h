/// Brief description: solver mainloop for CellWave simulator, Li Rinzel model w/ 2 buffers
///               
///  endogeneous & exogeneous buffer, IP3 diffusion, IP3 dependence in Ca dynamics
///

#ifndef SOLVER_2_BUFFER_H
#define SOLVER_2_BUFFER_H

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
#include "Reaction2Buffer.h"

namespace CellWave {
  //class Reaction2Buffer;

  class Solver2Buffer : public GenericSolver {
  public:
    static std::string solverType() 
      { 
	return Reaction2Buffer::rxnType();
      };

    enum { cc=0, pc=1, hc=2,
	   b1c=3, b2c=4};  // hardwired components. must match RxnSlepchenko...

    Solver2Buffer( Info &data_ ) 
      : GenericSolver( data_ ), ip3InfluxBoundaryID( 3 ), useIP3InfluxBC( false )
      {
	this->setSolverName(  Reaction2Buffer::rxnType()  );	
      }
    ~Solver2Buffer()
      {
	//default
      };

    bool readParameterFile( CellWave::ParameterReader &param);
    void setup();
    void setup( CellWave::ParameterReader &param);    
    void initialData();
    void solve(); 
    //void finish();

    int getIP3InfluxBoundaryID() { return ip3InfluxBoundaryID; }
    int getIP3InfluxBoundaryGridID() { return ip3InfluxGridNumber; }

    virtual real getTime();               // in base class
    virtual real getTotalWallTime();      // in base class
    virtual int  getNumberOfTimeSteps(); // in base class

    //.. IP3 volume source terms for IP3 influx in 2D
    void inline addIP3VolumeSources( const int igrid,
				     const realArray &xArray, 
				     const realArray &qArray, 
				     realArray &rhsArray);
    
    //.. solver data  
    CellWave::Reaction2Buffer *pChem;

    int    ip3InfluxBoundaryID;
    bool   useIP3InfluxBC;
    int    ip3InfluxGridNumber;

    double xInflux_p;
    double yInflux_p;
    double zInflux_p;
    double influxRadius_p;
    double influxRate_p;

    // -- IN BASE CLASS
    //CompositeGrid           cg;
    //Interpolant            interpolant;
    //CompositeGridOperators operators;
    //Ogshow show;
    //realCompositeGridFunction q;
    //GenericReaction *pChem;
    //real tcomp;   /// computational time
    //int  ktime;   /// timestep number
    //int physicalBoundaryID, fluxBoundaryID;
    

  };

//---------------------------INLINE FUNCTIONS-----------------------------------------------------

void inline Solver2Buffer::
addIP3VolumeSources( const int igrid,
		     const realArray &xArray, 
		     const realArray &qArray, 
		     realArray &rhsArray)
{
  //..Setup
  assert( pChem != NULL );
  CellWave::GenericReaction &chem = *pChem; // shorthand
  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  assert( (0<= igrid)  && (igrid< cg.numberOfComponentGrids()));
  MappedGrid & mg = cg[igrid];
  const IntegerArray & d   = mg.dimension();
  const IntegerArray & gir = mg.gridIndexRange();
  const int nd             = cg.numberOfDimensions();
  const int ncomp          = chem.getNumberOfSpecies();
  
  
}




}; //end namespace CellWave
#endif
