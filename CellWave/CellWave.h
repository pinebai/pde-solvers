///
///  CellWave -- simulating biochemical waves in cells
///  
///    includes main headers required for CellWave
///

#ifndef CELLWAVE_H
#define CELLWAVE_H "CellWave.h"

#include "GenericReaction.h"
//#include "ReactionFactory.h"
//#include "GenericSolver.h"
//#include "SolverFactory.h"

#include "ParameterReader.h"
#include "Info.h"
#include "DPrintf.h"

#define dprintf XXX

//#include "ExplicitSolver.h"

namespace CellWave {

  void start( int argc, char **argv);
  void finish();
  
  enum DebugPrintNumber {
    BroadcastPrint=0, 
    PrintOut=1,
    LogPrint=2, 
    DebugPrint=3, 
    DetailedDebugPrint=4,
    MinMaxDebugPrint=5,
    DebugReaction=6,
    DebugSolver=7
  };

}; //end namespace CellWave

#endif
