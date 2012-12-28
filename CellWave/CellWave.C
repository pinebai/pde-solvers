
#include "CellWave.h"
#include <stdio.h>
#include <stdarg.h>

namespace CellWave {

void start( int argc, char **argv)
{
  //..initialize debug print facility
  DPrintfCreate();
  DPrintfSetDefault( PrintOut );

  DPrintfOpen(LogPrint,            "OUT/runcellwave.log");
  DPrintfOpen(DebugPrint,          "OUT/debug.out");
  DPrintfOpen(DetailedDebugPrint,  "OUT/debugDetailed.out");
  DPrintfOpen(MinMaxDebugPrint,    "OUT/debugMinMax.out");
  DPrintfOpen(DebugReaction,       "OUT/debugReaction.out");
  DPrintfOpen(DebugSolver,         "OUT/debugSolver.out");

  DPrintf(LogPrint,  "...................CellWave run log............\n");
  DPrintf(DebugPrint,"--debug output--\n");
  DPrintf(DetailedDebugPrint,"--detailed debug output--\n");
  DPrintf(MinMaxDebugPrint,"--Min/Max debug output--\n");
  DPrintf(DebugReaction,"-- debug **reactions** output--\n");
  DPrintf(DebugSolver,"--debug **solver** output--\n");

}

void finish()
{
  CellWave::DPrintfDestruct();
}

}



