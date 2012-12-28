/// Brief description: Factory class that creates solver instances
///

#include <math.h>
#include <stdio.h>
#include <iostream.h>

#include "CellWave.h"
#include "GenericSolver.h"
#include "SolverFactory.h"

//..known solvers: add your new solver here
#include "SolverLiRinzel.h"
#include "Solver2Buffer.h"
#include "SolverSlepchenko2Buffer.h"

using namespace CellWave;

//.. constructor
SolverFactory::
SolverFactory()
{
  //..do nothing
}

SolverFactory::
~SolverFactory()
{
  //..do nothing
}

GenericSolver* SolverFactory::
getSolver( const std::string &name,  CellWave::Info &data )
{
  //..NEW SOLVERS MUST BE ADDED HERE
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  GenericSolver *pSolver = NULL;

  DPrintf(PRINT,"SolverFactory::getSolver <%s>\n", name.c_str());
  if ( name == CellWave::SolverLiRinzel::solverType() ) {
    DPrintf(PRINT," --- SolverFactory --> inst. LiRinzel\n");

    pSolver = new CellWave::SolverLiRinzel(data);
  } 
  else if (name == CellWave::SolverSlepchenko2Buffer::solverType()){
    DPrintf(PRINT," --- SolverFactory --> inst. Slepchenko\n");
    pSolver = new CellWave::SolverSlepchenko2Buffer(data);
  }
  else if (name == CellWave::Solver2Buffer::solverType()){
    DPrintf(PRINT," --- SolverFactory --> inst. 2Buffer\n");
    pSolver = new CellWave::Solver2Buffer(data);
  }
  if (pSolver==NULL) {
    DPrintf(BPRINT,"ERROR::cellWaveMain not able to instantiate solver '%s'\n",
	    name.c_str());
    throw "error, no solver instantiated";
  }
  assert(pSolver!= NULL );
  DPrintf(PRINT,">> Solver %s instantiated <<\n",name.c_str());
  DPrintf(CellWave::DebugSolver,">> Solver %s instantiated <<\n",
	  name.c_str());

  return( pSolver );
}
