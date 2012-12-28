#include <assert.h>
#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "CellWave.h"
#include "GenericSolver.h"
#include "SolverFactory.h"

//#include "getDiffusionDT.h"
//#include "SolverLiRinzel.h"
//#include "ReactionLiRinzelWagner.h"

//#include "SolverSlepchenko2Buffer.h"
//#include "ReactionSlepchenko2Buffer.h"

int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);
  CellWave::start(argc,argv);

  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  DPrintf(BPRINT," --------------------------------------------------------------------- \n");
  DPrintf(BPRINT,"   CellWave -- simulation of biochemical waves in cells\n");
  DPrintf(BPRINT," --------------------------------------------------------------------- \n");

  CellWave::Info data;
  CellWave::ParameterReader *pParameters = NULL;

  //*
  //*..Set parameters
  //*
  if (argc>1) {
    assert( argv[1] != NULL );
    std::string paramFileName = argv[1];
    pParameters = new CellWave::ParameterReader( paramFileName );
  }

  if( pParameters == NULL ) {
    assert( argv[0] != NULL );
    DPrintf(BPRINT,"**CellWave ERROR:: couldn't find/open the parameter file, exiting**\n");
    DPrintf(BPRINT,"\n  usage: %s <parameter file.par>\n", argv[0]);
    exit(-1);
  }
  DPrintf(BPRINT,"    ***  parameter file = %s\n", argv[1]);
  DPrintf(BPRINT," --------------------------------------------------------------------- \n");

  CellWave::ParameterReader &param = *pParameters;
  param.get( "name of grid file",    data.nameOfOGFile,  "" );
  param.get( "name of show file",    data.nameOfShowFile,"");
  param.get( "maximum timestep",     data.timeStepSize,  0.1);
  param.get( "number of timesteps",  data.numberOfTimeSteps, 1);
  param.get( "save frequency",	     data.saveEveryNthFrame, 10);
  param.get( "log frequency",	     data.logEveryNthFrame,  
	                             data.saveEveryNthFrame);
  param.get( "flush frequency",      data.flushFrequency, 10);

  DPrintf(PRINT,".. Grid file=%s, Show file=%s\n", 
	 data.nameOfOGFile.c_str(), data.nameOfShowFile.c_str());
  DPrintf(PRINT,"..     dt=%8.4e, num. steps=%d, save every %d frame, Tmax=%8.4e\n",
	 data.timeStepSize, data.numberOfTimeSteps, 
	 data.saveEveryNthFrame, data.maximumTime());
  
  //CellWave::ExplicitSolver solver( data );
  
  //* .... Solver Factory code
  CellWave::SolverFactory solverFactory;
  CellWave::GenericSolver *pSolver =NULL;
  std::string paramFileType;
  paramFileType = param.get("parameter file type");

  pSolver = solverFactory.getSolver( paramFileType, data );
  
  assert(pSolver!= NULL );
  DPrintf(CellWave::DebugSolver,">> CellWaveMain--Solver %s instantiated <<\n",
	  pSolver->getSolverName().c_str());

  //* .... Execute solver
  pSolver->setup(param);
  pSolver->printReactionParameters( CellWave::DebugReaction );

  pSolver->initialData();
  pSolver->solve();
  pSolver->finish();
  
  //
  //*..Finalize and output log data
  //

  DPrintf(PRINT, "----------------------------------------------------------------\n");
  DPrintf(PRINT,"Computation completed after %d steps at computational T=%g\n",
	  pSolver->getNumberOfTimeSteps (), pSolver->getTime());
  DPrintf(PRINT,"    Total CPU time %g sec, average CPU time/step = %g sec/step\n",
	  pSolver->getTotalWallTime(), pSolver->getTotalWallTime()/pSolver->getNumberOfTimeSteps());

  DPrintf(PRINT,"\n--> view the results by running:\n");
  DPrintf(PRINT,"     plotStuffx %s\n\n",data.nameOfShowFile.c_str());

  delete pSolver;
  CellWave::finish();
  Overture::finish();          
  return 0;
    
}
