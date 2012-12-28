#include <assert.h>
#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "CellWave.h"
//#include "getDiffusionDT.h"
#include "SolverLiRinzel.h"
#include "ReactionLiRinzelWagner.h"

int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture
  CellWave::start(argc,argv);

  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  DPrintf(BPRINT," --------------------------------------------------------------------- \n");
  DPrintf(BPRINT,"   CellWave -- simulation of biochemical waves in cells\n");
  DPrintf(BPRINT," --------------------------------------------------------------------- \n");

  CellWave::Info data;
  CellWave::ParameterReader *pParameters = NULL;

  //
  //..Set parameters
  //
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

  CellWave::ParameterReader &param = *pParameters;
  param.get( "name of grid file",    data.nameOfOGFile,  "" );
  param.get( "name of show file",    data.nameOfShowFile,"");
  param.get( "maximum timestep",     data.timeStepSize,  0.1);
  param.get( "number of timesteps",  data.numberOfTimeSteps, 1);
  param.get( "save frequency",	     data.saveEveryNthFrame, 10);

  DPrintf(PRINT,".. Grid file=%s, Show file=%s\n", 
	 data.nameOfOGFile.c_str(), data.nameOfShowFile.c_str());
  DPrintf(PRINT,"..     dt=%8.4e, num. steps=%d, save every %d frame, Tmax=%8.4e\n",
	 data.timeStepSize, data.numberOfTimeSteps, 
	 data.saveEveryNthFrame, data.maximumTime());
  
  //CellWave::ExplicitSolver solver( data );
  
  //
  //..Solver
  //
  CellWave::SolverLiRinzel solver( data );
  solver.setup(param);
  solver.initialData();
  solver.solve();
  
  //
  //..Finalize and output log data
  //
  Overture::finish();          

  DPrintf(PRINT, "----------------------------------------------------------------\n");
  DPrintf(PRINT,"Computation completed after %d steps at computational T=%g\n",
	  solver.getNumberOfTimeSteps (), solver.getTime());
  DPrintf(PRINT,"    Total CPU time %g sec, average CPU time/step = %g sec/step\n",
	  solver.getTotalWallTime(), solver.getTotalWallTime()/solver.getNumberOfTimeSteps());

  DPrintf(PRINT,"\n--> view the results by running:\n");
  DPrintf(PRINT,"     plotStuffx %s\n\n",data.nameOfShowFile.c_str());

  CellWave::finish();
  return 0;
    
}
