#include <iostream>
#include <stdio.h>
#include <math.h>
//#include "Reactions.h"
#include "ParameterReader.h"
#include "CellWave.h"
#include "ReactionSlepchenko2Buffer.h"

using namespace CellWave;

int main(int argc, char **argv)
{
  
  try {
    CellWave::start(argc,argv);
    const int BPRINT=CellWave::BroadcastPrint;
    const int PRINT=CellWave::PrintOut;
    
    DPrintf(BPRINT," --------------------------------------------------------------------- \n");
    DPrintf(BPRINT,"  .. odeSlepchenko2Buffer ..\n");
    DPrintf(BPRINT,"      simulating a well mixed LiRinzel rxn w/ 2 buffers\n");
    DPrintf(BPRINT," --------------------------------------------------------------------- \n");
    if(argc <= 1) {
      std::cerr <<"Usage: odeSlepchenko2Buffer <parameter file name>\n";
      throw "no parameter file";
    }
    
    CellWave::Info data;
    CellWave::ParameterReader *pParameters = NULL;
    //
    //..Set parameters
    //
    if (argc>1) {
      assert( argv[1] != NULL );
      std::string paramFileName = argv[1];
      pParameters = new ParameterReader( paramFileName );
    }

    if( pParameters == NULL ) {
      assert( argv[0] != NULL );
      DPrintf(BPRINT,"**CellWave ERROR:: couldn't find/open the parameter file, exiting**\n");
      DPrintf(BPRINT,"\n  usage: %s <parameter file.par>\n", argv[0]);
      exit(-1);
    }
    
    ParameterReader &param = *pParameters;
    int logEvery;
    param.get( "name of grid file",    data.nameOfOGFile,  "" );
    param.get( "name of show file",    data.nameOfShowFile,"");
    param.get( "maximum timestep",     data.timeStepSize,  0.1);
    param.get( "number of timesteps",  data.numberOfTimeSteps, 1);
    param.get( "save frequency",       data.saveEveryNthFrame, 10);
    
    DPrintf(PRINT,".. Grid file=%s, Show file=%s\n", 
	    data.nameOfOGFile.c_str(), data.nameOfShowFile.c_str());
    DPrintf(PRINT,"..     dt=%8.4e, num. steps=%d, save every %d frame, Tmax=%8.4e\n",
	    data.timeStepSize, data.numberOfTimeSteps, 
	    data.saveEveryNthFrame, data.maximumTime());

    //..setup rxns

    ReactionSlepchenko2Buffer chem;
    bool ok=  chem.readParameterFile( param );

    //..setup the ode solver
    double x=0., y=0., z=0.;
    double p,c,h, b1, b2; // pn, cn, hn;
    double pFlux, cFlux, hFlux, b1Flux, b2Flux;
    double pFlux_prev, cFlux_prev, hFlux_prev, b1Flux_prev, b2Flux_prev;
    int     nsteps = data.numberOfTimeSteps;
    double  dt     = data.timeStepSize;
    double  tcomp  = 0.;
    
    chem.setInitialData( x,y,z, c, h, p, b1, b2);

    //..output file
    std::string outFileName;
    param.get("ode output filename", outFileName, "ode.dat");
    const char *cFileName = outFileName.c_str();
    FILE *fOutput = fopen( cFileName, "w");
    if (fOutput == NULL) {
      printf("ERROR -- file could not be opened\n");
      printf("         filename = <%s>\n", cFileName); 
      throw "file error";
    }
    fflush(0);
    printf("\n");
    printf("-- initial data--");
    printf("c= %8.3e, h= %8.3e, p= %8.3e, b1= %8.3e, b2= %8.3e\n", 
	   c, h, p, b1, b2);
    
    //
    //.. prev step for AB2 = first step is Euler
    //
    // react.computeFlux( 0., c, p, h,
    //	       cFlux_prev,pFlux_prev,hFlux_prev);
    
    for (int kstep=0; kstep <nsteps; ++ kstep) {
      tcomp = kstep*dt;
      
      chem.computeRHS( c,h,p, b1,b2, cFlux, hFlux, pFlux, b1Flux, b2Flux );
      if ( kstep % data.saveEveryNthFrame == 0 ) {
	fprintf(fOutput,
	 " %16.8f  %16.8f  %16.8f %16.8f %16.8f  %16.8f %16.8f %16.8f %16.8f\n", 
	 tcomp,c,h,b1,b2, cFlux, hFlux, b1Flux, b2Flux );
      }
      DPrintf(DebugSolver,
	      " %16.8f  %16.8f  %16.8f %16.8f %16.8f\n",tcomp,c,h,b1,b2);
      if(kstep== 0) {
	cFlux_prev  = cFlux;  // First step = FORWARD EULER
	hFlux_prev  = hFlux;
	b1Flux_prev = b1Flux;
	b2Flux_prev = b2Flux;
      }

      c  +=  .5*dt*( 3.*cFlux   - cFlux_prev);
      h  +=  .5*dt*( 3.*hFlux   - hFlux_prev);
      b1 +=  .5*dt*( 3.*b1Flux  - b1Flux_prev);
      b2 +=  .5*dt*( 3.*b2Flux  - b2Flux_prev);
      //p +=  .5*( 3.*pFlux - pFlux_prev);
      
      cFlux_prev = cFlux;
      hFlux_prev = hFlux;
      b1Flux_prev = b1Flux;
      b2Flux_prev = b2Flux;
      
    }; //end for 
    fclose(fOutput);
    printf("-- final   data--");
    printf("c= %8.3e, h= %8.3e, p= %8.3e, b1= %8.3e, b2= %8.3e\n", 
	   c, h, p, b1, b2);
  } catch ( ... ) {
    std::cout << "Exiting, error caught.\n";
  }
}






