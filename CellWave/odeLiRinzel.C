#include <iostream>
#include <stdio.h>
#include <math.h>
//#include "Reactions.h"
#include "ParameterReader.h"
#include "CellWave.h"
#include "ReactionLiRinzelWagner.h"

using namespace CellWave;

int main(int argc, char **argv)
{
  
  try {
    CellWave::start(argc,argv);
    const int BPRINT=CellWave::BroadcastPrint;
    const int PRINT=CellWave::PrintOut;
    
    DPrintf(BPRINT," --------------------------------------------------------------------- \n");
    DPrintf(BPRINT,"   odeLiRinzel -- simulating a well mixed LiRinzel rxn\n");
    DPrintf(BPRINT," --------------------------------------------------------------------- \n");
    if(argc <= 1) {
      std::cerr <<"Usage: odeLiRinzel <parameter file name>\n";
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

    //..setup rxns

    ReactionLiRinzelWagner chem;
    bool ok=  chem.readParameterFile( param );

    //..setup the ode solver
    double x=0., y=0., z=0.;
    double p,c,h; // pn, cn, hn;
    double pFlux, cFlux, hFlux;
    double pFlux_prev, cFlux_prev, hFlux_prev;  
    int     nsteps = data.numberOfTimeSteps;
    double  dt     = data.timeStepSize;
    double  tcomp  = 0.;
    
    chem.setInitialData( x,y,z, c, h, p);

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
    printf("c= %8.3e, h= %8.3e, p= %8.3e\n", c, h, p);
    
    //
    //.. prev step for AB2 = first step is Euler
    //
    // react.computeFlux( 0., c, p, h,
    //	       cFlux_prev,pFlux_prev,hFlux_prev);
    
    for (int kstep=0; kstep <nsteps; ++ kstep) {
      tcomp = kstep*dt;
      
      chem.computeRHS( c,h,p, cFlux, hFlux, pFlux );
      fprintf(fOutput,
            " %16.8f  %16.8f  %16.8f %16.8f %16.8f\n", 
            tcomp,c,h, cFlux, hFlux );
      
      cFlux_prev = cFlux;  //FORWARD EULER
      hFlux_prev = hFlux;
      
      c +=  .5*dt*( 3.*cFlux - cFlux_prev);
      h +=  .5*dt*( 3.*hFlux - hFlux_prev);
      //p +=  .5*( 3.*pFlux - pFlux_prev);
      
      cFlux_prev = cFlux;
      hFlux_prev = hFlux;
      
    }; //end for 
    fclose(fOutput);
    printf("-- final   data--");
    printf("c= %8.3e, h= %8.3e, p= %8.3e\n", c, h, p);
  } catch ( ... ) {
    std::cout << "Exiting, error caught.\n";
  }
}






