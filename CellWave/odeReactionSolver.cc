#include <iostream>
#include <stdio.h>
#include <math.h>
#include "OXCWave/Reactions.h"
#include "FlowToolkit/ParameterReader.h"

using namespace OXCWave;
using FlowToolkit::ParameterReader;

void 
changeParameters(FlowToolkit::ParameterReader &par,
		 OXCWave::Reactions           &Q);
bool 
isParameterFileCompatible( ParameterReader &par);

void 
changeParameters(FlowToolkit::ParameterReader &par,
		 OXCWave::Reactions           &Q)
{
  par.get("nu_L",    Q.nu_L,      5e-4);
  par.get("d_I",     Q.d_I,       0.025); // micro M
  par.get("k_p",     Q.k_p,       0.4);   // micro M
  par.get("nu_P",    Q.nu_P,      0.1);   // micro M
  par.get("I_s",     Q.I_s,       0.12);  // micro M
  par.get("C_er",    Q.C_er,      10.);   // micro M
  par.get("tau_0",   Q.tau_0,     4);     // sec
  par.get("eta",     Q.eta,       1.);    // 1/ micro M

  //tau_0      = 25;     // sec  --- DIFFERENT from Wagner
  par.get("d_act",   Q.d_act,     1.2);   // micro M

  par.get("diffCalcium", Q.diffCalcium, 300.);  // micro m^2/sec
  //diffCalcium =0.; //TRY THIS
  
  par.get("lambda", Q.lambda,      112.5); // 1/sec
  par.get("d_inh",  Q.d_inh,       1.5);   // micro M
  par.get("beta",   Q.beta,        0.053); // nondimensional
  par.get("n_h",    Q.n_h,         1.); // hill coeff

  // ..parameters NOT USED BY WAGNER et al
  par.get("diffIP3",Q.diffIP3,     0.); //300.;   // IP3 diffusion
  par.get("k_i",    Q.k_i,         0.); // 1/sec, IP3 degradation, 0=synthetic IP3S

  //..initial data
  par.get("calcium_0", Q.calcium_0,  0.1153); // micro M
  par.get("ip3_0",     Q.ip3_0,      0.12);   // micro M
  par.get("h_0",       Q.h_0,        0.93);   // nondimensional
}

bool 
isParameterFileCompatible( ParameterReader &params )
{
  std::string typeString="";
  params.get("parameter file type",typeString);
  bool ok= ( typeString== "oxcwave_dimensional_reactions");
  if( !ok ) {
    std::cout <<" --ERROR: incompatible parameter file, "
	      <<"type="<<typeString << std::endl;
  }
  return( ok );
}

int main(int argc, char **argv)
{
  try {
    if(argc <= 1) {
      std::cerr <<"Usage: testReactions.x <parameter file name>\n";
      throw "no parameter file";
    }
    FlowToolkit::ParameterReader params(argv[1]);
    if( !isParameterFileCompatible( params ) ) {
      throw "Unknown parameter file type";
    }
    
    OXCWave::Reactions react;
    
    double p,c,h; // pn, cn, hn;
    double pFlux, cFlux, hFlux;
    double pFlux_prev, cFlux_prev, hFlux_prev;  
    int     nsteps = 200;
    double  dt     = 0.1;
    double  tcomp  = 0.;
    
    const Reactions::IP3DistributionType 
      ip3dist=  Reactions::HomogeneousIP3;
    
    //..changes to parameters
    changeParameters(params, react);
    params.get("number of timesteps", nsteps, nsteps);
    params.get("dt",                  dt,     dt);

    printf("--parameters before rescaling\n");
    react.outputParameters();    //..output parameters
    double cascale=1., c0=1.;
    params.get("calcium rescale", cascale, cascale);
    react.invariantCalciumScaling( cascale,c0,c );
    printf("--parameters after rescaling\n");
    react.outputParameters();
    
    react.evaluateInitialData( ip3dist, c, h, p);
    //..output file
    std::string outFileName;
    params.get("ode output filename", outFileName, "ode.dat");
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
      
      react.computeFlux( tcomp, c,p,h,
			 cFlux, pFlux, hFlux );
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






