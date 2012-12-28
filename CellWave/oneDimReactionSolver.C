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
    
    double p,c,h, pmax, cmax, hmax; // pn, cn, hn;
    double *pFlux, *cFlux, *hFlux;
    double *pFlux_prev, *cFlux_prev, *hFlux_prev;  

    int     nsteps = 200;
    double  dt     = 0.1;
    double  tcomp  = 0.;
    
    const Reactions::IP3DistributionType 
      ip3dist=  Reactions::HomogeneousIP3;

    int ndim=100;
    double xbegin=-500., xend=500.;
    double dx=1;
    params.get("npoints oned", ndim, ndim);
    params.get("xbegin oned", xbegin, xbegin );
    params.get("xend oned", xend, xend );
    
    //..changes to parameters
    changeParameters(params, react);
    params.get("number of timesteps", nsteps, nsteps);
    params.get("dt",                  dt,     dt);

    printf("--parameters--\n");
    react.outputParameters();    //..output parameters
    //double cascale=1., c0=1.;
    //params.get("calcium rescale", cascale, cascale);
    //react.invariantCalciumScaling( cascale,c0,c );
    //printf("--parameters after rescaling\n");
    //react.outputParameters();

    //..setup grid & initial data
    double *xcoord, *cvec, *hvec, *pvec;
    xcoord = new double[ ndim+1 ];
    cvec =   new double[ ndim+1 ];
    hvec =   new double[ ndim+1 ];
    pvec =   new double[ ndim+1 ];

    cFlux =   new double[ ndim+1 ];
    hFlux =   new double[ ndim+1 ];
    pFlux =   new double[ ndim+1 ];

    cFlux_prev =   new double[ ndim+1 ];
    hFlux_prev =   new double[ ndim+1 ];
    pFlux_prev =   new double[ ndim+1 ];
 
    dx = (xend - xbegin)/(1.* ndim);
    for(int i=0; i<= ndim; ++i) {
      xcoord[i] = xbegin + dx*double(i);
      react.evaluateInitialData( ip3dist, 
	       cvec[i], hvec[i], pvec[i], xcoord[i]);
    }

    //..output file
    std::string outFileName;
    int saveEvery=1;
    params.get("pde output filename", 
	       outFileName, "ode.dat");
    params.get("save every", saveEvery, saveEvery);
    const char *cFileName = outFileName.c_str();
    FILE *fOutput = fopen( cFileName, "w");
    
    
    //printf("\n");
    //printf("-- initial data--");
    //printf("c= %8.3e, h= %8.3e, p= %8.3e\n", c, h, p);
    
    //
    //.. prev step for AB2 = first step is Euler
    //
    //react.computeFlux( 0., c, p, h,
    //cFlux_prev,pFlux_prev,hFlux_prev);
    
    double diff= react.diffCalcium * react.beta;
    for (int kstep=0; kstep <nsteps; ++ kstep) {
      tcomp = kstep*dt;
      if (kstep % saveEvery == 0 ) {
	for(int i=0; i<=ndim; ++i ) {
	  fprintf(fOutput, " %16.8f %16.8f %16.8f\n", 
		  xcoord[i],cvec[i], hvec[i]);
	  //fprintf(fOutput,
	  //      " %16.8f  %16.8f  %16.8f %16.8f %16.8f\n", 
	  //      tcomp,c,h, cFlux, hFlux );
	}
	fprintf(fOutput," \n\n");
      }
      double dq= (diff/(dx*dx));
      for(int i=1; i< ndim; ++i ) {
	react.computeFlux( tcomp, cvec[i],pvec[i],hvec[i],
			   cFlux[i], pFlux[i], hFlux[i] );
	cFlux[i] += dq*( cvec[i+1]-2.*cvec[i]+cvec[i-1] ); 
      }; //end for
      react.computeFlux( tcomp, cvec[0],pvec[0],hvec[0],
			 cFlux[0], pFlux[0], hFlux[0] );
      react.computeFlux( tcomp, cvec[ndim],pvec[ndim],hvec[ndim],
			 cFlux[ndim], pFlux[ndim], hFlux[ndim] );
      cFlux[0]    += dq*( 2*cvec[1] -2*cvec[0] );
      cFlux[ndim] += dq*( 2*cvec[ndim-1] -2*cvec[ndim] );

      cmax=0.; pmax=0.; hmax=0.;
      for(int i=0; i<=ndim; ++i ) {
	cFlux_prev[i] = cFlux[i];  //FORWARD EULER
	hFlux_prev[i] = hFlux[i];

	cvec[i] +=  .5*dt*( 3.*cFlux[i] - cFlux_prev[i]);
	hvec[i] +=  .5*dt*( 3.*hFlux[i] - hFlux_prev[i]);
	//p +=  .5*( 3.*pFlux - pFlux_prev);
	
	if(cvec[i]>cmax) cmax=cvec[i];
	if(hvec[i]>hmax) hmax=hvec[i];
	cFlux_prev[i] = cFlux[i];
	hFlux_prev[i] = hFlux[i];
      }
      printf("%8.3f, cmax= %8.3f, hmax= %8.3f\n",
	     tcomp, cmax, hmax);
    }

    fclose(fOutput);
    printf("-- final   data--");
    printf("c= %8.3e, h= %8.3e, p= %8.3e\n", c, h, p);
  } catch ( ... ) {
    std::cout << "Exiting, error caught.\n";
  }
}






