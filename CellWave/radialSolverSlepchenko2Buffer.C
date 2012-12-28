#include <iostream>
#include <stdio.h>
#include <math.h>

#include "CellWave.h"
//#include "getDiffusionDT.h"
#include "SolverLiRinzel.h"
#include "ReactionLiRinzelWagner.h"

#include "SolverSlepchenko2Buffer.h"
#include "ReactionSlepchenko2Buffer.h"


using namespace CellWave;
using CellWave::ParameterReader;

bool 
isParameterFileCompatible( ParameterReader &params )
{
  std::string typeString="";
  params.get("parameter file type",typeString);
  bool ok= ( typeString== "Slepchenko2Buffer");
  if( !ok ) {
    std::cout <<" --ERROR: incompatible parameter file, "
	      <<"type="<<typeString << std::endl;
  }
  return( ok );
}

int main(int argc, char **argv)
{
  try {
    CellWave::start(argc,argv);
    const int BPRINT=CellWave::BroadcastPrint;
    const int PRINT=CellWave::PrintOut;
    
    DPrintf(BPRINT," --------------------------------------------------------------------- \n");
    DPrintf(BPRINT,"  .. radialSolverSlepchenko2Buffer ..\n");
    DPrintf(BPRINT,"      simulating Ca waves in a LiRinzel rxn w/ 2 buffers\n");
    DPrintf(BPRINT," --------------------------------------------------------------------- \n");
    if(argc <= 1) {
      std::cerr <<"Usage: radialSolverSlepchenko2Buffer <parameter file name>\n";
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
    param.get( "save frequency",       data.saveEveryNthFrame, 10);
    
    DPrintf(PRINT,".. Grid file=%s, Show file=%s\n", 
	    data.nameOfOGFile.c_str(), data.nameOfShowFile.c_str());
    DPrintf(PRINT,"..     dt=%8.4e, num. steps=%d, save every %d frame, Tmax=%8.4e\n",
	    data.timeStepSize, data.numberOfTimeSteps, 
	    data.saveEveryNthFrame, data.maximumTime());

    if( !isParameterFileCompatible( param ) ) {
      throw "Unknown parameter file type";
    }
    
    //..setup rxns

    ReactionSlepchenko2Buffer chem;
    bool ok=  chem.readParameterFile( param );
    enum { cc=0, pc=1, hc=2,
	   b1c=3, b2c=4,nComponents=5};// hardwired components. must match RxnSlepchenko...

    //..setup the radial 1-dim wave solver
    double x=0., y=0., z=0.;
    //double p,c,h, b1, b2; // pn, cn, hn;
    //double pFlux, cFlux, hFlux, b1Flux, b2Flux;
    //double pFlux_prev, cFlux_prev, hFlux_prev, b1Flux_prev, b2Flux_prev;
    double minval[nComponents], maxval[nComponents];

    int     nsteps = data.numberOfTimeSteps;
    double  dt     = data.timeStepSize;
    double  tcomp  = 0.;
    
    int    ndim=100;
    int    numberOfSpaceDimensions;
    double rInner= 0., rOuter=1000.;
    double dx=1;

    numberOfSpaceDimensions = chem.getNumberOfDimensions();

    param.get("npoints radial", ndim, ndim);
    param.get("r inner radial", rInner, 0.);
    param.get("r outer radial", rOuter, 1000.);

    //..changes to parameters
    param.get("number of timesteps", nsteps, nsteps);
    param.get("dt",                  dt,     dt);
 
    //..setup grid & initial data
    double *rcoord = new double[ ndim+1 ];
    double *rinv   = new double[ ndim+1 ];
    double *cvec =   new double[ ndim+1 ];
    double *hvec =   new double[ ndim+1 ];
    double *pvec =   new double[ ndim+1 ];
    double *b1vec=   new double[ ndim+1 ];
    double *b2vec=   new double[ ndim+1 ];

    double *cFlux =   new double[ ndim+1 ];
    double *hFlux =   new double[ ndim+1 ];
    double *pFlux =   new double[ ndim+1 ];
    double *b1Flux =  new double[ ndim+1 ];
    double *b2Flux =  new double[ ndim+1 ];

    double *cFlux_prev  =   new double[ ndim+1 ];
    double *hFlux_prev  =   new double[ ndim+1 ];
    double *pFlux_prev  =   new double[ ndim+1 ];
    double *b1Flux_prev =   new double[ ndim+1 ];
    double *b2Flux_prev =   new double[ ndim+1 ];
 
    printf("INITIAL DATA:  r, p\n");
    double dr = (rOuter - rInner)/(1.* ndim);
    for(int i=0; i<= ndim; ++i) {
      rcoord[i] = rInner + dr*double(i);
      double x,y,z;
      x=rcoord[i];  y=0.;  z=0.;
      chem.setInitialData( x,y,z, cvec[i], hvec[i], pvec[i], b1vec[i], b2vec[i]);
      printf("  %8.3f   %8.3f \n", rcoord[i], pvec[i]);
    }
    printf("-------------------------------\n");

    rinv[0]=0.;
    for( int i=1; i<=ndim; ++i ) rinv[i] = 1./rcoord[i];

    //..output file
    std::string outFileName;
    int saveEvery=1;
    param.get("pde output filename", 
	       outFileName, "pde.dat");
    param.get("save every", saveEvery, saveEvery);
    const char *cFileName = outFileName.c_str();
    FILE *fOutput = fopen( cFileName, "w");

    //
    //.. prev step for AB2 = first step is Euler
    //
    //react.computeFlux( 0., c, p, h,
    //cFlux_prev,pFlux_prev,hFlux_prev);
    
    //double diff= react.diffCalcium * react.beta;
    double visc[nComponents];
    visc[cc]  = chem.getDiffusionCoefficient( cc );
    visc[pc]  = chem.getDiffusionCoefficient( pc );
    visc[hc]  = chem.getDiffusionCoefficient( hc );
    visc[b1c] = chem.getDiffusionCoefficient( b1c );
    visc[b2c] = chem.getDiffusionCoefficient( b2c );
    double diff[nComponents][ndim];
    
    for (int kstep=0; kstep <nsteps; ++ kstep) {
      tcomp = kstep*dt;
      if (kstep % saveEvery == 0 ) {
	fprintf(fOutput, " 0 %16.8f ", tcomp);
	for(int i=0; i<=ndim; ++i ) fprintf(fOutput, " %16.8f ", rcoord[i]); // 0=rcoord
	fprintf(fOutput, "\n 1 %16.8f ", tcomp);
	for(int i=0; i<=ndim; ++i ) fprintf(fOutput, " %16.8f ", cvec[i]);   // 1=cvec
	fprintf(fOutput, "\n 2 %16.8f ", tcomp);
	for(int i=0; i<=ndim; ++i ) fprintf(fOutput, " %16.8f ", hvec[i]);   // 2=hvec
	fprintf(fOutput, "\n 3 %16.8f ", tcomp);
	for(int i=0; i<=ndim; ++i ) fprintf(fOutput, " %16.8f ", pvec[i]);   // 3=pvec
	fprintf(fOutput, "\n 4 %16.8f ", tcomp);
	for(int i=0; i<=ndim; ++i ) fprintf(fOutput, " %16.8f ", b1vec[i]);   // 4=b1vec
	fprintf(fOutput, "\n 5 %16.8f ", tcomp);
	for(int i=0; i<=ndim; ++i ) fprintf(fOutput, " %16.8f ", b2vec[i]);   // 5=b2vec
	fprintf(fOutput,"\n");
      }
      for(int i=0; i<= ndim; ++i ) {
	chem.computeRHS( cvec[i],    pvec[i],  hvec[i],  b1vec[i],  b2vec[i],
			   cFlux[i], pFlux[i], hFlux[i], b1Flux[i], b2Flux[i] );
      } //end for;

      double drinv2= 1./(dr*dr);
      for(int i=1; i<  ndim; ++i ) {
	double nr=(1./(2.*dr))*(double(numberOfSpaceDimensions)-1.)*rinv[i];

	diff[cc][i]  = drinv2*( cvec[i+1]  -2.*cvec[i]  +cvec[i-1] )   +nr*( cvec[i+1]-cvec[i-1]);
	diff[pc][i]  = drinv2*( pvec[i+1]  -2.*pvec[i]  +pvec[i-1] )   +nr*( pvec[i+1]-pvec[i-1]);
	diff[hc][i]  = drinv2*( hvec[i+1]  -2.*hvec[i]  +hvec[i-1] )   +nr*( hvec[i+1]-hvec[i-1]);
	diff[b1c][i] = drinv2*( b1vec[i+1] -2.*b1vec[i] +b1vec[i-1] )  +nr*( b1vec[i+1]-b1vec[i-1]);
	diff[b2c][i] = drinv2*( b2vec[i+1] -2.*b2vec[i] +b2vec[i-1] )  +nr*( b2vec[i+1]-b2vec[i-1]);
      }; //end for

      diff[cc][0] = drinv2*( 2*cvec[1]-2*cvec[0] );   diff[cc][ndim]  = drinv2*(2*cvec[ndim-1]-2*cvec[ndim]);
      diff[pc][0] = drinv2*( 2*pvec[1]-2*pvec[0] );   diff[pc][ndim]  = drinv2*(2*pvec[ndim-1]-2*pvec[ndim]);
      diff[hc][0] = drinv2*( 2*hvec[1]-2*hvec[0] );   diff[hc][ndim]  = drinv2*(2*hvec[ndim-1]-2*hvec[ndim]);
      diff[b1c][0]= drinv2*( 2*b1vec[1]-2*b1vec[0] ); diff[b1c][ndim] = drinv2*(2*b1vec[ndim-1]-2*b1vec[ndim]);
      diff[b2c][0]= drinv2*( 2*b2vec[1]-2*b2vec[0] ); diff[b2c][ndim] = drinv2*(2*b2vec[ndim-1]-2*b2vec[ndim]);

      for(int i=0; i<= ndim; ++i ) {
	cFlux[i]   +=  visc[cc] *diff[cc] [i];
	pFlux[i]   +=  visc[pc] *diff[pc] [i];
	hFlux[i]   +=  visc[hc] *diff[hc] [i];
	b1Flux[i]  +=  visc[b1c]*diff[b1c][i];
	b2Flux[i]  +=  visc[b2c]*diff[b2c][i];
      }

      for( int j=0; j<nComponents; ++j ) { maxval[j]=0.; };
      for(int i=0; i<=ndim; ++i ) {
	cFlux_prev[i] = cFlux[i];  //FORWARD EULER for first step
        pFlux_prev[i] = pFlux[i];
	hFlux_prev[i] = hFlux[i];
	b1Flux_prev[i] = b1Flux[i];
	b2Flux_prev[i] = b2Flux[i];

	cvec[i]  +=  .5*dt*( 3.*cFlux[i] - cFlux_prev[i]);
	pvec[i]  +=  .5*dt*( 3.*pFlux[i] - pFlux_prev[i]);
	hvec[i]  +=  .5*dt*( 3.*hFlux[i] - hFlux_prev[i]);
	b1vec[i] +=  .5*dt*( 3.*b1Flux[i] - b1Flux_prev[i]);
	b2vec[i] +=  .5*dt*( 3.*b2Flux[i] - b2Flux_prev[i]);
	
	if(cvec[i]>maxval[cc])   maxval[cc] =cvec[i];
	if(pvec[i]>maxval[pc])   maxval[pc] =pvec[i];
	if(hvec[i]>maxval[hc])   maxval[hc] =hvec[i];
	if(b1vec[i]>maxval[b1c]) maxval[b1c]=b1vec[i];
	if(b2vec[i]>maxval[b2c]) maxval[b2c]=b2vec[i];
	cFlux_prev[i] = cFlux[i];
	pFlux_prev[i] = pFlux[i];
	hFlux_prev[i] = hFlux[i];
	b1Flux_prev[i] = b1Flux[i];
	b2Flux_prev[i] = b2Flux[i];
      }
      printf("%8.3f, cmax= %8.3f, pmax= %8.3f, hmax= %8.3f, b1max= %8.3f, b2max= %8.3f\n",
	     tcomp, maxval[cc], maxval[pc], maxval[hc], maxval[b1c], maxval[b2c]);
    }
    fclose(fOutput);
  } catch ( ... ) {
    std::cout << "Exiting, error caught.\n";
  }
}






