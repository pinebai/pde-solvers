#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "OGTrigFunction.h"
#include "OGPolyFunction.h"
#include "PlotStuff.h"
#include "interpolatePoints.h"

#include "CellWave.h"

#include <math.h>
#include "GenericReactionMacros.h"

//
// ..Gaussian Kernel routines
//
struct GaussianKernelData {
  int    numberOfDimensions;
  double timeOffset;
  double viscosity;
  double x0, y0, z0;
  double totalMass;

};

GaussianKernelData 
setGaussianKernelData( double timeOffset, 
		       double viscosity, 
		       double x0, 
		       double y0, 
		       double z0,
		       double totalMass,
		       int numberOfDimensions)
{
  GaussianKernelData kernelData;
  kernelData.timeOffset = timeOffset;
  kernelData.viscosity  = viscosity;
  kernelData.x0         = x0;
  kernelData.y0         = y0;
  kernelData.z0         = z0;
  kernelData.totalMass  = totalMass;
  kernelData.numberOfDimensions = numberOfDimensions;

  return( kernelData );
}

void
printGaussianKernelData( int output, const GaussianKernelData &kernelData )
{
  using CellWave::DPrintf;

  DPrintf(output,"GaussianKernelData:\n");
  DPrintf(output,"   timeOffset = %f\n", kernelData.timeOffset );
  DPrintf(output,"   viscosity  = %f\n", kernelData.viscosity );
  DPrintf(output,"   x0         = %f\n", kernelData.x0 );
  DPrintf(output,"   y0         = %f\n",  kernelData.y0 );     
  DPrintf(output,"   z0         = %f\n",  kernelData.z0 );     
  DPrintf(output,"   mass       = %f\n",  kernelData.totalMass );     
  DPrintf(output,"   numberOfDimensions= %d\n",  kernelData.numberOfDimensions );
}

inline double 
evaluateGaussianKernel( const GaussianKernelData &kernelData, double tcomp, double xc, double yc, double zc )
{
  const double t=tcomp + kernelData.timeOffset;
  const double x=xc    - kernelData.x0;
  const double y=yc    - kernelData.y0;
  const double z=zc    - kernelData.z0;
  const double mu   = kernelData.viscosity;
  const int    ndim = kernelData.numberOfDimensions;  
  const double totalMass = kernelData.totalMass;
wwwwww  const double r2= x*x + y*y + z*z;
  const double pi= 4.*atan( 1.);

  double q;
  
  q = totalMass*exp( - r2/(4. * mu * t ))*pow( 4.* mu * pi * t, -ndim/2.);

  //  printf(" Gaussian=%8.3f, t=%8.3f, mu=%8.3f, r2=%8.3f, ndim/2.=%8.3f\n",
  //           q, t,mu,r2,ndim/2.);
  return( q );
}

void 
heatKernelFortranArray(const GaussianKernelData &kernelData,
		       const double &tcomp,
		       const int&ncomp,
		       const int &nd1a, const int &nd1b,
		       const int &nd2a, const int &nd2b,
		       const int &nd3a, const int &nd3b,
		       const int &n1a,  const int &n1b,
		       const int &n2a,  const int &n2b,
		       const int &n3a,  const int &n3b,
		       const double &xfirst,
		       double &qfirst )
{
  const double *xArray = &xfirst;
  double       *qArray = &qfirst;
  const int nd=kernelData.numberOfDimensions;

  if (nd==1 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int qc=0;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      //const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      //const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      const double y=0.,  z=0.;

      double &q      = qArray[ GRIDINDEX( i,j,k, qc) ];
      
      q = evaluateGaussianKernel( kernelData, tcomp, x,y,z);
    }
  }
  else if (nd==2 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int qc=0;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      //const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      const double z=0.;

      double &q      = qArray[ GRIDINDEX( i,j,k, qc) ];
      
      q = evaluateGaussianKernel( kernelData, tcomp, x,y,z);
    }
  }
  else   if (nd==3 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int qc=0;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      
      double &q      = qArray[ GRIDINDEX( i,j,k, qc) ];
      
      q = evaluateGaussianKernel( kernelData, tcomp, x,y,z);
    }
  }
  else printf("--unknown dimension = %d -- error\n", nd);
}


void
heatKernelGridFunction( GaussianKernelData &kernelData, 
			double tcomp, 
			MappedGrid        &mg, 
			realArray               &heatKernel )
{
  realArray & x = mg.vertex();  // array of vertices
  const IntegerArray & d    =  mg.dimension();
  const IntegerArray & gir  =  mg.gridIndexRange();
  const int nd=mg.numberOfDimensions();
  const int ncomp=1; //only set first component
  
  kernelData.numberOfDimensions = nd;
  // call a fortran function to compute du/dt
  // (This function does not currently solve the convection diffusion equation)
  //mySolver( t,dt,a,b,nu,nd, d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2),
  //           gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
  //          *x.getDataPointer(),*ug.getDataPointer(), *dudtg.getDataPointer() );
  heatKernelFortranArray( kernelData, tcomp,
			  ncomp,
			  d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2),
			  gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
			  *x.getDataPointer(), *heatKernel.getDataPointer() );
			  			  
}

int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture
  CellWave::start(argc,argv);  // initialize CellWave

  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;
  const int LOGPRINT=CellWave::LogPrint;

  try {

  DPrintf(BPRINT," ---------------------------------------------------------------------------- \n");
  DPrintf(BPRINT,"Solve the heat eq  u.t = viscosity*( u.xx + u.yy ) on an Overlapping grid \n");
  DPrintf(BPRINT,"Save results in a show file, use plotStuff to view this file                  \n");
  DPrintf(BPRINT,"  ... added TZ to test the accuracy ...                                       \n");
  DPrintf(BPRINT," ---------------------------------------------------------------------------- \n");

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
    throw "error";
  }

  //
  //..Read parameter values
  //
  CellWave::ParameterReader &param = *pParameters;

  std::string paramFileType="";
  std::string acceptableParamFileType="test diffusion";
  param.get("parameter file type", paramFileType, "");

  if( paramFileType != acceptableParamFileType ) {
    DPrintf(BPRINT,"ERROR: parameter file type '%s' unknown, exiting.");
    throw "error";    
  }

  std::string tempstr="";
  aString nameOfOGFile, nameOfShowFile;
  double timeStepSize=0.01;
  double viscosity=0.01;
  int numberOfTimeSteps=100;
  int saveEveryNthFrame=1;

  param.get( "name of grid file",    tempstr,  "" );        nameOfOGFile   = tempstr.c_str();
  param.get( "name of show file",    tempstr,"");           nameOfShowFile = tempstr.c_str();
  param.get( "maximum timestep",     timeStepSize,  0.1);
  param.get( "number of timesteps",  numberOfTimeSteps, 1);
  param.get( "save frequency",	     saveEveryNthFrame, 10);

  DPrintf(PRINT,".. Grid file=%s, Show file=%s\n", 
	    nameOfOGFile.c_str(), nameOfShowFile.c_str());
  DPrintf(PRINT,"..     dt=%8.4e, num. steps=%d, save every %d frame, Tmax=%8.4e\n",
	   timeStepSize,      numberOfTimeSteps, 
	   saveEveryNthFrame, timeStepSize*numberOfTimeSteps);

  param.get( "timestep size", timeStepSize, timeStepSize );
  param.get( "viscosity", viscosity, viscosity );
  param.get( "number of timesteps", numberOfTimeSteps, numberOfTimeSteps);
  param.get( "saveEveryNthFrame", saveEveryNthFrame, saveEveryNthFrame );

  bool openGraphicsWindow=TRUE;
  std::string isInteractive;
  param.get("interactive", isInteractive, "yes");
  openGraphicsWindow =  (isInteractive == "yes");
  DPrintf(PRINT,".. open graphics window (interactive) =");
  if (openGraphicsWindow) {
    DPrintf(PRINT,"yes\n");
  }
  else {
    DPrintf(PRINT,"no\n");
  }

  bool outputText=TRUE;
  std::string isOutputting;
  param.get("output", isOutputting, "yes");
  outputText =  (isOutputting == "yes");
  DPrintf(PRINT,".. isOutputting                       =");
  if ( outputText ) {
    DPrintf(PRINT,"yes\n");
  }
  else {
    DPrintf(PRINT,"no\n");
  }


  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update(MappedGrid::THEvertex | MappedGrid::THEmask);      // build vertices and mask

  //..Read & setup Gaussian Kernel data
  GaussianKernelData kernelData;
  { 
    double timeOffset;
    double x0, y0, z0;
    double totalMass;
    param.get( "time offset", timeOffset, 0.01);
    param.get( "x offset",   x0,  0.);
    param.get( "y offset",   y0,  0.);
    param.get( "z offset",   z0,  0.);
    param.get( "total mass", totalMass, 1.);
        
    kernelData=setGaussianKernelData( timeOffset,
				      viscosity,
				      x0, y0, z0,
				      totalMass,
				      cg.numberOfDimensions());
  }
  printGaussianKernelData( PRINT, kernelData);

  //..setup parameter values
  real tcomp=0., dt= timeStepSize;

  Interpolant interpolant(cg);                                 // Make an interpolant

  Ogshow show( nameOfShowFile );                               // create a show file
  show.saveGeneralComment("Diffusion Equation Point Spread");    // save a general comment in the show file
  show.setFlushFrequency(10);                                  // flush file every 10 frames

  // create the grid function for the solution
  Range all;
  const int numberOfFieldComponents=1;

  enum{ uComponent=0, uExactComponent=1, uErrorComponent=2, 
	  uRHSComponent=3, uForcingComponent=4, uRHSErrorComponent=5, numberOfOutputFields};
  realCompositeGridFunction q(cg,all,all,all,numberOfOutputFields);
  q.setName("u");
  q.setName("u", uComponent);
  q.setName("exact", uExactComponent);
  q.setName("pointwise error", uErrorComponent);
  q.setName("rhs", uRHSComponent );
  q.setName("forcing term", uForcingComponent);
  q.setName("rhs error", uRHSErrorComponent);

  //..create operators  & use those with the grid function
  CompositeGridOperators operators(cg);                        // operators for a CompositeGrid
  q.setOperators(operators);                                 
  // operators.setOrderOfAccuracy(4);                          // for fourth order

  const int nExactFields = 1;
  realCompositeGridFunction heatKernel( cg, all,all,all,nExactFields);
  Index I1,I2,I3, Id1,Id2,Id3;
  for( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
    MappedGrid &mg = cg[ig];
    getIndex( mg.indexRange(), I1,I2,I3);
    getIndex( mg.dimension(),  Id1,Id2,Id3);
    
    //..get heat kernel =exact solution
    realArray  &heatArray  = heatKernel[ig];
    heatKernelGridFunction( kernelData, tcomp, mg, heatArray) ;
    
    q[ig](Id1,Id2,Id3, uExactComponent) = heatArray(Id1,Id2,Id3);  

    //..set initial data
    q[ig](Id1,Id2,Id3, uComponent)      = heatArray(Id1,Id2,Id3); 
  }


  //..start timestepping (and setup graphics)
  PlotStuff ps(openGraphicsWindow,"testDiffusion");  // create a PlotStuff object
  PlotStuffParameters psp;                         // This object is used to change plotting parameters

  if( openGraphicsWindow) {
    PlotIt::contour(ps,q,psp);
    ps.redraw(TRUE);
    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,FALSE);
  }

  Index Ib1,Ib2,Ib3;
  double maxerr=0.;
  int kLastStep=0; // used to output last frame after tstep loop
  for( int ktime=0; ktime<numberOfTimeSteps; ++ktime ) {

    if( ktime % saveEveryNthFrame==0 ) {
      if(outputText) {
	DPrintf(PRINT,".. t=%8.3f,  max| error | = %16.8e\n", tcomp, maxerr );
      }

      DPrintf(LOGPRINT,".. saving frame at step %d, t=%e\n",ktime,tcomp);
      char buffer[80];
      show.startFrame();                                                // start a new frame
      show.saveComment(0,sPrintF(buffer,"Diffusion test: Solution %d ",ktime));  // comment 0 (shown on plot)
      show.saveComment(1,sPrintF(buffer,"  t = %e ",tcomp));                  // comment 1 (shown on plot)
      show.saveSolution( q );                                         // save the current grid function
    }

    tcomp += dt;
    maxerr=0.;
    for( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid &mg = cg[ig];
      getIndex( mg.indexRange(), I1,I2,I3);
      getIndex( mg.dimension(),  Id1,Id2,Id3);

      //..get heat kernel =exact solution
      realArray  &heatArray  = heatKernel[ig];
      heatKernelGridFunction( kernelData, tcomp, mg, heatArray) ;

      q[ig](Id1,Id2,Id3, uExactComponent) = heatArray(Id1,Id2,Id3);

      //..solve heat equation
      q[ig](I1,I2,I3, uRHSComponent) = viscosity*q[ig].laplacian()(I1,I2,I3,
								   uComponent);
      q[ig](I1,I2,I3, uComponent) += dt* q[ig](I1,I2,I3, uRHSComponent);

      //....boundary conditions:
      for( int side=Start; side<=End; side++ )
	for( int axis=axis1; axis<cg.numberOfDimensions(); axis++ ) {
	  if( mg.boundaryCondition()(side,axis) > 0 ) {
	    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	    q[ig](Ib1,Ib2,Ib3,uComponent)
	      =q[ig](Ib1,Ib2,Ib3, uExactComponent);
      }
    }

      //..compute error
      q[ig](I1,I2,I3, uErrorComponent) = 
	q[ig](I1,I2,I3, uComponent) - q[ig](I1,I2,I3, uExactComponent);
      
      double err=max( abs( q[ig](I1,I2,I3,uErrorComponent)));
      maxerr = max( maxerr, err );
      //DPrintf(PRINT,".. t=%8.3f,  max| error | = %16.8e\n", tcomp, maxerr );
    }
    kLastStep = ktime;
  }

  double maxSol=0., maxExact=0.;
  for( int ig=0; ig<cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid &mg = cg[ig];
      getIndex( mg.indexRange(), I1,I2,I3);
      maxSol  =max(maxSol,    max( abs( q[ig](I1,I2,I3,uComponent))));
      maxExact=max(maxExact,  max( abs( q[ig](I1,I2,I3,uExactComponent))));
  }

  //..print final info
  DPrintf(PRINT,"\nFINAL t=%8.3f:\n", tcomp);
  DPrintf(PRINT,"   MAX |ERROR|      = %16.8e\n\n", maxerr );
  DPrintf(PRINT,"   ..where max |Solution| = %16.8e\n", maxSol);
  DPrintf(PRINT,"         & max |Exact|    = %16.8e\n", maxExact);

  //..output final frame
  show.startFrame();                                                // start a new frame
  char buffer[80];
  show.saveComment(0,sPrintF(buffer,"Diffusion test: Solution %d ",kLastStep+1));  // comment 0 (shown on plot)
  show.saveComment(1,sPrintF(buffer,"  t = %e ",tcomp));                  // comment 1 (shown on plot)
  show.saveSolution( q );                                         // save the current grid function
  
  if (openGraphicsWindow) {
    PlotIt::contour(ps,q,psp);
    ps.redraw(TRUE);
    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,FALSE);
    //psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,TRUE);
  }

  } //end try
  catch (...) {
    DPrintf(BPRINT, "..fatal exception caught, finishing & exiting..\n");    
  }
  CellWave::finish();
  Overture::finish();          
  return 0;
}


