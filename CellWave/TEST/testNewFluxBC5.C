//
//  test interpolate -- 3D version that works, for implementing flux bc: du/dn= F[u1- u2]
//   
//   ** multicomponent tests!!  cellwave needs multiple components, 
//                                but we've only done 1 component before.
//
#include <iostream>
#include <string>

#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "OGTrigFunction.h"
#include "OGPolyFunction.h"

#include "FluxBC.h"
#include "getDiffusionDT.h"

#include "Display.h"

namespace {

  int localDebug=1; // 1+2+4+8+16;

  enum forcingTypes {
    noForcing = 0,
    polyForcing,
    trigForcing
  }; 

}



int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  printf(" ---------------------------------------------------------------------------- \n");
  printf("   try solving u_t = u_xx + u_yy with flux jump bcs.\n");
  printf(" ---------------------------------------------------------------------------- \n");

  //aString nameOfOGFile, nameOfShowFile, nameOfFluxBCType;
  //std::string nameOfOGFile="Grids/matchingAnnuli20um.hdf", 
  std::string nameOfOGFile="Grids/without-cell33-hirose2d.hdf",
    nameOfShowFile="OUT/nflux5.show", 
    nameOfFluxBCType ="interpolateFlux";
  int initialGrid=0;
  const int physicalBoundaryID=1, fluxBoundaryID=2; // for Overture grids

  if (argc>1) {
    nameOfOGFile       =  argv[1];
  }
  else {
    std::cerr << "usage: "<<argv[0]
	      << " <grid.hdf> <output.show>  <initial source grid #> <fluxBCType> "
	      << endl;
    std::cout << ".. using defaults \n";
  }

  if (argc>2) {
    nameOfShowFile     =  argv[2];
  }
  if (argc>3) {
    initialGrid = atoi( argv[3] );
    printf("[initial grid # = %d]\n", initialGrid );
  }

  if (argc>4) {
    printf("[setting flux bc from args]\n"); fflush(0);
    nameOfFluxBCType   =  argv[4];
  }

  // parameter values
  real t=0, dt=.0005;                                    // initialize time and time step
  real Tmax=100;

  //int saveEvery        =100,  logEvery=20;
  int saveEvery        =100,  logEvery=100;
  bool alwaysEuler = false;

  //real a=1., b=1., viscosity=.1;   // initialize parameters
  real a=1., b=1., viscosity=300;    // initialize parameters
  //real reactionFactor=10.;
  real reactionFactor= 10.;

  real fluxCoeff=1.;
  FluxBC fluxBC;
  FluxBC::FluxBCType bcType= FluxBC::noFluxBC;
  ::forcingTypes forcingOption= ::noForcing;          //::trigForcing;

  //forcingOption = ::trigForcing;

  fluxBC.setDebug(0); //monster amounts of dbg info
  //fluxBC.setDebug(31); //monster amounts of dbg info

  fluxBC.setFluxCoefficient( fluxCoeff );
  bcType = fluxBC.convertToFluxBCType( nameOfFluxBCType );
  
  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile.c_str() );
  cg.update();
  cg.update( CompositeGrid::THEvertex | CompositeGrid::THEvertexBoundaryNormal );
  CompositeGridOperators operators(cg);  // operators for a CompositeGrid
  Interpolant interpolant(cg);           // Make an interpolant

  fluxBC.updateToMatchGrid( cg, fluxBoundaryID );
  fluxBC.setInterpolant( interpolant );
  fluxBC.setOperators( operators );

  const double cfl=0.5; const int nSpecies=1;
  realArray  viscArray(1); viscArray=viscosity;
  real dtDiffusion = getDiffusionDT(cfl, nSpecies,   viscArray, cg );

  if (dtDiffusion<dt) dt= dtDiffusion;
  printf("--DT limits: from diffusion dt=%f, set to dt=%f\n", dtDiffusion, dt);

  //int numberOfTimeSteps=11;
  //int numberOfTimeSteps=10000;
  int numberOfTimeSteps = int(Tmax/dt);

  Ogshow show( nameOfShowFile.c_str() );                          // create a show file
  show.saveGeneralComment("test internal flux/jump bc");  // save a general comment in the show file
  show.setFlushFrequency(1);                              // flush file every N frames

  // create the grid function for the solution
  Range all;

  enum{ u0Component=0, u1Component=1, u2Component=2,      uExactComponent=3, 
	  u0ErrorComponent=4, 	  u1ErrorComponent=5, 	  u2ErrorComponent=6, 
	  u0RHSComponent=7,	  u1RHSComponent=8,	  u2RHSComponent=9,
	  uForcingComponent=10, 
	  u0RHSErrorComponent=11, u1RHSErrorComponent=12, u2RHSErrorComponent=13, 
	  numberOfOutputFields=14};
  const int numberOfFieldComponents=uExactComponent;
  realCompositeGridFunction q(cg,all,all,all,numberOfOutputFields);
  realCompositeGridFunction rhsp(cg,all,all,all,numberOfFieldComponents);
  realCompositeGridFunction unext(cg,all,all,all,numberOfFieldComponents);
  q.setName("u");
  q.setName("u0", u0Component);
  q.setName("u1", u1Component);
  q.setName("u2", u2Component);
  q.setName("exact", uExactComponent);
  q.setName("pointwise error 0", u0ErrorComponent);
  q.setName("pointwise error 1", u1ErrorComponent);
  q.setName("pointwise error 2", u2ErrorComponent);
  q.setName("rhs0", u0RHSComponent );
  q.setName("rhs1", u1RHSComponent );
  q.setName("rhs2", u2RHSComponent );
  q.setName("forcing term", uForcingComponent);
  q.setName("rhs0 error", u0RHSErrorComponent);
  q.setName("rhs1 error", u1RHSErrorComponent);
  q.setName("rhs2 error", u2RHSErrorComponent);
  q.setOperators( operators );
  q = 0.;

  realCompositeGridFunction u(cg,all,all,all,numberOfFieldComponents);
  u.setName("u");                                 // name the grid function
  u.setName("u0", u0Component);                   // name the grid function
  u.setName("u1", u1Component);                   // name the grid function
  u.setName("u2", u2Component);                   // name the grid function

  // create the exact solution (TZ solution)
  OGFunction *pExactSolution=NULL;
  if ( forcingOption == ::polyForcing ) {
    int degreeOfSpacePolynomial  = 2;
    int degreeOfTimePolynomial   = 1;
    pExactSolution = new OGPolyFunction( degreeOfSpacePolynomial, 
					 cg.numberOfDimensions(),
					 degreeOfTimePolynomial);
  }
  else if ( forcingOption == ::trigForcing ) {
    //real fx= 1., fy=1., fz = 1., ft=1.; // note: fz is not used in 2D
    //real fx= 2., fy=2., fz = 2., ft=2.; // note: fz is not used in 2D
    real fx= 1/20., fy=1/20., fz = 1/20., ft=1/20.; // note: fz is not used in 2D
    // defines cos(fx*pi*x)*cos(fy*pi*y)*cos(fz*pi*z)*cos(ft*pi*t)
    pExactSolution = new OGTrigFunction( fx,fy,fz,ft);
  }
  else if (forcingOption != ::noForcing ) {
    std::cerr << "Unknown forcing option ="<< forcingOption
	      << std::endl;
    forcingOption = ::noForcing;
  }
  //OGFunction & exact = *pExactSolution; // make a reference for readibility

  if ( pExactSolution != NULL ) {
    fluxBC.useTwilightZoneFlow( *pExactSolution );
  }
  fluxBC.setupInterpolation(numberOfOutputFields);

  u.setOperators(operators);                                 
  unext.setOperators(operators);
  // operators.setOrderOfAccuracy(4);                          // for fourth order

  // tell the operators to use TZ solution for BCs, if the forcing option >1 ?
  if ( forcingOption >0 ) {
    operators.setTwilightZoneFlow ( true );
    operators.setTwilightZoneFlowFunction( *pExactSolution );
  }

  // auxiliary data
  Index I1,I2,I3, Ib1, Ib2, Ib3;

  // set initial data
  if ( forcingOption >0 ) {
    OGFunction & exact = *pExactSolution; // make a reference for readibility
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.dimension(), I1,I2,I3);
      for(int j=0; j<numberOfFieldComponents; ++j ) {
	u[ig](I1,I2,I3, j+u0Component) = exact( mg, I1,I2,I3, 0, 0.) + 3;
	q[ig](I1,I2,I3, j+u0Component) = u[ig](I1,I2,I3,j+u0Component);
	q[ig](I1,I2,I3, j+u0ErrorComponent) = 0.;
      }

      q[ig](I1,I2,I3, uExactComponent) = exact( mg, I1,I2,I3, 0, 0.) +3;

    } //end for ig
  } else {
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.dimension(), I1,I2,I3);
      for(int j=0; j<numberOfFieldComponents; ++j ) {
	u[ig](I1,I2,I3, j+u0Component) =  1. *(ig==initialGrid);//double(ig);//0.;
	q[ig](I1,I2,I3, j+u0Component) = u[ig](I1,I2,I3, j+u0Component);
	q[ig](I1,I2,I3, j+u0RHSComponent) = -9.;
	q[ig](I1,I2,I3, j+u0ErrorComponent) = 0.;
      }
      q[ig](I1,I2,I3, uExactComponent) = 0.;
    }
  }

  if(localDebug&16) {
    Display disp;
    printf("*main*** show q\n");
    disp.display(q[0], "q[0]");
    disp.display(q[1], "q[1]");
  }
    
  char buffer[80];                       // buffer for sprintf

  printf("**********Time stepper: begin, dt=%f, Tmax=%f, nsteps=%d***\n",dt,Tmax,numberOfTimeSteps);

  for( int istep=0; istep<numberOfTimeSteps; istep++ )                    // take some time steps
  {
    double timerInterior=0.;
    double timerBCs=     0.;
    double timerInterpCode=0.;
    double timer00=getCPU();

    bool   saveThisStep = istep % saveEvery ==0;
    bool   logThisStep  = istep % logEvery  ==0;

    //u+=dt*( -a*u.x() - b*u.y() + viscosity*(u.xx() + u.yy())); // take a time step with Euler's method
    //u+=dt*viscosity*(u.xx() + u.yy()); // take a time step with Euler's method
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.indexRange(), I1,I2,I3);

      //..HEAT EQ
      //q[ig](I1,I2,I3, uRHSComponent)     =   viscosity*(u[ig].laplacian()(I1,I2,I3));
      for(int j=0; j<numberOfFieldComponents; ++j ) { 
	q[ig](I1,I2,I3, j+u0RHSComponent)   =   viscosity*(q[ig].laplacian()(I1,I2,I3,j+u0Component));
	q[ig](I1,I2,I3, j+u0RHSComponent)  +=  reactionFactor*q[ig](I1,I2,I3,j+u0Component)*(1- q[ig](I1,I2,I3,j+u0Component));

	if (forcingOption > 0) {
	  OGFunction & exact = *pExactSolution; // make a reference for readibility
	  q[ig](I1,I2,I3, j+u0RHSErrorComponent) = abs( q[ig](I1,I2,I3, j+u0RHSComponent) 
						     -  viscosity* exact.laplacian(mg, I1,I2,I3,0,t));
	  q[ig](I1,I2,I3, uForcingComponent) =  ( exact.t( mg, I1,I2,I3, 0, t) 
						  -viscosity*exact.laplacian(mg, I1,I2,I3,0,t));
	  q[ig](I1,I2,I3, j+u0RHSComponent)     -=  reactionFactor*exact(mg,I1,I2,I3,0,t)*(1- exact(mg,I1,I2,I3,0,t));
	  q[ig](I1,I2,I3, j+u0RHSComponent) += q[ig](I1,I2,I3, uForcingComponent);
	}

	//u[ig](I1,I2,I3) += dt* q[ig](I1,I2,I3, uRHSComponent); // EULER update
	if( alwaysEuler || istep==0) { 
	  if( (ig==0) && logThisStep ) printf("..euler..");
	  rhsp[ig](I1,I2,I3,j+u0Component) = q[ig](I1,I2,I3,j+u0RHSComponent); 
	}
	else {
	  if (ig==0 && logThisStep ) printf("..ab2..");
	}
	unext[ig](I1,I2,I3,j+u0Component) = u[ig](I1,I2,I3,j+u0Component) 
	        + dt*( (3./2.)*q[ig](I1,I2,I3, j+u0RHSComponent) 
	   	   -(1./2.)  *rhsp[ig](I1,I2,I3, j+u0Component ) );
	rhsp[ig](I1,I2,I3, j+u0Component) = q[ig](I1,I2,I3, j+u0RHSComponent );
      } //end for j
    } //end for ig
    double timer01=getCPU(); timerInterior=timer01-timer00;
    t+=dt;

    if( saveThisStep )  // save solution every 'saveEvery' steps
    {
      show.startFrame();                                         // start a new frame
      show.saveComment(0,sPrintF(buffer,"Here is solution %i",istep));   // comment 0 (shown on plot)
      show.saveComment(1,sPrintF(buffer,"  t=%f ",t));               // comment 1 (shown on plot)
      show.saveSolution( q );                                        // save the current grid function
    }
    
    unext.interpolate();                         // interpolate
    unext.periodicUpdate();
    //.. extrapolate ghostlines so flux/jump bcs ok later, then copy to u
    //.... unext remains a copy of u with ghostline data extrapolated
    unext.applyBoundaryCondition(0,BCTypes::extrapolate,BCTypes::allBoundaries); 
    for(int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      getIndex( cg[ig].dimension(), I1,I2,I3);
      for(int j=0; j<numberOfFieldComponents; ++j ) {
	u[ig](I1,I2,I3, j+u0Component ) = unext[ig]( I1,I2,I3, j+u0Component );
      }
    }

    Range allComponents(0,numberOfFieldComponents-1);
    u.applyBoundaryCondition(allComponents,BCTypes::extrapolate, physicalBoundaryID,0.);    //extrap. ghostlines in u
    u.applyBoundaryCondition(allComponents,BCTypes::extrapolate,     fluxBoundaryID,0.);    
    if (forcingOption == ::noForcing ) { // set all physical boundaries to No-Flux
      u.applyBoundaryCondition(allComponents,BCTypes::neumann,   physicalBoundaryID,0.);
      u.applyBoundaryCondition(allComponents,BCTypes::neumann,   fluxBoundaryID,0.);    

    }
    else if (forcingOption == ::trigForcing ) {
      u.applyBoundaryCondition(allComponents,BCTypes::neumann,   physicalBoundaryID,0.);
      u.applyBoundaryCondition(allComponents,BCTypes::neumann,   fluxBoundaryID,0.);    

    }

    //u.periodicUpdate();
    //u.finishBoundaryConditions();

    double timer02=getCPU();  timerBCs=timer02-timer01;

    //..copy solution to q:
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.dimension(), I1,I2,I3);
      for(int j=0; j<numberOfFieldComponents; ++j ) {	
	//u[ig](I1,I2,I3) =  1. *(ig==initialGrid);//double(ig);//0.;
	q[ig](I1,I2,I3, j+u0Component) = u[ig](I1,I2,I3, j+u0Component);
	//q[ig](I1,I2,I3, uExactComponent) = 0.;
	//q[ig](I1,I2,I3, uErrorComponent) = 0.;
      }
    }
    q.periodicUpdate();
    q.finishBoundaryConditions();
    for(int j=0; j<numberOfFieldComponents; ++j ) {	
      fluxBC.applyBoundaryCondition( q, j);
    }
    q.finishBoundaryConditions();

    if(localDebug&16) {
      Display disp;
      printf("*main*** show q\n");
      disp.display(q[0], "q[0]");
      disp.display(q[1], "q[1]");
    }

    double totalTime=timerBCs+timerInterior;
    if( logThisStep ) {
      printf(" step %5i time= %8.3e, interp %f\% (%f sec), interior %f\% (%f sec) ...\n", 
	     istep,t,
	     100.*timerBCs/totalTime, timerBCs, 
	     100.*timerInterior/totalTime, timerInterior);
      fflush(0);
    }
  }
  Overture::finish();          
  return 0;
}

