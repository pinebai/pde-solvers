//
// testFluxInterpolate -- driver for u_n = C [ u ] bc classes
//

#include "Overture.h"
#include "Ogshow.h"
#include "CompositeGridOperators.h"

#include <string>
#include <iostream>
#include <stdio.h>

#include "ParameterReader.h"
#include "getDiffusionDT.h"
//#include "FluxBC.h"
#include "FastFluxBC.h"

namespace {

  int localDebug=1; // 1+2+4+8;

  enum forcingTypes {
    noForcing = 0,
    polyForcing,
    trigForcing
  }; 

  enum fluxBCType {
    noFluxBC          = 0,
    interpolateFluxBC = 1
  }; 

  class TestFluxData {
  public:
    TestFluxData() :
        physicalBoundaryID(1), fluxBoundaryID(2),
        initialGrid(14),       nameOfFluxBCType("interpolateFlux"),
	nsteps(10), saveEvery(1)
    { };

    ~TestFluxData() {};

    void processCommandLine( int argc, char **argv );
    aString getNameOfOGFile() { return( nameOfOGFile.c_str()); };
    void getNameOfOGFile( aString     &name ) { name = nameOfOGFile.c_str(); };
    void getNameOfOGFile( std::string &name ) { name = nameOfOGFile; };

    aString getNameOfShowFile() { return( nameOfShowFile.c_str()); };
    void getNameOfShowFile( aString     &name ) { name = nameOfShowFile.c_str(); };
    void getNameOfShowFile( std::string &name ) { name = nameOfShowFile; };

    aString getNameOfFluxBCType() { return( nameOfFluxBCType.c_str()); };
    void getNameOfFluxBCType( aString     &name ) { name = nameOfFluxBCType.c_str(); };
    void getNameOfFluxBCType( std::string &name ) { name = nameOfFluxBCType; };

    int  getNumberOfSteps () { return nsteps; };
    int  getSaveEvery() { return saveEvery; };

    void print() { 
      std::cout << "The 'TestFluxData':\n";
      std::cout << "  nameOfOGFile         = "<< nameOfOGFile       << std::endl;
      std::cout << "  nameOfShowFile       = "<< nameOfShowFile     << std::endl;
      std::cout << "  nameOfFluxBCType     = "<< nameOfFluxBCType   << std::endl;
      std::cout << "  physicalBC ID        = "<< physicalBoundaryID << std::endl;
      std::cout << "  fluxBC ID            = "<< fluxBoundaryID     << std::endl;
      std::cout << "  initialGrid number   = "<< initialGrid        << std::endl;
      std::cout << "  nsteps               = "<< nsteps             << std::endl;
      std::cout << std::endl;
    }

  private: //do not use the copy constructor--> hence private
    TestFluxData( TestFluxData &X);
    TestFluxData & operator=( TestFluxData &X);
  public:
    std::string nameOfOGFile, nameOfShowFile, nameOfFluxBCType;
    const int physicalBoundaryID, fluxBoundaryID; // for Overture grids
    int initialGrid;
    int nsteps, saveEvery;
  };

  void 
  TestFluxData::processCommandLine( int argc, char **argv )
  {
    if (argc>1) {
      nameOfOGFile       =  argv[1];
      nameOfShowFile     =  argv[2];
    }
    else {
      std::cerr << "usage: "<<argv[0]
		<< " <grid.hdf> <output.show> <initial source grid #> <# of steps> "
		<< endl;
      exit(1);
    }
    if (argc>3) {
      initialGrid = atoi( argv[3] );
      printf("[initial grid # = %d]\n", initialGrid );
    }

    if (argc>4) {
      nsteps = atoi( argv[4] );
      printf("[number of steps = %d]\n", nsteps );
    }
    if (argc>5) {
      saveEvery = atoi( argv[5] );
      printf("[save every %d steps]\n", saveEvery);
    }
  }
  
  aString makeAString( std::string instring ) {
    return( instring.c_str() );
  }

  void saveFrame( Ogshow &show, realCompositeGridFunction &q, int istep, real tcomp )
  {
    char buffer[80];
    show.startFrame();
    show.saveComment(0,sPrintF(buffer,"Here is solution %i",istep));   // comment 0 (shown on plot)
    show.saveComment(1,sPrintF(buffer,"  t=%f ",tcomp));               // comment 1 (shown on plot)
    show.saveSolution( q );                                            // save the current grid function
  }

};



int
main( int argc, char **argv)
{
  Overture::start( argc, argv );
  printf("---------------------------------------------------------------\n");
  printf("       test flux bc classes \n");
  printf("---------------------------------------------------------------\n");

  ::TestFluxData fluxData;
  fluxData.processCommandLine( argc, argv );
  fluxData.print();

  //--------------------SETUP OVERTURE
  // parameter values
  real t=0, dt=.0005;                                    // initialize time and time step
  //real a=1., b=1., viscosity=.1;                      // initialize parameters
  real a=1., b=1., viscosity=10;                        // initialize parameters
  real fluxCoeff=10.;
  //FluxBC fluxBC;
  FastFluxBC fluxBC;
  FluxBC::FluxBCType bcType= FluxBC::noFluxBC;
  ::forcingTypes forcingOption= ::noForcing;          //::trigForcing;
  
  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg, fluxData.getNameOfOGFile() );
  cg.update();
  cg.update( CompositeGrid::THEvertex | CompositeGrid::THEvertexBoundaryNormal );
  CompositeGridOperators operators(cg);  // operators for a CompositeGrid
  Interpolant interpolant(cg);           // Make an interpolant

  fluxBC.updateToMatchGrid( cg, fluxData.fluxBoundaryID );
  fluxBC.setInterpolant( interpolant );
  fluxBC.setOperators( operators );

  //--------------------SETUP SOLVER
  const double cfl=0.5; const int nSpecies=1;
  realArray  viscArray(1); viscArray=viscosity;
  real dtDiffusion = getDiffusionDT(cfl, nSpecies,   viscArray, cg );

  printf("--DT limits: from diffusion dt=%f, preset dt=%f\n", dtDiffusion, dt);
  dt = cfl * dtDiffusion;


  Ogshow show( fluxData.getNameOfShowFile() );                          // create a show file
  show.saveGeneralComment("test internal flux/jump bc");
  show.setFlushFrequency(5);                            
  

  // create the grid function for the solution
  Range all;
  const int numberOfFieldComponents=1;
  realCompositeGridFunction u(cg,all,all,all,numberOfFieldComponents);
  u.setOperators( operators );
  u.setName("u");                                              // name the grid function

  enum{ uComponent=0, uExactComponent=1, uErrorComponent=2, 
	  uRHSComponent=3, uForcingComponent=4, uRHSErrorComponent=5, numberOfOutputFields};
  realCompositeGridFunction q(cg,all,all,all,numberOfOutputFields);
  q.setOperators( operators );
  realCompositeGridFunction unext(cg,all,all,all,1);
  unext.setOperators( operators );
  q.setName("u");
  q.setName("u", uComponent);
  q.setName("exact", uExactComponent);
  q.setName("pointwise error", uErrorComponent);
  q.setName("rhs", uRHSComponent );
  q.setName("forcing term", uForcingComponent);
  q.setName("rhs error", uRHSErrorComponent);


  // auxiliary data
  Index I1,I2,I3, Ib1, Ib2, Ib3;

  // set initial data
  for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
    MappedGrid & mg = cg[ig];
    getIndex( mg.dimension(), I1,I2,I3);
    u[ig](I1,I2,I3) =  1. *(ig==fluxData.initialGrid);//double(ig);//0.;
    q[ig](I1,I2,I3, uComponent) = u[ig](I1,I2,I3);
    q[ig](I1,I2,I3, uExactComponent) = 0.;
    q[ig](I1,I2,I3, uErrorComponent) = 0.;
  }
 
  //--------------SOLVER TIMESTEPPING LOOP

  char buffer[80];                       // buffer for sprintf
  int saveEvery        = fluxData.getSaveEvery();
  //int numberOfTimeSteps=11;
  int numberOfTimeSteps= fluxData.getNumberOfSteps();

  for( int istep=0; istep<numberOfTimeSteps; istep++ )                    // take some time steps
  {
    double timerInterior=0.;
    double timerBCs=     0.;
    double timerInterpCode=0.;
    double timer00=getCPU();

    if( istep % saveEvery == 0 ) ::saveFrame(show, q, istep, t );

    //u+=dt*( -a*u.x() - b*u.y() + viscosity*(u.xx() + u.yy())); // take a time step with Euler's method
    //u+=dt*viscosity*(u.xx() + u.yy()); // take a time step with Euler's method
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.indexRange(), I1,I2,I3);

      //..HEAT EQ
      q[ig](I1,I2,I3, uRHSComponent)     =   viscosity*(u[ig].laplacian()(I1,I2,I3));
      unext[ig](I1,I2,I3) = u[ig](I1,I2,I3) +dt* q[ig](I1,I2,I3, uRHSComponent); // EULER update

    } //end for ig
    double timer01=getCPU(); timerInterior=timer01-timer00;
    t+=dt;
    
    unext.interpolate();                                           // interpolate
    unext.periodicUpdate();
    //.. extrapolate ghostlines so flux/jump bcs ok later, then copy to u
    //.... unext remains a copy of u with ghostline data extrapolated
    unext.applyBoundaryCondition(0,BCTypes::extrapolate,BCTypes::allBoundaries); 
    for(int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      getIndex( cg[ig].dimension(), I1,I2,I3);
      u[ig](I1,I2,I3) = unext[ig]( I1,I2,I3 );
    }

    const int iSolComponent=0;
    fluxBC.applyBoundaryCondition( unext, u, iSolComponent );

    if (forcingOption == ::noForcing ) { // set all physical boundaries to No-Flux
      u.applyBoundaryCondition(0,BCTypes::neumann, fluxData.physicalBoundaryID,0.);
    }

    u.periodicUpdate();
    u.finishBoundaryConditions();

    double timer02=getCPU();  timerBCs=timer02-timer01;

    //..copy solution to q:
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.dimension(), I1,I2,I3);
      //u[ig](I1,I2,I3) =  1. *(ig==initialGrid);//double(ig);//0.;
      q[ig](I1,I2,I3, uComponent) = u[ig](I1,I2,I3);
      //q[ig](I1,I2,I3, uExactComponent) = 0.;
      //q[ig](I1,I2,I3, uErrorComponent) = 0.;
    }
    double totalTime=timerBCs+timerInterior;
    if( istep % saveEvery == 0 ) {
      printf(" step %5i time= %12.3e, interp %6.2f\% (%10e sec), interior %6.2f\% (%10e sec) ...\n", 
	     istep,t,
	     100.*timerBCs/totalTime, timerBCs, 
	     100.*timerInterior/totalTime, timerInterior);
    }
    fflush(0);
  } //end for istep


  Overture::finish();

}

