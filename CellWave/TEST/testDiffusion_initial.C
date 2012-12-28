#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "OGTrigFunction.h"
#include "OGPolyFunction.h"

namespace {
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
  printf("Solve the heat eq  u.t = viscosity*( u.xx + u.yy ) on an Overlapping grid \n");
  printf("Save results in a show file, use plotStuff to view this file                  \n");
  printf("  ... added TZ to test the accuracy ...                                       \n");
  printf(" ---------------------------------------------------------------------------- \n");

  aString nameOfOGFile, nameOfShowFile, nameOfForcingType;

  if (argc>3) {
    nameOfOGFile       =  argv[1];
    nameOfShowFile     =  argv[2];
    nameOfForcingType  =  argv[3];
  }
  else {
    std::cerr << "usage: "<<argv[0]
	      << " <grid.hdf> <output.show> "
	      << endl;
    exit(1);
  }

  // parameter values
  real t=0, dt=.0001;                                           // initialize time and time step
  real a=1., b=1., viscosity=.1;                               // initialize parameters
  ::forcingTypes forcingOption=::polyForcing;

  if ( nameOfForcingType == "noForcing" ) {
    forcingOption = ::noForcing;
  }
  else if ( nameOfForcingType == "polyForcing" ) {
    forcingOption = ::polyForcing;
  }
  else if ( nameOfForcingType == "trigForcing" ) {
    forcingOption = ::trigForcing;
  }
  else {
    std::cerr << "Unknown forcing type = ["
	      << nameOfForcingType << "]" << std::endl;
    forcingOption = ::noForcing;
  }

  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update();

  Interpolant interpolant(cg);                                 // Make an interpolant

  Ogshow show( nameOfShowFile );                               // create a show file
  show.saveGeneralComment("Diffusion Equation Point Spread");    // save a general comment in the show file
  //  show.setFlushFrequency(10);                                  // flush file every 10 frames
    
  // create the grid function for the solution
  Range all;
  const int numberOfFieldComponents=1;
  realCompositeGridFunction u(cg,all,all,all,numberOfFieldComponents);
  u.setName("u");                                              // name the grid function

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

  // create operators  & use those with the grid function
  CompositeGridOperators operators(cg);                        // operators for a CompositeGrid
  u.setOperators(operators);                                 
  // operators.setOrderOfAccuracy(4);                          // for fourth order

  // auxiliary data
  Index I1,I2,I3, Ib1, Ib2, Ib3;

  // set initial data
  for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
    MappedGrid & mg = cg[ig];
    getIndex( mg.dimension(), I1,I2,I3);
    u[ig](I1,I2,I3) = exact( mg, I1,I2,I3, 0, 0.);
    q[ig](I1,I2,I3, uComponent) = u[ig](I1,I2,I3);
    q[ig](I1,I2,I3, uErrorComponent) = 0.;
    q[ig](I1,I2,I3, uExactComponent) = exact( mg, I1,I2,I3, 0, 0.);
  } //end for ig

  char buffer[80];                                           // buffer for sprintf
  int saveEvery        =1;
  int numberOfTimeSteps=5;
  for( int istep=0; istep<numberOfTimeSteps; istep++ )       // take some time steps
  {
    if( istep % saveEvery == 0 )  // save solution every 'saveEvery' steps
    {
      show.startFrame();                                         // start a new frame
      show.saveComment(0,sPrintF(buffer,"Here is solution %i",istep));   // comment 0 (shown on plot)
      show.saveComment(1,sPrintF(buffer,"  t=%f ",t));               // comment 1 (shown on plot)
      show.saveSolution( q );                                        // save the current grid function
    }
    //u+=dt*( -a*u.x() - b*u.y() + viscosity*(u.xx() + u.yy())); // take a time step with Euler's method
    //u+=dt*viscosity*(u.xx() + u.yy()); // take a time step with Euler's method
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.indexRange(), I1,I2,I3);

      //..HEAT EQ
      //q[ig](I1,I2,I3, uRHSComponent)     =  ( viscosity*(u[ig].xx()(I1,I2,I3) + u[ig].yy()(I1,I2,I3)));
      q[ig](I1,I2,I3, uRHSComponent)     =  viscosity*u[ig].laplacian()(I1,I2,I3);
      if (forcingOption > 0) {
	OGFunction & exact = *pExactSolution; // make a reference for readibility
	q[ig](I1,I2,I3, uRHSErrorComponent) = abs( q[ig](I1,I2,I3, uRHSComponent) 
			      -  viscosity*(  exact.xx(mg, I1,I2,I3,0,t)+exact.yy(mg,I1,I2,I3,0,t)) );
	q[ig](I1,I2,I3, uForcingComponent) =  ( exact.t( mg, I1,I2,I3, 0, t) 
						   -viscosity*( exact.xx(mg, I1,I2,I3,0,t)+exact.yy(mg,I1,I2,I3,0,t)));
	q[ig](I1,I2,I3, uRHSComponent) += q[ig](I1,I2,I3, uForcingComponent);
      }
      u[ig](I1,I2,I3) += dt* q[ig](I1,I2,I3, uRHSComponent); // EULER update

#if 0
      //..REACTION EQ
      q[ig](I1,I2,I3, uRHSComponent)     =  -a*u[ig](I1,I2,I3);
      if (forcingOption > 0) {
	q[ig](I1,I2,I3, uRHSErrorComponent) = abs( q[ig](I1,I2,I3, uRHSComponent) 
						   - (-a*exact(mg, I1,I2,I3,0,t)));
	q[ig](I1,I2,I3, uForcingComponent) =   exact.t( mg, I1,I2,I3, 0, t) - ( -a*exact(mg,I1,I2,I3,0,t));
	q[ig](I1,I2,I3, uRHSComponent) += q[ig](I1,I2,I3, uForcingComponent);
      }
      u[ig](I1,I2,I3) += dt* q[ig](I1,I2,I3, uRHSComponent); // EULER update
#endif
    } //end for ig

    t+=dt;
    u.interpolate();                                           // interpolate
    u.periodicUpdate();
    // apply a dirichlet BC on all boundaries:
#define ForBoundary(side,axis)   for( axis=0; axis<mg.numberOfDimensions(); axis++ ) \
                                 for( side=0; side<=1; side++ )
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      int axis, side;

      Index Ib1,Ib2,Ib3;
      Index Ig1,Ig2,Ig3;
      if (forcingOption>0 ) {
	OGFunction & exact = *pExactSolution; // make a reference for readibility
	
	ForBoundary(side,axis)   {
	  if( mg.boundaryCondition()(side,axis) > 0  )   {
	    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	    
	    u[ig](Ib1, Ib2, Ib3) = exact(mg, Ib1, Ib2, Ib3, 0, t);
	  }    
	} // end ForBoundary
      }
      else {
	cout << "setting neumann bc's\n";
	u.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,1.);
      }
    } // end for ig
#undef  ForBoundary
    // u.applyBoundaryCondition(0,BCTypes::extrapolate,BCTypes::allBoundaries,0.); // for 4th order

    u.finishBoundaryConditions();
  
    real maxErr = 0.;
    real maxErrGrid[ cg.numberOfComponentGrids() ];
    if (forcingOption >0 ) {
      OGFunction & exact = *pExactSolution; // make a reference for readibility
      for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
	MappedGrid & mg = cg[ig];
	Index Ig1,Ig2,Ig3;
	getIndex( mg.gridIndexRange(), Ig1,Ig2,Ig3, -1); // interior point errors only!!
	getIndex( mg.dimension(), I1,I2,I3);
	q[ig](I1,I2,I3, uComponent)      = u[ig](I1,I2,I3);
	q[ig](I1,I2,I3, uExactComponent) = exact(mg, I1,I2,I3, 0, t);

	q[ig](I1,I2,I3, uErrorComponent) = 0.;
	q[ig](Ig1,Ig2,Ig3, uErrorComponent) = u[ig](Ig1,Ig2,Ig3) - exact(mg,Ig1,Ig2,Ig3,0,t);
	where( mg.mask()(Ig1,Ig2,Ig3) > 0 ) {
	  maxErr = max( maxErr, max(abs( q[ig](Ig1,Ig2,Ig3,uErrorComponent))));
	}
      } //end for ig
      printf(" step %5i time= %8.3f max error = %10.4e\n", istep,t,  maxErr );
    }
    else {
      for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
	MappedGrid & mg = cg[ig];
	Index Ig1,Ig2,Ig3;
	getIndex( mg.gridIndexRange(), Ig1,Ig2,Ig3, -1); // interior point errors only!!
	getIndex( mg.dimension(), I1,I2,I3);
	q[ig](I1,I2,I3, uComponent)      = u[ig](I1,I2,I3);
	q[ig](I1,I2,I3, uExactComponent) = 0.;

	q[ig](I1,I2,I3, uErrorComponent) = 0.;
	q[ig](Ig1,Ig2,Ig3, uErrorComponent) = 0.;
	where( mg.mask()(Ig1,Ig2,Ig3) > 0 ) {
	  maxErr = max( maxErr, max(abs( q[ig](Ig1,Ig2,Ig3,uErrorComponent))));
	}
      } //end for ig
      printf(" step %5i time= %8.3f max error = %10.4e\n", istep,t,  maxErr );
    }
  }
  Overture::finish();          
  return 0;
}

