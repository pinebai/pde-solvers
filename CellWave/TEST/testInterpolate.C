#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "OGTrigFunction.h"
#include "OGPolyFunction.h"

#include "interpolateFluxBoundary.h"

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
}



int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  printf(" ---------------------------------------------------------------------------- \n");
  printf("   try solving u_t = u_xx + u_yy with flux jump bcs.\n");
  printf(" ---------------------------------------------------------------------------- \n");

  aString nameOfOGFile, nameOfShowFile, nameOfFluxBCType;
  const int physicalBoundaryID=1, fluxBoundaryID=2; // for Overture grids

  if (argc>2) {
    nameOfOGFile       =  argv[1];
    nameOfShowFile     =  argv[2];
    nameOfFluxBCType   =  argv[3];
  }
  else {
    std::cerr << "usage: "<<argv[0]
	      << " <grid.hdf> <output.show> <fluxBCType>"
	      << endl;
    exit(1);
  }

  // parameter values
  real t=0, dt=.0001;                                     // initialize time and time step
  real a=1., b=1., viscosity=.1;                          // initialize parameters
  ::fluxBCType fluxBC= ::noFluxBC;
  ::forcingTypes forcingOption= ::noForcing;          //::trigForcing;

  if ( nameOfFluxBCType == "noFlux" ) {
    fluxBC = ::noFluxBC;
  }
  else if ( nameOfFluxBCType == "interpolateFlux" ) {
    fluxBC = ::interpolateFluxBC;
  }
  else {
    std::cerr << "Unknown fluxBC type = ["
	      << nameOfFluxBCType << "]" << std::endl;
    fluxBC = ::noFluxBC;
  }

  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update();
  cg.update( CompositeGrid::THEvertex | CompositeGrid::THEvertexBoundaryNormal );

  Interpolant interpolant(cg);                            // Make an interpolant

  Ogshow show( nameOfShowFile );                          // create a show file
  show.saveGeneralComment("test internal flux/jump bc");  // save a general comment in the show file
  show.setFlushFrequency(50);                         // flush file every 10 frames

  // create the grid function for the solution
  Range all;
  const int numberOfFieldComponents=1;
  realCompositeGridFunction u(cg,all,all,all,numberOfFieldComponents);
  u.setName("u");                                              // name the grid function

  enum{ uComponent=0, uExactComponent=1, uErrorComponent=2, 
	  uRHSComponent=3, uForcingComponent=4, uRHSErrorComponent=5, numberOfOutputFields};
  realCompositeGridFunction q(cg,all,all,all,numberOfOutputFields);
  realCompositeGridFunction unext(cg,all,all,all,1);
  q.setName("u");
  q.setName("u", uComponent);
  q.setName("exact", uExactComponent);
  q.setName("pointwise error", uErrorComponent);
  q.setName("rhs", uRHSComponent );
  q.setName("forcing term", uForcingComponent);
  q.setName("rhs error", uRHSErrorComponent);

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
    real fx= 2., fy=2., fz = 2., ft=2.; // note: fz is not used in 2D
    // defines cos(fx*pi*x)*cos(fy*pi*y)*cos(fz*pi*z)*cos(ft*pi*t)
    pExactSolution = new OGTrigFunction( fx,fy,fz,ft);
  }
  else if (forcingOption != ::noForcing ) {
    std::cerr << "Unknown forcing option ="<< forcingOption
	      << std::endl;
    forcingOption = ::noForcing;
  }
  //OGFunction & exact = *pExactSolution; // make a reference for readibility

  // create operators  & use those with the grid function
  CompositeGridOperators operators(cg);                        // operators for a CompositeGrid
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
      u[ig](I1,I2,I3) = exact( mg, I1,I2,I3, 0, 0.);
      q[ig](I1,I2,I3, uComponent) = u[ig](I1,I2,I3);
      q[ig](I1,I2,I3, uExactComponent) = exact( mg, I1,I2,I3, 0, 0.);
      q[ig](I1,I2,I3, uErrorComponent) = 0.;

    } //end for ig
  } else {
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      MappedGrid & mg = cg[ig];
      getIndex( mg.dimension(), I1,I2,I3);
      u[ig](I1,I2,I3) =  double(ig);//0.;
      q[ig](I1,I2,I3, uComponent) = u[ig](I1,I2,I3);
      q[ig](I1,I2,I3, uExactComponent) = 0.;
      q[ig](I1,I2,I3, uErrorComponent) = 0.;
    }
  }
      
  char buffer[80];                                             // buffer for sprintf
  int saveEvery        =10;
  int numberOfTimeSteps=300;
  //int numberOfTimeSteps=2;
  for( int istep=0; istep<numberOfTimeSteps; istep++ )                    // take some time steps
  {
    if (localDebug &2) printf("\n\n");
    if (localDebug &1)
      printf("..............step %d..............\n",istep);
    if (localDebug &2) printf("\n");

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
      q[ig](I1,I2,I3, uRHSComponent)     =   viscosity*(u[ig].laplacian()(I1,I2,I3));
      if (forcingOption > 0) {
	OGFunction & exact = *pExactSolution; // make a reference for readibility
	q[ig](I1,I2,I3, uRHSErrorComponent) = abs( q[ig](I1,I2,I3, uRHSComponent) 
						   -  viscosity* exact.laplacian(mg, I1,I2,I3,0,t));
	q[ig](I1,I2,I3, uForcingComponent) =  ( exact.t( mg, I1,I2,I3, 0, t) 
						-viscosity*exact.laplacian(mg, I1,I2,I3,0,t));
	q[ig](I1,I2,I3, uRHSComponent) += q[ig](I1,I2,I3, uForcingComponent);
      }
      //u[ig](I1,I2,I3) += dt* q[ig](I1,I2,I3, uRHSComponent); // EULER update
      unext[ig](I1,I2,I3) = u[ig](I1,I2,I3) +dt* q[ig](I1,I2,I3, uRHSComponent); // EULER update

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
    
    unext.interpolate();                                           // interpolate
    unext.periodicUpdate();
    //.. extrapolate ghostlines so flux/jump bcs ok later, then copy to u
    //.... unext remains a copy of u with ghostline data extrapolated
    unext.applyBoundaryCondition(0,BCTypes::extrapolate,BCTypes::allBoundaries); 
    for(int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      getIndex( cg[ig].dimension(), I1,I2,I3);
      u[ig](I1,I2,I3) = unext[ig]( I1,I2,I3 );
    }


    // apply exact bcs
#define ForBoundary(side,axis)   for( axis=0; axis<mg.numberOfDimensions(); axis++ ) \
                                 for( side=0; side<=1; side++ )
    for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
      //for ( int ig= cg.numberOfComponentGrids()-1; ig>=0; --ig ) { //inverted, first grid 1, then grid 0
      if (localDebug &2 )
      printf("--BC's for grid %d --\n", ig);
      MappedGrid & mg = cg[ig];
      MappedGridOperators &opmg = operators[ig];
      int axis, side;

      Index Ib1,Ib2,Ib3;
      Index Ig1,Ig2,Ig3;
      if (forcingOption>0 ) {
	OGFunction & exact = *pExactSolution; // make a reference for readibility
	
	ForBoundary(side,axis)   {
	  //const int physicalBoundaryID=1, fluxBoundaryID=2; // for Overture grids
	  if( mg.boundaryCondition()(side,axis) == physicalBoundaryID  ) { // physical boundaries
	    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	    
	    u[ig](Ib1, Ib2, Ib3) = exact(mg, Ib1, Ib2, Ib3, 0, t);
	  }    
	  else if( mg.boundaryCondition()(side,axis) == fluxBoundaryID  ) { // flux/jump bc (=concentration channel bc)
	    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	    
	    u[ig](Ib1, Ib2, Ib3) = exact(mg, Ib1, Ib2, Ib3, 0, t);          // **FIXME -- change this to du/dn = F[ u ]
	  }    
	} // end ForBoundary
      } 
      else { // else forcingOption == 0 (noForcing)
	
	realMappedGridFunction &umg=u[ig];
	assert( cg.numberOfComponentGrids() == 2 );
	const int igOther = 1-ig; //works only for 2 grids
	ForBoundary(side,axis)   {
	  //const int physicalBoundaryID=1, fluxBoundaryID=2; // for Overture grids

	  if( mg.boundaryCondition()(side,axis) == fluxBoundaryID  ) { // flux/jump bc
	    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	    const realArray &zz=mg.vertex();
	    int nAxes        =cg.numberOfDimensions();
	    int zaxis        =axis3;
	    if (nAxes==2) zaxis=axis2; // to make sure we don't seg.fault
	    const realArray &xb  = mg.vertex()(Ib1,Ib2,Ib3, axis1);
	    const realArray &yb  = mg.vertex()(Ib1,Ib2,Ib3, axis2);
	    const realArray &zb  = mg.vertex()(Ib1,Ib2,Ib3, zaxis);

	    int axisLength[3]={xb.getLength(axis1), xb.getLength(axis2), xb.getLength(axis3)};
	    int nInterpPoints= axisLength[axis1]*axisLength[axis2];
	    if (cg.numberOfDimensions()==3) nInterpPoints = nInterpPoints* axisLength[axis3];
	    
	    realArray xyInt( Ib1,Ib2,Ib3, nAxes );
	    Index All;
	    xyInt(All, All, All, axis1) = xb;
	    xyInt(All, All, All, axis2) = yb;
	    if( nAxes == 3) {
	      xyInt(All, All, All, axis3) = zb;
	    }
	    xyInt.reshape(nInterpPoints, nAxes);
	    	    
	    realArray       &uf     = u[ig];
	    realArray       &uOther = unext[igOther];
	    //realArray jump(Ib1,Ib2,Ib3), dummy(Ig1,Ig2,Ig3);
	    //realArray jump(nInterpPoints), dummy(nInterpPoints); 
	    realArray jump(nInterpPoints);

	    if(localDebug &8)   uOther.display(" u other ");
	    interpolatePointsPF(ig,xyInt, unext, jump);  //.. interpolate values from unext...
	    jump.reshape(Ib1,Ib2,Ib3);
	    if(localDebug &4)   jump.display("value"); 
	    jump = 100*( jump - uf(Ib1,Ib2,Ib3));
	    if(localDebug &4)    jump.display("jump");
	    //printf("========================================================================\n");

	    //.. ibc chosen to set a specified boundary (side,axis)
	    const int ibc=BCTypes::boundary1+side+2*axis;// ... take jumps from unext --> impose as flux in u
	    opmg.applyBoundaryCondition( umg, 0, BCTypes::neumann, ibc, jump ); 

	  }
	} // end ForBoundary
      }
    } // end for ig
#undef  ForBoundary

    if (forcingOption == ::noForcing ) { // set all physical boundaries to No-Flux
      u.applyBoundaryCondition(0,BCTypes::neumann, physicalBoundaryID,0.);
    }

    u.periodicUpdate();
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

