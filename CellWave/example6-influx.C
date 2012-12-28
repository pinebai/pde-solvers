#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#if 0
//.. looping calls to array interface, for initial data & rhs
void setInfluxBC( const double &tcomp,
		  const int&nd,    const int &ncomp,
		  const int &nd1a, const int &nd1b,
		  const int &nd2a, const int &nd2b,
		  const int &nd3a, const int &nd3b,
		  const int &n1a,  const int &n1b,
		  const int &n2a,  const int &n2b,
		  const int &n3a,  const int &n3b,
		  const double &xfirst, const double &qfirst, 
		  double &rhsfirst)
{
  const double *xArray   = &xfirst;
  const double *qArray   = &qfirst;
  double       *rhsArray = &rhsfirst;

  if( nd==3 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int cc=0, pc=1,hc=2;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      
      const double &p      = qArray[   GRIDINDEX( i,j,k, pc) ];
      double &rhs_p        = rhsArray[ GRIDINDEX( i,j,k, pc) ];
      
      //call single point rxn code -- depends on impl
      //computeRHS( c,h,p, b1, b2,   rhs_c, rhs_h, rhs_p, rhs_b1, rhs_b2,   x,y,z);
    }
  }

}
#endif

int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  printf(" ---------------------------------------------------------------------------- \n");
  printf("Solve: u.t + a*u.x + b*u.y = viscosity*( u.xx + u.yy ) on an Overlapping grid \n");
  printf("Save results in a show file, use plotStuff to view this file                  \n");
  printf(" ---------------------------------------------------------------------------- \n");

  //aString nameOfOGFile="Grids/cell20um2d_40_influx.hdf", nameOfShowFile="TEST/out-influx.show";
  //aString nameOfOGFile="Grids/box41-ip3Influx.hdf", nameOfShowFile="TEST/out-influx.show";
  aString nameOfOGFile="Grids/without-3d-influx.hdf", nameOfShowFile="TEST/OUT2/out-influx.show";
  cout << "example6>> Enter the name of the (old) overlapping grid file:" << endl;
  //cin >> nameOfOGFile;
  cout << nameOfOGFile << "\n";

  cout << "example6>> Enter the name of the (new) show file (blank for none):" << endl;
  //cin >> nameOfShowFile;
  cout << nameOfShowFile << "\n";

  real t=0, dt=          0.001;//0004;                              // initialize time and time step
  real a=0., b=0., viscosity=10.;                               // initialize parameters
  const int ip3InfluxBoundaryID=3;
  
  double ip3FluxMagnitude=10000.;
  //double x0_ip3=0., y0_ip3=12.;
  //double x0_ip3=10., y0_ip3=12.;
  double x0_ip3=45., y0_ip3= 45.;
  double ip3Width=3.;
    
  double pi=4.*atan(1.0);

  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update();
  cg.update( MappedGrid::THEvertex | MappedGrid::THEvertexBoundaryNormal );

  Interpolant interpolant(cg);                                 // Make an interpolant

  Ogshow show( nameOfShowFile );                               // create a show file
  show.saveGeneralComment("Convection Diffusion Equation");    // save a general comment in the show file
  show.setFlushFrequency(1);                                  // flush file every 10 frames
    
  CompositeGridOperators operators(cg);                        // operators for a CompositeGrid
  // operators.setOrderOfAccuracy(4);                          // for fourth order

  Range all;
  realCompositeGridFunction u(cg,all,all,all,2);               // create a grid function 
  u.setOperators(operators);                                 
  u.setName("solution");                                       // name the grid function
  u.setName("u",0);
  u.setName("rhs",1);

  //u=0.;                                                        // initial condition
  for (int igrid = 0; igrid < cg.numberOfComponentGrids(); igrid++) {
    Index I1,I2,I3;
    MappedGrid &mg              = cg[igrid];
    getIndex( mg.gridIndexRange(), I1,I2,I3);
    realArray & xArray = mg.vertex();  // array of vertices
	
#define XCV xArray(I1,I2,I3,axis1)
#define YCV xArray(I1,I2,I3,axis2)
    //u[igrid](I1,I2,I3,0) = xArray(I1,I2,I3,axis2);
    //u[igrid](I1,I2,I3,0) = ip3FluxMagnitude*exp( -( pow( XCV-x0_ip3,2.) + pow(YCV-y0_ip3,2.))/(4.*ip3Width) ) 
    //     /( 4.0*pi*ip3Width ); 
    //u[igrid](I1,I2,I3,0) =0.;
  }
  u=0.;

  char buffer[80];                                             // buffer for sprintf
  int numberOfTimeSteps=200000;
  int saveEvery=1;
  const int uc=0, fc=1;
  Index I1,I2,I3,If1,If2,If3, Ib1,Ib2,Ib3, Ig1,Ig2,Ig3;
  for( int i=0; i<numberOfTimeSteps; i++ )                    // take some time steps
  {
    //printf("--step %d, t=%f\n", i, t);
    if( i % saveEvery == 0 )  // save solution every 10 steps
    {
      printf("--step %d, t=%f\n", i, t);
      show.startFrame();                                         // start a new frame
      show.saveComment(0,sPrintF(buffer,"Here is solution %i",i));   // comment 0 (shown on plot)
      show.saveComment(1,sPrintF(buffer,"  t=%e ",t));               // comment 1 (shown on plot)
      show.saveSolution( u );                                        // save the current grid function
    }
    for( int igrid=0; igrid < cg.numberOfComponentGrids(); ++igrid) {
      //u+=dt*( -a*u.x() - b*u.y() + viscosity*(u.xx() + u.yy())); // take a time step with Euler's method
      MappedGrid &mg              = cg[igrid];
      realMappedGridFunction &umg = u[igrid];
      getIndex( mg.indexRange(), I1,I2,I3);
      umg(I1,I2,I3,fc) = viscosity*umg.laplacian()(I1,I2,I3,uc);
      umg(I1,I2,I3,uc) += dt*umg(I1,I2,I3,fc);
      //u[igrid](I1,I2,I3)
    }
    t+=dt;
    u.interpolate();                                           // interpolate
    // apply a dirichlet BC on all boundaries:
    //u.applyBoundaryCondition(0,BCTypes::dirichlet,BCTypes::allBoundaries,0.);
    u.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries, 0. ); //noflux
    for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid ) {
      MappedGrid &mg              = cg[igrid];
      realMappedGridFunction &umg = u[igrid];
      const int Start=0, End=1; // limits for 'side'
      for( int axis=0; axis< mg.numberOfDimensions(); ++axis) {
	for( int side=Start; side<=End; ++side ) {

	  int  ibc= mg.boundaryCondition()(side,axis);
	  if( ibc == ip3InfluxBoundaryID ) {
	    //printf ("---------------- setting influx bc on grid %d\n",igrid);
	    
	    const int thisBCIndex=BCTypes::boundary1+side+2*axis;
	    realArray & normal = mg.vertexBoundaryNormal(side,axis);
	    //normal.display("NORMAL VECTOR");
	    
	    realArray & uArray     = u[igrid];
	    realArray & xArray = mg.vertex();  // array of vertices
	    getBoundaryIndex(mg.gridIndexRange(), side, axis, Ib1,Ib2,Ib3);
	    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	    realArray ip3Flux(Ib1,Ib2,Ib3);
#define XC xArray(Ib1,Ib2,Ib3,axis1)
#define YC xArray(Ib1,Ib2,Ib3,axis2)
	    
	    ip3Flux =  ip3FluxMagnitude*exp( -( pow( XC-x0_ip3,2.) + pow(YC-y0_ip3,2.))/(4.*ip3Width) )  /( 4.0*pi*ip3Width ); 
	    //ip3Flux = XC;
	    double fluxMin= min(abs(ip3Flux)), fluxMax= max(abs(ip3Flux));
	    printf( " t=%f, grid=%2d, side=%d, axis=%d  ", t,igrid,side,axis);
	    printf( " ip3Width = %f, ip3FluxMagnitude= %f x0= %f, y0= %f, fmin= %f, fmax= %f\n", 
		    ip3Width, ip3FluxMagnitude,x0_ip3, y0_ip3,fluxMin,fluxMax);
	    
	    const int component=0; 
	    umg.applyBoundaryCondition(component, BCTypes::neumann, thisBCIndex, ip3Flux);
	    //umg.applyBoundaryCondition(component, BCTypes::dirichlet, thisBCIndex, ip3Flux);
	    //realArray z=umg(Ig1,Ig2,Ig3); z.display("u on the ghostline");
	  }    
	}      
      }//end for side
    }//end for axis
    
    // u.applyBoundaryCondition(0,BCTypes::extrapolate,BCTypes::allBoundaries,0.); // for 4th order
    u.finishBoundaryConditions();
  }

  Overture::finish();          
  return 0;
    
}
