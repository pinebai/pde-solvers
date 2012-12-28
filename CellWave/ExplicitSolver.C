#include "SolverLiRinzel.h"

using namespace CellWave;

void SolverLiRinzel::
setup()
{
  GenericSolver::genericSetup();
}

void SolverLiRinzel::
initialData()
{
  GenericSolver::genericInitialData();
}

void SolverLiRinzel::
solve()
{
  //..Setup
  assert( pChem != NULL );
  CellWave::GenericReaction &chem = *pChem; // shorthand

  //....set parameters -- IN Solver subclass
  //FIXME: parameters should be set from a parameter file.
  //FIXME: ALSO, hookup getDT to compute timestep
  //FIXME: what about the timestep from the reaction?? 

  real tcomp=0, dt=.1;    // works with D=378 um^2/sec, dx_min= 5 um   //FIXME:

  //dt = timeStepSize;                      // from command line     //FIXME
  //numberOfTimeSteps = totalNumberOfSteps; // from command line     //FIXME
  //saveEvery         = saveEveryNthFrame;  // from command line     //FIXME

  realCompositeGridFunction pconcentration, cconcentration, hconcentration;
  realCompositeGridFunction rhs_p, rhs_c, rhs_h; // components, and their RHS
  cconcentration.link(q,Range(cc,cc));  q.setName("Ca2+ (cytosol)",0); //FIXME
  pconcentration.link(q,Range(pc,pc));  q.setName("IP3",1);//FIXME
  hconcentration.link(q,Range(hc,hc));  q.setName("h Inhibitor",2);//FIXME

  rhs_c.link(f,Range(cc,cc));  f.setName("RHS for Ca2+",0);//FIXME
  rhs_p.link(f,Range(pc,pc));  f.setName("RHS for IP3",1);//FIXME
  rhs_h.link(f,Range(hc,hc));  f.setName("RHS for h Inhibitor",2);//FIXME

  //.... viscosity
  realArray viscosity( chem.getNumberOfSpecies());
  viscosity( cc ) =  chem.getDiffusionCoefficient( int(cc) );
  viscosity( hc ) =  chem.getDiffusionCoefficient( int(hc) );
  viscosity( pc ) =  chem.getDiffusionCoefficient( int(pc) );

  real cfl =0.5;
  real dtDiffusion = getDiffusionDT(cfl, chem.getNumberOfSpecies(), 
				    viscosity, cg );
  //real dtDiffusion;
  dt = dtDiffusion; // or maybe dt = min(dt, dtDiffusion)
  printf("*** dt= %8.4e, visc(p)= %8.4e, visc(c)= %8.4e\n",
	 dt, viscosity(pc), viscosity( cc));

  //..start timestepping
  real wallTime=0.;
  int iWallTimeSteps=0;
  real totalWallTime=0.;
  char buffer[80];

  const int &saveEvery = data.saveEveryNthFrame;

  { //....debug output
    real dp= chem.getDiffusionCoefficient( int(pc) );
    real dc= chem.getDiffusionCoefficient( int(cc) );
    real lengthScale = chem.getMaximumLengthScale();
    
    std::string offon[] = {"off", "on"};
    const char *isCaDiffusive 
      = offon[ int(chem.isDiffusive(int(cc))) ].c_str();
    const char *isIP3Diffusive 
      = offon[ int(chem.isDiffusive(int(pc))) ].c_str();
    
    printf("  nondim diffusion [c]= %8.4g(%s), [p]= %8.4g(%s),L=%8.4g \n",
	   dt*dp/lengthScale, isIP3Diffusive,
	   dt*dc/lengthScale, isCaDiffusive,  chem.getMaximumLengthScale());
  } //....end debug output
  
  for( int istep=0; istep<=data.numberOfTimeSteps+1; istep++ )
  {
    iWallTimeSteps++;
    real wallTimeThis=getCPU();
    Index I1,I2,I3;

    if( istep % saveEvery == 0 )  // save solution every 'saveEvery' steps
    {
      double cmax=0., pmax=0., hmax=0.;
      for (int igrid=0; igrid < cg.numberOfComponentGrids(); ++igrid) {
	//where() ... add 'where mask>0' to get correct ip3 max 

	getIndex(cg[igrid].indexRange(),I1,I2,I3);
        where( cg[igrid].mask()(I1,I2,I3)!=0 ) {

	  cmax = max(cmax, max(abs(cconcentration[igrid](I1,I2,I3))) );
	  pmax = max(pmax, max(abs(pconcentration[igrid](I1,I2,I3))) );
	  hmax = max(hmax, max(abs(hconcentration[igrid](I1,I2,I3))) );
	}
      }
      
      real avWallTime=wallTime/iWallTimeSteps;      
      
      printf("..step %5i (t=%6.2f): Max in c=%8.3e, p=%8.3e, h=%8.3e, CPU time/step %8.4e\n", 
	     istep, tcomp, cmax, pmax, hmax, avWallTime);
      wallTime =0.;     // reset
      iWallTimeSteps=0;

      show.startFrame();                                               // start a new frame
      show.saveComment(0,sPrintF(buffer,"Here is solution %i",istep)); // comment 0 (shown on plot)
      show.saveComment(1,sPrintF(buffer,"  t=%e ",tcomp));             // comment 1 (shown on plot)
      show.saveSolution( q );                                          // save the current grid function
    }

    //
    //..Perform one time step, update each grid
    //
    f = 0.;
    for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
      //Index I1, I2,I3;
      getIndex( cg[igrid].gridIndexRange(), I1,I2,I3);
      
      MappedGrid & mg = cg[igrid];
      realArray & qArray   = q[igrid];
      realArray & rhsArray = f[igrid];
      //realArray & dudtg = dudt[igrid];
      realArray & xArray = mg.vertex();  // array of vertices
      
      const IntegerArray & d   = mg.dimension();
      const IntegerArray & gir = mg.gridIndexRange();
      const int nd             = cg.numberOfDimensions();
      const int ncomp          = chem.getNumberOfSpecies();

      //  (1) Chemistry: (Laplace terms imposed in step (2))
      //        dc/dt = beta*lam( (ip3 rel. )* (ca rel.) h^3 )(C_er - c) -nu_P *( )) + D*Laplace c
      //        dp/dt = Di*Laplace p - k_i *p 
      //      .... call explicit loop here: arguments are the info for a Fortran style loop
      chem.callRHSLoop(tcomp, 
		       nd, ncomp,
		       d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2), 
		       gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2), 
		       *xArray.getDataPointer(),*qArray.getDataPointer(),
		       *rhsArray.getDataPointer() );      

      //  (2) Diffusion
      real dp= chem.getDiffusionCoefficient( int(pc) );
      real dc= chem.getDiffusionCoefficient( int(cc) );
      realArray tempLaplace(I1,I2,I3);

      //   .... IP3 diffusion:
      if( chem.isDiffusive( int(pc) )) {
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,pc);
	rhsArray(I1,I2,I3,pc) += dp*tempLaplace(I1,I2,I3); 
      }

      //    .... Ca2+ diffusion:
      if( chem.isDiffusive( int(pc) )) {
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,cc);
	rhsArray(I1,I2,I3,cc) += dc*tempLaplace(I1,I2,I3); 
      }

      //    .... h diffusion:
      //         --none--

    }; //end for igrid

    //  (3) Update solution on next time (Adams-Bashforth, 2nd order, fully explicit) //FIXME -- adaptive timestepping
    if ( istep == 0 ) fp=f; // first step euler.
    //qn = q + .5*dt*( 3.*f - fp );
    //  (4) Update previous time levels
    //q=qn; 
    //fp=f;

    for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
      Index I1,I2,I3;
      const int lastcomp=hc;
      Range active(0, lastcomp);
      //getIndex( cg[igrid].gridIndexRange(), I1,I2,I3);
      
      qn[igrid](I1,I2,I3,active) = q[igrid](I1,I2,I3,active)
	+ .5*dt*( 3.*f[igrid](I1,I2,I3, active) - fp[igrid](I1,I2,I3, active) );

      //  (4) Update previous time levels      
      fp[igrid]=f[igrid];
      q[igrid](I1,I2,I3,active) = qn[igrid](I1,I2,I3,active);

      q[igrid](I1,I2,I3, cc+3) =  rhs_c[igrid](I1,I2,I3);
      q[igrid](I1,I2,I3, pc+3) =  rhs_p[igrid](I1,I2,I3);
      q[igrid](I1,I2,I3, hc+3) =  rhs_h[igrid](I1,I2,I3);

    }
    //  (5) Set bc's: internal (interpolate), and physical

    q.periodicUpdate();
    q.interpolate();                                           // interpolate

    // apply a neumann BC on all boundaries -- h does not get a bc
    //cconcentration.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    //pconcentration.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    q.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    q.finishBoundaryConditions();

    // (6) update time to t_(n+1)
    tcomp +=dt;

    wallTimeThis = getCPU() - wallTimeThis;
    wallTime += wallTimeThis;
    totalWallTime += wallTimeThis;
  }
}

//
// .. utilities
//

real SolverLiRinzel::
getTime()
{
  return( tcomp );
}

real SolverLiRinzel::
getTotalWallTime()
{
  return( 0. );
}

int  SolverLiRinzel::
getNumberOfTimeSteps()
{
  return( ktime );
}

