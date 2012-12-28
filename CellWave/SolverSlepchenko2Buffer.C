#include "CellWave.h"
#include "SolverSlepchenko2Buffer.h"
#include "ReactionSlepchenko2Buffer.h"

using namespace CellWave;

bool  SolverSlepchenko2Buffer::
readParameterFile( CellWave::ParameterReader &param)
{
  throw "error"; //should not call this
  assert( pChem != NULL );
  ReactionSlepchenko2Buffer &chem = *pChem; // shorthand

  bool ok=  chem.readParameterFile( param );
  if (ok) {
    GenericSolver::readParameterFile( param );
    GenericSolver::updateModel();
  }
  return ok;
}

#if 0
void SolverSlepchenko2Buffer::
setup()
{
  CellWave::ReactionFactory rxnFactory;
  SolverSlepchenko2Buffer::pChem  
    = (ReactionSlepchenko2Buffer*) rxnFactory.getReaction("Slepchenko2Buffer"); //FIXME: static_cast<>
  GenericSolver::genericSetup( SolverSlepchenko2Buffer::pChem,  );
  GenericSolver::updateModel();
}
#endif

void SolverSlepchenko2Buffer::
setup(CellWave::ParameterReader &param)
{
  CellWave::ReactionFactory rxnFactory; // here -- should get the name from the param file!!
  SolverSlepchenko2Buffer::pChem  
    = (ReactionSlepchenko2Buffer*) rxnFactory.getReaction("Slepchenko2Buffer");
  GenericSolver::genericSetup( SolverSlepchenko2Buffer::pChem, param );
  //bool fileOK = this->readParameterFile( param );

}

//void SolverSlepchenko2Buffer::
//finish() -- in GenericSolver (base class)

void SolverSlepchenko2Buffer::
initialData()
{
  GenericSolver::genericInitialData();
}

void SolverSlepchenko2Buffer::
solve()
{
  //..Setup
  assert( pChem != NULL );
  CellWave::GenericReaction &chem = *pChem; // shorthand
  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  //FIXME: ALSO, hookup getDT to compute timestep
  //FIXME: what about the timestep from the reaction?? 

  tcomp=0;
  //.... viscosity
  realArray viscosity( chem.getNumberOfSpecies());
  viscosity( cc ) =  chem.getDiffusionCoefficient( int(cc) );
  viscosity( hc ) =  chem.getDiffusionCoefficient( int(hc) );
  viscosity( pc ) =  chem.getDiffusionCoefficient( int(pc) );
  viscosity( b1c ) =  chem.getDiffusionCoefficient( int(b1c) );
  viscosity( b2c ) =  chem.getDiffusionCoefficient( int(b2c) );

  real cfl =0.5;
  real dtDiffusion = getDiffusionDT(cfl, chem.getNumberOfSpecies(), 
				    viscosity, cg );
  //real dtDiffusion;
  real dt = dtDiffusion; // or maybe dt = min(dt, dtDiffusion)
  if (dt > data.timeStepSize ) {
    dt = data.timeStepSize;
  }
  DPrintf(DebugSolver,"*** dt= %8.4e, visc(p)= %8.4e, visc(c)= %8.4e, ",
	  dt, viscosity(pc), viscosity( cc));
  DPrintf(DebugSolver,"visc(b1)= %8.4e, visc(b2)= %8.4e\n",
	  viscosity(b1c), viscosity( b2c));
  
  //..start timestepping
  GenericSolver::wallTime=0.;
  GenericSolver::iWallTimeSteps=0;
  GenericSolver::totalWallTime=0.;

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
    
    DPrintf(DebugSolver,"  nondim diffusion [p]= %8.4g(%s), [c]= %8.4g(%s),L=%8.4g \n",
	   dt*dp/lengthScale, isIP3Diffusive,
	   dt*dc/lengthScale, isCaDiffusive,  chem.getMaximumLengthScale());
  } //....end debug output
  
  for( ktime=0; ktime<=data.numberOfTimeSteps+1; ++ktime )
  {
    iWallTimeSteps++;
    real wallTimeThis=getCPU();
    Index I1,I2,I3; // interior
    Index If1,If2,If3; // full range

    GenericSolver::outputStepStart(ktime,tcomp);

    if( ktime % saveEvery == 0 )  // save solution every 'saveEvery' steps
    {
      double cmax=0., pmax=0., hmax=0.;
      double cmin=1e100, pmin=1e100, hmin=1e100;
      double b1min=1e100, b2min=1e100;
      double b1max=0., b2max=0.;

      for (int igrid=0; igrid < cg.numberOfComponentGrids(); ++igrid) {
	//where() ... add 'where mask>0' to get correct ip3 max 

	getIndex(cg[igrid].indexRange(),I1,I2,I3);
        where( cg[igrid].mask()(I1,I2,I3)!=0 ) {

	  cmax = max(cmax, max(abs(q[igrid](I1,I2,I3,cc))) );
	  pmax = max(pmax, max(abs(q[igrid](I1,I2,I3,pc))) );
	  hmax = max(hmax, max(abs(q[igrid](I1,I2,I3,hc))) );
	  b1max= max(b1max, max(abs(q[igrid](I1,I2,I3,b1c))));
	  b2max= max(b2max, max(abs(q[igrid](I1,I2,I3,b2c))));

	  cmin = min(cmin, min(abs(q[igrid](I1,I2,I3,cc))) );
	  pmin = min(pmin, min(abs(q[igrid](I1,I2,I3,pc))) );
	  hmin = min(hmin, min(abs(q[igrid](I1,I2,I3,hc))) );
	  b1min= min(b1min, min(abs(q[igrid](I1,I2,I3,b1c))));
	  b2min= min(b2min, min(abs(q[igrid](I1,I2,I3,b2c))));
	}
      }
      
      real avWallTime=wallTime/iWallTimeSteps;      
      
      DPrintf(DebugSolver,
	      "..step %3i (t=%4.2f):c=[%6.3g,%6.3g] p=[%6.3g,%6.3g] h=[%6.3g,%6.3g] ",
	     ktime, tcomp, cmin,cmax, pmin,pmax, hmin,hmax );
      DPrintf(DebugSolver,
	      "b1=[%6.3g,%6.3g] b2=[%6.3g,%6.3g]\n",
	      b1min,b1max,b2min,b2max);

      DPrintf(PRINT,
	      "..step %3i (t=%4.2f):c=[%6.3g,%6.3g] p=[%6.3g,%6.3g] h=[%6.3g,%6.3g] ",
	     ktime, tcomp, cmin,cmax, pmin,pmax, hmin,hmax );
      DPrintf(PRINT,
	      "b1=[%8.3g,%8.3g] b2=[%8.3g,%8.3g]\n",
	      b1min,b1max,b2min,b2max);

      wallTime =0.;     // reset
      iWallTimeSteps=0;

      show.startFrame();                                               // start a new frame
      show.saveComment(0,sPrintF(buffer,"Here is solution %i",ktime)); // comment 0 (shown on plot)
      show.saveComment(1,sPrintF(buffer,"  t=%e ",tcomp));             // comment 1 (shown on plot)
      show.saveSolution( q );                                          // save the current grid function
      show.endFrame();
    }

    //
    //..Perform one time step, update each grid
    //
    f = 0.;
    for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
      //Index I1, I2,I3;
      getIndex( cg[igrid].gridIndexRange(), I1,I2,I3);
      getIndex( cg[igrid].dimension(),      If1,If2,If3);
      
      MappedGrid & mg = cg[igrid];
      realArray & qArray     = q[igrid];
      realArray & rhsArray = f[igrid];
      realArray & rhsArray_p = fp[igrid];
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

      //rhsArray.display("rhsArray");
      if( nucleusClass.getNumberOfNuclei() > 0 ) {
	realArray &maskArray = GenericSolver::getNucleusMaskArray( igrid );
	for (int ispecies=0; ispecies< ncomp; ++ispecies ) {
	  rhsArray(If1,If2,If3,ispecies) 
	    = maskArray(If1,If2,If3)*rhsArray(If1,If2,If3,ispecies);
	}
      }

      //  (2) Diffusion
      real dp= chem.getDiffusionCoefficient( int(pc) );
      real dc= chem.getDiffusionCoefficient( int(cc) );
      real db2=chem.getDiffusionCoefficient( int(b2c));
      realArray tempLaplace(I1,I2,I3);

      //   .... IP3 diffusion:
      bool alwaysDiffusion=true;
      if( alwaysDiffusion || chem.isDiffusive( int(pc) )) {
	DPrintf(DebugSolver,"diffuse IP3, nu=%8.3g\n", dp);
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,pc);
	rhsArray(I1,I2,I3,pc) += dp*tempLaplace(I1,I2,I3); //FIXME
      }

      //    .... Ca2+ diffusion:
      if( alwaysDiffusion || chem.isDiffusive( int(cc) )) {
	DPrintf(DebugSolver,"diffuse Ca, nu=%8.3g\n", dc);
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,cc);
	rhsArray(I1,I2,I3,cc) += dc*tempLaplace(I1,I2,I3);   //FIXME
      }

      //    .... h diffusion:
      //         --none--

      //    .... b1 diffusion:
      //         --none--

      //    .... b2 diffusion:
      if( alwaysDiffusion || chem.isDiffusive( int(b2c) )) {
	DPrintf(DebugSolver,"diffuse Buffer2, nu=%8.3g\n", db2);
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,b2c);
	rhsArray(I1,I2,I3,b2c) += db2*tempLaplace(I1,I2,I3);   //FIXME
      }
      
      //  (3) Update solution on next time (Adams-Bashforth, 2nd order, fully explicit) //FIXME -- adaptive timestepping
      if( ktime == 0 )  rhsArray_p = rhsArray;
	
      //..AB2 step
      DPrintf(DebugSolver,"dt=%8.3g\n",dt);
      for( int ic=0; ic<chem.getNumberOfSpecies(); ++ic ) {
	qArray(I1,I2,I3,ic) =  qArray(I1,I2,I3,ic) 
	  + .5*dt*( 3.*rhsArray(I1,I2,I3,ic) - rhsArray_p(I1,I2,I3,ic));
      }

      //  (4) Update previous time levels      
      rhsArray_p = rhsArray;

      //  (4B) DEBUG: copy RHS into solution array
      //qArray(I1,I2,I3, cc+3) =  rhsArray(I1,I2,I3,cc);
      // qArray(I1,I2,I3, pc+3) =  rhsArray(I1,I2,I3,pc);
      //qArray(I1,I2,I3, hc+3) =  rhsArray(I1,I2,I3,hc);

    }; //end for igrid

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


real SolverSlepchenko2Buffer::
getTime()
{
  return( GenericSolver::getTime() );
}

real SolverSlepchenko2Buffer:: 
getTotalWallTime()
{
  return( GenericSolver::getTotalWallTime() ); 
}

int  SolverSlepchenko2Buffer::
getNumberOfTimeSteps()
{
  return( GenericSolver::getNumberOfTimeSteps() );
}

