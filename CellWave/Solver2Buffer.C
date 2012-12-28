#include "CellWave.h"
#include "Solver2Buffer.h"
#include "Reaction2Buffer.h"
#include "MappedGridOperators.h"

using namespace CellWave;

bool  Solver2Buffer::
readParameterFile( CellWave::ParameterReader &param)
{
  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  //..IP3 top flux parameters
  int itemp;
  param.get( "use ip3 influx boundary", itemp, 0 ); useIP3InfluxBC= itemp;
  param.get( "ip3 influx grid number", ip3InfluxGridNumber, -1);
  param.get( "xInflux_p", xInflux_p, 0. );
  param.get( "yInflux_p", yInflux_p, 0. );
  param.get( "zInflux_p", zInflux_p, 0. );
  param.get( "influxRadius_p", influxRadius_p, 3.);
  param.get( "influxRate_p",   influxRate_p,   10.);

  if (useIP3InfluxBC) {
    DPrintf( DebugSolver,"--Solver2Buffer: using IP3 influx--\n");
    DPrintf( DebugSolver,"   ip3 influx grid number  = %d\n", ip3InfluxGridNumber );
    DPrintf( DebugSolver,"   xInflux_p               = %16.8e\n", xInflux_p );
    DPrintf( DebugSolver,"   yInflux_p               = %16.8e\n", yInflux_p );
    DPrintf( DebugSolver,"   zInflux_p               = %16.8e\n", zInflux_p );
    DPrintf( DebugSolver,"   influxRadius_p          = %16.8e\n", influxRadius_p);
    DPrintf( DebugSolver,"   influxRate_p            = %16.8e\n",   influxRate_p);
  }
  else {
    DPrintf( DebugSolver,"--Solver2Buffer: NOT using IP3 influx ** NO IP3 INFLUX --n");
  }
  bool ok = true;
  return ok;
}

#if 0
void Solver2Buffer::
setup()
{
  CellWave::ReactionFactory rxnFactory;
  Solver2Buffer::pChem  
    = (Reaction2Buffer*) rxnFactory.getReaction("2Buffer"); //FIXME: static_cast<>
  GenericSolver::genericSetup( Solver2Buffer::pChem,  );
  GenericSolver::updateModel();
}
#endif

void Solver2Buffer::
setup(CellWave::ParameterReader &param)
{
  CellWave::ReactionFactory rxnFactory; // here -- should get the name from the param file!!
  Solver2Buffer::pChem  
    = (Reaction2Buffer*) rxnFactory.getReaction("2Buffer");
  GenericSolver::genericSetup( Solver2Buffer::pChem, param );
  this->readParameterFile( param );

  useDebugShow=false; // in the solver routine

}

//void Solver2Buffer::
//finish() -- in GenericSolver (base class)

void Solver2Buffer::
initialData()
{
  GenericSolver::genericInitialData();
}

void Solver2Buffer::
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
  
  // ..start timestepping
  GenericSolver::wallTime=0.;
  GenericSolver::iWallTimeSteps=0;
  GenericSolver::totalWallTime=0.;

  char buffer[80];

  const int &saveEvery = data.saveEveryNthFrame;
  const int &logEvery  = data.logEveryNthFrame;

  Ogshow debugshow;
  if ( useDebugShow ) {
    debugshow.open("OUT/debugShow.show");
    debugshow.saveGeneralComment("CellWave: RHS info");
    debugshow.setFlushFrequency( data.flushFrequency );
  }

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

    //probes.collectData( ktime, tcomp, q);
    collectProbeData( ktime, tcomp, q );

    const bool saveThis= (ktime % saveEvery == 0);
    const bool logThis = (ktime % logEvery  == 0);

    if( logThis  ) 
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
	      "..step %3i (t=%8.2e):c=[%6.3g,%6.3g] p=[%6.3g,%6.3g] h=[%6.3g,%6.3g] ",
	     ktime, tcomp, cmin,cmax, pmin,pmax, hmin,hmax );
      DPrintf(DebugSolver,
	      "b1=[%6.3g,%6.3g] b2=[%6.3g,%6.3g]\n",
	      b1min,b1max,b2min,b2max);

      DPrintf(PRINT,
	      "..step %3i (t=%8.2e):c=[%6.3g,%6.3g] p=[%6.3g,%6.3g] h=[%6.3g,%6.3g] ",
	     ktime, tcomp, cmin,cmax, pmin,pmax, hmin,hmax );
      DPrintf(PRINT,
	      "b1=[%8.3g,%8.3g] b2=[%8.3g,%8.3g]\n",
	      b1min,b1max,b2min,b2max);

      wallTime =0.;     // reset
      iWallTimeSteps=0;
    }

    if( saveThis ) {
      show.startFrame();                                              
      show.saveComment(0,sPrintF(buffer,"Step %i",ktime));
      show.saveComment(1,sPrintF(buffer,"  t=%e ",tcomp));           
      show.saveSolution( q );                                        
      show.endFrame();
    }

    //
    //..Perform one time step, update each grid
    //
    f = 0.;
    for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
      //Index I1, I2,I3;
      //getIndex( cg[igrid].gridIndexRange(), I1,I2,I3);
      getIndex( cg[igrid].indexRange(), I1,I2,I3);
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

      //  (1a) Chemistry: (Laplace terms imposed in step (2))
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
	    = maskArray(If1,If2,If3)*rhsArray(If1,If2,If3,ispecies); //no rxn in nucleus, also for p
	}
      }

      //  (1b) IP3 Source terms for 2D model
      if (  (nd==2)  && (ip3InfluxGridNumber==igrid) && useIP3InfluxBC ) {
	addIP3VolumeSources( igrid, xArray, qArray, rhsArray );
      }
      
      //  (2) Diffusion
      real dp= chem.getDiffusionCoefficient( int(pc) );
      real dc= chem.getDiffusionCoefficient( int(cc) );
      real db2=chem.getDiffusionCoefficient( int(b2c));
      realArray tempLaplace(If1,If2,If3);

      //   .... IP3 diffusion --> update pnext
      bool alwaysDiffusion=true;
      if( alwaysDiffusion || chem.isDiffusive( int(pc) )) {
	if(ktime == 0) {
	  DPrintf(DebugSolver,"diffuse IP3, nu=%8.3g\n", dp);
	}
	tempLaplace =0.;
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,pc);
	rhsArray(If1,If2,If3,pc) += dp*tempLaplace(If1,If2,If3);
      }

      //    .... Ca2+ diffusion:
      if( alwaysDiffusion || chem.isDiffusive( int(cc) )) {
	if(ktime == 0) {
	  DPrintf(DebugSolver,"diffuse Ca, nu=%8.3g\n", dc);
	}
	tempLaplace =0.;
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,cc);
	rhsArray(If1,If2,If3,cc) += dc*tempLaplace(If1,If2,If3);   //FIXME
      }

      //    .... h diffusion:
      //         --none--

      //    .... b1 diffusion:
      //         --none--

      //    .... b2 diffusion:
      if( alwaysDiffusion || chem.isDiffusive( int(b2c) )) {
	if(ktime == 0) {
	  DPrintf(DebugSolver,"diffuse Buffer2, nu=%8.3g\n", db2);
	}
	tempLaplace = 0.;
	operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				     qArray, tempLaplace, I1,I2,I3,b2c);
	rhsArray(If1,If2,If3,b2c) += db2*tempLaplace(If1,If2,If3);   //FIXME
      }
      
      //  (3) Update solution on next time (Adams-Bashforth, 2nd order, fully explicit) //FIXME -- adaptive timestepping
      if( ktime == 0 )  rhsArray_p = rhsArray;
	
      //..AB2 step
      if(ktime == 0) DPrintf(DebugSolver,"dt=%8.3g\n",dt);
      for( int ic=0 ; ic<chem.getNumberOfSpecies(); ++ic ) { //all species
	qArray(I1,I2,I3,ic) =  qArray(I1,I2,I3,ic) 
	  + .5*dt*( 3.*rhsArray(I1,I2,I3,ic) - rhsArray_p(I1,I2,I3,ic));
      }

      // //  (3x) Update previous time levels       //later, use rhsArray_p fo interp
      // rhsArray_p = rhsArray;

      //  (3xx) DEBUG: copy RHS into solution array
      //qArray(I1,I2,I3, cc+3) =  rhsArray(I1,I2,I3,cc);
      // qArray(I1,I2,I3, pc+3) =  rhsArray(I1,I2,I3,pc);
      //qArray(I1,I2,I3, hc+3) =  rhsArray(I1,I2,I3,hc);

    }; //end for igrid

    //  (4) Set bc's: internal (interpolate), and physical

    q.periodicUpdate();
    q.interpolate();

    // (4a): apply a neumann BC on all boundaries -- h does not get a bc
    //cconcentration.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    //pconcentration.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    //q.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);

    Range comps(0, chem.getNumberOfSpecies()-1);
    q.applyBoundaryCondition(comps,BCTypes::extrapolate,BCTypes::allBoundaries,0.);
    q.periodicUpdate();
    q.finishBoundaryConditions(); //..not finished yet, but sets the corners correctly

    // (4b): apply flux BC

    const int numSpecies=chem.getNumberOfSpecies();
    fp.dataCopy(q);

    //..impose flux & noflux bcs:
    //     first impose noflux on internal flux bdries
    //     then apply flux/jump bcs on the internal flux bdries where
    //       flux coeff>0
    //
    #if 0
    for ( int ic=0; ic< chem.getNumberOfSpecies(); ++ic ) {
      const int fluxBC_ID=fluxBCData. getFluxBoundaryID();
      q.applyBoundaryCondition(ic,BCTypes::neumann,fluxBC_ID,0.);     
      if ( chem.hasFluxBC( ic ) ) {
	double fluxCoeff= chem.getFluxBCCoefficient( ic );
	setFluxBCCoefficient( fluxCoeff );
	applyFluxBoundaryConditions( q, ic );
      }
      applyNoFluxBoundaryConditions( q, ic);
    }
    #endif
    q.applyBoundaryCondition(comps,BCTypes::neumann,
                             BCTypes::allBoundaries,0.);


    // (4c): apply top inflow BC on IP3
    if ( useIP3InfluxBC ) {
      for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
	if ( igrid == getIP3InfluxBoundaryGridID())  {	  
	  MappedGrid & mg = cg[igrid];
	  MappedGridOperators &mgop = operators[igrid];
	  const double pi = 4.0*atan(1.0);

	  int side=-1, axis=-1; // find the top
	  const int Start=0, End=1; // limits for 'side'
	  bool foundTheTop=false;
	  for( axis=0; axis< mg.numberOfDimensions(); ++axis) {
	    for( side=Start; side<=End; ++side ) {
	      const int ibc= mg.boundaryCondition()(side,axis);
	      //DPrintf(DetailedDebugPrint, "> looking for top, grid %i, side %d, axis %d\n",
	      //igrid, side, axis);
	      if( ibc == ip3InfluxBoundaryID ) {
		foundTheTop=true;
		//DPrintf(DetailedDebugPrint, ">>> FOUND THE TOP. grd %i, sd %d, ax %d\n",
		//      igrid, side, axis);
		goto doneLookingForTop;
	      }
	    }// end for side
	  } //end for axis
	doneLookingForTop:
	  if( foundTheTop ) {
	    Index Ib1, Ib2, Ib3, Ig1,Ig2,Ig3;
	    MappedGrid & mg = cg[igrid];
	    DPrintf(DebugSolver, "   IP3 influx surf found: grid %d, side %d, axis %d, bc=%d\n",
		    igrid, side, axis, mg.boundaryCondition()(side,axis) ); fflush(0);
	    realMappedGridFunction &qmg = q[igrid];  
	    getBoundaryIndex(mg.gridIndexRange(), side, axis, Ib1,Ib2,Ib3);
	    getGhostIndex(   mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	    realArray & qArray     = q[igrid];
	    realArray & rhsArray = f[igrid];
	    realArray & rhsArray_p = fp[igrid];
	    realArray & xArray = mg.vertex();  // array of vertices
	    
	    //CREATE TEMPORARY STORAGE FOR INFLUX DATA
	    realArray ip3Flux(Ib1,Ib2,Ib3);
#define XC xArray(Ib1,Ib2,Ib3,axis1)
#define YC xArray(Ib1,Ib2,Ib3,axis2)
	    
	    ip3Flux =  influxRate_p*exp( -( pow( XC-xInflux_p,2.) 
	          + pow(YC-yInflux_p,2.))/(4.*influxRadius_p) 
	           )
	                     /( 4.0*pi*influxRadius_p );
	    //ip3Flux =  influxRate_p*exp( -( pow( XC-xInflux_p,2.))/(4.*influxRadius_p)) 
	    //             /( 4.0*pi*influxRadius_p );
	    //ip3Flux =  XC;
	    
	    double fluxMin= min(abs(ip3Flux)), fluxMax= max(abs(ip3Flux));
	    const int thisBCIndex=BCTypes::boundary1+side+2*axis;
	    //qmg.applyBoundaryCondition( pc, BCTypes::neumann, thisBCIndex, ip3Flux);
	    //qmg.applyBoundaryCondition( pc, BCTypes::dirichlet, thisBCIndex, ip3Flux);
	    qmg(Ib1,Ib2,Ib3,pc) = ip3Flux;
	    //mgop.applyBoundaryCondition( qmg, pc, BCTypes::dirichlet,thisBCIndex,ip3Flux);//DOES NOT WORK!!!
	    double resMin= min(abs( qmg(Ib1,Ib2,Ib3,pc))), resMax= max(abs(qmg(Ib1,Ib2,Ib3,pc)));
	    DPrintf(DebugSolver, "    setting IP3 influx: grid %d, fluxMin=%8.4e, fluxMax=%8.4e, pmin=%8.4e, pmax=%8.4e\n",
		    igrid, fluxMin, fluxMax, resMin, resMax);
	    DPrintf(DebugSolver, "    the boundary is side %d, axis %d, boundary1=%d, boundary5=%d, ibc=%d\n",
		    side, axis, BCTypes::boundary1, BCTypes::boundary5, thisBCIndex);
	  }
	  else {
	    DPrintf(DebugSolver, "   IP3 influx error. Top surf NOT FOUND, grid %\n",
		    igrid );
	  }
	  //endif foundTheTop
	}
      } //end for igrid
    } //end if useIP3InfluxBC

    q.finishBoundaryConditions(); 

    // (5): Update previous time levels
    fp.dataCopy( f );

    //debug: save rhs
    if ( useDebugShow && (ktime % saveEvery ==0 )) {
      debugshow.startFrame(); 
      debugshow.saveComment(0,sPrintF(buffer,"rhs(c,p,h,b1,b2) step %i",ktime));
      debugshow.saveComment(1,sPrintF(buffer,"  t=%e ",tcomp));           
      debugshow.saveSolution( f );                                        
      debugshow.endFrame();
    }

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


real Solver2Buffer::
getTime()
{
  return( GenericSolver::getTime() );
}

real Solver2Buffer:: 
getTotalWallTime()
{
  return( GenericSolver::getTotalWallTime() ); 
}

int  Solver2Buffer::
getNumberOfTimeSteps()
{
  return( GenericSolver::getNumberOfTimeSteps() );
}

