/// Brief description: common code for Solver classes in CellWave
///

#include <assert.h>
#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

//#include "CellWave.h"
#include "GenericReaction.h"
#include "ReactionFactory.h"
#include "Info.h"
#include "getDiffusionDT.h"

#include "GenericSolver.h"
#include "CellWave.h"

using namespace CellWave;

//.. virtual functions -- defined in derived classes
bool CellWave::GenericSolver::
readParameterFile( CellWave::ParameterReader &param)
{
  nucleusClass.readParameterFile( param );
  probes.readParameterFile( param );
  
  return true; //what's the return value here, then?
}


//.. Base class services

void GenericSolver::
genericSetup( GenericReaction *pChem_, CellWave::ParameterReader &param )
{
  //CalciumChemistry chem(cg.numberOfDimensions());  // initializes parameters automatically (by constructor)
  pChem = pChem_;
  assert( pChem != NULL );
  CellWave::GenericReaction &chem = *pChem; // shorthand
 
  // create and read in a CompositeGrid
  //CompositeGrid cg;
  //aString nameOfGrid = data.nameOfOGFile.c_str();
  int ierr=getFromADataBase( cg, data.nameOfOGFile.c_str() );

  DPrintf(DebugSolver,"GenericSolver -- ierr for getFromADataBase = %d\n",
	  ierr);
  if ( ierr == 1) { // not found
    DPrintf(BroadcastPrint,"ERROR -- grid file '%s' not found, exiting.\n",
	    data.nameOfOGFile.c_str());
    throw "error";
  }

  //cg.update();
  cg.update(MappedGrid::THEvertex | MappedGrid::THEcenter 
	    | MappedGrid::THEvertexBoundaryNormal);
  //Interpolant interpolant(cg);
  interpolant.updateToMatchGrid( cg );

  // ..create a show file
  // ..then save a general comment in the show file
  //Ogshow show( nameOfShowFile );                     
  show.open( data.nameOfShowFile.c_str() );
  show.saveGeneralComment("CellWave for Calcium wave modeling");
  show.setFlushFrequency( data.flushFrequency );   
  operators.updateToMatchGrid( cg );
  
  //..Read parameters
  bool ok=  chem.readParameterFile( param );
  if( ok ) {
    GenericSolver::readParameterFile( param );
    GenericSolver::updateModel();
  }
  else {
    DPrintf(BroadcastPrint, "ERROR -- GenericSolver, unable to read params for chem.\n");
    throw "error";
  }

  //..FluxBCs
  int noFluxBC;
  param.get("no flux bc",noFluxBC, 0); //BY DEFAULT, flux BC's are on
  if (noFluxBC) {
    turnOffFluxBoundaryConditions();
    DPrintf(BroadcastPrint,"--FluxBCs turned off");
    DPrintf(BroadcastPrint,"(GenericSolver::genericSetup)\n");
  }  
  else {
    turnOnFluxBoundaryConditions();
    DPrintf(BroadcastPrint,"--FluxBCs turned **on**");
    DPrintf(BroadcastPrint,"(GenericSolver::genericSetup)\n");
  }
  setupFluxBoundaryConditions( cg, interpolant, operators );

  //..Probes
  if ( probes.isOK() ) {
    probes.outputHeader(false);
    probes.initializeProbeValues( chem.getNumberOfSpecies() );
    probes.openOutputFile( "$Id: GenericSolver.C,v 1.23 2004/02/10 03:53:57 pfast Exp $" );
  }

}

void GenericSolver::
genericFinish()
{
  probes.closeOutputFile();
  show.close();
}

void GenericSolver::
collectProbeData( int ktime, double tcomp, realCompositeGridFunction &q)
{
  probes.collectData( ktime, tcomp, q );
}

void GenericSolver::
updateModel()
{
  assert( pChem != NULL );
  CellWave::GenericReaction &chem = *pChem; // shorthand

  //....create grid functions
  Range all;
  //q.updateToMatchGrid(cg,all,all,all,2* chem.getNumberOfSpecies()); //debug
  q.updateToMatchGrid(cg,all,all,all,chem.getNumberOfSpecies());
  q.setOperators( operators );
  q.setName("all");              

  //const int numberOfIP3Species=1; // ip3 is one species, no ip3 buffers etc
  //ip3Next.updateToMatchGrid(cg,all,all,all, numberOfIP3Species);
  //ip3Next.setOperators( operators );
  //ip3Next.setName("IP3 at time t(n+1)");

  nucleusGF.updateToMatchGrid( cg );
  nucleusClass.updateToMatchGrid( cg, nucleusGF);
  nucleusClass.evaluateGridFunction();

  int numSpecies=chem.getNumberOfSpecies();
  std::string compName =chem.getTitle();
  printf("GenericSolver::updateModel-->%s\n", compName.c_str());
  q.setName( compName.c_str() );
  for (int j=0; j<numSpecies; ++j)  {
    compName = chem.getLongComponentName( j );
    q.setName( compName.c_str(), j);
  }

  f. updateToMatchGrid(cg,all,all,all,chem.getNumberOfSpecies()); 
  fp.updateToMatchGrid(cg,all,all,all,chem.getNumberOfSpecies());
  fp.setOperators( operators );
  
};

void GenericSolver::
genericInitialData()
{
  assert( pChem != NULL );
  CellWave::GenericReaction &chem = *pChem; // shorthand
      
  //..set initial data!! --> allow hotspot in p,c,h (xc, yc, width)=exponential
  q= 0.;
  
  for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
    //Index I1, I2,I3;
    //getIndex( cg[igrid].indexRange(), I1,I2,I3);

    MappedGrid & mg = cg[igrid];
    realArray & qArray = q[igrid];
    //realArray & dudtg = dudt[igrid];
    realArray & xArray = mg.vertex();  // array of vertices

    const IntegerArray & d = mg.dimension();
    const IntegerArray & gir= mg.gridIndexRange();
    const int nd   =cg.numberOfDimensions();
    const int ncomp= chem.getNumberOfSpecies();
    double tcomp=0.; //dummy time here

    //..explicit loop here: arguments are the info for a Fortran style loop
    chem.callInitialDataLoop(tcomp, 
			     nd, ncomp,
			     d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2), 
			     gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2), 
			     *xArray.getDataPointer(),*qArray.getDataPointer() );
 
  }; //end for igrid
  q.periodicUpdate();
  q.interpolate();
}

void GenericSolver::
outputStepStart(int ktime_, double tcomp_) 
{
  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;
  
  //DPrintf(PRINT             ,"--------- Step #%d, tcomp= %g -----------------\n",
  //	      ktime_, tcomp_ );
  DPrintf(DebugPrint        ,"--------- Step #%d, tcomp= %g -----------------\n",
	  ktime_, tcomp_ );
  DPrintf(DetailedDebugPrint,"--------- Step #%d, tcomp= %g -----------------\n",
	  ktime_, tcomp_ );
  DPrintf(DebugSolver       ,"--------- Step #%d, tcomp= %g -----------------\n",
	  ktime_, tcomp_ );
  DPrintf(DebugReaction     ,"--------- Step #%d, tcomp= %g -----------------\n",
	  ktime_, tcomp_ );
}

void GenericSolver::
printReactionParameters( int iPrintChannel )
{
  if ( pChem!= NULL ) pChem->printParameters( iPrintChannel );
}

//
// ...FLUX BC CODE
//

void GenericSolver::
setupFluxBoundaryConditions( CompositeGrid          &cg,
			     Interpolant            &interpolant,
			     CompositeGridOperators &operators)
{
  assert( pChem!= NULL );
  if ( fluxBCData.useFluxBC ) { 
    fluxBCData.bc.updateToMatchGrid( cg, fluxBCData.fluxBoundaryID);
    fluxBCData.bc.setInterpolant( interpolant );
    fluxBCData.bc.setOperators(   operators );
    fluxBCData.bc.setupInterpolation( pChem->getNumberOfSpecies() ); //precomputes the stencils
  } 
  else {
    DPrintf(CellWave::DebugSolver,"**Warning(GenericSolver):");
    DPrintf(CellWave::DebugSolver,"    FluxBC not setup(no flux bc=true)**\n");
  }
}

void GenericSolver::
applyFluxBoundaryConditions( realCompositeGridFunction &u,
			     int iSolutionComponent)
{
  if ( fluxBCData.useFluxBC ) { 
    fluxBCData.bc.applyBoundaryCondition( u, iSolutionComponent);
  } 
  else {
    DPrintf(CellWave::DebugSolver,"**Warning(GenericSolver):");
    DPrintf(CellWave::DebugSolver,"    FluxBC not used(no flux bc=true)**\n"); 
  }
}

void GenericSolver::
applyNoFluxBoundaryConditions(realCompositeGridFunction &u,
			      int iSolutionComponent)
{
  const int noFluxID = fluxBCData.getPhysicalBoundaryID(); 
  u.applyBoundaryCondition( iSolutionComponent, BCTypes::neumann, 
			    noFluxID,  0. );
}

void GenericSolver::
setFluxBCCoefficient( double newFlux ) 
{
  fluxBCData.bc.setFluxCoefficient( newFlux );
}
