#include "Overture.h"
#include <float.h>

#include "FluxBC.h"

//

#include "Overture.h"
#include "CompositeGridOperators.h"
#include "interpolateFluxBoundary.h"

#include "Display.h"

#define FLUXBC_ForBoundary(side,axis) \
    for( axis=0; axis<mg.numberOfDimensions(); axis++ ) \
      for( side=0; side<=1; side++ )   


FluxBC::FluxBC( ) 
{ 
  //.. only base class constructor

  //debug=31; //**FIXME -- lots of debug output, please set debug=1
  debug=1;
}

FluxBC::
~FluxBC( ) 
{
  //.. only base class destructor
} 

void FluxBC::
updateToMatchGrid( CompositeGrid &cg, int bcID /* =0*/ )    //TODO: precompute the interp. stencils HERE!!
{ 
  this->setCompositeGrid( cg ); 
}

void FluxBC::
setCompositeGrid( CompositeGrid &cg ) 
{ 
  pCg = &cg; 
  if( pCg == NULL ) printf("ERROR FluxBC::setCompositeGrid has a problem -- pCg == NULL\n");
  //interpolatePoints ....
}

void  FluxBC::
setupInterpolation(const int numberOfOutputComponents /*=1 */)
  // sets up interpolation of bc id = idFluxBoundary
{
  //printf("FluxBC::applyBoundaryCondition called...\n");

  assert( pCg     != NULL ); CompositeGrid &cg          = *pCg;
  assert( pInterp != NULL ); Interpolant   &interp      = *pInterp;
  assert( pOp     != NULL ); CompositeGridOperators &op = *pOp;

  printf("***FluxBC::setupInterpolation: numComponents %d\n",
	 numberOfOutputComponents);

  //..interpolate in a two step process: 
  //   Collect data:
  //    *Step 1: setupInterpolation:
  //    (1) allocate the start/stop index arrays [ngrids][naxis][nsides]
  //    (2) loop through grids/edges; compute # points & collect xyz
  //    (3) form long xyzInterp array & set pointers to it

  //    *Step 2: applyBoundaryCondition -- see the next subroutine
  
  timerInterpSetupCode = getCPU();

  const int nGrids= cg.numberOfComponentGrids(); //shorthand
  const int nAxes = cg.numberOfDimensions();
  const int nSides= 2;

  // (1) allocate 
  edgeStart.redim(nGrids,nAxes,nSides);
  edgeEnd.redim(nGrids,nAxes,nSides);
  isInterpolatedEdge.redim(nGrids,nAxes,nSides);
  edgeCoords.redim(nGrids,nAxes,nSides);
  int totalNumberOfInterpolationPoints=0; //local counter, global count in xyzInterpolate

#define checkdim(X,n0,n1,n2) ( (n0<=X.size(0)) && (n1<=X.size(1)) && (n2<=X.size(2)))


  // (2) collect local xyz edge coordinates
  for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
    MappedGrid & mg = cg[ig];
    
    Index Ib1,Ib2,Ib3;
    Index Ig1,Ig2,Ig3;
    Index All;
    
    for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) {
      for( int side=0; side<=1; side++ ) {
	edgeStart(ig,axis,side) = -1;
	edgeEnd(ig,axis,side)   = -1;

	isInterpolatedEdge(ig, axis,side) = false;
	if( mg.boundaryCondition()(side,axis) == idFluxBoundary  ) {
	  isInterpolatedEdge(ig,axis,side) = true;
	  getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	  //getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	  //const realArray &zz=mg.vertex();
	  int zaxis        =axis3;
	  if (nAxes==2) zaxis=axis2; // to make sure we don't seg.fault
	  const realArray &xb  = mg.vertex()(Ib1,Ib2,Ib3, axis1);
	  const realArray &yb  = mg.vertex()(Ib1,Ib2,Ib3, axis2);
	  const realArray &zb  = mg.vertex()(Ib1,Ib2,Ib3, zaxis);
	  
	  int axisLength[3]={xb.getLength(axis1),xb.getLength(axis2),xb.getLength(axis3)};
	  int nInterpPoints= axisLength[axis1]*axisLength[axis2];
	  if (cg.numberOfDimensions()==3) nInterpPoints = nInterpPoints* axisLength[axis3];
	  
	  //realArray xyInt( Ib1,Ib2,Ib3, nAxes );
	  edgeCoords(ig,axis,side).redim( Ib1,Ib2,Ib3, nAxes );
	  realArray &xyz = edgeCoords(ig,axis,side); //shorthand
	  xyz(All, All, All, axis1) = xb;
	  xyz(All, All, All, axis2) = yb;
	  if( nAxes == 3) {
	    xyz(All, All, All, axis3) = zb;
	  }
	  xyz.reshape(nInterpPoints, nAxes);
	  totalNumberOfInterpolationPoints += nInterpPoints;

	  //printf("..On (grid %d, axis %d,side %d); ",  ig,axis,side);
	  //printf("passing in xyz array of dimensions [%d, %d], ntotal=%d\n",
	  //	 edgeCoords(ig,axis,side).getLength(0), 
	  //	 edgeCoords(ig,axis,side).getLength(1), totalNumberOfInterpolationPoints);
	} // end if ibc==fluxBoundaryID
      }// end for side
    }//end for axis
  }//end for ig

  // (3) form long xyzInterp array & set start/end indices
  xyzInterpolate.redim(totalNumberOfInterpolationPoints, nAxes );
  //const int nComponents=1; //set bcs one component at a time for now.
  valuesInterpolate.redim(totalNumberOfInterpolationPoints, 
			  numberOfOutputComponents );
  valuesInterpolate=-981.;
  ignoreGrid.redim(totalNumberOfInterpolationPoints);
  xyzInterpolate=-99;
  ignoreGrid=-1;
  int iNextStart=0;
  for (int ig=0; ig< nGrids; ++ig ) {
    for( int axis=0; axis<nAxes; ++axis ) {
      for( int side=0; side<nSides; ++side ) {
	if( isInterpolatedEdge(ig,axis,side) ) {
	  //printf("--grid %3d, axis %d, side %d:\n", ig,axis,side);
	  realArray &xyz          = edgeCoords(ig,axis,side); //shorthand
	  const int nInterpPoints = xyz.getLength(0);	  
	  const int iThisEnd      = iNextStart+nInterpPoints-1;
	  edgeStart(ig,axis,side) = iNextStart;
	  edgeEnd(ig,axis,side)   = iThisEnd;

	  Range dims(0,nAxes-1);
	  Range subset(iNextStart,iThisEnd);

#if 0  //debug
	  std::string namexyz=          "xyz              ";
	  std::string namexyzreshape=   "xyz reshaped     ";
	  std::string namexyzInterp=    "xyzInterpolate   ";
	  std::string namexyzInterpSub= "xyzInterp(sub)   ";
	  std::string subsetname=       "range subset     ";
	  std::string dimsname=         "range dims       ";
	  printBaseBound(namexyz, xyz);
	  printBaseBound(subsetname, subset);
	  printBaseBound(dimsname, dims);
#endif
	  xyz.reshape(subset,dims);

#if 0 //debug
	  printBaseBound(namexyzreshape, xyz);
	  printBaseBound(namexyzInterp, xyzInterpolate);fflush(0);
	  printBaseBound(namexyzInterpSub, xyzInterpolate(subset,dims));fflush(0);
#endif

	  xyzInterpolate(subset, dims) = xyz;  //interp. points for this edge
	  ignoreGrid(subset)           = ig;   //do not interp 'xyz' from grid =ig

	  iNextStart += nInterpPoints;
	}
      }
    }
  }
  interpolator. buildInterpolationInfo(  xyzInterpolate, cg, ignoreGrid );
  timerInterpSetupCode = getCPU() - timerInterpSetupCode;

  //..CHECK -- if some 'ignoreGrid(index)' items were not set==> skipped some = problem
  printf("Checking 'ignoreGrid'...\n");
  bool igok=true;
  for(int i=0; i<totalNumberOfInterpolationPoints; ++i ) {
    int ig=ignoreGrid(i);
    if( ig<0 ) {
      printf("..bug: ignoreGrid[%5d] = %3d\n", i, ig);
      igok=false;
    }
  }	 
  printf("--> ignoreGrid is "); 
  if( igok ) printf("OK\n");    //covered all items
  else       printf("NOT OK\n");//skipped some=problem

}


void  FluxBC::
applyBoundaryCondition( realCompositeGridFunction &q,
			int ic) 
{
  assert( pCg     != NULL ); CompositeGrid &cg          = *pCg;
  assert( pInterp != NULL ); Interpolant   &interp      = *pInterp;
  assert( pOp     != NULL ); CompositeGridOperators &op = *pOp;

  Display display; //David's display class

  realArray &qTemp = q[0];
  const int nComponents         = qTemp.getLength(3); // indices ( 0, 1, 2, --3-- ), 3=component
  const int nInterpComponents   =valuesInterpolate.getLength(1);
  assert( nComponents <= nInterpComponents );
  if(debug&4)printf("***FluxBC::applyBC: nComponents %d, nInterpComponents %d, fluxCoeff %g\n",
		    nComponents, nInterpComponents,  getFluxCoefficient());

  //printf("FluxBC::applyBoundaryCondition called...\n");

  timerInterpCode = getCPU();
  interpolator.interpolatePoints(q,valuesInterpolate);

  if(debug&8) {
    printf("-------------------DISPLAY: valuesInterpolate in FluxBC::applyBoundaryCondition %d  -------\n",
	   ic);
    display.display(valuesInterpolate,"valuesInterpolate");
  }

  for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
    MappedGrid & mg = cg[ig];
    MappedGridOperators    &opmg = op[ig];
    realMappedGridFunction &q_mg = q[ig];
    Index Ib1,Ib2,Ib3;
    Index Ig1,Ig2,Ig3;

    Index I1,I2,I3;
    getIndex( mg.dimension(), I1,I2,I3);
    
    //..flux/jump bc
    const int nAxes=cg.numberOfDimensions();
    int axis, side;
    for( axis=0; axis<mg.numberOfDimensions(); axis++ ) {
      for( side=0; side<=1; side++ ) {
	
	if( mg.boundaryCondition()(side,axis) == idFluxBoundary  ) {

	  //printf("FluxBC: grid, side, axis(%i,%i,%i) is a flux boundary (bc=%i)\n",
	  //        ig,side,axis,idFluxBoundary);
	  getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	  getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	  realArray       &qArray     = q[ig];

	  int nInterpPoints= edgeEnd(ig,axis,side) - edgeStart(ig,axis,side)+1;
	  Range dims(0,nAxes-1);
	  //Range subset(edgeStart(ig,axis,side,edgeEnd(ig, axis,side)));
	  Range subset(edgeStart(ig,axis,side),edgeEnd(ig, axis,side));

	  realArray jump(nInterpPoints);
	  jump.reshape(subset);

	  jump(subset) = valuesInterpolate(subset,ic);
	  jump.reshape(Ib1,Ib2,Ib3, Range(ic,ic));

	  if(debug&8) {
	    printf("..grid %d, side %d, axis %d, component %d, flux coeff %g-- JUMP\n", 
		   ig, side, axis, ic, getFluxCoefficient() );
	    display.display(jump,"jump: values for other grids");
	  }
	  jump = getFluxCoefficient()*( jump - qArray(Ib1,Ib2,Ib3,ic));

	  //if(debug&8  || (ic==1)) {
	  if(debug&8) {
	    printf("....Flux*[ jump ]\n");
	    display.display( jump, "flux coeff*[ other - this value]" );
	  }

	  //.. ibc chosen to set a specified boundary (side,axis)
	  // ... take jumps from unext --> impose as flux in u
	  const int ibc=BCTypes::boundary1+side+2*axis;
	  opmg.applyBoundaryCondition( q_mg, Range(ic,ic), BCTypes::neumann, ibc, jump ); 

	  if(debug&8)  display.display(jump,   "jump");
	  if(debug&16) display.display(qArray(I1,I2,I3,ic), "q function");

	} // end if ibc==fluxBoundaryID


      } 
    }
  }
  timerInterpCode = getCPU()-timerInterpCode;

}
