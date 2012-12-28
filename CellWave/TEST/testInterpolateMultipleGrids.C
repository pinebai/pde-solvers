//
// ..test interpolating flux boundaries on multiple component grids
//

#include "Overture.h"
#include <stdio.h>
#include <string>
#include "InterpolatePoints.h"

#include "ArraySimple.h"

void printBaseBound( std::string name, doubleArray X)
{
  int xbase[2];
  int xbound[2];

  for( int i=0; i<=1; ++i) {
    xbase[i]  = X.getBase(i);
    xbound[i] = X.getBound(i);
  }

  printf("array '%s' bounds: [%d:%d, %d:%d]\n",  name.c_str(),
	 xbase[0],xbound[0], xbase[1], xbound[1]);
};

void printBaseBound( std::string name, IntegerArray X)
{
  int xbase[2];
  int xbound[2];

  for( int i=0; i<=1; ++i) {
    xbase[i]  = X.getBase(i);
    xbound[i] = X.getBound(i);
  }

  printf("array '%s' bounds: [%d:%d, %d:%d]\n",  name.c_str(),
	 xbase[0],xbound[0], xbase[1], xbound[1]);
};

void printBaseBound( std::string name, Range X)
{
  int xbase[2];
  int xbound[2];

  xbase[0]  = X.getBase();
  xbound[0] = X.getBound();
  
  printf("range '%s' bounds: [%d:%d]       \n",  name.c_str(),
	 xbase[0],xbound[0]);
};



int main(int argc, char **argv)
{
  Overture::start( argc, argv );

  CompositeGrid cg;
  std::string gridName="Grids/interp20umMatching.hdf";
  //int ignoreGrid=-1;
  
  if( argc>1 ) gridName = argv[1];
  else { 
    printf("usage: %s <grid name.hdf> \nUsing default grid %s\n",
	   argv[0], gridName.c_str());
  }
  //if( argc>2 ) npoints=atoi(argv[2]);
  //if( argc>3 ) ignoreGrid=atoi(argv[3]);

  getFromADataBase( cg, gridName.c_str() );
  cg.update();
  int ndim   = cg.numberOfDimensions();

  realCompositeGridFunction f(cg);
  f=1.;

  printf("..Interpolating from %d dimensional grid '%s'\n",
	 ndim, gridName.c_str());

  const int idFluxBoundary=2;

  //..interpolate in a two step process: 
  //   Collect data:
  //    * 
  //    (1) allocate the start/stop index arrays [ngrids][naxis][nsides]
  //    (2) loop through grids/edges; compute # points & collect xyz
  //    (3) form long xyzInterp array & set pointers to it

  const int nGrids= cg.numberOfComponentGrids();
  const int nAxes = cg.numberOfDimensions();
  const int nSides= 2;

  // (1) allocate
  //ArraySimple<int>  edgeStart(nGrids,nAxes,nSides),edgeEnd(nGrids,nAxes,nSides);
  //ArraySimple<bool> isInterpolatedEdge(nGrids,nAxes,nSides);
  //ArraySimple<realArray> edgeCoords(nGrids,nAxes,nSides); //collect the xyz here
  ArraySimple<int>   edgeStart, edgeEnd;
  ArraySimple<bool>  isInterpolatedEdge;
  ArraySimple<realArray> edgeCoords;

  edgeStart.redim(nGrids,nAxes,nSides);
  edgeEnd.redim(nGrids,nAxes,nSides);
  isInterpolatedEdge.redim(nGrids,nAxes,nSides);
  edgeCoords.redim(nGrids,nAxes,nSides);

#define checkdim(X,n0,n1,n2) ( (n0<=X.size(0)) && (n1<=X.size(1)) && (n2<=X.size(2)))

  realArray     xyzInterpolate;
  IntegerArray  ignoreGrid;
  int totalNumberOfInterpolationPoints=0;

  // (2) collect xyz
  for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
    MappedGrid & mg = cg[ig];
    //MappedGridOperators      &opmg = op[ig];
    //realMappedGridFunction &qto_mg = qto[ig];
    
    Index Ib1,Ib2,Ib3;
    Index Ig1,Ig2,Ig3;
    Index All;
    
    //..flux/jump bc
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

	  printf("..On (grid %d, axis %d,side %d); passing in xyz array of dimensions [%d, %d], ntotal=%d\n",
		 ig,axis,side,
		 edgeCoords(ig,axis,side).getLength(0), 
		 edgeCoords(ig,axis,side).getLength(1), totalNumberOfInterpolationPoints);
	}
      }
    }
  }    
  // (3) form long xyzInterp array & set start/end indices
  xyzInterpolate.redim(totalNumberOfInterpolationPoints, nAxes );
  ignoreGrid.redim(totalNumberOfInterpolationPoints);
  xyzInterpolate=-99;
  ignoreGrid=-1;
  int iNextStart=0;
  for (int ig=0; ig< nGrids; ++ig ) {
    for( int axis=0; axis<nAxes; ++axis ) {
      for( int side=0; side<nSides; ++side ) {
	if( isInterpolatedEdge(ig,axis,side) ) {
	  printf("--grid %3d, axis %d, side %d:\n", ig,axis,side);
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

  //..CHECK
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
  if( igok ) printf("OK\n");
  else       printf("NOT OK\n");

  // (4) Build interpolation array
  InterpolatePoints interp;
  printf("..before 'buildInterpolationInfo\n");fflush(0);
  interp.buildInterpolationInfo(xyzInterpolate, cg, ignoreGrid );
  printf("..after  'buildInterpolationInfo\n");fflush(0);
  
  //..output the interp info
  IntegerArray indexValues, interpoleeGrid;
  interp.getInterpolationInfo(cg, indexValues, interpoleeGrid);
  //indexValues.display("IndexValues");
  //interpoleeGrid.display("InterpoleeGrid");

  // (5) try interpolating
  realArray vals(totalNumberOfInterpolationPoints);
  vals=0;
  interp.interpolatePoints(f, vals);

  // ..output
  printf("..Interpolation stencils (i,j,k) and the interpolee grid:\n");
  for( int i=0; i<totalNumberOfInterpolationPoints; ++i ) {
    printf("  (%4d,%4d,%4d), grid %3d, data = %10e\n", 
	   indexValues(i,0), indexValues(i,1), indexValues(i,2),
	   interpoleeGrid(i), vals(i)) ;
  }



  Overture::finish();
}
