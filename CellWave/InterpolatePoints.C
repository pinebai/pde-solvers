#include "Overture.h"
#include "display.h"
#include "InterpolatePoints.h"


// The macro MODR shifts a point back into the main periodic region
#define NRM(axis)  ( indexRange(End,axis)-indexRange(Start,axis)+1 )
#define MODR(i,axis)  ( \
  ( (i-indexRange(Start,axis)+NRM(axis)) % NRM(axis)) \
      +indexRange(Start,axis) \
                           )

static int localDebug=0;   // 1+2+4+8;

InterpolatePoints::
InterpolatePoints()
{
  indirection=NULL;
  interpolationLocation=NULL;
  interpolationLocationPlus=NULL;
  interpolationCoordinates=NULL;

}

InterpolatePoints::
~InterpolatePoints()
{
  delete [] indirection;
  delete [] interpolationLocation;
  delete [] interpolationLocationPlus;
  delete [] interpolationCoordinates;

}

const IntegerArray & InterpolatePoints::
getStatus() const
// ==================================================================================
// /Description:
//    Return the status array for the last interpolation. The values in status are from the
// InterpolationStatusEnum.
// ==================================================================================
{
  return status;
}


int InterpolatePoints::
getInterpolationInfo(CompositeGrid & cg, IntegerArray & indexValues,IntegerArray & interpoleeGrid) const
// ==================================================================================
// /Description:
//   Return the index values and interpoleeGrid for the last interpolation.
// ==================================================================================
{
  const int numberOfDimensions = cg.numberOfDimensions();
  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  int totalNumberOfInterpolationPoints=sum(numberOfInterpolationPoints);
  
  if( totalNumberOfInterpolationPoints==0 )
  {
    indexValues.redim(0);
    interpoleeGrid.redim(0);
    return 1;
  }
  
  assert( indirection!=NULL );

  indexValues.redim(totalNumberOfInterpolationPoints,3);
  indexValues=0;
  interpoleeGrid.redim(totalNumberOfInterpolationPoints);
  interpoleeGrid=0;
  
  // We must check the highest priority grid first since this was the order the points were generated.
  // A point may be extrpolated on a higher priodirt grid but then interpolated on a lower priority grid.
  int grid;
  for( grid=numberOfComponentGrids-1; grid>=0; grid-- )  // check highest priority grid first
  {
    const int num=numberOfInterpolationPoints(grid);
    if( num>0 )
    {
      IntegerArray & ia = indirection[grid];
      IntegerArray & ip = interpolationLocation[grid];
      const int *iap = ia.Array_Descriptor.Array_View_Pointer0;
#define IA(i0) iap[i0]
      const int *ipp = ip.Array_Descriptor.Array_View_Pointer1;
      const int ipDim0=ip.getRawDataSize(0);
#define IP(i0,i1) ipp[i0+ipDim0*(i1)]

      for( int i=0; i<num; i++ )
      {
	for( int axis=0; axis<numberOfDimensions; axis++ )
	  indexValues(IA(i),axis)=IP(i,axis);
	interpoleeGrid(IA(i))=grid;
      }
      
    }
  }
  
  return 0;

}

#undef IP
#undef IA


int InterpolatePoints::
buildInterpolationInfo(const realArray & positionToInterpolate, CompositeGrid & cg, 
		       int ignoreGrid)
// =================================================================================================
// /Description:
//    Build the interpolation location arrays that can be used to interpolate at some specified points.
//    Will not interpolate from component grid {\tt ignoreGrid} if specified.
//
// /Return value: if negative then the absolute value of the return value is the number of points 
//   not assigned.
// =================================================================================================
{
  int numberOfPointsToInterpolate=positionToInterpolate.getLength(0);
  IntegerArray ignoreGridArray(numberOfPointsToInterpolate);
  ignoreGridArray= ignoreGrid;
  return( buildInterpolationInfo( positionToInterpolate, cg, ignoreGridArray ));
}

int InterpolatePoints::
buildInterpolationInfo(const realArray & positionToInterpolate, 
		       CompositeGrid & cg, 
		       IntegerArray &ignoreGrid/*=nullIntArray()*/ )
// =================================================================================================
// /Description:
//    Build the interpolation location arrays that can be used to interpolate at some specified points.
//    Will not interpolate from component grid {\tt ignoreGrid} if specified.
//
// /Return value: if negative then the absolute value of the return value is the number of points 
//   not assigned.
// =================================================================================================
{
  const int debug=0;

  cg.update(MappedGrid::THEboundingBox);

  const int numberOfDimensions = cg.numberOfDimensions();

  const int numberOfComponentGrids = cg.numberOfComponentGrids();
  Range Axes=numberOfDimensions;

  int grid, axis;
  for( grid=0; grid<numberOfComponentGrids; grid++ )
  {
    if( !cg[grid].isAllVertexCentered() && !cg[grid].isAllCellCentered() )
    {
      cout << "interpolatePoints:ERROR: grids must be either vertex or cell centered, no mongrels! \n";
      return 1;
    }      
  }

  int numberOfPointsToInterpolate=positionToInterpolate.getLength(0);
  bool interpolateFromAllGrids=true; //if =false, we will ignore grid 'ignoreGrid(i)' for point i
  bool isIgnoreArrayCompatible= (ignoreGrid.getBase(0)  <= positionToInterpolate.getBase( 0 ))
                             && (ignoreGrid.getBound(0) >= positionToInterpolate.getBound( 0 ));
  if( isIgnoreArrayCompatible ) interpolateFromAllGrids=false; //if ignore is ok, we'll use it

  const real epsi=1.e-3;
  int extrap,pointWasExtrapolated;
  int returnValue=0;  // 0=ok, >0 error, <0 some points extrapolated

  delete [] indirection;
  delete [] interpolationLocation;
  delete [] interpolationLocationPlus;
  delete [] interpolationCoordinates;

  indirection = new IntegerArray [numberOfComponentGrids];
  interpolationLocation = new IntegerArray[numberOfComponentGrids];
  interpolationLocationPlus = new IntegerArray[numberOfComponentGrids];
  interpolationCoordinates = new RealArray[numberOfComponentGrids];

  numberOfInterpolationPoints.redim(numberOfComponentGrids);
  numberOfInterpolationPoints=0;
  int *numberOfInterpolationPointsp = numberOfInterpolationPoints.Array_Descriptor.Array_View_Pointer0;
#define NUMBEROFINTERPOLATIONPOINTS(i0) numberOfInterpolationPointsp[i0]
    int *pIgnoreGrid = ignoreGrid.Array_Descriptor.Array_View_Pointer0;
#define IGNOREGRID(i0) pIgnoreGrid[i0]


  // *****************************************************
  // ******** Find an interpolation stencil to use *******
  // *****************************************************

  Range R=numberOfPointsToInterpolate;
  

  IntegerArray ia0(R);
  status.redim(R);
  int *statusp = status.Array_Descriptor.Array_View_Pointer0;
#define STATUS(i0) statusp[i0]
    int *ia0p = ia0.Array_Descriptor.Array_View_Pointer0;
#define IA0(i0) ia0p[i0]

  status=notInterpolated;

  const realArray & x =positionToInterpolate;
  real bb[2][3];

  grid=numberOfComponentGrids-1;   // check this grid first
  for( grid=numberOfComponentGrids-1; grid>=0; grid-- )  // check highest priority grid first
  { 
    MappedGrid & mg = cg[grid];
    Mapping & mapping = mg.mapping().getMapping();
    const intArray & mask = mg.mask();
    
    const RealArray & gridSpacing = mg.gridSpacing();
    const IntegerArray & indexRange = mg.indexRange();
    const IntegerArray & dimension  = mg.dimension();
    const IntegerArray & isPeriodic = mg.isPeriodic();
    const real shift = (bool)mg.isAllVertexCentered() ? 0. : .5; // shift position for cell centered grids

    const int *dimensionp = dimension.Array_Descriptor.Array_View_Pointer1;
    const int dimensionDim0=dimension.getRawDataSize(0);
#define DIMENSION(i0,i1) dimensionp[i0+dimensionDim0*(i1)]
    const int *indexRangep = indexRange.Array_Descriptor.Array_View_Pointer1;
    const int indexRangeDim0=indexRange.getRawDataSize(0);
#define INDEXRANGE(i0,i1) indexRangep[i0+indexRangeDim0*(i1)]
    const real *gridSpacingp = gridSpacing.Array_Descriptor.Array_View_Pointer0;
#define GRIDSPACING(i0) gridSpacingp[i0]


    // get the bounding box for this grid --- *** increase bounding box a bit ??

    const RealArray & boundingBox = mg.boundingBox();
    // boundingBox.display("Here is the boundingBox");

    // increase the size of the bounding box to allow for interp. from ghost point
    real scale=0.;
    for( axis=0; axis<numberOfDimensions; axis++ )
      scale=max(scale,boundingBox(1,axis)-boundingBox(0,axis));

    const real delta=scale*.25;
    for( axis=0; axis<numberOfDimensions; axis++ )
    {
      bb[0][axis]=boundingBox(0,axis)-delta;
      bb[1][axis]=boundingBox(1,axis)+delta;
    }
    
    const real *xp = x.Array_Descriptor.Array_View_Pointer1;
    const int xDim0=x.getRawDataSize(0);
#define X(i0,i1) xp[i0+xDim0*(i1)]

    // make a list of points inside the bounding box
    int j=0;
    //.. check notInterpolated and extrapolated, but not on 'ignoreGrid'
    //.... in case ignoreGrid is incompatible-->interpolateFromAllGrids==true,  **pf
    //     and the OR statement will be shortcircuited. In that case            **pf
    //     IGNOREGRID(i) won't get evaluated                                    **pf
#define CHECK_INTERPOLATION_POINT(i) ( (STATUS(i)!=interpolated) \
                                     &&(interpolateFromAllGrids || (grid!=IGNOREGRID(i))))
    if( numberOfDimensions==2 )
    {
      for( int i=0; i<numberOfPointsToInterpolate; i++ )
      {
	if( CHECK_INTERPOLATION_POINT( i ) )
	{
	  real x0=X(i,0);
	  real y0=X(i,1);
	  if( x0>=bb[0][0] && x0<=bb[1][0] &&
	      y0>=bb[0][1] && y0<=bb[1][1] ) 
	  {
	    IA0(j)=i;
	    j++;
	  }
	}
      }
    }
    else
    {
      for( int i=0; i<numberOfPointsToInterpolate; i++ )
      {
	if( CHECK_INTERPOLATION_POINT(i) )
	{
	  real x0=X(i,0);
	  real y0=X(i,1);
	  real z0=X(i,2);
	  if( x0>=bb[0][0] && x0<=bb[1][0] &&
	      y0>=bb[0][1] && y0<=bb[1][1] &&
	      z0>=bb[0][2] && z0<=bb[1][2] ) 
	  {
	    IA0(j)=i;
	    j++;
	  }
	}
      }
    }
#undef CHECK_INTERPOLATION_POINT
    
    if(j>0 )
    {
      // attempt to interpolate points from this grid.

      int numberToCheck=j;
      Range I=numberToCheck;
      realArray ra(I,numberOfDimensions),xa(I,numberOfDimensions);

      real *rap = ra.Array_Descriptor.Array_View_Pointer1;
      const int raDim0=ra.getRawDataSize(0);
#define RA(i0,i1) rap[i0+raDim0*(i1)]

    real *xap = xa.Array_Descriptor.Array_View_Pointer1;
    const int xaDim0=xa.getRawDataSize(0);
#define XA(i0,i1) xap[i0+xaDim0*(i1)]

      IntegerArray & ia = indirection[grid];
      ia.redim(I);

      int *iap = ia.Array_Descriptor.Array_View_Pointer0;
#define IA(i0) iap[i0]

      ia(I)=ia0(I);
      if( numberOfDimensions==2 )
      {
	for( int i=0; i<numberToCheck; i++ )
	{
	  XA(i,0)=X(IA(i),0);
	  XA(i,1)=X(IA(i),1);
	}
      }
      else
      {
	for( int i=0; i<numberToCheck; i++ )
	{
	  XA(i,0)=X(IA(i),0);
	  XA(i,1)=X(IA(i),1);
	  XA(i,2)=X(IA(i),2);
	}
      }
      
      ra=-1;
      mapping.inverseMap(xa,ra);

      // display(ra,"ra after inversion");
      
      // compress possible points based on unit square coordinates
      const real offset0=.5+gridSpacing(0)*2.01;  // allow points to be within 2 ghost lines
      const real offset1=.5+gridSpacing(1)*2.01;  // allow points to be within 2 ghost lines
      const real offset2=.5+gridSpacing(2)*2.01;  // allow points to be within 2 ghost lines
      j=0;
      int i;
      if( numberOfDimensions==2 )
      {
	for( i=0; i<numberToCheck; i++ )
	{
	  if( fabs(RA(i,0)-.5)<offset0 && fabs(RA(i,1)-.5)<offset1 )
	  {
	    if( i!=j )
	    {
	      IA(j)=IA(i);
	      RA(j,0)=RA(i,0);
	      RA(j,1)=RA(i,1);
	    }
	    j++;
	  }
	}
      }
      else
      {
	for( i=0; i<numberToCheck; i++ )
	{
	  if( fabs(RA(i,0)-.5)<offset0 && fabs(RA(i,1)-.5)<offset1 && fabs(RA(i,2)-.5)<offset2 )
	  {
	    if( i!=j )
	    {
	      IA(j)=IA(i);
	      RA(j,0)=RA(i,0);
	      RA(j,1)=RA(i,1);
	      RA(j,2)=RA(i,2);
	    }
	    j++;
	  }
	}
      }
      
      numberToCheck=j;
      if(numberToCheck==0 )
        continue;
      
      I=numberToCheck;
      
      IntegerArray & ip  = interpolationLocation[grid];
      IntegerArray & ip1 = interpolationLocationPlus[grid];

      ip.redim(I,numberOfDimensions); ip1.redim(I,numberOfDimensions);
      int *ipp = ip.Array_Descriptor.Array_View_Pointer1;
      const int ipDim0=ip.getRawDataSize(0);
#define IP(i0,i1) ipp[i0+ipDim0*(i1)]
      int *ip1p = ip1.Array_Descriptor.Array_View_Pointer1;
      const int ip1Dim0=ip1.getRawDataSize(0);
#define IP1(i0,i1) ip1p[i0+ip1Dim0*(i1)]

      RealArray & dra = interpolationCoordinates[grid];
      dra.redim(I,numberOfDimensions);
      real *drap = dra.Array_Descriptor.Array_View_Pointer1;
      const int draDim0=dra.getRawDataSize(0);
#define DRA(i0,i1) drap[i0+draDim0*(i1)]

      RealArray dr(I,numberOfDimensions);
      real *drp = dr.Array_Descriptor.Array_View_Pointer1;
      const int drDim0=dr.getRawDataSize(0);
#define DR(i0,i1) drp[i0+drDim0*(i1)]

      for( axis=0; axis<numberOfDimensions; axis++ )
      { 
	for( i=0; i<numberToCheck; i++ )
	{
          real rr = RA(i,axis)/GRIDSPACING(axis)+INDEXRANGE(0,axis);
	  IP(i,axis)=rr>=0. ? int(rr+.5) : int(rr-.5); // closest point to r
//  	  IP(i,axis)=rr>=0. ? int(rr) : int(rr-1.); // closest point <= to r

	  IP(i,axis)=min(DIMENSION(End,axis)-1,max(DIMENSION(Start,axis)+1,IP(i,axis)));  // may no longer be < r
	  DR(i,axis)=rr-IP(i,axis)-shift;
	}
      }

      // dra(I,Axes)=min(fabs(dr(I,Axes)),1.);
      dra(I,Axes)=fabs(dr(I,Axes));
      // dra(I,Axes)=dr(I,Axes);

      //...........only use 4 points if dra bigger than epsilon, otherwise just use 2 points (ip1==ip),
      //    this lets us  interpolate near interpolation boundaries
      for( axis=0; axis<numberOfDimensions; axis++ )
      {
	for( i=0; i<numberToCheck; i++ )
	    IP1(i,axis)=IP(i,axis)+( DR(i,axis)>0. ? 1 : -1);

	if( isPeriodic(axis) )    // ........periodic wrap
	  for( i=0; i<numberToCheck; i++ )
	    IP1(i,axis)=MODR(IP1(i,axis),axis);
      }

      // define the valid subset of points for interpolation   *wdh* 021015
      // We allow interpolation from ghost points too; the number of ghost points
      // allowed depends on the discretization width
      const int *gridIndexRangep = &mg.gridIndexRange(0,0);
#define gridIndexRange(side,axis) gridIndexRangep[(side)+2*(axis)]
      int iRangep[6];
#define iRange(side,axis) iRangep[(side)+2*(axis)]
      const int extra=mg.discretizationWidth(0)/2;
      for( axis=0; axis<numberOfDimensions; axis++ )
      {
	iRange(0,axis)=gridIndexRange(0,axis)-extra;
	iRange(1,axis)=gridIndexRange(1,axis)+extra;
      }


      // compress interpolation info into arrays ia,ip,ip1,
      j=0;
      if( numberOfDimensions==2 )
      {
	const int * maskp = mask.Array_Descriptor.Array_View_Pointer1;
	const int maskDim0=mask.getRawDataSize(0);
#define MASK(i0,i1) maskp[i0+maskDim0*(i1)]

	for( i=0; i<numberToCheck; i++ )
	{
	  if( MASK(IP(i,0),IP (i,1))!=0 && MASK(IP1(i,0),IP (i,1))!=0 &&
	      MASK(IP(i,0),IP1(i,1))!=0 && MASK(IP1(i,0),IP1(i,1))!=0 )
	  {
	    if( i!=j )
	    {
	       IP(j,0)= IP(i,0);  IP(j,1)= IP(i,1); 
	      IP1(j,0)=IP1(i,0); IP1(j,1)=IP1(i,1); 
	      DRA(j,0)=DRA(i,0); DRA(j,1)=DRA(i,1); 

	      IA(j)=IA(i);
	    }
//           printf(" ia=%i x=(%8.2e,%8.2e) r=(%8.2e,%8.2e) from grid=%i dra=(%e,%e) ip=(%i,%i) ip1=(%i,%i) \n",
//                  ia(j),x(ia(j),0),x(ia(j),1),ra(i,0),ra(i,1),grid,
//                  dra(j,0),dra(j,1),ip(j,0),ip(j,1),ip1(j,0),ip1(j,1));


	    if(min(IP(j,0),IP1(j,0))>=iRange(0,0) && max(IP(j,0),IP1(j,0))<=iRange(1,0) &&
	       min(IP(j,1),IP1(j,1))>=iRange(0,1) && max(IP(j,1),IP1(j,1))<=iRange(1,1) )
	    {
	      STATUS(IA(j)) = interpolated;
	      j++;
	    }
	    else 
	    {
              if( STATUS(IA(j)) != extrapolated )  // if already extrapolated, use that value.
	      {
		STATUS(IA(j)) = extrapolated;
		j++;
	      }
	    }
	  }
	}
      }
      else // 3D
      {
	const int * maskp = mask.Array_Descriptor.Array_View_Pointer2;
	const int maskDim0=mask.getRawDataSize(0);
	const int maskDim1=mask.getRawDataSize(1);
#undef MASK
#define MASK(i0,i1,i2) maskp[i0+maskDim0*(i1+maskDim1*(i2))]

	for( i=0; i<numberToCheck; i++ )
	{
	  if( MASK(IP(i,0),IP (i,1),IP (i,2))!=0 && MASK(IP1(i,0),IP (i,1),IP (i,2))!=0 &&
	      MASK(IP(i,0),IP1(i,1),IP (i,2))!=0 && MASK(IP1(i,0),IP1(i,1),IP (i,2))!=0 &&
              MASK(IP(i,0),IP (i,1),IP1(i,2))!=0 && MASK(IP1(i,0),IP (i,1),IP1(i,2))!=0 &&
	      MASK(IP(i,0),IP1(i,1),IP1(i,2))!=0 && MASK(IP1(i,0),IP1(i,1),IP1(i,2))!=0 )
	  {
	    if( i!=j )
	    {
	       IP(j,0)= IP(i,0);  IP(j,1)= IP(i,1);  IP(j,2)= IP(i,2);   
	      IP1(j,0)=IP1(i,0); IP1(j,1)=IP1(i,1); IP1(j,2)=IP1(i,2); 
	      DRA(j,0)=DRA(i,0); DRA(j,1)=DRA(i,1); DRA(j,2)=DRA(i,2);  

	      IA(j)=IA(i);
	    }

	    if(min(IP(j,0),IP1(j,0))>=iRange(0,0) && max(IP(j,0),IP1(j,0))<=iRange(1,0) &&
	       min(IP(j,1),IP1(j,1))>=iRange(0,1) && max(IP(j,1),IP1(j,1))<=iRange(1,1) &&
	       min(IP(j,2),IP1(j,2))>=iRange(0,2) && max(IP(j,2),IP1(j,2))<=iRange(1,2) )
	    {
	      STATUS(IA(j)) = interpolated;
	      j++;
	    }
	    else 
	    {
              if( STATUS(IA(j)) != extrapolated )  // if already extrapolated, use that value.
	      {
		STATUS(IA(j)) = extrapolated;
		j++;
	      }
	    }
	  }
	}
	
      }  // end else 3d
      NUMBEROFINTERPOLATIONPOINTS(grid)=j;

    }
  }
  
  int numberInterpolated=0;
  int numberExtrapolated=0;
  int j=0;
  for( int i=0; i<numberOfPointsToInterpolate; i++ )
  {
    if( STATUS(i)==interpolated )
    {
      numberInterpolated++;
    }
    else
    {
      if( STATUS(i)==extrapolated ) 
	numberExtrapolated++;
	
      IA0(j)=i;
      j++;
    }
  }
  int numNotAssigned=numberOfPointsToInterpolate-numberInterpolated-numberExtrapolated;
    
  if( debug )
    printf("buildInterpolationInfo: total interpolated=%i extrapolated=%i not assigned=%i\n",
	   numberInterpolated,numberExtrapolated,numNotAssigned);
  if( numNotAssigned>0 )  
    printf("buildInterpolationInfo: WARNING: %i points not assigned!\n",numNotAssigned);
  
  returnValue=-numNotAssigned;
  return returnValue;
}

#undef DIMENSION
#undef INDEXRANGE
#undef GRIDSPACING
#undef gridIndexRange
#undef iRange
#undef X
#undef XA
#undef STATUS
#undef IA0
#undef IA
#undef RA
#undef IP
#undef IP1
#undef DR
#undef DRA


//\begin{>interpolatePointsInclude.tex}{}
int InterpolatePoints::
interpolatePoints(const realCompositeGridFunction & u,
		  realArray & uInterpolated, 
		  const Range & R0/* =nullRange */,           
		  const Range & R1/* =nullRange */,
		  const Range & R2/* =nullRange */,
		  const Range & R3/* =nullRange */,
		  const Range & R4/* =nullRange */ )
//=======================================================================================================
//  /Description:
//    Given some points in space, determine the values of a grid function u. If interpolation
//    is not possible then extrapolate from the nearest grid point. The extrapolation is zero-order
//    so that the value is just set equal to the value from the boundary.
//  /u (input): interpolate values from this grid function
//  /uInterpolated (output): uInterpolated(0:numberOfPointsToInterpolate-1,R0,R1,R2,R3,R4) : interpolated
//      values
//  /R0,R1,...,R4 (input): interpolate these components of the grid function. R0 is the range of values for
//     the first component of u, R1 the values for the second component, etc. By default all components
//      of u are interpolated.
// ==========================================================================================================
{
  // ***************************************************
  // ****** Interpolate points given the stencil *******
  // ***************************************************
  CompositeGrid & cg = *u.getCompositeGrid();

  const int numberOfDimensions = cg.numberOfDimensions();
  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  // determine component ranges to use:
  Range Ra[5] = {R0,R1,R2,R3,R4};  
  int i;
  for( i=0; i<5; i++ )
  {
    if( Ra[i].length()<=0 ) //     if( Ra[i]==nullRange )
      Ra[i] = Range(u.getComponentBase(i),u.getComponentBound(i));  
    else if( Ra[i].getBase()<u.getComponentBase(i) || Ra[i].getBound()>u.getComponentBound(i) )
    {
      cout << "interpolatePoints:ERROR: the component Range R" << i << " is out of range! \n";
      printf("R%i =(%i,%i) but the dimensions for component %i of u are (%i,%i) \n",i,
	     Ra[i].getBase(),Ra[i].getBound(),i,u.getComponentBase(i),u.getComponentBound(i));
      Overture::abort("error");
    }
    else if( i<3 && (Ra[i].getBase()<uInterpolated.getBase(i+1) || Ra[i].getBound()>uInterpolated.getBound(i+1)) )
    {
      cout << "interpolatePoints:ERROR: the component Range R" << i << " is out of range! \n";
      printf("R%i =(%i,%i) but the dimensions for index %i of uInterpolated are (%i,%i) \n",i,
	     Ra[i].getBase(),Ra[i].getBound(),i+1,uInterpolated.getBase(i+1),uInterpolated.getBound(i+1));
      Overture::abort("error");
    }
  }


  // We must check the highest priority grid first since this was the order the points were generated.
  // A point may be extrpolated on a higher priodirt grid but then interpolated on a lower priority grid.
  int grid;
  for( grid=numberOfComponentGrids-1; grid>=0; grid-- )  // check highest priority grid first
  {
    // interpolate from this grid.

    const int num=numberOfInterpolationPoints(grid);
    if( num>0 )
    {
      // printf("----interpolatePointsNew: interp %i points from grid %i\n",num,grid);

      IntegerArray & ia = indirection[grid];
      IntegerArray & ip = interpolationLocation[grid];
      IntegerArray & ip1= interpolationLocationPlus[grid];
      RealArray & dra = interpolationCoordinates[grid];
      const IntegerArray & gid = cg[grid].gridIndexRange();
      
      // display(ia,"ia");
      // display(ip,"ip");
      // display(dra,"dra");
      

      const realArray & ug = u[grid];


      const int *iap = ia.Array_Descriptor.Array_View_Pointer0;
#define IA(i0) iap[i0]

      const int *ipp = ip.Array_Descriptor.Array_View_Pointer1;
      const int ipDim0=ip.getRawDataSize(0);
#define IP(i0,i1) ipp[i0+ipDim0*(i1)]
      const int *ip1p = ip1.Array_Descriptor.Array_View_Pointer1;
      const int ip1Dim0=ip1.getRawDataSize(0);
#define IP1(i0,i1) ip1p[i0+ip1Dim0*(i1)]

      const real *drap = dra.Array_Descriptor.Array_View_Pointer1;
      const int draDim0=dra.getRawDataSize(0);
#define DRA(i0,i1) drap[i0+draDim0*(i1)]

      real *uInterpolatedp = uInterpolated.Array_Descriptor.Array_View_Pointer1;
      const int uInterpolatedDim0=uInterpolated.getRawDataSize(0);
#define UINTERPOLATED(i0,i1) uInterpolatedp[i0+uInterpolatedDim0*(i1)]


      // ...........Bi-Linear Interpolation:
      if( numberOfDimensions==2 )
      {
	const real *ugp = ug.Array_Descriptor.Array_View_Pointer2;
	const int ugDim0=ug.getRawDataSize(0);
	const int ugDim1=ug.getRawDataSize(1);
#define UG(i0,i1,i2) ugp[i0+ugDim0*(i1+ugDim1*(i2))]

	for( int c0=Ra[0].getBase(); c0<=Ra[0].getBound(); c0++)  // *** add more components ****
	{
	  for( int i=0; i<num; i++ )
	  {
	    UINTERPOLATED(IA(i),c0)= 
	      (1.-DRA(i,1))*(
		(1.-DRA(i,0))*UG(IP (i,0),IP(i,1),c0)
		   +DRA(i,0) *UG(IP1(i,0),IP(i,1),c0))
	      + DRA(i,1) *(
		(1.-DRA(i,0))*UG( IP(i,0),IP1(i,1),c0)
		   +DRA(i,0) *UG(IP1(i,0),IP1(i,1),c0));

//  	    if( c0==3 )
//  	    {
//  	      printf(" grid=%i i=%i ia=%i dra=(%7.1e,%7.1e) ip=(%i,%i) ip1=(%i,%i) gid=%i,%i "
//                       "u=(%9.3e,%9.3e,%9.3e,%9.3e) uI=%9.3e\n",
//  		     grid,i,ia(i),dra(i,0),dra(i,1),
//  		     ip(i,0),ip(i,1),ip1(i,0),ip1(i,1),
//                       gid(0,0),gid(1,0),
//                       UG(IP (i,0),IP(i,1),c0),UG(IP1(i,0),IP(i,1),c0),
//                       UG( IP(i,0),IP1(i,1),c0),UG(IP1(i,0),IP1(i,1),c0),
//                       uInterpolated(ia(i),c0));
//  	    }

	  }
	}
      }
      else // 3D
      {
	const real *ugp = ug.Array_Descriptor.Array_View_Pointer3;
	const int ugDim0=ug.getRawDataSize(0);
	const int ugDim1=ug.getRawDataSize(1);
	const int ugDim2=ug.getRawDataSize(2);
#undef UG
#define UG(i0,i1,i2,i3) ugp[i0+ugDim0*(i1+ugDim1*(i2+ugDim2*(i3)))]

	for( int c0=Ra[0].getBase(); c0<=Ra[0].getBound(); c0++)  // *** add more components ****
	{
	  for( int i=0; i<num; i++ )
	  {
	    UINTERPOLATED(IA(i),c0)= 
              (1.-DRA(i,2))*(
  	      (1.-DRA(i,1))*(
		(1.-DRA(i,0))*UG(IP (i,0),IP(i,1),IP(i,2),c0)
		   +DRA(i,0) *UG(IP1(i,0),IP(i,1),IP(i,2),c0))
	      + DRA(i,1) *(
		(1.-DRA(i,0))*UG( IP(i,0),IP1(i,1),IP(i,2),c0)
		   +DRA(i,0) *UG(IP1(i,0),IP1(i,1),IP(i,2),c0))
                           )
              + DRA(i,2)*(
  	      (1.-DRA(i,1))*(
		(1.-DRA(i,0))*UG(IP (i,0),IP(i,1),IP1(i,2),c0)
		   +DRA(i,0) *UG(IP1(i,0),IP(i,1),IP1(i,2),c0))
	      + DRA(i,1) *(
		(1.-DRA(i,0))*UG( IP(i,0),IP1(i,1),IP1(i,2),c0)
		   +DRA(i,0) *UG(IP1(i,0),IP1(i,1),IP1(i,2),c0))
		         );
	  }
	}
      }
    }
  }
  
  return 0;

}

#undef IP
#undef IP1
#undef UG
#undef IA
#undef UINTERPOLATED

//\begin{>interpolatePointsInclude.tex}{}
int InterpolatePoints::
interpolationCoefficients(const realArray & positionToInterpolate,
			  const CompositeGrid &cg,
			  realArray & uInterpolationCoeff )

{
  //=======================================================================================================
  //  /Description:
  //    Return the coefficients for the interpolation of a grid function u at some points in space. (kkc)
  //    If interpolation
  //    is not possible then extrapolate from the nearest grid point. The extrapolation is zero-order
  //    so that the value is just set equal to the value from the boundary.
  //  /cg (input): interpolate values from this grid 
  //  /uInterpolationCoeff (output): uInterpolationCoeff(0:numberOfPointsToInterpolate-1, 2^numberOfDimensions)
  //      interpolation coefficients
  //\end{interpolatePointsInclude.tex}  
  // ==========================================================================================================


  // 030228 kkc, most of this code is from interpolatePoints, the only difference is that the coefficients are stored
  //             instead of the interpolated value

  // note the extra index (the second in the list) to u, this entry holds the index to the coefficent for a 
  //      node in the parametric element.  they are ordered by index in parameter space. Hence,
  //      in a 3D grid, the coefficient at (ir1,ir2,ir3) = ir1 + ( 2*(ir2 + 2*ir3)), so (0,0,1) would be at
  //      index 4.

  //CompositeGrid & cg = *u.getCompositeGrid();

  const int numberOfDimensions = cg.numberOfDimensions();
  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  // determine component ranges to use:
//   Range Ra[5] = {R0,R1,R2,R3,R4};  
//   int i;
//   for( i=0; i<5; i++ )
//   {
//     if( Ra[i].length()<=0 ) //     if( Ra[i]==nullRange )
//       Ra[i] = Range(u.getComponentBase(i),u.getComponentBound(i));  
//     else if( Ra[i].getBase()<u.getComponentBase(i) || Ra[i].getBound()>u.getComponentBound(i) )
//     {
//       cout << "interpolationCoeffictions:ERROR: the component Range R" << i << " is out of range! \n";
//       printf("R%i =(%i,%i) but the dimensions for component %i of u are (%i,%i) \n",i,
// 	     Ra[i].getBase(),Ra[i].getBound(),i,u.getComponentBase(i),u.getComponentBound(i));
//       Overture::abort("error");
//     }
//     else if( i<3 && (Ra[i].getBase()<uInterpolated.getBase(i+1) || Ra[i].getBound()>uInterpolated.getBound(i+1)) )
//     {
//       cout << "interpolationCoefficients:ERROR: the component Range R" << i << " is out of range! \n";
//       printf("R%i =(%i,%i) but the dimensions for index %i of uInterpolated are (%i,%i) \n",i,
// 	     Ra[i].getBase(),Ra[i].getBound(),i+1,uInterpolated.getBase(i+1),uInterpolated.getBound(i+1));
//       Overture::abort("error");
//     }
//   }


  // We must check the highest priority grid first since this was the order the points were generated.
  // A point may be extrpolated on a higher priodirt grid but then interpolated on a lower priority grid.
  int grid;
  for( grid=numberOfComponentGrids-1; grid>=0; grid-- )  // check highest priority grid first
  {
    // interpolate from this grid.

    const int num=numberOfInterpolationPoints(grid);
    if( num>0 )
    {
      // printf("----interpolatePointsNew: interp %i points from grid %i\n",num,grid);

      IntegerArray & ia = indirection[grid];
      IntegerArray & ip = interpolationLocation[grid];
      IntegerArray & ip1= interpolationLocationPlus[grid];
      RealArray & dra = interpolationCoordinates[grid];
      const IntegerArray & gid = cg[grid].gridIndexRange();
      
      // display(ia,"ia");
      // display(ip,"ip");
      // display(dra,"dra");
      

      //      const realArray & ug = u[grid];


      const int *iap = ia.Array_Descriptor.Array_View_Pointer0;
#define IA(i0) iap[i0]

      const int *ipp = ip.Array_Descriptor.Array_View_Pointer1;
      const int ipDim0=ip.getRawDataSize(0);
#define IP(i0,i1) ipp[i0+ipDim0*(i1)]
      const int *ip1p = ip1.Array_Descriptor.Array_View_Pointer1;
      const int ip1Dim0=ip1.getRawDataSize(0);
#define IP1(i0,i1) ip1p[i0+ip1Dim0*(i1)]

      const real *drap = dra.Array_Descriptor.Array_View_Pointer1;
      const int draDim0=dra.getRawDataSize(0);
#define DRA(i0,i1) drap[i0+draDim0*(i1)]

      real *uInterpolatedp = uInterpolationCoeff.Array_Descriptor.Array_View_Pointer1;
      const int uInterpolatedDim0=uInterpolationCoeff.getRawDataSize(0);
      // kkc#define UINTERPOLATED(i0,i1) uInterpolatedp[i0+uInterpolatedDim0*(i1)]
      int nD = numberOfDimensions==2 ? 4 : 8;
#define UINTERPOLATED(i0,ix,iy,iz) uInterpolatedp[i0+ uInterpolatedDim0*( (ix)+2*((iy) +2*(iz)))]


      // ...........Bi-Linear Interpolation:
      if( numberOfDimensions==2 )
      {
// 	const real *ugp = ug.Array_Descriptor.Array_View_Pointer2;
// 	const int ugDim0=ug.getRawDataSize(0);
// 	const int ugDim1=ug.getRawDataSize(1);
#define UG(i0,i1,i2) ugp[i0+ugDim0*(i1+ugDim1*(i2))]

	//	for( int c0=Ra[0].getBase(); c0<=Ra[0].getBound(); c0++)  // *** add more components ****
	//	{
	int c0=0;
	  for( int i=0; i<num; i++ )
	  {
// 	    UINTERPOLATED(IA(i),c0)= 
// 	      (1.-DRA(i,1))*(
// 		(1.-DRA(i,0))*UG(IP (i,0),IP(i,1),c0)
// 		   +DRA(i,0) *UG(IP1(i,0),IP(i,1),c0))
// 	      + DRA(i,1) *(
// 		(1.-DRA(i,0))*UG( IP(i,0),IP1(i,1),c0)
// 		   +DRA(i,0) *UG(IP1(i,0),IP1(i,1),c0));

	    UINTERPOLATED(IA(i),0,0,0)= 
	      (1.-DRA(i,1))*(1.-DRA(i,0));

	    UINTERPOLATED(IA(i),0,1,0)= 
	      DRA(i,1) *(1.-DRA(i,0));

	    UINTERPOLATED(IA(i),1,0,0)= 
	      (1.-DRA(i,1))*DRA(i,0);

	    UINTERPOLATED(IA(i),1,1,0)= 
	      DRA(i,1) * DRA(i,0);

//  	    if( c0==3 )
//  	    {
//  	      printf(" grid=%i i=%i ia=%i dra=(%7.1e,%7.1e) ip=(%i,%i) ip1=(%i,%i) gid=%i,%i "
//                       "u=(%9.3e,%9.3e,%9.3e,%9.3e) uI=%9.3e\n",
//  		     grid,i,ia(i),dra(i,0),dra(i,1),
//  		     ip(i,0),ip(i,1),ip1(i,0),ip1(i,1),
//                       gid(0,0),gid(1,0),
//                       UG(IP (i,0),IP(i,1),c0),UG(IP1(i,0),IP(i,1),c0),
//                       UG( IP(i,0),IP1(i,1),c0),UG(IP1(i,0),IP1(i,1),c0),
//                       uInterpolated(ia(i),c0));
//  	    }

	  }
	  //	}
      }
      else // 3D
      {
// 	const real *ugp = ug.Array_Descriptor.Array_View_Pointer3;
// 	const int ugDim0=ug.getRawDataSize(0);
// 	const int ugDim1=ug.getRawDataSize(1);
// 	const int ugDim2=ug.getRawDataSize(2);
#undef UG
#define UG(i0,i1,i2,i3) ugp[i0+ugDim0*(i1+ugDim1*(i2+ugDim2*(i3)))]

	//	for( int c0=Ra[0].getBase(); c0<=Ra[0].getBound(); c0++)  // *** add more components ****
	//	{
	int c0=0;
	  for( int i=0; i<num; i++ )
	  {
// 	    UINTERPOLATED(IA(i),c0)= 
//               (1.-DRA(i,2))*(
//   	      (1.-DRA(i,1))*(
// 		(1.-DRA(i,0))*UG(IP (i,0),IP(i,1),IP(i,2),c0)
// 		   +DRA(i,0) *UG(IP1(i,0),IP(i,1),IP(i,2),c0))
// 	      + DRA(i,1) *(
// 		(1.-DRA(i,0))*UG( IP(i,0),IP1(i,1),IP(i,2),c0)
// 		   +DRA(i,0) *UG(IP1(i,0),IP1(i,1),IP(i,2),c0))
//                            )
//               + DRA(i,2)*(
//   	      (1.-DRA(i,1))*(
// 		(1.-DRA(i,0))*UG(IP (i,0),IP(i,1),IP1(i,2),c0)
// 		   +DRA(i,0) *UG(IP1(i,0),IP(i,1),IP1(i,2),c0))
// 	      + DRA(i,1) *(
// 		(1.-DRA(i,0))*UG( IP(i,0),IP1(i,1),IP1(i,2),c0)
// 		   +DRA(i,0) *UG(IP1(i,0),IP1(i,1),IP1(i,2),c0))
// 		         );

	    UINTERPOLATED(IA(i),0,0,0)= (1.-DRA(i,2))*(1.-DRA(i,1))*(1.-DRA(i,0));
	    UINTERPOLATED(IA(i),0,1,0)= (1.-DRA(i,2))*DRA(i,1) *(1.-DRA(i,0));
	    UINTERPOLATED(IA(i),1,0,0)= (1.-DRA(i,2))*(1.-DRA(i,1))*DRA(i,0);
	    UINTERPOLATED(IA(i),1,1,0)= (1.-DRA(i,2))*DRA(i,1)*DRA(i,0);

	    UINTERPOLATED(IA(i),0,0,1)= DRA(i,2)*(1.-DRA(i,1))*(1.-DRA(i,0));
	    UINTERPOLATED(IA(i),0,1,1)= DRA(i,2)*DRA(i,1) *(1.-DRA(i,0));
	    UINTERPOLATED(IA(i),1,0,1)= DRA(i,2)*(1.-DRA(i,1))*DRA(i,0);
	    UINTERPOLATED(IA(i),1,1,1)= DRA(i,2)*DRA(i,1)*DRA(i,0);

	  }
	  //}
      }
    }
  }
  
  return 0;

}

#undef IP
#undef IP1
#undef UG
#undef IA
#undef UINTERPOLATED


//\begin{>interpolatePointsInclude.tex}{}
int InterpolatePoints::
interpolatePoints(const realArray & positionToInterpolate,
		     const realCompositeGridFunction & u,
		     realArray & uInterpolated, 
		     const Range & R0/* =nullRange */,           
		     const Range & R1/* =nullRange */,
		     const Range & R2/* =nullRange */,
		     const Range & R3/* =nullRange */,
		     const Range & R4/* =nullRange */ )
//=======================================================================================================
//  /Description:
//    Given some points in space, determine the values of a grid function u. If interpolation
//    is not possible then extrapolate from the nearest grid point. The extrapolation is zero-order
//    so that the value is just set equal to the value from the boundary.
//  /positionToInterpolate (input):
//     positionToInterpolate(0:numberOfPointsToInterpolate-1,0:numberOfDimensions-1) : (x,y[,z]) positions
//          to interpolate. The first dimension of this array determines how many points to interpolate.
//  /u (input): interpolate values from this grid function
//  /uInterpolated (output): uInterpolated(0:numberOfPointsToInterpolate-1,R0,R1,R2,R3,R4) : interpolated
//      values
//  /R0,R1,...,R4 (input): interpolate these components of the grid function. R0 is the range of values for
//     the first component of u, R1 the values for the second component, etc. By default all components
//      of u are interpolated.
//  /indexGuess (input/ouput): indexGuess(0:numberOfPointsToInterpolate-1,0:numberOfDimensions-1) : 
//    (i1,i2[,i3]) values for initial 
//        guess for searches. Not required by default.
//  /interpoleeGrid(.) (input/output): interpoleeGrid(0:numberOfPointsToInterpolate-1) : try
//        this grid first. Not required by default. 
//  /wasInterpolated(.) (output) : If provided as an argument, on output wasInterpolated(i)=TRUE if the point
//     was successfully interpolated, or wasInterpolated(i)=FALSE if the point was extrapolated.
//  /Errors:  This routine in principle should always be able to interpolate or extrapolate.
//  /Return Values:
//    \begin{itemize}
//      \item 0 = success
//      \item 1 = error, unable to interpolate (this should never happen)
//      \item -N = could not interpolate N points, but could extrapolate -- extrapolation was performed
//         from the nearest grid point.
//    \end{itemize}
//  /Author: WDH
//\end{interpolatePointsInclude.tex}  
// =======================================================================================================
{

  int returnValue=0;
  CompositeGrid & cg = *u.getCompositeGrid();
  
  returnValue=buildInterpolationInfo(positionToInterpolate,cg );

  interpolatePoints(u,uInterpolated,R0,R1,R2,R3,R4);
  
  return returnValue;
  
}

#undef NRM
#undef MODR

//\begin{>interpolateAllPointsInclude.tex}{}
int InterpolatePoints::
interpolateAllPoints(const realCompositeGridFunction & uFrom,
		     realCompositeGridFunction & uTo )
//==============================================================================
//
// /Description:
//     Interpolate all values on one CompositeGridFunction, {\ff uTo},  
//   from the values of another CompositeGridFunction,
//   {\ff uFrom}. Values on {\ff uTo} are extrapolated if they lie outside the region covered by {\ff uFrom}.
//   This routine calls the {\ff interpolatePoints} function.
// /uFrom (input):
//      Use these values to interpolate from.
// /uTo (output):
//      Fill in all values on this grid (including ghost-points).
// /Errors:  This routine in principle should always be able to interpolate or extrapolate all
//   values.
// /Return Values:
//     \begin{itemize}
//       \item 0 = success
//       \item 1 = error, unable to interpolate 
//       \item -N = could not interpolate N points, but could extrapolate -- extrapolation was performed
//          from the nearest grid point.
//     \end{itemize}
//
// /Author: WDH
//
//\end{interpolateAllPointsInclude.tex}
//==============================================================================
{
  CompositeGrid & cgTo= (CompositeGrid&) *uTo.gridCollection;
  int numberOfExtrapolatedPoints=0;
  for( int grid=0; grid<cgTo.numberOfComponentGrids(); grid++)
  {
    numberOfExtrapolatedPoints+=interpolateAllPoints(uFrom,uTo[grid]);
  }
  return numberOfExtrapolatedPoints;

/* ---
  
  Range C0 = Range(uTo.getComponentBase(0),uTo.getComponentBound(0));


  // Index I1,I2,I3;
  int numberOfExtrapolatedPoints=0;
  int grid;
  for( grid=0; grid<cgTo.numberOfComponentGrids(); grid++)
  {
    // make a list of points to interpolate. No need to interpolate points with mask==0
    const intArray & mask = cgTo[grid].mask();
    const realArray & center = cgTo[grid].center();
    
    intArray ia;
    ia = (mask!=0).indexMap();
    
    if( ia.getLength(0)>0 )
    {
      Range I=ia.getLength(0);
      realArray x(I,cgTo.numberOfDimensions()), uInterpolated(I,C0);
    
      const int i3=center.getBase(2);
      if( cgTo.numberOfDimensions()==2 )
      {
	for( int axis=0; axis<2; axis++ )
	  x(I,axis)=center(ia(I,0),ia(I,1),i3,axis);
      }
      else
      {
	for( int axis=0; axis<3; axis++ )
	  x(I,axis)=center(ia(I,0),ia(I,1),ia(I,2),axis);
      }
    
      int num=interpolatePoints(x,uFrom,uInterpolated);

      realArray & u = uTo[grid];
      if( cgTo.numberOfDimensions()==2 )
      {
	for( int c0=C0.getBase(); c0<=C0.getBound(); c0++ )  // *** could avoid this copy if right shape
	  u(ia(I,0),ia(I,1),i3,c0)=uInterpolated(I,c0);
      }
      else
      {
	for( int c0=C0.getBase(); c0<=C0.getBound(); c0++ )  // *** could avoid this copy if right shape
	  u(ia(I,0),ia(I,1),ia(I,2),c0)=uInterpolated(I,c0);
      }

      // printf("interpolatePoints: number of extrapolated points on grid %i = %i\n",grid,num);
      numberOfExtrapolatedPoints-=num;

    }
    
  }
  return numberOfExtrapolatedPoints;
#else
  cout << "interpolateAllPoints:Error: not implemented for P++ yet \n";
  Overture::abort("error");
  return 0;
#endif
---- */

}



//\begin{>>interpolateAllPointsInclude.tex}{}
int InterpolatePoints::
interpolateAllPoints(const realCompositeGridFunction & uFrom,
                     realMappedGridFunction & uTo )
//==============================================================================
//
// /Description:
//     Interpolate all values on a realMappedGridFunction, {\ff uTo},  
//   from the values of another CompositeGridFunction,
//   {\ff uFrom}. Values on {\ff uTo} are extrapolated if they lie outside the region covered by {\ff uFrom}.
//   This routine calls the {\ff interpolatePoints} function.
// /uFrom (input):
//      Use these values to interpolate from.
// /uTo (output):
//      Fill in all values on this grid (including ghost-points).
// /Errors:  This routine in principle should always be able to interpolate or extrapolate all
//   values.
// /Return Values:
//     \begin{itemize}
//       \item 0 = success
//       \item 1 = error, unable to interpolate 
//       \item -N = could not interpolate N points, but could extrapolate -- extrapolation was performed
//          from the nearest grid point.
//     \end{itemize}
//
// /Author: WDH
//
//\end{interpolateAllPointsInclude.tex}
//==============================================================================
{
  int numberOfExtrapolatedPoints=0;
  
  MappedGrid & mg= *uTo.getMappedGrid();
  
  // make a list of points to interpolate. No need to interpolate points with mask==0
  const intArray & mask = mg.mask();
    
  intArray ia;
  ia = (mask!=0).indexMap();

  const int *iap = ia.Array_Descriptor.Array_View_Pointer1;
  const int iaDim0=ia.getRawDataSize(0);
#define IA(i0,i1) iap[i0+iaDim0*(i1)]
    
  const int numToInterpolate=ia.getLength(0);
  if( numToInterpolate>0 )
  {
    Range I=ia.getLength(0);
    Range C0 = Range(uTo.getComponentBase(0),uTo.getComponentBound(0));
    realArray x(I,mg.numberOfDimensions()), uInterpolated(I,C0);
    
    real *xp = x.Array_Descriptor.Array_View_Pointer1;
    const int xDim0=x.getRawDataSize(0);
#define X(i0,i1) xp[i0+xDim0*(i1)]
    real *uInterpolatedp = uInterpolated.Array_Descriptor.Array_View_Pointer1;
    const int uInterpolatedDim0=uInterpolated.getRawDataSize(0);
#define UINTERPOLATED(i0,i1) uInterpolatedp[i0+uInterpolatedDim0*(i1)]


    if( mg.isRectangular() )
    {
	real dx[3],xab[2][3];
	mg.getRectangularGridParameters( dx, xab );

	const int i0a=mg.gridIndexRange(0,0);
	const int i1a=mg.gridIndexRange(0,1);
	const int i2a=mg.gridIndexRange(0,2);

        const real xa=xab[0][0], dx0=dx[0];
        const real ya=xab[0][1], dy0=dx[1];
        const real za=xab[0][2], dz0=dx[2];
	
#define COORD0(i0,i1,i2) (xa+dx0*(i0-i0a))
#define COORD1(i0,i1,i2) (ya+dy0*(i1-i1a))
#define COORD2(i0,i1,i2) (za+dz0*(i2-i2a))

      if( mg.numberOfDimensions()==2 )
      {
        for( int i=0; i<numToInterpolate; i++ )
	{
	  X(i,0)=COORD0(IA(i,0),IA(i,1),0);
	  X(i,1)=COORD1(IA(i,0),IA(i,1),0);
	}
      }
      else
      {
        for( int i=0; i<numToInterpolate; i++ )
	{
	  X(i,0)=COORD0(IA(i,0),IA(i,1),IA(i,2));
	  X(i,1)=COORD1(IA(i,0),IA(i,1),IA(i,2));
	  X(i,2)=COORD2(IA(i,0),IA(i,1),IA(i,2));
	}
      }
    }
    else
    {
      mg.update(MappedGrid::THEcenter);
      
      const realArray & center = mg.center();
      const real *centerp = center.Array_Descriptor.Array_View_Pointer3;
      const int centerDim0=center.getRawDataSize(0);
      const int centerDim1=center.getRawDataSize(1);
      const int centerDim2=center.getRawDataSize(2);
#define CENTER(i0,i1,i2,i3) centerp[i0+centerDim0*(i1+centerDim1*(i2+centerDim2*(i3)))]
      const int i3=center.getBase(2);
      if( mg.numberOfDimensions()==2 )
      {
        for( int i=0; i<numToInterpolate; i++ )
	{
	  X(i,0)=CENTER(IA(i,0),IA(i,1),i3,0);
	  X(i,1)=CENTER(IA(i,0),IA(i,1),i3,1);
	}
      }
      else
      {
        for( int i=0; i<numToInterpolate; i++ )
	{
	  X(i,0)=CENTER(IA(i,0),IA(i,1),IA(i,2),0);
	  X(i,1)=CENTER(IA(i,0),IA(i,1),IA(i,2),1);
	  X(i,2)=CENTER(IA(i,0),IA(i,1),IA(i,2),2);
	}
      }
    }
    
    
    int num=interpolatePoints(x,uFrom,uInterpolated);

    realArray & u = uTo;
    real *up = u.Array_Descriptor.Array_View_Pointer3;
    const int uDim0=u.getRawDataSize(0);
    const int uDim1=u.getRawDataSize(1);
    const int uDim2=u.getRawDataSize(2);
#undef U
#define U(i0,i1,i2,i3) up[i0+uDim0*(i1+uDim1*(i2+uDim2*(i3)))]

    if( mg.numberOfDimensions()==2 )
    {
      const int i3=u.getBase(2);
      for( int c0=C0.getBase(); c0<=C0.getBound(); c0++ )  // *** could avoid this copy if right shape
	for( int i=0; i<numToInterpolate; i++ )
	  U(IA(i,0),IA(i,1),i3,c0)=UINTERPOLATED(i,c0);
    }
    else
    {
      for( int c0=C0.getBase(); c0<=C0.getBound(); c0++ )  // *** could avoid this copy if right shape
	for( int i=0; i<numToInterpolate; i++ )
	  U(IA(i,0),IA(i,1),IA(i,2),c0)=UINTERPOLATED(i,c0);
    }

    // printf("interpolatePoints: number of extrapolated points on grid %i = %i\n",grid,num);
    numberOfExtrapolatedPoints-=num;

  }
    
  return numberOfExtrapolatedPoints;

}

#undef CENTER
#undef IA
#undef X
#undef UINTERPOLATED
#undef COORD0
#undef COORD1
#undef COORD2
#undef U
