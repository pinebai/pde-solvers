#include "Overture.h"
#include <float.h>

#include "FluxBCBase.h"

#include "CompositeGridOperators.h"
#include "interpolateFluxBoundary.h"

#define FLUXBC_ForBoundary(side,axis) \
    for( axis=0; axis<mg.numberOfDimensions(); axis++ ) \
      for( side=0; side<=1; side++ )   


FluxBCBase::FluxBCBase( ) 
{ 
  //interpDebug=7;  // 1+2+4+8;
  interpDebug=0;  // 1+2+4+8;
  debug      =0;

  pCg = NULL; pInterp=NULL;     
  idFluxBoundary = 2;
  fluxBCType = interpolateFluxBC;
  fluxCoefficient = 1.;

  doNotUseTwilightZoneFlow(); //by default

}

FluxBCBase::
~FluxBCBase( ) 
{
} 

void FluxBCBase::
updateToMatchGrid( CompositeGrid &cg, int bcID /* =0*/ )    //TODO: precompute the interp. stencils HERE!!
{ 
  this->setCompositeGrid( cg ); 
}

void FluxBCBase::
setCompositeGrid( CompositeGrid &cg ) 
{ 
  pCg = &cg; 
  if( pCg == NULL ) printf("ERROR FluxBC::setCompositeGrid has a problem -- pCg == NULL\n");
  //interpolatePoints ....
};

void FluxBCBase:: 
setInterpolant( Interpolant &interpolant )
{ 
  pInterp = &interpolant; 
}

void FluxBCBase:: 
setOperators( CompositeGridOperators &operators )
{ 
  pOp = & operators; 
}

FluxBCBase::FluxBCType FluxBCBase::
setFluxBCType( std::string &name ) 
{
  fluxBCType =  this->convertToFluxBCType( name );
  return( this->getFluxBCType() );
}
  
FluxBCBase::FluxBCType  FluxBCBase::
convertToFluxBCType( std::string &name ) 
{
  FluxBCType ft;
  
  if ( name == "noFlux" ) {
    ft = FluxBCBase::noFluxBC;
  }
  else if ( name == "interpolateFlux" ) {
    ft = FluxBCBase::interpolateFluxBC;
  }
  else {
    std::cerr << "Unknown bcType type = ["
	      << name << "]" << std::endl;
    ft = FluxBCBase::noFluxBC;
  }
  return ( ft );
}

FluxBCBase::FluxBCType  FluxBCBase::
getFluxBCType() 
{ 
  return fluxBCType; 
};

void  FluxBCBase::
setFluxBoundaryID( int id )  //..then call 'update to match grid'
{ 
  idFluxBoundary = id; 
};

int   FluxBCBase::getFluxBoundaryID() 
{ 
  return idFluxBoundary; 
};

void    FluxBCBase::setFluxCoefficient( double f ) 
{ 
  if(debug&4) printf("--flux coefficient = %g\n", f );
  fluxCoefficient = f; 
};

double  FluxBCBase::getFluxCoefficient( ) 
{ 
  return fluxCoefficient; 
};

void     FluxBCBase::setDebug( int debug_ ) 
{ 
  debug= debug_; 
};

int      FluxBCBase::
getDebug() 
{ 
  return debug; 
};

double  FluxBCBase::
getInterpSetupTimer() 
{ 
  return timerInterpSetupCode; 
};

double  FluxBCBase::
getInterpTimer() 
{ 
  return timerInterpCode; 
};

void  FluxBCBase::
useTwilightZoneFlow( OGFunction &tzSolution)
{
  pExactSolution = &tzSolution; //TZ always assumed a continuous fcn
  useTZ=true;
}

void  FluxBCBase::
doNotUseTwilightZoneFlow()
{
  pExactSolution = NULL;
  useTZ=false;
}

bool  FluxBCBase::
isUsingTwilightZoneFlow()
{
  return useTZ;
}

void  FluxBCBase::
applyBoundaryCondition( realCompositeGridFunction &qfrom,
			realCompositeGridFunction &qto,
			int ic, 
			double tcomp /*=0.*/ ) 
{
  assert( pCg     != NULL ); CompositeGrid &cg          = *pCg;
  assert( pInterp != NULL ); Interpolant   &interp      = *pInterp;
  assert( pOp     != NULL ); CompositeGridOperators &op = *pOp;

  //printf("FluxBCBase::applyBoundaryCondition called...\n");

  timerInterpCode = 0.; // read the timer after this call with getInterpTimer()
  for ( int ig=0; ig< cg.numberOfComponentGrids(); ++ig ) {
    MappedGrid & mg = cg[ig];
    MappedGridOperators      &opmg = op[ig];
    realMappedGridFunction &qto_mg = qto[ig];
    
    Index Ib1,Ib2,Ib3;
    Index Ig1,Ig2,Ig3;
    
    //..flux/jump bc
    int axis, side;
    for( axis=0; axis<mg.numberOfDimensions(); axis++ ) {
      for( side=0; side<=1; side++ ) {
	
	if( mg.boundaryCondition()(side,axis) == idFluxBoundary  ) {
	  //printf("FluxBCBase: grid, side, axis(%i,%i,%i) is a flux boundary (bc=%i)\n",
	  //        ig,side,axis,idFluxBoundary);
	  getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
	  getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
	  const realArray &zz=mg.vertex();
	  int nAxes        =cg.numberOfDimensions();
	  int zaxis        =axis3;
	  if (nAxes==2) zaxis=axis2; // to make sure we don't seg.fault
	  const realArray &xb  = mg.vertex()(Ib1,Ib2,Ib3, axis1);
	  const realArray &yb  = mg.vertex()(Ib1,Ib2,Ib3, axis2);
	  const realArray &zb  = mg.vertex()(Ib1,Ib2,Ib3, zaxis);
	  
	  int axisLength[3]={xb.getLength(axis1),xb.getLength(axis2),xb.getLength(axis3)};
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
	  
	  realArray       &qf     = qto[ig];
	  realArray jump(nInterpPoints);

	  double timerTemp=getCPU();
	  interpolatePoints(ig,xyInt, qfrom, jump);//.. interpolate values from qfrom... TIMEHOG!!
	  timerInterpCode+= getCPU()-timerTemp;
	  jump.reshape(Ib1,Ib2,Ib3);
	  if(debug &4)   jump.display("value"); 
	  const double fluxCoeff=getFluxCoefficient();
	  jump = getFluxCoefficient()*( jump - qf(Ib1,Ib2,Ib3));
	  if(debug &4)    jump.display("jump");
	  
	  //.. ibc chosen to set a specified boundary (side,axis)
	  // ... take jumps from unext --> impose as flux in u
	  const int ibc=BCTypes::boundary1+side+2*axis;
	  //..Add Twilight zone correction if used
	  if( isUsingTwilightZoneFlow()) {
	    //.. add stuff to jump, on grid=ig, side, axis
	    OGFunction & exact = *pExactSolution; // make a reference for readibility
	    realArray & normal = mg.vertexBoundaryNormal(side,axis);    
	    jump +=  normal(Ib1,Ib2,Ib3,axis1)*exact.x(mg, Ib1, Ib2, Ib3, 0, tcomp)
 	            + normal(Ib1,Ib2,Ib3,axis2)*exact.y(mg, Ib1, Ib2, Ib3, 0, tcomp);
	  }
	  opmg.applyBoundaryCondition( qto_mg, Range(ic,ic), BCTypes::neumann, ibc, jump ); 
	} // end if ibc==fluxBoundaryID
      } 
    }
  }
}


//
// --- FluxBCBase utilities
// 

#ifdef GETLENGTH
#define GET_LENGTH dimension
#else
#define GET_LENGTH getLength
#endif

// The macro MODR shifts a point back into the main periodic region
#define NRM(axis,grid)  ( cg[grid].indexRange()(End,axis)-cg[grid].indexRange()(Start,axis)+1 )
#define MODR(i,axis,grid)  ( \
  ( (i-cg[grid].indexRange()(Start,axis)+NRM(axis,grid)) % NRM(axis,grid)) \
      +cg[grid].indexRange()(Start,axis) \
                           )

//\begin{>interpolatePointsInclude.tex}{}
int FluxBCBase::
interpolatePoints(int iThisGrid,
		    const realArray & positionToInterpolate,
		    const realCompositeGridFunction & u,
		    realArray & uInterpolated, 
		    intArray & indexGuess=Overture::nullIntegerDistributedArray(),
		    intArray & interpoleeGrid=Overture::nullIntegerDistributedArray(),
		    intArray & wasInterpolated=Overture::nullIntegerDistributedArray())
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

  const Range & R0=nullRange;
  const Range & R1=nullRange;
  const Range & R2=nullRange;
  const Range & R3=nullRange;
  const Range & R4=nullRange;

  CompositeGrid & cg = (CompositeGrid&) *u.gridCollection;
  cg.update(MappedGrid::THEcenter);


  int grid, axis;
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    if( cg[grid].mapping().mapPointer==NULL )
    {
      cout << "interpolatePoints:ERROR: grids must have mappings! \n";
      exit(1);
    }
    if( !cg[grid].isAllVertexCentered() && !cg[grid].isAllCellCentered() )
    {
      cout << "interpolatePoints:ERROR: grids must be either vertex or cell centered, no mongrels! \n";
      return 1;
    }      
  }
  // determine component ranges to use:
  Range R[5] = {R0,R1,R2,R3,R4};  
  int i;
  for( i=0; i<5; i++ )
  {
    if( R[i].length()<=0 ) //     if( R[i]==nullRange )
      R[i] = Range(u.getComponentBase(i),u.getComponentBound(i));  
    else if( R[i].getBase()<u.getComponentBase(i) || R[i].getBound()>u.getComponentBound(i) )
    {
      cout << "interpolatePoints:ERROR: the component Range R" << i << " is out of range! \n";
      printf("R%i =(%i,%i) but the dimensions for component %i of u are (%i,%i) \n",i,
	     R[i].getBase(),R[i].getBound(),i,u.getComponentBase(i),u.getComponentBound(i));
      Overture::abort("error");
    }
    else if( i<3 && (R[i].getBase()<uInterpolated.getBase(i+1) || R[i].getBound()>uInterpolated.getBound(i+1)) )
    {
      cout << "interpolatePoints:ERROR: the component Range R" << i << " is out of range! \n";
      printf("R%i =(%i,%i) but the dimensions for index %i of uInterpolated are (%i,%i) \n",i,
	     R[i].getBase(),R[i].getBound(),i+1,uInterpolated.getBase(i+1),uInterpolated.getBound(i+1));
      Overture::abort("error");
    }
  }

  int numberOfPointsToInterpolate=positionToInterpolate.GET_LENGTH(0);

  if( uInterpolated.GET_LENGTH(0) < numberOfPointsToInterpolate )
  {
    cout << "interpolatePoints::ERROR: there is not enough space to hold all the interpolated values\n";
    cout << "numberOfPointsToInterpolate=positionToInterpolate.getLength(0)=" << numberOfPointsToInterpolate
         << ", but uInterpolated.getLength(0) = " << uInterpolated.GET_LENGTH(0) << endl;
    return 1;
  }

  const real epsi=1.e-3;
  int extrap,pointWasExtrapolated;
  int returnValue=0;  // 0=ok, >0 error, <0 some points extrapolated
  char buff[100];


  intArray ip(3),ip1(3);
  realArray r(1,3),x(1,3),dr(3),dra(3); x=0.;
  Range Axes(0,cg.numberOfDimensions()-1);
  grid=0;  

  for( int ipt=0; ipt<numberOfPointsToInterpolate; ipt++ )
  {
    if( interpoleeGrid.GET_LENGTH(0) >= numberOfPointsToInterpolate )
      grid=min(cg.numberOfComponentGrids()-1,max(0,interpoleeGrid(ipt)));  // here is the first grid we check
    else
      grid=min(cg.numberOfComponentGrids()-1,max(0,grid));


    x(0,Axes)=positionToInterpolate(ipt,Axes);

    if( interpDebug & 8 ) {
      printf("******Inverting point ipt=%i, x=(%f,%f,%f), this grid=%d\n",ipt,x(0,0),x(0,1),x(0,2), iThisGrid);

    }

    real minimumDistance=REAL_MAX;
    real distance;

    // *****************************************************************************
    // *** Loop through the grids until we find a point we can interpolate from ***
    // *****************************************************************************
    pointWasExtrapolated=-1;  // do we have to extrapolate to get a value? 0=ok, 1=extrap, -1=error
    for( int gridn=0; gridn<cg.numberOfComponentGrids(); gridn++ )
    {
      extrap=FALSE;
      if( gridn>0 ) 
        grid = (grid+1) % cg.numberOfComponentGrids();  // here is the next grid to try;
      
      if( cg[grid].getGridType()==GenericGrid::unstructuredGrid )
      {
	continue;
      }

      if (grid == iThisGrid ) continue; // don't interpolate from same grid... this won't work

      
      const IntegerArray & indexRange = cg[grid].indexRange();
      if( cg.numberOfDimensions()==2 )
      {
        ip(2)=ip1(2)=cg[grid].indexRange()(Start,axis3);
      }
      cg[grid].mapping().inverseMap(x,r);   // Invert the mapping
      if( interpDebug & 4 )
      {
        r.display(sPrintF(buff,"Here is the inverse r on grid =%i ",grid));
        realArray z(1,3);
	cg[grid].mapping().map(r,z);
	z.display("Here is z from map(r,z)");
      }
      if( r(0,0)==ApproximateGlobalInverse::bogus )
      {
	if( cg[grid].mapping().mapPointer->approximateGlobalInverse )
	{
	  cg[grid].mapping().mapPointer->approximateGlobalInverse->findNearestGridPoint(0,0,x,r);
          if( interpDebug & 4 )
            r.display("Nearest grid point computed:");
	}
      }
      real shift = (bool)cg[grid].isAllVertexCentered() ? 0. : .5; // shift position for cell centered grids
      for( axis=0; axis<cg.numberOfDimensions(); axis++ )
      { // ip = closest point on the grid
	ip(axis)=int(r(0,axis)/cg[grid].gridSpacing()(axis)+cg[grid].indexRange()(0,axis)); 
	ip(axis)=min(cg[grid].indexRange()(End,axis),max(cg[grid].indexRange()(Start,axis),ip(axis)));
	dr(axis)=r(0,axis)/cg[grid].gridSpacing()(axis)+cg[grid].indexRange()(0,axis)-ip(axis)-shift;
      }
      
      dra(Axes)=min(fabs(dr(Axes)),1.);
      //...........only use 4 points if dra bigger than epsilon, otherwise just use 2 points,
      //    this lets us  interpolate near interpolation boundaries
      ip1(Axes)=ip(Axes);
      for( axis=0; axis<cg.numberOfDimensions(); axis++ )
      {
        if( dra(axis)>epsi )
          ip1(axis)+= dr(axis)>0. ? 1 : -1;
        if( cg[grid].isPeriodic()(axis) )    // ........periodic wrap
	  ip1(axis)=MODR(ip1(axis),axis,grid);
      }
      if( interpDebug & 4 )
      {
	ip.display("here is ip");
	ip1.display("Here is ip1");
	dr.display("here is dr");
      }
      
      //.............Unable to interpolate if outside the current grid, but
      //             extrapolate (to zero order) if this is the closest point
      //             so far
      if( cg[grid].mask()(ip(0),ip(1),ip(2))!=0  &&
          (ip1(axis1)<indexRange(Start,axis1) || ip1(axis1)>indexRange(End,axis1) ||
	   ip1(axis2)<indexRange(Start,axis2) || ip1(axis2)>indexRange(End,axis2) ||
	   ip1(axis3)<indexRange(Start,axis3) || ip1(axis3)>indexRange(End,axis3)) )
      {
        // distance= distance between point and center(ip)
	distance= cg.numberOfDimensions()==2 ? 
	     fabs(x(0,0)-cg[grid].center()(ip(0),ip(1),ip(2),0))
	    +fabs(x(0,1)-cg[grid].center()(ip(0),ip(1),ip(2),1))
	      :
	   fabs(x(0,0)-cg[grid].center()(ip(0),ip(1),ip(2),0))
	  +fabs(x(0,1)-cg[grid].center()(ip(0),ip(1),ip(2),1))
	  +fabs(x(0,2)-cg[grid].center()(ip(0),ip(1),ip(2),2));
		  
        if( distance<minimumDistance )
	{
          extrap=TRUE;
          for( axis=0; axis<cg.numberOfDimensions(); axis++ )
            if(ip1(axis)<indexRange(Start,axis) || ip1(axis)>indexRange(End,axis))
	      ip1(axis)=ip(axis);
	}
        else
         continue;    //  don't extrapolate using this point, already have a better guess
      }

      if( cg.numberOfDimensions()==2 )
      {
	//  ... (check to see whether all marked interpolation points are valid)...
/* ----
	if( (int)cg[grid].isAllVertexCentered() )
	{
	  if(cg[grid].mask()(ip(0) ,ip(1))==0 || cg[grid].mask()(ip(0) ,ip1(1))==0  ||
	     cg[grid].mask()(ip1(0),ip(1))==0 || cg[grid].mask()(ip1(0),ip1(1))==0 )
	  {
            if( interpDebug & 8 )
  	      cout << "*** unable to interpolate on this grid, mask=0\n";
	    continue ;  //       ....Unable to interpolate, try another grid
	  }
	}
	else if( cg[grid].mask()(ip(0) ,ip(1))==0 )  // check this for cell centered
	  continue;
--- */
	if(cg[grid].mask()(ip(0) ,ip(1))==0 || cg[grid].mask()(ip(0) ,ip1(1))==0  ||
	   cg[grid].mask()(ip1(0),ip(1))==0 || cg[grid].mask()(ip1(0),ip1(1))==0 )
	{
	  if( interpDebug & 8 )
	    cout << "*** unable to interpolate on this grid, mask=0\n";
	  continue ;  //       ....Unable to interpolate, try another grid
	}

        if( extrap )  // if we can extrapolate then keep track of the distance of extrapolation
          minimumDistance=distance;

	// ...........Bi-Linear Interpolation:
        if( u.positionOfComponent(0)==3 )
	{
	  if( interpDebug &4) printf("\n\n");
	  for( int c1=R[1].getBase(); c1<=R[1].getBound(); c1++)  // add more components *******
	  {
	    for( int c0=R[0].getBase(); c0<=R[0].getBound(); c0++)
	    {
	      uInterpolated(ipt,c0,c1)= 
		(1.-dra(1))*(
			     (1.-dra(0))*u[grid]( ip(0),ip(1),ip(2),c0)
			        +dra(0) *u[grid](ip1(0),ip(1),ip(2),c0))
		  + dra(1) *(
			     (1.-dra(0))*u[grid]( ip(0),ip1(1),ip(2),c0)
			        +dra(0)* u[grid](ip1(0),ip1(1),ip(2),c0));

	      if( interpDebug &4) printf(">> value of u at interpolated point = %16.8f\n",uInterpolated(ipt,c0,c1));
	    }
	  }
	  if( interpDebug &4) printf("\n\n");
	}
	else
	{
	  // *** All this junk is so that we can interpolate the general case ****
	  int i00[8], i01[8], i10[8], i11[8];
	  i00[u.positionOfCoordinate(0)]=ip(0);   i10[u.positionOfCoordinate(0)]=ip1(0);
	  i00[u.positionOfCoordinate(1)]=ip(1);   i10[u.positionOfCoordinate(1)]= ip(1);
	  i00[u.positionOfCoordinate(2)]=ip(2);   i10[u.positionOfCoordinate(2)]=ip(2);

	  i01[u.positionOfCoordinate(0)]= ip(0);  i11[u.positionOfCoordinate(0)]=ip1(0);
	  i01[u.positionOfCoordinate(1)]=ip1(1);  i11[u.positionOfCoordinate(1)]=ip1(1);
	  i01[u.positionOfCoordinate(2)]=ip(2);   i11[u.positionOfCoordinate(2)]=ip(2);

	  for( int c1=R[1].getBase(); c1<=R[1].getBound(); c1++)  // add more components *******
	  {
	    i00[u.positionOfComponent(1)]=c1;
	    i10[u.positionOfComponent(1)]=c1;
	    i01[u.positionOfComponent(1)]=c1;
	    i11[u.positionOfComponent(1)]=c1;
	    for( int c0=R[0].getBase(); c0<=R[0].getBound(); c0++)
	    {
	      i00[u.positionOfComponent(0)]=c0;
	      i10[u.positionOfComponent(0)]=c0;
	      i01[u.positionOfComponent(0)]=c0;
	      i11[u.positionOfComponent(0)]=c0;
	      uInterpolated(ipt,c0,c1)= 
		(1.-dra(1))*(
			     (1.-dra(0))*u[grid](i00[0],i00[1],i00[2],i00[3])
			     +dra(0)*u[grid](i10[0],i10[1],i10[2],i10[3]))
		  + dra(1) *(
			     (1.-dra(0))*u[grid](i01[0],i01[1],i01[2],i01[3])
			     +dra(0)*u[grid](i11[0],i11[1],i11[2],i11[3]));
	    }
	  }
	}
      }
      else
      {
//	if( (int)cg[grid].isAllVertexCentered() )
//	{
	  if(cg[grid].mask()( ip(0),ip(1), ip(2))==0 || cg[grid].mask()( ip(0),ip1(1), ip(2))==0  ||
	     cg[grid].mask()(ip1(0),ip(1), ip(2))==0 || cg[grid].mask()(ip1(0),ip1(1), ip(2))==0  ||
	     cg[grid].mask()( ip(0),ip(1),ip1(2))==0 || cg[grid].mask()( ip(0),ip1(1),ip1(2))==0  ||
	     cg[grid].mask()(ip1(0),ip(1),ip1(2))==0 || cg[grid].mask()(ip1(0),ip1(1),ip1(2))==0  )
	  {
            if( interpDebug & 8 )
  	      cout << "*** unable to interpolate on this grid, mask=0\n";
	    continue ;  //       ....Unable to interpolate, try another grid
	  }
//	}
//	else if( cg[grid].mask()(ip(0),ip(1),ip(2))==0 )  // check this for cell centered
//	  continue;
        if( u.positionOfComponent(0)==3 )
	{
	  for( int c1=R[1].getBase(); c1<=R[1].getBound(); c1++)  // add more components *******
	  {
	    for( int c0=R[0].getBase(); c0<=R[0].getBound(); c0++)
	    {
	      uInterpolated(ipt,c0,c1)= 
               (1-dra(2))*(
		(1.-dra(1))*(
			     (1.-dra(0))*u[grid](ip(0),ip(1),ip(2),c0)
			         +dra(0)*u[grid](ip1(0),ip(1),ip(2),c0))
		  + dra(1) *(
			     (1.-dra(0))*u[grid](ip(0),ip1(1),ip(2),c0)
			         +dra(0)*u[grid](ip1(0),ip1(1),ip(2),c0))
                          )
                 +dra(2)*(
		(1.-dra(1))*(
			     (1.-dra(0))*u[grid](ip(0),ip(1),ip1(2),c0)
			         +dra(0)*u[grid](ip1(0),ip(1),ip1(2),c0))
		  + dra(1) *(
			     (1.-dra(0))*u[grid](ip(0),ip1(1),ip1(2),c0)
			         +dra(0)*u[grid](ip1(0),ip1(1),ip1(2),c0))
                          );
	    }
	  }
	}
	else
	{
          cout << "interpolatePoints:error: positionOfComponent=" << u.positionOfComponent(0)
	    << " not implemented yet. Talk to WDH \n";
	  return 1;
	}
      }
      // return the values used:
      if( indexGuess.GET_LENGTH(0) >= numberOfPointsToInterpolate )
        for( axis=0; axis<cg.numberOfDimensions(); axis++ )
          indexGuess(ipt,axis)=ip(axis);

      if( interpoleeGrid.GET_LENGTH(0) >= numberOfPointsToInterpolate )
        interpoleeGrid(ipt)=grid;   //  !extrap ? k: -k;

      if( wasInterpolated.GET_LENGTH(0) >= numberOfPointsToInterpolate )
        wasInterpolated(ipt)=!extrap;
      
      if( extrap )
	pointWasExtrapolated=1;
      else
      {
        pointWasExtrapolated=0;
        break;   // point has been successfully interpolated, try next point
      }
    }

    if( pointWasExtrapolated==-1 )
    {
      cout << "interpolatePoints::ERROR: unable to interpolate or extrapolate a point! \n";
      cout << "This error should not occur\n";
      printf(" ipt=%i, x=(%f,%f,%f) \n",ipt,x(0,0),x(0,1),x(0,2));
    }
    else if( pointWasExtrapolated==1 )
      returnValue--;   // count number of extrapolated points
  } // end for ipt... over the interpolated points
  if( interpDebug & 4 )
    printf(" Closest point = (%i,%i), grid = %i \n",ip(0),ip(1),grid);
    

  return returnValue;
}

#undef NRM
#undef MODR
#undef GET_LENGTH
