//
// max timestep for a multicomponent diffusion eqn. 
//
#include <assert.h>

#include "Overture.h"  
#include "MappedGridOperators.h"

#include "getDiffusionDT.h"

void
getSpacing(
	   realArray &minSpacing,
	   realArray &maxSpacing,
	   CompositeGrid & cg,
	   const real alpha0 = -2.,
	   const real beta0  = 1. )
//======================================================================================
// /Description:
//  compute the minimum & maximum grid spacing in x,y,z on each grid
//  /minSpacing( igrid, axis) (output) : min. spacing on grid=igrid & axis
//  /maxSpacing( igrid, axis) (output) : min. spacing on grid=igrid & axis
//  /cg (input) : composite grid
//  /alpha0, beta0 (input) : parameters defining the ellipse for the stability region  
// ====================================================================================
{
  enum {xaxis=0, yaxis=1, zaxis=2};
  const int maxDimensions =3;
  int numberOfDimensions=cg.numberOfDimensions();
  int numberOfGrids     =cg.numberOfGrids();
  minSpacing.redim( numberOfGrids, maxDimensions );
  maxSpacing.redim( numberOfGrids, maxDimensions );
   
  for (int ig=0; ig< cg.numberOfGrids(); ++ig) {
    MappedGrid             &mg   = cg[ig];
    MappedGridOperators    op(mg);

    if( mg.isRectangular() )  {
      // ***** rectangular grid *****
      real dx[3];
      mg.getDeltaX(dx);
      for (int iaxis=0; iaxis<maxDimensions; ++iaxis ) {
	minSpacing( ig, iaxis ) = dx[ iaxis ];
	maxSpacing( ig, iaxis ) = dx[ iaxis ];
      }
    }
    else   {
      // ***** non-rectangular grid *****
      mg.update(MappedGrid::THEinverseVertexDerivative);  // make sure the jacobian derivatives are built
      // define an alias:
      realMappedGridFunction & rx = mg.inverseVertexDerivative();
      // we need to compute the derivatives of rx:
      rx.setOperators(op);
      // Get Index's for the interior+boundary points
      Index I1,I2,I3;
      getIndex( mg.indexRange(),I1,I2,I3);
      
      realArray a1,b1,nu11,nu12,nu22;
      const real a=0., b=0.;
      // a1 = a*r1.x + b*r1.y + nu ( r1.xx + r1.yy )     b1 = a*r2.x + b*r2.y + nu ( r2.xx + r2.yy )
      a1   = a*rx(I1,I2,I3,0,0) + b*rx(I1,I2,I3,0,1)
	- nu*( rx.x(I1,I2,I3,0,0)(I1,I2,I3,0,0)  + rx.y(I1,I2,I3,0,1)(I1,I2,I3,0,1) );
      b1   = a*rx(I1,I2,I3,1,0) + b*rx(I1,I2,I3,1,1)
	- nu*( rx.x(I1,I2,I3,1,0)(I1,I2,I3,1,0)  + rx.y(I1,I2,I3,1,1)(I1,I2,I3,1,1) );
      // nu11 = nu*( r1.x*r1.x + r1.y*r1.y )
      // nu12 = nu*( r1.x*r2.x + r1.y*r2.y )*2 
      // nu22 = nu*( r2.x*r2.x + r2.y*r2.y ) 
      nu11 = nu*( rx(I1,I2,I3,0,0)*rx(I1,I2,I3,0,0) + rx(I1,I2,I3,0,1)*rx(I1,I2,I3,0,1) );
      nu12 = nu*( rx(I1,I2,I3,0,0)*rx(I1,I2,I3,1,0) + rx(I1,I2,I3,0,1)*rx(I1,I2,I3,1,1) )*2.;
      nu22 = nu*( rx(I1,I2,I3,1,0)*rx(I1,I2,I3,1,0) + rx(I1,I2,I3,1,1)*rx(I1,I2,I3,1,1) );
      
      // Grid spacings on unit square:
      real dr1 = mg.gridSpacing()(axis1);
      real dr2 = mg.gridSpacing()(axis2);
      dt = cfl * min( 
	    pow(
		pow( abs(a1)*(1./(beta0*dr1))+abs(b1)*(1./beta0*dr2) , 2.)
		+pow(   nu11 *(4./(alpha0*dr1*dr1)) 
			+abs(nu12)*(1./(alpha0*dr1*dr2))
			+nu22 *(4./(alpha0*dr2*dr2)) , 2.)
		,-.5) 
	    );
    } //end is.rectangular
  }//end for ig (grids)
}
