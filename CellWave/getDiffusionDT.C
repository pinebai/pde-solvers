//
// max timestep for a multicomponent diffusion eqn. 
//
#include <assert.h>

#include "Overture.h"  
#include "MappedGridOperators.h"

#include "getDiffusionDT.h"

real
getDiffusionDT(const real & cfl, 
	       int numberOfComponents,
	       realArray &nuArray,
	       CompositeGrid & cg,
	       const real alpha0, /*= -2.,*/
	       const real beta0  /*= 1. */ )
{
  realArray dtArray;
  getDiffusionDT( cfl, numberOfComponents, nuArray, dtArray,
		  cg, alpha0, beta0 );
  return( min(dtArray) );
}

void
getDiffusionDT(const real & cfl, 
	       int numberOfComponents,
	       realArray &nuArray,
	       realArray &dtArray,
	       CompositeGrid & cg,
	       const real alpha0, /*= -2.,*/
	       const real beta0  /* = 1. */)
{
  realArray dtGridArray;
  getDiffusionDT( cfl, numberOfComponents, nuArray, dtArray,
		  dtGridArray, cg, alpha0, beta0 );
  //  dtGridArray.display("dtGridArray");
}


void
getDiffusionDT(const real & cfl, 
	       int numberOfComponents,
	       realArray &nuArray,
	       realArray &dtArray,
	       realArray &dtGridArray,
	       CompositeGrid & cg,
	       const real alpha0, /*= -2.,*/
	       const real beta0  /*= 1.*/ )
//======================================================================================
// /Description:
//  Determine the time step for the convection diffusion equation
//      u_t + a u_x + b u_y = nu( u_xx + u_yy )
//  discretized with the mapping method. Scale the maximum allowable time
//  step by the factor cfl. the stability region is assumed to lie within the 
//  ellipse (x/alpha0)^2 + (y/beta0)^2 = 1
// /cfl (input): Scale the time step by this factor (cfl=1 should normally be stable)
// /a (input) : coefficient of u_x
// /b (input) : coefficient of u_y
// /nu (input) : coefficient of (u_xx+u_yy)
// /alpha0, beta0 (input) : parameters defining the ellipse for the stability region  
// ====================================================================================
{
  dtArray.redim(numberOfComponents);
  dtGridArray.redim( cg.numberOfGrids(), numberOfComponents);

  dtArray= 1e100;      // no timestep limit
  dtGridArray= 1e100;  // no timestep limit

  // .. stability is evaluated from 
  //     nu*laplace = purely real (negative) eigenvalues
  for (int ig=0; ig< cg.numberOfGrids(); ++ig) {
    MappedGrid             &mg   = cg[ig];
    MappedGridOperators    op(mg);

    for (int ic=0; ic < numberOfComponents; ++ic) { // loop through different diffusions
      const real nu = nuArray(ic);
      if ( fabs(nu) < 1e-16 ) {
	//dtArray(ic)= 1e100; // no timestep limit
	break;
      }
      real dt;
      if( mg.isRectangular() )  {
	// ***** rectangular grid *****
	real dx[3];
	mg.getDeltaX(dx);
	dt = cfl / fabs( nu *(4./(alpha0*dx[0]*dx[0])+4./(alpha0*dx[1]*dx[1])));
	//dt = cfl * pow( pow( fabs(a)*(1./(beta0*dx[0]))+fabs(b)*(1./beta0*dx[1]) , 2.)
	//	     +pow( nu *(4./(alpha0*dx[0]*dx[0])+4./(alpha0*dx[1]*dx[1])) , 2.)
	//      ,-.5);
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
      if (ig>0) {
	dtArray(ic) = min( dt, dtArray(ic));
      }
      else {
	dtArray(ic) = dt;
      }
      dtGridArray( ig, ic ) = dt;
      //printf(".. grid %5d, getDiffusionDT;  dt=  %24.16e\n",ig, dt );
    } //end for ic (components)
  }//end for ig (grids)
}
