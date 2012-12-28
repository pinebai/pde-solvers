//
// -- representing the nucleus in CellWave
//

#include <stdio.h>
#include "Nucleus.h"
#include "GenericReactionMacros.h"

using namespace CellWave;

Nucleus::Nucleus()
{
  nucleusShape=NoNucleus; //default =no nucleus
  x0=0., y0=0., z0=0.;
  radius=10.;
  x1=10., y1=10., z1=10.;
  boundaryThickness =1e-5; //boundaryThickness =0.; //effectively zero thickness
}

Nucleus::
~Nucleus()
{
  //..do nothing
}

void Nucleus::
setup()
  // setup: precompute necessary quantities for fast evaluation
  //        of the mask
{
  //.. ? precompute stuff for the box version
}

void Nucleus::
setID( int id_ )
{
  id = id_;
}

void Nucleus::
setCenter( double x, double y, double z)
{
  x0=x, y0=y, z0=z;
}

void Nucleus::
setRadius( double radius_ )
{
  radius = radius_;
} 

void Nucleus::
setBoundaryThickness( double thickness_ )
{
  boundaryThickness = thickness_;
}

void Nucleus::
setCorners( double x0_, double y0_, double z0_,
	    double x1_, double y1_, double z1_)
{
  x0=x0_, y0=y0_, z0=z0_;
  x1=x1_, y1=y1_, z1=z1_;
}  

int Nucleus::
getID() const
{
  return(id);
}

void Nucleus::
getCenter( double &x, double &y, double &z) const
{
  x=x0, y=y0, z=z0;
}

double Nucleus::
getRadius() const
{
  return(radius);
} 

double Nucleus::
getBoundaryThickness() const
{
  return( boundaryThickness );
}

void Nucleus::
getCorners( double &x0_, double &y0_, double &z0_,
	    double &x1_, double &y1_, double &z1_) const 
{
  x0_=x0, y0_=y0, z0_=z0;
  x1_=x1, y1_=y1, z1_=z1;
}  

std::string Nucleus::
getNucleusShapeName()
{
  std::string name="unknown";
  if      ( nucleusShape == NoNucleus ) {
    name="No Nucleus";
  }
  else if ( nucleusShape == SphericalNucleus ) {
    name="Spherical Nucleus";
  }
  else if ( nucleusShape == BoxNucleus ) {
    name="Box Nucleus";
  }

  return( name );

}

double Nucleus::
getMask(       double x, double y, double z)
{
  double mask=1.;
  if ( nucleusShape == NoNucleus ) {
    mask=1.;
  }
  else if (nucleusShape == SphericalNucleus ) {
    mask = getSphericalMask( x,y,z );
  }
  else if (nucleusShape == BoxNucleus ) {
    mask = getBoxMask( x, y, z );
  }
  else {
    printf("--ERROR -- unknown nucleus shape %d in cell %d\n",
	   int( nucleusShape ), getID());
  }
}

void Nucleus::
getMaskArray( const double &tcomp,
                   const int &nd,    const int &ncomp,
	  	   const int &nd1a, const int &nd1b,
		   const int &nd2a, const int &nd2b,
		   const int &nd3a, const int &nd3b,
		   const int &n1a,  const int &n1b,
		   const int &n2a,  const int &n2b,
		   const int &n3a,  const int &n3b,
		   const double &xfirst, double &qfirst )
{
  const double *xArray = &xfirst;
  double       *qArray = &qfirst;
  const int xaxis=0, yaxis=1, zaxis=2;
  const int imask=0;

  // nd     = number of dimensions
  // ncomp  = number of components
  // nd?a, nd?b for ?=1,2,3 gives _array_ dimensions
  // n?a, n?b   for ?=1,2,3 gives _loop bounds_
  // ForGrid,GRIDINDEX:  defined in GenericReactionImpl.h

  //  printf("ReactionSlepchenko2Buffer--initial data %d\n", 
  //	 int(ip3Distribution));

  if ( nucleusShape == NoNucleus ) {
    ForGrid(i,j,k) {
      qArray[ GRIDINDEX( i,j,k, imask ) ] = 1.;
    }    
  }
  else if( nd ==3 ) {
    if ( nucleusShape == SphericalNucleus ) {
      ForGrid(i,j,k) {
	
	const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
	const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
	const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
	
	double &mask      = qArray[ GRIDINDEX( i,j,k, imask) ];
	mask = getSphericalMask( x,y,z);
      }
    }
    else if ( nucleusShape == BoxNucleus ) {
      ForGrid(i,j,k) {
	
	const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
	const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
	const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
	
	double &mask      = qArray[ GRIDINDEX( i,j,k, imask) ];
	mask = getBoxMask( x,y,z);
      }
    }
  }
  else  if (nd == 2) {
    if( nucleusShape == SphericalNucleus ) {
      ForGrid(i,j,k) {
	const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
	const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
	const double z= 0.; // z must be zero 2D initial data
	
	double &mask      = qArray[ GRIDINDEX( i,j,k, imask) ];
	mask = getSphericalMask( x,y,z); 
      }
    }
    else if ( nucleusShape == BoxNucleus ) {
      ForGrid(i,j,k) {
	const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
	const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
	const double z= 0.; // z must be zero 2D initial data
	
	double &mask      = qArray[ GRIDINDEX( i,j,k, imask) ];
	mask = getBoxMask( x,y,z); 
      }      
    }
  }
}
  
