#ifndef NUCLEUS_H
#define NUCLEUS_H
//
// -- representing the nucleus in CellWave
//    -- one/mapped grid  (= one 'Nucleus' instance per
//       physical cell)

#include <math.h>
#include <string>

namespace CellWave {

class Nucleus 
{
 public:
  //..public data
  enum NucleusShape { NoNucleus=0, SphericalNucleus=1, BoxNucleus=2 };

  //..public functions
  Nucleus();
  ~Nucleus();

  void setup();

  void setID( int id_ );
  void setCenter( double x, double y, double z);
  void setRadius( double radius_ );
  void setBoundaryThickness( double thickness );
  void setCorners( double x0_, double y0_, double z0_,
		   double x1_, double y1_, double z1_);
  void setShape( NucleusShape shape) { nucleusShape = shape; };

  int  getID() const;
  void getCenter( double &x, double &y, double &z) const;
  double getRadius( ) const;
  double getBoundaryThickness() const;
  void getCorners( double &x0_, double &y0_, double &z0_,
		   double &x1_, double &y1_, double &z1_) const;
  NucleusShape getShape() { return nucleusShape; };
  
  
  double getMask(       double x, double y, double z);
  void   getMaskArray( const double &tcomp,
                   const int &nd,    const int &ncomp,
	  	   const int &nd1a, const int &nd1b,
		   const int &nd2a, const int &nd2b,
		   const int &nd3a, const int &nd3b,
		   const int &n1a,  const int &n1b,
		   const int &n2a,  const int &n2b,
		   const int &n3a,  const int &n3b,
		   const double &x, double &q );
  
  inline double getSphericalMask( double x, double y, double z ); // fast pointwise eval
  inline double getBoxMask( double x, double y, double z);        // fast pointwise eval

  std::string getNucleusShapeName();
  
 public: // ..data

  int id; // cell ID to which this nucleus belongs. Could be different from mappedGrid id
  NucleusShape nucleusShape;
  double boundaryThickness; // sharp=0., otherwise give lenght>0 for blending bdry
  double x0, y0, z0; 
  double x1, y1, z1; 

  //..spherical
  // center = x0,y0,z0
  double radius;

  //..box
  // bottom left  corner= x0, y0, z0
  // top    right corner= x1, y1, z1

};
 
//------inline fcns
inline double  Nucleus:: 
getSphericalMask( double x, double y, double z )
{
#define SQ(q) ( (q)*(q))
  double mask=0.;
  double r=sqrt( SQ( x-x0 ) + SQ( y-y0 ) + SQ( z-z0));
  mask = 0.5*( tanh(2.*( r-radius )/boundaryThickness ) + 1. );
#undef SQ  
  return( mask );
}

inline double  Nucleus::
getBoxMask( double x, double y, double z )
{
  // is x,y,z inside the box?
  double mask;
  bool isInside;
  isInside = ( x0<= x) && (x <=x1 );
  isInside = isInside && ( (y0 <= y ) && (y <= y1 ));
  isInside = isInside && ( (z0 <= z ) && (z <= z1 ));

  if (isInside) {
    mask=0.;
  }
  else {
    mask=1.;
  }
  return( mask );
}


} // end namespace CellWave
#endif

