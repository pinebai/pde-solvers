// Brief description: base class for Ca/IP3 dynamics
///
///   CellWave::GenericCalciumIP3Reaction  -- base class for Ca2+/IP3 dynamics
///     part of CellWave
///
#ifndef CELLWAVE_GENERIC_CALCIUM_IP3_REACTION_H
#define CELLWAVE_GENERIC_CALCIUM_IP3_REACTION_H

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "GenericReaction.h"
#include "CellWave.h"

namespace CellWave {

  class GenericCalciumIP3Reaction : public GenericReaction {
  public:
    enum { cc=0, pc=1, hc=2, numberOfComponents};
    
    GenericCalciumIP3Reaction() 
      { //default
      };
    ~GenericCalciumIP3Reaction()
      { //default
      };

  private: // do not use copy constructors for now
    GenericCalciumIP3Reaction( const GenericCalciumIP3Reaction &X) 
      { };
    GenericCalciumIP3Reaction & operator=( const GenericCalciumIP3Reaction &X)
      { };
  public:

    double  getMaximumLengthScale() { return maxLengthScale; };
    void    setMaximumLengthScale( double len) { maxLengthScale = len; };

    //.. set parameters
    /* virtual */ bool
      readParameterFile( CellWave::ParameterReader &param);

    int getNumberOfDimensions( ) { return numberOfDimensions; };

    virtual void printParameters(int iPrintChannel);

    //.. looping calls to array interface, for initial data & rhs
    /* virtual */ void 
       callInitialDataLoop( const double &tcomp,
			   const int&nd,    const int &ncomp,
			   const int &nd1a, const int &nd1b,
			   const int &nd2a, const int &nd2b,
			   const int &nd3a, const int &nd3b,
			   const int &n1a,  const int &n1b,
			   const int &n2a,  const int &n2b,
			   const int &n3a,  const int &n3b,
			   const double &x, double &q );

    inline void setInitialData( const double &x, 
				const double &y, 
				const double &z,
				double &c,       
				double &h,       
				double &p);

    void initializeParameters();

    //.. data
  private:
    double maxLengthScale;
  public:
    ///..initial data values
    double calcium_0;   /// initial calcium2+ level
    double ip3_0;       /// initial IP3 level
    double h_0;         /// initial inhibitor level
    
    int    numberOfDimensions; /// used for spatially dependent initial data

    ///..initial blob data
    double xBlob_c;
    double yBlob_c;
    double zBlob_c;
    double blobWidth_c;
    double blobMin_c, blobMax_c;
    
    double xBlob_p;
    double yBlob_p;
    double zBlob_p;
    double blobWidth_p;
    double blobMin_p, blobMax_p;

    //..data for Exact Gaussian initial data
    double gaussianTotalConcentration;
    double gaussianTimeOffset;
    //  --> center is given by xBlob_p, yBlob_p, zBlob_p

    /// for inhomog distribution of IP3, Wagner et al, figure 4.
    enum IP3DistributionType {HomogeneousIP3=0, BlobIP3=1, BoxBlobIP3=2, 
			      RadialIP3=3, RadialWithGradientIP3=4, ExactGaussianIP3=5};
    enum CaDistributionType {HomogeneousCa=0, BlobCa=1 };
    
    IP3DistributionType ip3Distribution;
    CaDistributionType  caDistribution;
    double ip3Dist_I_s, ip3Dist_I_h, ip3Dist_r_c, ip3Dist_I_w; 
    double ip3Dist_xstar, ip3Dist_Iprime_h, ip3Dist_Iprime_w;

    double ip3scale;
    double ip3diffusion;
    
    //..box IP3 data
    double ip3BoxMax;
    double ip3Box_x0, ip3Box_y0, ip3Box_z0;
    double ip3Box_x1, ip3Box_y1, ip3Box_z1;

    void setBoxCorners( double x0, double y0, double z0, 
			double x1, double y1, double z1) 
      {
	ip3Box_x0= x0; ip3Box_y0= y0; ip3Box_z0= z0;
	ip3Box_x1= x1; ip3Box_y1= y1; ip3Box_z1= z1;
	DPrintf(DebugReaction,"Calcium Rxn, corners:\n");
	DPrintf(DebugReaction," (x0, y0, z0) = ( %f, %f, %f )\n",
		ip3Box_x0, ip3Box_y0, ip3Box_z0,
		ip3Box_x1, ip3Box_y1, ip3Box_z1 );

	DPrintf(DebugReaction," (x1, y1, z1) = ( %f, %f, %f )\n",
		ip3Box_x1, ip3Box_y1, ip3Box_z1 );
      };
    
  };//end class GenericCalciumIP3Reaction


  ///
  /// ..INLINE functions
  ///
  
  /// -- initial data ----------------------------------------
  
  inline void GenericCalciumIP3Reaction::
  setInitialData( const double &xin, 
		  const double &yin, 
		  const double &zin,
		  double &c,       
		  double &h,       
		  double &p)
    {
      double z =zin;
      double y =yin;
      double x =xin;
      if (numberOfDimensions <=2) {
	z=0.;
      }
      if (numberOfDimensions <=1) {
	y=0.;
      }

      //..H initial data
      h = h_0;

      //..Ca2+ initial data
      c = calcium_0;   ///..default

      if ( caDistribution == HomogeneousCa ) {
	c = calcium_0;   ///.. initialize to homogeneous c
      }   
      else if ( caDistribution == BlobCa ) {
	const double r2 =  pow( x-xBlob_c,2.)+pow( y- yBlob_c,2.)+pow( z-zBlob_c,2.);
	c = blobMin_c + ( blobMax_c -  blobMin_c)
	* exp(- r2/pow( blobWidth_c,2.) );
      }

      //..IP3 initial data
      if ( ip3Distribution == HomogeneousIP3 ) {
	p = ip3_0;
      }
      else if ( ip3Distribution == BlobIP3 ) { 
	const double r2 =  pow( x-xBlob_p,2.)+pow( y- yBlob_p,2.)+pow( z-zBlob_p,2.);
	p  = blobMin_p + ( blobMax_p - blobMin_p) * exp(- r2/pow( blobWidth_p,2.) );
      }
      else if ( ip3Distribution == BoxBlobIP3 ) { 
	// is x,y,z inside the box?
	double mask;
	bool isInside;
	isInside = ( ip3Box_x0<= x ) && (x <= ip3Box_x1 );
	isInside = isInside && ( (ip3Box_y0 <= y ) && (y <= ip3Box_y1 ));
	isInside = isInside && ( (ip3Box_z0 <= z ) && (z <= ip3Box_z1 ));
	DPrintf(DebugReaction,"x"); 
	if (isInside) {	  p= ip3BoxMax; }
	else          {	  p= ip3_0;	}
      }
      else if ( ( ip3Distribution == RadialIP3)
		||  (ip3Distribution == RadialWithGradientIP3)) {
	///fig 4, Wagner et al
	const double r =  sqrt( pow( x,2.)+pow( y,2.)+pow( z,2.) );
      	p= ip3Dist_I_s*( 1.+ip3Dist_I_h * exp( ( r/ ip3Dist_r_c    -1.) / ip3Dist_I_w));
      } 
      else if ( (ip3Distribution == ExactGaussianIP3) && (ip3diffusion>0) ) {
	const double pi       =  4.*atan(1.);
	const double halfNdim = double( numberOfDimensions )/2.;
	const double width    =  4.*pi*ip3diffusion*(-gaussianTimeOffset);
	const double r2       =  pow( x-xBlob_p,2.)+pow( y- yBlob_p,2.)+pow( z-zBlob_p,2.);
	p  = gaussianTotalConcentration*pow(width,-halfNdim)* exp(- r2/width );
      }
      else throw "error: unknown IP3 distribution type\n";
      
      if ( ip3Distribution == RadialWithGradientIP3) { 
	///fig 4, Wagner et al, second equation
	if ( x + ip3Dist_xstar <0 ) {
	  
	  const double r =  sqrt( pow( x,2.)+pow( y,2.)+pow( z,2.) ); /// note: z=0 in 2D
	  p += ip3Dist_Iprime_h *  exp( (r/ip3Dist_r_c -1.)/ip3Dist_I_w)
	    *  exp( - pow( y / ip3Dist_r_c* ip3Dist_Iprime_w , 4.));
	}///end if x + x* <0

      }///end if ip3Distribution == RadialWithGrad...
      //printf(" ( %e, %e, %e ) ", c,p,h);
      if ( p >100. ) {
	p=100.; // clip off at [IP3] =100    //  FIXME:   kludge for max IP3
      }
      p = ip3scale * p;
      //printf("   IP3:  x=%8.4g  y=%8.4g  z=%8.4g,  p= %8.4g\n", x,y,z,p);
    }
  

}; //end namespace CellWave
#endif
