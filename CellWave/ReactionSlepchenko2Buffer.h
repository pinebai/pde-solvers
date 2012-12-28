/// Brief description: Li-Rinzel rxn mechanism, as used by Wagner et al
///
///   CellWave::ReactionSlepchenko2Buffer  implements
///     the Li-Rinzel model with two Ca buffers (mobile=b1, and immobile=b2)
///     see Slepchenko, Schaff, and Choi, J. Comp. Phys. (2000), vol 162, p. 186-218
///   Eq. (5.1), p. 210 & Eqs. (4.2)--(4.3), p. 204, and IP3 diffusion + degradation
///
///     but currently no IP3 effects in Ca2+ dynamics. (J_0 is indep of IP3)
///
///      c_t   = D_c Laplace( c )   + f( c, h) + R1 + R2
///      h_t   = g( c,h)
///      p_t   = D_p Laplace( p )   - k_i* p
///      b1_t  = -R1
///      b2_t  = D_b1 Laplace( b2 )  - R2
///
///   where from the Li-Rinzel model
///     f(c,h) = J_0 ( c*h / ( c+d_act))^3 - V_m * c^2/( c^2 + K_m^2) + L
///     g(c,h) = k_on*( d_inh - (d_inh +c)*h),
///
///   and we have for the immobile Ca2+ buffer
///     R1 = - k1_on* (b1_tot - b1)*c + k1_off*b1 
///   and for the mobile Ca2+ buffer
///     R2 = - k2_on* (b2_tot - b2)*c + k2_off*b2
///

#ifndef REACTION_SLEPCHENKO2BUFFER_H 
#define REACTION_SLEPCHENKO2BUFFER_H "ReactionSlepchenko2Buffer.h"

#include <math.h>
#include <stdio.h>
#include <string>
#include <iostream>

#include "GenericCalciumIP3Reaction.h"
#include "ParameterReader.h"

namespace CellWave {
  
  class ReactionSlepchenko2Buffer : public GenericCalciumIP3Reaction {
  public:
    enum { cc=0, pc=1, hc=2,
	   b1c= 3, b2c=4, numberOfComponents};

    ReactionSlepchenko2Buffer();
    ~ReactionSlepchenko2Buffer();

  public: 
    static std::string rxnType() { return "Slepchenko2Buffer"; }
    void initializeParameters();

    virtual bool
      readParameterFile( CellWave::ParameterReader &param);

    virtual double getDiffusionCoefficient(const int component); 
    virtual bool   isDiffusive( const int component);

    virtual std::string getTitle();
    virtual std::string getLongComponentName( const int component );
    virtual std::string getShortComponentName( const int component );

    virtual int    getNumberOfSpecies();
    //virtual double  getLengthScale() { return maxLengthScale; };

    ///.. array loop interface, for initial data & rhs
    virtual void callInitialDataLoop( const double &tcomp,
                   const int&nd,    const int &ncomp,
	  	   const int &nd1a, const int &nd1b,
		   const int &nd2a, const int &nd2b,
		   const int &nd3a, const int &nd3b,
		   const int &n1a,  const int &n1b,
		   const int &n2a,  const int &n2b,
		   const int &n3a,  const int &n3b,
		   const double &x, double &q );

    virtual void callRHSLoop( const double &tcomp,
		   const int&nd,    const int &ncomp,
		   const int &nd1a, const int &nd1b,
		   const int &nd2a, const int &nd2b,
		   const int &nd3a, const int &nd3b,
		   const int &n1a,  const int &n1b,
		   const int &n2a,  const int &n2b,
		   const int &n3a,  const int &n3b,
		   const double &x, const double &q, 
		   double &rhs);

    /// Inline function:
    ///.. scalar calls to initial data, and rhs
    inline void setInitialData( const double &x, 
				const double &y, 
				const double &z,
				double &c,       
				double &h,       
				double &p,
				double &b1,
				double &b2);
    
    inline void computeRHS( const double c,  /// IN:       
			    const double h,       
			    const double p,
			    const double b1,
			    const double b2,
			    double & rhs_c,  /// OUT:
			    double & rhs_h,
			    double & rhs_p,
			    double & rhs_b1,
			    double & rhs_b2,
			    ///                 OPTIONAL input
			    const double x=0.,  
			    const double y=0., 
			    const double z=0. );


    //..DATA..................................................................
    //...now in GenericCalciumIP3 Rxn
    //double maxLengthScale; // for getting an estimate of domain size

    //..PARAMETERS............................................................
    ///....components
    //enum { cc=0, pc=1, hc=2,
    //       b1c=3, b2c=4, numberOfComponents};

    ///..parameters, as Table III, p. 210, in Slepchenko (see also text on p.210)
    double J_0;
    double d_act;
    double d_inh;
    double k_on;
    double V_m;
    double K_p;
    double leak;   // = L in Table III, Slepchenko et al.
    double b1_tot;
    double b2_tot;
    double K1;
    double K2;
    double k1_on, k1_off;
    double k2_on, k2_off;

    //..diffusions: these are D_c, D_p, and D_b1, respectively
    double diffCalcium;
    double diffIP3;
    double diffB2;
    
    double k_i; // IP3 degradation

    //..disable buffers at will
    int b1Mask;
    int b2Mask;

    ///..initial data values
    double b1_0;        /// initial immobile buffer level
    double b2_0;        /// initial mobile buffer level

    //-------generic initial data for c,p,h --> in GenericCa2+/IP3 rxn class
#if 0    
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
    
    
    /// for inhomog distribution of IP3, Wagner et al, figure 4.
    //enum IP3DistributionType {HomogeneousIP3=0, BlobIP3=1, 
    //RadialIP3=2, RadialWithGradientIP3=3};
    //enum CaDistributionType {HomogeneousCa=0, BlobCa=1 };
    //IP3DistributionType ip3Distribution;
    //CaDistributionType  caDistribution;
    double ip3Dist_I_s, ip3Dist_I_h, ip3Dist_r_c, ip3Dist_I_w; 
    double ip3Dist_xstar, ip3Dist_Iprime_h, ip3Dist_Iprime_w;
#endif
    
  }; // end class ReactionSlepchenko2Buffer
  
  ///
  /// ..INLINE functions
  ///
  
  /// -- initial data ----------------------------------------
  
  inline void ReactionSlepchenko2Buffer::
  setInitialData( const double &xin, 
		  const double &yin, 
		  const double &zin,
		  double &c,       
		  double &h,       
		  double &p,
		  double &b1,
		  double &b2)
    {
      /// call slow version, 'if' inside innermost loop:
      ///.. if _slow version is a problem, write one outer loop
      ///   for each type of initial data and call only one
      ///   kind of initial data in each innermost loop. not worth the trouble currently.
      double z =zin;
      double y =yin;
      double x =xin;
      if (numberOfDimensions <=2) {
	z=0.;
      }
      if (numberOfDimensions <=1) {
	y=0.;
      }

      GenericCalciumIP3Reaction::setInitialData(  x,y,z, c,h,p );

      //..buffer initial data
      b1= b1_0;
      b2= b2_0;

    }
  
  /// -- reaction RHS evaluator ------------------------------
  
  inline void ReactionSlepchenko2Buffer::
  computeRHS( const double c,  /// IN:       
	      const double h,       
	      const double p,
	      const double b1,
	      const double b2,
	      double & rhs_c,  /// OUT:
	      double & rhs_h,
	      double & rhs_p,
	      double & rhs_b1,
	      double & rhs_b2,
	      ///                 OPTIONAL input
	      const double x=0., /// can omit x,y,z if not needed for reaction terms 
	      const double y=0., 
	      const double z=0. )
  {
    double JFlux = J_0*pow( c*h/( c+d_act), 3.); 
    double JSerca= -V_m* pow(c,2.)/( pow(c,2.) + pow( K_p, 2.));
    double R1 = b1Mask*( -k1_on*( b1_tot - b1)*c +k1_off*c);
    double R2 = b2Mask*(-k2_on*( b2_tot - b2)*c +k2_off*c);
      
    rhs_c  =  JFlux + JSerca + leak + R1 + R2;       // Ca2+
    rhs_p  = - k_i *p;                               // IP3
    rhs_h  =  k_on*( d_inh - (c+ d_inh)* h);         // gating variable  

    rhs_b1 = -R1;
    rhs_b2 = -R2;
  }

}; //end namespace CellWave

#endif

