/// Brief description: Li-Rinzel rxn mechanism, as used by Wagner et al
///
///   CellWave::ReactionLiRinzelWagner  implements
///    the Li-Rinzel mechanism as used by Wagner et al, Biophys. J. '98
///    Note that IP3 could diffuse, we have the diffusion coeff. here.
///    The Ca2+ buffering is done through _effective buffering_,
///    which assumes rapid equilibria for the mobile & immobile buffers,
///    and that the buffer conc. is small(?)
///
#ifndef CELLWAVE_REACTION_LI_RINZEL_WAGNER_H
#define CELLWAVE_REACTION_LI_RINZEL_WAGNER_H "ReactionLiRinzelWagner.h"

#include <math.h>
#include <stdio.h>
#include <string>
#include <iostream>

#include "GenericCalciumIP3Reaction.h"
//#include "GenericReaction.h"
#include "ParameterReader.h"

#include "CellWave.h"

namespace CellWave {
  
  class ReactionLiRinzelWagner : public GenericCalciumIP3Reaction {
  public:
    enum { cc=0, pc=1, hc=2,
	   bcfast=3, bcslow=4, numberOfComponents};

    ReactionLiRinzelWagner();
    ~ReactionLiRinzelWagner();

  public: 
    static std::string rxnType() { return "LiRinzelWagner"; }
    void initializeParameters();

    /* virtual */ bool
      readParameterFile( CellWave::ParameterReader &param);

    /* virtual */ double getDiffusionCoefficient(const int component); 
    /* virtual */ bool   isDiffusive( const int component);

    /* virtual */ std::string getTitle();
    /* virtual */ std::string getLongComponentName( const int component );
    /* virtual */ std::string getShortComponentName( const int component );

    /* virtual */ int    getNumberOfSpecies();

    ///.. array loop interface, for initial data & rhs
    /* virtual */ void callInitialDataLoop( const double &tcomp,
                   const int&nd,    const int &ncomp,
	  	   const int &nd1a, const int &nd1b,
		   const int &nd2a, const int &nd2b,
		   const int &nd3a, const int &nd3b,
		   const int &n1a,  const int &n1b,
		   const int &n2a,  const int &n2b,
		   const int &n3a,  const int &n3b,
		   const double &x, double &q );

    /* virtual */ void callRHSLoop( const double &tcomp,
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
				double &p);
    
    inline void setInitialData_slow( const double &x, 
				     const double &y, 
				     const double &z,
				     double &c,       
				     double &h,       
				     double &p);

    inline void computeRHS( const double c,  /// IN:       
			    const double h,       
			    const double p,
			    double & rhs_c,  /// OUT:
			    double & rhs_h,
			    double & rhs_p,
			    ///                 OPTIONAL input
			    const double x=0.,  
			    const double y=0., 
			    const double z=0. );


    //..DATA..................................................................


    //..PARAMETERS............................................................
    ///....components
    //enum { cc=0, pc=1, hc=2,
    //bcfast=3, bcslow=4, numberOfComponents};
    
    ///..parameters, as Table 1, in Wagner, Pearson, Keizer
    double nu_L;       
    double d_I;
    double k_p;        
    double nu_P;       
    double I_s;        
    double C_er;       
    
    double tau_0;       /// inhibition relax time
    double d_act;      
    
    double diffCalcium; /// =D, micro m^2/sec     diffusion coeff/calcium
    
    double lambda;     
    double d_inh;      
    double beta;       
    
    /// .... params not used by Wagner et al
    double diffIP3;     /// =0                  diffusion coeff/IP3
    double k_i;         /// =0 or 0.2 1/sec,    IP3 degradation

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
    enum IP3DistributionType {HomogeneousIP3=0, BlobIP3=1, 
			      RadialIP3=2, RadialWithGradientIP3=3};
    enum CaDistributionType {HomogeneousCa=0, BlobCa=1 };
    
    IP3DistributionType ip3Distribution;
    CaDistributionType  caDistribution;
    double ip3Dist_I_s, ip3Dist_I_h, ip3Dist_r_c, ip3Dist_I_w; 
    double ip3Dist_xstar, ip3Dist_Iprime_h, ip3Dist_Iprime_w;
    
#endif

  }; // end class ReactionLiRinzelWagner
  

  ///
  /// ..INLINE functions
  ///
  
  /// -- initial data ----------------------------------------
  
  inline void ReactionLiRinzelWagner::
  setInitialData( const double &x, 
		  const double &y, 
		  const double &z,
		  double &c,       
		  double &h,       
		  double &p)
    {
      GenericCalciumIP3Reaction::setInitialData( x,y,z, c,h,p); //generic Ca2+/IP3 initial data
      //.. modifications from LiRinzel HERE.
      //    --none
    }
  
  inline void ReactionLiRinzelWagner::  /// initial version: if statement in inner loop
  setInitialData_slow( const double &x, 
		       const double &yin, 
		       const double &zin,
		       double &c,       
		       double &h,       
		       double &p)
    {
      DPrintf(DebugReaction,"-------ERROR should not be calling 'setInitialData_slow' in RxnLiRinzel\n");
#if 0
      double z =zin;
      double y =yin;
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

      if ( caDistribution == ReactionLiRinzelWagner::HomogeneousCa ) {
	c = calcium_0;   ///.. initialize to homogeneous c
      }   
      else if ( caDistribution == ReactionLiRinzelWagner::BlobCa ) {
	c = blobMin_c + ( blobMax_c -  blobMin_c)
	* exp( 
	      -( pow( x- xBlob_c,2) + pow( y- yBlob_c,2))/pow( blobWidth_c,2.) 
	      );
      }

      //..IP3 initial data
      if ( ip3Distribution == ReactionLiRinzelWagner::HomogeneousIP3 ) {
	p = ip3_0;
      }
      else if ( ip3Distribution == ReactionLiRinzelWagner::BlobIP3 ) { 
	const double r2 =  pow( x-xBlob_p,2.)+pow( y- yBlob_p,2.)+pow( z-zBlob_p,2.);
	p  = blobMin_p + ( blobMax_p - blobMin_p) * exp(- r2/pow( blobWidth_p,2.) );
      }
      else if ( ( ip3Distribution == ReactionLiRinzelWagner::RadialIP3)
		||  (ip3Distribution == ReactionLiRinzelWagner::RadialWithGradientIP3)) {
	///fig 4, Wagner et al
	const double r =  sqrt( pow( x,2.)+pow( y,2.)+pow( z,2.) );
      	p= ip3Dist_I_s*( 1.+ip3Dist_I_h * exp( ( r/ ip3Dist_r_c    -1.) / ip3Dist_I_w));
      } 
      else throw "error: unknown IP3 distribution type\n";
      
      if ( ip3Distribution == ReactionLiRinzelWagner::RadialWithGradientIP3) { 
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
      //printf("   IP3:  x=%8.4g  y=%8.4g  z=%8.4g,  p= %8.4g\n", x,y,z,p);

#endif
    }
  
  /// -- reaction RHS evaluator ------------------------------
  
  inline void ReactionLiRinzelWagner::
  computeRHS( const double c,  /// IN:       
	      const double h,       
	      const double p,
	      double & rhs_c,  /// OUT:
	      double & rhs_h,
	      double & rhs_p,
	      ///                 OPTIONAL input
	      const double x=0., /// can omit x,y,z if not needed for reaction terms 
	      const double y=0., 
	      const double z=0. )
  {
    double JFlux = (nu_L + pow( p/( p+ d_I),3.)*pow(c/(c+d_act),3.)*pow(h,3.))*(C_er-c);
    double JSerca= -nu_P* pow(c,2.)/( pow(c,2.) + ( k_p *k_p));
    
    rhs_c = beta*lambda*(  JFlux + JSerca );  // Ca2+
    rhs_p  = - k_i *p;                        // IP3
    rhs_h =  ( d_inh - (c+ d_inh)* h)/ tau_0; // gating variable  
    
    ///rhs_h = ( d_inh/(c+ d_inh) - h)/tau_0;  // non-standard, causes pulses
  }

}; //end namespace CellWave

#endif

