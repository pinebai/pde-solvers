/// CellWave::ReactionLiRinzelWagner -- Li-Rinzel reaction with Wagner et al notation
///
///    Li-Rinzel reaction 
///


#include <math.h>
#include <stdio.h>
#include <iostream>

#include "GenericReaction.h"
#include "ParameterReader.h"
#include "ReactionLiRinzelWagner.h"
#include "GenericReactionMacros.h"
#include "CellWave.h"

using namespace CellWave;

ReactionLiRinzelWagner::
ReactionLiRinzelWagner() {
  this->setReactionName( this->rxnType() );
  initializeParameters();
}

ReactionLiRinzelWagner::
~ReactionLiRinzelWagner() 
{ 
 ///default
};

int  ReactionLiRinzelWagner::
getNumberOfSpecies()
{
  return( 3 ); // this class can only have 3 species.
}
 
double ReactionLiRinzelWagner::
getDiffusionCoefficient(const int component)
{
  double diff=0.;
  if      ( cc == component ) {
    diff = beta*diffCalcium;
  }
  else if ( pc == component ) {
    diff = diffIP3;
  }
  else if ( hc == component ) {
    diff =0.;
  }
  else {
    throw "ReactionLiRinzelWagner:: ERROR, unknown component";
  }
  return(  diff );
}


bool ReactionLiRinzelWagner::
isDiffusive( const int component)
{
  bool diffusionFlag=false;
  double diff = getDiffusionCoefficient( component);

  diffusionFlag = fabs( diff ) > getMinimumDiffusionLimit();
#if 0
  printf(" ReactionLiRinzelWagner::isDiffusive q[%3d]: ",component);
  printf(" nu=%6.2g, nu_min=%5.3g ... diffusive? ", diff, getMinimumDiffusionLimit());
  if( diffusionFlag) {
    printf("yes\n");
  }
  else {
    printf("no\n");
  }
#endif

  return(  diffusionFlag );
}

std::string  ReactionLiRinzelWagner::
getTitle()
{
  std::string name="Li-Rinzel-Wagner model";
  return( name );
}

std::string  ReactionLiRinzelWagner::
getLongComponentName( const int component )
{
  std::string name;
  if( component == cc ) {
    name="Calcium";
  } 
  else if( component == pc ) {
    name="IP3";
  }
  else if( component == hc ) {
    name="IP3R Gating";
  }
  else {
    name="unknown";
  }
  return( name );
}


std::string  ReactionLiRinzelWagner::
getShortComponentName( const int component )
{
  std::string name;
  if( component == cc ) {
    name="c";
  } 
  else if( component == pc ) {
    name="p";
  }
  else if( component == hc ) {
    name="h";
  }
  else {
    name="unknown";
  }
  return( name );
}

//
// .. Loop code
//
void ReactionLiRinzelWagner::
callInitialDataLoop( const double &tcomp,
		     const int&nd, const int&ncomp,
		     const int &nd1a, const int &nd1b,
		     const int &nd2a, const int &nd2b,
		     const int &nd3a, const int &nd3b,
		     const int &n1a,  const int &n1b,
		     const int &n2a,  const int &n2b,
		     const int &n3a,  const int &n3b,
		     const double &xfirst,
		     double &qfirst )
{
  const double *xArray = &xfirst;
  double       *qArray = &qfirst;
  // nd     = number of dimensions
  // ncomp  = number of components
  // nd?a, nd?b for ?=1,2,3 gives _array_ dimensions
  // n?a, n?b   for ?=1,2,3 gives _loop bounds_
  // ForGrid,GRIDINDEX:  defined in GenericReactionImpl.h

  //  printf("ReactionLiRinzelWagner--initial data %d\n", 
  //	 int(ip3Distribution));
  
# define RADIUS( x,y,z ) (sqrt( pow(x,2.) +pow(y,2.) +pow(z,2.)))
  double maxLengthScale=0.;

  if( nd ==3 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int cc=0, pc=1,hc=2;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      
      double &c      = qArray[ GRIDINDEX( i,j,k, cc) ];
      double &p      = qArray[ GRIDINDEX( i,j,k, pc) ];
      double &h      = qArray[ GRIDINDEX( i,j,k, hc) ];
      
      //call single point rxn code -- depends on impl    
      setInitialData( x,y,z, c,h,p);
      double r = RADIUS( x,y,z );
      if ( maxLengthScale < r ) maxLengthScale = r;
    }
  }
  else  if (nd == 2) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int cc=0, pc=1,hc=2;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];

      const double z= 0.; // z must be zero 2D initial data
    
      double &c      = qArray[ GRIDINDEX( i,j,k, cc) ];
      double &p      = qArray[ GRIDINDEX( i,j,k, pc) ];
      double &h      = qArray[ GRIDINDEX( i,j,k, hc) ];
      
      //call single point rxn code -- depends on impl    
      setInitialData( x,y,z, c,h,p);
      double r = RADIUS( x,y, 0. );
      if ( maxLengthScale < r ) maxLengthScale = r;
    }
  }
# undef RADIUS
  GenericCalciumIP3Reaction::setMaximumLengthScale( maxLengthScale );

}

//.. looping calls to array interface, for initial data & rhs
void  ReactionLiRinzelWagner::
callRHSLoop( const double &tcomp,
	     const int&nd,    const int &ncomp,
	     const int &nd1a, const int &nd1b,
	     const int &nd2a, const int &nd2b,
	     const int &nd3a, const int &nd3b,
	     const int &n1a,  const int &n1b,
	     const int &n2a,  const int &n2b,
	     const int &n3a,  const int &n3b,
	     const double &xfirst, const double &qfirst, 
	     double &rhsfirst)
{
  const double *xArray   = &xfirst;
  const double *qArray   = &qfirst;
  double       *rhsArray = &rhsfirst;
  // nd     = number of dimensions
  // ncomp  = number of components
  // nd?a, nd?b for ?=1,2,3 gives _array_ dimensions
  // n?a, n?b   for ?=1,2,3 gives _loop bounds_
  // ForGrid,GRIDINDEX:  defined in GenericReactionImpl.h

  DPrintf(DebugReaction,"ip3_0 = %8.3g\n", ip3_0);

  if( nd==3 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int cc=0, pc=1,hc=2;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      
      const double &c      = qArray[ GRIDINDEX( i,j,k, cc) ];
      const double &p      = qArray[ GRIDINDEX( i,j,k, pc) ];
      const double &h      = qArray[ GRIDINDEX( i,j,k, hc) ];
      
      double &rhs_c        = rhsArray[ GRIDINDEX( i,j,k, cc) ];
      double &rhs_p        = rhsArray[ GRIDINDEX( i,j,k, pc) ];
      double &rhs_h        = rhsArray[ GRIDINDEX( i,j,k, hc) ];
      
      //call single point rxn code -- depends on impl
      computeRHS( c,h,p,   rhs_c, rhs_h, rhs_p,   x,y,z);
    }
  }
  else if( nd==2 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int cc=0, pc=1,hc=2;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      const double z=0.;
      
      const double &c      = qArray[ GRIDINDEX( i,j,k, cc) ];
      const double &p      = qArray[ GRIDINDEX( i,j,k, pc) ];
      const double &h      = qArray[ GRIDINDEX( i,j,k, hc) ];
      
      double &rhs_c        = rhsArray[ GRIDINDEX( i,j,k, cc) ];
      double &rhs_p        = rhsArray[ GRIDINDEX( i,j,k, pc) ];
      double &rhs_h        = rhsArray[ GRIDINDEX( i,j,k, hc) ];
      
      //call single point rxn code -- depends on impl
      computeRHS( c,h,p,   rhs_c, rhs_h, rhs_p,   x,y,z);
    }
  }

}




//
// .. internal member functions
//
void ReactionLiRinzelWagner::
initializeParameters()
{

  //.. parameters
  GenericCalciumIP3Reaction::initializeParameters();

  DPrintf(DebugReaction,"RXN LiRinzelWagner -- initializeParameters.\n");  
  /// ..set parameters from Table 1, in Wagner Pearson, Keizer
  nu_L       = 5e-4;  /// nondimensional
  d_I        = 0.025; /// micro M
  k_p        = 0.4;   /// micro M
  nu_P       = 0.1;   /// micro M
  I_s        = 0.12;  /// micro M
  C_er       = 10.;   /// micro M
  tau_0      = 4;     /// sec
  d_act      = 1.2;   /// micro M
  diffCalcium= 300.;  /// micro m^2/sec
  ///diffCalcium =0.; ///TRY THIS
  lambda     = 112.5; /// 1/sec
  d_inh      = 1.5;   /// micro M
  beta       = 0.053; /// nondimensional
  
  /// ..parameters NOT USED BY WAGNER et al
  diffIP3    = 0.;
  k_i        = 0.; ///2.;   /// 1/sec, IP3 degradation, 0=synthetic IP3S

  //diffIP3    = 300.;      /// IP3 diffusion
  //k_i        = 0.2; ///2.;   /// 1/sec, IP3 degradation, 0=synthetic IP3S
  
}; 

bool ReactionLiRinzelWagner::
readParameterFile( CellWave::ParameterReader &par )
{
  DPrintf(DebugReaction,"LiRinzel -- read param file\n");
  std::string typeString="";
  par.get("parameter file type",typeString);
  bool ok= ( typeString== rxnType() );
  if( !ok ) {
    DPrintf(BroadcastPrint," --ERROR: incompatible parameter file, type= %s\n",
	    typeString.c_str());
    return( ok );
  }

  GenericCalciumIP3Reaction::readParameterFile( par );

  par.get("nu_L",    nu_L,      5e-4);
  par.get("d_I",     d_I,       0.025); // micro M
  par.get("k_p",     k_p,       0.4);   // micro M
  par.get("nu_P",    nu_P,      0.1);   // micro M
  par.get("I_s",     I_s,       0.12);  // micro M
  par.get("C_er",    C_er,      10.);   // micro M
  par.get("tau_0",   tau_0,     4);     // sec
  // par.get("eta",     eta,       1.);    // 1/ micro M

  //tau_0      = 25;     // sec  --- DIFFERENT from Wagner
  par.get("d_act",   d_act,     1.2);   // micro M

  par.get("diffCalcium", diffCalcium, 300.);  // micro m^2/sec
  //diffCalcium =0.; //TRY THIS

  par.get("lambda", lambda,      112.5); // 1/sec
  par.get("d_inh",  d_inh,       1.5);   // micro M
  par.get("beta",   beta,        0.053); // nondimensional

  // ..parameters NOT USED BY WAGNER et al
  par.get("diffIP3",diffIP3,     0.); //300.;   // IP3 diffusion
  par.get("k_i",    k_i,         0.); // 1/sec, IP3 degradation, 0=synthetic IP3S

#if 0
  //..initial data
  par.get("calcium_0", calcium_0,  0.1153); // micro M
  par.get("ip3_0",     ip3_0,      0.24);   // micro M             FIXME -- not as Wagner et al
  par.get("h_0",       h_0,        0.93);   // nondimensional

  ///..initial blob data
  //caDistribution = BlobCa;
  //ip3Distribution = HomogeneousIP3;

  caDistribution  = HomogeneousCa;
  ip3Distribution = RadialWithGradientIP3;

  /// .... for Ca2+
  par.get("blobMin_c",    blobMin_c,    calcium_0 );
  par.get("blobMax_c",    blobMax_c,    2.);
  par.get("xBlob_c",      xBlob_c,     -200.);
  par.get("yBlob_c",      yBlob_c,      0.);
  par.get("zBlob_c",      zBlob_c,      0.);
  par.get("blobWidth_c",  blobWidth_c,  60.);
  
  /// .... for  IP3
  par.get("blobMin_p",    blobMin_p,    ip3_0);
  par.get("blobMax_p",    blobMax_p,    20.);
  par.get("xBlob_p",      xBlob_p,     -100.); 
  par.get("yBlob_p",      yBlob_p,      100.);
  par.get("zBlob_p",      zBlob_p,      0.);
  par.get("blobWidth_p",  blobWidth_p,  200);
#endif

  return( ok );
}








