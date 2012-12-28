/// CellWave::ReactionSlepchenko2Buffer -- Li-Rinzel reaction w/ 2 buffers, as in Slepchencko
///
///


#include <math.h>
#include <stdio.h>
#include <iostream>

#include "GenericReaction.h"
#include "ParameterReader.h"
#include "ReactionSlepchenko2Buffer.h"
#include "GenericReactionMacros.h"
#include "CellWave.h"

using namespace CellWave;

ReactionSlepchenko2Buffer::
ReactionSlepchenko2Buffer() {
  this->setReactionName( this->rxnType() );
  initializeParameters();
}

ReactionSlepchenko2Buffer::
~ReactionSlepchenko2Buffer() 
{ 
 ///default
};

int  ReactionSlepchenko2Buffer::
getNumberOfSpecies()
{
  return( 5 ); // this class can only have 5 species.
}
 
double ReactionSlepchenko2Buffer::
getDiffusionCoefficient(const int component)
{
  double diff=0.;
  if      ( cc == component ) {
    diff = diffCalcium;
  }
  else if ( pc == component ) {
    diff = diffIP3;
  }
  else if ( hc == component ) {
    diff =0.;
  }
  else if ( b1c == component ) {
    diff =0.;
  }
  else if ( b2c == component ) {
    diff = diffB2;
  }
  else {
    throw "ReactionSlepchenko2Buffer:: ERROR, unknown component";
  }
  return(  diff );
}


bool ReactionSlepchenko2Buffer::
isDiffusive( const int component)
{
  bool diffusionFlag=false;
  double diff = getDiffusionCoefficient( component);

  diffusionFlag = fabs( diff ) > getMinimumDiffusionLimit();
  return(  diffusionFlag );
}

std::string  ReactionSlepchenko2Buffer::
getTitle()
{
  std::string name="Slepchenko-2-buffer model";
  return( name );
}

std::string  ReactionSlepchenko2Buffer::
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
  else if( component == b1c ) {
    name="Immobile Buffer (b1)";
  }
  else if( component == b2c ) {
    name="Mobile Buffer   (b2)";
  }
  else {
    name="unknown";
  }
  return( name );
}


std::string  ReactionSlepchenko2Buffer::
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
  else if( component == b1c ) {
    name="b1";
  }
  else if( component == b2c ) {
    name="b2";
  }
  else {
    name="unknown";
  }
  return( name );
}

//
// .. Loop code
//
void ReactionSlepchenko2Buffer::
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

  //  printf("ReactionSlepchenko2Buffer--initial data %d\n", 
  //	 int(ip3Distribution));

#define RADIUS( x,y,z ) (sqrt( pow(x,2.) +pow(y,2.) +pow(z,2.)))
  double maxLengthScale = 0.;
  
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
      double &b1     = qArray[ GRIDINDEX( i,j,k, b1c) ];
      double &b2     = qArray[ GRIDINDEX( i,j,k, b2c) ];
      
      //call single point rxn code -- depends on impl    
      setInitialData( x,y,z, c,h,p, b1, b2);
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
      double &b1     = qArray[ GRIDINDEX( i,j,k, b1c) ];
      double &b2     = qArray[ GRIDINDEX( i,j,k, b2c) ];
      
      //call single point rxn code -- depends on impl    
      setInitialData( x,y,z, c,h,p, b1, b2);
      double r = RADIUS( x,y, 0. );
      if ( maxLengthScale < r ) maxLengthScale = r;
    }
  }
# undef RADIUS
  GenericCalciumIP3Reaction::setMaximumLengthScale( maxLengthScale );
}

//.. looping calls to array interface, for initial data & rhs
void  ReactionSlepchenko2Buffer::
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
      const double &b1     = qArray[ GRIDINDEX( i,j,k, b1c) ];
      const double &b2     = qArray[ GRIDINDEX( i,j,k, b2c) ];
      
      double &rhs_c        = rhsArray[ GRIDINDEX( i,j,k, cc) ];
      double &rhs_p        = rhsArray[ GRIDINDEX( i,j,k, pc) ];
      double &rhs_h        = rhsArray[ GRIDINDEX( i,j,k, hc) ];
      double &rhs_b1        = rhsArray[ GRIDINDEX( i,j,k, b1c) ];
      double &rhs_b2        = rhsArray[ GRIDINDEX( i,j,k, b2c) ];
      
      //call single point rxn code -- depends on impl
      computeRHS( c,h,p, b1, b2,   rhs_c, rhs_h, rhs_p, rhs_b1, rhs_b2,   x,y,z);
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
      const double &b1     = qArray[ GRIDINDEX( i,j,k, b1c) ];
      const double &b2     = qArray[ GRIDINDEX( i,j,k, b2c) ];
      
      double &rhs_c        = rhsArray[ GRIDINDEX( i,j,k, cc) ];
      double &rhs_p        = rhsArray[ GRIDINDEX( i,j,k, pc) ];
      double &rhs_h        = rhsArray[ GRIDINDEX( i,j,k, hc) ];
      double &rhs_b1        = rhsArray[ GRIDINDEX( i,j,k, b1c) ];
      double &rhs_b2        = rhsArray[ GRIDINDEX( i,j,k, b2c) ];
      
      //call single point rxn code -- depends on impl
      computeRHS( c,h,p, b1, b2,  rhs_c, rhs_h, rhs_p, rhs_b1, rhs_b2,  x,y,z);
    }
  }
}

//
// .. internal member functions
//
void ReactionSlepchenko2Buffer::
initializeParameters()
{
  GenericCalciumIP3Reaction::initializeParameters();
  DPrintf(DebugReaction,"RXN Slepchenko2Buffer -- initializeParameters.\n");  
    
  //.. data
  setMaximumLengthScale(0.);

  //.. parameters
  b1Mask = 1;  // set these to =0 to remove buffer terms from the eqns
  b2Mask = 1;

  b1_0 = 0.;
  b2_0 = 0.;
}; 

bool ReactionSlepchenko2Buffer::
readParameterFile( CellWave::ParameterReader &par )
{
  DPrintf(DebugReaction,"Slepchenko2Buffer -- read param file\n");
  std::string typeString="";
  par.get("parameter file type",typeString);
  //bool ok= ( typeString== "Slepchenko2Buffer" );
  bool ok= ( typeString== rxnType() );
  if( !ok ) {
    DPrintf(BroadcastPrint," --ERROR: incompatible parameter file, type= %s\n",
	    typeString.c_str());
    return( ok );
  }

  GenericCalciumIP3Reaction::readParameterFile( par );
    
  //..initial data -- in addition to GenericCa2+/IP3 stuff
  par.get("b1_0",      b1_0,       0.1);    // initial immobile buffer level
  par.get("b2_0",      b2_0,       0.1);    // initial mobile buffer level

  ///..parameters, as Table III, p. 210, in Slepchenko (see also text on p.210)
  par.get("J_0",       J_0,      1000);
  par.get("d_act",     d_act,    0.7);
  par.get("d_inh",     d_inh,    0.6);
  par.get("k_on",     k_on,     2.0);
  par.get("V_m",       V_m,      10.);
  par.get("K_p",       K_p,      0.25);
  par.get("leak",      leak,     1.51e-2);   // = L in Table III, Slepchenko et al.
  par.get("b1 total",  b1_tot,   200.);
  par.get("b2 total",  b2_tot,   9.5);  // default value from Fig. 13 in Slepchenko 

  par.get("K1",        K1,       10.);
  par.get("K2",        K2,       0.24);
  par.get("k1_on",     k1_on,    0.05); 
  k1_off = K1* k1_on;                    // only set K_j, and k_j_on:  k_j_off is computed 
  par.get("k2_on",     k2_on,    0.05);
  k2_off = K2* k2_on;

  //..diffusions: these are D_c, D_p, and D_b1, respectively
  par.get("calcium diffusion", diffCalcium, 300.);  // calcium diffusion
  par.get("IP3 diffusion",     diffIP3,     0.);     // IP3 diffusion: also GenericCalciumIP3Reaction.C
  par.get("b2 diffusion",      diffB2,      50.);    // mobile buffer diffusion
  par.get("k_i",               k_i,         0.); // 1/sec, IP3 degradation, 0=synthetic IP3S

  //..disable buffers at will
  par.get("b1Mask",            b1Mask,      1);  // set to =0 to remove buffer 1
  par.get("b2Mask",            b2Mask,      1);  // set to =0 to remove buffer 2


  return( ok );
}








