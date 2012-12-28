/// CellWave::Reaction2Buffer -- Li-Rinzel reaction w/ 2 buffers, as in Slepchencko
///
///


#include <math.h>
#include <stdio.h>
#include <iostream>

#include "GenericReaction.h"
#include "ParameterReader.h"
#include "Reaction2Buffer.h"
#include "GenericReactionMacros.h"
#include "CellWave.h"

using namespace CellWave;

Reaction2Buffer::
Reaction2Buffer() {
  this->setReactionName( this->rxnType() );
  initializeParameters();
}

Reaction2Buffer::
~Reaction2Buffer() 
{ 
 ///default
};

int  Reaction2Buffer::
getNumberOfSpecies()
{
  return( 5 ); // this class can only have 5 species.
}
 
double Reaction2Buffer::
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
    throw "Reaction2Buffer:: ERROR, unknown component";
  }
  return(  diff );
}


bool Reaction2Buffer::
isDiffusive( const int component)
{
  bool diffusionFlag=false;
  double diff = getDiffusionCoefficient( component);

  diffusionFlag = fabs( diff ) > getMinimumDiffusionLimit();
  return(  diffusionFlag );
}

double Reaction2Buffer::
getFluxBCCoefficient( const int component )
{
  double flux=0.;
  if      ( cc == component ) {
    flux = fluxBCCoefficientCalcium;
  }
  else if ( pc == component ) {
    flux = fluxBCCoefficientIP3;
  }
  else if ( hc == component ) {
    flux =0.;
  }
  else if ( b1c == component ) {
    flux = fluxBCCoefficientB1;
  }
  else if ( b2c == component ) {
    flux = fluxBCCoefficientB2;
  }
  else {
    throw "Reaction2Buffer:: ERROR, unknown component";
  }
  
  if( isDiffusive( component ) ) {
    const double diff=getDiffusionCoefficient( component );
    flux = flux/diff;
  } else {
    flux = 0.;
  }

  return( flux );
}

bool Reaction2Buffer::
hasFluxBC( const int component )
{
  bool fluxFlag=false;
  const double flux= getFluxBCCoefficient( component );
  const double diff= getDiffusionCoefficient( component );
  const double minDiff = 

  fluxFlag = fabs( flux ) > getMinimumDiffusionLimit();
  fluxFlag = fluxFlag && isDiffusive( component);
  // only species that diffuse can have fluxBCCoeff= flux/diffusion 
  return( fluxFlag );
}


std::string  Reaction2Buffer::
getTitle()
{
  std::string name="2 Buffer model: Li-Rinzel + 2x explicit buff + IP3";
  return( name );
}

std::string  Reaction2Buffer::
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


std::string  Reaction2Buffer::
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
void Reaction2Buffer::
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
  // ForGrid,GRIDINDEX:  defined in GenericReactionMacros.h

  //  printf("Reaction2Buffer--initial data %d\n", 
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
void  Reaction2Buffer::
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
  if (useRapidBuffering) {
    callRHSLoop_explicitBuffering( tcomp,nd, ncomp,nd1a,nd1b,
			  nd2a,nd2b, nd3a, nd3b, n1a, n1b, n2a,n2b,
			  n3a, n3b, xfirst, qfirst, rhsfirst);
  } else {
    callRHSLoop_rba( tcomp,nd, ncomp,nd1a,nd1b,
		     nd2a,nd2b, nd3a, nd3b, n1a, n1b, n2a,n2b,
		     n3a, n3b, xfirst, qfirst, rhsfirst);    
  }
}

//.. looping calls to array interface, for initial data & rhs/ full buffering
void  Reaction2Buffer::
callRHSLoop_explicitBuffering( const double &tcomp,
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
  // ForGrid,GRIDINDEX:  defined in GenericReactionMacros.h

  //DPrintf(DebugReaction,"ip3_0 = %8.3g\n", ip3_0);

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


//.. looping calls to array interface, for initial data & rhs/ with RBA
void  Reaction2Buffer::
callRHSLoop_rba( const double &tcomp,
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
  // ForGrid,GRIDINDEX:  defined in GenericReactionMacros.h

  //DPrintf(DebugReaction,"ip3_0 = %8.3g\n", ip3_0);

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
// ..for computing the rapid buffer species
//
//.. looping calls to array interface, for initial data & rhs/ full buffering
void  Reaction2Buffer::
computeRapidBuffers( const double &tcomp,
	     const int&nd,    const int &ncomp,
	     const int &nd1a, const int &nd1b,
	     const int &nd2a, const int &nd2b,
	     const int &nd3a, const int &nd3b,
	     const int &n1a,  const int &n1b,
	     const int &n2a,  const int &n2b,
	     const int &n3a,  const int &n3b,
	     const double &xfirst, const double &qfirst)
{
  const double *xArray   = &xfirst;
  const double *qArray   = &qfirst;
  // nd     = number of dimensions
  // ncomp  = number of components
  // nd?a, nd?b for ?=1,2,3 gives _array_ dimensions
  // n?a, n?b   for ?=1,2,3 gives _loop bounds_
  // ForGrid,GRIDINDEX:  defined in GenericReactionMacros.h

  //DPrintf(DebugReaction,"ip3_0 = %8.3g\n", ip3_0);

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
      
      //call single point rxn code -- depends on impl
      //computeRHS( c,h,p, b1, b2,   rhs_c, rhs_h, rhs_p, rhs_b1, rhs_b2,   x,y,z);
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
      
      //call single point rxn code -- depends on impl
      //computeRHS( c,h,p, b1, b2,  rhs_c, rhs_h, rhs_p, rhs_b1, rhs_b2,  x,y,z);
    }
  }
}

//
// .. internal member functions
//
void Reaction2Buffer::
initializeParameters()
{
  GenericCalciumIP3Reaction::initializeParameters();
  //DPrintf(DebugReaction,"RXN 2Buffer -- initializeParameters.\n");  
    
  //.. data
  setMaximumLengthScale(0.);

  //.. parameters
  b1Mask = 1;  // set these to =0 to remove buffer terms from the eqns
  b2Mask = 1;

  b1_0 = 0.;
  b2_0 = 0.;

  useRapidBuffering = false;
}; 

bool Reaction2Buffer::
readParameterFile( CellWave::ParameterReader &par )
{
  DPrintf(DebugReaction,"2Buffer -- read param file\n");
  std::string typeString="";
  par.get("parameter file type",typeString);
  //bool ok= ( typeString== "2Buffer" );
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

  //..for some parameters, not yet clear what their standard values are
  par.get("nu_L",     nu_L,      1.51e-2);  //FIX
  par.get("nu_c",     nu_c,         0);     //FIX
  par.get("d_act",     d_act,    0.7);
  par.get("d_inh",     d_inh,    0.6);
  par.get("d_I",       d_I,         0);     //FIX
  par.get("nu_m",     nu_m,      10.);
  par.get("K_p",       K_p,      0.25);
  par.get("C_er",     C_er,        10);
  par.get("k_on",     k_on,     2.0);

  par.get("gamma",   gamma,        0.);
  par.get("b1 total",  B1_tot,   200.);
  par.get("b2 total",  B2_tot,   9.5);  // default value from Fig. 13 in Slepchenko 

  par.get("K1",        K1,       10.);
  par.get("K2",        K2,       0.24);
  par.get("k1_on",     k1_on,    0.05); 
  k1_off = K1* k1_on;                    // only set K_j, and k_j_on:  k_j_off is computed 
  par.get("k2_on",     k2_on,    0.05);
  k2_off = K2* k2_on;

  par.get("H_flux",    H_flux,   2.);

  //..diffusions: these are D_c, D_p, and D_b1, respectively
  par.get("calcium diffusion", diffCalcium, 300.);  // calcium diffusion
  par.get("IP3 diffusion",     diffIP3,     0.);     // IP3 diffusion: also GenericCalciumIP3Reaction.C
  par.get("b2 diffusion",      diffB2,      50.);    // mobile buffer diffusion

  //..fluxBC Coefficient A: du/dn = A [ u ]
  par.get("flux bc coefficient calcium", fluxBCCoefficientCalcium,0.);
  par.get("flux bc coefficient IP3",     fluxBCCoefficientIP3,    0.);
  par.get("flux bc coefficient b1",      fluxBCCoefficientB1,     0.);
  par.get("flux bc coefficient b2",      fluxBCCoefficientB2,     0.);


  //..disable buffers at will
  par.get("b1Mask",            b1Mask,      1);  // set to =0 to remove buffer 1
  par.get("b2Mask",            b2Mask,      1);  // set to =0 to remove buffer 2


  return( ok );
}



void Reaction2Buffer::
printParameters(int iPrintChannel)
{
  const int ip=iPrintChannel;

  DPrintf(ip, "Parameters for Reaction2Buffer:\n");
  DPrintf(ip, "  'parameter' <'default'>= 'value'\n");
  DPrintf(ip, "====================================\n");

  DPrintf(ip, "GenericCa2+/IP3 parameters:\n\n");
  GenericCalciumIP3Reaction::printParameters(ip);
  
  DPrintf(ip, "\nSpecific parameters for 2Buffer:\n\n");
  printOneParameter(ip,"b1_0",      b1_0,       0.1);  
  printOneParameter(ip,"b2_0",      b2_0,       0.1);    // initial mobile buffer level

  printOneParameter(ip,"nu_L",     nu_L,      1.51e-2);
  printOneParameter(ip,"nu_c",     nu_c,         0);   
  printOneParameter(ip,"d_act",     d_act,    0.7);
  printOneParameter(ip,"d_inh",     d_inh,    0.6);
  printOneParameter(ip,"d_I",       d_I,         0);   
  printOneParameter(ip,"nu_m",     nu_m,      10.);
  printOneParameter(ip,"K_p",       K_p,      0.25);
  printOneParameter(ip,"C_er",     C_er,        10);
  printOneParameter(ip,"k_on",     k_on,     2.0);

  printOneParameter(ip,"gamma",   gamma,        0.);
  printOneParameter(ip,"b1 total",  B1_tot,   200.);
  printOneParameter(ip,"b2 total",  B2_tot,   9.5); 

  printOneParameter(ip,"K1",        K1,       10.);
  printOneParameter(ip,"K2",        K2,       0.24);
  printOneParameter(ip,"k1_on",     k1_on,    0.05); 
  printOneParameter(ip,"k2_on",     k2_on,    0.05);
  printOneParameter(ip,"H_flux",    H_flux,   2.);

  //..diffusions: these are D_c, D_p, and D_b1, respectively
  printOneParameter(ip,"calcium diffusion", diffCalcium, 300.);
  printOneParameter(ip,"IP3 diffusion",     diffIP3,     0.); 
  printOneParameter(ip,"b2 diffusion",      diffB2,      50.);

  //..fluxBC Coefficient A: du/dn = A [ u ]
  printOneParameter(ip,"flux bc coefficient calcium", fluxBCCoefficientCalcium,0.);
  printOneParameter(ip,"flux bc coefficient IP3",     fluxBCCoefficientIP3,  0.);
  printOneParameter(ip,"flux bc coefficient b1",      fluxBCCoefficientB1,   0.);
  printOneParameter(ip,"flux bc coefficient b2",      fluxBCCoefficientB2,  0.);

  //..disable buffers at will
  printOneParameter(ip,"b1Mask",            b1Mask,      1); 
  printOneParameter(ip,"b2Mask",            b2Mask,      1); 

  //..for comparison with Slepchenko et al.
  DPrintf(ip, "\nComparison with Slepchenko et al: (derived qty's)\n");
  double leak_ = nu_L*(C_er - calcium_0);
  double J0_   = nu_c*(C_er - calcium_0)*pow( ip3_0/(ip3_0+d_I), 3.);

  printOneParameter(ip,"LEAK",              leak_,     0.0151);
  printOneParameter(ip,"J0",                J0_,       1000);
  
  //..temp-- debug info on screen:
  DPrintf(CellWave::PrintOut, "\nComparison with Slepchenko et al: (derived qty's)\n");
  printOneParameter(PrintOut,"LEAK",              leak_,     0.0151);
  printOneParameter(PrintOut,"J0",                J0_,       1000);
  DPrintf(CellWave::PrintOut, "\n");


} 
