/// Brief description: CellWave::GenericReaction is the rxn base class 
///
///   CellWave::GenericReaction
///     part of CellWave
///     -- base class for all rxn mechanisms
///     -- to write your own rxn, derive from this, copy/modify the loop code

#include "GenericReaction.h"
#include "GenericReactionMacros.h"
#include "CellWave.h"

#include <stdlib.h>

//
// .. Parameters
//
void         CellWave::GenericReaction::
setReactionName( const std::string &name)
{
  reactionName = name;
}

std::string  CellWave::GenericReaction::
getReactionName()
{
  return( reactionName );
}

double  CellWave::GenericReaction::
getDiffusionCoefficient( const int component)
{
  return( -1.); // shouldn't call me
}

bool  CellWave::GenericReaction::
isDiffusive( const int component)
{
  return( false);
}

double CellWave::GenericReaction::
getFluxBCCoefficient( const int component )
{
  return( 0. );
}

bool CellWave::GenericReaction::
hasFluxBC( const int component )
{
  return( false );
}


std::string CellWave::GenericReaction::
getTitle( )
{
  std::string name="GenericReaction";
  return ( name);
}

std::string CellWave::GenericReaction::
getLongComponentName( const int component)
{
  std::string name="N/A";
  return( name );
}

std::string CellWave::GenericReaction::
getShortComponentName( const int component)
{
  std::string name="N/A";
  return( name );
}


int  CellWave::GenericReaction::
getNumberOfSpecies()
{
  return( 0 );
}

bool CellWave::GenericReaction::
readParameterFile( CellWave::ParameterReader &param)
{
  // do nothing
  DPrintf( BroadcastPrint, "ERROR -- GenericReaction.readParameterFile called\n");
  throw "error";
  return( false ); //error=false
}

void CellWave::GenericReaction::
printParameters( int iPrintChannel )
{
  return; // base class does nothing
}

void CellWave::GenericReaction::
printOneParameter(int ip, std::string name, double param, double def)
{

  char buf[80], bufdef[80];
  sprintf(buf, "%16.8e", param );
  sprintf(bufdef, "%16.8e", def   );
  std::string outpar(buf);
  std::string outdef(bufdef);

  this->printOneParameter(ip,name,outpar,outdef);
}

void CellWave::GenericReaction::
printOneParameter(int ip, std::string name, int param, int def)
{
  char buf[80], bufdef[80];
  sprintf(buf, "%16d", param );
  sprintf(bufdef, "%16d", def   );
  std::string outpar(buf);
  std::string outdef(bufdef);

  this->printOneParameter(ip,name,outpar,outdef);
}

void CellWave::GenericReaction::
printOneParameter(int ip, std::string name, double param, int def)
{
  this->printOneParameter(ip,name,param,1.0*def);
}

void CellWave::GenericReaction::
printOneParameter(int ip, std::string name,
		  std::string param, std::string def)
{
  DPrintf( ip, "  %20s <%s> = %s \n", 
	   name.c_str(), def.c_str(),  param.c_str()); 
}

//
//.. LOOP/Array interface
//

//.. looping calls to array interface, for initial data & rhs
void CellWave::GenericReaction::
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

  ForGrid(i,j,k) {
    const int xaxis=0, yaxis=1, zaxis=2;
    const int cc=0, pc=1,hc=2;
    
    const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
    const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
    const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
    
    double &c      = qArray[ GRIDINDEX( i,j,k, cc) ];
    double &p      = qArray[ GRIDINDEX( i,j,k, pc) ];
    double &h      = qArray[ GRIDINDEX( i,j,k, hc) ];

    printf("ERROR: you should not call GenericReaction::callInitialDataLoop\n");
    exit(-1);
    //call single point rxn code -- depends on impl    
    //setInitialData( x,y,z, c,h,p);
  }
}

//.. looping calls to array interface, for initial data & rhs
void CellWave::GenericReaction::
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
    
    printf("ERROR: you should not call GenericReaction::callInitialDataLoop\n");
    exit(-1);
    //call single point rxn code -- depends on impl
    //computeRHS( c,h,p,   rhs_c, rhs_h, rhs_p,   x,y,z);
  }
}

void CellWave::GenericReaction::
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
  printf("--------should not be calling GenericReaction ----\n");
}



