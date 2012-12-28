//
// NucleusGridFunction 
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>

#include "CellWave.h"
#include "NucleusGridFunction.h"

using namespace CellWave;

NucleusGridFunction::
NucleusGridFunction()
{
  pCG        = NULL;
  pNucleusGF = NULL;
  nucleusBoundaryThickness=20.;

}

NucleusGridFunction::
~NucleusGridFunction()
{
  pCG = NULL;
}

void NucleusGridFunction::
updateToMatchGrid( CompositeGrid &cg, 
		   doubleCompositeGridFunction &nucleusGF )
{
  nucleusGF.updateToMatchGrid( cg );
  
  pCG        = &cg;
  pNucleusGF = &nucleusGF;
}

void NucleusGridFunction::
readParameterFile( ParameterReader &params )
{
  //read:
  // nucleus type=none, file, box, sphere
  // nucleus file= <file.cwn>
  // nucleus corners= x0 y0 z0 x1 y1 z1
  // nucleus center= x0 y0 z0
  // nucleus radius= r0
  // nucleus boundary thickness= t0

  const std::string keyBoundaryThickness="nucleus boundary thickness";
  double thickness;
  params.get( "nucleus boundary thickness", thickness, -1. );
  if ( thickness>=0. ) nucleusBoundaryThickness=thickness;

  typedef boost::tokenizer<>::iterator TokIterator;
  typedef std::vector<double>          DoubleVector;

  bool noNucleus=false;
  std::string nucleusType;
  params.get( "nucleus type", nucleusType, "");
  if ( (nucleusType == "box" )   || (nucleusType == "cubic")) {
    Nucleus::NucleusShape nucleusShape=Nucleus::BoxNucleus;
    std::string doubleList="";
    params.get("nucleus corners",  doubleList ); 
    DPrintf(DebugSolver,"..NucleusGridFunction: nucleus type=box.\n");
    DPrintf(DebugSolver,"....  nucleus corners = '%s'\n", doubleList.c_str());

    boost::tokenizer<> tok(doubleList);
    DoubleVector       corners;

    DPrintf(DebugSolver,"    corners are = ");
    corners.clear();
    for (TokIterator it=tok.begin(); it!=tok.end(); ++it) {
      double q;
      sscanf(it->c_str(), "%lf", &q);
      DPrintf(DebugSolver," %lf ", q);
      corners.push_back( q );
    }
    DPrintf(DebugSolver,"\n");
    //..create the nucleic info
    if( corners.size() >5 ) { 
      Nucleus thisNuc;
      thisNuc.setID(1);
      thisNuc.setBoundaryThickness( nucleusBoundaryThickness );
      thisNuc.setCorners(corners[0], corners[1], corners[2],
			 corners[3], corners[4], corners[5]);
      thisNuc.setShape( nucleusShape );

      nucleus.push_back( thisNuc );
    }
  }
  else if ( (nucleusType == "sphere") || ( nucleusType == "spherical")) {
    Nucleus::NucleusShape nucleusShape=Nucleus::SphericalNucleus;
    double radius;
    params.get("nucleus radius",radius, 20.);
    std::string doubleList="";
    params.get("nucleus center",  doubleList ); 
    DPrintf(DebugSolver,"..NucleusGridFunction: nucleus type=sphere not supported yet:\n");
    DPrintf(DebugSolver,"..  nucleus center = '%s'\n", doubleList.c_str());

    boost::tokenizer<> tok(doubleList);
    DoubleVector       center;

    DPrintf(DebugSolver,"    center is = ");
    center.clear();
    for (TokIterator it=tok.begin(); it!=tok.end(); ++it) {
      double q;
      sscanf(it->c_str(), "%lf", &q);
      DPrintf(DebugSolver," %lf ", q);
      center.push_back( q );
    }
    DPrintf(DebugSolver,"\n");
    //..create the nucleic info
    double x0=0., y0=0., z0=0.;
    if( center.size() >1) {
      x0=center[0];
      y0=center[1];
    }
    if( center.size()>2) {
      z0=center[2];
    }

    Nucleus thisNuc;
    thisNuc.setID(1);
    
    thisNuc.setBoundaryThickness( nucleusBoundaryThickness );
    thisNuc.setCenter( x0, y0, z0);
    thisNuc.setRadius( radius);
    thisNuc.setShape( nucleusShape );
    nucleus.push_back( thisNuc );

  }
  else if ( nucleusType == "file" ) {
    std::string nucleusFileName="";
    params.get("nucleus file",nucleusFileName);
    DPrintf(DebugSolver,"nuclei read in from file='%s'\n", nucleusFileName.c_str());
    readCellNucleusFile( nucleusFileName );
  }
  else { // all other options mean--> no nucleus
    noNucleus=true;
  }
}

bool NucleusGridFunction::
readCellNucleusFile( const std::string cn_filename,
			    const std::string gridFileName/* = ""*/ )
{
  bool okFlag=false;
  clearGrid2NucleusMap();
  clearNucleus2GridMap();

  FILE *fp=fopen( cn_filename.c_str(),"r");
  if( !fp ) {
    okFlag=false;
    return( okFlag );
  }
  DPrintf(DebugSolver,"reading file '%s'...\n",cn_filename.c_str());
  DPrintf(DebugSolver,"------------------------------\n");
  
  const int bufferLength=1000;
  char buffer[bufferLength];
  
  while( fgets( buffer, bufferLength, fp)) {
    const int lineLength=strlen(buffer);

    typedef std::vector<std::string> StringVector;
    StringVector tokens;

    if( lineLength>0 ) { 
      if( buffer[0]=='#' ) {
	buffer[lineLength-1]=0;
	DPrintf(DebugSolver,"comment< %s >\n", buffer);
      }
      else {
	using namespace std;
	Nucleus thisNuc;
	buffer[lineLength-1]=0;
	//printf("regular< %s >\n", buffer);
	string line(buffer);
	//cout << "<"<<line<<">\n";

	typedef boost::tokenizer<>::iterator           TokIterator;
	typedef boost::char_delimiters_separator<char> TokSeparator;
	//see boost::tokenizer 'char_delimiters_separator' docs
	//  sep( returnable=false, returned="", separators=0)
	// --> separators=0 means anything for which iswhitespace() 
	//     is true is a separator
	TokSeparator sep(false,"",0);
	boost::tokenizer< TokSeparator> tok(line, sep);
	for (TokIterator it=tok.begin(); it!=tok.end(); ++it) {
	  //DPrintf(DebugSolver,"<%s> ", it->c_str());
	  tokens.push_back( *it );
	}
	//..INPUT FILE FORMAT FOR .cwn
	// format: <nucleus #> <radius> <x y z of center> <grid ID(s)>
	// lines with '#' in column 1 are comments
	int nID;      const int idIndex=0; 
	double rad;   const int idRadius=1;
	double x,y,z; const int idX=2, idY=3, idZ=4;
	//std::vector<int> gridIDs;

	sscanf(tokens[idIndex].c_str(),  "%d",   &nID);
	sscanf(tokens[idRadius].c_str(), "%lg",  &rad);
	sscanf(tokens[idX].c_str(),      "%lg",  &x);
	sscanf(tokens[idY].c_str(),      "%lg",  &y);
	sscanf(tokens[idZ].c_str(),      "%lg",  &z);
	DPrintf(DebugSolver,"  #tokens=%d, ztoken=%s -- ", 
	       tokens.size(), tokens[idZ].c_str());

	DPrintf(DebugSolver,"id=%d, R=%f, x=%f, y=%f, z=%f,",nID,rad,x,y,z);
	DPrintf(DebugSolver,"\n");
	thisNuc.setID( nID);
	thisNuc.setCenter(x,y,z);
	thisNuc.setRadius(rad);
	
	Nucleus::NucleusShape nucleusShape=Nucleus::SphericalNucleus;
	thisNuc.setShape( nucleusShape );
	thisNuc.setBoundaryThickness( nucleusBoundaryThickness );

	nucleus.push_back( thisNuc );

	DPrintf(DebugSolver,"gridIDs for nucleus # %d=",nID);
	for( int i=idZ+1; i<tokens.size(); ++i ) {
	  int gID=-1;
	  sscanf(tokens[i].c_str(), "%d", &gID);
	  DPrintf(DebugSolver," %d ",gID);
	  grid2NucleusMap.insert( std::make_pair(gID,nID));
	  nucleus2GridMap.insert( std::make_pair(nID,gID));
	}		 
	DPrintf(DebugSolver,"\n");
      }
    };

  }
  DPrintf(DebugSolver,"-------------done-------------\n");
  fclose(fp);
  okFlag=true;

  return( okFlag ); 
  
  
}

void NucleusGridFunction::
evaluateGridFunction( double tcomp  /* = 0. */ ) //CONTINUE FROM HERE
{
  assert( pCG != NULL );
  int grid;
  Index If1,If2,If3; // full range
  const int ncomp=1; // number of components in the mask

  DPrintf(DebugSolver,">evaluateGridFunction:\n");
  for( grid=0; grid< pCG->numberOfComponentGrids(); grid++ )  {
    DPrintf(DebugSolver,"..grid %d\n", grid);
    MappedGrid & mg =   (*pCG)[grid];
    getIndex(mg.dimension(),     If1,If2,If3);    // INDICES = by dimension   
    realArray & mask  =   (*pNucleusGF)[grid];
    realArray & x     =   mg.vertex();  // array of vertices

    realArray temp(If1,If2,If3, ncomp);

    const IntegerArray & d = mg.dimension();
    const IntegerArray & gir= mg.gridIndexRange();
    const int nd= pCG->numberOfDimensions();

    //..Loop over the nuclei for this grid
    mask = 1.; // mask is ==1 for points in cytosol, ==0 inside the nucleus

    IDVector nuclei;
    //getGridNuclei( grid, nuclei); //optimized
    getAllNuclei( nuclei );
    for ( int j=0; j<nuclei.size(); ++j ) {
      Nucleus &thisNuc = getNucleus( nuclei[ j ] );
      std::string nucleusType;
      
      DPrintf(DebugSolver,".... on grid %d, evaluate nucleus %d (%s)\n", 
	     grid, thisNuc.getID(), thisNuc.getNucleusShapeName().c_str());

      thisNuc.getMaskArray(tcomp,nd, ncomp, d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2),
			   gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
			   *x.getDataPointer(),*temp.getDataPointer() );
      mask = mask*temp; // collect masks from the nuclei to one mask for this grid
    } // end for j (over the nuclei for this grid) 
  } // end for grid
}

void NucleusGridFunction::
maskAGridFunction( doubleCompositeGridFunction &u)
{
  //-->ASSUME u is compatible with *pNucleusGridFunction.<-- add a check for this?
  assert( pNucleusGF != NULL );

  u = u* ( *pNucleusGF ); // apply the mask on u
  
}

void NucleusGridFunction::
printInfo()
{
  DPrintf(DebugSolver,"--List of Nuclei:\n");
  for(int i=0; i<nucleus.size(); ++i ) {
    Nucleus &nucl= nucleus[i];
    int id;
    double x,y,z,rad;
    
    id= nucl.getID();
    rad=nucl.getRadius();
    nucl.getCenter( x,y,z);
    
    DPrintf(DebugSolver," id=%d, x=%f, y=%f, z=%f, r=%f", id,x,y,z,rad);

    IDVector nucleusGrids;
    getNucleusGrids( id, nucleusGrids);

    DPrintf(DebugSolver,", grids= ");
    for( int i=0; i<nucleusGrids.size(); ++i ) {
      DPrintf(DebugSolver," %d ",nucleusGrids[i]);
    }
    DPrintf(DebugSolver,"\n");

  }
}

bool NucleusGridFunction::
checkConsistency()
{
  assert( pCG != NULL );
  assert( pNucleusGF != NULL );

  const int nGrids = pCG->numberOfGrids();

  //..iterate through the nuclei in grid2NucleusMap & nucleus2GridMap 
  //  to check the numbers stay within bounds.
 
}

int NucleusGridFunction::
getGridNuclei( int gridID, IDVector &nuclei)
{
  int numberOfNuclei=0;
  nuclei.clear();
  IterateIntMap pos;
  //DPrintf(DebugSolver,"..getGridNuclei( gridID = %d ), nuclei= ",gridID );
  for (pos =grid2NucleusMap.lower_bound( gridID );
       pos!=grid2NucleusMap.upper_bound( gridID );
       ++pos ) {

    //DPrintf(DebugSolver," %d ", pos->second);
    nuclei.push_back( pos->second );
    numberOfNuclei++;
  }
  //DPrintf(DebugSolver,"\n");

  return( numberOfNuclei );
}

int NucleusGridFunction::
getNucleusGrids( int nucleusID, IDVector &gridIDs)
{
  int numberOfNucleusGrids=0;
  gridIDs.clear();
  IterateIntMap pos;
  //  DPrintf(DebugSolver,"..getNucleusGrids( nID = %d ), grids= ",nucleusID);
  for (pos =nucleus2GridMap.lower_bound( nucleusID );
       pos!=nucleus2GridMap.upper_bound( nucleusID );
       ++pos ) {
    
    //DPrintf(DebugSolver," %d ", pos->second );
    gridIDs.push_back( pos->second );
    numberOfNucleusGrids++;
  }
  //DPrintf(DebugSolver,"\n");
  
  return( numberOfNucleusGrids );
}

int NucleusGridFunction::
getAllNuclei( IDVector &nuclei )
{
  nuclei.clear();
  DPrintf(DebugSolver,"..getAllNuclei -- # of nuclei= %d\n", getNumberOfNuclei() );
  fflush(0);
  for(int i=1; i<= getNumberOfNuclei(); ++i ) {
    nuclei.push_back( i );
  }
  return ( getNumberOfNuclei());
};
