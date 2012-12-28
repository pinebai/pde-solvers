//
// CellWave::Probes -- probe dynamic chemical concentrations during simulations
//  $Id: Probes.C,v 1.5 2003/05/21 17:56:15 pfast Exp $
//
#include <assert.h>
#include <time.h>
#include <stdio.h>

#include <boost/tokenizer.hpp>

#include "Overture.h"
#include "interpolatePoints.h"
#include "CellWave.h"
#include "Probes.h"

using namespace CellWave;
    
Probes::Probes( int numberOfSpecies_ /* = 0 */)
{
  useProbes = false;
  probeHowOften = 1; //1= every step  
  numberOfSpecies = numberOfSpecies_;
}

Probes::~Probes()
{
}

bool  Probes::
readParameterFile( CellWave::ParameterReader &params )
{
  DPrintf(DebugSolver,"..entering Probes::readParameterFile\n");
  params.get( "probe location filename", probeLocationFileName, "" );
  if ( probeLocationFileName == "") {
    useProbes=false;
  }
  else {
    useProbes=true;
  }
  params.get( "probe output filename", probeOutputFileName, "" );
  if ( probeOutputFileName == "" ) {
    useProbes=false;
  }
  //else --> useProbes=true if it was set to that before
  this->readProbeLocations();

  if( useProbes ) {
    parameterFileName = params.getParameterFileName();
    params.get("probe frequency", probeHowOften, 1 );
  }
  DPrintf( DebugSolver, "Probes initialized:\n");
  DPrintf( DebugSolver, "  use probes = %d\n", int( useProbes ));
  DPrintf( DebugSolver, "  probe locations from '%s'\n", 
	   probeLocationFileName.c_str());
  DPrintf( DebugSolver, "  probe output file is '%s'\n",
	   probeOutputFileName.c_str());
  DPrintf( DebugSolver, "  probe output frequency is %d\n", 
	   probeHowOften);

}

void Probes::
initializeProbeValues( int numberOfSpecies_ ) 
{
  int nProbes = numberOfProbes();
  numberOfSpecies = numberOfSpecies_;
  probeCurrentValue.redim( nProbes, numberOfSpecies);
}

void Probes::
readProbeLocations()
{
  DPrintf(DebugSolver,"..entering Probes::readProbeLocations\n");
  if( useProbes ) {
    FILE *fp;
    fp = fopen( probeLocationFileName.c_str(), "r");
    if ( fp == NULL ) {
      DPrintf( BroadcastPrint, "CellWave::Probes ERROR!\n");
      DPrintf( BroadcastPrint, " Cannot open probe location file '%s',",
	       probeLocationFileName.c_str());
      DPrintf( BroadcastPrint, " probes are turned off.");
      useProbes= false;
    }
    else { // probe input file ok, read the probe xyz
      this->clearProbeLocations();
      std::vector< OneProbe > probeVector;

      DPrintf( DebugSolver, "Probes: reading locations from '%s'\n",
	       probeLocationFileName.c_str());

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
	    //..INPUT FILE FORMAT FOR .probexyz
	    // format: <probe number> <x y z of center>
	    // lines with '#' in column 1 are comments
	    int pID;      const int idIndex=0;
	    double x,y,z; const int idX=1, idY=2, idZ=3;
	    //std::vector<int> gridIDs;
	    
	    sscanf(tokens[idIndex].c_str(),  "%d",   &pID);
	    sscanf(tokens[idX].c_str(),      "%lg",  &x);
	    sscanf(tokens[idY].c_str(),      "%lg",  &y);
	    if( tokens.size() >= 4 ) {
	      sscanf(tokens[idZ].c_str(),      "%lg",  &z);
	    }
	    else z=0.; 

	    DPrintf(DebugSolver,"id=%d, x=%f, y=%f, z=%f",pID,x,y,z);
	    DPrintf(DebugSolver,"\n");
	    
	    OneProbe probe( pID, x,y,z );
	    probeVector.push_back( probe );
	  }
	}
      } // while fgets
      fclose(fp);
      int nProbes= probeVector.size();
      probeLocation.redim( nProbes, 3 );
      for( int i=0; i< nProbes; ++i ) {
	assert( i < probeVector.size() ); 
	assert( probeVector[i].id == i+1 );
	assert( i < probeLocation.getLength( axis1 ));

	OneProbe & p = probeVector[i];

	probeLocation(i, axis1) = p.x;
	probeLocation(i, axis2) = p.y;
	probeLocation(i, axis3) = p.z;
      }
      initializeProbeValues( numberOfSpecies );
      
    } //endif fp==NULL
  } //endif useProbes
}

void Probes::
clearProbeLocations()
{
  probeLocation.redim(0);
}

void Probes::
openOutputFile(std::string comment ) // //write the probe locations, in comment lines
{
  DPrintf(DebugSolver,"..entering Probes::openOutputFile\n");
  if( !isOK() ) {  
    DPrintf(DebugSolver,"..quick exit from Probes::openOutputFile, not initialized\n");
    return;
  }

  //open file, then write a header with the probe locations
  //##CellWave probe file
  //##"comment"
  //#simulationStarted=
  //#paramFile=
  //#probeLocationFile=
  //#outputShowFile=
  //#outputComponents= c h p
  //#numberOfProbes=
  //#probeList=follows
  //## probe number  x  y  z
  //#   1            1. 2. 3.
  //##end of probe list
  //#probeOutput=follows
  //## probe number time <components>
  fProbe = fopen( probeOutputFileName.c_str(), "w");
  if( fProbe == NULL ) {
    DPrintf(BroadcastPrint,"Probes ERROR -- unable to open output file '%s'\n",
	    probeOutputFileName.c_str());
    throw "error";
  }
  const bool useHeader= !noHeader; 
  if ( useHeader && (fProbe != NULL) ) {

    fprintf( fProbe, "##CellWave  probe output file\n");
    fprintf( fProbe, "##%s\n", comment.c_str());
    
    time_t date;
    time(&date ); //from <time.h>
    fprintf( fProbe, "#simulationStarted=  %s", asctime(localtime( &date) ));
    fprintf( fProbe, "#parameterFile=      %s\n", parameterFileName.c_str());
    fprintf( fProbe, "#probeLocationFile=  %s\n", probeLocationFileName.c_str());
    fprintf( fProbe, "#outputShowFile=     -N/A-\n" );
    fprintf( fProbe, "#outputComponents=   %d\n", numberOfSpecies);
    fprintf( fProbe, "#numberOfProbes=     %d\n", numberOfProbes());
    fprintf( fProbe, "#probeList=notAvailable\n");
    //fprintf( fProbe, "#probeList=follows\n");
    //call output probe list:fprintf( fProbe, "## probe number  x  y  z\n");
    //call output probe list: fprintf( fProbe, "#   1            1. 2. 3.
    //call output probe listfprintf( fProbe, "##end of probe list
    fprintf( fProbe, "#probeOutput=follows\n");
    fprintf( fProbe, "## time probe1.component1 probe1.component2 ...\n");
  } else if (!useHeader) {
    DPrintf( DebugSolver, "CellWave::Probes  Warning -- no header output to probe sequence\n");
  } 
  else {
    DPrintf( BroadcastPrint, "CellWave::Probes  Warning.\n");
    DPrintf( BroadcastPrint, "  -->no probe output file specified--\n");
    DPrintf( BroadcastPrint, "     probes are output to showfile\n");
  }
  fflush( 0 );
}

void Probes::
collectData( int ktime, double tcomp, realCompositeGridFunction &u )
{
  DPrintf(DebugSolver,"..entering Probes::collectData(%d, %f)\n",
	  ktime, tcomp);

  if( useProbes && (ktime % probeHowOften ==0 )) {
    DPrintf(DebugSolver,".... collecting probes ....\n");
    int nProbes = probeLocation.getLength( axis1);
    int nSpecies= probeCurrentValue.getLength( axis2);
    assert( probeLocation.getLength( axis2 )>=3 );
    assert( nProbes == probeCurrentValue.getLength( axis1 ) );
    assert( nSpecies >= numberOfSpecies );

    this->probeTime= tcomp;
    interpolatePoints( probeLocation, u, probeCurrentValue); //from Overture
    //probeCurrentValue.display();
    this->writeData();
  }
}

void Probes::
writeData()
{
  DPrintf(DebugSolver,"..entering Probes::writeData\n");

  //assumes data has been collected by 'collectData'
  if( useProbes ) {
    int nProbes = probeLocation.getLength( axis1);
    int nSpecies= probeCurrentValue.getLength( axis2);
    assert( probeLocation.getLength( axis2 )>=3 );
    assert( nProbes == probeCurrentValue.getLength( axis1 ) );
    assert( nSpecies >= numberOfSpecies );
    assert( fProbe != NULL );
    
    fprintf(fProbe,"%16.8e ", probeTime );
    for(int iprobe=0; iprobe< nProbes; ++iprobe ) {
      for(int i=0; i<nSpecies; ++i ) {
	fprintf(fProbe," %16.8e ",probeCurrentValue( iprobe, i ));
      }
    }
    fprintf(fProbe,"\n");
    fflush(0 );
  }
}

void Probes::
closeOutputFile()
{
  DPrintf(DebugSolver,"..entering Probes::closeOutputFile\n");
  if( !isOK() ) {  
    DPrintf(DebugSolver,"..quick exit from Probes::closeOutputFile, not initialized\n");
    return;
  }

  assert( fProbe != NULL );
  fclose( fProbe );
  fProbe = NULL;
}

realArray & Probes::
getProbeValueSequence()
{
  return probeValueSequence;
}


realArray & Probes::
getProbeValues()
{
  return probeCurrentValue;
}
