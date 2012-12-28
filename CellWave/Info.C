///
/// Brief description: information/data class for CellWave
///   
///  CellWave::Info contains the main inputs for the simulator
///

#include "Info.h"

using namespace CellWave;

Info::
Info()
  : timeStepSize ( 0.001 ),
    numberOfTimeSteps( 100 ),
    saveEveryNthFrame ( 10 )
{
  // ..default
}
    
Info::
~Info()
{
  //..default
}

double Info::
maximumTime() 
{ 
  return ( timeStepSize*numberOfTimeSteps); 
};

