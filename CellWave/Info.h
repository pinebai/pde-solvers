///
/// Brief description: information/data class for CellWave
///   
///  CellWaveInfo contains the main inputs for the simulator
///
#ifndef INFO_H
#define INFO_H "Info.h"

#include <iostream>
#include <string>

namespace CellWave {
  
  class Info {
  public:
    Info();
    ~Info();
  private: // don't use the copy constructors
    Info( const Info &X) {};
    Info &operator=( const Info &X ) { };

  public:
    double maximumTime();

    // ..data
    int  numberOfSpecies; // number of reacting species
    double *viscosity;      // diffusion const of the species (no cross) 

    //aString nameOfOGFile, nameOfShowFile;
    std::string nameOfOGFile, nameOfShowFile;

    int  saveEveryNthFrame;// =10;
    int  logEveryNthFrame; // =10;
    int  flushFrequency;  //=10
    double timeStepSize;//= 0.001;
    int  numberOfTimeSteps;//=100;

  }; //end class Info
}; //end namespace 
#endif
