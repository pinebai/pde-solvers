#ifndef PROBES_H
#define PROBES_H
//
//  CellWave::Probes -- place measurement probes inside the silicon cell to record
//                      local chemical concentrations
//    $Id: Probes.h,v 1.4 2003/05/21 17:56:15 pfast Exp $
//

#include "Overture.h"
#include <string>
#include <vector>

namespace CellWave {

  struct OneProbe {
    OneProbe( int id_, double x_, double y_, double z_ =0. )
    { id=id_; x=x_; y= y_; z=z_; };

    int id;
    double x;
    double y;
    double z;
  };

  class Probes {
  public:
    Probes( int numberOfSpecies_=0 );
    ~Probes();

    bool isOK() { return useProbes; };

    bool readParameterFile( CellWave::ParameterReader &param);
    void outputHeader( bool useHeader=true) { noHeader= (!useHeader);} ;

    void initializeProbeValues( int numberOfSpecies_ );
    void readProbeLocations();

    void openOutputFile(std::string comment=""); //write the probe locations, in comment lines
    void collectData( int ktime, double tcomp, realCompositeGridFunction &u );
    void writeData();
    //void writeAll();
    void closeOutputFile();

    int  numberOfProbes() { return ( probeLocation.getLength( axis1 ));};
    void clearProbeLocations();

    realArray & getProbeValueSequence();
    realArray & getProbeValues();

    //..data
    std::string probeOutputFileName, probeLocationFileName;
    std::string parameterFileName;
    FILE *fProbe; //output file

    bool useProbes;
    bool noHeader;
    int  probeHowOften;
    int numberOfSpecies;      //how many chemicals are we probing

    realArray probeLocation; //(numberOfProbes, 3);
    double    probeTime;
    realArray probeCurrentValue;    //(numberOfProbes, 3);
    realArray probeValueSequence; //for saving to showfile
  };


}; //end of namespace CellWave
#endif
