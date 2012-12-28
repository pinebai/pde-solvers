
#include <stdio.h>
#include <assert.h>

#include "Overture.h"
#include "CellWave.h"

#include "Probes.h"

int main(int argc, char **argv)
{
  Overture::start(argc,argv);  // initialize Overture
  CellWave::start(argc,argv);
  
  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  DPrintf(BPRINT," ------------------------------------------------------------ \n");
  DPrintf(BPRINT," Test CellWave::Probes -- probe a grid function \n");
  DPrintf(BPRINT," ------------------------------------------------------------ \n");

  if (argc<2) {
    printf("  usage: %s <parameter file>\n", argv[0]);
    exit(-1);
  }
  std::string parameterFileName=argv[1];
  CellWave::ParameterReader param( parameterFileName );

  std::string overtureGridFileName;
  param.get( "name of grid file", overtureGridFileName, "Grids/oocyte_1000_um.hdf");
  
  
  // create and read in a CompositeGrid
  aString nameOfOGFile = overtureGridFileName.c_str();
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update();

  double tcomp = 0.;
  int nSpecies=3;
  Range all;
  realCompositeGridFunction u(cg,all,all,all,nSpecies);
  u=0.;                                   // initialize to zero
  Index I1,I2,I3;                         // A++ Index object

  CellWave::Probes probes (nSpecies);
  probes.outputHeader(false);
  probes.readParameterFile( param );
  //probes.readProbeLocations();
  probes.openOutputFile( "$Id: testProbes.C,v 1.2 2003/02/06 06:15:25 pfast Exp $" );
    
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    double pi = 4.*atan(1.);
    getIndex(cg[grid].indexRange(),I1,I2,I3);                

    u[grid](I1,I2,I3,axis1)=cg[grid].vertex()(I1,I2,I3,axis1);
    u[grid](I1,I2,I3,axis2)=cg[grid].vertex()(I1,I2,I3,axis2);
    if (cg.numberOfDimensions()>2) {
      u[grid](I1,I2,I3,axis3)=cg[grid].vertex()(I1,I2,I3,axis3);
    }
    else {
      u[grid](I1,I2,I3,axis3)=0.;
    }

    //    where( cg[grid].mask()(I1,I2,I3) > 0 ) {
    //  u[grid](I1,I2,I3)=sin(2.*pi*cg[grid].vertex()(I1,I2,I3,axis1)/500.)
    //                   *cos(2.*pi*cg[grid].vertex()(I1,I2,I3,axis2)/500.);  
    // }
    //for(int i=0; i< nSpecies; ++i ) {
    //  u[grid](I1,I2,I3,i) = 1.0*i+1.;
    //}
  }    
  Interpolant interpolant(cg);      // Make an interpolant
  interpolant.interpolate(u);       // interpolate

  int ktime=0;
  probes.collectData( ktime, tcomp, u );
  probes.closeOutputFile();

  CellWave::finish();
  Overture::finish();          
  return 0;  

}
