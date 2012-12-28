//
// convert a CellWave parameter file to Matlab parameter file
//    for computing nullClines & IP3 open probability
//

#include <iostream>
#include <string>
#include <stdio.h>

#include "CellWave.h"
#include "Reaction2Buffer.h"

int main(int argc, char **argv)
{
  CellWave::start(argc,argv);

  using CellWave::DPrintf;
  const int BPRINT=CellWave::BroadcastPrint;
  const int PRINT=CellWave::PrintOut;

  CellWave::ParameterReader *pParameters = NULL;

  //*
  //*..Set parameters
  //*
  std::string paramFileName = "";
  std::string outputFileName="matlabCellWaveParameters.m";
  if (argc>1) {
    assert( argv[1] != NULL );
    paramFileName = argv[1];
    pParameters = new CellWave::ParameterReader( paramFileName );
  } else {
    printf("usage: %s <CW param file> <optional: matlab output file>\n",argv[0]);
    return -1;
  }
  if (argc>2) {
    assert( argv[2] != NULL );
    outputFileName=argv[2];
    //outputFileName += ".m";
  }
  assert( pParameters != NULL );
  CellWave::ParameterReader &param = *pParameters;

  CellWave::Reaction2Buffer chem;
  chem.readParameterFile( param );

  //
  // -- generate matlab output file
  //
  FILE *fp = fopen( outputFileName.c_str(), "w+");
  if( fp == NULL ) {
    printf("ERROR -- could not open '%s' for writing.\n", 
	   outputFileName.c_str());
    return -1;
  }

  //..write: header
  fprintf( fp, "function params = parameters2Buffer\n");
  fprintf( fp, "%% function params = parameters2Buffer\n");
  fprintf( fp, "%% \n");
  fprintf( fp, "%% \n");
  fprintf( fp, "%% set  parameters for Li-Rinzel 2 Buffer model\n");
  fprintf( fp, "%% \n");
  fprintf( fp, "%% .. automatically generated from a CellWave param file:\n");
  fprintf( fp, "%%    original file= '%s'\n", paramFileName.c_str() );
  fprintf( fp, "\n");

  //..write: parameters
  fprintf( fp, "%24s = %e; \n","params.nu_L",    chem.nu_L);
  fprintf( fp, "%24s = %e; \n","params.nu_c",    chem.nu_c);
  fprintf( fp, "%24s = %e; \n","params.d_act",   chem.d_act);
  fprintf( fp, "%24s = %e; \n","params.d_inh",   chem.d_inh);
  fprintf( fp, "%24s = %e; \n","params.d_I",     chem.d_I);
  fprintf( fp, "%24s = %e; \n","params.V_m",     chem.nu_m);
  fprintf( fp, "%24s = %e; \n","params.k_p",     chem.K_p);
  fprintf( fp, "%24s = %e; \n","params.C_er",    chem.C_er);
  fprintf( fp, "%24s = %e; \n","params.k_on",    chem.k_on);
  fprintf( fp, "\n");

  fprintf( fp, "%24s = %e; \n","params.gamma",   chem.gamma );
  fprintf( fp, "%24s = %e; \n","params.b1tot",   chem.B1_tot);
  fprintf( fp, "%24s = %e; \n","params.b2tot",   chem.B2_tot);
  fprintf( fp, "%24s = %e; \n","params.K1",      chem.K1); 
  fprintf( fp, "%24s = %e; \n","params.K2",      chem.K2); 
  fprintf( fp, "%24s = %e; \n","params.H",       chem.H_flux);

  fprintf( fp, "%24s = %e; \n","params.c_0",     chem.calcium_0);
  fprintf( fp, "%24s = %e; \n","params.ip3_0",   chem.ip3_0);
  fprintf( fp, "%24s = %e; \n","params.h_0",     chem.h_0);
  fprintf( fp, "%24s = %e; \n","params.b1_0",    chem.b1_0);
  fprintf( fp, "%24s = %e; \n","params.b2_0",    chem.b2_0);
  fprintf( fp,"\n");

  CellWave::finish();
}
