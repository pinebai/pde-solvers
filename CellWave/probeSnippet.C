#include "interpolatePoints.h"

  FILE *fProbe = fopen("oxcp.dat", "w+");
    // (7) dump out solution values at a few probes
    int numberOfProbes = 10;
    realArray probeLocation(numberOfProbes, 3);
    realArray probeValue(numberOfProbes,3);
    real dx=1000./(numberOfProbes+1);
    for( int i=0; i< numberOfProbes; ++i ) {
      probeLocation( i, axis1) = -500. + (i+1)*dx;
      probeLocation( i, axis2) = 0.;
      probeLocation( i, axis3) = 0.;
    }
    interpolatePoints( probeLocation, q, probeValue);

    fprintf(fProbe, "%8.4e  ",tcomp);
    for( int i=0; i< numberOfProbes; ++i) {
      fprintf(fProbe, "%8.4e %8.4e %8.4e ", 
              probeValue(i,0), 
              probeValue(i,1),
              probeValue(i,2));
    }
    fprintf(fProbe, "\n");
    fflush(0);

