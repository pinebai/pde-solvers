#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "OGTrigFunction.h"
#include "OGPolyFunction.h"
#include "PlotStuff.h"
#include "interpolatePoints.h"

#include "CellWave.h"

#include <math.h>
#include "GenericReactionMacros.h"

#include "heatKernel.h"

GaussianKernelData 
setGaussianKernelData( double timeOffset, 
		       double viscosity, 
		       double x0, 
		       double y0, 
		       double z0,
		       double totalMass,
		       int numberOfDimensions)
{
  GaussianKernelData kernelData;
  kernelData.timeOffset = timeOffset;
  kernelData.viscosity  = viscosity;
  kernelData.x0         = x0;
  kernelData.y0         = y0;
  kernelData.z0         = z0;
  kernelData.totalMass  = totalMass;
  kernelData.numberOfDimensions = numberOfDimensions;

  return( kernelData );
}

void
printGaussianKernelData( int output, const GaussianKernelData &kernelData )
{
  using CellWave::DPrintf;

  DPrintf(output,"GaussianKernelData:\n");
  DPrintf(output,"   timeOffset = %f\n", kernelData.timeOffset );
  DPrintf(output,"   viscosity  = %f\n", kernelData.viscosity );
  DPrintf(output,"   x0         = %f\n", kernelData.x0 );
  DPrintf(output,"   y0         = %f\n",  kernelData.y0 );     
  DPrintf(output,"   z0         = %f\n",  kernelData.z0 );     
  DPrintf(output,"   mass       = %f\n",  kernelData.totalMass );     
  DPrintf(output,"   numberOfDimensions= %d\n",  kernelData.numberOfDimensions );
}

void 
heatKernelFortranArray(const GaussianKernelData &kernelData,
		       const double &tcomp,
		       const int&ncomp,
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
  const int nd=kernelData.numberOfDimensions;

  if (nd==1 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int qc=0;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      //const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      //const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      const double y=0.,  z=0.;

      double &q      = qArray[ GRIDINDEX( i,j,k, qc) ];
      
      q = evaluateGaussianKernel( kernelData, tcomp, x,y,z);
    }
  }
  else if (nd==2 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int qc=0;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      //const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      const double z=0.;

      double &q      = qArray[ GRIDINDEX( i,j,k, qc) ];
      
      q = evaluateGaussianKernel( kernelData, tcomp, x,y,z);
    }
  }
  else   if (nd==3 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int qc=0;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      
      double &q      = qArray[ GRIDINDEX( i,j,k, qc) ];
      
      q = evaluateGaussianKernel( kernelData, tcomp, x,y,z);
    }
  }
  else printf("--unknown dimension = %d -- error\n", nd);
}


void
heatKernelGridFunction( GaussianKernelData &kernelData, 
			double tcomp, 
			MappedGrid        &mg, 
			realArray               &heatKernel )
{
  realArray & x = mg.vertex();  // array of vertices
  const IntegerArray & d    =  mg.dimension();
  const IntegerArray & gir  =  mg.gridIndexRange();
  const int nd=mg.numberOfDimensions();
  const int ncomp=1; //only set first component
  
  kernelData.numberOfDimensions = nd;
  // call a fortran function to compute du/dt
  // (This function does not currently solve the convection diffusion equation)
  //mySolver( t,dt,a,b,nu,nd, d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2),
  //           gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
  //          *x.getDataPointer(),*ug.getDataPointer(), *dudtg.getDataPointer() );
  heatKernelFortranArray( kernelData, tcomp,
			  ncomp,
			  d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2),
			  gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
			  *x.getDataPointer(), *heatKernel.getDataPointer() );
			  			  
}

