#ifndef HEAT_KERNEL_H
#define HEAT_KERNEL_H

#include "Overture.h"
#include <math.h>

//
// ..Gaussian Kernel routines
//
struct GaussianKernelData {
  int    numberOfDimensions;
  double timeOffset;
  double viscosity;
  double x0, y0, z0;
  double totalMass;

};

GaussianKernelData 
setGaussianKernelData( double timeOffset, 
		       double viscosity, 
		       double x0, 
		       double y0, 
		       double z0,
		       double totalMass,
		       int numberOfDimensions);

void
printGaussianKernelData( int output, const GaussianKernelData &kernelData );

inline double 
evaluateGaussianKernel( const GaussianKernelData &kernelData, double tcomp, double xc, double yc, double zc )
{
  const double t=tcomp + kernelData.timeOffset;
  const double x=xc    - kernelData.x0;
  const double y=yc    - kernelData.y0;
  const double z=zc    - kernelData.z0;
  const double mu   = kernelData.viscosity;
  const int    ndim = kernelData.numberOfDimensions;  
  const double totalMass = kernelData.totalMass;
  const double r2= x*x + y*y + z*z;
  const double pi= 4.*atan( 1.);

  double q;
  
  q = totalMass*exp( - r2/(4. * mu * t ))*pow( 4.* mu * pi * t, -ndim/2.);

  //  printf(" Gaussian=%8.3f, t=%8.3f, mu=%8.3f, r2=%8.3f, ndim/2.=%8.3f\n",
  //           q, t,mu,r2,ndim/2.);
  return( q );
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
		       double &qfirst );

void
heatKernelGridFunction( GaussianKernelData &kernelData, 
			double              tcomp, 
			MappedGrid         &mg, 
			realArray          &heatKernel );

#endif

