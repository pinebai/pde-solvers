//
//  CellWave::GenericCalciumIP3
//
//

#include "CellWave.h"
#include "GenericCalciumIP3Reaction.h"
#include "GenericReactionMacros.h"
#include <boost/tokenizer.hpp>
#include <vector>

using namespace CellWave;

//
//.. LOOP/Array interface
//

//.. looping calls to array interface, for initial data & rhs
void GenericCalciumIP3Reaction::
callInitialDataLoop( const double &tcomp,
		     const int&nd, const int&ncomp,
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
  // nd     = number of dimensions
  // ncomp  = number of components
  // nd?a, nd?b for ?=1,2,3 gives _array_ dimensions
  // n?a, n?b   for ?=1,2,3 gives _loop bounds_
  // ForGrid,GRIDINDEX:  defined in GenericReactionImpl.h

  //  printf("ReactionLiRinzelWagner--initial data %d\n", 
  //	 int(ip3Distribution));

#define RADIUS( x,y,z ) (sqrt( pow(x,2.) +pow(y,2.) +pow(z,2.)))
  maxLengthScale = 0.;
  
  if( nd ==3 ) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int cc=0, pc=1,hc=2;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];
      const double &z= xArray[ GRIDINDEX( i,j,k, zaxis) ];
      
      double &c      = qArray[ GRIDINDEX( i,j,k, cc) ];
      double &p      = qArray[ GRIDINDEX( i,j,k, pc) ];
      double &h      = qArray[ GRIDINDEX( i,j,k, hc) ];
      
      //call single point rxn code -- depends on impl    
      setInitialData( x,y,z, c,h,p);
      double r = RADIUS( x,y,z );
      if ( maxLengthScale < r ) maxLengthScale = r;
    }
  }
  else  if (nd == 2) {
    ForGrid(i,j,k) {
      const int xaxis=0, yaxis=1, zaxis=2;
      const int cc=0, pc=1,hc=2;
      
      const double &x= xArray[ GRIDINDEX( i,j,k, xaxis) ];
      const double &y= xArray[ GRIDINDEX( i,j,k, yaxis) ];

      const double z= 0.; // z must be zero 2D initial data
    
      double &c      = qArray[ GRIDINDEX( i,j,k, cc) ];
      double &p      = qArray[ GRIDINDEX( i,j,k, pc) ];
      double &h      = qArray[ GRIDINDEX( i,j,k, hc) ];
      
      //call single point rxn code -- depends on impl    
      setInitialData( x,y,z, c,h,p);
      double r = RADIUS( x,y, 0. );
      if ( maxLengthScale < r ) maxLengthScale = r;
    }
  }
# undef RADIUS
}

void GenericCalciumIP3Reaction::
initializeParameters()
{
  DPrintf(DebugReaction,"GenericCa2+/IP3 rxn -- initializeParameters.\n");

  ///..initial data
  calcium_0  = 0.1153;   /// micro M
  ip3_0      = 0.12;     /// micro M
  h_0        = 0.93;     /// nondimensional

  //numberOfDimensions = numDimensions_;
  numberOfDimensions = 3; // default dimension =3

  ///..initial blob data
  /// .... for Ca2+
  blobMin_c = calcium_0; 
  blobMax_c = 2.;
  xBlob_c = -500.;
  yBlob_c = 0.;
  zBlob_c = 0.;
  blobWidth_c = 60;
  
  /// .... for  IP3
  blobMin_p = ip3_0;
  blobMax_p = 20.;    /// CHANGED
  xBlob_p = -100.; 
  yBlob_p =100.;
  zBlob_p =0.;
  blobWidth_p = 200;
  
  /// .... for inhomog. distribution of IP3, fig. 4 in Wagner et al
  ip3Dist_I_s = 0.12;  /// u Molar
  ip3Dist_I_h = 1.0;    
  ///ip3Dist_I_h  =0.5; /// less convex wave
  ip3Dist_I_w = 0.015;
  ip3Dist_r_c = 500.; /// u meter
  
  ip3Dist_xstar    = 167;   /// u meter
  ip3Dist_Iprime_h = 0.84;  /// u Molar
  ip3Dist_Iprime_w = 0.8; 

  ip3scale = 1.;
  
  ///ip3Distribution = HomogeneousIP3;
  ///ip3Distribution = RadialIP3;
  //ip3Distribution = RadialWithGradientIP3;
  //ip3Distribution  = BlobIP3;
  
  if (ip3Distribution == RadialWithGradientIP3) {
    /// CONCAVE WAVE PARAMS -- FIGURE 6 in Wagner et al -- redefine some params
    //blobMax_c = calcium_0;  ///2.;
    //diffCalcium = 377.36; /// Derr = 20um^2/sec
    //lambda      = 250.;   /// 1/sec 
  }
}; 

bool GenericCalciumIP3Reaction::
readParameterFile( ParameterReader &par )
{
  bool ok=true;
  DPrintf(DebugReaction,"GenericCalciumIP3Rxn -- read param file\n");
  std::string typeString="";

  //..initial data
  par.get("number of dimensions", numberOfDimensions, 3); //number of spatial dimensions
  par.get("calcium_0", calcium_0,  0.1153); // micro M
  par.get("ip3_0",     ip3_0,      0.24);   // micro M             FIXME -- not as Wagner et al
  par.get("h_0",       h_0,        0.93);   // nondimensional

  ///..initial blob data
  std::string typeName="";
  //caDistribution = BlobCa;
  //ip3Distribution = HomogeneousIP3;

  caDistribution  = HomogeneousCa;
  ip3Distribution = RadialWithGradientIP3;

  par.get("calcium initial data type",typeName,"");
  if( typeName == "homogeneous" ) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: ca distribution = HOMOGENEOUS.\n");
    caDistribution= HomogeneousCa;
  }
  else   if( (typeName == "blob") || (typeName == "spherical blob" )) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: ca distribution = BLOB.\n");
    caDistribution=BlobCa;
  }
  else {
    DPrintf(DebugReaction,"Warning -- undefined 'calcium initial data type', set to homogeneous.\n");
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: ca distribution = UNKNOWN.\n");
    DPrintf(PrintOut,"Warning -- undefined 'calcium initial data type', set to homogeneous.\n");
    caDistribution= HomogeneousCa;
  }

  par.get("IP3 initial data type",typeName,"");
  if( typeName == "homogeneous" ) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: IP3 distribution = HOMOGENEOUS.\n");
    ip3Distribution= HomogeneousIP3;
  }
  else   if( (typeName == "blob") || (typeName == "spherical blob") ) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: IP3 distribution = BLOB.\n");
    ip3Distribution= BlobIP3;
  }
  else   if( typeName == "box blob" ) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: IP3 distribution = BOX BLOB.\n");
    ip3Distribution= BoxBlobIP3;
  }
  else   if( typeName == "radial" ) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: IP3 distribution = RADIAL.\n");
    ip3Distribution= RadialIP3;
  }
  else   if( typeName == "radial with gradient" ) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: IP3 distribution = RADIAL w/ GRAD.\n");
    ip3Distribution= RadialWithGradientIP3;
  }
  else   if( (typeName == "gaussian") || (typeName == "exact gaussian") ) {
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: IP3 distribution = EXACT GAUSSIAN.\n");
    ip3Distribution= ExactGaussianIP3;
  }
  else {
    DPrintf(DebugReaction,"Warning -- undefined 'IP3 initial data type', set to homogeneous.\n");
    DPrintf(DebugReaction,"GenericCa2+/IP3 rxn: IP3 distribution = UNKNOWN.\n");
    DPrintf(PrintOut,"Warning -- undefined 'IP3 initial data type', set to homogeneous.\n");
    ip3Distribution= HomogeneousIP3;
  }
  
  /// .... for Ca2+
  par.get("blobMin_c",    blobMin_c,    calcium_0 );
  par.get("blobMax_c",    blobMax_c,    2.);
  par.get("xBlob_c",      xBlob_c,     -200.);
  par.get("yBlob_c",      yBlob_c,      0.);
  par.get("zBlob_c",      zBlob_c,      0.);
  par.get("blobWidth_c",  blobWidth_c,  60.);
  
  /// .... for  IP3
  par.get("blobMin_p",    blobMin_p,    ip3_0);
  par.get("blobMax_p",    blobMax_p,    20.);
  par.get("xBlob_p",      xBlob_p,     -100.); 
  par.get("yBlob_p",      yBlob_p,      100.);
  par.get("zBlob_p",      zBlob_p,      0.);
  par.get("blobWidth_p",  blobWidth_p,  200);

  par.get("gaussian total concentration", gaussianTotalConcentration,  1. );
  par.get("gaussian time offset",         gaussianTimeOffset,          -0.1);

  par.get("IP3 diffusion", ip3diffusion,  0.); // IP3 diffusion: also RxnSlepchenko & RxnLiRinzel 

  /// .... for inhomog. distribution of IP3, fig. 4 in Wagner et al
  par.get("ip3Dist_I_s",ip3Dist_I_s,0.12);
  par.get("ip3Dist_I_h",ip3Dist_I_h,1.0);
  par.get("ip3Dist_I_w",ip3Dist_I_w,0.015);
  par.get("ip3Dist_r_c",ip3Dist_r_c,500.);
  
  par.get("ip3Dist_xstar", ip3Dist_xstar,167);
  par.get("ip3Dist_Iprime_h",ip3Dist_Iprime_h, 0.84);
  par.get("ip3Dist_Iprime_w",ip3Dist_Iprime_w,0.8);
  par.get("ip3scale", ip3scale, 1.0);
  if (ip3Distribution == BoxBlobIP3 ) {
    std::string cornerString;
    par.get("IP3 box corners", cornerString,"");
    par.get("IP3 box maximum", ip3BoxMax, ip3_0);

    //..using Boost 1.30
    typedef boost::tokenizer<boost::char_separator<char> >            Tokenizer;
    typedef std::vector<double>                                       DoubleVector;
    typedef boost::tokenizer<boost::char_separator<char> >::iterator  TokIterator;

    boost::char_separator<char> sep(", ");
    Tokenizer tok(cornerString, sep);
    DoubleVector corners;

#if 0    
    //..old: for Boost 1.26
    typedef boost::tokenizer<>::iterator TokIterator;
    typedef std::vector<double>          DoubleVector;
    
    boost::tokenizer<>   tok(cornerString);
    std::vector<double>  corners;
#endif   

    DPrintf(DebugReaction,"    IP3 Box Blob\n");
    DPrintf(DebugReaction,"      IP3 min= %f, IP3 max=%f\n",
	    ip3_0, ip3BoxMax);
    DPrintf(DebugReaction,"Tokenizing corners '%s'\n", cornerString.c_str());
    DPrintf(DebugReaction,"      corners are = ");
    corners.clear();
    for (TokIterator it=tok.begin(); it!=tok.end(); ++it) {
      double q;
      printf("");
      sscanf(it->c_str(), "%lf", &q);
      DPrintf(DebugReaction," %g ['%s'];  ", q, it->c_str());
      corners.push_back( q );
    }
    DPrintf(DebugReaction,"\n");
    
    assert( corners.size() > 5 );
    setBoxCorners(corners[0], corners[1], corners[2],
		  corners[3], corners[4], corners[5]);
  }

  return( ok );
}

void GenericCalciumIP3Reaction::
printParameters(int iPrintChannel)
{
  const int ip=iPrintChannel;

  //..initial data
  printOneParameter(ip,"number of dimensions", numberOfDimensions, 3); 
  printOneParameter(ip,"calcium_0", calcium_0,  0.1153); 
  printOneParameter(ip,"ip3_0",     ip3_0,      0.24);   
  printOneParameter(ip,"h_0",       h_0,        0.93);   

  //par.get("calcium initial data type",typeName,"");
  std::string typeName="";
  if(     caDistribution== HomogeneousCa ) {
    typeName = "homogeneous";
  }
  else   if(     caDistribution==BlobCa ){
    typeName = "spherical blob";
  }
  printOneParameter(ip,"calcium initial data type", typeName, "homogeneous");

  if( ip3Distribution== HomogeneousIP3 ) {
    typeName = "homogeneous";
  }
  else   if(     ip3Distribution== BlobIP3 ) {
    typeName = "spherical blob";
  }
  else   if(     ip3Distribution== BoxBlobIP3 ) {
    typeName = "box blob";
  }
  else   if(     ip3Distribution== RadialIP3 ) {
    typeName = "radial";
  }
  else   if(     ip3Distribution== RadialWithGradientIP3 ) {
    typeName = "radial with gradient";
  }
  else   if(     ip3Distribution== ExactGaussianIP3 ) {
    typeName = "exact gaussian";
  }
  printOneParameter(ip,"IP3 initial data type",typeName,"homogeneous");
  
  printOneParameter(ip,"blobMin_c",    blobMin_c,    calcium_0 );
  printOneParameter(ip,"blobMax_c",    blobMax_c,    2.);
  printOneParameter(ip,"xBlob_c",      xBlob_c,     -200.);
  printOneParameter(ip,"yBlob_c",      yBlob_c,      0.);
  printOneParameter(ip,"zBlob_c",      zBlob_c,      0.);
  printOneParameter(ip,"blobWidth_c",  blobWidth_c,  60.);
  
  printOneParameter(ip,"blobMin_p",    blobMin_p,    ip3_0);
  printOneParameter(ip,"blobMax_p",    blobMax_p,    20.);
  printOneParameter(ip,"xBlob_p",      xBlob_p,     -100.); 
  printOneParameter(ip,"yBlob_p",      yBlob_p,      100.);
  printOneParameter(ip,"zBlob_p",      zBlob_p,      0.);
  printOneParameter(ip,"blobWidth_p",  blobWidth_p,  200);

  printOneParameter(ip,"gaussian total concentration", gaussianTotalConcentration,  1. );
  printOneParameter(ip,"gaussian time offset",         gaussianTimeOffset,          -0.1);

  printOneParameter(ip,"IP3 diffusion", ip3diffusion,  0.);

  printOneParameter(ip,"ip3Dist_I_s",ip3Dist_I_s,0.12);
  printOneParameter(ip,"ip3Dist_I_h",ip3Dist_I_h,1.0);
  printOneParameter(ip,"ip3Dist_I_w",ip3Dist_I_w,0.015);
  printOneParameter(ip,"ip3Dist_r_c",ip3Dist_r_c,500.);
  
  printOneParameter(ip,"ip3Dist_xstar", ip3Dist_xstar,167);
  printOneParameter(ip,"ip3Dist_Iprime_h",ip3Dist_Iprime_h, 0.84);
  printOneParameter(ip,"ip3Dist_Iprime_w",ip3Dist_Iprime_w,0.8);
  printOneParameter(ip,"ip3scale", ip3scale, 1.0);

}
