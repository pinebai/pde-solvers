//
// test using CVode in Overture
//

#include "Overture.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"  /* definitions of realtype, integertype           */
#include "cvode.h"          /* main CVODE header file                         */
#include "iterativ.h"       /* contains the enum for types of preconditioning */
#include "cvspgmr.h"        /* use CVSPGMR linear solver each internal step   */
#include "smalldense.h"     /* use generic DENSE solver for preconditioning   */
#include "nvector_serial.h" /* definitions of type N_Vector, macro NV_DATA_S  */
#include "sundialsmath.h"   /* contains SQR macro                             */

struct CVodeInfo {
  M_Env machEnv;
  realtype abstol, reltol, t, tout, ropt[OPT_SIZE];
  realtype floor;
  long int iopt[OPT_SIZE];
  N_Vector y;
  void *cvode_mem;
  int iout, flag;
  int neq;
};

struct ProblemConstants {
  int NUM_SPECIES;
  realtype KH;
  realtype VEL;
  realtype KV0;
  realtype Q1; 
  realtype Q2 ;
  realtype C3 ;
  realtype A3 ;
  realtype A4 ;
  realtype C1_SCALE; 
  realtype C2_SCALE ;
  
  realtype T0       ;
  int      NOUT      ;
  realtype TWOHR    ;
  realtype HALFDAY  ;
  realtype PI       ;
  
  realtype XMIN      ;
  realtype XMAX      ;
  realtype ZMIN     ;
  realtype ZMAX     ;
  realtype XMID     ;
  realtype ZMID     ;

  int     MX        ;
  int     MZ        ;
  int     NSMX      ;
  int     MM        ;
};

struct UserData {
  CVodeInfo         *ci;
  ProblemConstants  *pconst;
}

void setupCVodeInfo( int neq, CVodeInfo &ci)
{
  ci.neq=neq;
  /* Initialize serial machine environment */
  ci.machEnv = M_EnvInit_Serial(ci.neq);
  ci.reltol = 1e-5;
  ci.floor  = 100.0;
  ci.abstol = ci.reltol*ci.floor;
  ci.y = N_VNew(ci.neq, ci.machEnv);
}

void setupProblemConstants( ProblemConstants &ci )
{
  ci.NUM_SPECIES =2            ;/* number of species         */
  ci.KH          =4.0e-6       ;/* horizontal diffusivity Kh */
  ci.VEL         =0.001        ;/* advection velocity V      */
  ci.KV0         =1.0e-8       ;/* coefficient in Kv(z)      */
  ci.Q1          =1.63e-16     ;/* coefficients q1, q2, c3   */ 
  ci.Q2          =4.66e-16;
  ci.C3          =3.7e16;
  ci.A3          =22.62        ;/* coefficient in expression for q3(t) */
  ci.A4          =7.601        ;/* coefficient in expression for q4(t) */
  ci.C1_SCALE    =1.0e6        ;/* coefficients in initial profiles    */
  ci.C2_SCALE    =1.0e1;
  
  ci.T0          =0.0          ;/* initial time */
  ci.NOUT        =12           ;/* number of output times */
  ci.TWOHR       =7200.0       ;/* number of seconds in two hours  */
  ci.HALFDAY     =4.32e4       ;/* number of seconds in a half day */
  ci.PI          =4.*atan(1.); ;/* pi */ 
  
  ci.XMIN        = 0.0         ;/* grid boundaries in x  */
  ci.XMAX        =20.0         ; 
  ci.ZMIN        =30.0         ;/* grid boundaries in z  */
  ci.ZMAX        =50.0;
  ci.XMID        =10.0         ;/* grid midpoints in x,z */          
  ci.ZMID        =40.0;
  
  ci.MX          =400          ;  /* MX = number of x mesh points */
  ci.MZ          =400          ;  /* MZ = number of z mesh points */
  ci.NSMX        =800          ;  /* NSMX = NUM_SPECIES*MX */
  ci.MM          = ci.MX*ci.MZ; /* MM = MX*MZ */
}

void setupCVode( CVodeInfo &ci, ProblemConstants &cc, UserData &data  )
{
  cvode_mem = CVodeMalloc(ci.neq, f, cc.T0, ci.y, BDF, NEWTON, SS, &ci.reltol,
                          &ci.abstol, data, NULL, FALSE, ci.iopt, ci.ropt, ci.machEnv);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed."); return(1); }
}

int main( int argc, char **argv)
{
  Overture::start( argc, argv );

  CVodeInfo cvodeInfo;
  ProblemConstants pconst;

  setupProblemConstants( pconst );

  /* Allocate memory, and set problem data, initial values, tolerances */ 
  UserData data;
  data = AllocUserData();
  InitUserData(data,cvodeInfo,pconst);
  SetInitialProfiles(y, data->dx, data->dz);

  /* Call CVodeMalloc to initialize CVODE: 

     NEQ     is the problem size = number of equations
     f       is the user's right hand side function in y'=f(t,y)
     T0      is the initial time
     y       is the initial dependent variable vector
     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration
     SS      specifies scalar relative and absolute tolerances
     &reltol and &abstol are pointers to the scalar tolerances
     data    is the pointer to the user-defined block of coefficients
     FALSE   indicates there are no optional inputs in iopt and ropt
     iopt    and ropt arrays communicate optional integer and real input/output

     A pointer to CVODE problem memory is returned and stored in cvode_mem.  */

  setupCVode( cvodeInfo, pc, data );

  Overture::finish();
}
