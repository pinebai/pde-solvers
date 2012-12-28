//..Overture headers
#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "display.h"

//..CVode headers
#include "sundialstypes.h"
#include "cvode.h"
#include "cvdense.h"
#include "nvector_serial.h"
#include "dense.h"
 
#include "iterativ.h"       /* contains the enum for types of preconditioning */
#include "cvspgmr.h"        /* use CVSPGMR linear solver each internal step   */
#include "cvbandpre.h"      /* band preconditioner function prototypes        */
#include "sundialsmath.h"   /* contains SQR macro                             */

//void (*CVodeRHSFunction)(integertype N, 
//			 realtype t, 
//			 N_Vector y, 
//			 N_Vector ydot,
//			 void *f_data);

class CVodeData {
public:
  CVodeData() { setBDFSolver(); setNewtonSolver(); };
  ~CVodeData() { };

  void setBDFSolver()     { solverType = BDF; };
  void setAdamsSolver()   { solverType = ADAMS; };
  void setNewtonSolver()  { nonlinearSolverType = NEWTON; };
  void setFunctionalIterationSolver() { nonlinearSolverType = FUNCTIONAL; };

  long int  getSolver() { return( solverType );};
  long int  getNonlinearSolver() {return( nonlinearSolverType); };

  void setNumberOfEquations( int neq ) { numberOfEquations =neq;};
  int  getNumberOfEquations() { return numberOfEquations; };

  int       numberOfEquations; //NEQ;
  realtype  absoluteTolerance; //ATOL 
  realtype  relativeTolerance; //RTOL

  N_Vector  solutionVector;

  M_Env     machEnv;
  void     *cvode_mem;
  realtype ropt[OPT_SIZE], reltol, tcomp, tout, dtout;
  long int iopt[OPT_SIZE];
  N_Vector abstol;

  int       iout;
  int       solverFlag;

  realtype  time0;
  long int  solverType;          // BDF or ADAMS
  long int  nonlinearSolverType; //
};

class CVodeOvertureGrid { // ..data structure passed into RHS/Jac functions. Contains all CVode info
  //struct CVodeData {
  // int i;
  //};
public:
  CVodeOvertureGrid();
  ~CVodeOvertureGrid();

  void initialize(int numberOfComponents, MappedGrid &mg, realArray &uSolution );

  void setRHSFunction(   void (*rhsFunction)(integertype N, 
		      realtype t, 
		      N_Vector y, 
		      N_Vector ydot,
					     void *f_data));

  void setMappedGrid(MappedGrid &mg );

  void setGridFunction( int numberOfComponents, realArray &uSolution  );

  real *getArrayDataPointer( const realArray &u );

  void computeOvertureDimensions( );

  void setGridNumber( int ig ) { igrid = ig; };
  int  getGridNumber() { return igrid; };
  
  //..Overture data
  int     igrid;
  int     nd1a, nd1b, nd2a,nd2b,nd3a,nd3b; //dimensions
  int     n1a, n1b, n2a, n2b, n3a, n3b, ndimensions, ncomponents;
  int     getNumberOfGridPoints() 
    { return (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)*ncomponents; };

  MappedGrid *pMg;
  
  //..CVode data
  CVodeData cvodeData;
  void (*rhsFunction)(integertype N, 
		      realtype t, 
		      N_Vector y, 
		      N_Vector ydot,
		      void *f_data);

};

class CVodeOverture 
{
public:
  CVodeOverture( );
  ~CVodeOverture();
  void initialize( CompositeGrid &cg, realCompositeGridFunction &q);
  
  //..data

  typedef CVodeOvertureGrid*   PCVodeOvertureGrid;
  int     numberOfGrids, numberOfComponents;

  CVodeOvertureGrid *cvodeGrid;
};


