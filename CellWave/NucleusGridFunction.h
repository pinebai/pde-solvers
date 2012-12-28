#ifndef NUCLEUS_GRID_FUNCTION_H
#define NUCLEUS_GRID_FUNCTION_H

//
// NucleusGridFunction -- more than one nucleus/grid
//  .. contain maps for
//     1)  gridID --> nucleusIndices
//     2)  nucleusIndex --> gridIDs 
//   ( 3)  cellID --> nucleusIndices, for later )

#include <vector>
#include <map>

//..overture
#include "Overture.h"
#include "CompositeGridFunction.h"

//..cellwave
#include "Nucleus.h"
#include "ParameterReader.h"


namespace CellWave {

class NucleusGridFunction {
 public:
  //..public datatypes
  typedef std::vector<Nucleus>      NucleusVector;
  typedef std::vector<int>          IDVector;
  typedef std::multimap<int, int>   GridNucleusMap; // from gridID        --> nucleus Index
  typedef std::multimap<int, int>   NucleusGridMap; // from nucleus Index --> gridID
  typedef std::multimap<int, int>::iterator IterateIntMap;
  
  //..public interface
  NucleusGridFunction();
  ~NucleusGridFunction();
  
  void updateToMatchGrid( CompositeGrid &cg,
			  doubleCompositeGridFunction &nucleusGF);
  void readParameterFile( ParameterReader &params );

  bool readCellNucleusFile( const std::string cn_filename,
			    const std::string gridFileName = "" );
  void evaluateGridFunction( double tcomp =0.);
  void maskAGridFunction( doubleCompositeGridFunction &u);
  void maskAGridFunction( int gridNumber, realArray &ug, Index &I1, Index &I2, Index &I3);
  void printInfo();
  bool checkConsistency();
  
  int getNumberOfNuclei() { return( nucleus.size()); };
  int getGridNuclei( int gridID, IDVector &nuclei);
  int getNucleusGrids( int nucleusID, IDVector &gridIDs);
  int getAllNuclei( IDVector &nuclei);

  void clearGrid2NucleusMap() { grid2NucleusMap.clear();  };
  void clearNucleus2GridMap() { nucleus2GridMap.clear();  };

  Nucleus &getNucleus( int id ) {return nucleus[id-1]; };
  void     setNucleus( int nucleusID, const Nucleus &nucleus);

  //..data
  
  NucleusVector  nucleus;           // the list of nucleai
  GridNucleusMap grid2NucleusMap;   // map a grid # to nucleus Indices
  NucleusGridMap nucleus2GridMap;   // map the nucleus Index to grid #

  CompositeGrid               *pCG;
  doubleCompositeGridFunction *pNucleusGF;

  double nucleusBoundaryThickness;
};

} //end namespace CellWave

#endif
