create mappings 
  annulus 
    start and end angles 
    0 .3
    inner and outer radii 
    5. 10.
    boundary conditions 
    1 1 1 2 
    mappingName 
    bottom 
    exit 
  annulus 
    start and end angles 
    0 0.3 
    inner and outer radii 
    10. 20. 
    lines 
    21 13 
    boundary conditions 
    1 1 2 2 
    mappingName 
    middle
    exit 
  annulus 
    start and end angles 
    0 0.3 
    inner and outer radii 
    20. 30. 
    lines 
    46 15
    boundary conditions 
    1 1 2 1 
    mappingName 
    top
    exit 
  view mappings 
    bottom 
    middle
    top 
    exit 
  exit this menu
generate an overlapping grid 
  bottom 
  middle
  top 
  change parameters 
    prevent hole cutting 
      all 
      all 
      done 
    prevent interpolation
      all
      all
      done
    exit
  compute overlap
  exit
save a grid
matchingAnnuli3x.hdf
matchingAnnuli
exit
