create mappings 
  annulus 
    start and end angles 
    0 .4 
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
    1 1 2 1 
    exit 
  view mappings 
    exit 
  change a mapping 
  Annulus 
    mappingName 
    top 
    exit 
  view mappings 
    bottom 
    top 
    exit 
  change a mapping 
  bottom 
    start and end angles 
    0 0.3 
    exit 
  exit this menu 
generate an overlapping grid 
  bottom 
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
matchingAnnuli20um.hdf
matchingAnnuli
exit
