create mappings 
  rectangle 
    set corners 
    0 1 0 .5 
    lines 
    11 5 
    mappingName 
    upperGrid 
    exit 
  rectangle 
    set corners 
    0 1 -.5 0 
    lines 
    21 11 
    mappingName 
    lowerGrid 
    exit 
  exit this menu 
* 
generate an overlapping grid 
  upperGrid 
  lowerGrid 
  change parameters 
    mixed boundary 
      upperGrid 
        bottom (side=0,axis=1) 
        lowerGrid 
          specify mixed boundary points 
          5 7  -1 -1 
          done 
      * 
      lowerGrid 
        top    (side=1,axis=1) 
        upperGrid 
          specify mixed boundary points 
          10 14 11 11 
          done 
      done 
    non-cutting boundary points 
      upperGrid 
        bottom (side=0,axis=1) 
      0 10 0 0 0 0 
      lowerGrid 
        top    (side=1,axis=1) 
      0 20 10 10 0 0 
      done 
      exit 
    compute overlap 
    exit
  save a grid
  billMatchingSquares.hdf
  matchingSquares
  exit
