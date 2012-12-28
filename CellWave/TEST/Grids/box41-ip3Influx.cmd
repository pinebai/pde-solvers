create mappings
  box
    set corners
*    -1 1 -1 1 -1 1
    0 20 0 20 0 20
    lines
    41 41 41 
    boundary conditions 
    1 1 1 1 1 3
*    1 2 3 4 5 6
    exit
  exit this menu
generate an overlapping grid
  box
  compute overlap
  exit
save a grid
box41-ip3Influx.hdf
box
exit
