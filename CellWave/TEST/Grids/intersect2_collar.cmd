create mappings
  annulus
    inner and outer radii
    250 500
    lines
    60 20
    lines
    60 15
    lines
    80 15
    centre for annulus
    -200 -100
    boundary conditions
    -1 -1 0 1
    mappingName
    circle1
    exit
  annulus
    inner and outer radii
    250 500
    lines
    60 15
    boundary conditions
    -1 -1 0 1
   centre for annulus
    100 0
    mappingName
    circle2
    exit
  rectangle
    set corners
    -500 400 -400 300
    lines
    90 70
    mapping parameters
    lines 80 50
    close mapping dialog
    mapping parameters
    Boundary Condition: left    0
    Boundary Condition: right   0
    Boundary Condition: bottom  0
    Boundary Condition: top     0
    close mapping dialog
    exit

*
*
*
  exit this menu
generate an overlapping grid
  square
  circle1
  circle2
  change parameters
    prevent hole cutting
      circle1
      all
      circle2
      all
      done
    ghost points
      all
      2 2 2 2
    exit
  compute overlap
  set view:0 0.16242 -0.485669 0 11.0175 1 0 0 0 1 0 0 0 1
  reset:0
  exit
save a grid
intersect2.hdf
twocell
exit
