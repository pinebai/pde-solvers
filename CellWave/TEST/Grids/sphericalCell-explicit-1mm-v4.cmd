*
* command file to create a sphere in a box
*
*  time to make: 594s new: 3.5
*     cpu=2s (ov15 sun-ultra optimized)
*
create mappings
* first make a sphere
Sphere
    inner radius
    250
    outer radius
    500
exit
*
* now make a mapping for the north pole
*
reparameterize
  orthographic
    specify sa,sb
      2.5 2.5
  exit
  lines
     200 200 60
*     150 150 45
*    100 100 30
*    50 50 15
*    15 15 5
*  boundary conditions
*    0 0 0 0 1 0
*  share
*    0 0 0 0 1 0
  boundary conditions
    0 0 0 0 0 1
  share
    0 0 0 0 0 1
  mappingName
    north-pole
exit
*
* now make a mapping for the south pole
*
reparameterize
  orthographic
    choose north or south pole
      -1
    specify sa,sb
      2.5 2.5
  exit
  lines
    200 200 60
*   150 150 45
*   100 100 30
*   50 50 15
*    15 15 5
*  boundary conditions
*    0 0 0 0 1 0
*  share
*    0 0 0 0 1 0
  boundary conditions
    0 0 0 0 0 1
  share
    0 0 0 0 0 1
  mappingName
    south-pole
exit
*
* Here is the box
*
Box
  set corners
     -300 300 -300 300 -300 300
*    -.6 .6 -.6 .6 -.6 .6
*    -2 2 -2 2 -2 2 
  lines
   120 120 120
*   90 90 90
*   60 60 60
*   30 30 30
*    21 21 21
  boundary conditions
    0 0 0 0 0 0
  mappingName
    box
  exit
exit
*
generate an overlapping grid
  box
  north-pole
  south-pole
  done
  change parameters
    interpolation type
      explicit for all grids
    ghost points
      all
      2 2 2 2 2 2
  exit
  compute overlap
exit
save an overlapping grid
sphericalCell-explicit-v4.hdf
sphericalCell
exit
