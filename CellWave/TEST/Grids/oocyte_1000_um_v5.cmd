*
* simple 2D grid for 1mm diameter X. Laevis oocyte 
*   LENGTHS ARE IN micro m=10^-6 meters.
*
*  min spacing: 
*    annulus, dtheta= 4.2, dr=4.2
*    square,  dx= 5 
*
create mappings
*
rectangle
  set corners
    -500. 500. -500. 500.
**    -2. 2. -2. 2.
  lines
      400 400
*      200 200
*     100 100
*    50 50 
  boundary conditions
   0 0 0 0 
  mappingName
  square
exit
*
Annulus
  lines
     640 112
*    320 56
*    160 28
*    80 14
*  centre
*    0. 1.
    inner and outer radii
    250.  500.
  boundary conditions
    -1 -1 0 1 
exit
*
exit
generate an overlapping grid
    square
    Annulus
  done
  change parameters
    * choose implicit or explicit interpolation
    * interpolation type
    *   implicit for all grids
    ghost points
      all
      2 2 2 2 2 2
  exit
*  display intermediate results
  compute overlap
  exit
*
save an overlapping grid
oocyte5.hdf
oocyte
exit

*oocyte_1000_um.hdf

