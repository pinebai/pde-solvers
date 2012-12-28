*
* simple 2D grid for 1mm diameter X. Laevis oocyte 
*   LENGTHS ARE IN micro m=10^-6 meters.
*
*  min spacing: 
*    annulus, dtheta= 19, dr=17.8
*    square,  dx= 20
*
*  min dt, mu=beta*D=16 um^2/s
*  --> 
create mappings
*
rectangle
  set corners
    -500. 500. -500. 500.
**    -2. 2. -2. 2.
  lines
    50 50 
  boundary conditions
   0 0 0 0 
  mappingName
  square
exit
*
Annulus
  lines
    80 14
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
oocyte.hdf
oocyte
exit

*oocyte_1000_um.hdf

