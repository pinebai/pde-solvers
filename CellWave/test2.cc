#include <stdio.h>

int main( int argc, char **argv)
{
  int side, axis, i=0;
  for( axis=0; axis<3;  ++axis) {
    for( side=0; side<=1; ++side) {
      i++;
      printf(" i= %2d,  axis= %d, side= %d\n", i, axis, side);
      if( i==3 ) goto done;
    }
  }
 done:
  printf("********done -- side= %d, axis= %d, i=%d \n", side, axis, i );
  return 0;
}
