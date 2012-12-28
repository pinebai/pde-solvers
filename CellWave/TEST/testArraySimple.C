#include "Overture.h"
#include <iostream>
#include <stdio.h>

#include "ArraySimple.h"

const char* oktext( bool isok)
{
  if( isok ) return "ok";
  else       return "not ok";
}

int main(int argc, char **argv)
{
  Overture::start( argc, argv );

  printf("---test 1 ---\n");
  ArraySimple<int> start(10,3,2);
  ArraySimple<int> stop(10,3,2);

  for (int i=0;i<10;++i ) {
    for( int j=0; j<3; ++j ) {
      for( int k=0; k<2; ++k ) {
	start(i,j,k) = i*j*k;
	stop(i,j,k) =  i*j*k+1;
      }
    }
  }

  for (int i=0;i<10;++i ) {
    for( int j=0; j<3; ++j ) {
      for( int k=0; k<2; ++k ) {
	printf(" i,j,k=(%d,%d,%d), start=%d, stop=%d.\n",
	       i,j,k,start(i,j,k), stop(i,j,k));
      }
    }
  }

  printf("---test 2 ---\n");
  const int n1=100, n2=30, n3=20;
#if 0  //--this version is ok
  ArraySimple<int> *pTest;
  pTest = new ArraySimple<int>(n1,n2,n3);
  ArraySimple<int> &test = *pTest;
#endif

  //--this version is also ok
  ArraySimple<int> test;
  test.redim(n1,n2,n3);
  for (int i=0;i<n1;++i ) {
    for( int j=0; j<n2; ++j ) {
      for( int k=0; k<n3; ++k ) {
	test(i,j,k) = i*j*k;
      }
    }
  }  

  for (int i=0;i<n1;++i ) {
    for( int j=0; j<n2; ++j ) {
      for( int k=0; k<n3; ++k ) {
	int it_is=test(i,j,k);
        int it_should_be= i*j*k;
	bool isok= ( it_is == it_should_be );
	printf(" test= %8d  should be = %8d  (%s)\n",test(i,j,k),i*j*k,oktext(isok) );
      }
    }
  }  

}
