//
// testing array subsets
//
#include "Overture.h"

int main(int argc, char **argv)
{
  Overture::start( argc, argv);
  IntegerArray aa(10);
  aa=0;
  Range sub(3,6);
  IntegerArray bb(sub);
  bb=1;

  aa(sub) = bb(sub);

  aa.display("aa");
  Overture::finish();
}
