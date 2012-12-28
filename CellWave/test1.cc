#include <string>
#include <stdio.h>

int main( int argc, char **argv) 
{
  std::string za="filename";
  za += ".m";

  printf(" new name='%s'\n", za.c_str());

  return 1;
}
