#include <string>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  std::string name="test parameter";
  double param=7.19;
  double def=13.2;

  char buf[80], bufdef[80];
  sprintf(buf, "%16.8e", param );
  sprintf(bufdef, "%16.8e", def   );
  std::string outpar(buf);
  std::string outdef(bufdef);

  printf("  %16s <%16s> = %16s \n", 
	 name.c_str(), outdef.c_str(), outpar.c_str());
    
}
