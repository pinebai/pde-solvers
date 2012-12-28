#include <stdio.h>
#include <time.h>

int main(int argc, char **argv)
{
  time_t date;

  //printf("short date= %s\n", asctime(localtime(time(NULL))));
  time(&date );

  printf("dateOfSimulation=%s\n",asctime(localtime( &date) ));
  

}
