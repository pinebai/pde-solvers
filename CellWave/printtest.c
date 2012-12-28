/*
 *  test indirect printing 
 *
 */

#include <stdio.h>
#include <stdarg.h>

void dprintf(char *fmt, ...)
{
  va_list ap; /* pointer to each unnamed arg in turn */
  char *p, *aval;

  const int nlenCBuf = 65535;
  char cBuf[nlenCBuf];
  
  va_start( ap, fmt);
  vsprintf( cBuf, fmt, ap );
  
  printf("-----------Processed print, output follows:----------\n");
  printf("%s", cBuf );
  printf("-----------done--------------------------------------\n");
  
}

int main(int argc, char **argv)
{
  dprintf("Hello there %i, this is <%8.4g> what do you say?\n",7, 3.14);
}
