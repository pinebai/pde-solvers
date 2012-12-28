#include "DPrintf.h"

int main(int argc, char **argv)
{
  CellWave::DPrintfCreate();

  /*CellWave::DPrintfStop(1);*/
  CellWave::DPrintfOpen(1, "debug.out");
  
  CellWave::DPrintf(0, "Hello world\n");
  CellWave::DPrintf(1, "debugging to stream %d", 1);
  CellWave::DPrintf(1, "... this should be good.\n");

  CellWave::DPrintfDestruct();
}
