/// Brief description: debug print facility from CellWave
///
///

#ifndef DPRINTF_H
#define DPRINTF_H

namespace CellWave {
  void DPrintfCreate();
  void DPrintfDestruct();

  //  void DPrintfOpenSwill( int portnumber );
  //  void DPrintfCloseSwill();

  void DPrintfSetDefault( int def_ );
  int  DPrintfOpen( int tag, char *fname );
  void DPrintfClose( int tag );
  void DPrintfStart( int dlevel );
  void DPrintfStop( int dlevel );
  void DPrintfStopDebug();
  int  DPrintfGetStatus( int dlevel );
  void DPrintf( int dlevel, char *fmt, ...);

  void DPrintfOutputFileHandle(int tag);
  int  DPrintfGetFileHandle(int tag);

}; 

#endif
