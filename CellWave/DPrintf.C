/// Brief description: debug print facility (from CellWave)
///

#include <stdio.h>
#include <stdarg.h>

namespace CellWave {
# define MAX_DPRINTF_LEVELS 20
# define MAX_DPRINTF_LENGTH 65535
  static int DPrintf_default;
  static int DPrintf_outputflag[MAX_DPRINTF_LEVELS];
  static FILE* DPrintf_file[MAX_DPRINTF_LEVELS];

  void DPrintfCreate() {
    int i;
    for ( i=0; i< MAX_DPRINTF_LEVELS; ++i) {
      DPrintf_outputflag[i] = 1; /* enabled */
      DPrintf_file[i] = stdout;
    }
    DPrintf_default = 0;
  }

  void DPrintfDestruct() {
    int i;
    for ( i=0; i< MAX_DPRINTF_LEVELS; ++i) {
      FILE*  fp=DPrintf_file[i];
      if( fp != stdout ) fclose( fp );
    }
    DPrintf_default = -1;
  }

  void DPrintfOutputFileHandle(int tag) {
    printf("  DPrintf::file_handle[%d] = %d\n",
	   tag, int(DPrintf_file[tag]));
    fflush(0);
  }

  int  DPrintfGetFileHandle(int tag) {
    return( int( DPrintf_file[tag]));
  }

  void DPrintfSetDefault( int def_ ) {
    DPrintf_default = def_;
  };

  int DPrintfOpen( int tag, char *fname )  {
    int ierr;
    FILE *fp;
    
    fp = fopen( fname, "w+");
    ierr = (fp!=0);
    if (ierr != 0 ) DPrintf_file[tag] = fp;
    return( ierr );
  }

  void DPrintfClose( int dlevel ) {
    if ( ( 0<= dlevel ) && (dlevel < MAX_DPRINTF_LEVELS) ) {
      FILE *fp=DPrintf_file[ dlevel ];
      if (fp != stdout ) {
	fclose(fp );
	DPrintf_file[ dlevel ] = stdout;
      }
    }
  }

  void DPrintfStart( int dlevel ) {
    if ( ( 0<= dlevel ) && (dlevel < MAX_DPRINTF_LEVELS) ) {
      DPrintf_outputflag[ dlevel ] = 1;
    }
  }

  void DPrintfStop( int dlevel ) {
    if ( ( 0<= dlevel ) && (dlevel < MAX_DPRINTF_LEVELS) ) {
      DPrintf_outputflag[ dlevel ] = 0;
    }
  }

  void DPrintfStopDebug() {
    int i;
    for ( i=2; i<MAX_DPRINTF_LEVELS; ++i) {
      DPrintfStop( i);
    }
  }

  int  DPrintfGetStatus( int dlevel ) {
    if ( ( 0<= dlevel ) && (dlevel < MAX_DPRINTF_LEVELS) ){
      return( DPrintf_outputflag[ dlevel ] );
    }
  }  

  void DPrintf( int dlevel, char *fmt, ...)  {
    va_list ap; /* pointer to each unnamed arg in turn */

    //fprintf( DPrintf_file[dlevel], "@");

    if ( 0 == DPrintf_outputflag[ dlevel ] ) return; /* no output */
    
    const int nlenCBuf = MAX_DPRINTF_LENGTH;
    static char cBuf[nlenCBuf];
    
    va_start( ap, fmt);
    vsprintf( cBuf, fmt, ap );
    va_end( ap );

    if ( ( 0<= dlevel ) && (dlevel < MAX_DPRINTF_LEVELS) ) {
      //fprintf( DPrintf_file[dlevel],"DPRINT       (%d):",dlevel); //XXX DEBUG CODE
      fprintf( DPrintf_file[dlevel],"%s", cBuf );
    }
    if (dlevel == 0) {
      int i;
      for ( i =1; i<MAX_DPRINTF_LEVELS; ++i ) {
	if (DPrintf_file[i] !=stdout ) {
	  if ( 0 != DPrintf_outputflag[ i ] ) {
	    //fprintf( DPrintf_file[i],"BROADCAST(%d):",dlevel); //XXX DEBUG CODE
	    fprintf( DPrintf_file[i], "%s", cBuf); /* broadcast to files*/
	  }
	}
      } // end for i
    }// end if dlevel==0
    fflush(0);
  };

}; // end namespace FlowToolkit

