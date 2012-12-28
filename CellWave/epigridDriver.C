#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "epigridder.h"

// EpiGrid file format .epig
// #<vertices>
//    nvertices
//    1 x y
//    2 x y
//    ...
// #<edges>
//    nedges
//    1 nverts vert1 vert2 ... vertn
//    ...
// #<cells>
//    ncells
//    1 nedges edge1 ... edgen
// #<tficells>
//    nTFICells 
//    1 leftedge rightedge
// 
int main(int argc, char **argv)
{
  Edge    edge;
  TFICell cell;
  FILE *fpIn, *fpOut;
  const int BUFSIZE=160;
  char buf[BUFSIZE];
  char nametag[BUFSIZE];

#define MAXVERTICES  3000
#define MAXEDGES     1000
#define MAXTFICELLS  300 

  const int nxpoints=20, nypoints=20;
  //#define MAXCELLS     300

  Vertex   vertices[ MAXVERTICES ];
  Edge     edges[ MAXEDGES ];
  TFICell  tficells[ MAXTFICELLS ];
  //Cell     cells[ MAXCELLS ];

  int      numberOfVertices;
  int      numberOfEdges;
  int      numberOfTFICells;
  int      numberOfCells;

  if (argc <= 2 ) {
    printf("usage: %s <cell description.epig> <overture output file.cmd>\n", argv[0]);
    exit(-1);
  }

  fpIn = fopen( argv[1], "r");
  if (fpIn == NULL ) {
    printf("Couldn't open %s for reading\n", argv[1]);
  }

  fgets( buf,     BUFSIZE, fpIn ); //comment header
  fgets( nametag, BUFSIZE, fpIn ); //<vertices>
  fgets( buf,     BUFSIZE, fpIn );
  sscanf( buf, "%d", &numberOfVertices);

  printf("* numberOfVertices= %d  %s",  numberOfVertices, nametag);
  for( int i=0; i< numberOfVertices; ++i ) {
    inputVertex( fpIn, vertices[ i ]);
  }

  fgets( nametag, BUFSIZE, fpIn ); //<edges>
  fgets( buf,     BUFSIZE, fpIn );
  sscanf( buf, "%d", &numberOfEdges);

  printf("* numberOfEdges= %d  %s",  numberOfEdges, nametag);
  for( int i=0; i< numberOfEdges; ++i ) {
    inputEdge( fpIn, edges[ i ]);
  }

  fgets( nametag, BUFSIZE, fpIn ); //<cells>
  fgets( buf,     BUFSIZE, fpIn );
  sscanf( buf, "%d", &numberOfCells);

  printf("* numberOfCells= %d  %s", numberOfCells, nametag);
  for( int i=0; i< numberOfCells; ++i ) {
    fgets( buf, BUFSIZE, fpIn ); // skip the cell definitions
  }

  fgets( nametag, BUFSIZE, fpIn ); //<tficells>
  fgets( buf,     BUFSIZE, fpIn );
  sscanf( buf, "%d", &numberOfTFICells);

  printf("* numberOfTFICells= %d  %s",  numberOfTFICells, nametag);
  for( int i=0; i< numberOfTFICells; ++i ) {
    inputTFICell( fpIn, nxpoints, nypoints, tficells[ i ]);
  }
  printf("--done reading--\n");
  fflush(0);
  fclose( fpIn );

  fpOut = fopen(argv[2], "w+");
  
  fprintf(fpOut,"** test epigridder **\n");
  fprintf(fpOut,"create mappings\n");

  for(int i=0; i< numberOfEdges; ++i) {
    outputEdge( fpOut, vertices, edges[i]);
  }

  for(int i=0; i< numberOfTFICells; ++i) {
    outputTFICell( fpOut, vertices, edges, tficells[i]);
  }

  fflush(0);
  fclose( fpOut );
}
