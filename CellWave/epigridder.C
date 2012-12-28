//
// EpiGridder
//

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "epigridder.h"

void inputVertex(FILE *fp, Vertex &vert)
{
  const int BUFSIZE=160;
  const int MAXTOK=10;
  char buf[BUFSIZE];
  char *pbuf;

  assert( fp != NULL );
  fgets( buf, BUFSIZE, fp );
  //sscanf(buf, "%d %g; %gl", &vert.id, &vert.x, &vert.y);

  sscanf( strtok( buf, " "),   "%d",   &vert.id);
  sscanf( strtok( NULL, " "),  "%lf",  &vert.x);
  sscanf( strtok( NULL, " "),  "%lf",  &vert.y);
  
  //fscanf(fp, "%d %g; %gl", &vert.id, &vert.x, &vert.y);
  printf(" > vertex %d,  x= %g,  y= %g \n", 
	 vert.id, vert.x, vert.y);
}

void inputEdge(FILE *fp, Edge &edge)
{
  const int BUFSIZE=160;
  const int MAXTOK=10;
  char buf[BUFSIZE];
  char *pbuf;

  assert( fp != NULL );
  fgets( buf, BUFSIZE, fp );

  sscanf( strtok( buf, " "),   "%d",   &edge.id);
  sscanf( strtok( NULL, " "),  "%d",   &edge.npoints);
  printf(" > edge %d: vertices= ", edge.id);
  for (int i=0; i< edge.npoints; ++i ) {
    sscanf( strtok( NULL, " "),  "%d",  &edge.vertID[i]);
    printf("%d ", edge.vertID[i]);
  }
  printf("\n");
  //fscanf(fp, "%d %g; %gl", &vert.id, &vert.x, &vert.y);
}

void inputTFICell(FILE *fp, int nxpoints, int nypoints, TFICell &cell)
{
  const int BUFSIZE=160;
  const int MAXTOK=10;
  char buf[BUFSIZE];
  char *pbuf;

  cell.nxpoints = nxpoints;
  cell.nypoints = nypoints;

  assert( fp != NULL );
  fgets( buf, BUFSIZE, fp );

  sscanf( strtok( buf, " "),   "%d",   &cell.id);
  //sscanf( strtok( NULL, " "),  "%d",   &cell.npoints);
  printf(" > TFI cell %d: edges= ", cell.id, 2);
  sscanf( strtok( NULL, " "),  "%d",  &cell.edgeLeft);
  sscanf( strtok( NULL, " "),  "%d",  &cell.edgeRight);
  printf("%d ", cell.edgeLeft);
  printf("%d ", cell.edgeRight);
  printf("\n");
  //fscanf(fp, "%d %g; %gl", &vert.id, &vert.x, &vert.y);
}

void outputEdge(FILE *fp, const Vertex *verts, const Edge &edge)
{
  assert( fp != NULL);
  
  fprintf(fp, "***   EDGE %d\n", edge.id );
  fprintf(fp, " spline\n");
  fprintf(fp, "   enter spline points\n");
  fprintf(fp, "   %d\n", edge.npoints);
  for( int i=0; i< edge.npoints; ++i ) {
    int vid=edge.vertID[i]-1;
    fprintf(fp, "   %8.3e %8.3e\n", verts[vid].x, verts[vid].y);
  }
  fprintf(fp, "   shape preserving (toggle)\n");
  fprintf(fp, "   mappingName\n");
  fprintf(fp, "   edge%d\n", edge.id);
  fprintf(fp, "   exit\n");
  fprintf(fp, "*\n");

};

void outputTFICell(FILE *fp, const Vertex *verts, 
		   const Edge *edges,
		   const TFICell &cell)
{
  assert( fp != NULL );

  fprintf(fp, "*** CELL %d\n", cell.id);
  fprintf(fp, " tfi\n");
  fprintf(fp, "   choose left curve   (r_1=0)\n");
  fprintf(fp, "     edge%d\n",  cell.edgeLeft);
  fprintf(fp, "   choose right curve  (r_1=1)\n");
  fprintf(fp, "     edge%d\n",  cell.edgeRight);
  fprintf(fp, "   lines\n");
  fprintf(fp, "   %d %d\n", cell.nxpoints,  cell.nypoints);
  fprintf(fp, "   mappingName\n");
  fprintf(fp, "   cell%d\n", cell.id );	  
  fprintf(fp, "   exit\n");
  fprintf(fp, "   *\n");
}

