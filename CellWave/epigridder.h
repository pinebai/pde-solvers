//
// EpiGridder -- build a multigrid epithelial grid 
//
#

const int MAXEDGEPOINTS=10;

struct Vertex {
  int    id;
  double x;
  double y;
};

struct Edge {
  int id;
  int npoints;
  //double x[MAXEDGEPOINTS];
  //double y[MAXEDGEPOINTS];
  int vertID[MAXEDGEPOINTS];
};

struct TFICell {
  int id;
  int edgeLeft;
  int edgeRight;
  int nxpoints;
  int nypoints;
};


void inputVertex(FILE *fp, Vertex &vert);
void inputEdge(FILE *fp, Edge &edge);
void inputTFICell(FILE *fp, int nxpoints, int nypoints, TFICell &cell);

void outputEdge(FILE *fp,const Vertex *verts, const Edge &edge);
void outputTFICell(FILE *fp,  const Vertex *verts, 
		   const Edge *edges, const TFICell &cell);
