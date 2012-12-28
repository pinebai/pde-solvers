///
/// Brief description: Loop macros for CellWave::GenericChemistry
///

/// For loop over a 3D grid. Assumes defined:
///		     const int&nd, const int&ncomp,
///		     const int &nd1a, const int &nd1b,
///		     const int &nd2a, const int &nd2b,
///		     const int &nd3a, const int &nd3b,
///		     const int &n1a,  const int &n1b,
///		     const int &n2a,  const int &n2b,
///		     const int &n3a,  const int &n3b
#define ForGrid(i,j,k) for (int k=nd3a; k<=nd3b; ++k)    \
                         for (int j=nd2a; j<=nd2b; ++j)  \
                           for (int i=nd1a; i<=nd1b; ++i)

//const int jlen=1+nd1b-nd1a;
//const int klen=jlen*(1+nd2b-nd2a);
//const int dlen=klen*(1+nd3b-nd3a);
//#define INDEX(i,j,k,ic)  ( (i-nd1a)+jlen*(j-nd2a)+klen*(k-nd3a)+dlen*ic )

#define GRIDINDEX(i,j,k,ic) ( (i-nd1a)+(1+nd1b-nd1a)*(j-nd2a)+    \
                            (1+nd1b-nd1a)*(1+nd2b-nd2a)*(k-nd3a)  \
                            + (1+nd1b-nd1a)*(1+nd2b-nd2a)*(1+nd3b-nd3a)*ic )
