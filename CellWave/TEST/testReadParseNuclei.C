//
// read & parse nucleus data files (.cwn)
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>

int main(int argc, char **argv)
{
  const std::string cn_filename="Grids/testNuclei.cwn";
  bool okFlag=false;

  FILE *fp=fopen( cn_filename.c_str(),"r");
  if( !fp ) {
    okFlag=false;
    printf("error: CW nucleus file '%s' not found, exiting.\n",
	   cn_filename.c_str());
    return( int(okFlag) );
  }
  printf("reading file '%s'...\n",cn_filename.c_str());
  printf("------------------------------\n");

  const int bufferLength=1000;
  char buffer[bufferLength];
  
  while( fgets( buffer, bufferLength, fp)) {
    const int lineLength=strlen(buffer);

    typedef std::vector<std::string> StringVector;
    StringVector tokens;

    if( lineLength>0 ) { 
      if( buffer[0]=='#' ) {
	buffer[lineLength-1]=0;
	printf("comment< %s >\n", buffer);
      }
      else {
	using namespace std;
	buffer[lineLength-1]=0;
	//printf("regular< %s >\n", buffer);
	string line(buffer);
	//cout << "<"<<line<<">\n";

	typedef boost::tokenizer<>::iterator TokIterator;
	boost::tokenizer<> tok(line);
	for (TokIterator it=tok.begin(); it!=tok.end(); ++it) {
	  printf(" <%s> ", it->c_str());
	  tokens.push_back( *it );
	}
	printf("   number of tokens=%d\n", tokens.size());
	
      }
    };

  }
  printf("-------------done-------------\n");
  fclose(fp);
  okFlag=true;

  return( int(okFlag ));

}
