#include <iostream>
#include <stdio.h>
#include <vector>
#include <boost/tokenizer.hpp>


int main(int argc, char **argv)
{
  std::string cornerString="   0.1     0.1      0.1   20 20 10";
  
  typedef boost::tokenizer<boost::char_separator<char> >            Tokenizer;
  typedef std::vector<double>                                       DoubleVector;
  typedef boost::tokenizer<boost::char_separator<char> >::iterator  TokIterator;

  boost::char_separator<char> sep(", ");
  Tokenizer tok(cornerString, sep);
  DoubleVector corners;
  corners.clear();

  for (TokIterator it=tok.begin(); it!=tok.end(); ++it ) {
    double q;
    printf("");
    sscanf(it->c_str(), "%lf", &q);
    //DPrintf(DebugReaction," %lf  ['%s'];  ", q, it->c_str());
    corners.push_back( q );
  }
  printf("the corners are:\n   ");
  printf("%e", corners[0]);
  for(int i=1; i<corners.size(); ++i ) {
    printf(", %e", corners[i]);
  }
  printf("\n");
}
