#include <stdio.h>
#include <iostream>
#include <map>
#include <string>

int main(int argc, char **argv)
{
  typedef   std::multimap<int,std::string> IntStringMap;
  typedef   std::multimap<int,std::string>::iterator Iter;

  IntStringMap mm;

  mm.insert(std::make_pair(1,"hello"));
  mm.insert(std::make_pair(1,"world"));

  Iter pos;
  const int key=1;
  std::cout <<">>: ";
  for( pos=mm.lower_bound(key); 
       pos!=mm.upper_bound(key);
       ++pos) {
    std::cout << pos->second << " ";
  }
  std::cout << std::endl;

}
