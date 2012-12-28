// -*- c++ -*-
#ifndef CELLWAVE_PARAMETER_READER_H
#define CELLWAVE_PARAMETER_READER_H "CellWave/ParameterReader.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>

#include <map>
#include <boost/any.hpp>
#include <boost/tokenizer.hpp>

//#define DEBUGOUT(x) {x}
#define DEBUGOUT(x) {}

namespace CellWave {

  class ParameterReader {
  public:
    const int maxBufferLength;
    const int maxParameterLength;
    std::string outputFileName;
    std::string parameterFileName;
    typedef std::map<std::string, boost::any> ParameterList;
    enum WarningCodes {WarningOK=0, WarningWrongType=-1, 
		       WarningSetFromDefault=1  };

    ParameterReader( int &argc, char ** &argv) 
      : maxBufferLength(80), outputFileName(""), maxParameterLength(40) 
    {
    }; 
    ParameterReader( std::string fname ) 
      : maxBufferLength(80), outputFileName(""), maxParameterLength(40) 
    {
      parameterFileName = fname;
      readFile( fname );
    };  

    ~ParameterReader() {}; // standard

    ParameterList q;

    ParameterList & getParameterList() {return this->q;};

    void setOutputFile( std::string outFile_) {
      outputFileName = outFile_;
      std::FILE *fp=std::fopen(outputFileName.c_str(),"w");
      std::fclose(fp); //zaps the file
    };

    std::string getParameterFileName()
    {
      return( parameterFileName );
    }
    
    bool findKey( const std::string &paramName, boost::any &z )
    {
      ParameterReader::ParameterList::iterator pos;
      pos = this->q.find(paramName);
      if (pos != q.end()) { // found
	z = pos->second;
	return true;
      }
      else return false;
    };
  
    bool
    get( const std::string &paramName,
		    int &iValue, int defaultValue=0)
    {
      WarningCodes error;
      return( this->get(paramName, iValue, defaultValue, error));
    };

    bool
    get( const std::string &paramName,
		    int &iValue, int defaultValue, WarningCodes &error)
    {
      DEBUGOUT(std::cout <<"get(int)\n";);      
      char buf[maxBufferLength];
      sprintf(buf,"%10i",defaultValue);
      std::string s, temp(buf);
      
      bool useDefault= get( paramName, s, temp, error);
      if (useDefault) {
	iValue=defaultValue;
      } else {
	std::sscanf( s.c_str(), "%i", &iValue);
      }
      return useDefault;
    };

    bool
    get(  const std::string &paramName,
		    bool &bValue, bool defaultValue=false)
    {
      WarningCodes error;
      return( this->get(paramName, bValue, defaultValue,error));
    }

    bool
    get(  const std::string &paramName,
	  bool &bValue, bool defaultValue, WarningCodes &error)
    {
      int temp, intDefault=defaultValue;
      get( paramName, temp, intDefault);
      bValue=temp;
    }

    bool
    get( const std::string &paramName,
		    double &dValue, double defaultValue=0)
    {
      WarningCodes error;
      return( this->get(paramName, dValue, defaultValue, error));
    }
  
    double
    get( const std::string &paramName,
		    double &dValue, double defaultValue, WarningCodes &error)
    {
      DEBUGOUT(std::cout <<"get(double)\n";);
      char buf[maxBufferLength];
      sprintf(buf,"%24.16g",defaultValue);
      std::string s, temp(buf);

      bool useDefault= get( paramName, s, temp, error);
      if (useDefault) {
	dValue=defaultValue;
      } else {
	std::sscanf( s.c_str(), "%lg", &dValue);
      }
      return useDefault;
    };   
    
    bool
    get( const std::string &paramName,
		    std::string &sValue, std::string defaultValue="")
    {
      WarningCodes error;
      return( this->get(paramName, sValue, defaultValue, error));
    };

    std::string 
    get( const std::string &paramName )
    {
      WarningCodes error;
      std::string temp, defaultValue="";
      this->get( paramName, temp, defaultValue);
      return (temp);
    }

    bool
    get( const char *paramName,
		    std::string &sValue, std::string defaultValue="")
    {
      std::string temp = paramName;
      WarningCodes error;
      return( this->get(paramName, sValue, defaultValue, error));
    };
 
    bool
    get( const std::string &paramName,
		    std::string &sValue, std::string defaultValue, 
		    WarningCodes &error)
    {
      DEBUGOUT(std::cout <<"get(std::string)\n";);
      boost::any z;
      sValue = defaultValue;
      bool useDefault=true;
      error = WarningOK;
      if (findKey( paramName, z)) {
	DEBUGOUT(std::cout << "..found "<<paramName<<"\n"; );
	try {
	  sValue = boost::any_cast<std::string>( z );
	  useDefault=false;
	  DEBUGOUT(std::cout <<"..set the value as "<<sValue<<"\n";);
	}
	catch ( ... ) {
	  error= WarningWrongType;
	  DEBUGOUT( std::cerr << "Warning:: parameter "
		    << paramName <<" has wrong type\n";);
	}
      }
      else error=WarningSetFromDefault;
      if (useDefault) {
	DEBUGOUT(std::cerr << "Warning:: parameter "
		 <<paramName<<" set from defaults\n";);
      }
      
      DEBUGOUT(std::cout << "GOT: "<<paramName<<"== "<< sValue<<std::endl;);
      if (outputFileName != "") {
	std::FILE *fp=std::fopen(outputFileName.c_str(),"a");
	char fmt[maxBufferLength];
	sprintf(fmt,"%%%is = %%s\n",maxParameterLength);
	//fprintf(fp,"%s = %s\n",paramName.c_str(), sValue.c_str());
	fprintf(fp,fmt,paramName.c_str(), sValue.c_str());
	std::fclose(fp);
      }
      return useDefault;
    }    

    int oldReadFile(std::string fileName)
    {
      using namespace std;
      ParameterList & q = this->getParameterList();

      int iLines=0;
      FILE *fp;
      fp = fopen(fileName.c_str(), "r");
      if (fp == NULL) return -1;
      const int maxLength = 1024;
      char *buf = new char[maxLength];
      string param;
      string value;
      for ( ;; ) {
	char* flag= fgets( buf, maxLength, fp);
	if (flag==NULL) break;
	string input(buf);
	param="";
	value="";
	
	parseALine( input, param, value);
	DEBUGOUT(std::cout << "adding key["<<param<<"]\n";)
	q[param] = value;
      }
    }

    int readFile(std::string fileName)
    {
      using namespace std;
      ParameterList & q = this->getParameterList();

      int iLines=0;      
      //FILE *fp;
      //fp = fopen(fileName.c_str(), "r");
      ifstream fs(fileName.c_str());
      if (!fs) return -1;
      string input, param, value;
      while( getline(fs, input )) {
	param=""; value="";
	
	parseALine( input, param, value);
	if ( param == "" ) {
	  DEBUGOUT(cout << "empty key["<<param<<"]\n";);
	} else {
	  DEBUGOUT(std::cout << "adding key["<<param<<"]\n";);
	  q[param] = value;
	}
      }
    }

    void
    eatSurroundingWhiteSpace(std::string &inout)
    {
      int last=-1;
      //std::cout << "BEFORE: ["<<inout.substr(0,inout.length())<<"]";
      //std::cout << ", length="<<inout.length()<<std::endl;
      for (last=inout.length()-1; last>=0; --last ) {
	int inv=inout[last];
	//cout << "  <i="<<last<<":"<< inout[last]<<":"<<inv<<">   ";
	if ( !isspace( inout[last])) break;
      }
      //cout <<"\n-- finally last="<<last<<", <"
	//   <<inout.substr(0,last+1)<<">"<<endl;
      if (last>=0) inout = inout.substr(0,last+1);

      int first;
      for (first=0; first<=last-1; ++first) {
	if( !isspace( inout[first])) break;
      }
      inout = inout.substr(first); // cut [0,first-1]
      //std::cout << "AFTER: ["<<inout<<"]\n";

    }
  
    bool
    parseALine(std::string input, std::string &param, std::string &value)
    {
      
      std::string str=input;
      param="";
      int idx = input.find("#");
      if (idx != std::string::npos) {
	if (idx>0) str = input.substr(0,idx-1);
	else str="";
      }
      typedef boost::tokenizer<boost::char_delimiters_separator<char> > Tok;
      const char dropped[] = "=";   
      boost::char_delimiters_separator<char> sep(false,"",dropped);
      Tok tok(str, sep);
      
      Tok::iterator tok_iter=tok.begin();
      if (tok_iter != tok.end()) {
	param = *tok_iter;
	if ( ++tok_iter!= tok.end()) {
	  value = *tok_iter;
	  //int len=value.length();
	  eatSurroundingWhiteSpace( value );
	};
      }
      
      //if (param!="") {
      //boost::char_delimiters_separator<char> sep2;
      //Tok nowhitespace(param, sep2);
      //Tok::iterator tok_iter2 = nowhitespace.begin();
      //param =*tok_iter2;
      //      }
      eatSurroundingWhiteSpace( param );
      if (param=="") value="";
      DEBUGOUT(std::cout <<param<< "= ["<<value<<"]"<<  std::endl;)
      return( param!="" );
    }
    
  }; //end class ParameterReader
}; //end namespace CellWave

#endif //CELLWAVE_PARAMETER_READER_H
