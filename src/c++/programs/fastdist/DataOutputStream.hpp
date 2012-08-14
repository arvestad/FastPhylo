#ifndef DATAOUTPUTSTREAM_HPP
#define DATAOUTPUTSTREAM_HPP

#include <cstdio>
#include <vector>
#include <string>
#include "Exception.hpp"
#include "DistanceMatrix.hpp"
#include "DistanceRow.hpp" //Mehmood's Changes here. Email: malagori@kth.se

//#include "FloatDistanceMatrix.hpp"
#include "Extrainfos.hpp"


class DataOutputStream
{
public:
  DataOutputStream( );
  DataOutputStream(char * filename );
  virtual ~DataOutputStream() {};
  virtual void print( StrDblMatrix & dm )=0;  //Mehmood's Changes here. Email: malagori@kth.se
  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos )=0;
  virtual void printEndRun()=0;
  virtual void printRow( StrFloRow & dm , std::string name, int row)=0 ;
  virtual void printHeader( size_t numNodes )=0 ;
  virtual void printBootstrapSpliter(size_t numNodes)=0 ;


  //Mehmood's addition here
  // FAST PRINTING OF FLOATS
  static const char ONEDIGIT[128];
  static const char TENDIGIT[128];

protected:
  FILE * fp;
  bool file_was_opened;

};

//void setXmlFlag(bool xmlFlag);
void printPHYLIPfast(const StrDblMatrix &dm, FILE *out, bool flag);

/*class PhylipDmOutputStream : public DataOutputStream
{
public:
  PhylipDmOutputStream(char * filename ) : DataOutputStream(filename) {};
  virtual ~PhylipDmOutputStream() {};
  virtual void print( StrDblMatrix & dm );
};
*/
#endif // DATAOUTPUTSTREAM_HPP
