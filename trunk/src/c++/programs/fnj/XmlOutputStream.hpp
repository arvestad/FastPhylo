#ifndef XMLOUTPUTSTREAM_HPP
#define XMLOUTPUTSTREAM_HPP

#include <cstdio>
#include "DataOutputStream.hpp"



class XmlOutputStream : public DataOutputStream
{
public:
  XmlOutputStream();
  XmlOutputStream(char * filename );
  virtual ~XmlOutputStream();
protected:
  virtual void print( tree2int_map & tree2count, bool noCounts, std::vector<std::string> & names, Extrainfos & extrainfos );

 /** printNewick unparses the xml string created before. Of course it would better if the SequenceTree would produce both a newick string
     and some xml representation of newick. 
**/
  virtual void printNewick(std::ostream * fp , std::string s );
  virtual void printNewickNode(std::ostream * fp , xmlNode * node);
};

#endif // XMLOUTPUTSTREAM_HPP
