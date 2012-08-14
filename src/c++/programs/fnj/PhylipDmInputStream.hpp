/*
 * PhylipDmInputStream.hpp
 *      Auther: Mehmood Alam Khan Email: malagori@kth.se
 */
#ifndef PHYLIPDMINPUTSTREAM_HPP
#define PHYLIPDMINPUTSTREAM_HPP

#include <cstdio>
#include "DataInputStream.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


class PhylipDmInputStream : public DataInputStream
{
public:
  PhylipDmInputStream(char * filename );
  virtual ~PhylipDmInputStream();
  virtual readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos );
  virtual readstatus readFloatDM( StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos )  {};
protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;

 
};


#endif // PHYLIPDMINPUTSTREAM_HPP
