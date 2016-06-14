#include "FastaInputStream.hpp"
#include <cstdio>


using namespace std;

FastaInputStream::~FastaInputStream() {
  if ( file_was_opened ) 
    fin.close();
}

FastaInputStream::FastaInputStream(char * filename)  
{ 
  file_was_opened = false;
  if ( filename == 0 )
   {
	fp = & std::cin;
   }
  else
    {
     fin.open(filename, ifstream::in );
     if ( ! fin.good() )
     {
       fin.close();
       fin.clear();
       THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
     }
     file_was_opened = true;
     fp = & fin;
    }
}



bool
FastaInputStream::read( std::vector<DNA_b128_String> &b128seqs, std::string & runId, std::vector<std::string> &names,  Extrainfos &extrainfos)
{
  std::vector<Sequence> seqs;
  if ( ! readSequences(seqs,runId,extrainfos) )
  {
  	return false;
  }
  names.clear();
  names.reserve(seqs.size());
  for( size_t i=0;i<seqs.size();i++) {
    names.push_back(seqs[i].name);
  }
  Sequences2DNA_b128(seqs,b128seqs);
  return true;
}

bool
FastaInputStream::readSeqName(std::vector<Sequence> &seqs, std::string &line, int linesRead) {

  seqs.resize(linesRead);
  std::string seqStr;
  Sequence &s = seqs[linesRead-1];
  string::size_type findPos = line.find_first_of("\t ",1);
  string::size_type nameEndPos;
  if ( findPos == string::npos ) {
    nameEndPos = line.size() -1;
  }
  else {
    nameEndPos = findPos-1;
  }


  s.name = line.substr(1, nameEndPos);

  return true;
}
bool
FastaInputStream::readSeq(std::vector<Sequence> &seqs, std::string &line, int linesRead) {

  seqs.resize(linesRead);

  Sequence &s = seqs[linesRead-1];
  std::string seqStr = s.seq;

  seqStr += line;

  if (!fp->eof() &&  (seqStr.size() == 0 || seqStr.find_first_not_of("acgtumrwsykvhdbnxACGTUMRWSYKVHDBNX -.?") != string::npos ) ) {
    THROW_EXCEPTION("Malformed Fasta format3\n");
    exit(EXIT_FAILURE);
  } else {
     Sequence &s = seqs[linesRead-1];
     s.seq = seqStr;
  }

  return true;
}

bool
FastaInputStream::readSequences(std::vector<Sequence> &seqs,  std::string & runId, Extrainfos &extrainfos) {
  std::string line;
  int linesRead = 0;
  while ( getline ( *fp, line ) ) {
      line.erase(line.find_last_not_of(" \n\r\t")+1);
      if ( line.size() > 0 && line[0] == '>' ) {
        ++linesRead;
        readSeqName( seqs , line ,linesRead );

      }else{
      	readSeq( seqs , line ,linesRead );
      }
  }
  return true;
   //  THROW_EXCEPTION("Malformed Fasta2 format\n");
     //exit(EXIT_FAILURE);
}

/*#include "FastaInputStream.hpp"
#include <cstdio>

using namespace std;

FastaInputStream::~FastaInputStream() {
  if ( file_was_opened )
    fin.close();
}

FastaInputStream::FastaInputStream(char * filename = 0 )
{
  file_was_opened = false;
  if ( filename == 0 )
    {
   fp = & std::cin;    }
  else
    {
     fin.open(filename, ifstream::in );
     if ( ! fin.good() )
     {
       fin.close();
       fin.clear();
       THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
     }
     file_was_opened = true;
     fp = & fin;
    }
}



bool
FastaInputStream::read( std::vector<DNA_b128_String> &b128seqs, std::string & runId, std::vector<std::string> &names,  Extrainfos &extrainfos) 
{
  std::vector<Sequence> seqs;
  if ( ! readSequences(seqs,runId,extrainfos) ) return false;
  names.clear();names.reserve(seqs.size());
  for( size_t i=0;i<seqs.size();i++) {
    names.push_back(seqs[i].name);
  }
  Sequences2DNA_b128(seqs,b128seqs);
  return true;
}

bool
FastaInputStream::readSeq(std::vector<Sequence> &seqs, std::string &line, int linesRead) {

  linesRead++;
  seqs.resize(linesRead);
  std::string seqStr;
  Sequence &s = seqs[linesRead-1];
  string::size_type findPos = line.find_first_of("\t ",1);
  string::size_type nameEndPos;  
  if ( findPos == string::npos ) {
    nameEndPos = line.size() -1;
  }
  else {
    nameEndPos = findPos-1;
  }


  s.name = line.substr(1, nameEndPos);

  bool readGreaterThan = false;
  while ( ! readGreaterThan && getline ( *fp, line ) ) {
    //       c != 'n' && line[0] != '\0' && ! fp->eof() && fp->good()
    if ( line.size() > 0 ) {
      if ( line[0] == '>' ) {
        readGreaterThan = true;
        readSeq( seqs, line, linesRead );        
      } else { 
	seqStr += line;
      }
    }
  }
  if ( seqStr.size() == 0 || seqStr.find_first_not_of("acgtumrwsykvhdbnxACGTUMRWSYKVHDBNX -.?") != string::npos ) {
    THROW_EXCEPTION("Malformed Fasta format\n");
    exit(EXIT_FAILURE);
  } else {
     Sequence &s = seqs[linesRead-1];
     s.seq = seqStr;
  }
  return true;
}


bool
FastaInputStream::readSequences(std::vector<Sequence> &seqs,  std::string & runId, Extrainfos &extrainfos) {

  std::string line;
  while ( getline ( *fp, line ) ) {
      if ( line.size() > 0 && line[0] == '>' ) {
        readSeq( seqs , line ,0 );        
        return true;
      }
  }
     THROW_EXCEPTION("Malformed Fasta format\n");
     exit(EXIT_FAILURE);
}
*/
