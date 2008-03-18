//--------------------------------------------------
//                                        
// File: Sequence.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Sequence.hpp,v 1.9 2006/12/25 18:40:39 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <string>
#include "stl_utils.hpp"
#include "Object.hpp"
#include "log_utils.hpp"
#include <fstream>
#include "file_utils.hpp"

//
// A simple class for handling named sequences.
// 
// This class allows for reading and writing of PHYLIP sequence files
// and also bootstrapping.
//
//


class Sequence : public Object{
public:

  static hashstr stringhasher;
  
  std::string name;
  std::string seq;
  
  Sequence();
  Sequence(std::string n, std::string s);
  Sequence(const Sequence &s);

  Sequence& operator=(const Sequence &s);

  virtual std::istream& objInitFromStream(std::istream &in);  

  //------------------------------------------
  //PRINTING
  //Ex. "HUMAN    agct-agct"
  virtual std::ostream& printOn(std::ostream& os) const; 
  //print only name
  virtual std::ostream& printShort(std::ostream& os) const;
  virtual void printWithoutGaps(std::ostream& os) const;


  //------------------------------------------
  // Hashes on the name if it exists otherwise the sequence
  virtual size_t hashCode() const;
  // uses the name if it exists otherwise the sequence
  virtual bool equals(const Object *o) const;

  //------------------------------------------
  //returns true if all characters in the sequences are in chars.
  bool onlyContains(std::string &chars);

  //----------------------------
  // READING PHYLIP SEQUENCEFILE
  // The vector is cleared of all sequences and new  sequences are added
  static void readSequences(std::vector<Sequence> &seqs, std::istream &in);
  static void printSequences(std::vector<Sequence> &seqs, std::ofstream &out);

  //-------------------------------------------
  //BOOTSTRAPPING
  //creates a bootstrapped data set from the sequences in seqs.
  static void bootstrapSequences(std::vector<Sequence> &seqs, std::vector<Sequence> &boot);

};



#endif // SEQUENCE_HPP








