//--------------------------------------------------
//                                        
// File: Sequence.cpp                              
//                             
// Author: Isaac Elias, Mehmood Alam Khan
// e-mail: isaac@nada.kth.se, malagori@kth.se
//                             
// cvs: $Id: Sequence.cpp,v 1.45 2006/12/26 11:57:36 isaac Exp $                                 
//
//--------------------------------------------------

#include "Sequence.hpp"
#include <string>
#include "log_utils.hpp"
#include <fstream>
#include "nucleotide.hpp"
#include "stl_utils.hpp"

using namespace std;

hashstr Sequence::stringhasher;

Sequence::Sequence(){
	name = "";
	seq = "";
}
Sequence::Sequence(std::string n, std::string s){
	name = n;
	seq = s;
}
Sequence::Sequence(const Sequence &s){
	name = s.name;
	seq = s.seq;
}

Sequence& Sequence::operator=(const Sequence &s){
	name = s.name;
	seq = s.seq;
	return *this;
}

std::ostream& Sequence::printOn(std::ostream& os) const{
	if ( name.size() == 0 && seq.size() == 0 )
		return os;
	if ( name.length() < 10 ){
		os  << std::setw(10) << std::left;
		os << name;
	}
	else
		os << name << " ";
	os << seq;
	return os;
}

std::ostream& Sequence::printShort(std::ostream& os) const{
	os << name;
	return os;
}

std::istream& Sequence::objInitFromStream(std::istream &in){
	name.clear();
	in >> name;
	seq.clear();
	in >> seq;
	return in;
}

void Sequence::printWithoutGaps(std::ostream& os) const{
	if ( name.size() == 0 && seq.size() == 0 )
		return;
	if ( name.length() < 10 ){
		os  << std::setw(10) << std::left;
		os << name;
	}
	else
		os << name << " ";

	for(size_t i=0 ; i<seq.length() ; i++ ){
		if(seq[i]!=' ')
			os<<seq[i];
	}

}

// Hashes on the name if it exists otherwise the sequence
size_t Sequence::hashCode() const {
	if ( name.size() == 0 )
		return stringhasher(seq.c_str());
	return stringhasher(name.c_str());
}
// uses the name if it exists otherwise the sequence
bool Sequence::equals(const Object *o) const {
	const Sequence *otherseq = (const Sequence *)o;
	if ( name.size() == 0 ){
		if ( otherseq->name.size() == 0 )
			return seq == otherseq->seq;
		else
			return false;
	}
	if ( otherseq->name.size() == 0 )
		return false;

	return name == otherseq->name;
}

bool Sequence::onlyContains(std::string &chars){
	for ( size_t i = 0 ; i < seq.size() ; i++ )
		if ( chars.find(seq[i]) != chars.npos )
			return false;

	return true;
}


// //-----------------------------------------------
// void
// Sequence::readSequences(std::vector<Sequence> &seqs, ifstream &fin){
//   int numSequences;
//   unsigned int seqlen;
//   char tmp[100];
//   fin.getline(tmp,100);
//   sscanf(tmp,"%d %d",&numSequences,&seqlen);

//   const size_t startindex = seqs.size();
//   seqs.resize(seqs.size()+numSequences);

//   //read the names and map the sequences onto the tree
//   for ( int i = 0 ; i < numSequences ; i++ ){
//     Sequence &s = seqs[startindex + i];
//     fin >> s.name;
//     s.seq.clear();
//     s.seq.reserve(seqlen+5);

//     while (1){
//       char c = fin.get();
//       nucleotide n = char2nucleotide(c);
//       if ( DNA_NOT_ALLOWED == n ){
//         if ( !isspace(c) ){
//           USER_ERROR("Bad character \'" << c << "\'");
//         }
//         else if ( c != '\n' )
//           continue;
//         //if '\n'
//         break;
//       }

//       s.seq.push_back(nucleotide2char(n));
//     }
//   }
//   //read remaining sequences

//   //The sequences aren't neccesarily on one line but my be spread out interleaving
//   //over several lines. Therefore we read until seqlen chars have been read.
//   while ( seqs[0].seq.length() < seqlen ){
//     for ( int i = 0 ; i < numSequences ; i++ ){
//       char c = fin.peek();
//       if ( !isspace(c) ){
//         //skip first 10 chars
//         for ( int j = 10 ; j != 0 ; j-- )
//           fin.get();
//       }
//       while ( c == '\n' )
//         c = fin.get();

//       Sequence &s = seqs[startindex+i];
//       while (1){
//         c = fin.get();
//         nucleotide n = char2nucleotide(c);
//         if ( DNA_NOT_ALLOWED == n ){
//           if ( !isspace(c) ){
//             USER_ERROR("Bad character \'" << c << "\'");
//           }
//           else if ( c != '\n' )//skip space
//             continue;
//           //if '\n'
//           break;
//         }

//         s.seq.append(1,nucleotide2char(n));
//       } 
//     }
//   }


//   // CHECK THAT ALL STRINGS HAVE THE SAME LENGTH
//   for ( int i = 0 ; i < numSequences ; i++ )
//     if ( seqs[startindex+i].seq.length() != seqlen ){
//       USER_ERROR("Sequence not of correct length: " << seqs[i].name << "    length is " << seqs[i].seq.length()); 
//     }
// }


void
Sequence::readSequences(std::vector<Sequence> &seqs, istream &fin){
	int numSequences;
	unsigned int seqlen;
	const int MAXLINE = 16384;
	char line[MAXLINE];
	do {//skip lines which does not contain two numbers.
		fin.getline(line,MAXLINE);
	}while(sscanf(line,"%d %d",&numSequences,&seqlen) != 2 );

	seqs.resize(numSequences);

	//read the names and sequences
	char tmpName[11];
	for ( int i = 0 ; i < numSequences ; i++ ){
		Sequence &s = seqs[i];
		s.seq.clear();
		s.seq.reserve(seqlen+1);

		fin.getline(tmpName,11);//reads atmost 10 chars
		//skip lines without 10 chars per line
		if( !fin.fail() ){
			if( fin.eof() ) THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
			i--;
			continue;
		}
		s.name.clear();
		appendUntil(s.name,tmpName,fin.gcount(), ' ');

		fin.clear();
		fin.getline(line,MAXLINE);

		while( fin.fail() && fin.gcount()==MAXLINE-1 ){//didn't read all the line
			appendAllNonChars(s.seq,line,fin.gcount(), ' ');
			fin.clear();
			fin.getline(line,MAXLINE);
		}
		if( !fin.fail()) {//we read it all including the newline char unless it ended with eof
			appendAllNonChars(s.seq,line, fin.gcount() - (fin.eof()? 0: 1), ' ');
		}
		else //fail
			THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
		//    PRINT_V(s.seq);
	}//end for loop

//Mehmood's Changes here'
		//The sequences aren't neccesarily on one line but my be spread out interleaving
		//over several lines. Therefore we read until seqlen chars have been read.
		bool whileTrue= true;
		// check if all the sequences has been completely read.
		for ( int i = 0 ; i < numSequences ; i++ ){
			if ( seqs[i].seq.length() == seqlen ){
				whileTrue= false;
			}
		}

		while ( whileTrue){
			for ( int i = 0 ; i < numSequences ; i++ ){
				Sequence &s = seqs[i];

				fin.getline(line,MAXLINE);

				std::string myStr;
				myStr=line;

				if (myStr.empty()){
					if( fin.eof() )
						THROW_EXCEPTION("Sequence not of correct length: " << seqs[i].name<< "    length is " << seqs[i].seq.length());
					i=i-1;
					continue;
				}

				if (!( fin.fail() && fin.gcount()==MAXLINE-1 )){
					appendAllNonChars(s.seq,line, fin.gcount() - (fin.eof()? 0: 1), ' ');

				}


				// check if any one of the sequences is read completely and has correct
				// seq.length but some don't have that break the loop with error msg:
				for ( int i = 0 ; i < numSequences ; i++ ){
						if ( seqs[i].seq.length() == seqlen ){
							whileTrue= false;
						}
					}

			}//end for loop
		}

// Mehmood's changes end here

//	//The sequences aren't neccesarily on one line but my be spread out interleaving
//	//over several lines. Therefore we read until seqlen chars have been read.
//	while ( seqs[0].seq.length() < seqlen ){
//		for ( int i = 0 ; i < numSequences ; i++ ){
//			//      PRINT_V(i);
//			Sequence &s = seqs[i];
//			// mehmood's changes
//			//      fin.getline(tmpName,11);//reads atmost 10 chars
//			//      //skip lines without 10 chars per line
//			//      if( !fin.fail() ){
//			//	if( fin.eof() ) THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
//			//	i--;
//			//	continue;
//			//      }
//			//// mehmood's changes
//			fin.clear();
//			fin.getline(line,MAXLINE);
//			//PRINT_V(line);
//			while( fin.fail() && fin.gcount()==MAXLINE-1 ){//didn't read all the line
//				appendAllNonChars(s.seq,line,fin.gcount(), ' ');
//				fin.clear();
//				fin.getline(line,MAXLINE);
//				//PRINT_V(line);
//			}
//			if( !fin.fail()) {//we read it all including the newline char unless it ended with eof
//				appendAllNonChars(s.seq,line,fin.gcount() - (fin.eof()? 0: 1), ' ');
//			}
//			else //fail
//				THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
//			//PRINT_V(s.seq);
//		}//end for loop
//	}

	// CHECK THAT ALL STRINGS HAVE THE SAME LENGTH
	for ( int i = 0 ; i < numSequences ; i++ ){
		if ( seqs[i].seq.length() != seqlen ){
			THROW_EXCEPTION("Sequence not of correct length: " << seqs[i].name << "    length is " << seqs[i].seq.length());
		}
	}
}



void 
Sequence::printSequences(std::vector<Sequence> &seqs, std::ofstream &out){
	out << seqs.size() << "\t " << seqs[0].seq.length() << endl;
	for ( size_t i = 0 ; i < seqs.size() ; i++ ){
		out << seqs[i] << endl;
	}
}

void 
Sequence::bootstrapSequences(std::vector<Sequence> &seqs, std::vector<Sequence> &boot){


	//ensure capacity in bootsequences
	boot.resize(seqs.size());

	const size_t seqlen = seqs[0].seq.length();
	for ( size_t i = 0 ; i < seqs.size() ; i++ ){
		boot[i].name = seqs[i].name;
		if ( boot[i].seq.length() < seqlen )
			boot[i].seq.assign(seqlen,' ');
		else if ( boot[i].seq.length() > seqlen )
			boot[i].seq.erase(seqlen);

		assert( boot[i].seq.length() == seqlen );
	}

	// Do the bootstrapping
	size_t pos=0;
	size_t seq;
	vector<int> samplePositions(seqlen);
	const size_t stride = 32;

	if( stride < seqlen){
		for( pos=0; pos<(seqlen-stride); pos+= stride)
			for( size_t i=0; i<stride; i++)
				samplePositions[pos+i] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));
		for (; pos<seqlen; pos++)
			samplePositions[pos] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));


		for( seq=0; seq<seqs.size(); seq++){
			string & b = boot[seq].seq;
			const string & s = seqs[seq].seq;

			for( pos=0; pos<seqlen-stride; pos+=stride){
				for( size_t i=0; i<stride; i++)
					b[pos+i] = s[samplePositions[pos+i]];
			}
			for (; pos<seqlen; pos++)
				b[pos] = s[samplePositions[pos]];
		}
	}
	else{//seqlen<stride
		for (pos=0; pos<seqlen; pos++)
			samplePositions[pos] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));

		for( seq=0; seq<seqs.size(); seq++){
			string & b = boot[seq].seq;
			const string & s = seqs[seq].seq;
			for (pos=0; pos<seqlen; pos++)
				b[pos] = s[samplePositions[pos]];
		}
	}
}










