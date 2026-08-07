//--------------------------------------------------
//                                        
// File: file_utils.cpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: file_utils.cpp,v 1.4 2006/09/10 10:18:26 isaac Exp $                                 
//
//--------------------------------------------------

#include "fastphylo/core/file_utils.hpp"
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/core/Exception.hpp"

using namespace std;

/* CHECK IF FILE EXISTS */
int
file_exists(const char *fname){
	FILE *ftmp = fopen(fname,"r");
	if ( ftmp != nullptr ){
		fclose(ftmp);
		return 1;
	}
	else { {
		return 0;
}
}
}



FILE *
open_write_file_interactive(const char *fname){
	// Binary mode, not "w"/"a" (text mode) - on Windows the CRT
	// translates every outgoing '\n' to "\r\n" in text mode, silently
	// corrupting output that is meant to be byte-identical across
	// platforms (this project's own code always writes '\n' explicitly,
	// so it never relied on that translation for correctness).
	const char *mode ="wb";
	if ( file_exists(fname) != 0 ){
		while ( true ){
			cout << "File exists: \"" << fname << "\"" << endl;
			cout << "Do you want to (w)rite or (a)ppend? "<<endl;
			string choice;
			cin >> choice;
			if ( choice == "w" ){
				mode = "wb";
				break;
			}
			else if ( choice == "a" ){
				mode = "ab";
				break;
			}
			cout << "bad choice" << endl;
		}
	}

	string tmp;
	FILE *ftmp = fopen(fname,mode);
	if ( ftmp == nullptr ){
		cout << "Couldn't open file \"" << fname << "\"" << endl;
		cout << "What file do you want to write to? " << endl;
		cin >> tmp;
		fname = tmp.c_str();
		return open_write_file(fname);
	}
	else { {
		return ftmp;
}
}
}



FILE *
open_write_file(const char *fname){
	// Binary mode - see open_write_file_interactive()'s comment above.
	const char *mode ="wb";

	FILE *ftmp = fopen(fname,mode);
	if ( ftmp == nullptr ){
		THROW_EXCEPTION("Couldn't open file \"" << fname << "\"");
	}
	else { {
		return ftmp;
}
}
}

void
open_write_stream_interactive(const char *fname, ofstream &of){
	// Binary mode - see open_write_file_interactive()'s comment above.
	ofstream::openmode mode = ofstream::out | ofstream::binary;
	if ( file_exists(fname) != 0 ){
		while ( true ){
			cout << "File exists: \"" << fname << "\"" << endl;
			cout << "Do you want to (w)rite or (a)ppend? "<<endl;
			string choice;
			cin >> choice;
			if ( choice == "w" ){
				mode = ofstream::out | ofstream::binary;
				break;
			}
			else if ( choice == "a" ){
				mode = ofstream::app | ofstream::binary;
				break;
			}
			cout << "bad choice" << endl;
		}
	}

	of.open(fname,mode);

	if ( !of.good() ){
		of.close();
		of.clear();
		USER_WARNING("can't open file " << string(fname));
		return;
	}
}

void
open_write_stream(const string fname, ofstream &of){
	open_write_stream(fname.c_str(),of);
}

void
open_write_stream(const char *fname, ofstream &of){
	// Binary mode - see open_write_file_interactive()'s comment above.
	ofstream::openmode mode = ofstream::out | ofstream::binary;

	of.open(fname,mode);

	if ( !of.good() ){
		of.close();
		of.clear();
		THROW_EXCEPTION("Can't open file " << fname);
	}
}


FILE *
open_read_file_interactive(const char *fname){
	string tmp;
	while ( true ){
		FILE *ftmp = fopen(fname,"r");
		if ( ftmp != nullptr ) {
			return ftmp;
}

		cout << "File doesn't exist: \"" << fname << "\"" << endl;
		cout << "What file do you want to read? " << endl;

		cin >> tmp;
		fname = tmp.c_str();
	}
}

FILE *
open_read_file(const char *fname){
	FILE *ftmp = fopen(fname,"r");
	if ( ftmp != nullptr ) {
		return ftmp;
}
	THROW_EXCEPTION("File doesn't exist: \"" << fname <<"\"");
}

// Added by Mehmood Khan Malagori; email: malagori@kth.se
std::ofstream *
open_write_binary(const char *fname){
	auto *ofs = new ofstream(fname, ios::binary);
	if ( !ofs->good() ){
		ofs->close();
		ofs->clear();
		THROW_EXCEPTION("Can't open file " << fname);
		ofs = nullptr;
	}

	return ofs;
}
// end open_write_binary()

void
open_read_stream_interactive(const char *fname, ifstream &fin){

	string tmp;
	while ( true ){
		fin.open(fname,ifstream::in);

		if ( fin.good() ) {
			return;
}
		fin.close();
		fin.clear();
		cout << "File doesn't exist: \"" << fname << "\"" << endl;
		cout << "What file do you want to read? " << endl;

		cin >> tmp;
		fname = tmp.c_str();
	}
}


void
open_read_stream(const string fname, ofstream &fin){
	open_read_stream(fname,fin);
}

void
open_read_stream(const char *fname, ifstream &fin){
	fin.open(fname,ifstream::in);

	if ( fin.good() ) {
		return;
}

	fin.close();
	fin.clear();

	THROW_EXCEPTION("File doesn't exist: \"" << fname << "\"");
}


void
skipWhiteSpace(FILE *f){
	char c = fgetc(f);
	while ( isspace(c) != 0 ) {
		c = fgetc(f);
}

	ungetc(c,f);
}

void
skipWhiteSpace(std::istream &in){
	in >> ws;
}


void
skipUntil(std::istream &in, char *chars){

	char c;
	in.get(c);
	while ( in.good() && strchr(chars,c) == nullptr ){
		in.get(c);
	}
	in.unget();
}

void
appendToken(std::istream &in, std::string &str){
	char c;
	in.get(c);
	while ( in.good() && (isspace(c) == 0) ){
		str.push_back(c);
		in.get(c);
	}
	in.unget();
}


void
appendUntil(std::istream &in, std::string &str,  char *chars){
	char c;
	in.get(c);
	while (  in.good() && strchr(chars,c) == nullptr ){
		str.push_back(c);
		in.get(c);
	}
	in.unget();


}
