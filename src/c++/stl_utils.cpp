
#include "stl_utils.hpp"
#include <string>
#include <iostream>
#include "file_utils.hpp"


using namespace std;



string operator + (const string & s, const int i) {
    // a 32bit int in decimal asciiz form is less than 12 characters
    char buf[12];
#ifdef WIN32
    _snprintf(buf, sizeof(buf), "%d", i);
#else
    snprintf(buf, sizeof(buf), "%d", i);
#endif
    string t = s;
    t.append(buf);
    return t;
}

string operator + (const int i, const string & s) {
    // a 32bit int in decimal asciiz form is less than 12 characters
    char buf[12];
#ifdef WIN32
    _snprintf(buf, sizeof(buf), "%d", i);
#else
    snprintf(buf, sizeof(buf), "%d", i);
#endif
    string t = s;
    t.insert(0,buf);
    return t;
}


string operator + (const string & s, const float i) {
    // a 32bit int in decimal asciiz form is less than 12 characters
    char buf[20];
#ifdef WIN32
    _snprintf(buf, sizeof(buf), "%f", i);
#else
    snprintf(buf, sizeof(buf), "%f", i);
#endif
    string t = s;
    t.append(buf);
    return t;
}

string operator + (const float i, const string & s) {
    // a 32bit int in decimal asciiz form is less than 12 characters
    char buf[20];
#ifdef WIN32
    _snprintf(buf, sizeof(buf), "%f", i);
#else
    snprintf(buf, sizeof(buf), "%f", i);
#endif
    string t = s;
    t.insert(0,buf);
    return t;
}


std::istream& operator>>(std::istream &in, std::string *str){
  in >> *str;
  return in;
}

std::ostream& operator<<(std::ostream & os, std::string *str){
  os << *str;
  return os;
}


void
print_map(str2int_map &m){

  cout << "map_size = " << m.size() << endl;
  str2int_map::iterator i = m.begin();
  for ( ; i != m.end() ; ++i ){
    cout << (*i).first << " -> " <<(*i).second << endl; 
  }
}




