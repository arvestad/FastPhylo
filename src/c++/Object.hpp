///////////////////////////////////////////////
//
// File: Object.hpp
//
// Author: Isaac Elias
//
// cvs: $Id: Object.hpp,v 1.10 2006/12/08 11:09:12 isaac Exp $
///////////////////////////////////////////////

#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <string>
#include <iostream>
#include "log_utils.hpp"

class Object
{
public:
    virtual ~Object(){};

  //PRINTING
  virtual std::ostream& printOn(std::ostream& os) const;
  //calls printOn
  virtual std::string toString() const;

  //empty init function
  virtual std::istream& objInitFromStream(std::istream &is){return is;}

  //overwrite
  virtual size_t hashCode() const { return (size_t) this;}
  virtual bool equals(const Object *o) const {return (void *)this == (void*)o;}

  //don't overwrite
  virtual bool equals(const Object &o) const {return this->equals(&o);}
  virtual bool operator==(const Object &o) const { return this->equals(&o);}
  virtual bool operator!=(const Object &o) const { return !(this->equals(&o));}
};



struct objeq
{
  bool operator()(const Object & o1, const Object & o2) const
  {
    return o1.equals(o2);
  }
};
struct objhash{
  size_t operator()(const Object &o) const{
    return o.hashCode();
  }
};

struct objeq_ptr
{
  bool operator()(const Object * o1, const Object * o2) const
  {
    return o1->equals(o2);
  }
};
struct objhash_ptr{
  size_t operator()(const Object *o) const{
    return o->hashCode();
  }
};

struct objeq_ref
{
  bool operator()(const Object &o1, const Object &o2) const
  {
    return o1.equals(o2);
  }
};
struct objhash_ref{
  size_t operator()(const Object &o) const{
    return o.hashCode();
  }
};



/*
 * Out stream operators
 */
std::ostream& operator<<(std::ostream & os,  const Object &obj);
std::ostream& operator<<(std::ostream & os,  const Object *obj);

/*
 * Calls the Object.initFromStream
 */
std::istream& operator>>(std::istream & in,  Object &obj);
std::istream& operator>>(std::istream & in,  Object *obj);


std::string operator + (const std::string & s,  const Object & obj);
std::string operator + (const std::string & s,  const Object * obj);
std::string operator + ( const Object & obj, const std::string & s);
std::string operator + ( const Object * obj, const std::string & s);

#endif // OBJECT_HPP









