///////////////////////////////////////////////
//
// File: Object.hpp
//
// Author: Isaac Elias
//
///////////////////////////////////////////////

#pragma once

#include <string>
#include <iostream>

class Object
{
  public:
    virtual ~Object()
    {
    }

    // PRINTING
    virtual std::ostream &printOn(std::ostream &os) const;
    // calls printOn
    virtual std::string toString() const;

    // empty init function
    virtual std::istream &objInitFromStream(std::istream &is)
    {
        return is;
    }

    // overwrite
    virtual size_t hashCode() const
    {
        return reinterpret_cast<size_t>(this);
    }
    virtual bool equals(const Object *o) const
    {
        return static_cast<const void *>(this) == static_cast<const void *>(o);
    }

    // don't overwrite
    virtual bool equals(const Object &o) const
    {
        return this->equals(&o);
    }
    virtual bool operator==(const Object &o) const
    {
        return this->equals(&o);
    }
    virtual bool operator!=(const Object &o) const
    {
        return !(this->equals(&o));
    }
};

struct objeq
{
    bool operator()(const Object &o1, const Object &o2) const
    {
        return o1.equals(o2);
    }
};
struct objhash
{
    size_t operator()(const Object &o) const
    {
        return o.hashCode();
    }
};

/*
 * Out stream operators
 */
std::ostream &operator<<(std::ostream &os, const Object &obj);
// Needed for DistanceMatrix<T*, ...>/Data_init<T*>::operator() to compile
// when T* is an Object-derived pointer (e.g. NJMatrix's TreeNode<...>*):
// the class template's virtual objInitFromStream()/printOn() overrides are
// instantiated for the vtable regardless of whether they're ever actually
// called, and Data_init<T*>/Data_printOn<T*> stream through this pointer
// overload via the derived-to-base implicit conversion.
std::ostream &operator<<(std::ostream &os, const Object *obj);

/*
 * Calls the Object.initFromStream
 */
std::istream &operator>>(std::istream &in, Object *obj);
