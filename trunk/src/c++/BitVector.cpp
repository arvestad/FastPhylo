//--------------------------------------------------
//                                        
// File: BitVector.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: BitVector.cpp,v 1.9 2006/12/08 11:09:12 isaac Exp $                                 
//
//--------------------------------------------------


#include "BitVector.hpp"
#include "log_utils.hpp"

using namespace std;


BitVector::BitVector(size_t numBits) : bits(numBits/32+1,0){
  this->numBits = numBits;
}

BitVector::BitVector(const BitVector &copyOf) : bits(copyOf.bits){
  this->numBits = copyOf.numBits;
}

BitVector& 
BitVector::operator=(const BitVector &bv){
  numBits = bv.numBits;
  bits = bv.bits;
  return *this;
}

std::ostream& 
BitVector::printOn(std::ostream& os) const{
 
  for(size_t i=0 ; i<numBits ; i++){
    if( i%32 == 0 )
      os << " ";
    os << getBit(i);
  }
  return os;
}

void 
BitVector::flippAllBits(){

  for(size_t i = 0 ; i<bits.size() ; i++){
    bits[i] = (bits[i] ^ 0xffFFffFFU);
  }

}


void 
BitVector::setAllBits(){

  for(size_t i = 0 ; i<bits.size() ; i++){
    bits[i] = 0xffFFffFFU;
  }
}
void 
BitVector::clearAllBits(){
  
  for(size_t i = 0 ; i<bits.size() ; i++){
    bits[i] = 0;
  }
}
void 
BitVector::bitwiseAnd(const BitVector &bv){
  size_t end = getLastCombinedHolderPosition(bv);
  for(size_t i = 0 ; i<end ; i++){
    bits[i] = bits[i] & bv.bits[i];
  }
}
void 
BitVector::bitwiseOr(const BitVector &bv){
  size_t end = getLastCombinedHolderPosition(bv);
  for(size_t i = 0 ; i<end ; i++){
    bits[i] = bits[i] | bv.bits[i];
  }

}
void 
BitVector::bitwiseXor(const BitVector &bv){
  size_t end = getLastCombinedHolderPosition(bv);
  for(size_t i = 0 ; i<end ; i++){
    bits[i] = bits[i] ^ bv.bits[i];
  }
}
void 
BitVector::bitwiseEqual(const BitVector &bv){
  size_t end = getLastCombinedHolderPosition(bv);
  for(size_t i = 0 ; i<end ; i++){
    bits[i] = (bits[i] ^ bv.bits[i]) ^ 0xffFFffFFU;
  }

}

size_t 
BitVector::hashCode() const{
  unsigned long long h = 0;
  
  for(size_t i = 0 ; i<bits.size()-1 ; i++){
    h = ((h << 32) | bits[i]) % 2654435761U;
  }

  h = ((h << 32) | getLastHolderWithClearedUnusedBitPositions()) % 2654435761U;

  //  cout <<( unsigned long) h << endl;
  return (unsigned long) h;
}

bool   
BitVector::equals(const Object *o) const{
  BitVector *bv = (BitVector*) o;

  if(bv==NULL) return false;
  if(numBits!=bv->numBits) return false;
  
  for(size_t i = 0 ; i<bits.size()-1 ; i++){
    if(bits[i]!=bv->bits[i]) return false;
  }
  if(getLastHolderWithClearedUnusedBitPositions()!=bv->getLastHolderWithClearedUnusedBitPositions()) 
    return false;

  return true;
}









