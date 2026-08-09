//--------------------------------------------------
//
// File: BitVector.cpp
//
// Author: Isaac Elias
// e-mail: isaac@nada.kth.se
//
//--------------------------------------------------

#include "fastphylo/core/BitVector.hpp"
#include "fastphylo/core/log_utils.hpp"

using namespace std;

BitVector::BitVector(size_t numBits) : bits((numBits / 32) + 1, 0)
{
    this->numBits = numBits;
}

BitVector::BitVector(const BitVector &copyOf) : bits(copyOf.bits)
{
    this->numBits = copyOf.numBits;
}

BitVector &BitVector::operator=(const BitVector &bv)
{
    numBits = bv.numBits;
    bits = bv.bits;
    return *this;
}

std::ostream &operator<<(std::ostream &os, const BitVector &bv)
{
    for (size_t i = 0; i < bv.getNumBits(); i++)
    {
        if (i % 32 == 0)
        {
            os << " ";
        }
        os << bv.getBit(i);
    }
    return os;
}

void BitVector::flippAllBits()
{

    for (auto &bit : bits)
    {
        bit = (bit ^ 0xffFFffFFU);
    }
}

void BitVector::setAllBits()
{

    for (auto &bit : bits)
    {
        bit = 0xffFFffFFU;
    }
}
void BitVector::clearAllBits()
{

    for (auto &bit : bits)
    {
        bit = 0;
    }
}
void BitVector::bitwiseAnd(const BitVector &bv)
{
    size_t end = getLastCombinedHolderPosition(bv);
    for (size_t i = 0; i < end; i++)
    {
        bits[i] = bits[i] & bv.bits[i];
    }
}
void BitVector::bitwiseOr(const BitVector &bv)
{
    size_t end = getLastCombinedHolderPosition(bv);
    for (size_t i = 0; i < end; i++)
    {
        bits[i] = bits[i] | bv.bits[i];
    }
}
void BitVector::bitwiseXor(const BitVector &bv)
{
    size_t end = getLastCombinedHolderPosition(bv);
    for (size_t i = 0; i < end; i++)
    {
        bits[i] = bits[i] ^ bv.bits[i];
    }
}
void BitVector::bitwiseEqual(const BitVector &bv)
{
    size_t end = getLastCombinedHolderPosition(bv);
    for (size_t i = 0; i < end; i++)
    {
        bits[i] = (bits[i] ^ bv.bits[i]) ^ 0xffFFffFFU;
    }
}
