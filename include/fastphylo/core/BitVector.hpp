///////////////////////////////////////////////
//
// File: BitVector.hpp
//
// Author: Isaac Elias
//
///////////////////////////////////////////////

#pragma once

#include <string>
#include <iostream>
#include "fastphylo/core/Object.hpp"
#include <vector>

#define BITHOLDER_POS(POS) ((POS) & 0xfffff)
#define BITVEC_POS(POS) ((POS) >> 5)

class BitVector : public Object
{

  public:
    BitVector(size_t numBits = 0);
    BitVector(const BitVector &copyOf);
    BitVector &operator=(const BitVector &bv);

    size_t getNumBits() const
    {
        return numBits;
    }
    void setNumBits(size_t i)
    {
        numBits = i;
        bits.resize(numBits / 32 + 1, 0);
    }
    int getBit(size_t pos) const
    {
        size_t bitholder = bits[BITVEC_POS(pos)];
        size_t intval = (bitholder >> BITHOLDER_POS(pos)) & 0x1;
        return intval;
    }
    void setBit(size_t pos)
    {
        size_t &bitholder = bits[BITVEC_POS(pos)];
        bitholder = bitholder | (0x1 << (BITHOLDER_POS(pos)));
    }
    void clearBit(size_t pos)
    {
        size_t &bitholder = bits[BITVEC_POS(pos)];
        bitholder = bitholder & (~(0x1 << (BITHOLDER_POS(pos))));
    }
    void setAllBits();
    void clearAllBits();
    void flippAllBits();
    void flippAllIfPositionIsCleared(size_t pos)
    {
        if (getBit(pos) == 0)
            flippAllBits();
    }

    // if different length the operation is performed only on the first
    // bits which the vectors have in common
    void bitwiseAnd(const BitVector &bv);
    void bitwiseOr(const BitVector &bv);
    void bitwiseXor(const BitVector &bv);
    // all bits in which the vectors agree are set
    void bitwiseEqual(const BitVector &bv);

    std::ostream &printOn(std::ostream &os) const override;
    std::istream &objInitFromStream(std::istream &is) override
    {
        return is;
    }

  private:
    std::vector<size_t> bits;
    size_t numBits;

    size_t getLastHolderWithClearedUnusedBitPositions() const
    {
        size_t rest = numBits % 32;
        size_t p = bits.size() - 1;
        return ((bits[p] << (32 - rest)) >> (32 - rest));
    }

    size_t getLastCombinedHolderPosition(const BitVector &bv) const
    {
        size_t end = bits.size();
        end = (end <= bv.bits.size() ? end : bv.bits.size());
        return end;
    }
};

// Found via ADL - what `os << someBitVector` actually binds to (an
// exact-match free function beats the inherited, conversion-requiring
// Object::operator<<). printOn() above forwards here.
std::ostream &operator<<(std::ostream &os, const BitVector &bv);
