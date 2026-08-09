//--------------------------------------------------
//
// File: DistanceMatrix_impl.hpp
//
// Author: Isaac Elias
// e-mail: isaac@nada.kth.se
//
//--------------------------------------------------
#pragma once

#include <iomanip>
#include <string>

DM_TEMPLATE void DISTANCEMATRIX::assureSize()
{
    identifiers.resize(size);
    D.resize(size);
    for (size_t i = 0; i < size; i++)
        D[i].resize(size);
}

DM_TEMPLATE
DISTANCEMATRIX::DistanceMatrix(size_t size)
{
    this->size = size;
    assureSize();
}

DM_TEMPLATE
DISTANCEMATRIX::DistanceMatrix(const DISTANCEMATRIX &dm)
{
    (*this) = dm;
}

DM_TEMPLATE DISTANCEMATRIX &DISTANCEMATRIX::operator=(const DISTANCEMATRIX &dm)
{
    size = dm.size;

    identifiers.resize(size);
    D.resize(size);

    // copy the upper triangle
    for (size_t i = 0; i < size; i++)
    {
        identifiers[i] = dm.identifiers[i];
        DISTVEC &thisvec = D[i];
        thisvec.resize(size);
        const DISTVEC &thatvec = dm.D[i];
        for (size_t j = i; j < size; j++)
            thisvec[j] = thatvec[j];
    }

    return *this;
}

DM_TEMPLATE void DISTANCEMATRIX::setDefaultValues(DistanceType &defval, Identifier &defid)
{

    for (size_t i = 0; i < size; i++)
        identifiers[i] = defid;

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j <= i; j++)
            D[j][i] = defval;
    }
}

DM_TEMPLATE std::ostream &operator<<(std::ostream &out, const DISTANCEMATRIX &dm)
{
    out << std::setw(5) << std::right << dm.size << std::endl;

    out.precision(6);
    for (size_t i = 0; i < dm.size; i++)
    {
        out << std::setw(10) << std::left;
        dm.idPrintOn(out, dm.identifiers[i]);
        out << "  ";
        size_t j = 0;
        for (; j < i; j++)
        {
            out << std::setw(10) << std::right;
            dm.distPrintOn(out, dm.D[j][i]);
            out << " ";
        }
        for (; j < dm.size; j++)
        {
            out << std::setw(10) << std::right;
            dm.distPrintOn(out, dm.D[i][j]);
            out << " ";
        }
        out << std::endl;
    }

    return out;
}

// SWAP AND REMOVE LAST ROW
DM_TEMPLATE void DISTANCEMATRIX::swapRowToLast(size_t row)
{

    assert(row < size);
    size_t lastRow = size - 1;
    if (row == lastRow)
        return;

    // switch the distances
    for (size_t i = 0; i < row; i++)
    {
        DistanceType tmp = D[i][row];
        D[i][row] = D[i][lastRow];
        D[i][lastRow] = tmp;
    }

    DistanceType tmp = D[row][row];
    D[row][row] = D[lastRow][lastRow];
    D[lastRow][lastRow] = tmp;

    for (size_t i = row + 1; i < lastRow; i++)
    {
        DistanceType tmp = D[row][i];
        D[row][i] = D[i][lastRow];
        D[i][lastRow] = tmp;
    }

    // Switch the Identifiers
    Identifier tmpI = identifiers[row];
    identifiers[row] = identifiers[lastRow];
    identifiers[lastRow] = tmpI;
}

DM_TEMPLATE void DISTANCEMATRIX::removeLastRow()
{
    size--;
}
