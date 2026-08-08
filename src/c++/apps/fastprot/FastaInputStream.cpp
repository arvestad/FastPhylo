#include "FastaInputStream.hpp"

using namespace std;

FastaInputStream::FastaInputStream(char *filename)
    : reader(filename, "abcdefghiklmnopqrstuvwyzxABCDEFGHIKLMNOPQRSTUVWYZX -.?")
{
}

bool FastaInputStream::read(vector<Sequence> &seqs, string &runId, vector<string> &names, Extrainfos &extrainfos)
{
    if (!readSequences(seqs, runId, extrainfos))
    {
        return false;
    }
    names.clear();
    names.reserve(seqs.size());
    for (auto &seq : seqs)
    {
        names.push_back(seq.name);
    }
    return true;
}

bool FastaInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos)
{
    return reader.readSequences(seqs, runId, extrainfos);
}
