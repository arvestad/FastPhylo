#!/usr/bin/env python3
"""Converts an aligned FASTA file to PHYLIP sequential format (names
padded/truncated to 10 chars) for driving PHYLIP's protdist, to
benchmark against fastprot on the same alignment."""
import sys


def main():
    fasta_path, phylip_path = sys.argv[1], sys.argv[2]
    names = []
    seqs = []
    seq = ""
    for line in open(fasta_path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seq:
                seqs.append(seq)
            names.append(line[1:].strip()[:10].ljust(10))
            seq = ""
        else:
            seq += line.strip()
    if seq:
        seqs.append(seq)

    lengths = set(len(s) for s in seqs)
    assert len(lengths) == 1, f"sequences not aligned (lengths={lengths})"
    length = lengths.pop()

    with open(phylip_path, "w") as f:
        f.write(f"  {len(seqs)}   {length}\n")
        for name, s in zip(names, seqs):
            f.write(f"{name}{s}\n")


if __name__ == "__main__":
    main()
