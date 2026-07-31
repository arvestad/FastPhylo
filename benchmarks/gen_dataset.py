#!/usr/bin/env python3
"""Deterministic synthetic protein FASTA generator for benchmarks/tests.

Generates n_seqs sequences of the given length, diverged from a shared
random ancestor at a divergence rate uniformly drawn per-sequence in
[0.05, 0.40] (representative of realistic protein family divergence -
see phase0_audit.md's profiling dataset, which used the same range).
Deterministic for a given (n_seqs, length, seed).
"""
import argparse
import random


def main():
    p = argparse.ArgumentParser()
    p.add_argument("n_seqs", type=int)
    p.add_argument("length", type=int)
    p.add_argument("outfile")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    random.seed(args.seed)
    aa = "ARNDCQEGHILKMFPSTWYV"
    base = ''.join(random.choice(aa) for _ in range(args.length))
    with open(args.outfile, "w") as f:
        for i in range(args.n_seqs):
            s = list(base)
            nmut = random.randint(int(args.length * 0.05), int(args.length * 0.40))
            for _ in range(nmut):
                pos = random.randrange(args.length)
                s[pos] = random.choice(aa)
            f.write(f">seq{i}\n{''.join(s)}\n")


if __name__ == "__main__":
    main()
