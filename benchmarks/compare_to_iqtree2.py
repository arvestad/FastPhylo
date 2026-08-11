#!/usr/bin/env python3
"""Cross-validates fastprot's embedded protein model matrices against
IQ-TREE2's (model/modelprotein.cpp, fetched to /tmp/iqtree_modelprotein.cpp).
Parses IQ-TREE2's lower-triangular exchangeability matrix + frequency
vector format, builds the actual rate matrix Q(i,j)=R(i,j)*pi(j) (off-
diagonal), Q(i,i)=-row sum, and compares element-by-element against
fastprot's dumped matrices (/tmp/fastprot_matrices.txt, from
dump_model_matrices.cpp), after normalizing for the fact that empirical
rate matrices are only defined up to an arbitrary overall time-scale
(one global multiplicative constant per model), which the two sources
need not share.
"""
import re

IQTREE_FILE = "/tmp/iqtree_modelprotein.cpp"
FASTPROT_FILE = "/tmp/fastprot_matrices.txt"

# fastprot's canonical order (ProtSeqCode.hpp): A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
# IQ-TREE2/PAML use the same order for their .dat-style embedded models.


def parse_iqtree_block(text, model_name):
    pattern = re.compile(r"model " + re.escape(model_name) + r"=\n(.*?);", re.DOTALL)
    m = pattern.search(text)
    if not m:
        return None
    lines = [l.strip() for l in m.group(1).strip().split("\n") if l.strip()]
    # First 19 lines: lower-triangular R (line k has k+1 values, k=1..19... actually
    # line 1 has 1 value (R(1,0)), line 2 has 2 (R(2,0),R(2,1)), ..., line 19 has 19 values)
    # Last line: 20 frequencies.
    r_lines = lines[:19]
    freq_line = lines[19]
    freqs = [float(x) for x in freq_line.split()]
    assert len(freqs) == 20, f"{model_name}: expected 20 freqs, got {len(freqs)}"

    R = [[0.0] * 20 for _ in range(20)]
    for i, line in enumerate(r_lines, start=1):
        vals = [float(x) for x in line.split()]
        assert len(vals) == i, f"{model_name} row {i}: expected {i} values, got {len(vals)}"
        for j, v in enumerate(vals):
            R[i][j] = v
            R[j][i] = v

    Q = [[0.0] * 20 for _ in range(20)]
    for i in range(20):
        row_sum = 0.0
        for j in range(20):
            if i != j:
                Q[i][j] = R[i][j] * freqs[j]
                row_sum += Q[i][j]
        Q[i][i] = -row_sum
    return Q, freqs


def parse_fastprot_dump(text, model_name):
    lines = text.strip().split("\n")
    idx = lines.index(f"{model_name}_Q")
    q_rows = lines[idx + 1 : idx + 21]
    Q = [[float(x) for x in row.split()] for row in q_rows]
    idx2 = lines.index(f"{model_name}_EQ")
    freqs = [float(x) for x in lines[idx2 + 1].split()]
    return Q, freqs


def compare_model(iqtree_text, fastprot_text, model_name):
    result = parse_iqtree_block(iqtree_text, model_name)
    if result is None:
        print(f"{model_name}: NOT FOUND in IQ-TREE2 source (skipping)")
        return
    Q_ref, freq_ref = result
    Q_fp, freq_fp = parse_fastprot_dump(fastprot_text, model_name)

    # Frequencies: compare directly (no scale ambiguity here).
    max_freq_diff = max(abs(a - b) for a, b in zip(freq_ref, freq_fp))

    # Rate matrix: empirical Q is only defined up to one arbitrary global
    # scale (the unit of time is arbitrary) - estimate the scale from the
    # ratio of off-diagonal sums, then check every element agrees with
    # that SAME single scale factor (not one scale per element - that
    # would hide exactly the kind of bug we're checking for).
    off_diag_ref = sum(abs(Q_ref[i][j]) for i in range(20) for j in range(20) if i != j)
    off_diag_fp = sum(abs(Q_fp[i][j]) for i in range(20) for j in range(20) if i != j)
    scale = off_diag_fp / off_diag_ref if off_diag_ref else float("nan")

    max_rel_diff = 0.0
    worst = None
    for i in range(20):
        for j in range(20):
            ref_scaled = Q_ref[i][j] * scale
            fp = Q_fp[i][j]
            denom = max(abs(ref_scaled), 1e-12)
            rel_diff = abs(fp - ref_scaled) / denom
            if rel_diff > max_rel_diff:
                max_rel_diff = rel_diff
                worst = (i, j, Q_ref[i][j], fp, scale)

    status = "OK" if max_rel_diff < 0.01 and max_freq_diff < 1e-3 else "MISMATCH"
    print(
        f"{model_name}: scale={scale:.6f}  max_freq_diff={max_freq_diff:.2e}  "
        f"max_rel_matrix_diff={max_rel_diff:.4%}  [{status}]"
    )
    if status == "MISMATCH":
        i, j, ref, fp, sc = worst
        print(f"  worst cell: ({i},{j})  iqtree*scale={ref*sc:.6f}  fastprot={fp:.6f}  (raw iqtree R*pi={ref:.6f})")


def main():
    iqtree_text = open(IQTREE_FILE).read()
    fastprot_text = open(FASTPROT_FILE).read()
    for model_name in ["WAG", "JTT", "DAYHOFF", "LG"]:
        compare_model(iqtree_text, fastprot_text, model_name)
    print("MVR: not in IQ-TREE2's model catalog - no reference available there")


if __name__ == "__main__":
    main()
