#! /bin/bash
#
# Differential test for the speed2026a protein-distance speedup
# (see plan.md, phase1_design.md): runs fastprot's ID/JC/JCK/JCSS models
# (the ones wired to the FASTPROT_NEW_PROT_CODE switch in
# ProtDistCalc.cpp - ED/ML are untouched, see phase0_audit.md) against
# each of its old and new count_id_dist() primitives on the same inputs,
# and fails loudly if the outputs differ. Extends RunExamples.sh's
# diff-against-a-known-good-output convention to an old-vs-new
# differential comparison instead of a fixed fixture.
#
# Datasets:
#   protein_seq.fasta   - existing tiny fixture (4 seqs, ~12 residues),
#                          used by RunExamples.sh examples 6/7.
#   globin_family.fasta - real biological data: 25 reviewed UniProt
#                          globin sequences, aligned with muscle 5.2
#                          (see phase0_audit.md / this script's git log
#                          for provenance), 167 columns.

if [ -z $FASTPHYLOPATH ]
then
    echo "Set FASTPHYLOPATH to the directory with the project binaries!"
    exit 1
fi

export PATH=${FASTPHYLOPATH}:$PATH

FAIL=0

function diff_check {
    # $1: description  $2: input format (fasta/phylip)  $3: input file  $4: model
    local desc="$1" fmt="$2" infile="$3" model="$4"
    local old_out=$(mktemp)
    local new_out=$(mktemp)

    fastprot -I "$fmt" "$infile" -D "$model" -O phylip -o "$old_out"
    FASTPROT_NEW_PROT_CODE=1 fastprot -I "$fmt" "$infile" -D "$model" -O phylip -o "$new_out"

    if diff -q "$old_out" "$new_out" > /dev/null; then
        echo "PASS: $desc, model=$model"
    else
        echo "FAIL: $desc, model=$model - old and new primitives disagree"
        diff "$old_out" "$new_out" | head -5
        FAIL=1
    fi
    rm -f "$old_out" "$new_out"
}

for model in ID JC JCK JCSS; do
    diff_check "tiny fixture (protein_seq.fasta)" fasta protein_seq.fasta "$model"
    diff_check "real globin family (globin_family.fasta)" fasta globin_family.fasta "$model"
done

echo
if [ $FAIL -eq 0 ]; then
    echo "RunProtDiff.sh: all old-vs-new comparisons matched"
else
    echo "RunProtDiff.sh: DISAGREEMENTS FOUND (see FAIL lines above)"
fi
exit $FAIL
