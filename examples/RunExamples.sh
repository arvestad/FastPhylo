#! /bin/bash

if [ -z $FASTPHYLOPATH ]
then
    echo "Set FASTPHYLOPATH to the directory with the project binaries!"
    exit -1
fi

export PATH=${FASTPHYLOPATH}:$PATH
#echo $PATH

function run_example {
    # Parameters:
    # $1: the number of the example
    # $2: the string with the shell command to run.
    # $3: the output file

    echo
    echo "Example "$1": "$2
    $2 > $3
    echo "Output in "$3
}


run_example 1  "fastdist -I phylip seq.phylip" ex1.out
run_example 2  "fastdist -I fasta seq.fasta" ex2.out
run_example 3  "fastdist -I xml -O xml seq.xml" ex3.out
#run_example 4  "cat seq.xml | fastdist -I xml -O xml" ex4.out

run_example 6  "fastprot -I phylip protein_seq.phylip -O xml" ex6.out
run_example 7  "fastprot -I fasta protein_seq.fasta" ex7.out
run_example 8  "fnj -r 2 -I phylip dm.phylip" ex8.out
run_example 9  "fnj -I xml dm.xml" ex9.out
# speed2026a: ID/JC/JCK/JCSS all run through the byte-encoded SIMD
# primitive as of Phase 6 (see phase0_audit.md/phase1_design.md); these
# fixtures are this project's durable regression coverage for that.
run_example 11 "fastprot -I fasta protein_seq.fasta -D ID" ex11.out
run_example 12 "fastprot -I fasta protein_seq.fasta -D JC" ex12.out
run_example 13 "fastprot -I fasta protein_seq.fasta -D JCK" ex13.out
run_example 14 "fastprot -I fasta protein_seq.fasta -D JCSS" ex14.out
run_example 15 "fastprot -I fasta globin_family.fasta -D JCK" ex15.out
#run_example 10 "cat seq.phylip | fastdist -I phylip -O phylip -b 3 -r 2 | fnj -I phylip -O xml -r 2 -d 4" ex10.out

echo
for i in ex*.out; do 
    diff -q $i expected_output/$i
done
