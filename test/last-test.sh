#! /bin/sh

# Exercise LAST programs, and compare the output to a reference
# output.  More tests should be added!

cd $(dirname $0)

# Make sure we use this version of LAST:
PATH=../src:../scripts:$PATH

dnaSeq=galGal3-M-32.fa
protSeq=Q2LCP8.fa
fastq=SRR001981-1k.fastq
mat=../examples/hoxd70.mat
gc=../examples/vertebrateMito.gc
seed=../examples/yass.seed
db=/tmp/last-test

{
    echo TEST 1  # spaced seeds, soft-masking, centroid alignment, matrix file
    lastdb -c -m110 $db $dnaSeq
    lastal -u1 -j5 -p $mat -x3400 -e2500 $db $dnaSeq
    echo

    echo TEST 2  # multiple volumes & query batches
    lastdb -m1 -s1 $db $dnaSeq
    lastal -f0 -i1 -w0 $db $dnaSeq
    echo

    echo TEST 3  # match-counting, with multiple query batches
    lastal -j0 -i1 -s0 $db $dnaSeq
    echo

    echo TEST 4  # FASTQ quality scores
    lastal -Q1 -e90 -a9 $db $fastq
    echo

    echo TEST 5  # translated alignment & genetic code file
    lastdb -p $db $protSeq
    lastal -F12 -pBL62 -e40 -G $gc $db $dnaSeq
    echo

    echo TEST 6  # subset seed file, soft-masking
    lastdb -c -u $seed $db $dnaSeq
    lastal -s0 -f0 -e18 $db $dnaSeq
    echo

    echo TEST 7  # asymmetric scoring matrix
    lastal -s0 -f0 -p asymmetric.mat -e2000 $db $dnaSeq
    echo

    echo TEST 8  # FASTQ-Illumina quality scores
    lastdb -m1111110 $db $dnaSeq
    lastal -Q3 -e110 $db illumina100.txt
    echo

    echo TEST 9  # PRB-format quality data
    lastal -Q4 -e90 $db mouse_tss_prb.txt
    echo

    echo TEST 10  # probabilistic alignment with quality scores
    lastal -Q1 -j6 -e90 -a9 $db $fastq
    echo

    echo TEST 11  # sparse index, generalized affine gap costs
    lastdb -w2 -c $db $dnaSeq
    lastal -r3 -q3 -a21 -c2 -e60 -f0 $db $dnaSeq
    echo

    echo TEST 12  # generalized affine gaps, frameshifts, tabular output
    lastdb -p -c $db $protSeq
    lastal -F12 -pBL62 -c2 -e40 -f0 $db $dnaSeq
    echo

    echo TEST 13  # gapless alignment, protein-protein alignment, seed freq
    lastal -j1 -f0 -e37 -m100 $db $protSeq
    echo

    echo TEST 14  # fastq-versus-fastq, seed freq
    lastdb -Q1 $db sd-ccs-100.fq
    lastal -Q1 -r1 -q2 -a1 -b1 -e44 -m100 -s0 $db sd-ccs-100.fq
    echo

    echo TEST 15  # incomplete sorting, lastal on one volume
    lastdb -i10 -s1 $db $dnaSeq
    lastal -Q1 -e90 -a9 -f0 ${db}0 $fastq
    echo

    echo TEST 16  # multiple seeds, transition constraints
    lastdb -c -m 11101T011T11,111001010010111 $db $dnaSeq
    lastal -s0 -f0 -e18 $db $dnaSeq
    echo

    echo TEST 17  # Iedera notation
    lastdb -c -m '#@#--##--#-#' $db $dnaSeq
    lastal -s0 -f0 -e18 $db $dnaSeq
    echo

    echo TEST 18  # overlap alignment, tabular output ending in gaps
    lastdb -m1111110 $db $dnaSeq
    lastal -T1 -Q1 -e60 -a9 -f0 $db $fastq
    echo

    echo TEST 19  # probabilistic overlap alignment
    lastal -T1 -Q1 -e60 -a9 -j4 $db $fastq
    echo

    echo TEST 20  # expected counts
    lastal -s0 -e18 -j7 $db $dnaSeq
    echo

    echo TEST 21  # named multi-seed, sparse query seeding
    lastdb -c -uMAM8 $db hg19-M.fa
    lastal -e34 -k128 -f0 $db galGal3-M-32.fa
    echo

    echo TEST 22  # named score matrix, sparse query seeding
    lastal -pHOXD70 -k128 -f0 $db galGal3-M-32.fa
    echo
} |
grep -v version |  # omit header lines with the LAST version number
diff last-test.out -

# Test: last-bisulfite, last-merge-batches, last-split, named seeds
lastdb -uBISF f hg19-M.fa
lastdb -uBISR r hg19-M.fa
../examples/last-bisulfite.sh f r bs100.fastq | grep -v '^#' | diff bs100.maf -
rm f.* r.*

./maf-convert-test.sh

# Test: lastdb, lastal, last-split, maf-sort, maf-join
cd ../examples
./multiMito.sh | diff multiMito.maf -

rm $db*
