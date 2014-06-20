#! /bin/sh

# Regression tests

cd $(dirname $0)

PATH=../src:$PATH

maf=SRR359290-1k.maf

{
    last-split -h

    last-split $maf

    last-split -c0 $maf

    last-split -t0 $maf

    last-split -M7.5 -S2 $maf

    last-split -m0.001 -s180 $maf

    last-split -n $maf

    last-split -d0 -m0.001 -s180 $maf
    last-split -d1 -m0.001 -s180 $maf
    last-split -d2 -m0.001 -s180 $maf

    grep -v '^q' $maf | last-split -m0.001 -s180
} |
diff last-split-test.out -
