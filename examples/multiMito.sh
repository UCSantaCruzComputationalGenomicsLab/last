#! /bin/sh

# This script demonstrates using LAST and maf-join.py to construct a
# multiple alignment of the human, mouse, chicken, and fugu
# mitochondrial genomes.

# You have to run this from inside the examples directory, else it
# won't find the files.

# Put the LAST programs and scripts into the command search path:
PATH=$PATH:../src:../scripts

# Make a LAST database of the human sequence:
lastdb -c humanMito humanMito.fa

# Align the mouse sequence to the human sequence:
# (last-split will use a score threshold of 19+6=25)
lastal -e19 -j4 humanMito mouseMito.fa | last-split | maf-sort.sh > hm.maf

# Align the chicken sequence to the human sequence:
lastal -e19 -j4 humanMito chickenMito.fa | last-split | maf-sort.sh > hc.maf

# Align the fugu sequence to the human sequence:
lastal -e19 -j4 humanMito fuguMito.fa | last-split | maf-sort.sh > hf.maf

# Join the pairwise alignments into a multiple alignment:
maf-join.py hm.maf hc.maf hf.maf

# Clean up the intermediate files that we made:
rm humanMito.??? h?.maf
