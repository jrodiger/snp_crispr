#!/bin/bash

# arguments
species=$1
input=$2
output=$3
pam=$4
option1=$5
option2=$6

# check if using alternate reference for specificity
if [ "$option1" == "-alt"  ] || [ "$option2" == "-alt"  ]; then
	alt="true"
else
	alt="false"
fi

# check if considering all variants with sgRNA sequence
if [ "$option1" == "-all"  ]; then
	all="-all"
else
	all=""
fi

# paths
wt="tmp/$output-wt.fasta"
snp="tmp/$output-snp.fasta"
blast_wt="tmp/$output-wt_blast.txt"
blast_snp="tmp/$output-snp_blast.txt"

# pipeline
python src/1-find_crispr_designs.py $species $input tmp/$output $pam $all
perl src/2-blast_crisprs.pl -s $species -f1 $wt -f2 $snp -o tmp/$output -pam $pam -alt $alt
perl src/3-calculate_scores.pl -b1 $blast_wt -b2 $blast_snp -out $output
python src/4-process_results.py $blast_wt $blast_snp $input tmp/results/$output
