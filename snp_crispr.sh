#!/bin/bash

# arguments
species=$1
input=$2
pam=$3
option1=$4
option2=$5
outputfilename="designs"

# parallel
using_parallel=false

# check for both argument
if [ "$option1" == "-both" ]; then
	rm error.log 2>/dev/null
	echo "./snp_crispr.sh $species $input $pam 1" > tmp/commands.txt
	echo "./snp_crispr.sh $species $input $pam -all 2" >> tmp/commands.txt
	parallel < tmp/commands.txt
	python src/5-combine_results.py
	cat error.log 2>/dev/null
	rm -f tmp/*
	exit 0
fi

# check for all argument + parallel
if [ "$option1" == "1" ]; then
	using_parallel=true
	outputfilename="part1"
	all=""
elif [ "$option1" == "-all" ]; then
	all="-all"
	if [ "$option2" == "2" ]; then
		using_parallel=true
		outputfilename="part2"
	fi
else
	all=""
fi

# paths
wt="tmp/$outputfilename-wt.fasta"
snp="tmp/$outputfilename-snp.fasta"
blast_wt="tmp/$outputfilename-wt_blast.txt"
blast_snp="tmp/$outputfilename-snp_blast.txt"

# don't delete previous error log twice in parallel
if [ "$using_parallel" == false ]; then
	rm error.log 2>/dev/null
fi

# run pipeline + check error log
python src/1-find_crispr_designs.py $species $input tmp/$outputfilename $pam $all
perl src/2-blast_crisprs.pl -s $species -f1 $wt -f2 $snp -o tmp/$outputfilename -pam $pam
perl src/3-calculate_scores.pl -b1 $blast_wt -b2 $blast_snp -out $outputfilename
python src/4-process_results.py $blast_wt $blast_snp results/$outputfilename

# don't print error log / delete tmp files twice when running in parallel
if [ "$using_parallel" == false ]; then
	cat error.log 2>/dev/null
	rm -f tmp/*
fi
