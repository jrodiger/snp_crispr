#!/bin/bash

# arguments
species=$1
input=$2
pam=$3
threads=$4
option1=$5
option2=$6

# remove previous logs
rm error.log 2>/dev/null
rm results/no_designs.csv 2>/dev/null
rm -rf tmp/*

# check if file is in VCF format + convert
if [ ${input: -4} == ".vcf" ]; then
	python src/convert_vcf.py $input
	input="tmp/$input"
fi

# split input by chromosome
python src/0-process_input.py $input $threads

# check if using alternate reference for specificity
if [ "$option1" == "-alt"  ] || [ "$option2" == "-alt"  ]; then
	alt="-alt"
else
	alt=""
fi

# run pipeline on each tmp chromosome input file
for file in tmp/input/*; do
	output="$(basename "$file")"
	if [ "$option1" == "-both" ]; then
		echo "./src/design.sh $species $file "$output"_1 $pam -all $alt" >> tmp/commands.txt
		echo "./src/design.sh $species $file "$output"_2 $pam $alt" >> tmp/commands.txt
	else
    	echo "./src/design.sh $species $file $output $pam $option1 $alt" >> tmp/commands.txt
	fi
done

# run pipeline on each chromosome + option in parallel
parallel --jobs $threads < tmp/commands.txt

# combine results
python src/5-combine_results.py $input

# print error log + delete tmp files
cat error.log 2>/dev/null
rm -rf tmp/*
