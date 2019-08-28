#!/bin/bash

# arguments
species=$1
input=$2
pam=$3
threads=$4
option=$5

# remove previous logs
rm error.log 2>/dev/null
rm results/no_designs.csv 2>/dev/null

# split input by chromosome
python src/0-process_input.py $input $threads

# run pipeline on each tmp chromosome input file
for file in tmp/input/*; do
	output="$(basename "$file")"
	if [ "$option" == "-both" ]; then
		echo "./src/design.sh $species $file "$output"_1 $pam -all" >> tmp/commands.txt
		echo "./src/design.sh $species $file "$output"_2 $pam" >> tmp/commands.txt
	else
    	echo "./src/design.sh $species $file $output $pam $option" >> tmp/commands.txt
	fi
done

# run pipeline on each chromosome + option in parallel
parallel --jobs $threads < tmp/commands.txt

# combine results
python src/5-combine_results.py $input

# print error log + delete tmp files
cat error.log 2>/dev/null
rm -rf tmp/*
