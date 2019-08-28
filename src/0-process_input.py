#!/usr/bin/env python3
import sys
import os

# Jonathan Rodiger - 2019

###################
## process input ##
###################

# create input file dir in tmp
# split user input by chromosome
# run pipeline on each chromosome input in parallel

user_input  = sys.argv[1]
num_threads = sys.argv[2]

chr2input = {}

os.makedirs('tmp/input')
os.makedirs('tmp/results')

with open(user_input, 'r') as f:
	header = next(f)
	for line in f:
		data = line.split(',')
		chromosome = data[1]
		if chromosome not in chr2input:
			chr2input[chromosome] = []
		chr2input[chromosome].append(line)

for chromosome in chr2input:
	with open('tmp/input/' + chromosome, 'w') as out:
		for line in chr2input[chromosome]:
			out.write(line)
