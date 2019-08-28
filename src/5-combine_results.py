#!/usr/bin/env python3
from Bio.Seq import Seq
import sys
import os

# Jonathan Rodiger - 2019

# arguments
input_file = sys.argv[1]

# dictionary forces result keys to be unique
results = {}
successful_designs = {}

# loop through all results files
for filename in os.listdir('tmp/results'):
	with open('tmp/results/' + filename, 'r') as f:
		header = next(f)
		for line in f:
			results[line] = None

# write combined results
with open('results/designs.csv', 'w') as out:
	out.write(header)
	for line in results:
		out.write(line)
		data = line.split(',')
		chromosome = data[1]
		start      = data[2]
		stop       = data[3]
		strand     = data[4]
		position   = data[5]
		# store succesfully designed variants
		for variant in data[6].split(';'):
			ref, variant = variant.split('>')
			successful_designs[','.join([chromosome, position, strand, ref, variant])] = None

# check for input variants w/ no designs
failed_designs = False
with open('results/no_designs.csv', 'a') as out:
	with open(input_file, 'r') as f:
		out.write('chromosome,position,strand,reference,variant\n')
		next(f)
		for line in f:
			data = line.split(',')
			chromosome = data[1]
			position   = data[2]
			strand     = data[3]
			ref        = data[4]
			variant    = data[5]
			input_row  = ','.join([chromosome, position, strand, ref, variant])
			# check if design is on opposite strand of input
			if strand == '+':
				alt_strand = '-'
			else:
				alt_strand = '+'
			alt_ref = str(Seq(ref).reverse_complement())
			alt_variant = str(Seq(variant).reverse_complement())
			alt_input_row = ','.join([chromosome, position, alt_strand, alt_ref, alt_variant])
			if input_row not in successful_designs and alt_input_row not in successful_designs:
				failed_designs = True
				out.write(input_row + '\n')
if not failed_designs:
	os.remove('results/no_designs.csv')
