#!/usr/bin/env python3
import os

# Jonathan Rodiger - 2019

# dictionary forces result keys to be unique
results = {}

with open('results/part1.csv', 'r') as f:
	header = next(f)
	for line in f:
		results[line] = None

with open('results/part2.csv', 'r') as f:
	next(f)
	for line in f:
		results[line] = None

with open('results/designs.csv', 'w') as out:
	out.write(header)
	for line in results:
		out.write(line)

os.remove('results/part1.csv')
os.remove('results/part2.csv')
