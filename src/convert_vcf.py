#!usr/bin/env python3
import sys

# vcf input file
vcf_input  = sys.argv[1]

with open('tmp/' + vcf_input, 'w') as out:
	out.write('gene_symbol,chromosome,position,strand,reference,variant,group(optional)\n')
	with open(vcf_input, 'r') as f:
		for line in f:
			if line[0] == '#':
				continue
			data = line.split('\t')
			chromosome = data[0]
			position = data[1]
			reference = data[3].upper()
			alternate = data[4].upper()

			# check chromosome id format (varies depending on file)
			if 'chr' in chromosome:
				chromosome = chromosome[3:]

			# check if multiple alternate alleles
			alternate = alternate.split(',')

			for alt in alternate:
				out.write(','.join(['',chromosome,position,'+',reference,alt,'']) + '\n')
