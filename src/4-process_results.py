#!/usr/bin/env python3
from Bio.Seq import Seq
import shutil
import sys
import os

# Jonathan Rodiger - 2019

# Check for off-targets in blast report w/ only mismatch in N of PAM (-NGG/-NAG)
def check_designs(filename):
	with open(filename, 'r') as f:
		for line in f:
			data = line.split('\t')
			if data[10] == '0' and data[11] == '0' and data[9][20] == 'X':
				global bad_designs
				global bad_design_list
				bad_designs = True
				bad_design_list.append(data[5])


# Returns designs w/ snp positions as lowercase
def lowercase_snps(wt_crispr, variant_crispr):
	wt_crispr = list(wt_crispr)
	variant_crispr = list(variant_crispr)
	# check sequences for mismatch/snp + make lowercase
	for i in range(0, len(wt_crispr)):
		if wt_crispr[i] != variant_crispr[i]:
			wt_crispr[i] = wt_crispr[i].lower()
			variant_crispr[i] = variant_crispr[i].lower()
	wt_crispr = ''.join(wt_crispr)
	variant_crispr = ''.join(variant_crispr)
	return wt_crispr, variant_crispr


# Returns designs w/ indel variant positions as lowercase
def lowercase_indels(ref, variant, start, stop, position, strand, wt_crispr, variant_crispr):
	# only convert to int for indels because snps can have multiple
	# positions for the same design seperated by semicolons
	start    = int(start)
	stop     = int(stop)
	position = int(position)
	# modify reverse complement, then switch back when done
	if strand == '-':
		tmp = start
		start = stop
		stop = tmp
		wt_crispr = str(Seq(wt_crispr).reverse_complement())
		variant_crispr = str(Seq(variant_crispr).reverse_complement())
	wt_crispr = list(wt_crispr)
	variant_crispr = list(variant_crispr)
	# wt crispr
	if len(ref) == 1:
		# check if reference position is within design
		if position - start >= 0:
			wt_crispr[position-start] = wt_crispr[position-start].lower()
	else:
		ref_start = position - start
		ref_stop = ref_start + len(ref)
		# check if range outside of design seq
		if ref_start < 0:
			ref_start = 0
		if ref_stop > len(wt_crispr):
			ref_stop = len(wt_crispr)
		# lowercase reference range
		for i in range(ref_start, ref_stop):
			wt_crispr[i] = wt_crispr[i].lower()
	# variant crispr
	variant_start = position - start
	variant_stop = variant_start + len(variant)
	# check if range outside of design seq
	if variant_start < 0:
		variant_start = 0
	if variant_stop > len(variant_crispr):
		variant_stop = len(variant_crispr)
	# lowercase variant range
	for i in range(variant_start, variant_stop):
		variant_crispr[i] = variant_crispr[i].lower()
	wt_crispr = ''.join(wt_crispr)
	variant_crispr = ''.join(variant_crispr)
	# take reverse complement to revert to original design
	if strand == '-':
		wt_crispr = str(Seq(wt_crispr).reverse_complement())
		variant_crispr = str(Seq(variant_crispr).reverse_complement())

	return wt_crispr, variant_crispr


# Returns true if indel doesn't overlap w/ 20mer
def indel_outside_20mer(variant_crispr):
	for i in range(0,20):
		if variant_crispr[i].islower():
			return False
	return True


# Returns distance from pam to closest variant
def distance_to_pam(variant_crispr):
	last_snp = 0
	for i in range(0, len(variant_crispr)):
		if variant_crispr[i].islower():
			last_snp = i
	dist_to_pam = 19 - last_snp
	return dist_to_pam


if __name__ == '__main__':
	wt_blast      = sys.argv[1]
	snp_blast     = sys.argv[2]
	final_summary = sys.argv[3] + '.csv'

	# check if designs were found
	if not os.path.exists(final_summary):
		print('No designs found')
		exit()

	# globals
	bad_designs       = False
	bad_design_list   = []
	summary_lines     = []
	bad_summary_lines = []

	# check for severe off targets
	check_designs(wt_blast)
	check_designs(snp_blast)

	# Remove bad designs if found
	if bad_designs == True:
		with open(final_summary, 'r') as f:
			for line in f:
				summary_lines.append(line)
				data = line.split(',')
				if data[7] in bad_design_list or data[8] in bad_design_list:
					bad_summary_lines.append(line)
		with open(final_summary, 'w') as out:
			for line in summary_lines:
				if line not in bad_summary_lines:
					out.write(line)

	summary_lines = []
	# Show variant positions as lowercase + annotate distance to PAM from variant
	with open(final_summary, 'r') as f:
		header = next(f).strip('\n')
		for line in f:
			data           = line.strip('\n').split(',')
			start          = data[2]
			stop           = data[3]
			strand         = data[4]
			position       = data[5]
			# only store first snp if multiple, used for indels
			ref, variant   = data[6].split(';')[0].split('>')
			wt_crispr      = data[7]
			variant_crispr = data[8]
			# snps
			if len(ref) == 1 and len(variant) == 1:
				wt_crispr, variant_crispr = lowercase_snps(wt_crispr, variant_crispr)
			# indels
			else:
				wt_crispr, variant_crispr = lowercase_indels(ref, variant, start, stop, 
												position, strand, wt_crispr, variant_crispr)
				# skip designs that only overlap w/ PAM
				if indel_outside_20mer(variant_crispr):
					continue
			data[7] = wt_crispr
			data[8] = variant_crispr
			dist_to_pam = distance_to_pam(variant_crispr)
			summary_lines.append(','.join(data) + ',' + str(dist_to_pam) + '\n')
	with open(final_summary, 'w') as out:
		out.write(header + ',dist_to_pam\n')
		for line in summary_lines:
			out.write(line)

	# remove tmp files
	shutil.rmtree('tmp')
	os.makedirs('tmp')
