#!/usr/bin/env python3
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
			data       = line.strip('\n').split(',')
			wt_crispr  = list(data[7])
			snp_crispr = list(data[8])
			for i in range(0, len(wt_crispr)):
				if wt_crispr[i] != snp_crispr[i]:
					wt_crispr[i] = wt_crispr[i].lower()
					snp_crispr[i] = snp_crispr[i].lower()
			data[7] = ''.join(wt_crispr)
			data[8] = ''.join(snp_crispr)
			snp_positions = data[5].split(';')
			last_snp = 0
			for pos in snp_positions:
				if abs(int(data[2]) - int(pos)) > last_snp:
					last_snp = abs(int(data[2]) - int(pos))
			dist_to_pam = 19 - last_snp
			summary_lines.append(','.join(data) + ',' + str(dist_to_pam) + '\n')
	with open(final_summary, 'w') as out:
		out.write(header + ',dist_to_pam\n')
		for line in summary_lines:
			out.write(line)

	# remove tmp files
	shutil.rmtree('tmp')
	os.makedirs('tmp')
