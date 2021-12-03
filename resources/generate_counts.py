# Author: Raj Chari
# Date: May 2nd, 2020
# Purpose: script to generate a counts file that can be used for subsequent analysis
# Input: sam alignment file, reference fasta file
# Output: strict counts file, mapped alignment count file, bar plot file, mapping summary file

import sys
import argparse
import re
import os
import glob
import subprocess
from Bio import SeqIO

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

from collections import defaultdict

def generate_library_output(sam_file, reference_file, strict_count_file, sra_count_file, bar_plot_file_strict, bar_plot_file_sra, mapping_summary_file):
	# global variables
	reference_db = defaultdict(str)
	spacer_db = defaultdict(str)
	strict_counts_db = defaultdict(int)
	sra_counts_db = defaultdict(int)
	total_reads = 0
	total_strict = 0
	total_sra = 0

	# go through reference file
	for record in SeqIO.parse(reference_file,'fasta'):
		name = str(record.id)
		# find first CACCG
		sequence = str(record.seq)
		index = sequence.find('CACCG')
		spacer = sequence[index+5:index+25]
		reference_db[name] = spacer
		if spacer not in spacer_db:
			spacer_db[spacer] = name
		else:
			print('Duplicate spacer!' + spacer)
	reference_file.close()

	# go through SAM file
	for line in sam_file:
		line = line.rstrip('\r\n')
		if line.startswith('@')==False:
			# add to total
			total_reads += 1

			# tally for sra
			parts = line.split('\t')
			if parts[2] != '*':
				cigar = parts[5]
				if cigar[2]=='M' and int(cigar[0:2]) >= 20:
					total_sra += 1
					if parts[2] not in sra_counts_db:
						sra_counts_db[parts[2]] = 1
					else:
						sra_counts_db[parts[2]] += 1

			# tally for strict
			seq = parts[9]
			spacer = seq[0:20]
			notindb = 0

			# check in spacer db
			if spacer in spacer_db:
				total_strict += 1
				refname = spacer_db[spacer]
				#print(ref)
				if refname not in strict_counts_db:
					strict_counts_db[refname] = 1
				else:
					strict_counts_db[refname] += 1
			else:
				notindb += 1

	# write the output to both
	num_zero_strict = 0
	num_zero_sra = 0
	for spacer_id in reference_db:
		# break spacer ID
		gene_symbol = 'N/A'
		if '-' in spacer_id:
			spacer_id_parts = spacer_id.split('-')
			gene_symbol = spacer_id_parts[0]

		# grab counts
		strict_count = 0
		sra_count = 0
		if spacer_id in strict_counts_db:
			strict_count = strict_counts_db[spacer_id]
		else:
			strict_counts_db[spacer_id] = 0

		if spacer_id in sra_counts_db:
			sra_count = sra_counts_db[spacer_id]
		else:
			sra_counts_db[spacer_id] = 0

		# write to file
		strict_count_file.write(spacer_id + ',' + gene_symbol + ',' + str(strict_count) + '\n')
		sra_count_file.write(spacer_id + ',' + gene_symbol + ',' + str(sra_count) + '\n')

		# add to zero count
		if strict_count==0:
			num_zero_strict += 1
		if sra_count==0:
			num_zero_sra += 1

	# total spacers
	total_spacers = len(reference_db.keys())

	# write mapping summary
	mapping_summary_file.write('Sample,Total_Reads,Mapped_Strict,Mapped_SRA\n')
	sample_name = os.path.basename(sam_file.name)
	sample_name = sample_name.replace('_bwamem.sam','')
	mapping_summary_file.write(sample_name + ',' + str(total_reads) + ',' + str(total_strict) + ',' + str(total_sra) + '\n')

	# plot the graphs
	plot_data(strict_counts_db, bar_plot_file_strict, sample_name, num_zero_strict, total_reads, total_strict, 'strict')
	plot_data(sra_counts_db, bar_plot_file_sra, sample_name, num_zero_sra, total_reads, total_sra, 'sra')

	# close file handles
	sam_file.close()
	strict_count_file.close()
	sra_count_file.close()
	mapping_summary_file.close()

def plot_data(counts_db, filename, sample, zeros, total_reads, total_mapped, graph_type):

	# x axis is the number of possible values
	x_locs = range(len(counts_db.keys()))

	# plot bar plot
	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

	# generate the CDF data
	sorted_counts = sorted(list(counts_db.values()),reverse=True)
	cdf_data = []
	rolling_cdf = 0.00
	for value in sorted_counts:
		proportion = value / total_reads
		rolling_cdf = rolling_cdf + proportion
		cdf_data.append(rolling_cdf)

	# percent mapped
	percent_mapped = (total_mapped / total_reads) * 100
	rounded = round(percent_mapped,2)

	# set the title
	fig.suptitle(sample + '- percent mapped - ' + graph_type + ': ' + str(rounded) + ', num zero: ' + str(zeros) + ' out of ' + str(len(counts_db.keys())) + ' guides', fontsize=10)

	ax1.plot(x_locs,cdf_data,color='b')
	ax1.set_ylabel('Cumulative proportion of total')
	ax2.set_xlabel('Guide sequence rank ordered by decreasing abundance')

	ax2.bar(x_locs,sorted_counts,color='b')
	ax2.set_ylabel('Number of reads')
	plt.tight_layout()

	# save the file
	plt.savefig(filename.name)	

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-a','--sam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-r','--reference_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-s','--strict_count_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-b','--sra_count_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-g','--bar_plot_file_strict',type=argparse.FileType('w'),required=True)
	parser.add_argument('-f','--bar_plot_file_sra',type=argparse.FileType('w'),required=True)
	parser.add_argument('-m','--mapping_summary_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	generate_library_output(opts.sam_file, opts.reference_file, opts.strict_count_file, opts.sra_count_file, opts.bar_plot_file_strict, opts.bar_plot_file_sra, opts.mapping_summary_file)
 
if __name__ == '__main__':
	main(sys.argv[1:])