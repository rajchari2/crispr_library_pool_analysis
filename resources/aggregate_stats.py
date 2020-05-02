# Author: Raj Chari
# Date: May 2nd, 2020
# Purpose: to aggregate all of the mapping statistics
# Input: sam alignment file, reference fasta file
# Output: strict counts file, mapped alignment count file, bar plot file, mapping summary file

import sys
import argparse
import re
import glob
import subprocess
import os
from Bio import SeqIO

from collections import defaultdict


def aggregate_results(mapping_file_list, strict_file_list, sra_file_list, mapping_stats_file, strict_matrix_file, sra_matrix_file):
	# go through list of mapping stats
	header_written = False
	for summary_file in mapping_file_list:
		infile = open(summary_file,'r')
		line_count = 0
		for line in infile:
			if line_count > 0:
				mapping_stats_file.write(line)
			elif header_written==False:
				mapping_stats_file.write(line)
				header_written = True
			line_count += 1
		infile.close()
	mapping_stats_file.close()

	# now go through and write the matrix for counts for strict
	sample_list_strict = []
	guide_list = defaultdict(str)
	strict_db = defaultdict(dict)
	for strict_file in strict_file_list:
		sample_name = os.path.basename(strict_file)
		sample_name = sample_name.replace('_strict_match_counts.csv','')
		sample_list_strict.append(sample_name)
		infile = open(strict_file,'r')
		for line in infile:
			line = line.rstrip('\r\n')
			parts = line.split(',')
			guide_key = parts[0] + ',' + parts[1]
			strict_db[sample_name][guide_key] = parts[2]
			guide_list[guide_key] = 'Y'
		infile.close()

	# write to file
	strict_matrix_file.write('Guide-name,Gene,' + ','.join(sample_list_strict) + '\n')
	for guide in guide_list:
		ltw = guide
		for sample in sample_list_strict:
			ltw = ltw + ',' + strict_db[sample][guide]
		strict_matrix_file.write(ltw + '\n')

	# now go through sra
	sample_list_sra = []
	sra_db = defaultdict(dict)
	for sra_file in sra_file_list:
		sample_name = os.path.basename(sra_file)
		sample_name = sample_name.replace('_sra_counts.csv','')
		sample_list_sra.append(sample_name)
		infile = open(sra_file,'r')
		for line in infile:
			line = line.rstrip('\r\n')
			parts = line.split(',')
			guide_key = parts[0] + ',' + parts[1]
			sra_db[sample_name][guide_key] = parts[2]
		infile.close()

	# write to file
	sra_matrix_file.write('Guide-name,Gene,' + ','.join(sample_list_sra) + '\n')
	for guide in guide_list:
		ltw = guide
		for sample in sample_list_sra:
			ltw = ltw + ',' + sra_db[sample][guide]
		sra_matrix_file.write(ltw + '\n')

	# close file handles
	mapping_stats_file.close()
	strict_matrix_file.close()
	sra_matrix_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--mapping_file_list',nargs='+',required=True)
	parser.add_argument('-a','--strict_file_list',nargs='+',required=True)
	parser.add_argument('-b','--sra_file_list',nargs='+',required=True)
	parser.add_argument('-m','--mapping_stats_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-s','--strict_matrix_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-c','--sra_matrix_file',type=argparse.FileType('w'),required=True)	
	opts = parser.parse_args(argv)
	aggregate_results(opts.mapping_file_list, opts.strict_file_list, opts.sra_file_list, opts.mapping_stats_file, opts.strict_matrix_file, opts.sra_matrix_file)

if __name__ == '__main__':
	main(sys.argv[1:])