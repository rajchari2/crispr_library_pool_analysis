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
from Bio import SeqIO

from collections import defaultdict


def aggregate_results(summary_file_list, mapping_stats_file):
	# go through input file
	header_written = False
	for summary_file in summary_file_list:
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


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--summary_file_list',nargs='+',required=True)
	parser.add_argument('-o','--mapping_stats_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	aggregate_results(opts.summary_file_list, opts.mapping_stats_file)

if __name__ == '__main__':
	main(sys.argv[1:])