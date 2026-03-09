# script to make one giant matrix

import sys
import argparse
import re
import glob
import subprocess
import os
from Bio import SeqIO

from collections import defaultdict

def aggregate_results(counts_file_list, output_prefix):
	# create two output files; one for raw and one for normalized
	output_raw = output_prefix + '_raw.csv'
	output_norm = output_prefix + '_norm.csv'
	o_raw = open(output_raw,'w')
	o_norm = open(output_norm,'w')

	# storage of counts / sample
	total_per_sample = defaultdict(int)

	data_table = defaultdict(dict)
	sample_list = []
	for infile in counts_file_list:
		infile = infile.rstrip('\r\n')
		ifile = open(infile,'r')
		sample_name = infile.replace('_counts_orf.csv','')
		sample_list.append(sample_name)
		for line in ifile:
			line = line.rstrip('\r\n')
			orf,bc,count = line.split(',')
			data_table[orf][sample_name] = int(count)
			if sample_name not in total_per_sample:
				total_per_sample[sample_name] = int(count)
			else:
				total_per_sample[sample_name] += int(count)
		ifile.close()
	counts_file_list.close()
	o_raw.write('ORF,' + ','.join(sample_list) + '\n')
	o_norm.write('ORF,' + ','.join(sample_list) + '\n')
	for orf in data_table:
		ltw_raw = orf
		ltw_norm = orf
		for sample in sample_list:
			if sample not in data_table[orf]:
				ltw_raw = ltw_raw + ',0'
				ltw_norm = ltw_norm + ',0'
			else:
				ltw_raw = ltw_raw + ',' + str(data_table[orf][sample])
				norm_value = int(data_table[orf][sample]) / total_per_sample[sample] * 100000
				ltw_norm = ltw_norm + ',' + str(norm_value)

		o_raw.write(ltw_raw + '\n')
		o_norm.write(ltw_norm + '\n')
	o_raw.close()
	o_norm.close()


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--counts_file_list',type=argparse.FileType('r'),required=True)	
	parser.add_argument('-o','--output_prefix',required=True)	
	opts = parser.parse_args(argv)
	aggregate_results(opts.counts_file_list, opts.output_prefix)

if __name__ == '__main__':
	main(sys.argv[1:])

