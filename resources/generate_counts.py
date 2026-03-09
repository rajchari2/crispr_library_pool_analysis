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
from Bio.Seq import Seq

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

from collections import defaultdict

def generate_library_output(input_file, reference_file, library_type, strict_count_file, bar_plot_file_strict, mapping_summary_file):
	# if the library type is single guide or ORF, the reference file is a fasta file
	sample_name = ''
	# deal with each type of library separately
	if library_type=='single_guide':
		# global variables
		reference_db = defaultdict(str)
		spacer_db = defaultdict(str)
		strict_counts_db = defaultdict(int)
		total_reads = 0
		total_strict = 0		

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
		for line in input_file:
			line = line.rstrip('\r\n')
			if line.startswith('@')==False:
				# split
				parts = line.split('\t')

				# add to total
				total_reads += 1

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
		for spacer_id in reference_db:
			# break spacer ID
			gene_symbol = 'N/A'
			if '-' in spacer_id:
				spacer_id_parts = spacer_id.split('-')
				gene_symbol = spacer_id_parts[0]

			# grab counts
			strict_count = 0
			if spacer_id in strict_counts_db:
				strict_count = strict_counts_db[spacer_id]
			else:
				strict_counts_db[spacer_id] = 0

			# write to file
			strict_count_file.write(spacer_id + ',' + gene_symbol + ',' + str(strict_count) + '\n')

			# add to zero count
			if strict_count==0:
				num_zero_strict += 1

		# total spacers
		total_spacers = len(reference_db.keys())

		# write mapping summary
		mapping_summary_file.write('Sample,Total_Reads,Mapped_Strict\n')
		sample_name = os.path.basename(input_file.name)
		sample_name = sample_name.replace('_bwamem.sam','')
		mapping_summary_file.write(sample_name + ',' + str(total_reads) + ',' + str(total_strict) + '\n')

		# plot the graphs
		plot_data(strict_counts_db, bar_plot_file_strict, sample_name, num_zero_strict, total_reads, total_strict, 'strict')

		# close file handles
		input_file.close()
		strict_count_file.close()
		mapping_summary_file.close()

	elif library_type == 'orf':
		reference_db = defaultdict(str)
		strict_counts_db = defaultdict(int)
		total_reads = 0
		total_strict = 0
		quality_count = 0
		barcode_length = 0
		num_zero_strict = 0
		# store the barcodes for each ORF
		for record in SeqIO.parse(reference_file,'fasta'):
			name = str(record.id)
			sequence = str(record.seq)
			barcode_length = len(sequence)
			reference_db[sequence] = name
		reference_file.close()

		# data file will be trimmed fastq where we want the first n bases of sequence
		for record in SeqIO.parse(input_file,'fastq'):
			total_reads += 1
			barcode = str(record.seq)[0:barcode_length]
			scores = record.letter_annotations["phred_quality"]
			passed = True
			# go through and make sure the score is > 20
			i = 0
			while i < barcode_length:
				if int(scores[i]) < 20:
					passed = False
					break
				i += 1
			if passed==True:
				quality_count += 1
				if barcode in reference_db:
					total_strict += 1
					orf = reference_db[barcode]
					if orf not in strict_counts_db:
						strict_counts_db[orf] = 1
					else:
						strict_counts_db[orf] += 1
		

		# go through and write
		for barcode in reference_db:
			orf_name = reference_db[barcode]
			if orf_name in strict_counts_db:
				count = strict_counts_db[orf_name]
			else:
				count = 0
			strict_count_file.write(orf_name + ',' + barcode + ',' + str(count) + '\n')
			# add to zero count
			if count==0:
				num_zero_strict += 1
				
		# write mapping summary
		mapping_summary_file.write('Sample,Total_Reads,Quality_Reads,Mapped_Strict\n')
		sample_name = os.path.basename(input_file.name)
		sample_name = sample_name.replace('_R1_trimmed.fastq','')
		mapping_summary_file.write(sample_name + ',' + str(total_reads) + ',' + str(quality_count) + ',' + str(total_strict) + '\n')

		# plot the graphs
		plot_data(strict_counts_db, bar_plot_file_strict, sample_name, num_zero_strict, quality_count, total_strict, 'strict')

		# close file handle
		input_file.close()
		strict_count_file.close()
		mapping_summary_file.close()

	elif library_type=='dual_guide':
		# global variables
		reference_db = defaultdict(str)
		# counts
		strict_counts_db = defaultdict(int)
		r1_counts_db = defaultdict(int)
		r2_counts_db = defaultdict(int)
		# read matches
		read_db = defaultdict(list)
		r1_db = defaultdict(list)
		r2_db = defaultdict(list)
		total_reads = 0
		total_passed_reads = 0
		r1_r2_matchcount = 0
		r1_match_count = 0
		r2_match_count = 0

		# reference R2 file
		r2_filename = input_file.name.replace('R1','R2')
		r2_handle = open(r2_filename,'r')

		# go through reference file
		for line in reference_file:
			line = line.rstrip('\r\n')
			g_name,sequence = line.split('\t')
			guide1 = sequence[12:32]
			guide2 = sequence[78:98]
			index = guide1 + '_' + guide2
			reference_db[index] = g_name
			# keep track of which guides belong to which pair
			r1_db[guide1].append(g_name)
			r2_db[guide2].append(g_name)
		reference_file.close()

		# go through and store sequence for each read, first do R1
		for record in SeqIO.parse(input_file,'fastq'):
			fastq_id = str(record.id)
			sequence = str(record.seq)
			guide = sequence[0:20]
			if fastq_id not in read_db:
				read_db[fastq_id] = []
			read_db[fastq_id].append(guide)
		input_file.close()

		total_reads = len(read_db.keys())

		# go through R2
		for record in SeqIO.parse(r2_handle,'fastq'):
			fastq_id = str(record.id)
			sequence = str(record.seq)
			guide_RC = sequence[0:20]
			guide = str(Seq(guide_RC).reverse_complement())
			if fastq_id in read_db:
				read_db[fastq_id].append(guide)
			else:
				total_reads += 1
		r2_handle.close()

		# go through read db and try to match
		for read in read_db:
			if len(read_db[read])==2:
				total_passed_reads += 1
				g1 = read_db[read][0]
				g2 = read_db[read][1]
				test_index = g1 + '_' + g2
				# first check for dual match
				if test_index in reference_db:
					r1_r2_matchcount += 1
					g_name = reference_db[test_index]
					if g_name not in strict_counts_db:
						strict_counts_db[g_name] = 1
					else:
						strict_counts_db[g_name] += 1
				# g1 matches in read 1 db and it's unique
				elif g1 in r1_db and len(r1_db[g1])==1:
					r1_match_count += 1
					# get the appropriate ID
					g_id = r1_db[g1][0]
					if g_id not in r1_counts_db:
						r1_counts_db[g_id] = 1
					else:
						r1_counts_db[g_id] += 1
				elif g2 in r2_db and len(r2_db[g2])==1:
					r2_match_count += 1
					# get the appropriate ID
					g_id = r2_db[g2][0]
					if g_id not in r2_counts_db:
						r2_counts_db[g_id] = 1
					else:
						r2_counts_db[g_id] += 1				

		# write mapping summary
		sample_name = os.path.basename(input_file.name)
		mapping_summary_file.write('Sample,Total_Reads,Total_Passed_Reads,MatchBoth,MatchG1,MatchG2\n')
		mapping_summary_file.write(sample_name + ',' + str(total_reads) + ',' + str(total_passed_reads) + ',' + str(r1_r2_matchcount) + ',' + str(r1_match_count) + ',' + str(r2_match_count) + '\n')
		mapping_summary_file.close()

		# write output of counts
		strict_count_file.write('Guide_name,Spacer1,Spacer2,Both,R1,R2\n')
		for index in reference_db:
			g_name = reference_db[index]
			# counts
			r1_r2_count = 0
			r1_count = 0
			r2_count = 0
			if g_name in strict_counts_db:
				r1_r2_count = strict_counts_db[g_name]
			if g_name in r1_counts_db:
				r1_count = r1_counts_db[g_name]
			if g_name in r2_counts_db:
				r2_count = r2_counts_db[g_name]
			spacer1,spacer2 = index.split('_')
			strict_count_file.write(g_name + ',' + spacer1 + ',' + spacer2 + ',' + str(r1_r2_count) + ',' + str(r1_count) + ',' + str(r2_count) + '\n')



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
	fig.suptitle(sample + '- percent mapped - ' + graph_type + ': ' + str(rounded) + ', num zero: ' + str(zeros) + ' out of ' + str(len(counts_db.keys())) + ' sequences', fontsize=10)

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
	parser.add_argument('-i','--input_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-r','--reference_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-s','--strict_count_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-g','--bar_plot_file_strict',type=argparse.FileType('w'),required=True)
	parser.add_argument('-m','--mapping_summary_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-l','--library_type',required=True)
	opts = parser.parse_args(argv)
	generate_library_output(opts.input_file, opts.reference_file, opts.library_type, opts.strict_count_file, opts.bar_plot_file_strict, opts.mapping_summary_file)
 
if __name__ == '__main__':
	main(sys.argv[1:])