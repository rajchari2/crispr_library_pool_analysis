# script to make a yaml file from a tab-delimited sample file for analysis of CRISPR library data
# pipeline produces two counts files -> one based on strict string mapping and one based on short read alignment
# updated to analyze dual guide libraries and barcoded ORF libraries



# input file: [tab-delimited - sample name {tab} reference file {tab} adapter
# output file: yaml file which can run in the pipeline

# import appropriate libraries
import yaml
import sys
import argparse
import re
import glob
import subprocess
from collections import defaultdict

def makeYAML (input_file,ngs_run,analysis_name,output_file,directory_name):
	# dictionary to hold yaml
	yaml_dict = {}
	yaml_dict['samples'] = []
	yaml_dict['ngs_run'] = ngs_run
	yaml_dict['analysis_name'] = analysis_name

	for line in input_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')

		# variables to store fields
		sample_name = parts[0]
		reference_file = parts[1]
		five_prime_adapter_R1 = parts[2]
		five_prime_adapter_R2 = parts[3]

		# go through and change the filenames
		path_to_check = directory_name + '/' + sample_name + '_*.fastq.gz'
		file_list = sorted(glob.glob(path_to_check))

		# and change each file
		mvCommand_R1 = 'mv ' + file_list[0] + ' ' + directory_name + '/' + sample_name + '_R1.fastq.gz'
		p = subprocess.Popen(mvCommand_R1,shell=True)
		p.communicate()
		if five_prime_adapter_R2 != 'N/A':
			mvCommand_R2 = 'mv ' + file_list[1] + ' ' + directory_name + '/' + sample_name + '_R2.fastq.gz'
			p = subprocess.Popen(mvCommand_R2,shell=True)
			p.communicate()
		# add to dictionary
		sample_dict = {}
		sample_dict[sample_name] = {}
		sample_dict[sample_name]['name'] = sample_name
		sample_dict[sample_name]['reference'] = reference_file
		sample_dict[sample_name]['five_prime_adapter_R1'] = five_prime_adapter_R1
		sample_dict[sample_name]['five_prime_adapter_R2'] = five_prime_adapter_R2

		# add to yaml dict
		yaml_dict['samples'].append(sample_dict)

	# write final yaml file
	with open(output_file, 'w') as yaml_file:
		yaml.dump(yaml_dict, yaml_file, default_flow_style=False)

    # close everything
	input_file.close()
	yaml_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-n','--ngs_run',required=True)
	parser.add_argument('-a','--analysis_name',required=True)
	parser.add_argument('-d','--data_directory',required=True)
	parser.add_argument('-o','--output_file',required=True)
	opts = parser.parse_args(argv)
	makeYAML(opts.input_file,opts.ngs_run,opts.analysis_name,opts.output_file,opts.data_directory)
 
if __name__ == '__main__':
	main(sys.argv[1:])