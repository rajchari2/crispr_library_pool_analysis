# test

import snakemake
import os
from collections import defaultdict

# get all of the samples
samples = config["samples"]
ngs_run = config["ngs_run"]
sample_list = []
reference_list = []
sample_to_reference = defaultdict(str)
sample_to_adapter = defaultdict(str)

# go through samples
for sample in samples:
	for sample_name in sample:
		sample_list.append(sample_name)
		ref_file = sample[sample_name]['reference']
		adapter = sample[sample_name]['adapter']
		if ref_file not in reference_list:
			reference_list.append(ref_file)
		sample_to_reference[sample_name] = ref_file
		sample_to_adapter[sample_name] = adapter

# functions to get variables
def get_reference(wildcards):
	return sample_to_reference[wildcards.sample]

def get_adapter(wildcards):
	return sample_to_adapter[wildcards.sample]

rule bwa_index:
	input:
		reference_file = '{reference}'
	output:
		bai = '{reference}.bwt'
	run:
		if os.path.exists(output.bai)==False:
			shell('bwa index {input.reference_file}')

rule build_indexes:
	input:
		ref_list = expand(rules.bwa_index.output.bai,reference=reference_list)
	output:
		index_file_created = 'processed/' + ngs_run + '_Indexes_Created.tab'
	run:
		with open(output.index_file_created, "w") as out:
			out.write('References indexed\n')

rule cutadapter:
	input:
		fastq = 'data_files/{sample}_R1.fastq.gz'
	params:
		adapter = get_adapter
	output:
		fastq_trimmed = 'data_files/{sample}_R1_trimmed.fastq'
	shell:
		'cutadapt -g {params.adapter} --discard-untrimmed --quiet -o {output.fastq_trimmed} {input.fastq}'


rule bwa_mem:
	input:
		reference_file = get_reference,
		fastq = rules.cutadapter.output.fastq_trimmed
	output:
		sam = 'processed/{sample}_bwamem.sam'
	shell:
		'bwa mem {input.reference_file} {input.fastq} > {output.sam}'

rule generate_counts:
	input:
		sam_file = rules.bwa_mem.output.sam,
		reference_file = get_reference
	output:
		strict_count_file = 'output/{sample}_strict_match_counts.csv',
		sra_count_file = 'output/{sample}_sra_counts.csv',
		bar_plot_file_strict = 'final_output/{sample}_library_representation_strict.png',
		bar_plot_file_sra = 'final_output/{sample}_library_representation_sra.png',
		mapping_summary_file = 'output/{sample}_mapping_summary.csv'
	shell:
		'python resources/generate_counts.py -a {input.sam_file} -r {input.reference_file} -s {output.strict_count_file} -b {output.sra_count_file} -g {output.bar_plot_file_strict} -f {output.bar_plot_file_sra} -m {output.mapping_summary_file}'

rule aggregate_mapping_stats:
	input:
		summary_list = sorted(expand(rules.generate_counts.output.mapping_summary_file,sample=sample_list)),
		strict_files = sorted(expand(rules.generate_counts.output.strict_count_file,sample=sample_list)),
		sra_files = sorted(expand(rules.generate_counts.output.sra_count_file,sample=sample_list))
	output:
		aggregate_file = 'final_output/' + ngs_run + '_library_mapping_stats.csv',
		strict_matrix = 'final_output/' + ngs_run + '_count_matrix_strict.csv',
		sra_matrix = 'final_output/' + ngs_run + '_count_matrix_sra.csv'
	shell:
		'python resources/aggregate_stats.py -i {input.summary_list} -a {input.strict_files} -b {input.sra_files}  -m {output.aggregate_file} -s {output.strict_matrix} -c {output.sra_matrix}'

rule process_libraries:
	input:
		report_file = rules.aggregate_mapping_stats.output.aggregate_file,
		strict_matrix_file = rules.aggregate_mapping_stats.output.strict_matrix,
		sra_matrix_file = rules.aggregate_mapping_stats.output.sra_matrix
