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
sample_to_five_prime_adapter_R1 = defaultdict(str)
sample_to_five_prime_adapter_R2 = defaultdict(str)
sample_to_library_type = defaultdict(str)

# go through samples
for sample in samples:
	for sample_name in sample:
		sample_list.append(sample_name)
		ref_file = sample[sample_name]['reference']
		five_prime_adapter_R1 = sample[sample_name]['five_prime_adapter_R1']
		five_prime_adapter_R2 = sample[sample_name]['five_prime_adapter_R2']
		if ref_file not in reference_list:
			reference_list.append(ref_file)
		sample_to_reference[sample_name] = ref_file
		sample_to_five_prime_adapter_R1[sample_name] = five_prime_adapter_R1
		sample_to_five_prime_adapter_R2[sample_name] = five_prime_adapter_R2
		
# functions to get variables
def get_reference(wildcards):
	return sample_to_reference[wildcards.sample]

def get_five_prime_adapter_R1(wildcards):
	return sample_to_five_prime_adapter_R1[wildcards.sample]

def get_five_prime_adapter_R2(wildcards):
	return sample_to_five_prime_adapter_R2[wildcards.sample]

rule cut_five_adapter_R1:
	input:
		fastq = 'data_files/{sample}_R1.fastq.gz'
	params:
		adapter = get_five_prime_adapter_R1
	output:
		fastq_trimmed_R1 = 'data_files/{sample}_R1_trimmed.fastq'
	shell:
		'cutadapt -g {params.adapter} --discard-untrimmed --quiet -o {output.fastq_trimmed_R1} {input.fastq}'

rule cut_five_adapter_R2:
	input:
		fastq = 'data_files/{sample}_R2.fastq.gz'
	params:
		adapter = get_five_prime_adapter_R2
	output:
		fastq_trimmed_R2 = 'data_files/{sample}_R2_trimmed.fastq'
	shell:
		'cutadapt -g {params.adapter} --discard-untrimmed --quiet -o {output.fastq_trimmed_R2} {input.fastq}'

rule bwa_mem:
	input:
		reference_file = get_reference,
		fastq = rules.cut_five_adapter_R1.output.fastq_trimmed_R1
	output:
		sam = 'processed/{sample}_bwamem.sam'
	shell:
		'bwa mem {input.reference_file} {input.fastq} > {output.sam}'

rule generate_counts_single_guide:
	input:
		sam_file = rules.bwa_mem.output.sam,
		reference_file = get_reference,
	output:
		strict_count_file = 'final_output/{sample}_strict_match_counts.csv',
		bar_plot_file_strict = 'final_output/{sample}_library_representation_strict.png',
		mapping_summary_file = 'final_output/{sample}_mapping_summary_sg.csv'
	shell:
		'python resources/generate_counts.py -i {input.sam_file} -r {input.reference_file} -s {output.strict_count_file} -g {output.bar_plot_file_strict} -m {output.mapping_summary_file} -l single_guide'

rule generate_counts_orf:
	input:
		fastq_r1 = rules.cut_five_adapter_R1.output.fastq_trimmed_R1,
		reference_file = get_reference
	output:
		orf_count_file = 'final_output/{sample}_counts_orf.csv',
		bar_plot_file_strict = 'final_output/{sample}_representation_orf_strict.png',
		mapping_summary_file = 'final_output/{sample}_mapping_summary_orf.csv'
	shell:
		'python resources/generate_counts.py -i {input.fastq_r1} -r {input.reference_file} -s {output.orf_count_file} -g {output.bar_plot_file_strict} -m {output.mapping_summary_file} -l orf'


rule generate_counts_dual_guide:
	input:
		fastq_r1 = rules.cut_five_adapter_R1.output.fastq_trimmed_R1,
		fastq_r2 = rules.cut_five_adapter_R2.output.fastq_trimmed_R2,
		reference_file = get_reference
	output:
		dg_count_file = 'final_output/{sample}_counts_dg.csv',
		bar_plot_file_strict = 'final_output/{sample}_representation_dg_strict.png',
		mapping_summary_file = 'final_output/{sample}_mapping_summary_dg.csv'
	shell:
		'python resources/generate_counts.py -i {input.fastq_r1} -r {input.reference_file} -s {output.dg_count_file} -g {output.bar_plot_file_strict} -m {output.mapping_summary_file} -l dual_guide'

rule process_single_guide_libraries:
	input:
		count_list = sorted(expand(rules.generate_counts_single_guide.output.strict_count_file,sample=sample_list)),
		png_list = sorted(expand(rules.generate_counts_single_guide.output.bar_plot_file_strict,sample=sample_list))

rule process_orf_libraries:
	input:
		count_list = sorted(expand(rules.generate_counts_orf.output.orf_count_file,sample=sample_list))

rule process_dual_guide_libraries:
	input:
		count_list = sorted(expand(rules.generate_counts_dual_guide.output.dg_count_file,sample=sample_list))
