import os
import pandas as pd

GUAP_FOLDER = config["GUAP_DIR"]
aligners = config["aligner"]
variant_callers = config["variant_caller"]
out_dir = config["working_dir"]
ref_bwa = config["reference_index"]
ref_bowtie2 = config["reference_index"]
ref_fasta = config["reference_fasta"]
known_variants = config["known_variants"]
bed_file = config["bed_file"]


bamfilename = ""

common_rules = config["common_rules"]
include: f'{common_rules}/common.smk'

def get_final_output(wildcards):
    final_input = []
    if config['skip_QC'] is False:
        final_input.extend(expand(
                f"QC/{{sample}}_{RS}{{R}}{TAIL}_fastqc.{{ext}}",
                ext = ["zip", "html"],
                R = [1, 2],
                sample = samples_names
        ))
        if config["trimmomatic"] is True:
            final_input.extend(expand(
                f"QC/{{sample}}_{RS}{{R}}_fastqc.{{ext}}",
                ext = ["zip", "html"],
                R = [1, 2],
                sample = samples
            ))
    final_input.extend(expand(
        "{sample}_{aligner}.stats",
         sample = samples, 
         aligner = aligners))

    final_input.extend(expand(
        "{variant_caller}/{sample}_{aligner}_picard.vcf", 
        sample = samples, 
        variant_caller = variant_callers,
        aligner = aligners)) 

    if variant_callers == "GATK":
        final_input.extend(expand("GATK/{sample}_{aligner}_picard.pdf", sample = samples, aligner = aligners))
    return final_input

include: 'variant_caller.smk'
include: 'bam_processing.smk'
include: 'QC.smk'
include: 'genome_align.smk'
include: 'utils.smk'




