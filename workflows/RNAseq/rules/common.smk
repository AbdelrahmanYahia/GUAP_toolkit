import os
import pandas as pd

common_rules = config["common_rules"]
include: f'{common_rules}/common.smk'

ref_fasta = config["reference_fasta"]
ref_gtf = config["gtf_file"]
ref_index = config["reference_index"]
aligner = config["aligner"]
quantifier = config["quantifier"]
bamfilename = "Aligned.sortedByCoord.out"
dirstr = "star/"
repattern = "/Aligned.sortedByCoord.out.bam"


if aligner == 'hisat2':
    bamfilename = 'aligned_sorted'
    dirstr = "hisat2/"
    repattern = "/aligned_sorted_bam"
elif aligner == 'kallisto':
    bamfilename = 'pseudoalignments'
    quantifier = "quant"

def get_align_input(wildcards):
    sample = wildcards.sample
    units = SampleTable.loc[sample]
    if config["trimmomatic"] is True:
            mydict = dict(
        zip(
            ["R1", "R2"],
                [
                    f"trimmomatic/{sample}_{RS}1.{EXT}",
                    f"trimmomatic/{sample}_{RS}2.{EXT}",
                ],
        )
    )
    else:
        source = PATH

        mydict = dict(
            zip(
                ["R1", "R2"],
                [
                    f"{config['input']}/{sample}_{units.sample_number}{lane}_{RS}1{TAIL}.{EXT}",
                    f"{config['input']}/{sample}_{units.sample_number}{lane}_{RS}2{TAIL}.{EXT}",
                ],
            )
        )
    return mydict

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


    if aligner in ["hisat2"] and quantifier in ["featurecounts", "htseq"]:
        final_input.extend(expand(
            "Counts/{aligner}_{quantifier}.counts", 
            aligner = aligner,
            quantifier = quantifier)
        )
    elif aligner in ["star"] and quantifier in ["featurecounts", "htseq"]:
        final_input.extend(expand(
            "Counts/{aligner}_{quantifier}.counts", 
            aligner = aligner,
            quantifier = ["quant", quantifier]))

    elif aligner in ["kallisto"]:
        final_input.extend(expand(
            "Counts/{aligner}_{quantifier}.counts", 
            aligner = aligner,
            quantifier = "quant")
        )

    if config["perform_DE"]:
        final_input.extend(expand(
            "Downstream/{aligner}_{quantifier}_tmp.txt", 
            aligner = aligner,
            quantifier = quantifier
        ))

    final_input.extend(expand("{aligner}/{sample}/{bamfilename}.stats", sample=samples,aligner=aligner, bamfilename=bamfilename))
    return final_input

include: 'utils.smk'
include: 'mapping_qc.smk'
include: 'pseudo_aligners.smk'
include: 'QC.smk'
include: 'quant.smk'
include: 'spliceaware_aligners.smk'
include: 'trim.smk'
include: 'DE.smk'
