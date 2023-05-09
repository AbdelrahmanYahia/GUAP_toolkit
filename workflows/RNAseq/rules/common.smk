import os
import pandas as pd
# get options
PATH = config["path"]
EXT = config["ext"]
RS = config["R"]
TAIL = f"_{config['tail']}"
sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file,index_col=0)
samples = list(SampleTable.index) # same as IDs
files_R1s = list(SampleTable.iloc[:, 0])
files_R2s = list(SampleTable.iloc[:, 8])
samples_IDs = list(SampleTable.iloc[:, 2])
samples_names = list(SampleTable.iloc[:, 1])
ALL_THREADS = config["threads"]
MEM = config["total_mem"]
GUAP_FOLDER = config["GUAP_DIR"]
samples_dir = PATH
ref_fasta = config["reference_fasta"]
ref_gtf = config["gtf_file"]
ref_index = config["reference_index"]
aligner = config["aligner"]
quantifier = config["quantifier"]
R = [1, 2]
R1_pattern = config["R1_pattern"]
R2_pattern = config["R2_pattern"]
bamfilename = "Aligned.sortedByCoord.out"
source = PATH
lane="_L001"
working_dir = config["working_dir"]
source_dir = config["GUAP_DIR"]
dirstr = "star/"
repattern = "/Aligned.sortedByCoord.out.bam"



if aligner == 'hisat2':
    bamfilename = 'aligned_sorted'
    dirstr = "hisat2/"
    repattern = "/aligned_sorted_bam"
elif aligner == 'kallisto':
    bamfilename = 'pseudoalignments'
    quantifier = "quant"


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

# rule decompress:
#     input: 
#         f"reads/{{sample}}_{RS}{{R}}{EXTT}"
#     output:
#         temp(f"reads/{{sample}}_{RS}{{R}}{EXT}")
#     shell:
#         "gunzip {input}"
        







