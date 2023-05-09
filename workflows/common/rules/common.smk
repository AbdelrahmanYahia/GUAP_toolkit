import os
import pandas as pd

# get options
PATH = config["path"]
EXTT = config["ext"]
RS = config["R"]
if config["decompress"]:
    EXT = EXTT.replace(".gz","")
else:
    EXT = EXTT

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
R = [1, 2]
R1_pattern = config["R1_pattern"]
R2_pattern = config["R2_pattern"]
source = PATH
lane="_L001"
working_dir = config["working_dir"]
source_dir = config["GUAP_DIR"]
common_rules = config["common_rules"]

include: f'{common_rules}/utils.smk'

rule decompress:
    input: 
        get_decompress_input
    output:
        temp(f"{source}/{{sample}}_{RS}{{R}}{TAIL}.{EXT}")
    shell:
        "gunzip {input}"
