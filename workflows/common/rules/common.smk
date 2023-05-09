PATH = config["path"]
EXTT = config["ext"]
EXT = EXTT.replace(".gz","")
RS = config["R"]
TAIL = config["tail"]
sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file,index_col=0)
SAMPLES = list(SampleTable.index)
SAMPLES_IDs = list(SampleTable.iloc[:, 0])
ALL_THREADS = config["threads"]
MEM = config["total_mem"]
GUAP_FOLDER = config["GUAP_DIR"]

def get_cutadapt_input():
    if config["trimmomatic"]:
        return("trimmomatic")
    else:
        return("reads")

def get_multiqc_input(wildcards):
    final_input = []
    if config["fastqc"]:
        final_input = expand(f"QC/{{sample}}_R{{R}}_fastqc.html",sample=SAMPLES_IDs, R=[1,2])
    if config["trimmomatic"]:
        final_input.extend(expand(f"trimmomatic/{{sample}}_R1{EXT}", sample=SAMPLES_IDs))
    if config["cutadapt"]:
        final_input.extend(expand(f"cutadapt/{{sample}}_R1{EXT}", sample=SAMPLES_IDs))   
    return final_input


rule getreads:
    output:
        R1 = (f"reads/{{sample}}_R1{EXTT}"),
        R2 = (f"reads/{{sample}}_R2{EXTT}")
    params: 
        input_path=PATH
    run:
        shell("cp {params.input_path}/{wildcards.sample}*_{RS}1{TAIL}{EXTT} {output.R1}"),
        shell("cp {params.input_path}/{wildcards.sample}*_{RS}2{TAIL}{EXTT} {output.R2}")

rule decompress:
    input: 
        f"reads/{{sample}}_R{{R}}{EXTT}"
    output:
        temp(f"reads/{{sample}}_R{{R}}{EXT}")
    shell:
        "gunzip {input}"

rule cleaning:
    input:
        "multiqc/multiqc_report.html"
    output:
        "log.txt"
    shell:
        """
        rm -r reads >> log.txt
        echo -e "\033[;34;1mALL DONE\033[;39;m"
        """

include: "QCandFilter.smk"
include: "aligners.smk"
