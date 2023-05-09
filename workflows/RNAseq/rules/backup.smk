import os
import pandas as pd


PATH = config["path"]
EXTT = config["ext"]
EXT = EXTT.replace(".gz","")
RS = config["R"]
TAIL = config["tail"]

if TAIL is None:
    TAIL = ""

sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file,index_col=0)
samples = list(SampleTable.index)
samples_IDs = list(SampleTable.iloc[:, 0])
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

if config["skip_QC"]:
    source = "reads"
else:
    source = "trimmomatic"

bam_ext = "Aligned.sortedByCoord.out.bam"


workdir: config["working_dir"]

def get_final_output(wildcards):
    final_input = []
    
    if aligner in ["hisat2", "star"]:
        final_input.extend(expand("Counts/{aligner}_{quantifier}.counts", 
            aligner = aligner,
            quantifier = quantifier)
        )
    else:
        final_input.extend(expand("Counts/{aligner}.counts", 
            aligner = aligner,
            quantifier = quantifier)
        )
    final_input.extend(expand("{aligner}/{sample}_{bam_ext}.stats", sample=samples_IDs,aligner=aligner, bam_ext=bam_ext))
    print(final_input)
    return final_input


workdir: config["working_dir"]


rule all:
    input: "multiqc/multiqc_report.html"


rule getreads:
    output:
        R1 = (f"reads/{{sample}}_R1{EXTT}"),
        R2 = (f"reads/{{sample}}_R2{EXTT}")
    params: 
        input_path=PATH
    run:
        shell("cp {params.input_path}/{wildcards.sample}*{RS}1{TAIL}{EXTT} {output.R1}"),
        shell("cp {params.input_path}/{wildcards.sample}*{RS}2{TAIL}{EXTT} {output.R2}")

rule decompress:
    input: 
        f"reads/{{sample}}_{RS}{{R}}{EXTT}"
    output:
        temp(f"reads/{{sample}}_{RS}{{R}}{EXT}")
    shell:
        "gunzip {input}"

rule Fastqc:
    input:
        f"{source}/{{sample}}_{RS}{{R}}{EXT}"
    log:
        f"logs/QC/QC_{{sample}}_{RS}{{R}}.log"
    output:
        zip=f"QC/{{sample}}_R{{R}}_fastqc.zip",
        html=f"QC/{{sample}}_R{{R}}_fastqc.html"
    benchmark: f"benchamrks/QC/{{sample}}_{RS}{{R}}fastqc.txt"
    threads: 2
    params:
        path="QC"
    
    shell:
        "fastqc {input} --threads {threads} -o {params.path} > {log} 2>&1"

rule trimmomatic:
    input:
        R1 = (f"reads/{{sample}}_{RS}1{EXT}"),
        R2 = (f"reads/{{sample}}_{RS}2{EXT}")

    output:
        log="logs/trimmomatic/{sample}.log",
        summary="logs/trimmomatic/{sample}.summary",
        nf1=f"trimmomatic/{{sample}}_R1{EXT}",
        nf2=f"trimmomatic/{{sample}}_R2{EXT}",
        nfu1=temp(f"trimmomatic/U/{{sample}}_R1_U{EXT}"),
        nfu2=temp(f"trimmomatic/U/{{sample}}_R2_U{EXT}")

    benchmark: "benchamrks/QC/{sample}_trim.txt"
    threads: config["trim_t"]
    params:
        size = config["slidingwindow_size"],
        quality = config["slidingwindow_quality"],
        extra = config['trimmomatic_extra_args'],
        minlen = config['trim_min_length']
    log: 
        "logs/trimmomatic/{sample}.txt"
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 -trimlog {output.log} \
                -summary {output.summary} {input.R1} {input.R2} {output.nf1} {output.nfu1} {output.nf2} {output.nfu2} \
                SLIDINGWINDOW:{params.size}:{params.quality} MINLEN:{params.minlen} > {log} 2>&1
        """

rule multiqc:
    input:
        get_final_output
    output:
        "multiqc/multiqc_report.html"
    shell:
        "multiqc . -o multiqc/"


rule mapping_qc:
    input:
       "{aligner}/{sample}_{bam_ext}"
    output:
        "{aligner}/{sample}_{bam_ext}.stats"
    shell:
        """
        samtools stats {input} > {output}
        """

rule hisat2:
    input: 
        R1 = f"{source}/{{sample}}_{RS}1{EXT}",
        R2 = f"{source}/{{sample}}_{RS}2{EXT}"
    threads: config["align_t"]
    params:
        index = ref_index,
        fa = ref_fasta,
        extra_args = config['aligner_extra_args']
    log:
        "logs/hisat2/{sample}.txt"
    output:
        "hisat2/{sample}.sam"
    benchmark: "benchamrks/hisat2/{sample}.txt"
    shell:
        """
        hisat2 -x {params.index} --threads {threads} \
            -1 {input.R1} -2 {input.R2} {params.extra_args} > {output} 2> {log}
        """

rule convert_sam:
    input: 
    	"hisat2/{sample}.sam"
    
    output: 
    	"hisat2/{sample}.bam"
    
    shell: 
    	"samtools view -hbo {output} {input}"

rule sort_bam:
    input:
        "hisat2/{sample}.bam"

    output:
        "hisat2/{sample}_sorted.bam"

    shell:
        "samtools sort {input} -o {output}"


rule STAR_align:
    input:
        R1 = f"{source}/{{sample}}_{RS}1{EXT}",
        R2 = f"{source}/{{sample}}_{RS}2{EXT}"

    output:
        "star/{sample}_Aligned.sortedByCoord.out.bam"

    threads: config["align_t"]
    params:
        index = ref_index,
        extra_args = config['aligner_extra_args']

    log: "logs/star/{sample}_STAR.log"

    benchmark: "benchamrks/star/{sample}.txt"

    shell:
        """
    STAR --genomeDir {params.index} \
        --outSAMunmapped Within \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD MC \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --runThreadN {threads} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFileNamePrefix "star/{wildcards.sample}_" \
        --readFilesIn {input.R1} {input.R2} {params.extra_args} > {log} 2>&1
        """


rule featureCounts:
    input:
        expand("{aligner}/{sample}_{bam_ext}", 
                aligner = aligner,
                sample = samples_IDs, 
                bam_ext = bam_ext)

    output: 
        "Counts/{aligner}_featurecounts.precounts"

    log: 
        "logs/counts/{aligner}_featurecounts.log"

    threads: config["quant_t"]
    benchmark: "benchamrks/counts/{aligner}_fc.txt"
    params:
        gtf = ref_gtf,
        feature = config["feature"],
        extra = config['quantifier_extra_args']
        
    shell:
        """
        featureCounts -T {threads} -g {params.feature} \
            -o {output} \
            -a {params.gtf} \
            {input} \
            > {log} 2>&1
        """

rule modify_Fcounts:
    input:
        "Counts/{aligner}_featurecounts.precounts"
    output:
        "Counts/{aligner}_featurecounts.counts"
    shell:
        "tail -n +2 {input} | cut -f 1,7- > {output}"


rule htseq:
    input:
        expand("{aligner}/{sample}_{bam_ext}", 
                sample = samples_IDs,
                bam_ext = bam_ext,
                aligner = aligner)

    output: 
        "Counts/{aligner}_htseq.counts"
    benchmark: "benchamrks/counts/{aligner}_htseq.txt"
    threads: config["quant_t"]

    params:
        gtf = ref_gtf,
        feature = config["feature"],
        extra = config['quantifier_extra_args']

    log: 
        "logs/counts/{aligner}_htseq.log"

    shell:
        """
        echo -e "geneid,{input}\n" > {input}
        htseq-count -t transcript {params.extra}\
                    -f bam -i {params.feature} \
                    -r pos {input} \
                    {params.gtf} >> {input} 2> {log}
        """
