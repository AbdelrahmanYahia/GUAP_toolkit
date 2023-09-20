import pandas as pd

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

working_dir = config["working_dir"]
sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file)

files_R1s = set(SampleTable.iloc[:, 1-1])
files_R2s = set(SampleTable.iloc[:, 8-1])
samples = set(SampleTable.iloc[:, 3-1]) # sample full name
samples_IDs = set(SampleTable.iloc[:, 4-1])


ALL_THREADS = config["threads"]
MEM = config["total_mem"]
GUAP_FOLDER = config["GUAP_DIR"]
R = [1, 2]
source = PATH

working_dir = config["working_dir"]
source_dir = config["GUAP_DIR"]
common_rules = config["common_rules"]

samples_dir = config["input"]

out_dir = config["output"]


samples = samples_IDs
ref = config["index"]
gff = config["gff"]
kraken2 = config["kraken2_index"]
indir  = PATH
logdir = "logs"


def get_final_output(wildcards):
    final_output = []
    
    final_output.extend(expand(
            "{dir}/{sample}_stats.txt",
            sample = samples, dir = ["01_map_to_ref", "05_assembly-mapping"]
    ))

    final_output.extend(expand(
            "04_contamination/{dir}/{sample}.report",
            sample = samples, dir = ["unmapped_reads", "filttered"]
    ))


    final_output.extend(expand(
            "06_call_variants/{sample}.filtterd.vcf.gz",
            sample = samples
    ))

    final_output.extend(expand(
            "03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly-QC",
            sample = samples
    ))

    final_output.extend(expand(
            "03_Assembly/RAW/{sample}_filterd-Assembly-QC",
            sample = samples
    ))

    return final_output

def get_log_file(sample, rule):
    return f"{logdir}/{sample}_{rule}.log"

rule QC:
    input:
        R1=f"{indir}/{{sample}}_R1.fastq.gz",
        R2=f"{indir}/{{sample}}_R2.fastq.gz"
    output:
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq",
        rep="00_filterred/{sample}_R2.html"
    threads: 12

    log: get_log_file("{sample}", "QC")

    shell:
        """
        fastp --in1 {input.R1} --in2 {input.R2} \
            --out1 {output.R1} --out2 {output.R2} \
            --thread {threads} -h {output.rep} > {log} 2>&1
        """

rule index_ref:
    input:
        ref
    output:
        directory('bowtie2_index/')
    threads: 12

    log: get_log_file("Bowtie_index", "index_ref")

    shell:
        """
        mkdir -p bowtie2_index
        bowtie2-build --threads {threads} {input} bowtie2_index/Org_ref > {log} 2>&1
        """

rule map_to_ref:
    input:
        rules.index_ref.output,
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq",
    output:
        "01_map_to_ref/{sample}.sam"
    threads: 12

    log: get_log_file("{sample}", "map_to_ref")

    shell:
        """
        bowtie2 --threads {threads} -x bowtie2_index/Org_ref \
                -1 {input.R1} \
                -2 {input.R2} \
                -S {output} > {log} 2>&1
        """ 

rule fixmate:
    input: "01_map_to_ref/{sample}.sam"
    output: "01_map_to_ref/{sample}_fixmate.sam"
    shell: 
        """
        samtools sort -n \
                -O sam {input} | samtools fixmate \
                -m -O bam - {output} 
        """
rule sort:
    input: "01_map_to_ref/{sample}_fixmate.sam"
    output: "01_map_to_ref/{sample}_fixmate.bam"
    shell: 
        """
        samtools sort -O bam {input} \
                -o {output} 
        """

rule remove_duplicate:
    input: "01_map_to_ref/{sample}_fixmate.bam"
    output: "01_map_to_ref/{sample}_dedup.bam"
    shell: "samtools markdup -r -S {input} {output}"


rule filter_varaints:
    input:
        "06_call_variants/{sample}.vcf.gz"
    output:
        "06_call_variants/{sample}.filtterd.vcf.gz"
    threads: 12
    params:
        ref=ref
    shell:
        """
        bcftools view {input} | vcfutils.pl varFilter - > {output}
        """

rule call_variants:
    input:
        "01_map_to_ref/{sample}_dedup.bam"
    output:
        "06_call_variants/{sample}.vcf.gz"
    threads: 12
    params:
        ref=ref
    shell:
        """
        bcftools mpileup --threads {threads}  -Ou \
            -f {params.ref} {input} | bcftools call --ploidy 1 \
            --threads {threads} -mv -o {output}
        """

rule mapping_stats:
    input: "01_map_to_ref/{sample}_dedup.bam"
    output: 
        sam="01_map_to_ref/{sample}_stats.txt",
        quali=directory("01_map_to_ref/{sample}_stats")
    threads: 12

    log: get_log_file("{sample}", "mapping_stats")

    shell:
        """
        samtools flagstat {input} > {output.sam} &
        qualimap bamqc -bam {input} -outdir {output.quali} > {log} 2>&1
        wait
        """

rule extract_mapped:
    input: "01_map_to_ref/{sample}_dedup.bam"
    output: 
        bam="01_map_to_ref/{sample}_concordant.bam",
        R1="02_Extract_reads/mapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/mapped_reads/{sample}_R2.fastq",
        U="02_Extract_reads/mapped_reads/{sample}_U.fastq"
    shell: 
        """
        samtools view -b -f 3 {input} > 01_map_to_ref/{wildcards.sample}_concordant.bam
        samtools fastq -1 {output.R1} \
            -2 {output.R2} \
            -0 {output.U} 01_map_to_ref/{wildcards.sample}_concordant.bam 
        """

rule extract_unmapped:
    input: "01_map_to_ref/{sample}_dedup.bam"
    output: 
        bam="01_map_to_ref/{sample}_unmapped.bam",
        R1="02_Extract_reads/unmapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/unmapped_reads/{sample}_R2.fastq",
        U="02_Extract_reads/unmapped_reads/{sample}_U.fastq"
    shell: 
        """
        samtools view -b -f 4 {input} > 01_map_to_ref/{wildcards.sample}_unmapped.bam
        samtools fastq -1 {output.R1} \
            -2 {output.R2} \
            -0 {output.U} 01_map_to_ref/{wildcards.sample}_unmapped.bam 
        """

rule assemble_RAW_reads:
    input:
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq"
    output:
        eldir=directory("03_Assembly/RAW/{sample}_filterd-Assembly"),
        asmpl="03_Assembly/RAW/{sample}_filterd-Assembly/assembly.fasta"
    threads: 12

    log: get_log_file("{sample}", "assemble_RAW_reads")
    shell:
        """
        unicycler -1 {input.R1} -2 {input.R2} \
              -o {output.eldir} \
              --threads {threads} > {log} 2>&1
        """

rule assmble_qc_RAW_reads:
    input:
        "03_Assembly/RAW/{sample}_filterd-Assembly/assembly.fasta"
    output:
        directory("03_Assembly/RAW/{sample}_filterd-Assembly-QC")
    shell:
        "quast -o {output} {input}"

rule assemble_mapped_reads:
    input:
        R1="02_Extract_reads/mapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/mapped_reads/{sample}_R2.fastq"
    output:
        eldir=directory("03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly"),
        asmpl="03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly/assembly.fasta"
    threads: 12

    log: get_log_file("{sample}", "assemble_mapped_reads")

    shell:
        """
    unicycler -1 {input.R1} -2 {input.R2} \
              -o {output.eldir} \
              --threads {threads}  > {log} 2>&1
        """

rule assmble_qc_mapped_reads:
    input:
        "03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly/assembly.fasta"
    output:
        directory("03_Assembly/mapped_reads/{sample}_mapped_reads-Assembly-QC")
    log: get_log_file("{sample}", "assmble_qc_mapped_reads")

    shell:
        "quast -o {output} {input} "


rule kraken_RAW_reads:
    input:
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq"
    output:
        report="04_contamination/filttered/{sample}.report",
        log="04_contamination/filttered/{sample}.out"
    params:
        db=kraken2
    threads: 12

    shell:
        """
        kraken2 --db {params.db} --threads {threads} \
            --report {output.report} \
            --paired {input.R1} {input.R2} > {output.log}
        """

rule kraken_unmapped_reads:
    input:
        R1="02_Extract_reads/unmapped_reads/{sample}_R1.fastq",
        R2="02_Extract_reads/unmapped_reads/{sample}_R2.fastq"
    output:
        report="04_contamination/unmapped_reads/{sample}.report",
        log="04_contamination/unmapped_reads/{sample}.out"
    params:
        db=kraken2
    threads: 12

    shell:
        """
        kraken2 --db {params.db} --threads {threads} \
            --report {output.report} \
            --paired {input.R1} {input.R2} > {output.log}
        """

rule index_ref_assembly:
    input:
        "03_Assembly/RAW/{sample}_filterd-Assembly/assembly.fasta"
    output:
        directory("03_Assembly/RAW/{sample}_filterd-Assembly/bowtie2_index")
    threads: 12

    log: get_log_file("{sample}", "index_ref_assembly")

    shell:
        """
        mkdir -p 03_Assembly/RAW/{wildcards.sample}_filterd-Assembly/bowtie2_index/
        bowtie2-build --threads {threads} {input} 03_Assembly/RAW/{wildcards.sample}_filterd-Assembly/bowtie2_index/assembly > {log} 2>&1
        """

rule map_to_assembly:
    input:
        rules.index_ref_assembly.output,
        R1="00_filterred/{sample}_R1.fastq",
        R2="00_filterred/{sample}_R2.fastq",
    output:
        "05_assembly-mapping/{sample}.sam"
    threads: 12

    log: get_log_file("{sample}", "map_to_assembly")

    shell:
        """
        bowtie2 --threads {threads} -x 03_Assembly/RAW/{wildcards.sample}_filterd-Assembly/bowtie2_index/assembly \
                -1 {input.R1} \
                -2 {input.R2} \
                -S {output} > {log} 2>&1
        """ 

rule fixmate_assembly:
    input: "05_assembly-mapping/{sample}.sam"
    output: "05_assembly-mapping/{sample}_fixmate.sam"
    shell: 
        """
        samtools sort -n \
                -O sam {input} | samtools fixmate \
                -m -O bam - {output}
        """

rule sort_assembly:
    input: "05_assembly-mapping/{sample}_fixmate.sam"
    output: "05_assembly-mapping/{sample}_fixmate.bam"
    shell: 
        """
        samtools sort -O bam {input} \
                -o {output} 
        """

rule remove_duplicate_assembly:
    input: "05_assembly-mapping/{sample}_fixmate.bam"
    output: "05_assembly-mapping/{sample}_dedup.bam"
    shell: "samtools markdup -r -S {input} {output}"


rule mapping_stats_assembly:
    input: "05_assembly-mapping/{sample}_dedup.bam"
    output: 
        sam="05_assembly-mapping/{sample}_stats.txt",
        quali=directory("05_assembly-mapping/{sample}_stats")
    threads: 12

    log: get_log_file("{sample}", "mapping_stats_assembly")

    shell:
        """
        samtools flagstat {input} > {output.sam} &
        qualimap bamqc -bam {input} -outdir {output.quali} > {log} 2>&1
        wait
        """

        
rule multiqc:
    input:
        get_final_output
    
    conda: "../env/wes_gatk.yml"

    benchmark: "benchamrks/Multiqc/report.txt"

    output:
        "multiqc/multiqc_report.html"
    threads: 12

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        "multiqc . -o multiqc/"

