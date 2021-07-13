include: "common.smk"

rule Fastqc:
    input:
        rawread="Data/samples/{sample}.fastq"
    log:
        "OUT/logs/QC/QC_{sample}.log"
    output:
        zip="OUT/rawQC/{sample}_fastqc.zip",
        html="OUT/rawQC/{sample}_fastqc.html"
    threads:5
    params:
        path="OUT/rawQC/"
    
    shell:
        "fastqc {input.rawread} --threads {threads} -o {params.path} 2> {log}"

rule cutadapt1:
    input:
        "Data/samples/{sample}.fastq"
    output:
        temp("OUT/cutadapt_temp/{sample}.fastq")
    log:
        "OUT/logs/cutAdapt/cut1_{sample}.log"
    params:
        adaptor=config["trimming"]["adaptor"]
    
    threads:5
    shell:
        "cutadapt -j {threads} -a {params.adaptor} -o {output} --minimum-length 23 {input} 2> {log}"

rule cutadapt:
    input:
        "OUT/cutadapt_temp/{sample}.fastq"
    output:
        "OUT/cutadapt/{sample}.fastq"
    log:
        "OUT/logs/cutAdapt/cutadapt_{sample}.log"
    threads:5
    shell:
        "cutadapt -j {threads} --minimum-length=12 --maximum-length=35 -u 4 -u -4 -o {output} {input} 2> {log}"