ALIGNER=config["aligner"]

rule aling:
    input:
       "OUT/cutadapt/{sample}.fastq"
    output:
        "OUT/miRNA_direct/mapped_reads/{sample}.sam"
    params:
        genome=config["indexes"]["bowtie_index_hum"],
        aligner=config["aligner"],

    
    log:
        "OUT/miRNA_direct/logs/bowtie_{sample}.log"
    
    threads: 10
    shell:
        "bowtie --threads {threads} -q -v 0 -k 10 -S -t {params.genome} {input}  -S {output} 2> {log}"
