include: "common.smk"

rule bowtie_aling:
    input:
       "OUT/cutadapt/{sample}.fastq"
    output:
        "OUT/miRNA_direct/mapped_reads/{sample}.sam"
    params:
        genome=config["indexes"]["bowtie_index_hum"]
    
    log:
        "OUT/miRNA_direct/logs/bowtie_{sample}.log"
    
    threads: 10
    shell:
        "bowtie --threads {threads} -q -v 0 -k 10 -S -t {params.genome} {input}  -S {output} 2> {log}"

rule Counts:
    input:
        insure_output = expand(rules.bowtie_aling.output,sample=SAMPLES)
        # allsams = glob_wildcards("mapped_reads/{sample}.sam")
    output:
        temp("OUT/miRNA_direct/counts/All.txt")
    params:
        anno=config["indexes"]["mirbase_annotation"]
    threads:8
    log:"OUT/miRNA_direct/logs/featureCounts.log"
    shell:
        "featureCounts -t miRNA -g Name -O -T {threads} -s 1 -M -a {params.anno} -o {output} {input} 2> {log}"
        "\n"
        

rule edit_counts:
    input:"OUT/miRNA_direct/counts/All.txt"
    output:"OUT/miRNA_direct/counts/final_counts.txt"
    shell:
        "tail -n +2 {input} | cut -f 1,7- > {output}"