include: "common.smk"

rule bowtie_aling:
    input:
       "OUT/cutadapt/{sample}.fastq"
    output:
        "OUT/miRNA_human/mapped_reads/{sample}.sam"
    params:
        genome=config["indexes"]["bowtie_index_hum"]
    
    threads: 20
    shell:
        """bowtie -p {threads} -v 0 -t -m 2 \
                    --best --strata -e 99999 \
                    --chunkmbs 2048 \
                    {params.genome} \
                    -q {input} \
                    -S {output}
        """

rule calcstats:
    input:
       rules.bowtie_aling.output
    output:
        "OUT/miRNA_human/mapped_reads/{sample}_human.stats"
    threads: 1
    shell:
        """
        samtools stats {input} > {output}
        """

rule Counts:
    input:
        insure_output = expand(rules.bowtie_aling.output,sample=SAMPLES)
        # allsams = glob_wildcards("mapped_reads/{sample}.sam")
    output:
        temp("OUT/miRNA_human/counts/All.txt")
    params:
        anno=config["indexes"]["mirbase_annotation"]
    threads:8
    log:"OUT/miRNA_direct/logs/featureCounts.log"
    shell:
        "featureCounts -t miRNA -g Name -O -T {threads} -s 1 -M -a {params.anno} -o {output} {input} 2> {log}"
        "\n"
        

rule edit_counts:
    input:"OUT/miRNA_human/counts/All.txt"
    output:"OUT/miRNA_human/counts/final_human_counts.txt"
    shell:
        "tail -n +2 {input} | cut -f 1,7- > {output}"