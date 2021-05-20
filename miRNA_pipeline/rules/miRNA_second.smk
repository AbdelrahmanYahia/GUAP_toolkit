include: "common.smk"

rule bowtie_mirBase:
    input:
       "OUT/cutadapt/{sample}.fastq"
    output:
        BAM="OUT/miRNA_second/miRBase/bam/{sample}.bam",
        UNM="OUT/miRNA_second/miRBase/unmapped/{sample}_unmapped.fastq"
    params:
        mirbase_ref=config["indexes"]["bowtie_index_mirbase"]
    
    log:
        "OUT/miRNA_second/logs/miRBase/bowtie_{sample}.log"
    
    threads: 5
    shell:
        "bowtie -n 0 -l 32 --norc --best --strata -m 1 --threads {threads} {params.mirbase_ref} {input} --un {output.UNM} -S | samtools sort --output-fmt BAM - > {output.BAM} 2> {log}"

rule miRBase_count:
    input:
       rules.bowtie_mirBase.output.BAM
    output:
        "OUT/miRNA_second/counts/{sample}.mirna.counts"
    log:
        "OUT/miRNA_second/logs/miRBase/samtools_{sample}.log"
    threads: 1
    shell:
        """
        samtools index {input} 2> {log}
        samtools idxstats {input} | cut -f1,3 - > {output} 2>> {log}
        """

rule bowtie_human:
    input:
       rules.bowtie_mirBase.output.UNM
    output:
        "OUT/miRNA_second/human/bam/{sample}.bam"
    params:
        hsa_ref=config["indexes"]["bowtie_index_hum"]    
    log:
        "OUT/miRNA_second/logs/human/bowtie_{sample}.log"
    
    threads: 5
    shell:
        "bowtie -n 1 -l 32 --norc --best --strata -m 1 --threads {threads} {params.hsa_ref} {input}  -S | samtools sort --output-fmt BAM - > {output} 2> {log}"

rule human_tag_bam:
    input:
       rules.bowtie_human.output
    output:
        "OUT/miRNA_second/human/tbam/{sample}_tagged.bam"
    params:
        bed_file=config["indexes"]["miRBED"]
    log:
        "OUT/miRNA_second/logs/human/samtools_{sample}.log"
    threads: 1
    shell:
        """
        samtools index {input} 2> {log}
        bedtools tag -i {input} -files {params.bed_file} -names -tag XQ > {output} 2>> {log}
        """

rule human_counting:
    input:
       rules.human_tag_bam.output
    output:
        "OUT/miRNA_second/counts/{sample}.human.counts"
    shell:
        "scripts/count_human.sh -i {input} -o {output}"

rule merge_counts:
    input:
       insure_output = expand(rules.human_counting.output,sample=SAMPLES),
       insure_output2 = expand(rules.miRBase_count.output,sample=SAMPLES)

    output:
        "all.csv"
    params:
        path="OUT/miRNA_second/counts/"
    shell:
        """
        python3 scripts/merge.py -i {params.path} -o {output}
        mv OUT/miRNA_second/counts/all.csv .
        """