include: "common.smk"

rule bowtie_aling_miR:
    input:
       "OUT/cutadapt/{sample}.fastq"
    output:
        "OUT/miRBase/mapped_reads/{sample}.bam"
    threads: 20
    shell:
        """bowtie -p {threads} -v 0 -t -m 2 \
                    --best --strata -e 99999 \
                    --chunkmbs 2048 \
                    /media/genomics/DB_Storage/db/miRNA_index/mature_idx \
                    -q {input} \
                    -S | samtools sort --output-fmt BAM - > {output}
        """

rule miRBase_count:
    input:
       rules.bowtie_aling_miR.output
    output:
        "OUT/miRBase/counts/{sample}.counts"
    threads: 1
    shell:
        """
        samtools index {input} 
        samtools idxstats {input} | cut -f1,3 - > {output}
        """

rule calcstats_miR:
    input:
       rules.bowtie_aling_miR.output
    output:
        "OUT/miRBase/mapped_reads/{sample}_miRBase.stats"
    threads: 1
    shell:
        """
        touch {output}
        samtools stats {input} > {output}
        """

