rule sort_and_convert_sam:
    input:
        "{aligner}/{sample}.sam"

    output:
        "{aligner}/{sample}.sorted.bam"
    
    shell:
        """
        samtools sort {input} -o {output}
        samtools index {output}
        """


rule mrk_duplicates_picard:
    input:
        "{aligner}/{sample}.sorted.bam"

    output:
        bam = "{aligner}/{sample}_picard.dedub.bam",
        matrix = "{aligner}/{sample}_picard.dedub.matrix"

    log: 
        "logs/{sample}_{aligner}_picard.dedub.log"

    shell:
        """
        picard MarkDuplicates -I {input} \
            -O {output.bam} \
            -M {output.matrix} > {log}
        """


rule bqsr_dedub_report:
    input: "{aligner}/{sample}_picard.dedub.bam"
    output: "GATK/{sample}_{aligner}_picard.report"
    params: 
        known_sites = known_variants,
        ref = ref_fasta
    threads: 1
    benchmark: "benchamrks/{sample}_{aligner}_GATK_pqsr.txt"
    shell:
        """
        gatk --java-options "-Xmx16G" BaseRecalibrator -R {params.ref} \
            -I {input} --known-sites {params.known_sites} \
            -O {output} 
        """
