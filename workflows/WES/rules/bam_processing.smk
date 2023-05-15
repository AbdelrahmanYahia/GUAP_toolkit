rule sort_bam:
    input:
        "{aligner}/{sample}.bam"

    output:
        out = "{aligner}/{sample}_sorted.bam"

    shell:
        "samtools sort {input} -o {output.out}"


rule mrk_duplicates_picard:
    input:
        "{aligner}/{sample}_sorted.bam"

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
        gatk --java-options "-Xmx8G" BaseRecalibrator -R {params.ref} \
            -I {input} --known-sites {params.known_sites} \
            -O {output} 
        """

rule BaseRecalibrator:
    input: 
        bam = "{aligner}/{sample}_picard.dedub.bam",
        report = "GATK/{sample}_{aligner}_picard.report"
    benchmark: "benchamrks/{sample}_{aligner}_GATK_apply_BQSR.txt"
    output: "GATK/{sample}_{aligner}_picard.pqsr.bam"
    threads: 1
    params: 
        known_sites = known_variants,
        ref = ref_fasta
    shell:
        """
        gatk --java-options "-Xmx8G" ApplyBQSR  -R {params.ref} \
            -I {input.bam} --emit-original-quals \
            -bqsr {input.report} -O {output} \
            --add-output-sam-program-record
        """


rule bqsr_calibrated_report:
    input: "GATK/{sample}_{aligner}_picard.pqsr.bam"
    output: "GATK/{sample}_{aligner}_picard_pqsr.report"
    params: 
        known_sites = known_variants,
        ref = ref_fasta
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx8G" BaseRecalibrator -R {params.ref} \
            -I {input} --known-sites {params.known_sites} \
            -O {output} 
        """
