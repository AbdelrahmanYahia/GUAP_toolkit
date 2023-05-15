
rule AnalyzeCovariates:
    input: 
        raw = "GATK/{sample}_{aligner}_picard.report", 
        bqsr = "GATK/{sample}_{aligner}_picard_pqsr.report"

    output: "GATK/{sample}_{aligner}_picard.pdf"
    threads: 1
    shell:
        """
        gatk AnalyzeCovariates -before {input.raw} \
            -after {input.bqsr} -plots {output}
        """

rule HaplotypeCaller:
    input: "GATK/{sample}_{aligner}_picard.pqsr.bam"
    output: "GATK/{sample}_{aligner}_picard.vcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file

    threads: config["threads_calling"]

    benchmark: "benchamrks/{sample}_{aligner}_GATK_HpTC.txt"

    shell:
        """
        gatk --java-options "-Xmx14G -XX:ParallelGCThreads={threads}" HaplotypeCaller -R {params.ref} \
            -L {params.bed} \
            -I {input} --native-pair-hmm-threads {threads} -O {output}
        """


rule mpileup:
    input:
        "{aligner}/{sample}_picard.dedub.bam"

    output:
        "mpileup/{sample}_{aligner}_picard.vcf"
    
    threads: config["threads_calling"]

    benchmark: "benchamrks/{sample}_{aligner}_mpileup.txt"

    params:
        fa = ref_fasta,
        bed = bed_file

    shell:
        """
        bcftools mpileup --threads 12 --regions-file {params.bed} \
            -Ou -f {params.fa} {input} |\
        bcftools call --threads 12 \
            -Ov -mv > {output}
        """

rule lofreq:
    input:
        "{aligner}/{sample}_picard.dedub.bam"

    output:
        "lofreq/{sample}_{aligner}_picard.vcf"
    
    threads: config["threads_calling"]

    benchmark: "benchamrks/{sample}_{aligner}_lofreq.txt"

    params:
        fa = ref_fasta,
        bed = bed_file

    shell:
        """
        lofreq call -f {params.fa} \
                --out {output} {input} \
                --call-indels 
        """

