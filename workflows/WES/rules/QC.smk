
include: f'{common_rules}/fastqc.smk'
include: f'{common_rules}/trimmomatic.smk'


rule QC_alignment:
    input:
        "{aligner}/{sample}_sorted.bam"

    output:
        cov = "{sample}_{aligner}.cov",
        stats = "{sample}_{aligner}.stats"

    shell:
        """
        samtools depth {input} | awk '{{sum+=$3}} END {{print "Average = ",sum/NR, "No of covered Nuc = ", NR}}' > {output.cov}
        samtools flagstat {input} > {output.stats}
        """

rule multiqc:
    input:
        get_final_output
    output:
        "multiqc/multiqc_report.html"
    shell:
        "multiqc . -o multiqc/"