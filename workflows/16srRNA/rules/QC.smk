include: f'{common_rules}/fastqc.smk'
include: f'{common_rules}/trimmomatic.smk'


rule multiqc:
    input:
        get_multiqc_input
    output:
        "multiqc/multiqc_report.html"
    shell:
        "multiqc . -o multiqc/"

