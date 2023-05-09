
rule Fastqc:
    input:
        f"reads/{{sample}}_R{{R}}{EXT}"
    log:
        "logs/QC/QC_{sample}_R{R}.log"
    output:
        zip=f"QC/{{sample}}_R{{R}}_fastqc.zip",
        html=f"QC/{{sample}}_R{{R}}_fastqc.html"

    threads: USE_THREADS
    params:
        path="QC"
    
    shell:
        "fastqc {input} --threads {threads} -o {params.path} >> {params.path}/log.txt 2> {log}"


rule trimmomatic:
    input:
        R1 = (f"reads/{{sample}}_R1{EXT}"),
        R2 = (f"reads/{{sample}}_R2{EXT}")

    output:
        log="logs/trimmomatic/{sample}.log",
        summary="logs/trimmomatic/{sample}.summary",
        nf1=f"trimmomatic/{{sample}}_R1{EXT}",
        nf2=f"trimmomatic/{{sample}}_R2{EXT}",
        nfu1=temp(f"trimmomatic/U/{{sample}}_R1_U{EXT}"),
        nfu2=temp(f"trimmomatic/U/{{sample}}_R2_U{EXT}")
    
    threads: USE_THREADS
    params:
        minlen=config["trim_min_length"]
    log: 
        "logs/trimmomatic/{sample}.txt"
    shell:
        """
            trimmomatic PE -threads {threads} -phred33 -trimlog {output.log} \
                    -summary {output.summary} {input.R1} {input.R2} {output.nf1} {output.nfu1} {output.nf2} {output.nfu2} \
                    SLIDINGWINDOW:4:10 MINLEN:{params.minlen} >> {log} 2>&1
        """

rule multiqc:
    input:
        get_multiqc_input
    output:
        "multiqc/multiqc_report.html"
    shell:
        "multiqc . -o multiqc/"

