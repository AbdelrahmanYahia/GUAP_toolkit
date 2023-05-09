
rule trimmomatic:
    input:
        unpack(get_raw_fasta)

    output:
        log="logs/trimmomatic/{sample}.log",
        summary="logs/trimmomatic/{sample}.summary",
        nf1 = f"trimmomatic/{{sample}}_{RS}1.{EXT}",
        nf2 = f"trimmomatic/{{sample}}_{RS}2.{EXT}",
        nfu1=temp(f"trimmomatic/U/{{sample}}_R1_U.{EXT}"),
        nfu2=temp(f"trimmomatic/U/{{sample}}_R2_U.{EXT}")

    benchmark: "benchamrks/QC/{sample}_trim.txt"
    threads: config["trim_t"]
    params:
        size = config["slidingwindow_size"],
        quality = config["slidingwindow_quality"],
        extra = config['trimmomatic_extra_args'],
        minlen = config['trim_min_length'],

    log: 
        "logs/trimmomatic/{sample}.txt"
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 -trimlog {output.log} \
                -summary {output.summary} {input.R1} {input.R2} {output.nf1} {output.nfu1} {output.nf2} {output.nfu2} \
                SLIDINGWINDOW:{params.size}:{params.quality} MINLEN:{params.minlen} > {log} 2>&1
        """