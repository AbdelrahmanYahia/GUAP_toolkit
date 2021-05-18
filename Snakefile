import os
import pandas as pd

configfile: "config.yaml"

include: "rules/common.smk"

sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file,index_col=0)
SAMPLES = list(SampleTable.index)

rule all:
    input:
        expand("OUT/rawQC/{sample}_fastqc.{extension}", sample=SAMPLES, extension=["zip","html"]),
        # expand("OUT/mapped_reads/{sample}.sam", sample=SAMPLES),
        "OUT/counts/All.txt"