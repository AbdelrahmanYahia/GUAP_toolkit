
rule mapping_qc:
    input:
       "{aligner}/{sample}/{bamfilename}.bam"
    output:
        "{aligner}/{sample}/{bamfilename}.stats"
    shell:
        """
        samtools stats {input} > {output}
        """
