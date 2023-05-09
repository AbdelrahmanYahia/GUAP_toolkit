
rule kallisto:
    input: 
        unpack(get_align_input)
    threads: config["align_t"]
    params:
        index = ref_index,
        extra_args = config['aligner_extra_args'],
    log:
        "logs/kallisto/{sample}.txt"
    output:
        outdir = directory("kallisto/{sample}"),
        bam="kallisto/{sample}/pseudoalignments.bam",
        counts="kallisto/{sample}/abundance.tsv"

    benchmark: "benchamrks/kallisto/{sample}.txt"
    shell:
        """
        kallisto quant -i {params.index} \
            -o {output.outdir} {input.R1} {input.R2} \
            -t {threads} \
            --pseudobam {params.extra_args} 2> {log}
        """

rule merge_kallisto_counts:
    input: 
        expand("kallisto/{sample}/abundance.tsv", sample=samples)
    threads: 1

    output:
        "Counts/kallisto_quant_ENST.counts"

    shell:
        """
        paste {input} | cut -f 1,4,9,14,19,24,29 > Counts/kallisto_counts.temp
        ls -1 {input} | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){{my $sample_name = $1; $sample_name =~ s/kallisto\///g; print "\t$sample_name"}}' | perl -ne 'print "Geneid$_\n"' > Counts/header.tsv
        cat Counts/header.tsv Counts/kallisto_counts.temp | grep -v "est_counts" > Counts/kallisto_quant.temp
        mv Counts/kallisto_quant.temp {output}
        rm -f Counts/header.tsv
        """

rule get_gene_symbol_kallisto:
    input: 
        "Counts/kallisto_quant_ENST.counts"    
    threads: 1
    params:
        db = "EnsDb",
        org = "Homo sapeins"
    output:
        "Counts/kallisto_quant.counts"
        
    shell:
        f"""
        Rscript {source_dir}/workflows/RNAseq/scripts/geneIDtoSymbol.R -i {{input}} -o {{output}} --transcript
        """