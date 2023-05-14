
rule hisat2:
    input: 
        unpack(get_align_input)
    threads: config["align_t"]
    params:
        index = ref_index,
        fa = ref_fasta,
        extra_args = config['aligner_extra_args']
    log:
        "logs/hisat2/{sample}.txt"
    output:
        "hisat2/{sample}/aligned.sam"
    benchmark: "benchamrks/hisat2/{sample}.txt"
    shell:
        """
        hisat2 -x {params.index} --threads {threads} \
            -1 {input.R1} -2 {input.R2} {params.extra_args} > {output} 2> {log}
        """

rule convert_sam:
    input: 
    	"hisat2/{sample}/aligned.sam"
    
    output: 
    	"hisat2/{sample}/aligned.bam"
    
    shell: 
    	"samtools view -hbo {output} {input}"

rule sort_bam:
    input:
        "hisat2/{sample}/aligned.bam"

    output:
        "hisat2/{sample}/aligned_sorted.bam"

    shell:
        "samtools sort {input} -o {output}"


rule STAR_align:
    input:
        unpack(get_align_input)

    output:
        "star/{sample}/Aligned.sortedByCoord.out.bam",
        "star/{sample}/ReadsPerGene.out.tab"

    threads: config["align_t"]
    params:
        index = ref_index,
        extra_args = config['aligner_extra_args']

    log: "logs/star/{sample}_STAR.log"

    benchmark: "benchamrks/star/{sample}.txt"

    shell:
        """
    STAR --genomeDir {params.index} \
        --outSAMunmapped Within \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD MC \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --runThreadN {threads} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFileNamePrefix "star/{wildcards.sample}/" \
        --readFilesIn {input.R1} {input.R2} {params.extra_args} > {log} 2>&1
        """


rule merge_star_counts:
    input: 
        ins_in = expand("star/{sample}/ReadsPerGene.out.tab", sample=samples)
        
    threads: 1
    params:
        strand = 1,
        in_dir = "star"
    output:
        temp("Counts/star_quant.temp")
    log:
        "logs/star/merging_counts.log"
    shell:
        f"""
        Rscript {source_dir}/workflows/RNAseq/scripts/merge_counts.R -i {{params.in_dir}} -o {{output}} -s {{params.strand}} > {{log}} 2>&1
        """

rule get_gene_symbol:
    input: 
        "Counts/star_quant.temp"    
    threads: 1
    params:
        db = "EnsDb",
        org = "Homo sapeins"
    output:
        "Counts/star_quant.counts"
    log:
        "logs/star/gene_to_symbl.log"
    shell:
        f"""
        Rscript {source_dir}/workflows/RNAseq/scripts/geneIDtoSymbol.R -i {{input}} -o {{output}} > {{log}} 2>&1
        """