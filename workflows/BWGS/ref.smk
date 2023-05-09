rule QC:
    input:
        R1="../sample/{sample}_R1.fastq",
        R2="../sample/{sample}_R2.fastq"
    output:
        R1="filterd/{sample}_R1.fastq",
        R2="filterd/{sample}_R2.fastq",
        rep="filterd/{sample}_R2.html"
    threads: 4
    shell:
        """
        fastp --in1 {input.R1} --in2 {input.R2} \
            --out1 {output.R1} --out2 {output.R2} \
            --thread {threads} -h {output.rep}
        """

rule index_ref:
    input:
        "../refs/Enterobacter_mori.fna.gz"
    output:
        "refs/bowtie2_index/Enterobacter_mori"
    threads: 4
    shell:
        "bowtie2-build --threads {threads} {input} {output}"

rule align:
    input:
        R1="filterd/{sample}_R1.fastq",
        R2="filterd/{sample}_R2.fastq",
        index="refs/bowtie2_index/Enterobacter_mori"
    output:
        "align/{sample}.sam"
    threads: 4
    shell:
        """
        bowtie2 --threads {threads} -x {input.index} \
                -1 {input.index} \
                -2 {input.R1} \
                -S {input.R2}
        """ 

rule fixmate:
    input: "align/{sample}.sam"
    output: "align/{sample}_fixmate.bam"
    shell: 
        """
        samtools sort -n \
                -O sam {input} | samtools fixmate \
                -m -O bam - {output}
        """

rule remove_duplicate:
    input: "align/{sample}_fixmate.bam"
    output: "align/{sample}_dedup.bam"
    shell: "samtools markdup -r -S {input} {output}"

rule mapping_stats:
    input: "align/{sample}_dedup.bam"
    output: 
        sam="align/{sample}_stats.txt",
        quali=directory("align/{sample}_stats")
    threads: 2
    shell:
        """
        samtools flagstat {input} > {output.sam} &
        qualimap bamqc -bam {input} -outdir {output.quali}
        wait
        """

rule extract_mapped:
    input: "align/{sample}_dedup.bam"
    output: 
        bam="align/{sample}_concordant.bam",
        R1="mapped_reads/{sample}_R1.fastq.gz",
        R2="mapped_reads/{sample}_R2.fastq.gz",
        U="mapped_reads/{sample}_U.fastq.gz"
    shell: 
        """
        samtools view -h -b -f 3 {input} > align/{sample}_concordant.bam
        samtools fastq -1 {output.R1} \
            -2 {output.R1} \
            -0 {output.U} align/{sample}_concordant.bam
        """

rule extract_unmapped:
    input: "align/{sample}_dedup.bam"
    output: 
        bam="align/{sample}_unmapped.bam",
        R1="unmapped_reads/{sample}_R1.fastq.gz",
        R2="unmapped_reads/{sample}_R2.fastq.gz",
        U="unmapped_reads/{sample}_U.fastq.gz"
    shell: 
        """
        samtools view -b -f 4 {input} > align/{sample}_unmapped.bam
        samtools fastq -1 {output.R1} \
            -2 {output.R1} \
            -0 {output.U} align/{sample}_unmapped.bam
        """

rule assemble_RAW_reads:
    input:
        R1="filterd/{sample}_R1.fastq",
        R2="filterd/{sample}_R2.fastq"
    output:
        eldir=directory("{sample}_filterd-Assembly"),
        asmpl="{sample}_filterd-Assembly/assembly.fasta"
    threads: 4
    shell:
        """
    unicycler -1 {input.R1} -2 {input.R2} \
              -o {output.eldir} \
              --threads {threads} --no_pilon --no_correct
        """

rule assmble_qc_RAW_reads:
    input:
        "{sample}_filterd-Assembly/assembly.fasta"
    output:
        directory("{sample}_filterd-Assembly/QC")
    shell:
        "quast -o {output} {input}"


rule assemble_mapped_reads:
    input:
        R1="mapped_reads/{sample}_R1.fastq.gz",
        R2="mapped_reads/{sample}_R2.fastq.gz"
    output:
        eldir=directory("{sample}_mapped_reads-Assembly"),
        asmpl="{sample}_mapped_reads-Assembly/assembly.fasta"
    threads: 4
    shell:
        """
    unicycler -1 {input.R1} -2 {input.R2} \
              -o {output.eldir} \
              --threads {threads} --no_pilon --no_correct
        """

rule assmble_qc_mapped_reads:
    input:
        "{sample}_mapped_reads-Assembly/assembly.fasta"
    output:
        directory("{sample}_mapped_reads-Assembly/QC")
    shell:
        "quast -o {output} {input}"


rule kraken_RAW_reads:
    input:
        R1="filterd/{sample}_R1.fastq",
        R2="filterd/{sample}_R2.fastq"
    output:
        report="kraken/filterd/{sample}.report",
        log="kraken/filterd/{sample}.out"
    params:
        db="/media/genomics/AlphaFold1/db/Kraken2/k2_pluspfp_20210127"
    threads: 4
    shell:
        """
        kraken2 --db {params.db} --threads {threads} \
            --report {output.report} \
            --paired {input.R1} {input.R2} > {output.log}
        """

rule kraken_unmapped_reads:
    input:
        R1="unmapped_reads/{sample}_R1.fastq.gz",
        R2="unmapped_reads/{sample}_R2.fastq.gz"
    output:
        report="kraken/unmapped_reads/{sample}.report",
        log="kraken/unmapped_reads/{sample}.out"
    params:
        db="/media/genomics/AlphaFold1/db/Kraken2/k2_pluspfp_20210127"
    threads: 4
    shell:
        """
        kraken2 --db {params.db} --threads {threads} \
            --report {output.report} \
            --paired {input.R1} {input.R2} > {output.log}
        """

rule index_ref_assembly:
    input:
        "{sample}_filterd-Assembly/assembly.fasta"
    output:
        "{sample}_filterd-Assembly/bowtie2_index/assembly"
    threads: 4
    shell:
        "bowtie2-build --threads {threads} {input} {output}"

rule align_assembly:
    input:
        R1="filterd/{sample}_R1.fastq",
        R2="filterd/{sample}_R2.fastq",
        index="{sample}_filterd-Assembly/bowtie2_index/assembly"
    output:
        "assembly-mapping/{sample}.sam"
    threads: 4
    shell:
        """
        bowtie2 --threads {threads} -x {input.index} \
                -1 {input.index} \
                -2 {input.R1} \
                -S {input.R2}
        """ 

rule fixmate_assembly:
    input: "assembly-mapping/{sample}.sam"
    output: "assembly-mapping/{sample}_fixmate.bam"
    shell: 
        """
        samtools sort -n \
                -O sam {input} | samtools fixmate \
                -m -O bam - {output}
        """

rule remove_duplicate_assembly:
    input: "assembly-mapping/{sample}_fixmate.bam"
    output: "assembly-mapping/{sample}_dedup.bam"
    shell: "samtools markdup -r -S {input} {output}"

rule mapping_stats_assembly:
    input: "assembly-mapping/{sample}_dedup.bam"
    output: 
        sam="assembly-mapping/{sample}_stats.txt",
        quali=directory("assembly-mapping/{sample}_stats")
    threads: 2
    shell:
        """
        samtools flagstat {input} > {output.sam} &
        qualimap bamqc -bam {input} -outdir {output.quali}
        wait
        """