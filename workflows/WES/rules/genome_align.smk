
rule bwa_align:
    input:
        unpack(get_align_input)

    output:
        "bwa/{sample}.sam"

    threads: config["threads_align"]
    params:
        index = ref_bwa,
        fa = ref_fasta

    log: "logs/{sample}_bwa.log"

    benchmark: "benchamrks/{sample}_bwa.txt"

    shell:
        """
        R1={input.R1}
        SM=$(basename $R1 | cut -d"_" -f1)   
        LB=$(basename $R1 | cut -d"_" -f1,2)  
        PL="Illumina"
        name=$(basename $R1 | cut -d'_' -f1)
        RGID=$(head -n1 $R1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
        PU=$RGID.$LB 
        bwa mem -t {threads} -M \
            -R "@RG\\tID:$RGID\\tSM:$SM\\tPL:$PL\\tLB:$LB\\tPU:$PU" {params.index} {input.R1} {input.R2} > {output} 2> {log}
        """

rule convert_sam:
    input: 
    	"bwa/{sample}.sam"
    
    output: 
    	"bwa/{sample}.bam"
    
    shell: 
    	"samtools view -hbo {output} {input}"


rule bowtie_align:
    input:
        unpack(get_align_input)

    output:
        "bowtie2/{sample}.noRG.sam"

    benchmark: "benchamrks/{sample}_bowtie2.txt"

    threads: config["threads_align"]
    params:
        index = ref_bowtie2

    log: "logs/{sample}_bowtie2.log"
    shell:
        """
        bowtie2 --threads {threads} \
            -x {params.index} -1 {input.R1} \
            -2 {input.R2} -S {output} > {log}
        """

rule add_read_group:
    input: "bowtie2/{sample}.noRG.sam"
    output: "bowtie2/{sample}.bam"
    shell:
        """
        PL="Illumina"

        picard AddOrReplaceReadGroups -I {input} -O bowtie2/{wildcards.sample}.sam \
            -SO coordinate -RGID {wildcards.sample} -LB {wildcards.sample} -PL $PL \
            -PU {wildcards.sample} -SM {wildcards.sample} \
            -CREATE_INDEX true
        
        samtools view -hbo {output} bowtie2/{wildcards.sample}.sam
        """

