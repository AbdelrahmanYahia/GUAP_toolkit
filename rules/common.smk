sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file,index_col=0)
SAMPLES = list(SampleTable.index)

rule Fastqc:
	input:
		rawread="Data/samples/{sample}.fastq"
	output:
		zip="OUT/rawQC/{sample}_fastqc.zip",
		html="OUT/rawQC/{sample}_fastqc.html"
	threads:1
	params:
		path="OUT/rawQC/"
	shell:
		"fastqc {input.rawread} --threads {threads} -o {params.path}"

rule cutadapt1:
    input:
        "Data/samples/{sample}.fastq"
    output:
        temp("OUT/cutadapt_temp/{sample}.fastq")
    params:
        adaptor=config["trimming"]["adaptor"]
    
    threads:1
    shell:
        "cutadapt -j {threads} -a {params.adaptor} -o {output} --minimum-length 23 {input}"

rule cutadapt:
    input:
        "OUT/cutadapt_temp/{sample}.fastq"
    output:
        "OUT/cutadapt/{sample}.fastq"        
    threads:1
    shell:
        "cutadapt -j {threads} --minimum-length=12 --maximum-length=35 -u 4 -u -4 -o {output} {input}"

rule bowtie_aling:
    input:
       "OUT/cutadapt/{sample}.fastq"
    output:
        "OUT/mapped_reads/{sample}.sam"
    params:
        genome=config["ref"]["hsa_genome"]
        
    threads: 4
    shell:
        "bowtie --threads {threads} -q -v 0 -k 10 -S -t {params.genome} {input}  -S {output}"

rule Counts:
    input:
        insure_output = expand(rules.bowtie_aling.output,sample=SAMPLES)
        # allsams = glob_wildcards("mapped_reads/{sample}.sam")
    output:
        "OUT/counts/All.txt"
    params:
        anno=config["ref"]["miRBase_annotation"]
    threads:8
    shell:
        "featureCounts -t miRNA -g Name -O -T {threads} -s 1 -M -a {params.anno} -o {output} {input}"

