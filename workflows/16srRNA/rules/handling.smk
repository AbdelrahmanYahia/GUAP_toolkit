CDIR = get_cutadapt_input()
PREF = get_rename_prefix()
AID = get_analysis_input_dir()



rule rename:
    input:
        seqs=f"QIIME2/{PREF}-rep-seqs.qza",
        tab=f"QIIME2/{PREF}-table.qza"
    output:
        seqs="QIIME2/rep-seqs.qza",
        tab="QIIME2/table.qza"
    shell:
        """
        mv {input.seqs} {output.seqs}
        mv {input.tab} {output.tab}
        """


rule remove_primers:
    input:
        nf1=f"{CDIR}/{{sample}}_R1{EXT}",
        nf2=f"{CDIR}/{{sample}}_R2{EXT}"
    output:
        of1=f"cutadapt/{{sample}}_R1{EXT}",
        of2=f"cutadapt/{{sample}}_R2{EXT}"
    log: 
        "logs/cutadapt/{sample}.txt"
    threads: USE_THREADS
    params:
        forprimer=config["forward_primer"],
        revprimer=config["reverse_primer"],
        m=config["min_length"]

    shell:
        """
        cutadapt -j {threads} -a ^{params.forprimer}... \
        -A ^{params.revprimer}... \
        -m {params.m} --discard-untrimmed \
        -o {output.of1} \
        -p {output.of2} \
        {input.nf1} {input.nf2}  >> {log} 2>&1
        """


rule get_qiime_s_table:
    input: 
        get_analysis_input(),
        mydir=f"{GUAP_FOLDER}/bin/16s/"
    params:
        indir=f"{AID}"
    output:
        "QIIME2/samples.tsv"
    shell:
        """
        Rscript {input.mydir}/scripts/generate_sample_table.R \
                    -i {params.indir} \
                    --r1-pattern '_R1.fastq' \
		            --r2-pattern '_R2.fastq' \
                    -o QIIME2
        """

rule getreads:
    output:
        R1 = (f"reads/{{sample}}_R1{EXTT}"),
        R2 = (f"reads/{{sample}}_R2{EXTT}")
    params: 
        input_path= PATH
    run:
        shell("cp {params.input_path}/{wildcards.sample}_*_{RS}1{TAIL}{EXTT} {output.R1}"),
        shell("cp {params.input_path}/{wildcards.sample}_*_{RS}2{TAIL}{EXTT} {output.R2}"),
        shell("mkdir -p QIIME2"),
        shell("mkdir -p QIIME2/phyloseq"),
        shell("mkdir -p QIIME2/classify"),
        shell("mkdir -p QIIME2/visualization"),
        shell("mkdir -p QIIME2/align")

rule decompress:
    input: 
        f"reads/{{sample}}_R{{R}}{EXTT}"
    output:
        temp(f"reads/{{sample}}_R{{R}}{EXT}")
    shell:
        "gunzip {input}"


rule cleaning:
    input:
        get_final_output
    output:
        temp("log.txt")
    params:
        dir=f"{AID}"
    shell:
        """
        rm -r {params.dir}/filtered >> log.txt
        rm -r reads >> log.txt
        echo -e "\033[;34;1mALL DONE\033[;39;m"
        """
