PREF = get_rename_prefix()

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
        unpack(get_align_input)
    output:
        of1 = f"cutadapt/{{sample}}_{RS}1.{EXT}",
        of2 = f"cutadapt/{{sample}}_{RS}2.{EXT}"
        
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
        {input.R1} {input.R2}  >> {log} 2>&1
        """

rule get_qiime_s_table:
    input: 
        get_analysis_input(),
        mydir=f"{GUAP_FOLDER}/workflows/16srRNA"
    params:
        indir=f"{AID}",
        R1=R1_pattern,
        R2=R2_pattern

    output:
        "QIIME2/samples.tsv"
    shell:
        """
        Rscript {input.mydir}/scripts/generate_sample_table.R \
                    -i {params.indir} \
                    --r1-pattern {params.R1} \
		            --r2-pattern {params.R2} \
                    -o QIIME2
        """


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
