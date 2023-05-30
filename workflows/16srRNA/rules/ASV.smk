
rule dada2_infer:
    input:
        get_analysis_input()
    output:
        seqs="DADA2/rep-seqs.fna",
        seqt="DADA2/seqtab-nochim.txt"
    params:
        indir=f"{AID}",
        mydir=f"{GUAP_FOLDER}/workflows/16srRNA",
        t=config["trunc_f"],
        T=config["trunc_r"],
        l=config["trim_l"],
        r=config["trim_r"],
        e=config["maxee_f"],
        E=config["maxee_r"],
        n=config["name"],
        m=config["min_overlap"],
        chim=config["chimera_method"],
        R1=R1_pattern,
        R2=R2_pattern

    threads: config["threads"]
    shell:
        """
            Rscript {params.mydir}/scripts/G16s.v0.9.R \
                -i {params.indir} -o "${{PWD}}/DADA2" \
                -p {threads} \
                -n '{params.n}' \
                -t {params.t} \
                -T {params.T} \
                -l {params.l} \
                -r {params.r} \
                -e {params.e} \
                --r1-pattern {params.R1} \
		        --r2-pattern {params.R2} \
                -m {params.m} --chimethod {params.chim} \
                -E {params.E} -s --use-exsisting
        """

rule QIIME_dada2_infer:
    input:
        "QIIME2/demux/paired-end-demux.qza"
    output:
        seqs="QIIME2/dada-rep-seqs.qza",
        table="QIIME2/dada-table.qza",
        stats="QIIME2/table-stats.qza"

    params:
        t=config["trunc_f"],
        T=config["trunc_r"],
        l=config["trim_l"],
        r=config["trim_r"],
        e=config["maxee_f"],
        E=config["maxee_r"],        
        chim=config["chimera_method"],
        m=config["min_overlap"]

    threads: config["threads"]
    shell:
        """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input} \
            --o-representative-sequences {output.seqs} \
            --o-table {output.table} \
            --o-denoising-stats {output.stats} \
            --p-n-threads {threads} \
            --p-trunc-len-f {params.t} \
            --p-trunc-len-r {params.T} \
            --p-trim-left-f {params.l} \
            --p-trim-left-r {params.r} \
            --p-max-ee-f {params.e} \
            --p-max-ee-r {params.E} \
            --p-chimera-method {params.chim} \
            --p-min-overlap {params.m}
        
        qiime tools export \
            --input-path QIIME2/table-stats.qza \
            --output-path stats
            
        qiime metadata tabulate \
            --m-input-file QIIME2/table-stats.qza \
            --o-visualization QIIME2/visualization/table-stats.qzv
        
        """


rule QIIME_deblur_infer:
    input:
        inp="QIIME2/demux/demux-filtered.qza"
    output:
        seqs="QIIME2/deblur-rep-seqs.qza",
        table="QIIME2/deblur-table.qza",
        stats="QIIME2/deblur/deblur-stats.qza",
        viz="QIIME2/visualization/deblur-stats.qzv"

    params:
        l=config["deblur_trim_length"]

    threads: config["threads"]
    shell:
        """
        qiime deblur denoise-16S --p-jobs-to-start {threads} \
            --i-demultiplexed-seqs {input.inp} \
            --p-trim-length {params.l} \
            --o-representative-sequences {output.seqs} \
            --o-table {output.table} \
            --p-sample-stats \
            --o-stats {output.stats}

        qiime tools export \
            --input-path QIIME2/deblur/deblur-stats.qza \
            --output-path stats &

        qiime deblur visualize-stats \
            --i-deblur-stats QIIME2/deblur/deblur-stats.qza \
            --o-visualization {output.viz} &
        
        wait
        """
