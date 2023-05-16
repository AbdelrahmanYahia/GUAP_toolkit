
rule align:
    input:
        "QIIME2/rep-seqs.qza"
    output:
        aligned="QIIME2/align/aligned-rep-seqs.qza",
        masked="QIIME2/align/masked-aligned-rep-seqs.qza",
        tree="QIIME2/align/unrooted-tree.qza",
        rooted="QIIME2/align/rooted-tree.qza"
    threads: USE_THREADS_QIIME
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input} \
            --o-alignment {output.aligned} \
            --o-masked-alignment {output.masked} \
            --o-tree {output.tree} \
            --o-rooted-tree {output.rooted} --p-n-threads {threads}
        """

rule export_phyloseq:
    input:
        table="QIIME2/table.qza",
        meta=config["metadata"],
        taxon="QIIME2/classify/taxonomy.qza",
        rooted="QIIME2/align/rooted-tree.qza"
    output:
        "ps.rds"
    threads: 1
    shell:
        """
        Rscript ~/GUAP/bin/16s/scripts/export_phyloseq.R \
            -w ${{PWD}} \
            -s {input.table}\
            -m {input.meta} \
            -t {input.rooted} \
            -x {input.taxon} \
            -o "{output}"
        """


rule summarize:
    input:
        table="QIIME2/table.qza",
        meta=config["metadata"],
        repseqs="QIIME2/rep-seqs.qza"

    output:
        table="QIIME2/visualization/table.qzv",
        seqs="QIIME2/visualization/rep-seqs.qzv"
    threads: 2
    shell:
        """
        qiime feature-table summarize \
            --i-table {input.table} \
            --o-visualization {output.table} \
            --m-sample-metadata-file {input.meta} &

        qiime feature-table tabulate-seqs \
            --i-data {input.repseqs} \
            --o-visualization {output.seqs}

        wait
        """

rule plot_bar:
    input:
        table="QIIME2/table.qza",
        meta=config["metadata"],
        taxon="QIIME2/classify/taxonomy.qza"

    output:
        "QIIME2/visualization/bar_plot.qzv"

    threads: 1
    shell:
        """
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxon} \
            --m-metadata-file {input.meta} \
            --o-visualization {output}

        qiime tools export --input-path {input.taxon} --output-path "./QIIME2/classify" 

        """


rule tabulate:
    input:
        "QIIME2/classify/taxonomy.qza"

    output:
        "QIIME2/visualization/taxonomy.qzv"

    threads: 1
    shell:
        """
        qiime metadata tabulate \
            --m-input-file {input} \
            --o-visualization {output}
        """

rule metrics:
    input:
        meta=config["metadata"],
        table="QIIME2/table.qza",
        rot="QIIME2/align/rooted-tree.qza"

    output:
        directory("QIIME2/core-metrics-results/")

    params:
        SD=config["sampling_depth"],
        DL=config["deblur"],
        QI=config["use_QIIME2"]
    
    threads: USE_THREADS_QIIME
    shell:
        """
        sampling_depth={params.SD}
        deblur={params.DL}
        use_QIIME2={params.QI}

        if [[ ${{sampling_depth}} == "None" || "null" ]] ; then
            if [ "${{deblur}}" == "True" ] || [ "${{deblur}}" == "true" ] ; then 
                sampling_depth=`cat stats/stats.csv | rev | cut -d "," -f3 | rev | tail -n +2 | sort -n | head -1`
                echo -e "\033[;33;1mWARNING: \033[;39;mUsing minimum Depth as sampling depth; \033[;33;1m${{sampling_depth}}\033[;39;m"
            elif [ "${{use_QIIME2}}" == "True" ] || [ "${{use_QIIME2}}" == "true" ] ; then 
                sampling_depth=`cat stats/stats.tsv | rev | cut  -f2 | rev | tail -n +3 | sort -n | head -1`
                echo -e "\033[;33;1mWARNING: \033[;39;mUsing minimum Depth as sampling depth; \033[;33;1m${{sampling_depth}}\033[;39;m"
            elif [ -d "DADA2" ] ; then
                sampling_depth=`cat ./DADA2/stats.csv | rev | cut -d "," -f1 | rev | tail -n +2 | sort -n | head -1`
                echo -e "\033[;33;1mWARNING: \033[;39;mUsing minimum Depth as sampling depth; \033[;33;1m${{sampling_depth}}\033[;39;m"
            else 
                echo -e "\033[;31;1mERROR: \033[;39;msomething worng! check params"
                exit 1
            fi
        fi
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny {input.rot} \
            --i-table {input.table} \
            --p-sampling-depth ${{sampling_depth}} \
            --m-metadata-file {input.meta} \
            --p-n-jobs-or-threads {threads} \
            --output-dir QIIME2/core-metrics-results

        cd QIIME2/core-metrics-results
        for i in *.qzv
        do
            mv $i ../visualization
        done
        cd ../../../

        """

rule export_tree_files:
    input:
        un="QIIME2/align/unrooted-tree.qza",
        rot="QIIME2/align/rooted-tree.qza"

    output:
        "QIIME2/unrooted-tree.nwk",
        "QIIME2/rooted-tree.nwk"

    threads: 1
    shell:
        """
        qiime tools export \
            --input-path {input.un} \
            --output-path QIIME2

        mv QIIME2/tree.nwk QIIME2/unrooted-tree.nwk

        qiime tools export \
            --input-path {input.rot} \
            --output-path QIIME2

        mv QIIME2/tree.nwk QIIME2/rooted-tree.nwk
        """
rule alpha_rarefaction:
    input:
        table="QIIME2/table.qza",
        meta=config["metadata"],
        rot="QIIME2/align/rooted-tree.qza"

    output:
        rare="QIIME2/visualization/alpha-rarefaction.qzv"        
    params:
        MD=config["max_depth"],
        DL=config["deblur"],
        QI=config["use_QIIME2"]
    
    threads: 1
    shell:
        """
        max_depth={params.MD}
        deblur={params.DL}
        use_QIIME2={params.QI}

        if [[ ${{max_depth}} == "None" || "null" ]] ; then
            if [ "${{deblur}}" == "True" ] || [ "${{deblur}}" == "true" ] ; then 
                max_depth=`cat stats/stats.csv | rev | cut -d "," -f3 | rev | tail -n +2 | sort -n | tail -1`
                echo -e "\033[;33;1mWARNING: \033[;39;mUsing maximum Depth as max-depth; \033[;33;1m${{max_depth}}\033[;39;m"
            elif [ "${{use_QIIME2}}" == "True" ] || [ "${{use_QIIME2}}" == "true" ] ; then 
                max_depth=`cat stats/stats.tsv | rev | cut  -f2 | rev | tail -n +3 | sort -n | tail -1`
                echo -e "\033[;33;1mWARNING: \033[;39;mUsing maximum Depth as max-depth; \033[;33;1m${{max_depth}}\033[;39;m"
            elif [ -d "DADA2" ] ; then
                max_depth=`cat ./DADA2/stats.csv | rev | cut -d "," -f1 | rev | tail -n +2 | sort -n | tail -1`
                echo -e "\033[;33;1mWARNING: \033[;39;mUsing maximum Depth as max-depth; \033[;33;1m${{max_depth}}\033[;39;m"
            else 
                echo -e "\033[;31;1mERROR: \033[;39;msomething worng! check params"
                exit 1
            fi
        fi
        qiime diversity alpha-rarefaction \
            --i-table {input.table} \
            --i-phylogeny {input.rot} \
            --m-metadata-file {input.meta} \
            --p-max-depth ${{max_depth}} \
            --o-visualization {output.rare}

        qiime tools export  \
            --input-path QIIME2/visualization/alpha-rarefaction.qzv \
            --output-path QIIME2/exports/alpha-rarefaction
        """

rule alpha_group_significace:
    input:
        rules.metrics.output
    params:
        inp="QIIME2/core-metrics-results/faith_pd_vector.qza",
        meta=config["metadata"]

    output:
        "QIIME2/visualization/faith-pd-group-significance.qzv"

    threads: 1
    shell:
        """
        qiime diversity alpha-group-significance \
        --i-alpha-diversity {params.inp} \
        --m-metadata-file {params.meta} \
        --o-visualization {output}

        qiime tools export  \
        --input-path QIIME2/visualization/faith-pd-group-significance.qzv \
        --output-path QIIME2/exports/faith-pd-group-significance 
        """

rule shannon_vector:
    input:
        rules.metrics.output
    params:
        inp="QIIME2/core-metrics-results/shannon_vector.qza",
        meta=config["metadata"]

    output:
        "QIIME2/visualization/shannon_vector-group-significance.qzv"

    threads: 1
    shell:
        """
        qiime diversity alpha-group-significance \
        --i-alpha-diversity {params.inp} \
        --m-metadata-file {params.meta} \
        --o-visualization {output}

        qiime tools export  \
        --input-path QIIME2/visualization/shannon_vector-group-significance.qzv \
        --output-path QIIME2/exports/shannon_vector-group-significance 
        """

rule evenness_group:
    input:
        rules.metrics.output
    params:
        inp="QIIME2/core-metrics-results/evenness_vector.qza",
        meta=config["metadata"]

    output:
        "QIIME2/visualization/evenness-group-significance.qzv"

    threads: 1
    shell:
        """
        qiime diversity alpha-group-significance \
        --i-alpha-diversity {params.inp} \
        --m-metadata-file {params.meta} \
        --o-visualization {output}

        qiime tools export  \
        --input-path QIIME2/visualization/evenness-group-significance.qzv \
        --output-path QIIME2/exports/evenness-group-significance 
        """

rule beta_group_significance:
    input:
        rules.metrics.output
    params:
        inp1="QIIME2/core-metrics-results/unweighted_unifrac_distance_matrix.qza",
        inp2="QIIME2/core-metrics-results/weighted_unifrac_distance_matrix.qza",
        meta=config["metadata"],
        cond=config["condition_name"]

    output:
        out1="QIIME2/visualization/unweighted-unifrac-condition-significance.qzv",
        out2="QIIME2/visualization/weighted-unifrac-condition-significance.qzv"
        
    threads: 2
    shell:
        """
    	qiime diversity beta-group-significance \
    	--i-distance-matrix {params.inp1} \
    	--m-metadata-file {params.meta} \
    	--m-metadata-column {params.cond} \
    	--o-visualization {output.out1} \
    	--p-pairwise &

    	qiime diversity beta-group-significance \
    	--i-distance-matrix {params.inp2} \
    	--m-metadata-file {params.meta} \
    	--m-metadata-column {params.cond} \
    	--o-visualization {output.out2} \
    	--p-pairwise 

        wait
        """

rule plot_beta_group:
    input:
        rules.metrics.output
    params:
        inp1="QIIME2/core-metrics-results/unweighted_unifrac_pcoa_results.qza",
        inp2="QIIME2/core-metrics-results/bray_curtis_pcoa_results.qza",
        meta=config["metadata"],

    output:
        viz1="QIIME2/visualization/unweighted-unifrac-emperor.qzv",
        viz2="QIIME2/visualization/bray-curtis-emperor.qzv"

    threads: 2
    shell:
        """
        qiime emperor plot \
        --i-pcoa {params.inp1} \
        --m-metadata-file {params.meta} \
        --p-custom-axes \
        --o-visualization {output.viz1} &

        qiime emperor plot \
        --i-pcoa {params.inp2} \
        --m-metadata-file {params.meta} \
        --o-visualization {output.viz2} 

        wait
        """



