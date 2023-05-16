
rule classify:
    input:
        get_classifeir_input
    output:
        "QIIME2/classify/taxonomy.qza"

    params:
        mydir=f"{GUAP_FOLDER}/bin/16s/snakemake",
        classifier=config["classifier"],
        classifier_threads=config["classifier_threads"],
        choose=config["choose_classifier"]

    threads: config["threads"]

    shell:
        """
        if [ {params.choose} == 'dada' ]; then 
            Rscript {params.mydir}/scripts/8_assign_taxa.R \
                -i {input} -o "${{PWD}}/DADA2" \
                --classifier {params.classifier} \
                -p {threads} -s

            qiime tools import \
                --type 'FeatureData[Taxonomy]' \
                --input-format HeaderlessTSVTaxonomyFormat \
                --input-path "${{PWD}}/DADA2/tax.txt" \
                --output-path {output}

        elif [ {params.choose} == 'qiime' ]; then 
            qiime feature-classifier classify-sklearn \
                --i-classifier {params.classifier} \
                --i-reads {input} --p-n-jobs {params.classifier_threads} \
                --o-classification {output} 
        
        else
            echo -e "ERROR"
            exit 1
        fi

        """
