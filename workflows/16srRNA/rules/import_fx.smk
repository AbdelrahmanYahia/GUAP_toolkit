
rule qiimeimport:
    input:
        "QIIME2/samples.tsv"
    output:
        reps="QIIME2/demux/paired-end-demux.qza",
        vis="QIIME2/visualization/demux.qzv"

    threads: 1
    shell:
        """
        qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path {input} \
        --output-path {output.reps} \
        --input-format PairedEndFastqManifestPhred33V2

        qiime demux summarize \
        --i-data QIIME2/demux/paired-end-demux.qza \
        --o-visualization {output.vis} 
        """

rule QIIME_Quality_filter:
    input:
        "QIIME2/demux/paired-end-demux.qza"
    output:
        seqs="QIIME2/demux/demux-filtered.qza",
        stats="QIIME2/demux/demux-filter-stats.qza",
        viz="QIIME2/visualization/demux-filter-stats.qzv"
    threads: 1
    shell:
        """
        qiime quality-filter q-score \
            --i-demux {input} \
            --o-filtered-sequences {output.seqs} \
            --o-filter-stats {output.stats}

        qiime metadata tabulate \
            --m-input-file QIIME2/demux/demux-filter-stats.qza \
            --o-visualization {output.viz}
        """

rule qiimeimportdada:
    input:
        seqs=rules.dada2_infer.output.seqs,
        seqt=rules.dada2_infer.output.seqt

    output:
        repseqs="QIIME2/r-rep-seqs.qza",
        table="QIIME2/r-table.qza"

    threads: 1
    shell:
        """
        qiime tools import \
            --input-path {input.seqs} \
            --type "FeatureData[Sequence]" \
            --output-path {output.repseqs}
        
        cd QIIME2/

        echo -n "#OTU Table" | cat - ../DADA2/seqtab-nochim.txt > biom-table.txt
        biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5

        qiime tools import \
            --input-path table.biom \
            --type "FeatureTable[Frequency]" \
            --input-format BIOMV210Format \
            --output-path r-table.qza
        
        rm biom-table.txt
        
        cd ../../
        """

