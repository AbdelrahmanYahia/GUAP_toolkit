qiime_import_direct(){
	echo -e "${YEL}###########      QIIME import    ###########${NC}"
	qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path samples.tsv \
	--output-path demux/paired-end-demux.qza \
	--input-format PairedEndFastqManifestPhred33V2
	error_cheker $? 0 "QIIME import"

	# demultiplexing
	qiime demux summarize \
	--i-data demux/paired-end-demux.qza \
	--o-visualization visualization/demux.qzv 
}

qiime_quality_filter(){
	echo -e "${YEL}###########    Quality filter    ###########${NC}"
	qiime quality-filter q-score \
		--i-demux demux/paired-end-demux.qza \
		--o-filtered-sequences demux/demux-filtered.qza \
		--o-filter-stats demux/demux-filter-stats.qza
}


qiime2_import_dada(){
	echo -e "${YEL}###########      QIIME import    ###########${NC}"
	qiime tools import \
		--input-path ../DADA2/rep-seqs.fna \
		--type "FeatureData[Sequence]" \
		--output-path rep-seqs.qza
	error_cheker $? 0 "qiime2 import dada"

	echo -n "#OTU Table" | cat - ../DADA2/seqtab-nochim.txt > biom-table.txt
	error_cheker $? 0 "export biom-table"

	biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5
	error_cheker $? 0 "convert biom table"

	qiime tools import \
	--input-path table.biom \
	--type "FeatureTable[Frequency]" \
	--input-format BIOMV210Format \
	--output-path table.qza
	error_cheker $? 0 "import dada biom to qiime2"

	if [ -f "${working_dir}/DADA2/tax.txt" ]; then 
	echo -e "${YEL}###########      QIIME importing classification    ###########${NC}"

		qiime tools import \
			--type 'FeatureData[Taxonomy]' \
			--input-format HeaderlessTSVTaxonomyFormat \
			--input-path "${working_dir}/DADA2/tax.txt" \
			--output-path classify/taxonomy.qza
		error_cheker $? 0 "import calssification of dada2 in qiime2"
	fi
	
}
