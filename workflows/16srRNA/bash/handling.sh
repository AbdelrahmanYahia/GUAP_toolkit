
export_phyloseq(){
	echo -e "${YEL}#######    Exporting phyloseq object   #####${NC}"
	Rscript ${GUAP_DIR}/bin/16s/scripts/export_phyloseq.R \
		-w "${OUTPUT}" \
		-s "${OUTPUT}/QIIME2/table.qza"\
		-m ${metadata} \
		-t "${OUTPUT}/QIIME2/align/rooted-tree.qza" \
		-x "${OUTPUT}/QIIME2/classify/taxonomy.qza" \
		-o "${OUTPUT}/ps.rds"
}


summrize(){
	echo -e "${YEL}#############     summarize     ############${NC}"
	
	qiime feature-table summarize \
	--i-table ./table.qza \
	--o-visualization visualization/table.qzv \
	--m-sample-metadata-file ${metadata} &

	qiime feature-table tabulate-seqs \
	--i-data ./rep-seqs.qza \
	--o-visualization visualization/rep-seqs.qzv
	wait
}

final_handling(){
	echo -e "${YEL}#############    biom to csv    ############${NC}"
	# Convert biom format to tsv format
	biom convert \
		-i ./phyloseq/feature-table.biom \
		-o ./phyloseq/otu_table.tsv \
		--to-tsv
    error_cheker $? 1 "biom to csv"


	cd phyloseq
	sed -i "1d" otu_table.tsv
	sed -i "s/#OTU ID//" otu_table.tsv
	cd ../
    error_cheker $? 1 "modify biom file "


	echo -e "${YEL}###### Export representative sequences #####${NC}"
	# Export representative sequences
	qiime tools export \
		--input-path rep-seqs.qza \
		--output-path phyloseq

    error_cheker $? 1 "export ASVs"


	echo -e "${YEL}##########     export tree files    ########${NC}"
	# Export tree files
	qiime tools export \
		--input-path align/unrooted-tree.qza \
		--output-path phyloseq

    error_cheker $? 0 "exporting tree files"

	cd phyloseq
	mv tree.nwk unrooted_tree.nwk
	cd ../


	echo -e "${YEL}############     rooted tree    ############${NC}"
	qiime tools export \
		--input-path align/rooted-tree.qza \
		--output-path phyloseq

	cd phyloseq
	mv tree.nwk rooted_tree.nwk
	cd ../
}
