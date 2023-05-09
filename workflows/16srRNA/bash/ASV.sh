
deblur_infer(){
	echo -e "${YEL}################     ASV     ###############${NC}"
	qiime deblur denoise-16S \
		--i-demultiplexed-seqs $1 \
		--p-trim-length ${deblur_trim_length} \
		--o-representative-sequences rep-seqs.qza \
		--o-table table.qza \
		--p-sample-stats \
		--p-jobs-to-start ${threads} \
		--o-stats deblur/deblur-stats.qza
	error_cheker $? 0 "deblur denoise"

	qiime tools export \
		--input-path deblur/deblur-stats.qza \
		--output-path stats

	qiime metadata tabulate \
		--m-input-file demux/demux-filter-stats.qza \
		--o-visualization visualization/demux-filter-stats.qzv 

	qiime deblur visualize-stats \
		--i-deblur-stats deblur/deblur-stats.qza \
		--o-visualization visualization/deblur-stats.qzv 
	
}

dada2_infer(){
	echo -e "${YEL}################     ASV     ###############${NC}"
	qiime dada2 denoise-paired \
	--p-n-threads ${threads} \
	--i-demultiplexed-seqs $1 \
	--p-trunc-len-r ${trunc_r} \
	--p-trunc-len-f ${trunc_f} \
	--p-trim-left-r ${trim_r} \
	--p-trim-left-f ${trim_l} \
	--p-max-ee-f ${maxee_f} \
	--p-max-ee-r ${maxee_r} \
	--o-representative-sequences rep-seqs.qza \
	--o-table table.qza \
	--o-denoising-stats table-stats.qza
	error_cheker $? 0 "DADA2 deniose"

	qiime tools export \
		--input-path table-stats.qza \
		--output-path stats

	qiime metadata tabulate \
	--m-input-file table-stats.qza \
	--o-visualization visualization/table-stats.qzv 

}

dada2_R(){
	if [[ ${choose_classifier} == 'dada' ]]; then
		classifydada="--classifier ${classifier} --assign-taxa"
	else
		classifydada=""
	fi

	echo -e "${YEL}########    DADA2 Rscript running   ########${NC}"
	if [ "${continue}" == "True" ] || [ "${continue}" == "true" ] ; then
		dadacont="--use-exsisting"
	else
		dadacont=""
	fi

	Rscript ${GUAP_DIR}/bin/16s/scripts/G16s.v0.9.R \
		-i ${INPUT} -o "${PWD}/DADA2"\
		-p ${threads} \
		-n ${name} \
		-t ${trunc_f} \
		-T ${trunc_r} \
		-l ${trim_l} \
		-r ${trim_r} \
		-e ${maxee_f} \
		--r1-pattern ${R1_pattern} \
		--r2-pattern ${R2_pattern} ${dadacont} \
		-m ${min_overlap} --chimethod ${chimera_method} \
		-E ${maxee_r} GUAP 16s -i indir -o outdir -m metadata.tsv -c classifier.qza \
         -tf 220 -tr 170 -l 10 -r 10 -ef 5 -er 7 --downstream \
         --export-figs --condition-name condition
