
align_(){
	echo -e "${YEL}###############     align     ##############${NC}"

	qiime phylogeny align-to-tree-mafft-fasttree \
		--i-sequences rep-seqs.qza \
		--o-alignment align/aligned-rep-seqs.qza \
		--o-masked-alignment align/masked-aligned-rep-seqs.qza \
		--o-tree align/unrooted-tree.qza \
		--o-rooted-tree align/rooted-tree.qza --p-n-threads ${threads}
}


taxa_plot(){
	echo -e "${YEL}#############     taxa plot     ############${NC}"
	
	qiime taxa barplot \
		--i-table table.qza \
		--i-taxonomy classify/taxonomy.qza \
		--m-metadata-file ${metadata} \
		--o-visualization visualization/bar_plot.qzv 
	error_cheker $? 1 "bar plot"
	qiime tools export \
	--input-path visualization/bar_plot.qzv \
	--output-path exports/barplot
	error_cheker $? 1 "export bar plot"

	qiime tools export --input-path classify/taxonomy.qza --output-path "./classify" 

}

export_biom(){
	echo -e "${YEL}###########    export biom     #############${NC}"

	# Export OTU table
	qiime tools export \
		--input-path table.qza \
		--output-path phyloseq 
}


core_metrics(){
	echo -e "${YEL}########  core-metrics-phylogenetic  #######${NC}"
	qiime diversity core-metrics-phylogenetic \
	--i-phylogeny align/rooted-tree.qza \
	--i-table table.qza \
	--p-sampling-depth ${sampling_depth} \
	--m-metadata-file ${metadata} \
	--p-n-jobs-or-threads ${threads} \
	--output-dir core-metrics-results
	error_cheker $? 0 "core metrics phylogenetic"

	for i in core-metrics-results/*.qzv ; do
		mv $i ../visualization ; done
}

alpha_rarefaction(){
	echo -e "${YEL}#######  diversity alpha-rarefaction  ######${NC}"
	qiime diversity alpha-rarefaction \
	--i-table table.qza \
	--i-phylogeny align/rooted-tree.qza \
	--m-metadata-file ${metadata} \
	--p-max-depth ${max_depth} \
	--o-visualization visualization/alpha-rarefaction.qzv
	error_cheker $? 1 "diversity alpha-rarefaction"
	qiime tools export  \
	--input-path visualization/alpha-rarefaction.qzv \
	--output-path exports/alpha-rarefaction 
}

shannon_vector(){
	echo -e "${YEL}#######  shannon_vector-group-significance  ######${NC}"

	qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/shannon_vector.qza \
	--m-metadata-file ${metadata} \
	--o-visualization visualization/shannon_vector-group-significance.qzv
	error_cheker $? 1 "shannon_vector group significance"
	qiime tools export  \
	--input-path visualization/shannon_vector-group-significance.qzv \
	--output-path exports/shannon_vector-group-significance 
}

alpha_group_significace(){
	echo -e "${YEL}########  alpha-group-significance  ########${NC}"
	qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
	--m-metadata-file ${metadata} \
	--o-visualization visualization/faith-pd-group-significance.qzv 
	error_cheker $? 1 "alpha group siginifcance"
	qiime tools export  \
	--input-path visualization/faith-pd-group-significance.qzv \
	--output-path exports/faith-pd-group-significance 

}

evenness_group(){
	echo -e "${YEL}########  evenness-group-significance  ######${NC}"
	qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/evenness_vector.qza \
	--m-metadata-file ${metadata} \
	--o-visualization visualization/evenness-group-significance.qzv 
	error_cheker $? 1 "evennes group significance"
	qiime tools export  \
	--input-path visualization/evenness-group-significance.qzv \
	--output-path exports/evenness-group-significance 
	error_cheker $?	1 "exporting evennesss group"
}

beta_group_significance(){
	echo -e "${YEL}##### unweighted-unifrac-subject-group #####${NC}"
	qiime diversity beta-group-significance \
	--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	--m-metadata-file ${metadata} \
	--m-metadata-column ${condition_name} \
	--o-visualization visualization/unweighted-unifrac-condition-significance.qzv \
	--p-pairwise 

	qiime diversity beta-group-significance \
	--i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
	--m-metadata-file ${metadata} \
	--m-metadata-column ${condition_name} \
	--o-visualization visualization/weighted-unifrac-condition-significance.qzv \
	--p-pairwise 
}

plot_beta_group(){
	echo -e "${YEL}###########      emperor plot     ##########${NC}"
	qiime emperor plot \
	--i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
	--m-metadata-file ${metadata} \
	--p-custom-axes \
	--o-visualization visualization/unweighted-unifrac-emperor.qzv 
	error_cheker $? 1 "emperor plot"

	qiime emperor plot \
	--i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
	--m-metadata-file ${metadata} \
	--o-visualization visualization/bray-curtis-emperor.qzv 
}

performe_downstream(){
	mkdir -p QIIME2/exports
	get_sampling_depth
	step_checker 0 core_metrics
	step_checker 1 alpha_rarefaction
	step_checker 1 shannon_vector
	step_checker 1 alpha_group_significace
	step_checker 1 evenness_group
	step_checker 1 plot_beta_group
	step_checker 1 beta_group_significance
}

