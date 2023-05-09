
train_classifier(){
	qiime feature-classifier extract-reads \
		--i-sequences  ${t_i_seqs}\
		--p-f-primer ${forward_primer} \
		--p-r-primer ${reverse_primer} \
		--p-min-length ${trainset_minlength} \
		--p-max-length ${trainset_maxlength} \
		--o-reads ref-seqs.qza --p-n-jobs {threads}
	error_cheker $? 0 "extract reads of classifier"
	    # Train classifier
	qiime feature-classifier fit-classifier-naive-bayes \
		--i-reference-reads ref-seqs.qza \
		--i-reference-taxonomy ${t_i_taxa} \
		--o-classifier Trained_classifier.qza \
	error_cheker $? 0 "fit classifier"
	classifier="Trained_classifier.qza"

}

classify_(){
	echo -e "${YEL}#############      classify     ############${NC}"
	qiime feature-classifier classify-sklearn \
		--i-classifier ${classifier} \
		--i-reads rep-seqs.qza --p-n-jobs ${classifier_threads} \
		--o-classification classify/taxonomy.qza
	error_cheker $? 0 "classifying "
	
	qiime metadata tabulate \
		--m-input-file classify/taxonomy.qza \
		--o-visualization visualization/taxonomy.qzv 
}

