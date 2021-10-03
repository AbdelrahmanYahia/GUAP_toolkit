#!/bin/bash
take_args(){
	read -p $'\x1b[1;31mmodify config params (y/n)?\e[0m' CHOICE
	case "$CHOICE" in 
	y|Y ) 
		read -e -p $'\x1b[1;34mTrain Classifier(y/n)? : \e[0m' trainset_train
		if [ "${trainset_train}" == "TRUE" ] || [ "${trainset_train}" == "True" ] || [ "${trainset_train}" == "T" ] || [ "${trainset_train}" == "true" ] 
		then
			read -p "Train set forward primer : " trainset_fprimer
			read -p "Train set reverse primer : " trainset_rprimer
			read -p "Train set min length : " trainset_minlength
			read -p "Train set max length : " trainset_maxlength
		else
			read -e -p "Classifier path: " indexes_classifier
		fi
		read -p "Use Deblur : (TRUE/FALSE)" deblur_use
		if [ "${deblur_use}" == "TRUE" ] || [ "${deblur_use}" == "True" ] || [ "${deblur_use}" == "T" ] || [ "${deblur_use}" == "true" ] 
		then
			read -p "Deblur Trim length: " deblur_trimlength
		fi
		read -p "DADA2 trunclength F : " dada_trunclength_f
		read -p "DADA2 trunclength R : " dada_trunclength_r
		read -p "DADA2 maxEE F : " dada_maxEE_f
		read -p "DADA2 maxEE R : " dada_maxEE_r
		;;
		n|N )  : ;;
		* ) echo -e   "[     ${RED}ERROR${NC}   ] invalid" 1>&2; exit 1;;
	esac
}

run_16s() {
	####################################
	mkdir -p ${OUTPUT}
	OUTPUT=`realpath ${OUTPUT}`
	metaname="$(basename -- $meta)"
	cp ${meta} ${OUTPUT}/${metaname}
	cd ${OUTPUT}
	if [ "$trimmomatic_use" == "True" ]
	then
	echo -e "${YEL}#########    Trimmomatic Running    ########${NC}"
	conda activate main
	trimmer_pe () { # takes R1 file 
		# prepare file names 
		f1=$1
		f2=$(echo $f1 | sed 's/R1/R2/')
		f1_base="$(basename -- $f1)"
		f2_base=$(echo $f1_base | sed 's/R1/R2/')
		logname=${f1_base%'.fastq'}'.log'
		summaryname=${f1_base%'.fastq'}'.summry'
		
		new_f1=${f1_base%'.fastq'}'.trim.fastq'
		newf1=${new_f1}
		newf_2=${f2_base%'.fastq'}'.trim.fastq'
		newf2=${newf_2}

		newf_1U=${f1_base%'.fastq'}'.trim.se.fastq'
		newf1U=${newf_1U}
		newf_2U=${f2_base%'.fastq'}'.trim.se.fastq'
		newf2U=${newf_2U}
		
		
		adap="$CONDA_PREFIX/share/trimmomatic-0.39-1/adapters"
		
		# performes trimmomatic on sample with 20 threads, no adaptor trimming 
		# to cut adaptor modify adap to proper adaptor and uncomment ILLUMINACLIP line 
		trimmomatic PE -threads ${threads} -phred33 -trimlog ${logname} \
				-summary ${summaryname} $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
				SLIDINGWINDOW:4:10 MINLEN:30 \
				ILLUMINACLIP:$adap/TruSeq3-PE-2.fa:2:30:10:1 >> trim/LOG.txt 2>> trim/ERRORs.txt

		## PE -> paired ended
		## SLIDINGWINDOW: Performs a sliding window trimming approach.
		## ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads
	}

	mkdir -p trim
	mkdir -p trim/logs

	for i in ${INPUT}/*R1*
	do  
		i_base="$(basename -- $i)"
		echo -e "${BLU}${i_base} trimming started${NC}"
		trimmer_pe $i 
	done
	echo -e "${GRE}Finalizing...${NC}"
	rm *.se*
	mv *.log trim/logs >> trim/LOG.txt 2>> trim/ERRORs.txt
	mv *.summry trim/logs >> trim/LOG.txt 2>> trim/ERRORs.txt
	mv *.trim* trim >> trim/LOG.txt 2>> trim/ERRORs.txt
	cd trim
	mv *.txt logs
	for i in *.trim.fastq
	do
		mv $i ${i%'.trim.fastq'}'.fastq'
	done
	cd .. 
	INPUT=`realpath ./trim`
	echo -e "${RED}New path is ${INPUT}${NC}"
	conda deactivate
	fi
	echo -e "${YEL}#######    Generating sample table    ######${NC}"
	Rscript ~/GUAP/bin/16s/bash/generate_sample_table.R -i ${INPUT} -o ${OUTPUT}
	mkdir -p QIIME2
	mkdir -p QIIME2/phyloseq
	mkdir -p QIIME2/classify
	mkdir -p QIIME2/visualization
	mkdir -p QIIME2/align
	mkdir -p QIIME2/exports
	mkdir -p QIIME2/exports/abundance
	mkdir -p QIIME2/taxonomy

	if [ "${deblur_use}" == "TRUE" ] || [ "${deblur_use}" == "True" ] || [ "${deblur_use}" == "T" ] || [ "${deblur_use}" == "true" ] 
	then
		conda activate qiime2-2021.4
		cd QIIME2
		mkdir -p deblur
		mkdir -p demux
		cd ../
		echo -e "${YEL}###########      QIIME import    ###########${NC}"
		qiime tools import \
		--type 'SampleData[PairedEndSequencesWithQuality]' \
		--input-path samples.tsv \
		--output-path QIIME2/demux/paired-end-demux.qza \
		--input-format PairedEndFastqManifestPhred33V2
			error_cheker $?

		cd QIIME2

		# demultiplexing
		qiime demux summarize \
		--i-data demux/paired-end-demux.qza \
		--o-visualization visualization/demux.qzv 
			error_cheker $?

		echo -e "${YEL}###########    Quality filter    ###########${NC}"

		qiime quality-filter q-score \
		--i-demux demux/paired-end-demux.qza \
		--o-filtered-sequences demux/demux-filtered.qza \
		--o-filter-stats demux/demux-filter-stats.qza
			error_cheker $?


		echo -e "${YEL}################     ASV     ###############${NC}"
		qiime deblur denoise-16S \
		--i-demultiplexed-seqs demux/demux-filtered.qza \
		--p-trim-length ${deblur_trimlength} \
		--o-representative-sequences rep-seqs.qza \
		--o-table table.qza \
		--p-sample-stats \
		--o-stats deblur/deblur-stats.qza
			error_cheker $?

		qiime metadata tabulate \
		--m-input-file demux/demux-filter-stats.qza \
		--o-visualization visualization/demux-filter-stats.qzv 
			error_cheker $?

		qiime deblur visualize-stats \
		--i-deblur-stats deblur/deblur-stats.qza \
		--o-visualization visualization/deblur-stats.qzv
			error_cheker $?
	else
		############  R part ###########
		echo -e "${YEL}########    DADA2 Rscript running   ########${NC}"
		Rscript ~/GUAP/bin/16s/bash/G16s.v0.9.R \
			-i ${INPUT} -o "${PWD}/DADA2"\
			-p ${threads} \
			-n ${name_of_file} \
			-t ${dada_trunclength_f} \
			-T ${dada_trunclength_r} \
			-l ${dada_trim_r} \
			-r ${dada_trim_l} \
			-e ${dada_maxEE_f} \
			-E ${dada_maxEE_r} -s
		error_cheker $?

		###############################

		conda activate qiime2-2021.4
		cd QIIME2
		echo -e "${YEL}###########      QIIME import    ###########${NC}"
		qiime tools import \
			--input-path ../DADA2/rep-seqs.fna \
			--type "FeatureData[Sequence]" \
			--output-path rep-seqs.qza
    	error_cheker $?

		echo -n "#OTU Table" | cat - ../DADA2/seqtab-nochim.txt > biom-table.txt
    	error_cheker $?

		biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5
    	error_cheker $?

		qiime tools import \
		--input-path table.biom \
		--type "FeatureTable[Frequency]" \
		--input-format BIOMV210Format \
		--output-path table.qza
	fi
    error_cheker $?

	###############################
	if [ "${trainset_train}" == "TRUE" ] || [ "${trainset_train}" == "True" ] || [ "${trainset_train}" == "T" ] || [ "${trainset_train}" == "true" ] 
	then
		conda activate qiime2
		qiime feature-classifier extract-reads \
			--i-sequences  ${indexes_trainclassifier_seqs}\
			--p-f-primer ${trainset_fprimer} \
			--p-r-primer ${trainset_rprimer} \
			--p-min-length ${trainset_minlength} \
			--p-max-length ${trainset_maxlength} \
			--o-reads ref-seqs.qza
		error_cheker $?
		qiime feature-classifier fit-classifier-naive-bayes \
			--i-reference-reads ref-seqs.qza \
			--i-reference-taxonomy ${indexes_trainclassifier_taxa} \
			--o-classifier ${indexes_classifier}
		error_cheker $?
	fi

	echo -e "${YEL}#############      classify     ############${NC}"
	qiime feature-classifier classify-sklearn \
		--i-classifier ${indexes_classifier} \
		--i-reads rep-seqs.qza --p-n-jobs ${threads} \
		--o-classification classify/taxonomy.qza 

	echo -e "${YEL}###############     align     ##############${NC}"
	qiime phylogeny align-to-tree-mafft-fasttree \
		--i-sequences rep-seqs.qza \
		--o-alignment align/aligned-rep-seqs.qza \
		--o-masked-alignment align/masked-aligned-rep-seqs.qza \
		--o-tree align/unrooted-tree.qza \
		--o-rooted-tree align/rooted-tree.qza --p-n-threads ${threads}
		error_cheker $?
	echo -e "${YEL}#######    Exporting phyloseq object   #####${NC}"
	Rscript ~/GUAP/bin/16s/bash/export_phyloseq.R \
		-w "${OUTPUT}" \
		-s "${OUTPUT}/QIIME2/table.qza"\
		-m "${OUTPUT}/sample-metadata.tsv" \
		-t "${OUTPUT}/QIIME2/align/rooted-tree.qza" \
		-x "${OUTPUT}/QIIME2/classify/taxonomy.qza" \
		-o "${OUTPUT}/ps.rds"
	echo -e "${YEL}##########    Generating figures   #########${NC}"
	Rscript ~/GUAP/bin/16s/bash/QIIME2_import_viz.R \
		-w "${OUTPUT}" \
		-i "${OUTPUT}/ps.rds"\
		-r ${core_met_phylo_samplingdepth} \
		-n ${name_of_file} \
		-c condition \
		-S 21-days \
		-C Before_treatment
	echo -e "${YEL}#############     summarize     ############${NC}"
	qiime feature-table summarize \
	--i-table ./table.qza \
	--o-visualization visualization/table.qzv \
	--m-sample-metadata-file ../sample-metadata.tsv

	qiime feature-table tabulate-seqs \
	--i-data ./rep-seqs.qza \
	--o-visualization visualization/rep-seqs.qzv
	echo -e "${YEL}#############     taxa plot     ############${NC}"
	qiime taxa barplot \
		--i-table table.qza \
		--i-taxonomy classify/taxonomy.qza \
		--m-metadata-file ../sample-metadata.tsv \
		--o-visualization visualization/bar_plot.qzv 

	qiime tools export --input-path classify/taxonomy.qza --output-path "./classify" 
	echo -e "${YEL}###############    tabulate   ##############${NC}"
		qiime metadata tabulate \
		--m-input-file classify/taxonomy.qza \
		--o-visualization visualization/taxonomy.qzv 
	echo -e "${YEL}###########    export biom     #############${NC}"

	# Export OTU table
	qiime tools export \
		--input-path table.qza \
		--output-path phyloseq 
    error_cheker $?

	echo -e "${YEL}#############    biom to csv    ############${NC}"
	# Convert biom format to tsv format
	biom convert \
		-i ./phyloseq/feature-table.biom \
		-o ./phyloseq/otu_table.tsv \
		--to-tsv
    error_cheker $?


	cd phyloseq
	sed -i "1d" otu_table.tsv
	sed -i "s/#OTU ID//" otu_table.tsv
	cd ../
    error_cheker $?

	echo -e "${YEL}###### Export representative sequences #####${NC}"
	# Export representative sequences
	qiime tools export \
		--input-path rep-seqs.qza \
		--output-path phyloseq

    error_cheker $?

	echo -e "${YEL}##########     export tree files    ########${NC}"
	# Export tree files
	qiime tools export \
		--input-path align/unrooted-tree.qza \
		--output-path phyloseq

    error_cheker $?

	cd phyloseq
	mv tree.nwk unrooted_tree.nwk
	cd ../
    error_cheker $?

	echo -e "${YEL}############     rooted tree    ############${NC}"
	qiime tools export \
		--input-path align/rooted-tree.qza \
		--output-path phyloseq
    error_cheker $?

	cd phyloseq
	mv tree.nwk rooted_tree.nwk
	cd ../
	
	if [ ${interactive} == "TRUE" ]; then
			read -p $'\x1b[1;34mmodify rarefacrtion params (y/n)?\e[0m' CHOICE
		case "$CHOICE" in 
		y|Y ) 
			read -p "Alpha rarefaction max depth : " alphararefaction_maxdepth
			read -p "Rarefaction sampling depth  : " core_met_phylo_samplingdepth
		;;
		n|N )  : ;;
		* ) echo -e   "[     ${RED}ERROR${NC}   ] invalid" 1>&2; exit 1;;
		esac
	fi
    error_cheker $?

	echo -e "${YEL}########  core-metrics-phylogenetic  #######${NC}"
	qiime diversity core-metrics-phylogenetic \
	--i-phylogeny align/rooted-tree.qza \
	--i-table table.qza \
	--p-sampling-depth ${core_met_phylo_samplingdepth} \
	--m-metadata-file ../sample-metadata.tsv \
	--output-dir core-metrics-results

    error_cheker $?

	echo -e "${YEL}#######  diversity alpha-rarefaction  ######${NC}"
	qiime diversity alpha-rarefaction \
	--i-table table.qza \
	--i-phylogeny align/rooted-tree.qza \
	--m-metadata-file ../sample-metadata.tsv \
	--p-max-depth ${alphararefaction_maxdepth} \
	--o-visualization visualization/alpha-rarefaction.qzv

	qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/shannon_vector.qza \
	--m-metadata-file ../sample-metadata.tsv \
	--o-visualization visualization/shannon_vector-group-significance.qzv

	echo -e "${YEL}########  alpha-group-significance  ########${NC}"
	qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
	--m-metadata-file ../sample-metadata.tsv \
	--o-visualization visualization/faith-pd-group-significance.qzv 

	echo -e "${YEL}########  evenness-group-significanc  ######${NC}"
	qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/evenness_vector.qza \
	--m-metadata-file ../sample-metadata.tsv \
	--o-visualization visualization/evenness-group-significance.qzv 

	echo -e "${YEL}##### unweighted-unifrac-subject-group #####${NC}"
	qiime diversity beta-group-significance \
	--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	--m-metadata-file ../sample-metadata.tsv \
	--m-metadata-column condition \
	--o-visualization visualization/unweighted-unifrac-condition-significance.qzv \
	--p-pairwise 

	echo -e "${YEL}###########      emperor plot     ##########${NC}"
	qiime emperor plot \
	--i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
	--m-metadata-file ../sample-metadata.tsv \
	--p-custom-axes \
	--o-visualization visualization/unweighted-unifrac-emperor.qzv 

	echo -e "${YEL}##########    alpha and beta    ############${NC}"
	qiime emperor plot \
	--i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
	--m-metadata-file ../sample-metadata.tsv \
	--o-visualization visualization/bray-curtis-emperor.qzv 

	echo -e "${YEL}#############     exports    ###############${NC}"

	qiime tools export  \
	--input-path visualization/alpha-rarefaction.qzv \
	--output-path exports/alpha-rarefaction 

	qiime tools export  \
	--input-path visualization/evenness-group-significance.qzv \
	--output-path exports/evenness-group-significance 

	qiime tools export  \
	--input-path visualization/faith-pd-group-significance.qzv \
	--output-path exports/faith-pd-group-significance 

	qiime tools export  \
	--input-path visualization/shannon_vector-group-significance.qzv \
	--output-path exports/shannon_vector-group-significance 
	#export the both absoulte and relative abundance 
	mkdir -p exports
	mkdir -p exports/abundance
	###############################
	echo -e "${YEL}########      absoulutes STARTED    ########${NC}"

	absolutes (){
		i=$1
		taxaa=$2
		echo -e "${YEL}##### running absoulutes on $i $taxaa ######${NC}"

		#absolute
		qiime taxa collapse  \
		--i-table table.qza \
		--i-taxonomy classify/taxonomy.qza  \
		--p-level $i  \
		--o-collapsed-table exports/abundance/level-$i-$taxaa-table.qza 
		error_cheker $?
		qiime tools export \
		--input-path exports/abundance/level-$i-$taxaa-table.qza \
		--output-path taxonomy
		error_cheker $?
		biom convert -i taxonomy/feature-table.biom -o taxonomy/abs-level-$i-$taxaa-table.tsv \
		--to-tsv
		#relative
		qiime feature-table relative-frequency \
		--i-table exports/abundance/level-$i-$taxaa-table.qza \
		--o-relative-frequency-table exports/abundance/rel-level-$i-$taxaa-table.qza
		error_cheker $?
		qiime tools export \
		--input-path exports/abundance/rel-level-$i-$taxaa-table.qza \
		--output-path taxonomy
		error_cheker $?
		biom convert -i taxonomy/feature-table.biom -o taxonomy/rel-level-$i-$taxaa-table.tsv \
		--to-tsv
		error_cheker $?
	}

	for i in {1..7} 
	do 
		if [ $i == "1" ]; then taxaa="Kingdom"; fi
		if [ $i == "2" ]; then taxaa="Phylum"; fi
		if [ $i == "3" ]; then taxaa="Class"; fi
		if [ $i == "4" ]; then taxaa="Order"; fi
		if [ $i == "5" ]; then taxaa="Family"; fi
		if [ $i == "6" ]; then taxaa="Genus"; fi
		if [ $i == "7" ]; then taxaa="Species"; fi

		absolutes $i $taxaa &
		
	done
	echo -e "${GRE}###################  Done  #################${NC}"
}