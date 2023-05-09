get_n_samples(){
	max_n_samples=$(( $threads / $1 ))
	if [[ $max_n_samples -eq 0 ]] ; then max_n_samples=1 ; fi
}

gunzip_with_echo(){
	# performs gunzipping on samples to decompress the samples 
	i="${INPUT}/${i}"
	i_base="$(basename -- $i)"
    echo -e "${BLU}gunzipping $i_base...${NC}"
    gunzip $1 
    echo -e "${GRE}done $i_base... ${NC}"
}

parallel_gunzip(){
	# uses GNU parallel to gunzip all samples in parallel
	find $1 -maxdepth 1 -name '*.gz' -print0 | xargs -0 -I {} -P $2 bash -c 'gunzip_with_echo "$@"' _ {}
}

run_fastp(){

	fastp_pe () { 
		f1="${INPUT}/$1"
		f2=$(echo $f1 | sed "s/_${R}1/_${R}2/")
		f1_base="$(basename -- $f1)"
		f2_base=$(echo $f1_base | sed "s/_${R}1/_${R}2/")
		logname_=${f1_base%"_${R}1${tail}${ext}"}'.html'
		logname="QC/fastp/${logname_}"
		new_f1=${f1_base%"$ext"}'.filterred.fastq'
		newf1="QC/fastp/${new_f1}"
		newf_2=${f2_base%"$ext"}'.filterred.fastq'
		newf2="QC/fastp/${newf_2}"
		
        fastp --in1 ${f1} --in2 ${f2} \
            --out1 ${newf1} --out2 ${newf2} \
            --thread 4 -h ${logname} > logs/QC/${f1_base%"${R}1${tail}${ext}"}.log 2> logs/QC/${f1_base%"${R}1${tail}${ext}"}.error
	}

	export R=$R
	export ext=$ext
	export INPUT=$INPUT
	export RED=$RED
	export GRE=$GRE
	export BLU=$BLU
	export YEL=$YEL
	export NC=$NC
	export -f fastp_pe

	mkdir -p QC/fastp
	mkdir -p logs/QC

	echo -e "${YEL}#########    Fastp Running    ########${NC}"
	get_n_samples 4
	parallel --bar -j $max_n_samples fastp_pe ::: ${R1s[@]}
	error_cheker $? 0 "Running Fastp"
	echo -e "${GRE}Finalizing...${NC}"
	cd QC/fastp
	for i in *.filterred.fastq
	do
		mv $i ${i%'.filterred.fastq'}'.fastq'
	done
	cd ../../
	rm fastp.json
}

QC(){

	if [[ "$skip_QC" == "False"  || "$skip_QC" == "false" ]] ; then
		step_checker 0 run_fastp
		INPUT="${working_dir}/QC/fastp"
	fi

}

index_ref(){
	echo -e "${YEL}#########    Indexing Refs    ########${NC}"
	inputr=$1
	outputr=$2
	mkdir -p ${OUTPUT}/bowtie2_index
	mkdir -p logs/index_ref/
	bowtie2-build --threads ${threads_index} ${inputr} ${OUTPUT}/bowtie2_index/${outputr} > logs/index_ref/${outputr}.log 2> logs/index_ref/${outputr}.error
	samtools faidx ${inputr}

}

align_pe () { 
	f1=$1
	f2=$2
	index_to_use=$3
	f1_base=$4
	outdirs=$5
	out=$6
	mkdir -p ${outdirs}
	echo -e "${YEL}#######  Aliging ${f1_base}   #####${NC}"
	bowtie2 --threads ${threads_align} -x ${index_to_use} \
			-1 ${f1} \
			-2 ${f2} \
			-S "${outdirs}/${out}" > logs/mapping/${f1_base%"${R}1${tail}${ext}"}.log 2> logs/mapping/${f1_base%"_${R}1${tail}${ext}"}.error
}

fix_mate (){
	in=$1
	out=$2
	echo -e "${YEL}Fix-mate${NC}"
    # fix mates and compres 
    samtools sort --threads ${threads_align} -n \
                  --threads ${threads_align} -O sam $in | samtools fixmate \
                  --threads ${threads_align} -m -O bam - $out
}

convert_to_bam(){
	echo -e "${YEL}convert to sorted bam${NC}"
    # convert to bam file and sort
    samtools sort --threads ${threads_align} -O bam $1 \
                  -o $2
}

mrk_dub(){
    echo -e "${YEL}remove duplicates ${NC}"
    # remove duplicates 
    samtools markdup --threads ${threads_align} -r -S $1 $2
}

mpping_stats(){
    echo -e "${YEL}Mapping statistics ${NC}"
    # Mapping statistics 
    samtools flagstat $1 > ${2}-mapping-stats.txt
    qualimap bamqc -bam $1 -outdir ${2}-alignment-stats > ${2}-qualimap.log 2>&1
}

extract_unmapped(){
    echo -e "${YEL}unmapped reads ${NC}"
    # unmapped reads 
    samtools view --threads ${threads_align} -b -f 4 $1 > ${1%".bam"}".unmapped.bam"
    samtools fastq --threads ${threads_align} -1 ${2}_R1.fastq.qz \
                   -2 ${2}_R2.fastq.qz \
                   -s ${2}_U.fastq.qz ${1%".bam"}".unmapped.bam"
}

extract_mapped(){
	samtools view --threads ${threads_align} -h -b -f 3 $1 > ${1%".bam"}".mapped.bam"
	samtools fastq --threads ${threads_align} -1 ${2}_R1.fastq.qz \
		-2 ${2}_R2.fastq.qz \
		-s ${2}_U.fastq.qz ${1%".bam"}".mapped.bam"
}

call_vcf(){
	echo -e "${YEL}call vcf ${NC}"
	bcftools mpileup --threads ${threads_align}  -Ou \
		-f $1 $2 | bcftools call --ploidy 1 \
		--threads ${threads_align} -mv -o $3
}

filter_vcf(){
	echo -e "${YEL}filter vcf ${NC}"
	bcftools view $1 | vcfutils.pl varFilter - > $2
}

assembly(){
	f1=$1
	f2=$(echo $f1 | sed "s/_${R}1/_${R}2/")
	outdir=$2
	mkdir -p $outdir
	unicycler -1 $f1 -2 $f2 \
		-o ${outdir} \
		--threads ${threads_assemble} --no_pilon --no_correct
}

align_workflow () {
	f1="$1"
	index_to_use="$2"
	outdirs="$3"
	index_fasta="$4"
	mkdir -p ${outdirs}/mapped
	mkdir -p ${outdirs}/unmapped
	mkdir -p logs/mapping
	f2=$(echo $f1 | sed "s/_${R}1/_${R}2/")
	f1_base="$(basename -- $f1)"
	file_base="${f1_base%"_${R}1${tail}${ext}"}"
	outsam="${file_base}.sam"
	step_checker 0 "align_pe ${f1} ${f2} ${index_to_use} ${f1_base} ${outdirs} ${outsam}"
	step_checker 0 "fix_mate ${outdirs}/${outsam} ${outdirs}/${file_base}.fixmate.sam"
	step_checker 0 "convert_to_bam  ${outdirs}/${file_base}.fixmate.sam  ${outdirs}/${file_base}.bam"
	step_checker 0 "mrk_dub ${outdirs}/${file_base}.bam ${outdirs}/${file_base}.dedub.bam"
	step_checker 0 "extract_mapped ${outdirs}/${file_base}.dedub.bam ${outdirs}/mapped/${file_base}"
	step_checker 0 "extract_unmapped ${outdirs}/${file_base}.dedub.bam ${outdirs}/unmapped/${file_base}"
	step_checker 1 "mpping_stats ${outdirs}/${file_base}.dedub.bam ${outdirs}/${file_base}"
	step_checker 0 "call_vcf ${index_fasta} ${outdirs}/${file_base}.dedub.bam  ${outdirs}/${file_base}.vcf"
	step_checker 0 "filter_vcf ${outdirs}/${file_base}.vcf ${outdirs}/${file_base}filttered.vcf"
	step_checker 1 "assembly ${outdirs}/unmapped/${file_base}_R1.fastq.qz ${outdirs}/Assembly/unmapped_reads"
	step_checker 1 "assembly ${outdirs}/mapped/${file_base}_R1.fastq.qz ${outdirs}/Assembly/mapped_reads"
	if [[ ${skip_check_contaminant} != true ]]; then
		step_checker 1 "de_contaminate ${kraken2_index} ${outdirs}/unmapped/${file_base}_R1.fastq.qz ${outdirs}/contaminant/unmapped_reads"
	fi
}


de_contaminate () {
	echo -e "${YEL}Kraken2 is running...${NC}"
	db="$1"
	f1="$2"
	outdirs="$3"
	f1_base="$(basename -- $f1)"
	file_base="${f1_base%"_${R}1${tail}${ext}"}"
	mkdir -p ${outdirs}
	f2=$(echo $f1 | sed "s/_${R}1/_${R}2/")

	kraken2 --db ${db} --threads ${threads_kraken} \
		--report "${outdirs}/${file_base}.report" \
		--paired $f1 $f2 > "${outdirs}/${file_base}.out"
}

prep_dirs(){
	IDs=$( tail -n +2 ${OUTPUT}/samples.tsv | while read -r uniquename id restofname extension R1 R2 tail PE path R; do echo "$id"; done )
	for i in ${IDs[@]} ; do
		mkdir -p ${OUTPUT}/${i}/Assembly/mapped_reads
		mkdir -p ${OUTPUT}/${i}/Assembly/raw_reads
		mkdir -p ${OUTPUT}/${i}/mapping
		mkdir -p ${OUTPUT}/${i}/contamination_check
	done
}

run_WGS(){

	####################################

	cd ${OUTPUT}
	prep_dirs
	if [ $compressed == "True" ]; then 
		parallel_gunzip ${INPUT} ${threads} >> log.txt
	fi
	
	ext=${ext%'.gz'}''
	R1s=$( tail -n +2 ${OUTPUT}/samples.tsv | while read -r uniquename id restofname extension R1 R2 tail PE path R; do echo "${R1%'.gz'}"; done ) 
	R2s=$( tail -n +2 ${OUTPUT}/samples.tsv | while read -r uniquename id restofname extension R1 R2 tail PE path R; do echo "${R2%'.gz'}"; done ) 
	
	step_checker 0 QC
	
	
	if [[ $index_path == "null" ]]; then
		if [[ org_name != null || org_id != null || org_txid != null || org_auto_download != false ]]; then 
			step_checker 0 "index_ref ${ref} ${ref_name}"
			index_path=`realpath bowtie2_index/${ref_name}`
		fi
	fi
	
	for i in ${R1s[@]}; do
		ID=${i%"${R1_pattern}"}''
		if [[ ${pre_check_contaminant} == true ]]; then
			step_checker 1 "de_contaminate ${kraken2_index} ${INPUT}/${i} ${ID}/contaminant/Raw_reads"
		fi
		step_checker 0 "align_workflow ${INPUT}/${i} ${index_path} ${ID}/mapping/raw_reads ${ref}"
		if [[ ${skip_assembly} != true || ${skip_assembly} != True ]];then 
			step_checker 0 "assembly ${INPUT}/${i} ${ID}/Assembly/raw_reads"
		fi

	done

}

