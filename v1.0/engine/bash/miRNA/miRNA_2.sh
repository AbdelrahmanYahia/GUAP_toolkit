source functions.sh

miRNAtra2_pipeline () { # takes 4 arguments: 1. aligner name 2. samples dir 3. number of threads 4. output dir
    # set variables
    sample_check $2
    pipline_value=$?
    aligner_name=$1
    files_dir=$2
    threads=$3
    outputdit=$4
    currentdir=`pwd`
    # creates working directories 
    mkdir -p $outputdit/Adaptor_QC
    mkdir -p $outputdit/counts
    mkdir -p $outputdit/Mapping_Genome
    touch $outputdit/LOG.txt
    touch $outputdit/ERRORs.txt

    # Single ended analysis
    SE_anaysis(){
        # calculating number of samples for progress bar 
        echo -e "${YEL}PERFORMING miRNATRA 2 DIRECT${NC}"
        cd $files_dir
        nu_of_R1=`ls * | wc -l`
        total_processes=$(( $nu_of_R1 * 100 ))
        current_value=$((1))
        echo -e "\n${YEL}Analysis Started: ${NC}"
        progress_bar $total_processes $(( $nu_of_R1 * 1 )) "${BRIGHT}${BLUE}STARTING ANALYSIS${NORMAL}"
        current_value=$((0))
        no_of_gz=`ls *.gz 2> /dev/null | wc -l `
        if [ $((no_of_gz)) -ne 0 ] # checks if samples are compressed 
        then 
            for i in *.fastq.gz;do
                gunzip $i
            done
        fi
        # starts analysis with cutting adaptors
        for target in *.fastq
        do             
            cutadapt_mirna $target >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv *.trim.* $outputdit/Adaptor_QC >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            current_value=$(( current_value + 35 ))
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ADAPTOR  TRIMMING${NORMAL}" 2> /dev/null
        done
        cd $outputdit/Adaptor_QC
        # aligning samples to human-genome 
        for i in *.trim.fastq
        do
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase  MAPPING${NORMAL}" 2> /dev/null
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            bowtie_mirSecond  $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv ${i%'.fastq'}'.hum.sam' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            current_value=$(( current_value + 35 )) 
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase MAPPING ${NORMAL}" 2> /dev/null
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done
        
        cd $outputdit/Mapping_Genome
        # renaming sam files for featurecounts  
        for i in *.hum.sam
        do 
            mv $i ${i%'.trim.hum.sam'}''
        done
        all_sams=()
        progress_bar $total_processes $current_value "${BRIGHT}${BLUE} COUNTING            ${NORMAL}" 2> /dev/null
        # creates an array with all sample names 
        for i in *
        do 
            all_sams+=$i
            all_sams+=" "
        done
        # performs featureCounts usning mirbase annotation file 
        featureCounts -t miRNA -g Name -O -T 8 -s 1 -M -a $mirbase_annotation -o "all.counts" $all_sams >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        error_cheker $?
        mv "all.counts" $outputdit/counts >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        mv "all.counts.summary" $outputdit/counts >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        cd $outputdit/counts
        # modifies output for final counts table 
        tail -n +2 all.counts | cut -f 1,7- > all.txt
        mv all.txt ../
        to_add=$(( $nu_of_R1 * 30))
        current_value=$(( current_value + to_add ))
        progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALL miRNA DONE${NORMAL}" 

    }
    # paired end analysis under development 
    PE_analysis(){
        echo -e "${RED}PAIRED END ANALYSIS IS UNDER DEVELOPMENT${NC}"
        exit
    }
    if [ $pipline_value == 1 ]
    then
        echo -e "${YEL}SE analysis will continue${NC}"
        if [ ${y} == "TRUE" ]; then :; else continue_ ; fi

        SE_anaysis
    elif [ $pipline_value == 2 ]
    then 
        echo -e "${YEL}PE analysis will continue${NC}"
        if [ ${y} == "TRUE" ]; then :; else continue_ ; fi

        PE_analysis
    else
        echo -e "${RED}ERORR${NC}"
        exit 
    fi
}
