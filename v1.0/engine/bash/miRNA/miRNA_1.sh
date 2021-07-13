source functions.sh

miRNAtra1_pipeline () { # takes 4 arguments: 1. aligner name 2. samples dir 3. number of threads 4. output dir
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
    mkdir -p $outputdit/sams
    mkdir -p $outputdit/Mapping_miRBase
    mkdir -p $outputdit/counts
    mkdir -p $outputdit/Mapping_Genome
    mkdir -p $outputdit/mirbas_unmapped
    touch $outputdit/LOG.txt
    touch $outputdit/ERRORs.txt

    # Single ended analysis
    SE_anaysis(){
        # calculating number of samples for progress bar 
        echo -e "${YEL}PERFORMING miRNATRA 1 DIRECT${NC}"
        cd $files_dir
        nu_of_R1=`ls * | wc -l`
        no_of_gz=`ls *.gz 2> /dev/null | wc -l `
        total_processes=$(( $nu_of_R1 * 100 ))
        current_value=$((1))
        echo -e "\n${YEL}Analysis Started: ${NC}"
        progress_bar $total_processes $(( $nu_of_R1 * 1 )) "${BRIGHT}${BLUE}STARTING ANALYSIS${NORMAL}"
        current_value=$((0))
        # starts analysis with cutting adaptors
        if [ $((no_of_gz)) -ne 0 ] # checks if samples are compressed 
        then 
            for i in *.fastq.gz;do
                gunzip $i
            done
        fi
        for target in *.fastq
        do             
            cutadapt_mirna $target >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv *.trim.* $outputdit/Adaptor_QC >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            current_value=$(( current_value + 15 ))
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ADAPTOR  TRIMMING${NORMAL}" 
        done
        cd $outputdit/Adaptor_QC
        # aligning samples to mirbase database 
        for i in *.trim.fastq
        do
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase  MAPPING${NORMAL}"
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            bowtie_mirbase_count  $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt # this function exports counts file and unmapped file
            error_cheker $?
            # moving each file to proper location
            mv ${i%'.fastq'}'.miRBase.sam' $outputdit/sams >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt 
            mv ${i%'.fastq'}'.miRBase.sorted.bam' $outputdit/Mapping_miRBase >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv ${i%'.fastq'}'.miRBase.sorted.bam.bai' $outputdit/Mapping_miRBase >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv ${i%'.fastq'}'_unmapped.fastq' $outputdit/mirbas_unmapped >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv ${i%'.fastq'}'.counts.txt' $outputdit/counts >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt

            current_value=$(( current_value + 15 )) 
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase MAPPING ${NORMAL}"
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done
        
        cd $outputdit/mirbas_unmapped
        # aligning unaligned reads to human genome 
        for i in *.fastq
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} human MAPPING  ${NORMAL}"
            bowtie_mirhuman $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt # this function exports a tagged bam file with miRNA locations to be counted 
            error_cheker $?
            mv ${i%'.fastq'}'.hum.sam' $outputdit/sams >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv ${i%'.fastq'}'.hum.sorted.bam' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv ${i%'.fastq'}'.hum.sorted.bam.bai' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv ${i%'.fastq'}'_hum.tagged.bam' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt            
            current_value=$(( current_value + 20 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            
        done
        
        cd $outputdit/Mapping_Genome
        # countinng tagged sites from mirna_list variable and tagged bam
        for i in *.tagged.bam
        do 
            sample_name=${i%'.trim_unmapped_hum.tagged.bam'}''
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} human COUNTING  ${NORMAL}"
            for j in ${mirna_list_array[@]}
            do
                j_count=`samtools view $i | grep $j | wc -l` 
                error_cheker $?
                echo -e "$j\t$j_count" >> ${sample_name}'.count.txt'
            done
            mv ${sample_name}'.count.txt' $outputdit/counts >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            current_value=$(( current_value + 45 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done
        progress_bar $total_processes $current_value "${BRIGHT}${BLUE}  FINISHING    ${NORMAL}"

        cd $outputdit/counts
        # uses python script ( merge.py ) to merge counts from mirbase and mapped to ref-genome 
        python3 $current_dir/merge.py -i `pwd` -o all.csv >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        mv all.csv ../ >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        to_add=$(( $nu_of_R1 * 5))
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
