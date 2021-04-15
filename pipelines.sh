# file containg piplines 
# opens functions script file
source functions.sh

# get working directory
current_dir=`pwd`

# metagenomics pipline 
kraken_pipeline_prog () { # takes 4 arguments: 1. aligner name 2. samples dir 3. number of threads 4. output dir 
    sample_check $2 True # checks sampels file names
    pipline_value=$? # exit value of function checks to choose single end or paired end
    aligner_name=$1
    files_dir=$2
    threads=$3
    outputdit=$4
    currentdir=`pwd`
    # creates directories 
    mkdir -p $outputdit/bams
    mkdir -p $outputdit/trimmed
    mkdir -p $outputdit/logs
    mkdir -p $outputdit/final_samples
    mkdir -p $outputdit/kraken
    mkdir -p $outputdit/Host_removed_samples
    mkdir -p $outputdit/korona
    touch $outputdit/LOG.txt
    # contiue with paired end analysis 
    PE_analysis(){
        echo -e "${YEL}PERFORMING KRAKEN USING ${aligner_name}${NC}"
        cd $files_dir
        nu_of_R1=`ls *R1* | wc -l` # gets number of R1 files
        total_processes=$(( $nu_of_R1 * 100 )) # assigns total progress value according to number of samples 
        current_value=$((1)) 
        echo -e "\n${YEL}Analysis Started: ${NC}"
        progress_bar $total_processes $(( $nu_of_R1 * 1 )) "${BRIGHT}${BLUE}TRIMMING SAMPLES${NORMAL}"
        current_value=$((0))

        for target in *R1* # starts with trimming of samples 
        do
            if [ "${i: -9}" == '.fastq.gz' ];then # uses trimmomatic function in functions file and move outputs to proper folders 
                trimmer_pe $target >> $outputdit/LOG.txt 2>&1 # redirects standard output and standard errors to log file
                mv *.log $outputdit/logs >> $outputdit/LOG.txt 2>&1
                mv *.summry $outputdit/logs >> $outputdit/LOG.txt 2>&1
                mv *.trim* $outputdit/trimmed >> $outputdit/LOG.txt 2>&1
                current_value=$(( current_value + 15 )) # updates current value 
                echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
                progress_bar $total_processes $current_value "${BRIGHT}${BLUE}TRIMMING SAMPLES${NORMAL}" # prints progress
            fi
        done
        
        cd $outputdit/trimmed
        mkdir -p SEs
        mv *.se.* SEs
        #  aligning samples to aligner acoording to user deffinition
        for i in *R1*
        do
            if [ "${i: -14}" == '.trim.fastq.gz' ];then # gets proper file 
                echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1 # reports file name to log file 
                progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}" # updates progres bar 
                # checks aligner name and use corresponding function
                if [ $aligner_name == "hisat2" ];then 
                    runhisat2 $threads $i >> $outputdit/LOG.txt 2>&1
                elif [ $aligner_name == "bwa" ]; then
                    bwa_rm_human $threads $i >> $outputdit/LOG.txt 2>&1
                elif [ $aligner_name == "bowtie2" ]; then
                    runbowtie2 $threads $i >> $outputdit/LOG.txt 2>&1
                else
                    echo -e   "[     ${RED}ERROR${NC}   ] invalid" 1>&2; exit 1
                fi
                mv ${i%'.fastq.gz'}'.sorted.bam' $outputdit/bams >> $outputdit/LOG.txt 2>&1 # moves files to proper directory
                current_value=$(( current_value + 25 )) # updates progress bar value
                progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}" # prints progress 
                echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
            fi
        done

        cd $outputdit/bams

        # extracts unmapped from bam files ( not human reads )
        for i in *.bam
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} REMOVING  HOST  ${NORMAL}"
            extract_unmapped $i ${i%'.trim.sorted.bam'}'_Host_Removed' ${i%'R1_001.trim.sorted.bam'} $threads >> $outputdit/LOG.txt 2>&1 
            rm ${i%'.trim.sorted.bam'}'_Host_Removed' >> $outputdit/LOG.txt 2>&1
            mv *.fastq.gz $outputdit/Host_removed_samples >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done

        cd $outputdit/Host_removed_samples
        
        # perfomrs kraken 
        for i in *R1.fastq.gz
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}PERFORMING KRAKEN${NORMAL}"
            kraken_ $threads ${i%'.fastq.gz'}'' $i >> $outputdit/LOG.txt 2>&1
            mv *.report *.out $outputdit/kraken >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 50 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done

        cd $outputdit/kraken

        # generates korona
        for i in *.out
        do
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}GENRATING KORONA${NORMAL}"
            kraken2korona ${i%'.out'}'.html' $i >> $outputdit/LOG.txt 2>&1
            mv ${i%'.out'}'.html' $outputdit/korona >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done
        progress_bar $total_processes $current_value "${BRIGHT}${GREEN}KRAKEN ANALYSIS DONE
        ${NORMAL}"
    }

    # performes single ended analysis
    SE_anaysis(){
        echo -e "${YEL}PERFORMING KRAKEN USING ${aligner_name}${NC}"
        cd $files_dir
        nu_of_R1=`ls * | wc -l`
        total_processes=$(( $nu_of_R1 * 100 ))
        current_value=$((1))
        echo -e "\n${YEL}Analysis Started: ${NC}"
        progress_bar $total_processes $(( $nu_of_R1 * 1 )) "${BRIGHT}${BLUE}TRIMMING SAMPLES${NORMAL}"
        current_value=$((0))
        # trimming samples 
        for target in *.fastq.gz
        do 
            trimmer $target >> $outputdit/LOG.txt 2>&1
            mv *.log $outputdit/logs >> $outputdit/LOG.txt 2>&1
            mv *.summry $outputdit/logs >> $outputdit/LOG.txt 2>&1
            mv *.trim* $outputdit/trimmed >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 15 ))
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}TRIMMING SAMPLES${NORMAL}" 
        done
        
        cd $outputdit/trimmed
        #  aligning samples to aligner acoording to user deffinition
        for i in *.trim.fastq.gz
        do
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}"
            # continue with aligner according to aligner name
            if [ $aligner_name == "hisat2" ]; then 
                run_SE_hisat2 $threads $i >> $outputdit/LOG.txt 2>&1
            elif [ $aligner_name == "bwa" ]; then 
                bwa_SE_rm_human $threads $i >> $outputdit/LOG.txt 2>&1
            elif [ $aligner_name == "bowtie2" ]; then 
                run_SE_bowtie2 $threads $i >> $outputdit/LOG.txt 2>&1
            else
                echo -e   "[     ${RED}ERROR${NC}   ] invalid" 1>&2; exit 1
            fi
            mv ${i%'.fastq.gz'}'.sorted.bam' $outputdit/bams >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 25 )) 
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}"
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done
        
        cd $outputdit/bams
        # extract unmapped reads
        for i in *.bam
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} REMOVING  HOST  ${NORMAL}"
            extract_SE_unmapped $i ${i%'.trim.sorted.bam'}'_Host_Removed.fastq' >> $outputdit/LOG.txt 2>&1
            mv *.fastq $outputdit/Host_removed_samples >> $outputdit/LOG.txt 2>&1           
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
            
        done

        cd $outputdit/Host_removed_samples
        # performes kraken 
        for i in *
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}PERFORMING KRAKEN${NORMAL}"
            kraken_SE $threads ${i%'.fastq'}'' $i >> $outputdit/LOG.txt 2>&1
            mv *.report *.out $outputdit/kraken >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 50 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done

        cd $outputdit/kraken
        # gnerates korona
        for i in *.out
        do
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}GENRATING KORONA${NORMAL}"
            kraken2korona ${i%'.out'}'.html' $i >> $outputdit/LOG.txt 2>&1
            mv ${i%'.out'}'.html' $outputdit/korona >> $outputdit/LOG.txt 2>&1            
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done
        progress_bar $total_processes $current_value "${BRIGHT}${GREEN}KRAKEN ANALYSIS DONE
        ${NORMAL}"
    }
    # checks sample_check function output to choose pipeline 
    if [ $pipline_value == 1 ]
    then
        echo -e "${YEL}SE analysis will continue${NC}"
        continue_
        SE_anaysis
    elif [ $pipline_value == 2 ]
    then 
        echo -e "${YEL}PE analysis will continue${NC}"
        continue_
        PE_analysis
    else
        echo -e "${RED}ERORR${NC}"
        exit 
    fi
}

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
    # Single ended analysis
    SE_anaysis(){
        # calculating number of samples for progress bar 
        echo -e "${YEL}PERFORMING miRNATRA 1 DIRECT${NC}"
        cd $files_dir
        nu_of_R1=`ls * | wc -l`
        total_processes=$(( $nu_of_R1 * 100 ))
        current_value=$((1))
        echo -e "\n${YEL}Analysis Started: ${NC}"
        progress_bar $total_processes $(( $nu_of_R1 * 1 )) "${BRIGHT}${BLUE}STARTING ANALYSIS${NORMAL}"
        current_value=$((0))
        # starts analysis with cutting adaptors
        for target in *.fastq
        do             
            cutadapt_mirna $target >> $outputdit/LOG.txt 2>&1
            mv *.trim.* $outputdit/Adaptor_QC >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 15 ))
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ADAPTOR  TRIMMING${NORMAL}" 
        done
        cd $outputdit/Adaptor_QC
        # aligning samples to mirbase database 
        for i in *.trim.fastq
        do
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase  MAPPING${NORMAL}"
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            bowtie_mirbase_count  $threads $i >> $outputdit/LOG.txt 2>&1 # this function exports counts file and unmapped file
            # moving each file to proper location
            mv ${i%'.fastq'}'.miRBase.sam' $outputdit/sams >> $outputdit/LOG.txt 2>&1 
            mv ${i%'.fastq'}'.miRBase.sorted.bam' $outputdit/Mapping_miRBase >> $outputdit/LOG.txt 2>&1
            mv ${i%'.fastq'}'.miRBase.sorted.bam.bai' $outputdit/Mapping_miRBase >> $outputdit/LOG.txt 2>&1
            mv ${i%'.fastq'}'_unmapped.fastq' $outputdit/mirbas_unmapped >> $outputdit/LOG.txt 2>&1
            mv ${i%'.fastq'}'.counts.txt' $outputdit/counts >> $outputdit/LOG.txt 2>&1

            current_value=$(( current_value + 15 )) 
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase MAPPING ${NORMAL}"
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done
        
        cd $outputdit/mirbas_unmapped
        # aligning unaligned reads to human genome 
        for i in *.fastq
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} human MAPPING  ${NORMAL}"
            bowtie_mirhuman $threads $i >> $outputdit/LOG.txt 2>&1 # this function exports a tagged bam file with miRNA locations to be counted 
            mv ${i%'.fastq'}'.hum.sam' $outputdit/sams >> $outputdit/LOG.txt 2>&1
            mv ${i%'.fastq'}'.hum.sorted.bam' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>&1
            mv ${i%'.fastq'}'.hum.sorted.bam.bai' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>&1
            mv ${i%'.fastq'}'_hum.tagged.bam' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>&1            
            current_value=$(( current_value + 20 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
            
        done
        
        cd $outputdit/Mapping_Genome
        # countinng tagged sites from mirna_list variable and tagged bam
        for i in *.tagged.bam
        do 
            sample_name=${i%'.trim_unmapped_hum.tagged.bam'}''
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} human COUNTING  ${NORMAL}"
            for j in ${mirna_list_array[@]}
            do
                j_count=`samtools view $i | grep $j | wc -l` 
                echo -e "$j\t$j_count" >> ${sample_name}'.count.txt'
            done
            mv ${sample_name}'.count.txt' $outputdit/counts >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 45 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done
        progress_bar $total_processes $current_value "${BRIGHT}${BLUE}  FINISHING    ${NORMAL}"

        cd $outputdit/counts
        # uses python script ( merge.py ) to merge counts from mirbase and mapped to ref-genome 
        python3 $current_dir/merge.py -i `pwd` -o all.csv >> $outputdit/LOG.txt 2>&1
        mv all.csv ../ >> $outputdit/LOG.txt 2>&1
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
        continue_
        SE_anaysis
    elif [ $pipline_value == 2 ]
    then 
        echo -e "${YEL}PE analysis will continue${NC}"
        continue_
        PE_analysis
    else
        echo -e "${RED}ERORR${NC}"
        exit 
    fi
}

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
        # starts analysis with cutting adaptors
        for target in *.fastq
        do             
            cutadapt_mirna $target >> $outputdit/LOG.txt 2>&1
            mv *.trim.* $outputdit/Adaptor_QC >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 35 ))
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ADAPTOR  TRIMMING${NORMAL}" 
        done
        cd $outputdit/Adaptor_QC
        # aligning samples to human-genome 
        for i in *.trim.fastq
        do
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase  MAPPING${NORMAL}"
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>&1
            bowtie_mirSecond  $threads $i >> $outputdit/LOG.txt 2>&1
            mv ${i%'.fastq'}'.hum.sam' $outputdit/Mapping_Genome >> $outputdit/LOG.txt 2>&1
            current_value=$(( current_value + 35 )) 
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} miRBase MAPPING ${NORMAL}"
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
        done
        
        cd $outputdit/Mapping_Genome
        # renaming sam files for featurecounts  
        for i in *.hum.sam
        do 
            mv $i ${i%'.trim.hum.sam'}''
        done
        all_sams=()
        progress_bar $total_processes $current_value "${BRIGHT}${BLUE} COUNTING            ${NORMAL}"
        # creates an array with all sample names 
        for i in *
        do 
            all_sams+=$i
            all_sams+=" "
        done
        # performs featureCounts usning mirbase annotation file 
        featureCounts -t miRNA -g Name -O -T 8 -s 1 -M -a $mirbase_annotation -o "all.counts" $all_sams >> $outputdit/LOG.txt 2>&1
        mv "all.counts" $outputdit/counts >> $outputdit/LOG.txt 2>&1
        mv "all.counts.summary" $outputdit/counts >> $outputdit/LOG.txt 2>&1
        echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>&1
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
        continue_
        SE_anaysis
    elif [ $pipline_value == 2 ]
    then 
        echo -e "${YEL}PE analysis will continue${NC}"
        continue_
        PE_analysis
    else
        echo -e "${RED}ERORR${NC}"
        exit 
    fi
}

