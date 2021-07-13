# file containg piplines 
# opens functions script file
source functions.sh

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
    touch $outputdit/ERRORs.txt

    # contiue with paired end analysis 
    PE_analysis(){
        echo -e "${YEL}PERFORMING KRAKEN USING ${aligner_name}${NC}"
        cd $files_dir
        nu_of_R1=`ls *R1* | wc -l` # gets number of R1 files
        if [ $(( $nu_of_R1 )) -eq 0 ];then
            exit
        fi
        total_processes=$(( $nu_of_R1 * 100 )) # assigns total progress value according to number of samples 
        current_value=$((1)) 
        echo -e "\n${YEL}Analysis Started: ${NC}"
        progress_bar $total_processes $(( $nu_of_R1 * 1 )) "${BRIGHT}${BLUE}TRIMMING SAMPLES${NORMAL}"
        current_value=$((0))

        for target in *R1* # starts with trimming of samples 
        do
            if [ "${i: -9}" == '.fastq.gz' ];then # uses trimmomatic function in functions file and move outputs to proper folders 
                trimmer_pe $target >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt # redirects standard output and standard errors to log file
                error_cheker $?
                mv *.log $outputdit/logs >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                error_cheker $?
                mv *.summry $outputdit/logs >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                error_cheker $?
                mv *.trim* $outputdit/trimmed >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                error_cheker $?
                current_value=$(( current_value + 15 )) # updates current value 
                echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
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
                echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt # reports file name to log file 
                progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}" # updates progres bar 
                # checks aligner name and use corresponding function
                if [ $aligner_name == "hisat2" ];then 
                    runhisat2 $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                    error_cheker $?
                elif [ $aligner_name == "bwa" ]; then
                    bwa_rm_human $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                    error_cheker $?
                elif [ $aligner_name == "bowtie2" ]; then
                    runbowtie2 $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                    error_cheker $?
                else
                    echo -e   "[     ${RED}ERROR${NC}   ] invalid" 1>&2; exit 1
                fi
                mv ${i%'.fastq.gz'}'.sorted.bam' $outputdit/bams >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt # moves files to proper directory
                error_cheker $?
                current_value=$(( current_value + 25 )) # updates progress bar value
                progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}" # prints progress 
                echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            fi
        done

        cd $outputdit/bams
        # extracts unmapped from bam files ( not human reads )
        for i in *.bam
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} REMOVING  HOST  ${NORMAL}"
            extract_unmapped $i ${i%'.trim.sorted.bam'}'_Host_Removed' ${i%'R1_001.trim.sorted.bam'} $threads >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt 
            error_cheker $?
            rm ${i%'.trim.sorted.bam'}'_Host_Removed' >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv *.fastq.gz $outputdit/Host_removed_samples >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done

        cd $outputdit/Host_removed_samples
        # perfomrs kraken 
        for i in *R1.fastq.gz
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}PERFORMING KRAKEN${NORMAL}"
            kraken_ $threads ${i%'.fastq.gz'}'' $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv *.report *.out $outputdit/kraken >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            current_value=$(( current_value + 50 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done

        cd $outputdit/kraken
        # generates korona
        for i in *.out
        do
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}GENRATING KORONA${NORMAL}"
            kraken2korona ${i%'.out'}'.html' $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv ${i%'.out'}'.html' $outputdit/korona >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
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
            trimmer $target >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv *.log $outputdit/logs >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv *.summry $outputdit/logs >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            mv *.trim* $outputdit/trimmed >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            current_value=$(( current_value + 15 ))
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}TRIMMING SAMPLES${NORMAL}" 
        done
        
        cd $outputdit/trimmed
        #  aligning samples to aligner acoording to user deffinition
        for i in *.trim.fastq.gz
        do
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}"
            # continue with aligner according to aligner name
            if [ $aligner_name == "hisat2" ]; then 
                run_SE_hisat2 $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                error_cheker $?
            elif [ $aligner_name == "bwa" ]; then 
                bwa_SE_rm_human $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                error_cheker $?
            elif [ $aligner_name == "bowtie2" ]; then 
                run_SE_bowtie2 $threads $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
                error_cheker $?
            else
                echo -e   "[     ${RED}ERROR${NC}   ] invalid" 1>&2; exit 1
            fi
            mv ${i%'.fastq.gz'}'.sorted.bam' $outputdit/bams >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            current_value=$(( current_value + 25 )) 
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}ALIGNING SAMPLES${NORMAL}"
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done
        
        cd $outputdit/bams
        # extract unmapped reads
        for i in *.bam
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE} REMOVING  HOST  ${NORMAL}"
            extract_SE_unmapped $i ${i%'.trim.sorted.bam'}'_Host_Removed.fastq' >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv *.fastq $outputdit/Host_removed_samples >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt           
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            
        done

        cd $outputdit/Host_removed_samples
        # performes kraken 
        for i in *
        do 
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}PERFORMING KRAKEN${NORMAL}"
            kraken_SE $threads ${i%'.fastq'}'' $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv *.report *.out $outputdit/kraken >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            current_value=$(( current_value + 50 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done

        cd $outputdit/kraken
        # gnerates korona
        for i in *.out
        do
            echo -e "SAMPLE IS : ${i}" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            progress_bar $total_processes $current_value "${BRIGHT}${BLUE}GENRATING KORONA${NORMAL}"
            kraken2korona ${i%'.out'}'.html' $i >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
            error_cheker $?
            mv ${i%'.out'}'.html' $outputdit/korona >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt            
            current_value=$(( current_value + 5 )) 
            echo -e "_________________________________________ \n" >> $outputdit/LOG.txt 2>> $outputdit/ERRORs.txt
        done
        progress_bar $total_processes $current_value "${BRIGHT}${GREEN}KRAKEN ANALYSIS DONE
        ${NORMAL}"
    }
    # checks sample_check function output to choose pipeline 
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
