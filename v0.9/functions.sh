#!/bin/bash
source indexes.sh

# progress bar ## takes three arguments 
# 1. total 2. current 3. message to print 
progress_bar() {
    BAR="#"
    MIN=1
    MAX=100
    TIMES=$1
    COUNT=$2
    MESSAGETOSHOW=$3
    PROGRESS=$(echo | awk -v t=$TIMES -v c=$COUNT -v max=$MAX '{ print int(c / t * max) }')
    if [ $(( PROGRESS )) -eq 0 ];then
        PROGRESS=1
    fi
    PROGRESS_BAR=""
    for _i in $(seq $MIN $PROGRESS); do PROGRESS_BAR="${PROGRESS_BAR}${BAR}"; done
    
    printf "\r[%-100s] %3d%% - %s " $PROGRESS_BAR $PROGRESS "$MESSAGETOSHOW"
}

# asks user to continue or not
continue_ (){

    read -p "Continue (y/n)?" CHOICE
    case "$CHOICE" in 
        y|Y ) 
            $@
            :
        ;;
        n|N ) echo -e "[      ${GRE}OK${NC}     ] Process stopped." 1>&2; exit 1;;
        * ) echo -e   "[     ${RED}ERROR${NC}   ] invalid" 1>&2; exit 1;;
    esac
}

# checks sample extension and paired or not 
sample_check(){
    dogz=$2 # if true will gunzip unzipped samlpes 
    counter=1
    return_value=0 # exit code 
    no_of_files=`ls ${1}/ | wc -l ` # gets number of files 
    echo -e "${YEL}NUMBER OF FILES = ${no_of_files}${NC}" 
    if [ $((no_of_files%2)) -ne 0 ] # checks if samples are not pairs 
    then 
        echo -e "${RED}SAMPLES NOT IN PAIRS${NC} ${YEL}CHOOSING SE PIPELINE${NC}"
        if [ `ls ${1}/*R1* 2> /dev/null | wc -l ` -ne 0 ];then
            for R1 in ${1}/*R1*
            do
                R2=$(echo $R1 | sed 's/R1/R2/')
                    if [ -f $R2 ]; then # if finds R2 returns erorr
                        echo -e "${RED}ERORR WITH SAMPLE NAMES (EC1)${NC}"
                        return_value=0
                        break
                    else
                        return_value=1
                    fi
            done
        else # if samples are paires 
            if [ `ls ${1}/*R2* 2> /dev/null | wc -l ` -ne 0 ];then 
                return_value=0
                echo -e "${RED}ERORR WITH SAMPLES NAMES (EC2)${NC}"
            else
                return_value=1
            fi
        fi
    else
        if [ `ls ${1}/*R1* 2> /dev/null | wc -l ` -ne 0 ] # checks presence of R2 
        then
            for R1 in ${1}/*R1*
            do
                R2=$(echo $R1 | sed 's/R1/R2/')
                if [ ! -f $R2 ]; then
                    return_value=1
                    break
                else
                    return_value=2
                fi
            done
        else
            echo -e "${RED}(EC3)${YEL} SAMPLES SEEMS TO BE SE ${NC}"
            return_value=1
        fi
    fi
    # checks files extensions 
    for i in ${1}/*
    do
        if [ "${i: -9}" == '.fastq.gz' ] # files should be with extension fastq.qz 
        then 
            echo -e "${i} ${GRE}PASSED${NC}"
            :
        elif [ ${i: -6} == ".fq.gz" ]
        then
            newname=`echo $i | sed 's/.fq.gz/.fastq.gz/'` # replaces fq.qz to fastq.qz
            mv $i $newname
            echo -e "${i} ${YEL}MODEFIED ${RED}(EC4)${NC}"
            :
        elif [ ${i: -6} == ".fastq" ] # if not compressed and analysis need compression it will compress file 
        then
            if [[ $dogz == 'True' ]]
            then
                gzip $i
                echo -e "${i} ${YEL}MODEFIED  ${RED}(EC5)${NC}"
                :
            else
                echo -e "${i} ${GRE}PASSED${NC}"
            fi
        elif [ ${i: -3} == ".fq" ]
        then
            newname=`echo $i | sed 's/.fq/.fastq/'`
            mv ${i} $newname
            if [ $dogz == 'True' ]
            then
                gzip ${i%'.fq'}'.fastq'
                echo -e "${i} ${YEL}MODEFIED  ${RED}(EC6)${NC}"
            else
                echo -e "${i} ${GRE}PASSED${NC}"
            fi
            :
        else
            echo -e "${RED}${i} DID NOT PASS (EC7)${NC}"
        fi
        counter=$(( counter + 1 ))
    done
    # checks for files with extension not fastq or fastq.gz 
    if [ `ls ${1}/ -I *'.fast'* | wc -w ` -ne 0 ];then 
        echo -e "${RED}FOUND FILES WITH EXTENSION NOT fastq.gz (EC 7)${NC}"
        return_value=0
    fi
    return $return_value
    
}

# checks sample extension and paired or not 
dir_check(){
    check_extension=$1 # if true will gunzip unzipped samlpes 
    counter=1
    return_value=0 # exit code 
    no_of_files=`ls | wc -l ` # gets number of files 
    if [ $((no_of_files)) -eq 0 ] # checks if samples are not pairs 
    then 
        return_value=0
        echo -e "${RED}ERROR${NC}"
        exit
    else
        if [ `ls *"${check_extension}" 2> /dev/null | wc -l ` -eq 0 ];then 
                return_value=0
                echo -e "${RED}ERROR${NC}"
                exit
        else
            return_value=1
        fi
    fi
    return $return_value 
}

# trimmomatic paired end 
trimmer_pe () { # takes R1 file 
    # prepare file names 
    f1=$1
    f2=$(echo $f1 | sed 's/R1/R2/')
    logname=${f1%'.fastq.gz'}'.log'
    summaryname=${f1%'.fastq.gz'}'.summry'
    
    new_f1=${f1%'.fastq.gz'}'.trim.fastq.gz'
    newf1=${new_f1}
    newf_2=${f2%'.fastq.gz'}'.trim.fastq.gz'
    newf2=${newf_2}

    newf_1U=${f1%'.fastq.gz'}'.trim.se.fastq.gz'
    newf1U=${newf_1U}
    newf_2U=${f2%'.fastq.gz'}'.trim.se.fastq.gz'
    newf2U=${newf_2U}
    
    adap="$CONDA_PREFIX/share/trimmomatic-0.39-1/adapters"
    
    # performes trimmomatic on sample with 20 threads, no adaptor trimming 
    # to cut adaptor modify adap to proper adaptor and uncomment ILLUMINACLIP line 
    trimmomatic PE -threads 20 -phred33 -trimlog ${logname} \
            -summary ${summaryname}  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
            SLIDINGWINDOW:4:10 MINLEN:30 \
            # ILLUMINACLIP:$adap/TruSeq3-PE-2.fa:2:30:10:1
    ## PE -> paired ended
    ## SLIDINGWINDOW: Performs a sliding window trimming approach.
    ## ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads
}

# Trimmomatic single ended
trimmer () { # takes R1 file 
    f1=$1
    # prepare file names 
    logname=${f1%'.fastq.gz'}'.log'
    summaryname=${f1%'.fastq.gz'}'.summry'
    
    new_f1=${f1%'fastq.gz'}'trim.fastq.gz'
    newf1=${new_f1}
 
    adap="$CONDA_PREFIX/share/trimmomatic-0.39-1/adapters"
    # performes trimmomatic on sample with 20 threads, no adaptor trimming 
    # to cut adaptor modify adap to proper adaptor and uncomment ILLUMINACLIP line
    trimmomatic SE -threads 20 -phred33 \
                -trimlog ${logname} \
                -summary ${summaryname}  $f1  $newf1 \
                SLIDINGWINDOW:4:10 MINLEN:30 \
                # ILLUMINACLIP:$adap/TruSeq3-SE.fa:2:30:10

}

# hisat2 for paried ended reads 
runhisat2 () { # takes argumetns : 1. threads 2. sample R1 
    R1=$2
    R2=$(echo $R1 | sed 's/R1/R2/')
    hisat2 -p $1 -x $hisat2_index -1 $R1 -2 $R2 | \
            samtools view -bS | \
            samtools sort -o ${R1%'.fastq.gz'}'.sorted.bam' - 
            # generates bam files directly 
}
# hisat2 for single ended reads 
run_SE_hisat2 () { # takes argumetns : 1. threads 2. sample R1 
    R1=$2
    hisat2 -p $1 -x $hisat2_index -U $R1 | \
            samtools view -bS | \
            samtools sort -o ${R1%'.fastq.gz'}'.sorted.bam' - 
            # generates bam files directly
}

# bwa for paired ended reads 
bwa_rm_human () { # takes argumetns : 1. threads 2. sample R1 
    R1=$2
    R2=$(echo $R1 | sed 's/R1/R2/')
    bwa mem -t $1 $bwa_index $R1 $R2 | \
            samtools view -bS | \
            samtools sort -o ${R1%'.fastq.gz'}'.sorted.bam' - 
}

# bwa for single ended reads 
bwa_SE_rm_human () { # takes argumetns : 1. threads 2. sample R1 
    R1=$2
    bwa mem -t $1 $bwa_index $R1 | \
            samtools view -bS | \
            samtools sort -o ${R1%'.fastq.gz'}'.sorted.bam' - 
}

# bwa for paired ended reads 
runbwa () { # takes argumetns : 1. threads 2. sample R1 
    R1=$2
    R2=$(echo $R1 | sed 's/R1/R2/')
    bwa mem -t $1 $3 $R1 $R2 | \
            samtools view -bS | \
            samtools sort -o ${R1%'.fastq.gz'}'.sorted.bam' - 
}

# bowtie2 for paired ended reads 
runbowtie2 (){ # takes argumetns : 1. threads 2. sample R1 
    R1=$2
    R2=$(echo $R1 | sed 's/R1/R2/')
    bowtie2 -p $1 --very-sensitive-local -q -x $bowtie2_index -1 $R1 -2 $R2 | \
            samtools view -bS | \
            samtools sort -o ${R1%'.fastq.gz'}'.sorted.bam' - 
}

# bowtie2 for single ended reads
run_SE_bowtie2 (){ # takes argumetns : 1. threads 2. sample R1 
    R1=$2
    bowtie2 -p $1 --very-sensitive-local -q -x $bowtie2_index -U $R1 | \
            samtools view -bS | \
            samtools sort -o ${R1%'.fastq.gz'}'.sorted.bam' - 
}

# extracts unmapped reads from paired ended aligned files
extract_unmapped () { # takes 3 arguments: 1. file in (BAM) 2. file out 3. file name (intermediate)
    filein=$1 
    fileout=$2
    filename=$3 
    samtools view -b -f 12 -F 256 $filein > $fileout
    samtools sort -n -m 5G -@ $4 $fileout -o $fileout.sorted.bam
    samtools fastq -@ $4 $fileout.sorted.bam \
            -1 "$filename"_HOST_REMOVED_R1.fastq.gz \
            -2 "$filename"_HOST_REMOVED_R2.fastq.gz \
            -0 /dev/null -s /dev/null -n
}

# extracts unmapped reads from single ended aligned files 
extract_SE_unmapped () { # takes 2 arguments 1. file in (BAM) 2. file out (fastq)
    filein=$1 
    fileout=$2
    samtools fastq -f 4 $filein > $fileout
}

# kraken paired end analysis 
kraken_ () { # takes 3 arguments 1. number of threads 2. report name 3. file in (fastq)
    numofT=$1
    reportname=$2'.report'
    R1=$3
    R2=$(echo $R1 | sed 's/R1/R2/')
    kraken2 --db $kraken_db --threads $numofT \
            --report $reportname \
            --paired $R1 $R2 > ${2%'.report'}'.out'
}

# kraken single end analysis 
kraken_SE () { # takes 3 arguments 1. number of threads 2. report name 3. file in (fastq)
    numofT=$1
    reportname=$2'.report'
    R1=$3
    kraken2 --db $kraken_db --threads $numofT \
            --report $reportname \
             $R1 > ${2%'.report'}'.out'

}

# generates korona figure 
kraken2korona () { # takes name of output file and kraken output (OUT) file
    output=$1
    samples=$2
    ktImportTaxonomy -q 2 -t 3 \
                    -o $output \
                    $samples 
}

# bowtie alignment against mirbase for miRNA analysis 
bowtie_mirbase_count (){ # takes number of threads and file in (R1)
    R1=$2
    bowtie -n 0 -l 32 --norc --best --strata -m 1 --threads $1 $bowtie_index_mirbase $R1 --un ${R1%'.fastq'}'_unmapped.fastq' -S ${R1%'.fastq'}'.miRBase.sam'
    samtools sort --output-fmt BAM ${R1%'.fastq'}'.miRBase.sam' > ${R1%'.fastq'}'.miRBase.sorted.bam'
    samtools index ${R1%'.fastq'}'.miRBase.sorted.bam'
    samtools idxstats ${R1%'.fastq'}'.miRBase.sorted.bam' | cut -f1,3 - >  ${R1%'.fastq'}'.counts.txt'
}

# bowtie alignemnt against human for miRNA analysis 
bowtie_mirhuman (){ # takes number of threads and file in (R1)
    R1=$2
    bowtie -n 1 -l 32 --norc --best --strata -m 1 --threads $1 $bowtie_index_hum $R1  -S ${R1%'.fastq'}'.hum.sam'
    samtools sort --output-fmt BAM ${R1%'.fastq'}'.hum.sam' > ${R1%'.fastq'}'.hum.sorted.bam'
    samtools index ${R1%'.fastq'}'.hum.sorted.bam'
    bedtools tag -i ${R1%'.fastq'}'.hum.sorted.bam' -files $miRBED -names -tag XQ > ${R1%'.fastq'}'_hum.tagged.bam'
}

# single ended cutadapt for miRNA analysis 
cutadapt_mirna(){
    target=$1
    cutadapt -a TGGAATTCTCGGGTGCCAAGG -o ${target%'.fastq'}'.temp.trim.fastq' --minimum-length 23 $target >> $outputdit/LOG.txt 2>&1
    cutadapt --minimum-length=12 --maximum-length=35 -u 4 -u -4 -o ${target%'.fastq'}'.trim.fastq' ${target%'.fastq'}'.temp.trim.fastq'
    rm *'.temp.'*
}

# bowtie alignemnt against human for miRNA analysis 
bowtie_mirSecond (){ # takes number of threads and file in (R1)
    R1=$2
    bowtie --threads $1 -q -v 0 -k 10 -S -t $bowtie_index_hum $R1  -S ${R1%'.fastq'}'.hum.sam'
}
error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}
