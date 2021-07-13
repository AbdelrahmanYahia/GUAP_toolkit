# file containg piplines 
# opens functions script file
source ../common/common_functions.sh
source ../common/indexes.sh

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
            --paired $R1 $R2 > ${2%'.report'}'.out' \
            --unclassified-out ${R1}"#.unclassified"
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
