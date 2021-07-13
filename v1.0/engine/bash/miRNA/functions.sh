source ../common/common_functions.sh
source ../common/indexes.sh

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