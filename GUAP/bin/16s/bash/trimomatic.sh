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
