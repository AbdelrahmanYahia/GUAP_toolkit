# configration for bash engine
kraken_db="/media/genomics/DB_Storage/db/Kraken2/k2_pluspfp_20210127"
bwa_index="/home/genomics/Documents/DR_Amany/new/BWA_Index/GRCh38_latest_genomic.fna"
hisat2_index="/media/genomics/old/Data/Genomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2_human_genome/genome"
bowtie2_index="/media/genomics/Data/old/Genomes/Homo_sapiens/UCSC/hg19/Sequence/bowtie2_index/genome"
bowtie_index_hum="/media/genomics/DB_Storage/db/miRNA/Index/human-Grch38-Ref/hg38"
bowtie_index_mirbase="/media/genomics/DB_Storage/db/miRNA/Index/miRBase/mat_hsa.fa"
mirbase_annotation="/media/genomics/DB_Storage/db/mirbase_hsa_gff/hsa.gff3"
miRBED="/media/genomics/DB_Storage/db/miRNA/Index/hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed"

##########################################################################################################






















#########################################################################################################
############################################# DON'T EDIT ################################################
#########################################################################################################

source Scripts/bash/functions.sh
source Scripts/bash/Kraken_direct.sh
source Scripts/bash/Kraken_rm_human.sh
source Scripts/bash/miRNA_1.sh
source Scripts/bash/miRNA_2.sh
source Scripts/bash/indexes.sh
