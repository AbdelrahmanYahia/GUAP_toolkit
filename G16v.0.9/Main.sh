eval "$(conda shell.bash hook)"
source parse_yml.sh
eval $(parse_yaml config.yml)
RED='\x1b[1;31m'
NC='\e[0m'
conda deactivate
############  R part ###########
echo -e "${RED}########    DADA2 Rscript running   ########${NC}"
Rscript G16s.v0.9.R \
        -i "${PWD}/samples" -o "${PWD}/DADA2"\
        -p ${threads} \
        -t ${dada_trunclength_f} \
        -T ${dada_trunclength_r} \
        -e ${dada_maxEE_f} \
        -E ${dada_maxEE_r} -s > R.out

###############################
echo -e "${RED}###########      QIIME import    ###########${NC}"
conda activate qiime2-2021.4

mkdir -p QIIME2
cd QIIME2
mkdir -p phyloseq
mkdir -p classify
mkdir -p visualization
mkdir -p align

qiime tools import \
--input-path ../DADA2/rep-seqs.fna \
--type 'FeatureData[Sequence]' \
--output-path rep-seqs.qza

echo -n "#OTU Table" | cat - ../DADA2/seqtab-nochim.txt > biom-table.txt

biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5

qiime tools import \
--input-path table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path table.qza
################################
source ../Process.sh 

echo -e "${RED}###########    export biom     #############${NC}"
# Export OTU table
qiime tools export \
        --input-path table.qza \
        --output-path phyloseq &

echo -e "${RED}#############    biom to csv    ############${NC}"
# Convert biom format to tsv format
biom convert \
        -i phyloseq/feature-table.biom \
        -o phyloseq/otu_table.tsv \
        --to-tsv


cd phyloseq
sed -i '1d' otu_table.tsv
sed -i 's/#OTU ID//' otu_table.tsv
cd ../


echo -e "${RED}###### Export representative sequences #####${NC}"
# Export representative sequences
qiime tools export \
        --input-path rep-seqs.qza \
        --output-path phyloseq


echo -e "${RED}##########     export tree files    ########${NC}"
# Export tree files
qiime tools export \
        --input-path align/unrooted-tree.qza \
        --output-path phyloseq


cd phyloseq
mv tree.nwk unrooted_tree.nwk
cd ../

echo -e "${RED}############     rooted tree    ############${NC}"
qiime tools export \
        --input-path align/rooted-tree.qza \
        --output-path phyloseq

cd phyloseq
mv tree.nwk rooted_tree.nwk
cd ../

source ../Downstream.sh

cd core-metrics-results
for i in *.qzv
do
  mv $i ../visualization &
done
wait
echo -e "${RED}###################  Done  #################${NC}"
