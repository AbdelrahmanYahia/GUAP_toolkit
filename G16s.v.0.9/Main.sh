eval "$(conda shell.bash hook)"
source parse_yml.sh
eval $(parse_yaml config.yml)
RED='\x1b[1;31m'
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
NC='\e[0m'
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

error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}


conda deactivate
############  R part ###########
echo -e "${YEL}########    DADA2 Rscript running   ########${NC}"
Rscript G16s.v0.9.R \
        -i "${PWD}/samples" -o "${PWD}/DADA2"\
        -p ${threads} \
        -t ${dada_trunclength_f} \
        -T ${dada_trunclength_r} \
        -e ${dada_maxEE_f} \
        -E ${dada_maxEE_r} -s > R.out

###############################
# echo -e "${YEL}###########      QIIME import    ###########${NC}"
# conda activate qiime2-2021.4

# mkdir -p QIIME2
# cd QIIME2
# mkdir -p phyloseq
# mkdir -p classify
# mkdir -p visualization
# mkdir -p align

# qiime tools import \
# --input-path ../DADA2/rep-seqs.fna \
# --type 'FeatureData[Sequence]' \
# --output-path rep-seqs.qza

# echo -n "#OTU Table" | cat - ../DADA2/seqtab-nochim.txt > biom-table.txt

# biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5

# qiime tools import \
# --input-path table.biom \
# --type 'FeatureTable[Frequency]' \
# --input-format BIOMV210Format \
# --output-path table.qza
# ################################
# echo -e "${YEL}#############     summarize     ############${NC}"
# qiime feature-table summarize \
#   --i-table table.qza \
#   --o-visualization visualization/table.qzv \
#   --m-sample-metadata-file ../sample-metadata.tsv &

# qiime feature-table tabulate-seqs \
#   --i-data rep-seqs.qza \
#   --o-visualization visualization/rep-seqs.qzv &

# source ../Process.sh 

# echo -e "${YEL}###########    export biom     #############${NC}"
# # Export OTU table
# qiime tools export \
#         --input-path table.qza \
#         --output-path phyloseq &

# echo -e "${YEL}#############    biom to csv    ############${NC}"
# # Convert biom format to tsv format
# biom convert \
#         -i ./phyloseq/feature-table.biom \
#         -o ./phyloseq/otu_table.tsv \
#         --to-tsv


# cd phyloseq
# sed -i '1d' otu_table.tsv
# sed -i 's/#OTU ID//' otu_table.tsv
# cd ../


# echo -e "${YEL}###### Export representative sequences #####${NC}"
# # Export representative sequences
# qiime tools export \
#         --input-path rep-seqs.qza \
#         --output-path phyloseq


# echo -e "${YEL}##########     export tree files    ########${NC}"
# # Export tree files
# qiime tools export \
#         --input-path align/unrooted-tree.qza \
#         --output-path phyloseq


# cd phyloseq
# mv tree.nwk unrooted_tree.nwk
# cd ../

# echo -e "${YEL}############     rooted tree    ############${NC}"
# qiime tools export \
#         --input-path align/rooted-tree.qza \
#         --output-path phyloseq

# cd phyloseq
# mv tree.nwk rooted_tree.nwk
# cd ../


# continue_
# source ../Downstream.sh
# cd core-metrics-results
# for i in *.qzv
# do
#   mv $i ../visualization &
# done
# wait
echo -e "${GRE}###################  Done  #################${NC}"
