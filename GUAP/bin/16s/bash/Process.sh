# 
if [ "${trainset_train}" == "TRUE" ] || [ "${trainset_train}" == "True" ] || [ "${trainset_train}" == "T" ] || [ "${trainset_train}" == "true" ] 
then
  conda activate qiime2
  source train_set.sh
fi

echo -e "${YEL}#############      classify     ############${NC}"
qiime feature-classifier classify-sklearn \
  --i-classifier ${indexes_classifier} \
  --i-reads rep-seqs.qza --p-n-jobs ${threads} \
  --o-classification classify/taxonomy.qza &

echo -e "${YEL}###############     align     ##############${NC}"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment align/aligned-rep-seqs.qza \
  --o-masked-alignment align/masked-aligned-rep-seqs.qza \
  --o-tree align/unrooted-tree.qza \
  --o-rooted-tree align/rooted-tree.qza --p-n-threads ${threads} &
wait

echo -e "${YEL}########  core-metrics-phylogenetic  #######${NC}"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny align/rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth ${core_met_phylo_samplingdepth} \
  --m-metadata-file ../sample-metadata.tsv \
  --output-dir core-metrics-results