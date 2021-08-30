
echo -e "${RED}#############      classify     ############${NC}"
qiime feature-classifier classify-sklearn \
  --i-classifier ../Silva_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification classify/taxonomy.qza &

echo -e "${RED}###############     align     ##############${NC}"
qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences rep-seqs.qza \
            --o-alignment align/aligned-rep-seqs.qza \
            --o-masked-alignment align/masked-aligned-rep-seqs.qza \
            --o-tree align/unrooted-tree.qza \
            --o-rooted-tree align/rooted-tree.qza --p-n-threads ${threads} &
wait

echo -e "${RED}########  core-metrics-phylogenetic  #######${NC}"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny align/rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth ${core_met_phylo_samplingdepth} \
  --m-metadata-file ../sample-metadata.tsv \
  --output-dir core-metrics-results 
