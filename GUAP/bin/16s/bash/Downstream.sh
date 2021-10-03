echo -e "${YEL}#############     summarize     ############${NC}"
qiime feature-table summarize \
  --i-table ./table.qza \
  --o-visualization visualization/table.qzv \
  --m-sample-metadata-file ../sample-metadata.tsv &

qiime feature-table tabulate-seqs \
  --i-data ./rep-seqs.qza \
  --o-visualization visualization/rep-seqs.qzv &

echo -e "${YEL}#############     taxa plot     ############${NC}"
qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy classify/taxonomy.qza \
    --m-metadata-file ../sample-metadata.tsv \
    --o-visualization visualization/bar_plot.qzv &

qiime tools export --input-path classify/taxonomy.qza --output-path "./classify/taxonomy" &

echo -e "${YEL}###############    tabulate   ##############${NC}"
qiime metadata tabulate \
  --m-input-file classify/taxonomy.qza \
  --o-visualization visualization/taxonomy.qzv &

echo -e "${YEL}#######  diversity alpha-rarefaction  ######${NC}"
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny align/rooted-tree.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --p-max-depth ${alphararefaction_maxdepth} \
  --o-visualization visualization/alpha-rarefaction.qzv &

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --o-visualization visualization/shannon_vector-group-significance.qzv &

wait
error_cheker $?


echo -e "${YEL}########  alpha-group-significance  ########${NC}"
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --o-visualization visualization/faith-pd-group-significance.qzv &

echo -e "${YEL}########  evenness-group-significanc  ######${NC}"
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --o-visualization visualization/evenness-group-significance.qzv &

echo -e "${YEL}##### unweighted-unifrac-subject-group #####${NC}"
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --m-metadata-column condition \
  --o-visualization visualization/unweighted-unifrac-condition-significance.qzv \
  --p-pairwise &

echo -e "${YEL}###########      emperor plot     ##########${NC}"
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --p-custom-axes \
  --o-visualization visualization/unweighted-unifrac-emperor-days-since-experiment-start.qzvm &

echo -e "${YEL}##########    alpha and beta    ############${NC}"
qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file ../sample-metadata.tsv \
  --o-visualization visualization/bray-curtis-emperor-days-since-experiment-start.qzv 

wait
error_cheker $?

echo -e "${YEL}#############     exports    ###############${NC}"

qiime tools export  \
  --input-path visualization/alpha-rarefaction.qzv \
  --output-path exports/alpha-rarefaction &

qiime tools export  \
  --input-path visualization/evenness-group-significance.qzv \
  --output-path exports/evenness-group-significance &

qiime tools export  \
  --input-path visualization/faith-pd-group-significance.qzv \
  --output-path exports/faith-pd-group-significance &

qiime tools export  \
  --input-path visualization/shannon_vector-group-significance.qzv \
  --output-path exports/shannon_vector-group-significance &
wait

error_cheker $?
