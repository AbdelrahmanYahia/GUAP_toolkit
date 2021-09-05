
echo -e "${YEL}###############     aldex2     #############${NC}"
qiime aldex2 aldex2 \
  --i-table table.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column ${DAT_metaDcolumn} \
  --output-dir DAT

echo -e "${YEL}############     effect-plot    ############${NC}"
qiime aldex2 effect-plot \
  --i-table DAT/differentials.qza \
  --o-visualization DAT/DATv

echo -e "${YEL}###########  extract-difference  ###########${NC}"
qiime aldex2 extract-differences \
  --i-table DAT/differentials.qza \
  --o-differentials DAT/sig \
  --p-sig-threshold 0.1 \
  --p-effect-threshold 0 \
  --p-difference-threshold 0

qiime tools export \
  --input-path DAT/sig.qza \
  --output-path differentially-expressed-features
