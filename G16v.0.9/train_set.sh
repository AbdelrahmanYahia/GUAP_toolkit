qiime feature-classifier extract-reads \
  --i-sequences /media/genomics/DB_Storage/db/16s-dbs/QIIME2_silva/silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-min-length 350 \
  --p-max-length 580 \
  --o-reads ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy /media/genomics/DB_Storage/db/16s-dbs/QIIME2_silva/silva-138-99-tax.qza \
  --o-classifier Silva_classifier.qza