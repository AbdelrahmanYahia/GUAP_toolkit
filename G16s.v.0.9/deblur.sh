mkdir -p deblur
mkdir -p demux
cd ../
echo -e "${YEL}###########      QIIME import    ###########${NC}"
# qiime tools import \
#   --type 'SampleData[PairedEndSequencesWithQuality]' \
#   --input-path samples.tsv \
#   --output-path QIIME2/demux/paired-end-demux.qza \
#   --input-format PairedEndFastqManifestPhred33V2

cd QIIME2

# # demultiplexing
# qiime demux summarize \
#   --i-data demux/paired-end-demux.qza \
#   --o-visualization visualization/demux.qzv &

echo -e "${YEL}###########    Quality filter    ###########${NC}"

qiime quality-filter q-score \
  --i-demux demux/paired-end-demux.qza \
  --o-filtered-sequences demux/demux-filtered.qza \
  --o-filter-stats demux/demux-filter-stats.qza


echo -e "${YEL}################     ASV     ###############${NC}"

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux/demux-filtered.qza \
  --p-trim-length ${deblur_trimlength} \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --p-sample-stats \
  --o-stats deblur/deblur-stats.qza

qiime metadata tabulate \
  --m-input-file demux/demux-filter-stats.qza \
  --o-visualization visualization/demux-filter-stats.qzv &

qiime deblur visualize-stats \
  --i-deblur-stats deblur/deblur-stats.qza \
  --o-visualization visualization/deblur-stats.qzv

wait