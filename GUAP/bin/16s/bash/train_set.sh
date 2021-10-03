qiime feature-classifier extract-reads \
  --i-sequences  ${indexes_trainclassifier_seqs}\
  --p-f-primer ${trainset_fprimer} \
  --p-r-primer ${trainset_rprimer} \
  --p-min-length ${trainset_minlength} \
  --p-max-length ${trainset_maxlength} \
  --o-reads ref-seqs.qza
error_cheker $?
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ${indexes_trainclassifier_taxa} \
  --o-classifier ${indexes_classifier}
error_cheker $?