help_message_WGS = """GUAP WGS analyis pipeline:

\033[;33;1musage:\033[;39;m GUAP WGS -i|--input PATH -o|--output PATH  

___________________________________________________________

\033[;33;1mbasic options:\033[;39;m

 -i, --input       DIR     input directory 
 -o, --output      DIR     output directory


\033[;33;1mconfigurations:\033[;39;m

 --skip-QC                     skipps fastqc step
 --skip-assembly               skipps assembly step
 --skip-assembly-QC            skipps aligning reads to assembly step
 --pre-check-contaminant       performs kraken2 on reads before mapping
 --skip-check-contaminant      skipps kraken2 of unmapped reads

\033[;36;1mDecontamination:\033[;39;m
 --kraken2-index       FILE    path to kraken2 index 

\033[;36;1mPerformance:\033[;39;m
 -t, --threads         INT     Number of total threads [if not supplied will use all available threads]
 --threads-index       INT     Number of threads to use during indexing ref [default = All given threads]
 --threads-assemble    INT     Number of threads to use per sample assembly [default = 4]
 --threads-align       INT     Number of threads to use per sample align [default = 4]

\033[;36;1mReference Indexing:\033[;39;m
 --index               FILE    fasta file of reference genome
 --index-name          STR     name of index to use
 --index-path          DIR     path to index bowtie2 files, use if pre-indexed 
 \033[;34;1mReference Downloading(use to download and index):\033[;39;m
 --org-name            STR     Scientific name of organism to search and download from NCBI
 --org-txid            STR     Organism NCBI taxid to search and download 
 --org-id              STR     Organism NCBI id to search and download
 --org-download-direct         Download first found record automaticaly 
 --org-auto-download   STR     Uses kraken2 to identify top abundant species and search NCBI and download ref.
    \033[;33;1m--kraken2-index [should be set when using --org-auto-download]\033[;39;m

\033[;36;1mOther configs:\033[;39;m
 --continue                    Continue on last stopped process
 --overwrite                   Overwrites analysis dir if found 
 
  
  this workflow is still under development thank you for testing for any bugs;
    please contact: Abdelrhman @ abdelrhman.yahia@57357.org
"""

help_message_WES = """GUAP WES analyis pipeline:

\033[;33;1musage:\033[;39;m GUAP WES -i|--input PATH -o|--output PATH --aligner bwa|bowtie2 \\ 
                --reference-fasta PATH --bed-file PATH

___________________________________________________________

\033[;33;1mbasic options:\033[;39;m
 -i, --input       DIR     input directory 
 -o, --output      DIR     output directory

\033[;33;1mconfigurations:\033[;39;m
 --skip-QC                       skipps fastqc step
 --aligner        bwa|bowtie2    choose aligner from [bwa|bowtie2]
 --variant-caller GATK | mpileup            choose variant caller [GATK|mpileup]
 --known-variants PATH         Path to known variants file (vcf) 
 --bed-file       PATH         Path to regeions file (bed)
 -n, --name       STR          name of the analysis

\033[;36;1mIndexing reference:\033[;39;m
 --index-reference               Index the reference using  choosen aligner 
 --reference-fasta     PATH      Path to reference fa file 
 --reference-index     PATH      Path and prefix of index files of proper aligner

\033[;36;1mPerformance:\033[;39;m
 -t, --threads         INT       Number of total threads [if not supplied will use all available threads]
 --threads-index       INT       Number of threads to use during indexing ref [default = All given threads]
 --threads-align       INT       Number of threads to use per sample align [default = 4]
 --threads-calling     INT       Number of threads to use per sample variant calling [default = 4]
 --bash                          uses bash scripts instead of snakemake (slower)

\033[;36;1mSnakemake options:\033[;39;m
 --snakemake-dag           Exports Rules DAG for the analysis
 --snakemake-dry-run       performs snakemake dry run

\033[;36;1mOther configs:\033[;39;m
  --verbose                print many output (most effective in snakemake mode)
 --continue                Continue on last stopped process
 --overwrite               Overwrites analysis dir if found 
 
  
\033[;36;1mExample run:\033[;39;m
guap WES \033[;35;1m-i\033[;39;m samples \033[;35;1m-o\033[;39;m out \033[;35;1m-n\033[;39;m test_1 \\
         \033[;35;1m--aligner\033[;39;m bwa \033[;35;1m--variant-caller\033[;39;m GATK \\
         \033[;35;1m--known-variants\033[;39;m exome.hg38.vcf.gz \\
         \033[;35;1m--reference-fasta\033[;39;m hg38.fa.gz \\
         \033[;35;1m--reference-index\033[;39;m bwa/hg38/hg38 \\
         \033[;35;1m--threads-align\033[;39;m 12 \033[;35;1m--threads-calling\033[;39;m 12 \\
         \033[;35;1m--bed-file\033[;39;m exon_coding_seq.bed \\
         \033[;35;1m--overwrite --snakemake-dry-run\033[;39;m

  this workflow is still under development thank you for testing for any bugs;
    please contact: Abdelrhman @ abdelrhman.yahia@57357.org
"""

help_message_RNA = """GUAP RNA analyis pipeline:

\033[;33;1musage:\033[;39;m guap RNA -i|--input PATH -o|--output PATH --aligner star|hisat2 \\ 
                --reference-fasta PATH --gtf-file PATH

___________________________________________________________

\033[;33;1mbasic options:\033[;39;m
 -i, --input       DIR     input directory 
 -o, --output      DIR     output directory

\033[;33;1mconfigurations:\033[;39;m
 --skip-QC         skipps fastqc step
 --aligner         choose aligner from [star|hisat2|salmon|kallisto|rsem]
 --quantifier      choose quantifier [featurecounts|htseq]
 --gtf-file        Path to regeions file (gff/gtf/gff3)
 -n, --name        name of the analysis

\033[;36;1mIndexing reference:\033[;39;m
 --reference-fasta     PATH      Path to reference fa file 
 --reference-index     PATH      Path and prefix of index files of proper aligner

\033[;36;1mPerformance:\033[;39;m
 -t, --threads         INT       Number of total threads [if not supplied will use all available threads]
 --threads-index       INT       Number of threads to use during indexing ref [default = All given threads]
 --threads-qc          INT       Number of threads to use per sample QC [default = 4]
 --bash                          uses bash scripts instead of snakemake (slower)

\033[;36;1mSnakemake options:\033[;39;m
 --snakemake-dag           Exports Rules DAG for the analysis
 --snakemake-dry-run       performs snakemake dry run

\033[;36;1mOther configs:\033[;39;m
  --verbose                print many output (most effective in snakemake mode)
 --continue                Continue on last stopped process
 --overwrite               Overwrites analysis dir if found 

  this workflow is still under development thank you for testing for any bugs;
    please contact: Abdelrhman @ abdelrhman.yahia@57357.org
"""

help_message_16s = """GUAP 16s analyis pipeline:

\033[;33;1musage:\033[;39;m GUAP 16s [-i|--input PATH] [-o|--output PATH] 
                [-m|--metadata FILE] [-c|--classifier FILE]

___________________________________________________________

\033[;33;1mbasic options:\033[;39;m

 -i, --input       DIR     input directory 
 -o, --output      DIR     output directory
 -m, --metadata    FILE    metadata file 
 -t, --threads     INT     Number of threads to use 
 -c, --classifier  FILE    classifier file (should be either
                           .qza file use with QIIME2 or 
                           fasta file for use with DADA2)


-cmt, --create-metadata-table 	create metadata empty 
				table from sample names in the input dir


\033[;33;1mconfigurations:\033[;39;m
\033[;36;1mbasic:\033[;39;m
 --bash                    uses bash scripts instead of snakemake
 --downstream              performs some downstream analysis in QIIME2 
   -d, --sampling-depth  INT   Sampling Depth for core-mertic phylogenetic
   -md, --max-depth      INT   Maximum Sampling Depth for alpha-rarefaction
 --export-figs             exports phyloseq basic figures 
    requires:                                   
  --condition-name STR     condition name in metadata file
 -n, --name        STR     name of the analysis
 --verbose                 print many output (most effective in snakemake mode)
 --continue                continue the analysis when re-run 
 --snakemake-dag           Exports Rules DAG for the analysis
 --snakemake-dry-run       performs snakemake dry run

\033[;36;1mQC:\033[;39;m
 --skip-QC                     skipps fastqc step

 --skip-trimmomatic            skipps trimmomatic step
 --trim-min-length      INT    trimmomatic minimum length of read

 --remove-primers              perform cutadapt to remove primers
                                  (requires: -fp [forward primer sequence], 
                                             -rp [reverse primer sequence])
 --min-length           INT    cutadapt min length [default=50]

\033[;36;1mASV generation:\033[;39;m
DADA2 
 -tf,--trunc-f          INT    trunclength Forward [default=0]
 -tr,--trunc-r          INT    trunclength Reverse [default=0]
 -l,--trim-l            INT    trunclength Forward [default=0]
 -r,--trim-r            INT    trunclength Reverse [default=0]
 -ef, --maxee-f         INT    maxEE Forward [default=4]
 -er,--maxee-r          INT    maxEE Reverse [default=5]
 -mo,--min-overlap      INT    minimum overlap to merge [default=10]
 -cm, --chimera-method  [consensus|pooled] chimera method [default=consensus]

 --deblur                       uses deblur to generate ASV table
 -L, --deblur-trim-length INT   Deblur trim length (required with --deblur)

 --use-QIIME2                   uses DADA2 in QIIME2 rather than R

\033[;36;1mClassification:\033[;39;m
 --choose-classifier    [dada|qiime]    Choose to classify taxanomy defualt: qiime
 --train    trainset for qiime2 naive classifier [can't be used with (dada classifier))]
   requires:
 --t-i-seqs            FILE    representivie sequences file (qza or FASTA)
 --t-i-taxa            FILE    taxonomy file (qza or csv)
 --trainset-minlength  INT     length of amblicons to be exclude if lower   
 --trainset-maxlength  INT     length of amblicons to be exclude if greater
 -fp, --forward-primer STR     forward primer sequence (used with cutadapt and train set)
 -rp, --reverse-primer STR     reverse primer sequence (used with cutadapt and train set)
                                   

\033[;33;1mexample usage:\033[;39;m 
GUAP 16s -i indir -o outdir -m metadata.tsv -c classifier.qza \\
         -tf 220 -tr 170 -l 10 -r 10 -ef 5 -er 7 --downstream \\
         --export-figs --condition-name condition

  this workflow is still under development thank you for testing for any bugs;
    please contact: Abdelrhman @ abdelrhman.yahia@57357.org
"""
