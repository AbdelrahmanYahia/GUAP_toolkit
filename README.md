# GUAP

**G**enomics **U**nit **A**nalysis **P**ipelines (aka. **GUAP**) is a collection of bash scripts that performs several pipelines usually conducted as a routine NGS work. The script works via command line interface. the user can give a samples directory containing all samples and choose the desired analysis pipeline and the program will perform all the upstream analysis. The supported pipelines are:

- Metagenomics pipeline ( metagenomes identification using kraken)
- miRNA analysis 1 ( standard genome alignment pipeline of miRNA)
- miRNA analysis 2 ( supervised miRNA analysis pipeline)

---

- List of used software:
    - Kraken v.1.1.1
    - KronaTools v.2.7.1 - ktImportTaxonomy
    - trimmomatic v.0.39
    - cutadapt v.2.8
    - hisat2  v.2.2.1
    - bwa v.0.7.17-r1188
    - bowtie v.1.2.3
    - bowtie2 v.2.4.1
    - samtools v.1.10
    - htslib v.1.10.2
    - bedtools v.2.27.1
    - python v.3.8.5
    - pandas v.1.2.2
- Scripts

    ```bash
    ├── functions.sh
    ├── GUAP
    ├── merge.py
    └── pipelines.sh
    ```

---

## Main workflow

GUAP takes several options as shown in the following block:

```bash
__usage="
GUAP v1.0
Usage: $(basename $0) [OPTIONS]

GUAP -i <inputdir> -o <outputdir> -c <analysis> -a <aligner> -t INT 

Options:
  -c <miRNA|kraken>                     Choose analysis                default = kraken
  -a <bwa|bowtie|bowtie2|hisat2>        Choose aligner                 default = bowtie2
  -i <str>                              Input directory path
  -o <str>                              Output directory path
  -t <int>                              Number of threads              default = 50
  -m <first|second>                     miRNA pipeline choise
  -e <bash|snakemake|nextflow>          use snakemake, nextflow or bash
  -h                                    Help message ( This message )
  -y                                    skip continue check
  -k                                    Direct flag
"
```


### Metagenomics pipeline using kraken:

Run the command:

`$ ./GUAP -i inputdir -o outputdir -a bowtie2 -c kraken -t num_of_threads` 

The script will first check the number of samples and if they are paired ended or not using the function `sample_check` from `functions.sh` according to exit code of the function the pipeline will continue with SE or PE analysis, the example provided upove is with single ended analysis.

The main analysis starts by trimming samples using `trimmomatic` then aliging samples to refference human genome according user defined aligner, unmapped reads will then be extract from sorted bam file using `samtools` the unmapped reads will then be used with `kraken` with k2_pluspfp database, korona figure will then be generated from output file using `ktImportTaxonomy` 



### miRNA tradiotinal pipeline

This pipeline is the traditional pipeline for miRNA analysis 

Run the command:

`$ ./GUAP -i inputdir -o outputdir -a bowtie -c mirna -t num_of_threads -m second` 

The script will first check the number of samples and if they are paired ended or not using the function `sample_check` from `functions.sh` according to exit code of the function the pipeline will continue with SE or PE analysis, the example provided upove is with single ended analysis. 

The main analysis starts by trimming samples using `cutadapt` then aliging samples to refference human genome using `bowtie` the resulted sam files of all samples will then be used to generate counts using `featureCounts` and file will be modified by the command `tail -n +2 all.counts | cut -f 1,7- > all.txt` .



### miRNA supervised pipeline

This pipeline is the subervised pipeline for miRNA analysis 

Run the command:

`$ ./GUAP -i inputdir -o outputdir -a bowtie -c mirna -t num_of_threads -m first`

The script will first check the number of samples and if they are paired ended or not using the function `sample_check` from `functions.sh` according to exit code of the function the pipeline will continue with SE or PE analysis, the example provided upove is with single ended analysis. 

The main analysis starts by trimming samples using `cutadapt` then aliging samples to samples using `bowtie` to mirbase database with options `bowtie -n 0 -l 32 --norc --best --strata -m 1 --threads $1 $bowtie_index_mirbase <input> --un <unmapped fastaname> -S <sam file>` then convert sam to sorted indexed bam file using `samtools` the bam file will then be xounted using `samtools idxstats | cut-f1,3` to counts file. The unmapped reads wil then be aligned using `bowtie` to human genome using the command `bowtie -n 1 -l 32 --norc --best --strata -m 1 --threads $1 $bowtie_index_hum <input> -S <output>` then the sam file will be converted to bam and tagged using the command `bedtools tag -i <input> -files $hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed -names -tag XQ > <tagged_bam>` the output will be a tagged bam file. to count the file a list of all miRNA names is looped and counted each one using the command `samtools view $i | grep $j | wc -l` which `i` is the sample name and `j` is miRNA name. 

All counts files from mirbase and mapped to genome are then merged together using `merge.py` python script. This script takes the names of the samples and creates a dictionary, with keys > sample names and values > dictionary with keys > mirna name and value > counts data. The data are filled by firstly looping on all mirbase files and generating counts then looping on all genome files and adding the value of each mirna name to there corresponding values. the resulted dictionary is then exported as a CSV file with all counts data. 

---

## Datasets used:

1. Kraken database: k2_pluspfp_20210127
2. bowtie human index
3. bowtie mirbase index
4. bwa human index (GRCH38)
5. bowtie2 human index (GRCH38)
6. hisat2 human index (GRCH38)
7. mirbase annotation: human gff3 file 
8. miRNA names list

