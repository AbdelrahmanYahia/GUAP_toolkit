#!/bin/bash

# color codes for command line 
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
CYN='\x1b[1;36m'
NC='\e[0m'

checker_array=()
__usage="
GUAP v1.0
Usage: $(basename $0) [OPTIONS]

${RED}GUAP WGS analyis pipeline:${NC}

\033[;33;1musage:\033[;39;m guap WGS -i|--input PATH -o|--output PATH  

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

${GRE}____________________________________________________________${NC}


${RED}GUAP WES analyis pipeline:${NC}

\033[;33;1musage:\033[;39;m guap WES -i|--input PATH -o|--output PATH --aligner bwa|bowtie2 \\ 
                --reference-fasta PATH --bed-file PATH


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

${GRE}___________________________________________________________${NC}

${RED}GUAP 16s analyis pipeline:${NC}

\033[;33;1musage:\033[;39;m guap 16s [-i|--input PATH] [-o|--output PATH] 
                [-m|--metadata FILE] [-c|--classifier FILE]


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
guap 16s -i indir -o outdir -m metadata.tsv -c classifier.qza \\
         -tf 220 -tr 170 -l 10 -r 10 -ef 5 -er 7 --downstream \\
         --export-figs --condition-name condition

${GRE}___________________________________________________________${NC}

${RED}GUAP DE:${NC}

Usage: guap DE -i [INPUT] -c [clinical] -S [sample] -C [control] -o [out-dir]
Performs DE on counts table

Options:
        -i CHARACTER, --input=CHARACTER
                input counts file

        -t, --tab
                input is a tab separated file

        -c CHARACTER, --clinical=CHARACTER
                input clincial file

        -S CHARACTER, --sample=CHARACTER
                Sample name in clinical

        -C CHARACTER, --control=CHARACTER
                Control name in clinical

        -d CHARACTER, --dirstr=CHARACTER
                Dir str to remove [default= OUT_miRNA_human_mapped_reads_]

        -r CHARACTER, --re=CHARACTER
                RE pattern to remove [default= _S([0-9]+)_L([0-9]+)_R(1|2)_([0-9]+).sam]

        -z NUMBER, --sample-zeros=NUMBER
                Allowed number of zeros in sample group [default= 1]

        -x NUMBER, --control-zeros=NUMBER
                Allowed number of zeros in control group [default= 1]

        -o CHARACTER, --out-dir=CHARACTER
                output file Directory [default= out]

        -n CHARACTER, --name=CHARACTER
                Name of analysis [default= GUAP-DE-out]

        -F, --Full-run
                Performe Full DE run [default]

        -p CHARACTER, --wd=CHARACTER
                input counts file

        -Q, --QC
                Perfom QC only

        -h, --help
                Show this help message and exit

${GRE}___________________________________________________________${NC}

${RED}GUAP RNA analyis pipeline:${NC}

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

${GRE}___________________________________________________________${NC}

${RED}GUAP 16s downstream GUI${NC}

Usage: guap 16s_downstream

${GRE}___________________________________________________________${NC}

${YEL}All availble workflows:${NC}
- WES
- WGS
- DE
- RNA
- 16s
- 16s_downstream
${GRE}___________________________________________________________${NC}

This workflow is still under development thank you for testing 
for any bugs please contact: Abdelrhman @ abdelrhman.yahia@57357.org
"

# function to print usage message and exits 
usage() {
    echo -e "$__usage" 1>&2; exit 1
}

# asks user to continue or not
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

step_checker(){
    tocheck=$2
    cont=$1
    checker_array=`cat "${OUTPUT}/.checker_array.txt"`
    if [ "${continue}" == "True" ] || [ "${continue}" == "true" ] ; then
        if [[ " ${checker_array[*]} " =~ " ${tocheck} " ]]; then
            :	
        else
            eval $(echo ${tocheck})
            if [[ $? -ne 0 ]]; then 
                echo -e "${RED}ERROR IN: ${BLE}${tocheck}${NC}"
                if [[ ${cont} -ne 1 ]]; then 
                    echo -e "${RED}GUAP STOPPED${NC}"
                    print_elapsed
                    exit 1
                else
                    :
                fi
            else
                checker_array+=${tocheck}
                checker_array+=" "
                echo -e "${checker_array[@]}" > "${OUTPUT}/.checker_array.txt"
            fi
        fi
    else
        eval $(echo ${tocheck})
        if [[ $? -ne 0 ]]; then 
            echo -e "${RED}ERROR IN: ${BLE}${tocheck}${NC}"
            print_elapsed
            if [[ ${cont} -ne 1 ]]; then 
                echo -e "${RED}GUAP STOPPED${NC}"
                exit 1
            else
                :
            fi
        else
            checker_array+=${tocheck}
            checker_array+=" "
            echo -e "${checker_array[@]}" > "${OUTPUT}/.checker_array.txt"
        fi
    fi
}

# checks sample extension and paired or not 

print_elapsed(){
  elapsed=$(( SECONDS - start_time ))
  fullelapsed=`eval "echo Elapsed time: $(date -ud "@$elapsed" +'%H:%M:%S')"`
  echo -e "${YEL}${fullelapsed}${NC}"
}
export -f print_elapsed

error_cheker(){
    lastexit=$1
    cont=$2
    name_=$3
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR IN: ${BLE}${name_}${NC}"
        print_elapsed
        if [[ ${cont} -ne 1 ]]; then 
            echo -e "${RED}GUAP STOPPED${NC}"
            exit 1
        else
            echo -e "$? passed ${name_}"
        fi
    fi
}


function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

get_sampling_depth(){
    if [[ ${sampling_depth} == "null" ]] ; then
		if [ "${deblur}" == "True" ] || [ "${deblur}" == "true" ] ; then 
			sampling_depth=`cat stats/stats.csv | rev | cut -d "," -f3 | rev | tail -n +2 | sort -n | head -1`
			echo -e "${YEL}WARNING: ${NC}Using minimum Depth as sampling depth; ${YEL}${sampling_depth}${NC}"
		elif [ "${use_QIIME2}" == "True" ] || [ "${use_QIIME2}" == "true" ] ; then 
			sampling_depth=`cat stats/stats.tsv | rev | cut  -f2 | rev | tail -n +3 | sort -n | head -1`
			echo -e "${YEL}WARNING: ${NC}Using minimum Depth as sampling depth; ${YEL}${sampling_depth}${NC}"
		elif [ -d "../DADA2" ] ; then
    		sampling_depth=`cat ../DADA2/stats.csv | rev | cut -d "," -f1 | rev | tail -n +2 | sort -n | head -1`
			echo -e "${YEL}WARNING: ${NC}Using minimum Depth as sampling depth; ${YEL}${sampling_depth}${NC}"
		else 
			echo -e "${RED}ERROR: ${NC}something worng! check params"
        	print_elapsed
        	exit
		fi
        if grep -Fq -- "sampling_depth: null" "${OUTPUT}/config.yaml"; then
            conf=`sed -r 's/^(\s*)(sampling_depth:\s*null\s*$)/\sampling_depth: '"$sampling_depth"'/' "${OUTPUT}/config.yaml"`
            echo -e "${conf}" > "${OUTPUT}/config.yaml"
        fi
	fi

	if [[ ${max_depth} == "null" ]] ; then
		if [ "${deblur}" == "True" ] || [ "${deblur}" == "true" ] ; then 
			max_depth=`cat stats/stats.csv | rev | cut -d "," -f3 | rev | tail -n +2 | sort -n | tail -1`
			echo -e "${YEL}WARNING: ${NC}Using Maximum Depth as max-depth; ${YEL}${max_depth}${NC}"
		elif [ "${use_QIIME2}" == "True" ] || [ "${use_QIIME2}" == "true" ] ; then 
			max_depth=`cat stats/stats.tsv | rev | cut  -f2 | rev | tail -n +3 | sort -n | tail -1`
			echo -e "${YEL}WARNING: ${NC}Using Maximum Depth as max-depth; ${YEL}${max_depth}${NC}"
		elif [ -d "../DADA2" ] ; then
    		max_depth=`cat ../DADA2/stats.csv | rev | cut -d "," -f1 | rev | tail -n +2 | sort -n | tail -1`
			echo -e "${YEL}WARNING: ${NC}Using Maximum Depth as max-depth; ${YEL}${max_depth}${NC}"
		else 
			echo -e "${RED}ERROR: ${NC}something worng! check params"
        	print_elapsed
        	exit
		fi
        if grep -Fq -- "max_depth: null" "${OUTPUT}/config.yaml"; then
            conf=`sed -r 's/^(\s*)(max_depth:\s*null\s*$)/\max_depth: '"$max_depth"'/' "${OUTPUT}/config.yaml"`
            echo -e "${conf}" > "${OUTPUT}/config.yaml"
        fi
	fi
}

download_ref(){
    if [[ $org_download_direct != true ]]; then 
        directdownload=""
    else 
        directdownload="--direct"
    fi
    if [[ $index == null && $index_name == null && $index_path == null ]] ; then
        if [[ $org_name == null ]] ; then
            if [[ $org_id == null ]]; then
                if [[ $org_txid == null ]]; then
                    if [[ $org_auto_download == false ]]; then 
                        error_cheker 1 0 "validating index and reference"
                    else 
                        :
                    fi
                else
                    mkdir -p ${working_dir}/user_downloaded_ref
                    python3 ${myfolder}/bin/WGS/scripts/NCBI_download.py --txid $org_txid -o ${working_dir}/user_downloaded_ref $directdownload
                    error_cheker $? 0 "Downloading reference"
                fi
            else
                mkdir -p ${working_dir}/user_downloaded_ref
                python3 ${myfolder}/bin/WGS/scripts/NCBI_download.py --id $org_id -o ${working_dir}/user_downloaded_ref $directdownload
                error_cheker $? 0 "Downloading reference"
            fi
        else
            mkdir -p ${working_dir}/user_downloaded_ref
            python3 ${myfolder}/bin/WGS/scripts/NCBI_download.py --name "$org_name" -o ${working_dir}/user_downloaded_ref $directdownload
            error_cheker $? 0 "Downloading reference"
        fi
    fi
    get_ref_name=`ls ${working_dir}/user_downloaded_ref`
    gunzip "${working_dir}/user_downloaded_ref/${get_ref_name}" 
    get_ref_name=${get_ref_name%".gz"}
    ref="${working_dir}/user_downloaded_ref/${get_ref_name}"
    ref_name=${get_ref_name%'.fna'}''
}


