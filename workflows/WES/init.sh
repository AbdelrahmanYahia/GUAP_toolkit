
help_message_WES="GUAP WES analyis pipeline:

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
"


if [[ $2 == "-h" || $2 == "--help" ]];then
    echo -e "${help_message_WES}"
else
    echo -e "${YEL}Activating GUAP conda ENV...${NC}"
    eval "$(conda shell.bash hook)"
    conda deactivate
    conda activate GUAP
    # parse rest of arguments and create samples sheet 
    python3 ${myfolder}/workflows/WES/WES.py "${@:2}"
    error_cheker $? 0 "validating input"
    # parse user options stored in config.yaml current env
    eval $(parse_yaml ${current}/config.yaml)
    error_cheker $? 0 "parsing yaml"
    mv ${current}/config.yaml ${working_dir}
    error_cheker $? 0 "moving yaml file to workdir"
    # store command used at working dir
    echo $@ > "${working_dir}/command.txt"
    echo -e "
================================
${CYN}GUAP settings:${NC}
================================
Workflow           :${CYN} WES${NC}
snakemake          :${CYN} ${snakemake}${NC}
aligner            :${CYN} ${aligner}${NC}
variant caller     :${CYN} ${variant_caller}${NC}
total threads      :${CYN} ${threads}${NC}
threads to index   :${CYN} ${threads_index}${NC}
threads to align   :${CYN} ${threads_align}${NC}
threads to calling :${CYN} ${threads_calling}${NC}
================================"
    if [[ $snakemake == "true" ]]; then
        cd ${myfolder}/workflows/WES/snakemake
        if [[ $verbose == "true" ]]; then
            qtag=""
        else
            qtag="-q"
        fi
        if [[ $snakemake_dry_run == "true" ]]; then 
            echo -e "${BLU}NOTE: Using snakemake workflow to perform dry-run...${NC}"
            if [[ $snakemake_dag == "true" ]]; then
                snakemake --configfile "${working_dir}/config.yaml" -n --rulegraph | dot -Tpng > "${working_dir}/${name}.png"
            else
                snakemake -j ${threads} --configfile "${working_dir}/config.yaml" -n ${qtag}
            fi
        elif [[ $snakemake_dag == "true" ]]; then
            echo -e "${BLU}NOTE: Using snakemake workflow to Generate DAG...${NC}"
            snakemake --configfile "${working_dir}/config.yaml" -n --rulegraph | dot -Tpng > "${working_dir}/${name}.png"
        else
            echo -e "${BLU}NOTE: Using snakemake workflow to perform the ananlysis...${NC}"
            snakemake -j ${threads} --configfile "${working_dir}/config.yaml" ${qtag}
        fi 
    elif [ $bash == "true" ]; then 
        ehco -e "${RED}BASH ${CYN}is ${YEL} currently ${GRN}under dev.${NC}"
        # source ${myfolder}/workflows/WGS/bash/Main.sh
    	# max_n_samples=$(( $threads / 4 ))
        # if [[ $max_n_samples -eq 0 ]] ; then max_n_samples=1 ; fi
        # INPUT=${path}
        # OUTPUT=${working_dir}
        # touch "${OUTPUT}/.checker_array.txt"
        
        # R1s=$( tail -n +2 ${OUTPUT}/samples.tsv | while read -r uniquename id restofname extension R1 R2 tail PE path R; do echo "$R1"; done ) 
        # R2s=$( tail -n +2 ${OUTPUT}/samples.tsv | while read -r uniquename id restofname extension R1 R2 tail PE path R; do echo "$R2"; done ) 
        
        # error_cheker $? 0 "running WES workflow"
    else
        error_cheker 1 0 "WES init"
    fi
fi
