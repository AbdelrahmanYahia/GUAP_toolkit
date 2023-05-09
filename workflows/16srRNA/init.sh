help_message="GUAP 16s analyis pipeline:

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
"

if [[ $2 == "-h" || $2 == "--help" ]];then
    echo -e "${help_message}"
else
    echo -e "${YEL}Activating GUAP conda ENV...${NC}"
    eval "$(conda shell.bash hook)"
    conda deactivate
    conda activate GUAP
    # parse rest of arguments and create samples sheet 
    python3 ${myfolder}/workflows/16s/16s.py "${@:2}"
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
Workflow           :${CYN} 16s${NC}
snakemake          :${CYN} ${snakemake}${NC}
deblur             :${CYN} ${deblur}${NC}
qiime              :${CYN} ${use_QIIME2}${NC}
classifier         :${CYN} ${choose_classifier}${NC}
export figures     :${CYN} ${export_figs}${NC}
trimmomatic        :${CYN} ${trimmomatic}${NC}
Total threads      :${CYN} ${threads}${NC}
================================"
    if [[ $snakemake == "true" && $bashdownstream != "true" ]]; then
        cd ${myfolder}/workflows/16s/snakemake
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
    elif [ $bashdownstream == "true" ]; then 
    source ${myfolder}/workflows/16s/bash/Main
    echo -e "${BLU}NOTE: Using Bash to perform Downstream only...${NC}"
    donwstreamonly
    error_cheker $? 0 "performing 16s downstream"
    
    elif [[ "${export_figs_only}" == "true" ]] ; then 
    echo -e "${BLU}NOTE: Generating Phyloseq figures only...${NC}"
    source ${myfolder}/workflows/16s/bash/Main
    INPUT=${path}
        OUTPUT=${working_dir}
    cd ${OUTPUT}
    step_checker export_ps_figs
    else
    source ${myfolder}/workflows/16s/bash/Main
    echo -e "${BLU}NOTE: Using Bash to perform the ananlysis...${NC}"
    if [ $downstream == "true" ] ; then 
        run_16s
        error_cheker $? 0 "running 16s workflow"
        performe_downstream
        error_cheker $? 0 "performing 16s downstream"
        step_checker export_ps_figs
    else
        run_16s
        error_cheker $? 0 "running 16s workflow"
    fi
    fi
fi
