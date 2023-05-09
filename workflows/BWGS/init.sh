
help_message_WGS="GUAP WGS analyis pipeline:

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
"

if [[ $2 == "-h" || $2 == "--help" ]];then
    echo -e "${help_message_WGS}"
else
    echo -e "${YEL}Activating GUAP conda ENV...${NC}"
    eval "$(conda shell.bash hook)"
    conda deactivate
    conda activate GUAP
    # parse rest of arguments and create samples sheet 
    python3 ${myfolder}/workflows/WGS/WGS.py "${@:2}"
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
Workflow           :${CYN} WGS metagenomics${NC}
snakemake          :${CYN} ${snakemake}${NC}
skip QC            :${CYN} ${skip_QC}${NC}
skip kraken        :${CYN} ${skip_check_contaminant}${NC}
skip-assembly      :${CYN} ${skip_assembly}${NC}
total threads      :${CYN} ${threads}${NC}
================================"
    source ${myfolder}/bin/WGS/bash/Main.sh
    	max_n_samples=$(( $threads / 4 ))
	if [[ $max_n_samples -eq 0 ]] ; then max_n_samples=1 ; fi
	INPUT=${path}
	OUTPUT=${working_dir}
	touch "${OUTPUT}/.checker_array.txt"
	
	R1s=$( tail -n +2 ${OUTPUT}/samples.tsv | while read -r uniquename id restofname extension R1 R2 tail PE path R; do echo "$R1"; done ) 
	R2s=$( tail -n +2 ${OUTPUT}/samples.tsv | while read -r uniquename id restofname extension R1 R2 tail PE path R; do echo "$R2"; done ) 
	
    download_ref
    error_cheker $? 0 "Downloading ref"


    echo -e "${BLU}NOTE: Using Bash to perform the analysis...${NC}"
    run_WGS
    error_cheker $? 0 "running WGS workflow"
fi