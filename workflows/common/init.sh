echo -e "${YEL}Activating GUAP conda ENV...${NC}"
eval "$(conda shell.bash hook)"
conda deactivate
conda activate GUAP
python3 ${myfolder}/bin/modules/modules.py "${@:2}"
error_cheker $?
if [[ $2 == "-h" || $2 == "--help" ]];then
    exit
else
    eval $(parse_yaml ${current}/config.yaml)
    error_cheker $?
    mv ${current}/config.yaml ${working_dir}
    error_cheker $?
    echo $@ > "${working_dir}/command.txt"
    
    if [ "${GUAP_args}" != "null" ]; then 
    source ${myfolder}/bin/modules/bash/guap_module.sh
    run_fx_pe
    
    elif [[ ${module} == "qc" || ${module} == "trimmomatic" || ${module} == "cutadapt" ]] ; then 
    cd ${myfolder}/bin/modules/snakemake
    snakemake -j ${threads} --configfile "${working_dir}/config.yaml" -q
    error_cheker $?
    else
    echo -e "${RED}Currently ${module} is not supported${NC}"
    
    fi
fi