#!/bin/bash
# isntall GUAP source code in home

YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'
getdir=${PWD}
# Create destination folder
DESTINATION="${HOME}/"
mkdir -p ${DESTINATION}

error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}

echo -e "${YEL}Unpacking package...${NC}"
# Find __ARCHIVE__ maker, read archive content and decompress it
ARCHIVE=$(awk '/^__ARCHIVE__/ {print NR + 1; exit 0; }' "${0}")
tail -n+${ARCHIVE} "${0}" | tar xpJ -C ${DESTINATION}
error_cheker $?
echo -e "${YEL}configuring package...${NC}"
cd ~/GUAP
error_cheker $?
chmod 777 GUAP
error_cheker $?
echo -e "${YEL}Updating environmental variables...${NC}"
echo '' >> ~/.bashrc
echo 'GUAP(){ ' >> ~/.bashrc
echo '~/GUAP/GUAP "$@"' >> ~/.bashrc
echo '}' >> ~/.bashrc
error_cheker $? 
echo 'export -f GUAP' >> ~/.bashrc 
error_cheker $?
cd install
echo -e "${YEL}Installing GUAP ENV${NC}"

conda env create -n GUAP --file GUAP.yml --quiet >> log.txt
error_cheker $?

echo -e "${YEL}Activating GUAP ENV${NC}"
eval "$(conda shell.bash hook)"
conda activate GUAP
error_cheker $?

echo -e "${YEL}Installing R packages${NC}"
touch log.txt
Rscript install_packages.R >> log.txt
error_cheker $?

echo -e "${YEL}Restart terminal and run 'GUAP [16s/DE] -h'${NC}"
echo -e "${GRE}ALL DONE...${NC}"

exit 0

























__ARCHIVE__
