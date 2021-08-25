#!/bin/bash
source indexes.sh

# progress bar ## takes three arguments 
# 1. total 2. current 3. message to print 
progress_bar() {
    BAR="#"
    MIN=1
    MAX=100
    TIMES=$1
    COUNT=$2
    MESSAGETOSHOW=$3
    PROGRESS=$(echo | awk -v t=$TIMES -v c=$COUNT -v max=$MAX '{ print int(c / t * max) }')
    if [ $(( PROGRESS )) -eq 0 ];then
        PROGRESS=1
    fi
    PROGRESS_BAR=""
    for _i in $(seq $MIN $PROGRESS); do PROGRESS_BAR="${PROGRESS_BAR}${BAR}"; done
    
    printf "\r[%-100s] %3d%% - %s " $PROGRESS_BAR $PROGRESS "$MESSAGETOSHOW"
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

# checks sample extension and paired or not 
sample_check(){
    dogz=$2 # if true will gunzip unzipped samlpes 
    counter=1
    return_value=0 # exit code 
    no_of_files=`ls ${1}/ | wc -l ` # gets number of files 
    echo -e "${YEL}NUMBER OF FILES = ${no_of_files}${NC}" 
    if [ $((no_of_files%2)) -ne 0 ] # checks if samples are not pairs 
    then 
        echo -e "${RED}SAMPLES NOT IN PAIRS${NC} ${YEL}CHOOSING SE PIPELINE${NC}"
        if [ `ls ${1}/*R1* 2> /dev/null | wc -l ` -ne 0 ];then
            for R1 in ${1}/*R1*
            do
                R2=$(echo $R1 | sed 's/R1/R2/')
                    if [ -f $R2 ]; then # if finds R2 returns erorr
                        echo -e "${RED}ERORR WITH SAMPLE NAMES (EC1)${NC}"
                        return_value=0
                        break
                    else
                        return_value=1
                    fi
            done
        else # if samples are paires 
            if [ `ls ${1}/*R2* 2> /dev/null | wc -l ` -ne 0 ];then 
                return_value=0
                echo -e "${RED}ERORR WITH SAMPLES NAMES (EC2)${NC}"
            else
                return_value=1
            fi
        fi
    else
        if [ `ls ${1}/*R1* 2> /dev/null | wc -l ` -ne 0 ] # checks presence of R2 
        then
            for R1 in ${1}/*R1*
            do
                R2=$(echo $R1 | sed 's/R1/R2/')
                if [ ! -f $R2 ]; then
                    return_value=1
                    break
                else
                    return_value=2
                fi
            done
        else
            echo -e "${RED}(EC3)${YEL} SAMPLES SEEMS TO BE SE ${NC}"
            return_value=1
        fi
    fi
    # checks files extensions 
    for i in ${1}/*
    do
        if [ "${i: -9}" == '.fastq.gz' ] # files should be with extension fastq.qz 
        then 
            echo -e "${i} ${GRE}PASSED${NC}"
            :
        elif [ ${i: -6} == ".fq.gz" ]
        then
            newname=`echo $i | sed 's/.fq.gz/.fastq.gz/'` # replaces fq.qz to fastq.qz
            mv $i $newname
            echo -e "${i} ${YEL}MODEFIED ${RED}(EC4)${NC}"
            :
        elif [ ${i: -6} == ".fastq" ] # if not compressed and analysis need compression it will compress file 
        then
            if [[ $dogz == 'True' ]]
            then
                gzip $i
                echo -e "${i} ${YEL}MODEFIED  ${RED}(EC5)${NC}"
                :
            else
                echo -e "${i} ${GRE}PASSED${NC}"
            fi
        elif [ ${i: -3} == ".fq" ]
        then
            newname=`echo $i | sed 's/.fq/.fastq/'`
            mv ${i} $newname
            if [ $dogz == 'True' ]
            then
                gzip ${i%'.fq'}'.fastq'
                echo -e "${i} ${YEL}MODEFIED  ${RED}(EC6)${NC}"
            else
                echo -e "${i} ${GRE}PASSED${NC}"
            fi
            :
        else
            echo -e "${RED}${i} DID NOT PASS (EC7)${NC}"
        fi
        counter=$(( counter + 1 ))
    done
    # checks for files with extension not fastq or fastq.gz 
    if [ `ls ${1}/ -I *'.fast'* | wc -w ` -ne 0 ];then 
        echo -e "${RED}FOUND FILES WITH EXTENSION NOT fastq.gz (EC 7)${NC}"
        return_value=0
    fi
    return $return_value
    
}

# checks sample extension and paired or not 
dir_check(){
    check_extension=$1 # if true will gunzip unzipped samlpes 
    counter=1
    return_value=0 # exit code 
    no_of_files=`ls | wc -l ` # gets number of files 
    if [ $((no_of_files)) -eq 0 ] # checks if samples are not pairs 
    then 
        return_value=0
        echo -e "${RED}ERROR${NC}"
        exit
    else
        if [ `ls *"${check_extension}" 2> /dev/null | wc -l ` -eq 0 ];then 
                return_value=0
                echo -e "${RED}ERROR${NC}"
                exit
        else
            return_value=1
        fi
    fi
    return $return_value 
}

error_cheker(){
    lastexit=$1
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR in exit value of last process${NC}"
        exit
    fi
}

read_config (){
    while read -r name value
    do
        echo "Content of $name is ${value//\"/}"
    done < $1
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
