#!/bin/bash

# colors
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'

# error check function
error_cheker(){
    lastexit=$1
    cont=$2
    name_=$3
    if [[ $(( lastexit )) -ne 0 ]];then
        echo -e "${RED}ERROR IN: ${BLE}${name_}${NC}"
        if [[ ${cont} -ne 1 ]]; then 
            echo -e "${RED}GUAP STOPPED${NC}" ; exit 1
        else
            echo -e "$? passed ${name_}"
        fi
    fi
}

# decalre help message
usage_message="
Usage: $(basename $0) -i fasta -o dir -g gtf

options:
    -i STR                  reference fasta
    -x                      create transcritome file using gtf 
    -g STR                  gtf annotation file
    -t STR                  transciptome file name to create (used with -x)
    -o PATH                 path to save the index 
    -a [aligner]            choose aligner from
                                [star|hisat2|kallisto|salmon]
    -t INT   (default=1)    number of threads
    -e STR                  extra args for aligner 
    -h                      print this message


"

# function that prints the help message and exits
print_usage(){
    echo -e "$usage_message" 1>&2; exit 1
}

# function to generate transcriptome from fa and gtf using gffread
generate_transcriptome_fa(){ # takes 1: ref fa, 2: gtf, 3: output file
    ref_fa=$1
    gtf=$2
    out=$3

    gffread $gtf -g $ref_fa -w $out
    error_cheker $? 0 "generating transcriptome from fasta"
    echo -e "${GRE}  $ref convert to trascriptome fa"

}

# main function for indexing the reference
index(){
    ref_fa=$1
    gtf=$2
    outpath=$3
    aligner=$4
    threads=$5
    extra=$6
    mkdir -p $outpath
	cd $outpath

	if [ $aligner == "star" ]
	then
		STAR --runMode genomeGenerate --runThreadN $threads --genomeDir $outpath --genomeFastaFiles $ref_fa --sjdbGTFfile $gtf $extra
	elif [ $aligner == "salmon" ]
	then
		salmon index -t $ref_fa -i $outpath $extra
	elif [ $aligner == "kallisto" ]
	then
		kallisto index -i $outpath $ref_fa $extra
	elif [ $aligner == "hisat2" ]
	then
		mkdir hisat2_index
		hisat2-build $ref_fa $outpath -p $threads $extra
	else 
		echo "${RED}ERROR!" ; exit 1
	fi
}

X=F

while getopts hxi:o:a:g:t:e: OPTION
do
	case "${OPTION}"
    	in
		i) input_path=${OPTARG}
            INPUT=$(realpath "$input_path");;
        o) output_path=${OPTARG}
            OUTPUT=$(realpath "$output_path");;
		a)
			a=${OPTARG}
                    	if [ "$a" == "hisat2" ] || [ "$a" == "star" ] || [ "$a" == "salmon" ] || [ "$a" == "kallisto" ] 
                    	then 
                        	:
                    	else
                        	echo -e "${RED}$a is not a valid -a argument${NC}"
                    	usage
                    	fi
                    	;;
		g) GTF_abs=${OPTARG}
            GTF=$(realpath "$GTF_abs");;
		h) usage;;
        t) threads=${OPTARG} ;;
		x) X="T";;
        e) extra_args=${OPTARG}
		\?) usage;;
    	esac
done

if [ "$F" == "T" ]; then 
    generate_transcriptome_fa $INPUT $GTF $new_ref
    INPUT=$new_ref
fi

index $INPUT $GTF $OUTPUT $a $threads $extra_args
error_cheker $? 0 "generating index for ${a} "