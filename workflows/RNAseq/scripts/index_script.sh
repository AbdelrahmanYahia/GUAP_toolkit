#!/bin/bash

# colors
YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'


__usage="
Usage: $(basename $0) [OPTIONS]
Options:
	-i <str>					reference genome Input file
	-o <str>					Output file
	-a <hisat2/star/salmon/kallsito>		Choose an aligner
	-g <str>					GTF file 
	-x 							create transcripts fa file from genome and use as reference		
	-O							transcriptome out file			
	-h							Help message
"


usage() {
	echo -e "$__usage" 1>&2; exit 1
}


generate_transcriptome_fa(){
	reference=$1
	gtf=$2
	out=$3

	gffread $gtf -g $reference -w $output

}

Index(){
	reference_dir=$1
	output_dir=$3
	aligner=$4
	gtf_dir=$5
	tool_name=$6
	tool_path=$7
	cd $output_dir
	if [ $aligner == "star" ]
	then
		STAR --runMode genomeGenerate --genomeDir star_index/ --genomeFastaFiles $reference_dir --sjdbGTFfile $gtf_dir --genomeSAindexNbases 11
	elif [ $aligner == "salmon" ]
	then
		salmon index -t $reference_TRA_dir -i salmon_index -k 31
	elif [ $aligner == "kallisto" ]
	then
		kallisto index -i kaliisto_index $reference_TRA_dir
	elif [ $aligner == "hisat2" ]
	then
		mkdir hisat2_index
		hisat2-build $reference_dir hisat2_index/hisat2
	else
		echo "${RED}ERROR!" ; exit 1
	fi
}


while getopts hxi:o:a:g:r:p: OPTION
do
	case "${OPTION}"
    	in
		i) INPUT=${OPTARG};;
        o) OUTPUT=${OPTARG};;
		a)
			a=${OPTARG}
                    	if [ "$a" == "hisat2" ] || [ "$a" == "star" ] || [ "$a" == "salmon" ] || [ "$a" == "kallisto" ] || [ "$a" == "RSEM" ]
                    	then 
                        	:
                    	else
                        	echo -e "${RED}$a is not a valid -a argument${NC}"
                    	usage
                    	fi
                    	;;
		g) GTF=${OPTARG};;
		O) new_ref=${OPTARG};;
		h) usage;;
		x) 
			generate_transcriptome_fa $INPUT $GTF $new_ref
			INPUT=$new_ref
			;;
		\?) usage;;
    	esac
done

Index $INPUT $OUTPUT $a $GTF $toolname $toolpath
