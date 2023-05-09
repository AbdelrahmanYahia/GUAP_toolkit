from ..workflows import *

class rRNA(WorkflowCli):
    name = '16s'
    help = '''rRNA 16s analysis workflow'''
    usage = "guap 16s [-i|--input PATH] [-o|--output PATH] [-m|--metadata FILE] [-c|--classifier FILE]"

    def add_arguments(self, parser):
        # basic configuration
        basic_conf = parser.add_argument_group(f'{CYN}basic config{NC}')
        basic_conf.add_argument('-i', '--input', help= "Input directory path", metavar='path',type=str)  
        basic_conf.add_argument('-o', '--output', help= "Output directory path", metavar='path',type=str)  
        basic_conf.add_argument('-t', '--threads', help= "Number of threads", type=int)

        # QC
        qc_conf = parser.add_argument_group(f'{CYN}QC{NC}')
        qc_conf.add_argument("--remove-primers", action="store_true", help="perform cutadapt to remove primers")
        qc_conf.add_argument("--min-length",  type=int,default=50, help='cutadapt (remove-primers) min length') 
        qc_conf.add_argument('-fp',"--forward-primer", help='forward primmer for cutadapt or train set ',required='--train' in sys.argv or '--remove-primers' in sys.argv)
        qc_conf.add_argument('-rp',"--reverse-primer", help='reverse primmer for cutadapt or train set', required= '--train' in sys.argv or '--remove-primers' in sys.argv)
        qc_conf.add_argument('--trimmomatic', dest='trimmomatic', action='store_true', help="Use trimmomatic [default]")
        qc_conf.add_argument("--trim-min-length",  type=int,default=30, help='trimmomatic min length') 
        qc_conf.add_argument('--skip-trimmomatic', dest='trimmomatic', action='store_false', help="Do NOT use trimmomatic")
        qc_conf.add_argument('--skip-QC', action='store_true', help="Skipp Fastqc step")

        # asv options
        asv = parser.add_argument_group(f'{CYN}ASV options{NC}')
        asv.add_argument("--deblur", action="store_true", help="Use deblur deniosing pipeline")
        asv.add_argument('-L',"--deblur-trim-length", default="-1", help='Deblur trim length, use with --deblur flag',type=int,metavar=int) 
        asv.add_argument('--use-QIIME2', action='store_true', help="Use QIIME2 to perform DADA2")
        asv.add_argument('-tf',"--trunc-f", default=0, help='DADA2 trunclength Forward', type=int) 
        asv.add_argument('-tr',"--trunc-r", default=0, help='DADA2 trunclength Reverse', type=int) 
        asv.add_argument('-l',"--trim-l", default=0, help='DADA2 trunclength Forward', type=int) 
        asv.add_argument('-r',"--trim-r", default=0, help='DADA2 trunclength Reverse', type=int) 
        asv.add_argument('-ef',"--maxee-f", default=4, help='DADA2 maxEE Forward', type=int) 
        asv.add_argument('-er',"--maxee-r", default=5, help='DADA2 maxEE Reverse', type=int)
        asv.add_argument('-mo',"--min-overlap", default=10, help='DADA2 minimum overlap to merge', type=int)
        asv.add_argument('-cm',"--chimera-method", default='consensus', help='DADA2 chimera method', type=str, choices=['consensus','pooled'])

        # classifier options 
        class_conf = parser.add_argument_group(f'{CYN}Classifier{NC}')
        class_conf.add_argument('--choose-classifier', help= "Choose to classify taxanomy [dada or qiime] defualt: qiime", default="qiime", choices=['dada', 'qiime'])
        class_conf.add_argument("-c", "--classifier", help='Classifier path', metavar='path',type=os.path.abspath)
        class_conf.add_argument('-ct', '--classifier-threads', help= "Number of threads used by classify-sklearn",default=2,type=int)

        train_setConf = parser.add_argument_group(f'{CYN}Train classifier{NC}')
        train_setConf.add_argument("--train", action="store_true",help="Train set for qiime2 naive classifier")
        train_setConf.add_argument("--t-i-seqs", help='Train set input sequences (qza)',required='--train' in sys.argv,type=os.path.abspath)
        train_setConf.add_argument("--t-i-taxa", help='Train set input taxa (qza)',required='--train' in sys.argv,type=os.path.abspath)
        train_setConf.add_argument("--trainset-minlength", help='Train set input sequences (qza)',required='--train' in sys.argv,type=os.path.abspath)
        train_setConf.add_argument("--trainset-maxlength", help='Train set input sequences (qza)',required='--train' in sys.argv,type=os.path.abspath)


        # downstream options
        downstream_conf = parser.add_argument_group(f'{CYN}Downstream{NC}')
        downstream_conf.add_argument("--downstream", action="store_true", help="Perform Downstream analysis")
        downstream_conf.add_argument('-m', '--metadata', help= "metadata file full path",required='--downstream' in sys.argv, metavar='path',type=os.path.abspath)  
        downstream_conf.add_argument("--bashdownstream", action="store_true",help="berform downstream analysis on analysis run by --bash argument (downstream only)")
        downstream_conf.add_argument('-d',"--sampling-depth",  type=int, help='Sampling Depth for core-mertic phylogenetic') 
        downstream_conf.add_argument('-md',"--max-depth",  type=int, help='Maximum Sampling Depth for alpha-rarefaction') 
        downstream_conf.add_argument( '--condition-name', help= "condition name in metadata file", required='--export-figs' in sys.argv, type=str)
        downstream_conf.add_argument('--export-figs', action='store_true', help="export phyloseq figures (alpha, beta and bar plots)")
        downstream_conf.add_argument('--export-figs-only', action='store_true', help="export phyloseq figures (alpha, beta and bar plots)")


        # performance
        snakemake_options = parser.add_argument_group(f'{CYN}snakemake{NC}')
        snakemake_options.add_argument('--bash', dest='snakemake', action='store_false', help="Use bash scripts instead of snakemake")
        snakemake_options.add_argument('--snakemake', dest='snakemake', action='store_true', help="Use snakemake workflow (ability to contiue) [default]")
        snakemake_options.add_argument('--snakemake-dry-run', action='store_true', help="performs snakemake dry run")
        snakemake_options.add_argument('--snakemake-dag', action='store_true', help="performs snakemake dry run and exports DAG")


        # other options
        other_conf = parser.add_argument_group(f'{CYN}Other{NC}')
        other_conf.add_argument('--continue', action='store_true', help="continue analysis when re-run")
        other_conf.add_argument('-n',"--name", default="/GUAP-16s-dada2", help='Name of files') 
        other_conf.add_argument('--version', action='version', help="Print version number")
        other_conf.add_argument('--verbose', action='store_true', help="verbose")
        other_conf.add_argument('--quit', dest='verbose', action='store_false', help="print many output")

        parser.set_defaults(snakemake=True)
        parser.set_defaults(trimmomatic=True)
        parser.set_defaults(verbose=False)

    def run(self, args):
        print(vars(args))


class WES(WorkflowCli):
    name = 'WES'
    help = '''guap WES -i/--input dir'''
    usage = ""

    def add_arguments(self, parser):
        # basic inputs
        parser.add_argument('-i', '--input', help= "Input directory path", metavar='path', type=os.path.abspath)  
        parser.add_argument('-o', '--output', help= "Output directory path", metavar='path', type=os.path.abspath)
        # performance
        parser.add_argument('-t', '--threads', help= "Number of total no. of threads", type=int)
        parser.add_argument('--threads-index', help= "Number of threads to use during indexing ref", type=int )
        parser.add_argument('--threads-align', help= "Number of threads to use per sample aliging", type=int, default = 4)
        parser.add_argument('--threads-calling', help= "Number of threads to use per sample variant caller", type=int, default = 4)
        # main process
        parser.add_argument('--aligner', help = "Choose aligner", choices=["bwa", "bowtie2"], type=str, default='bwa')
        parser.add_argument('--variant-caller', help = "Choose variant caller", choices=["GATK", "mpileup", "lofreq"], type=str, default='GATK')
        parser.add_argument('--bed-file', help='bed file path', metavar='path', type=os.path.abspath)
        parser.add_argument('--reference-fasta',metavar='path', type=os.path.abspath, help="path to reference fasta file")
        parser.add_argument('--reference-index',metavar='path', type=os.path.abspath, help="path to reference index")
        parser.add_argument('--known-variants',metavar='path', type=os.path.abspath, help="path to reference fasta file")
        # QC
        parser.add_argument('--skip-QC', action='store_true', help="Skip QC step")
        parser.add_argument('--index-reference', action='store_true', help="index reference")
        # snakemake options
        parser.add_argument('--snakemake-dry-run', action='store_true', help="performs snakemake dry run")
        parser.add_argument('--snakemake-dag', action='store_true', help="performs snakemake dry run and exports DAG")
        parser.add_argument('--snakemake', dest='snakemake', action='store_true', help="Use snakemake workflow (ability to contiue) [default]")
        parser.add_argument('--bash', dest='snakemake', action='store_false', help="Use bash scripts instead of snakemake")
        # other utils
        parser.add_argument('--continue', action='store_true', help="continue analysis when re-run")
        parser.add_argument('--overwrite', action='store_true', help="overwrite output dir if exsits")
        parser.add_argument('--verbose', action='store_true', help="verbose")
        parser.add_argument('--quit', dest='verbose', action='store_false', help="print many output")

        # defaults
        parser.set_defaults(verbose=False)
        parser.set_defaults(snakemake=True)


    def run(self, args):
        print(vars(args))

