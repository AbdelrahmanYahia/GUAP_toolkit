import sys
import argparse
import subprocess
from ..utils.globals import *
from ..utils.parse_input import parse_input_args


def smk_cmd(all_args):
    smk_cmd = ""

    if all_args['dry_run'] is True:
        smk_cmd += " -n"
        glogger.prnt_warning("performing dry run")
        if all_args['export_dag'] is True:
            smk_cmd += f" --rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'"
            glogger.prnt_warning("exporting dag")
        else:
            pass
    elif all_args['export_dag'] is True and all_args['dry_run'] != True:
        smk_cmd += f"--rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'"

    else:
        smk_cmd = smk_cmd + f"{all_args['smk_extra_args']}"


    snakemake_cmd = f"snakemake --snakefile '{GUAP_DIR}/workflows/RNAseq/Snakefile' --configfile '{all_args['working_dir']}/config.yaml' -j {all_args['threads']} {smk_cmd}"
    return snakemake_cmd

# Create a custom action to preserve quoted strings
class StoreQuotedString(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # If the value is surrounded by quotes, remove them
        if isinstance(values, str) and values.startswith('"') and values.endswith('"'):
            values = values[1:-1]
        setattr(namespace, self.dest, values)

# base class for subcommands (workflows)
class WorkflowCli:
    def __init__(self,subparser):
        self.parser = subparser.add_parser(self.name, help=self.help, usage=self.usage)
        self.add_arguments(self.parser)

    def add_arguments(self, parser):
        pass
    
    def run(self, args):
        pass

# main RNA subparser class
class RNA(WorkflowCli):
    name = 'RNA'
    help = '''guap RNA -i/--input dir'''
    usage = f"""{YEL}Basic Run Usage example:{NC}
guap RNA --input dir \\
        --output dir \\
        --aligner [star|hisat2|kallisto|salmon] \\
        --quantifier [htseq|featurecounts] \\
        --gtf-file file \\
        --reference-fasta file \\
        --reference-index dir
        """

    def add_arguments(self, parser):

        # basic configuration
        basic_conf = parser.add_argument_group(f'{CYN}basic config{NC}')

        basic_conf.add_argument(
            '-i', '--input', 
            help="Input directory path", 
            metavar='in path', 
            type=os.path.abspath, 
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)
        ) 

        basic_conf.add_argument(
            '-o', '--output', 
            help= "Output directory path", 
            metavar='out path', 
            type=os.path.abspath, 
            required = not any(arg in ["--print-last-run"] for arg in sys.argv)
        )  

        basic_conf.add_argument(
            '-L', '--last-run-params', 
            help= "Uses Last Run commands, Except -i and -o", 
            action='store_true'
        ) 

        # workflow configure
        workflow_conf = parser.add_argument_group(f'{CYN}Workflow configure{NC}')

        workflow_conf.add_argument(
            '-t', '--threads', 
            metavar = "N",
            help= "Number of total threads to use [default = all]", 
            type=int,

        )

        workflow_conf.add_argument(
            '--reference-fasta',
            metavar='path/to/file.fa',
            type=os.path.abspath,
            help="path to reference fasta file"
        )

        workflow_conf.add_argument(
            '--gtf-file',
            help='gtf file path',
            metavar='path/to/file.gtf',
            type=os.path.abspath
        )

        qc_conf = parser.add_argument_group(f'{BLU}QC configuration{NC}')

        qc_conf.add_argument(
            '--trimmomatic',
            dest='trimmomatic',
            action='store_true',
            help="Use trimmomatic"
        )

        qc_conf.add_argument(
            '--trim-t', 
            help= "Number of threads to use during trim step", 
            type=int ,
            metavar = "N",
            default= 4 
        )

        qc_conf.add_argument(
            "--trim-min-length", 
            type=int,
            default=30,
            metavar = "N",
            help='trimmomatic min length [default = 30]'
        )

        qc_conf.add_argument(
            "--slidingwindow-size", 
            type=int,
            default=4,
            metavar = "N",
            help='trimmomatic sliding window size [default = 4]'
        )

        qc_conf.add_argument(
            "--slidingwindow-quality", 
            type=int,
            default=10,
            metavar = "N",
            help='trimmomatic sliding window quality score [default = 10]'
        )

        qc_conf.add_argument(
            '--trimmomatic-extra-args',
             type=str,
            metavar="='-args'",
            help="A string value of extra args for trimmomatic (must be used with = with no spaces (--trimmomatic-extra-args='-arg1 -arg2'))",
            default=""
        )


        qc_conf.add_argument(
            '--skip-QC',
             action='store_true',
             help="Skipp Fastqc step"
        )

        qc_conf.set_defaults(trimmomatic=False)

        aligner_conf = parser.add_argument_group(f'{BLU}Aligner configuration{NC}')

        aligner_conf.add_argument(
            '--aligner', 
            help = "Choose aligner [default = star]",
            choices=["salmon", "kallisto", "hisat2", "star"], 
            metavar = "salmon|kallisto|hisat2|star",
            type=str
        )

        aligner_conf.add_argument(
            '--align-t', 
            help= "Number of threads to use during indexing ref [default = 4]", 
            type=int, 
            default= 4,
            metavar = "N",
        )

        aligner_conf.add_argument(
            '--aligner-extra-args', 
            help = "Extra arguments for aligner, use it with no spaces and add = ( --aligner-extra-args='-arg1 -arg2' ) ",
            type=str,
            metavar = "'-args'",
            default=""
        )

        aligner_conf.add_argument(
            '--reference-index',
            metavar='path/to/ref', 
            type=os.path.abspath,
            help="path to reference index"
        )

        quantifier_conf = parser.add_argument_group(f'{BLU}quantifier configuration{NC}')

        quantifier_conf.add_argument(
            '--quantifier',
            help = "Choose quantifier, if empty will skip quantification",
            choices=["htseq", "featurecounts"], 
            type=str,
            metavar = 'htseq|featurecounts'
        )

        quantifier_conf.add_argument(
            '--quant-t', 
            help= "Number of threads to use during quantification [default = 4]", 
            type=int , 
            default= 4
        )

        quantifier_conf.add_argument(
            '--feature', 
            help="feature to use in quantifing [default = gene_name]", 
            default= "gene_name"
        )

        quantifier_conf.add_argument(
            '--quantifier-extra-args', 
            metavar="='-args'",
            help="A string value of extra args for quatnifier (must be used with = with no spaces (--quantifier-extra-args='-arg1 -arg2'))",
            default="", type=str
        )

        # Snakemake Options
        snakemake_options = parser.add_argument_group(f'{CYN}Snakemake Options{NC}')

        snakemake_options.add_argument(
            '--dry-run', 
            action='store_true', 
            help="performs snakemake dry run"
        )

        snakemake_options.add_argument(
            '--export-dag', 
            action='store_true', 
            help="performs snakemake dry run and exports DAG"
        )

        snakemake_options.add_argument(
            "--smk-extra-args", 
            metavar="='-args'",
            help="A string value of extra args for snakemake(must be used with = with no spaces (--smk-extra-args='-arg1 -arg2'))",
            default="", type=str
        )

        # DE config
        de_config = parser.add_argument_group(f'{CYN}Differential expression{NC}')

        de_config.add_argument(
            '--perform-DE', 
            action='store_true', 
            help="perform Differential expression analysis"
        )

        de_config.add_argument(
            "--metadata", 
            help="Samples clinical/metadata file", 
            metavar='metadata file', 
            type=os.path.abspath, 
            required='--perform-DE' in sys.argv
        )

        de_config.add_argument(
            "--control-name", 
            help="Control name in clinical/metadata file", 
            metavar='str', 
            type=str,
            required='--perform-DE' in sys.argv
        )

        de_config.add_argument(
            "--sample-name", 
            help="Samples name in clinical/metadata file", 
            metavar='str', 
            type=str,
            required='--perform-DE' in sys.argv
        )

        de_config.add_argument(
            "--dirstr", 
            help="directory string to rempove from column names in counts file", 
            metavar='str', 
            type=str,
            default=""
        )

        de_config.add_argument(
            "--ID-extension-ptrn", 
            help="RE pattern to rempove from column names in counts file", 
            metavar='str', 
            type=str,
            default=""
        )

        de_config.add_argument(
            "--sample-zeros", 
            help="Number zeros allowed per gene in all samples [default = 1]", 
            metavar='N', 
            type=str,
            default=1
        )

        de_config.add_argument(
            "--control-zeros", 
            help="Number zeros allowed per gene in all controls [default = 1]", 
            metavar='N', 
            type=str,
            default = 1
        )
        # other options
        other_conf = parser.add_argument_group(f'{CYN}Other{NC}')

        other_conf.add_argument(
            '--continue', 
            action='store_true', 
            help="continue analysis when re-run"
        )

        other_conf.add_argument(
            '--overwrite', 
            action='store_true', 
            help="overwrite output dir if exsits"
        )

        other_conf.add_argument(
            '-n',"--name", 
            default=f"RNA_run[{os.environ['start_time']}]", 
            metavar = 'str',
            help=f"Name of files [ default = RNA_run[date time] ]"
        )

        other_conf.add_argument(
            '--verbose', 
            action='store_true', 
            help="verbose"
        )

        other_conf.add_argument(
            '--quit', 
            dest='verbose', 
            action='store_false', 
            help="print many output"
        )

        other_conf.add_argument(
            '--print-last-run', 
            action='store_true', 
            help="Prints last run on screen"
        )

        other_conf.set_defaults(verbose=False)


    def run(self, args):
        if args.print_last_run:
            with open(f"{GUAP_DIR}/.last_run.txt", 'r') as last_run:
                lines = last_run.readlines()
            last_command = lines[0]
            print(f"python3 {GUAP_DIR}guap.py {last_command}")
            exit()
        all_args = parse_input_args(args)
        snakemake_cmd = smk_cmd(all_args)
        try:
            if all_args['export_dag'] is True and all_args['dry_run'] != True:
                if all_args['continue']:
                    subprocess.run(f"snakemake --snakefile '{GUAP_DIR}/workflows/RNAseq/Snakefile' --configfile '{all_args['working_dir']}/config.yaml' -j {all_args['threads']} {all_args['smk_extra_args']}", shell=True)
                else:
                    subprocess.run(snakemake_cmd, shell=True)
                    subprocess.run(f"snakemake --snakefile '{GUAP_DIR}/workflows/RNAseq/Snakefile' --configfile '{all_args['working_dir']}/config.yaml' -j {all_args['threads']} {all_args['smk_extra_args']}", shell=True)

                print(f"{PRP}{runtime.elapsed()}{NC}")
            else:
                subprocess.run(f"{snakemake_cmd} -q -n --rulegraph | dot -Tpng > '{all_args['working_dir']}/{all_args['name']}.png'", shell=True)
                subprocess.run(snakemake_cmd, shell=True)
                print(f"{PRP}{runtime.elapsed()}{NC}")
        except:
            glogger.prnt_fatel("Error in snakemake dry run")
            print(f"{PRP}{runtime.elapsed()}{NC}") 



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


class RNA(WorkflowCli):
    name = 'RNA'
    help = '''guap RNA -i/--input dir'''
    usage = f"""{YEL}Basic Run Usage example:{NC}
guap RNA --input dir \\ 
        --output dir \\
        --aligner [star|hisat2|kallisto|salmon] \\
        --quantifier [htseq|featurecounts] \\
        --gtf-file file \\
        --reference-fasta file \\
        --reference-index dir \\
        """

    def add_arguments(self, parser):

        # basic configuration
        basic_conf = parser.add_argument_group(f'{CYN}basic config{NC}')

        basic_conf.add_argument(
            '-i', '--input', 
            help="Input directory path", 
            metavar='path', 
            type=os.path.abspath, 
            required = True
        ) 

        basic_conf.add_argument(
            '-o', '--output', 
            help= "Output directory path", 
            metavar='path', 
            type=os.path.abspath, 
            required = True
        )  

        # workflow configure
        workflow_conf = parser.add_argument_group(f'{CYN}Workflow configure{NC}')

        workflow_conf.add_argument(
            '-t', '--threads', 
            help= "Number of total threads to use", 
            type=int
        )

        workflow_conf.add_argument(
            '--reference-fasta',
            metavar='path',
            type=os.path.abspath,
            help="path to reference fasta file"
        )

        workflow_conf.add_argument(
            '--gtf-file',
            help='gtf file path',
            metavar='path',
            type=os.path.abspath
        )

        qc_conf = parser.add_argument_group(f'{BLU}QC configuration{NC}')

        qc_conf.add_argument(
            '--trimmomatic',
            dest='trimmomatic',
            action='store_true',
            help="Use trimmomatic"
        )

        qc_conf.add_argument(
            '--trim-t', 
            help= "Number of threads to use during trim step", 
            type=int , 
            default= 4 
        )

        qc_conf.add_argument(
            "--trim-min-length", 
            type=int,
            default=30,
            help='trimmomatic min length'
        )

        qc_conf.add_argument(
            "--slidingwindow-size", 
             type=int,
            default=4,
             help='trimmomatic sliding window size'
        )

        qc_conf.add_argument(
            "--slidingwindow-quality", 
            type=int,
            default=10,
            help='trimmomatic sliding window quality score'
        )

        qc_conf.add_argument(
            '--trimmomatic-extra-args',
             type=str,
             help="Extra arguments for aligner, don't forget to use ('')"
        )


        qc_conf.add_argument(
            '--skip-QC',
             action='store_true',
             help="Skipp Fastqc step"
        )

        qc_conf.set_defaults(trimmomatic=False)

        aligner_conf = parser.add_argument_group(f'{BLU}Aligner configuration{NC}')

        aligner_conf.add_argument(
            '--aligner', 
            help = "Choose aligner",
            choices=["salmon", "kallisto", "hisat2", "star"], 
            type=str
        )

        aligner_conf.add_argument(
            '--align-t', 
            help= "Number of threads to use during indexing ref", 
            type=int, 
            default= 4 
        )

        aligner_conf.add_argument(
            '--aligner-extra-args', 
            help = "Extra arguments for aligner, don't forget to use ('') ",
            type=str
        )

        aligner_conf.add_argument(
            '--reference-index',
            metavar='path', 
            type=os.path.abspath,
            help="path to reference index"
        )

        quantifier_conf = parser.add_argument_group(f'{BLU}quantifier configuration{NC}')

        quantifier_conf.add_argument(
            '--quantifier',
            help = "Choose quantifier",
            choices=["htseq", "featurecounts"], 
            type=str
        )

        quantifier_conf.add_argument(
            '--quant-t', 
            help= "Number of threads to use during indexing ref", 
            type=int , 
            default= 2
        )

        quantifier_conf.add_argument(
            '--feature', 
            help="feature to use in quantifing", 
            default= "gene_name"
        )

        quantifier_conf.add_argument(
            '--quantifier-extra-args', 
            help = "Extra arguments for quantifier, don't forget to use ('') ",
            type=str
        )
        


        # Snakemake Options
        snakemake_options = parser.add_argument_group(f'{CYN}Snakemake Options{NC}')

        snakemake_options.add_argument(
            '--snakemake-dry-run', 
            action='store_true', 
            help="performs snakemake dry run"
        )

        snakemake_options.add_argument(
            '--snakemake-dag', 
            action='store_true', 
            help="performs snakemake dry run and exports DAG"
        )

        snakemake_options.add_argument(
            "--smk-extra-args", 
            help="A string value of extra args for snakemake(don't forget the (''))"
        )


        # other options
        other_conf = parser.add_argument_group(f'{CYN}Other{NC}')

        other_conf.add_argument(
            '--continue', 
            action='store_true', 
            help="continue analysis when re-run"
        )

        other_conf.add_argument(
            '--overwrite', 
            action='store_true', 
            help="overwrite output dir if exsits"
        )

        other_conf.add_argument(
            '-n',"--name", 
            default=f"/RNA_run[{os.environ['start_time']}]", 
            help='Name of files'
        )

        other_conf.add_argument(
            '--verbose', 
            action='store_true', 
            help="verbose"
        )

        other_conf.add_argument(
            '--quit', 
            dest='verbose', 
            action='store_false', 
            help="print many output"
        )

        other_conf.set_defaults(verbose=False)


    def run(self, args):
        pass


class WGS(WorkflowCli):
    name = 'WGS'
    help = '''guap WGS -i/--input dir'''
    usage = ""

    def add_arguments(self, parser):
        parser.add_argument('-i', '--input', help= "Input directory path", metavar='path', type=str)  
        parser.add_argument('-o', '--output', help= "Output directory path", metavar='path', type=str)

        parser.add_argument('-t', '--threads', help= "Number of total no. of threads", type=int)
        parser.add_argument('--threads-index', help= "Number of threads to use during indexing ref", type=int )
        parser.add_argument('--threads-assemble', help= "Number of threads to use per sample assembly", type=int, default = 4)
        parser.add_argument('--threads-align', help= "Number of threads to use per sample aliging", type=int, default = 4)
        parser.add_argument('--threads-kraken', help= "Number of threads to use per sample aliging", type=int)

        parser.add_argument('--index', help= "fasta file to index", type=str, metavar='path')
        parser.add_argument('--index-name', help= "name of index to use", type=str)
        parser.add_argument('--index-path', help= "path to index files", type=str, metavar='path')
        parser.add_argument('--org-name', help= "",type=str)
        parser.add_argument('--org-txid', help= "",type=str)
        parser.add_argument('--org-id', help= "",type=str)
        parser.add_argument('--org-download-direct', help= "Download first found record automaticaly", action='store_true')
        parser.add_argument('--org-auto-download', help= "",action='store_true')

        parser.add_argument('--kraken2-index', help= "", type=str, metavar='path')

        parser.add_argument('--skip-QC', action='store_true', help="Skip  step")
        parser.add_argument('--skip-check-contaminant', action='store_true', help="Skip  step")
        parser.add_argument('--skip-assembly', action='store_true', help="Skip  step")
        parser.add_argument('--skip-assembly-QC', action='store_true', help="Skip  step")
        parser.add_argument('--pre-check-contaminant', action='store_true', help="Skip  step")

        parser.add_argument('--continue', action='store_true', help="continue analysis when re-run")
        parser.add_argument('--overwrite', action='store_true', help="overwrite output dir if exsits")

    def run(self, args):
        print(vars(args))