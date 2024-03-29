from ..workflows import *
import re

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
        super().run(args, "RNAseq")