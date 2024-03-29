from ..workflows import *

class WES(WorkflowCli):
    name = 'WES'
    help = '''guap WES -i/--input dir'''
    usage = f"""{YEL}Basic Run Usage example:{NC}
    guap WES -i indir -o outdir --bed-file file \
            --aligner bwa --reference-fasta fasta.fasta \
            --reference-index indexpath 
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
            '--threads', 
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
            '--bed-file', 
            help='bed file path', 
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
            help = "Choose aligner [default = bwa]",
            choices=["bwa", "bowtie2"], 
            metavar = "bwa|bowtie2",
            type=str,
            default='bwa'
        )

        aligner_conf.add_argument(
            '--threads-index', 
            help= "Number of threads to use during indexing ref [default = 4]", 
            type=int, 
            default= 4,
            metavar = "N",
        )

        aligner_conf.add_argument(
            '--threads-align', 
            help= "Number of threads to use during sample alignment [default = 4]", 
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

        variant_caller_conf = parser.add_argument_group(f'{BLU}Variant caller configuration{NC}')

        variant_caller_conf.add_argument(
                '--variant-caller', 
                help = "Choose variant caller", 
                choices=["GATK", "mpileup", "lofreq"], 
                type=str, 
                default='GATK'
        )

        variant_caller_conf.add_argument(
            '--known-variants',
            metavar='path', 
            type=os.path.abspath, 
            help="path to reference fasta file"
        )

        variant_caller_conf.add_argument(
            '--caller-extra-args', 
            help = "Extra arguments for caller, use it with no spaces and add = ( --caller-extra-args='-arg1 -arg2' ) ",
            type=str,
            metavar = "'-args'",
            default=""
        )

        variant_caller_conf.add_argument(
            '--threads-calling', 
            help= "Number of threads to use during variant calling [default = 4]", 
            type=int, 
            default= 4,
            metavar = "N",
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
            default=f"guap_run[{os.environ['start_time']}]", 
            metavar = 'str',
            help=f"Name of files [ default = guap_run[date time] ]"
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

    def run(self, args):
            super().run(args, "WES")