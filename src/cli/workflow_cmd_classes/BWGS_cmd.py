from ..workflows import *

class BWGS(WorkflowCli):
    name = 'BWGS'
    help = '''guap BWGS -i/--input dir'''
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
        parser.add_argument('--skip-check-contaminant', action='store_true', help="Skip check contaminant step")
        parser.add_argument('--skip-assembly', action='store_true', help="Skip assembly step")
        parser.add_argument('--skip-assembly-QC', action='store_true', help="Skip QC step")
        parser.add_argument('--pre-check-contaminant', action='store_true', help="Skip pre check contaminant step")

        parser.add_argument('--continue', action='store_true', help="continue analysis when re-run")
        parser.add_argument('--overwrite', action='store_true', help="overwrite output dir if exsits")

    def run(self, args):
        print(vars(args))