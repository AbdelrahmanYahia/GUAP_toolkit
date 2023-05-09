from ..workflows import *

class WES(WorkflowCli):
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

        parser.add_argument('-i', '--input', help= "Input directory path", metavar='path', type=os.path.abspath)  
        parser.add_argument('-o', '--output', help= "Output directory path", metavar='path', type=os.path.abspath)

        parser.add_argument('-t', '--threads', help= "Number of total no. of threads", type=int)
        parser.add_argument('--threads-index', help= "Number of threads to use during indexing ref", type=int )
        parser.add_argument('--threads-align', help= "Number of threads to use per sample aliging", type=int, default = 4)
        parser.add_argument('--threads-calling', help= "Number of threads to use per sample variant caller", type=int, default = 4)

        parser.add_argument('--aligner', help = "Choose aligner", choices=["bwa", "bowtie2"], type=str, default='bwa')
        parser.add_argument('--variant-caller', help = "Choose variant caller", choices=["GATK", "mpileup", "lofreq"], type=str, default='GATK')
        parser.add_argument('--bed-file', help='bed file path', metavar='path', type=os.path.abspath)
        parser.add_argument('--reference-fasta',metavar='path', type=os.path.abspath, help="path to reference fasta file")
        parser.add_argument('--reference-index',metavar='path', type=os.path.abspath, help="path to reference index")
        parser.add_argument('--known-variants',metavar='path', type=os.path.abspath, help="path to reference fasta file")

        parser.add_argument('--skip-QC', action='store_true', help="Skip QC step")
        parser.add_argument('--index-reference', action='store_true', help="index reference")

        parser.add_argument('--snakemake-dry-run', action='store_true', help="performs snakemake dry run")
        parser.add_argument('--snakemake-dag', action='store_true', help="performs snakemake dry run and exports DAG")
        parser.add_argument('--snakemake', dest='snakemake', action='store_true', help="Use snakemake workflow (ability to contiue) [default]")
        parser.add_argument('--bash', dest='snakemake', action='store_false', help="Use bash scripts instead of snakemake")

        parser.add_argument('--continue', action='store_true', help="continue analysis when re-run")
        parser.add_argument('--overwrite', action='store_true', help="overwrite output dir if exsits")
        parser.add_argument('--verbose', action='store_true', help="verbose")
        parser.add_argument('--quit', dest='verbose', action='store_false', help="print many output")
        parser.add_argument("-h", "--help", action='store_true')
        parser.add_argument('-n',"--name", default="GUAP-WES-Run", help='') 

        parser.set_defaults(verbose=False)
        parser.set_defaults(snakemake=True)

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

