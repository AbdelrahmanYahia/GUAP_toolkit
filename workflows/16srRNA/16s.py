import os
import sys
import argparse
from pathlib import Path
import psutil
import os
import sys
import argparse
from pathlib import Path
import yaml 
import pandas as pd
from collections import defaultdict
import re 
import logging
GUAP_DIR = (os.path.abspath(__file__)).replace("/workflows/16s/16s.py","")
sys.path.append(f"{GUAP_DIR}/bin")
from common_functions import *
from help_messages import *


logging.getLogger().setLevel(logging.INFO)
logging.basicConfig(format='%(levelname)s: %(message)s')
import multiprocessing
all_threads = multiprocessing.cpu_count()
if all_threads < 4:
    logging.error(f"\033[;31;1mYour System doesn't have enough threads to run the analyis")
    exit(1)
all_mem =int( (psutil.virtual_memory().total ) / 1000000000 )
if all_mem < 7:
    logging.error(f"\033[;31;1mYour System doesn't have enough memory to run the analyis")
    exit(1)

parser = argparse.ArgumentParser(description='16s analysis pipeline',
                                 prog='GUAP 16s',
                                 usage='GUAP 16s [-i|--input PATH] [-o|--output PATH] [-m|--metadata FILE] [-c|--classifier FILE]',
                                 add_help=False)
parser.version = '1.2'
parser.add_argument('-i', '--input', help= "Input directory path", metavar='path',type=str)  
parser.add_argument('-o', '--output', help= "Output directory path", metavar='path',type=str)  
parser.add_argument('-t', '--threads', help= "Number of threads", type=int)
parser.add_argument('-ct', '--classifier-threads', help= "Number of threads used by classify-sklearn",default=2,type=int)
parser.add_argument('-m', '--metadata', help= "metadata file full path",required='--downstream' in sys.argv, metavar='path',type=os.path.abspath)  
parser.add_argument('--version', action='version', help="Print version number")
parser.add_argument("--deblur", action="store_true", help="Use deblur deniosing pipeline")
parser.add_argument('-cmt',"--create-metadata-table", action="store_true", help="create metadata empty table from sample names in the input dir")
parser.add_argument('-L',"--deblur-trim-length", default="-1", help='Deblur trim length, use with --deblur flag',type=int,metavar=int) 
parser.add_argument('-tf',"--trunc-f", default=0, help='DADA2 trunclength Forward', type=int) 
parser.add_argument('-tr',"--trunc-r", default=0, help='DADA2 trunclength Reverse', type=int) 
parser.add_argument('-l',"--trim-l", default=0, help='DADA2 trunclength Forward', type=int) 
parser.add_argument('-r',"--trim-r", default=0, help='DADA2 trunclength Reverse', type=int) 
parser.add_argument('-ef',"--maxee-f", default=4, help='DADA2 maxEE Forward', type=int) 
parser.add_argument('-er',"--maxee-r", default=5, help='DADA2 maxEE Reverse', type=int)
parser.add_argument('-mo',"--min-overlap", default=10, help='DADA2 minimum overlap to merge', type=int)
parser.add_argument('-cm',"--chimera-method", default='consensus', help='DADA2 chimera method', type=str, choices=['consensus','pooled'])
parser.add_argument("--remove-primers", action="store_true", help="perform cutadapt to remove primers")
parser.add_argument("--min-length",  type=int,default=50, help='cutadapt (remove-primers) min length') 
parser.add_argument("--downstream", action="store_true", help="Perform Downstream analysis")
parser.add_argument('-d',"--sampling-depth",  type=int, help='Sampling Depth for core-mertic phylogenetic') 
parser.add_argument('-md',"--max-depth",  type=int, help='Maximum Sampling Depth for alpha-rarefaction') 
parser.add_argument("--bashdownstream", action="store_true",help="berform downstream analysis on analysis run by --bash argument (downstream only)")
parser.add_argument('-fp',"--forward-primer", help='forward primmer for cutadapt or train set ',required='--train' in sys.argv or '--remove-primers' in sys.argv)
parser.add_argument('-rp',"--reverse-primer", help='reverse primmer for cutadapt or train set', required= '--train' in sys.argv or '--remove-primers' in sys.argv)
parser.add_argument("--train", action="store_true",help="Train set for qiime2 naive classifier")
parser.add_argument("--t-i-seqs", help='Train set input sequences (qza)',required='--train' in sys.argv,type=os.path.abspath)
parser.add_argument("--t-i-taxa", help='Train set input taxa (qza)',required='--train' in sys.argv,type=os.path.abspath)
parser.add_argument("--trainset-minlength", help='Train set input sequences (qza)',required='--train' in sys.argv,type=os.path.abspath)
parser.add_argument("--trainset-maxlength", help='Train set input sequences (qza)',required='--train' in sys.argv,type=os.path.abspath)
parser.add_argument('--choose-classifier', help= "Choose to classify taxanomy [dada or qiime] defualt: qiime", default="qiime", choices=['dada', 'qiime'])
parser.add_argument("-c", "--classifier", help='Classifier path', metavar='path',type=os.path.abspath)
parser.add_argument('--trimmomatic', dest='trimmomatic', action='store_true', help="Use trimmomatic [default]")
parser.add_argument("--trim-min-length",  type=int,default=30, help='trimmomatic min length') 
parser.add_argument('--skip-trimmomatic', dest='trimmomatic', action='store_false', help="Do NOT use trimmomatic")
parser.add_argument('--skip-QC', action='store_true', help="Skipp Fastqc step")
parser.add_argument('--snakemake-dry-run', action='store_true', help="performs snakemake dry run")
parser.add_argument('--snakemake-dag', action='store_true', help="performs snakemake dry run and exports DAG")
parser.add_argument('-n',"--name", default="/GUAP-16s-dada2", help='Name of files') 
parser.add_argument('--snakemake', dest='snakemake', action='store_true', help="Use snakemake workflow (ability to contiue) [default]")
parser.add_argument('--use-QIIME2', action='store_true', help="Use QIIME2 to perform DADA2")
parser.add_argument('--verbose', action='store_true', help="verbose")
parser.add_argument('--quit', dest='verbose', action='store_false', help="print many output")
parser.add_argument('--export-figs', action='store_true', help="export phyloseq figures (alpha, beta and bar plots)")
parser.add_argument('--continue', action='store_true', help="continue analysis when re-run")
parser.add_argument('--export-figs-only', action='store_true', help="export phyloseq figures (alpha, beta and bar plots)")
parser.add_argument( '--condition-name', help= "condition name in metadata file", required='--export-figs' in sys.argv, type=str)
parser.add_argument('--bash', dest='snakemake', action='store_false', help="Use bash scripts instead of snakemake")
parser.set_defaults(snakemake=True)
parser.set_defaults(trimmomatic=True)
parser.set_defaults(verbose=False)
parser.add_argument("-h", "--help", action='store_true')

args = parser.parse_args()

def qiime_table_checker(df):
    if ("sample-id" or "id" or "sampleid" or "sample id") in df.columns:
        if "#q2:types" in (df.iloc[:1]).values.tolist()[0]:
            return True
        else:
            return False
    else:
        return False

# check if no arguments 
if not len(sys.argv) > 1:
    print(f"\033[;31;1mERROR: \033[;39;mNo Arguments supplied")
    print(help_message_16s)
    exit(1)

# check for help
if args.help is True:
    print(help_message_16s)
    exit(0)

# check n of threads and use all threads if not supplied 
if args.threads is None:
    logging.warning(f"\033[;33;1mYou did not specify any number of threads, I will use all ({all_threads}). \033[;39;m")
    args.threads = all_threads

elif args.threads < 4:
    logging.error(f"\033[;31;1mGUAP 16s can't use less than 4 threads, would you please change the number of threads?  \033[;39;m")
    new_threads = int(input("Number of threads: (4 and above, enter any thing else to exit) "))
    try:
        int(new_threads)
        if int(new_threads) < 4:
            logging.error(f"\033[;31;1mValue is smaller than 4. exiting...")
            exit(1)
        else:
            args.threads = new_threads
    except ValueError:
        logging.error(f"\033[;31;1mexiting...")
        exit(1)


outpath = os.path.abspath(args.output)
# create output path if doesn't exsit 
if not os.path.exists(outpath):
    print(f"\033[;33;1mcreating working directory @ {outpath} \033[;39;m")
    os.mkdir(outpath)


# create sample table 
samples = defaultdict(dict)
sample_names = set()
pattern = ""
path = os.path.abspath(args.input)
samples = defaultdict(dict)
allfiles = os.listdir(path)

for file in allfiles:
    if os.path.isfile(path+"/"+file):
        if "fastq" in file or "fq" in file:
            if full_match(file) is False:
                continue
            else:
                sample_name, filename, file_extension, matched, matched_groups = full_match(file)
                if "1" in matched_groups[4]:
                    f2 = file.replace("_R1","_R2")
                    f2 = f2.replace("_1","_2")
                    samples[matched_groups[1]]["Sample ID"] = matched_groups[2]
                    samples[matched_groups[1]]["Sample rest of name"] = matched_groups[3]
                    samples[matched_groups[1]]["Extension"] = file_extension
                    samples[matched_groups[1]]["R1"] = file
                    samples[matched_groups[1]]["R2"] = ""
                    if "_001." in file:
                        samples[matched_groups[1]]["tail"] = "_001"
                    else:
                        samples[matched_groups[1]]["tail"] = ""
                    if f2 in allfiles:
                        samples[matched_groups[1]]["R2"] = f2
                        samples[matched_groups[1]]["PE"] = True
                    else:
                        samples[matched_groups[1]]["R2"] = ""
                        samples[matched_groups[1]]["PE"] = False
                    samples[matched_groups[1]]["path"] = (path)
                    if "R" in matched_groups[4]:
                        samples[matched_groups[1]]["R"] = "R"
                    else:
                        samples[matched_groups[1]]["R"] = ""
                    
                elif "2" in matched_groups[4]:
                    continue
                else:
                    samples[matched_groups[1]]["Sample ID"] = matched_groups[2]
                    samples[matched_groups[1]]["Sample rest of name"] = matched_groups[3]
                    samples[matched_groups[1]]["Extension"] = file_extension
                    samples[matched_groups[1]]["R1"] = file
                    samples[matched_groups[1]]["R2"] = ""
                    samples[matched_groups[1]]["PE"] = False
                    if "R" in matched_groups[4]:
                        samples[matched_groups[1]]["R"] = "R"
                    else:
                        samples[matched_groups[1]]["R"] = ""

                    if "_001." in file:
                        samples[matched_groups[1]]["tail"] = "_001"
                    else:
                        samples[matched_groups[1]]["tail"] = ""
                    samples[matched_groups[1]]["path"] = (path)
        else:
            continue
m_samples = samples
samples= pd.DataFrame(samples).T
samples = samples[["Sample ID", "Sample rest of name", "Extension", "R1", "R2", "tail", "PE", "path", "R"]]
samples = samples.sort_values(by=['Sample ID'])

ext = str(check_extension(samples)[0])
tail = str(check_tail(samples)[0])
PE = bool(check_PE(samples))
R = str(check_R(samples)[0])
compressed = False
EXT = ext

# to perform gunzipping 
if ".gz" in ext:
    compressed = True
    EXT = ext.replace(".gz","")
    
# check for paired end files 
if PE is False:
    logging.error(f"\033[;31;1mSamples not in pairs... exiting...")
    exit(1)

# check if analysis run before and created sample table 
if os.path.exists(outpath+"/"+"samples.tsv"):
    logging.warning(f"\033[;33;1mFound an exsiting sample.tsv file in output directory, will not override.\033[;39;m")
else:
    samples.to_csv(outpath+"/"+"samples.tsv",sep='\t')    

# check metadata file 
if args.metadata is None or args.create_metadata_table:
    if args.create_metadata_table:
        pass
    else:
        create_meta = input(f"\033[;33;1mNo metadata file supplied,\033[;39;m do I create an empty one with sample IDs? (y/n) ")
        if create_meta == ( 'y' or 'Y'):
            pass
        else:
            exit(1)
    logging.warning(f"\033[;33;1mI will create one at '{args.output}', please re-run the analysis with -m {args.output}/sample-metada.tsv after modifing the file (check https://docs.qiime2.org/2021.11/tutorials/metadata/)\033[;39;m")
    header = pd.DataFrame({"1": ["#q2:types", "categorical"]}).T
    header.columns = ["sample-id", "condition"]

    new_samples = {}
    c = 2
    for sample in m_samples.values():
        new_samples[c] = [sample["Sample ID"], "BLANK"]
        c += 1
    samples_df = pd.DataFrame(new_samples).T
    samples_df.columns = ["sample-id", "condition"]
    samples_df = samples_df.sort_values(["sample-id"])
    metadata_file = pd.concat([header, samples_df])
    metadata_file.to_csv(outpath+"/"+"sample-metadata.tsv",sep='\t',index=False) 
    args.metadata = f"{outpath}/sample-metadata.tsv"
    exit(1)

else:
    if not os.path.isfile(args.metadata):
        logging.error(f"\033[;31;1m{args.metadata} Doesn't exsit!")
        exit(1)
    else:
        metadataname, metadataextension = os.path.splitext(args.metadata)
        if metadataextension == (".tsv"):
            the_file = pd.read_csv(args.metadata,sep="\t")
        elif metadataextension == (".csv"):
            the_file = pd.read_csv(args.metadata,sep="\t")
        else:
            logging.error(f"\033[;31;1m{args.metadata} Have strange extension (use: tsv or csv)")
            exit(1)

        rest_of_cols = the_file.columns[1:]

        if False in samples['Sample ID'].isin(the_file.T.iloc[0]).tolist():
            logging.error(f"\033[;31;1m Please check {args.metadata} sample IDs!\033[;39;m")
            print(f"\033[;31;1mNote: \033[;39;mYou can re-run your code WITHOUT suppling a metadata file and GUAP will ask to create an empty one for you with the samples IDs.")
            exit(1)        
        else:
            if qiime_table_checker(the_file):
                if metadataextension == (".csv"):
                    the_file.to_csv(outpath+"/"+"sample-metadata.tsv",sep='\t',index=False) 
                    args.metadata = "{outpath}/sample-metadata.tsv"
            else:
                print(f"\033[;31;1mNote: \033[;39;mModifing sample metadata file to be compatible with QIIME2")
                columns_names = rest_of_cols.to_list()
                columns_names.insert(0, 'sample-id')
                columns_names
                t = f'{(len(columns_names) -1 )* "categorical,"}'.split(",")[:-1]
                t.insert(0, "#q2:types")
                header_N = pd.DataFrame({1: t}).T
                header_N.columns = columns_names
                the_file.columns = columns_names
                last_file = pd.concat([header_N, the_file])
                last_file.to_csv(outpath+"/"+"sample-metadata.tsv",sep='\t',index=False) 
                args.metadata = f"{outpath}/sample-metadata.tsv"


# check if no classifier is selected 
if args.train is False and args.classifier is None:
    parser.error("\033[;31;1m--classifier is required when --train is not set.")
    sys.exit()

# check if proper classifier is selected with the analysis
if args.deblur is True and args.choose_classifier == "dada":
    parser.error("\033[;31;1m--choose-classifier dada is not currently supported with deblur.")
    sys.exit()
if args.use_QIIME2 is True and args.choose_classifier == "dada":
    parser.error("\033[;31;1m--choose-classifier dada is not currently supported with QIIME2 DADA2.")
    sys.exit()


# check if classifier exsits 
if not os.path.isfile(args.classifier):
    logging.error(f"\033[;31;1m{args.classifier} Doesn't exsit!")
    exit(1)


# create config file 
with open('config.yaml', 'w') as yaml_file:
    yaml.safe_dump(vars(args), yaml_file, default_flow_style=False, sort_keys=False)

with open('config.yaml', 'a') as yaml_file:
    yaml_file.writelines(f"path: {path}\n")
    yaml_file.writelines(f"working_dir: {outpath}\n")
    yaml_file.writelines(f"ext: {ext}\n")
    yaml_file.writelines(f"tail: {tail}\n")
    yaml_file.writelines(f"R: {R}\n")
    yaml_file.writelines(f"R1_pattern: _{R}1{tail}{EXT}\n")
    yaml_file.writelines(f"R2_pattern: _{R}2{tail}{EXT}\n")
    yaml_file.writelines(f"compressed: {compressed}\n")
    yaml_file.writelines(f"total_mem: {all_mem}\n")
    yaml_file.writelines(f"GUAP_DIR: {GUAP_DIR}")
