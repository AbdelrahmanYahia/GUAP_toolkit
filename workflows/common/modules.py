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

# create parser for command line options 
parser = argparse.ArgumentParser(description='Modules',
                                 prog='GUAP modules',
                                 usage='%(prog)s -i [input] -o [output] -m [module]',
                                 epilog='''This workflow is still under development
                                 thank you for testing for any bugs please contact:
                                 Abdelrhman @ abdelrhman.yahia@57357.org
                                 ''')
parser.version = '1.2'
parser.add_argument('-i', '--input', help= "Input directory path", required=True, metavar='path',type=str)  
parser.add_argument('-o', '--output', help= "Output directory path", required=True, metavar='path',type=str)  
parser.add_argument('-t', '--threads', help= "Number of threads", type=int)
parser.add_argument('-m', '--module', help= "module name", default="qc", choices=['qc', 'align', 'cutadapt', 'trimmomatic', 'GUAP'])

parser.add_argument('--pe', dest='pe', action='store_true', help="analysis is in paired end mode [default]")
parser.add_argument('--se', dest='pe', action='store_false', help="analysis in single end mode")
parser.set_defaults(pe=True)

parser.add_argument("--GUAP-args", type=str, help='custom code')


parser.add_argument("--fastqc", action="store_true")


parser.add_argument('--trimmomatic', dest='trimmomatic', action='store_true', help="Use trimmomatic")
parser.add_argument("--trim-min-length",  type=int,default=30, help='trimmomatic min length') 
parser.add_argument("--trim-leading",  type=int,default=30, help='trimmomatic LEADING') 
parser.add_argument("--trim-trailing",  type=int,default=30, help='trimmomatic min TRAILING') 
parser.add_argument("--trim-sliding",  type=int,default=10, help='trimmomatic slidingwindow') 
parser.add_argument("--trim-extra",  type=str, help='extra arguments for trimmomatic', default="") 

parser.add_argument("--cutadapt", action="store_true")
parser.add_argument('-fp',"--forward-primer", help='Cutadapt forward primmer', default="CCTACGGGNGGCWGCAG") 
parser.add_argument('-fr',"--reverse-primer", help='Cutadapt reverse primmer', default="GACTACHVGGGTATCTAATCC")
parser.add_argument('-ml',"--min-length",  type=int,default=150, help='Cutadapt min length') 
parser.add_argument("--cut-extra",  type=str, help='extra arguments for cutadapt', default="") 

parser.add_argument("--align", action="store_true", help="perform alignemnt")
parser.add_argument( '--index', help= "Index for aligner", required='--align' in sys.argv)
parser.add_argument('-a', '--aligner', help= "aligner name", default="bwa", choices=['bwa', 'bowtie2', 'bowtie', 'hisat2'])
parser.add_argument("--sam", action="store_true", help="keep sam file")
parser.add_argument("--align-extra",  type=str, help=f"extra arguments for [aligner]", default="") 

parser.add_argument('--version', action='version', help="Print version number")

args = parser.parse_args()
if args.module == "qc":
    args.fastqc = True

if args.module == "cutadapt":
    args.cutadapt = True

if args.module == "trimmomatic":
    args.trimmomatic = True

if args.module == "align":
    args.align = True

def full_match(file):
    filename, file_extension = os.path.splitext(file)
    if "fastq" in filename or "fq" in filename:
        filename , new_ext = os.path.splitext(filename)
        file_extension = new_ext + file_extension
    pattern = "(((\w+)(_S\d+_L\d+_))(R1|R2|r1|r2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(R1|R2|r1|r2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(1|2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_S\d+_L\d+)_)(R1|R2|r1|r2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(R1|R2|r1|r2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_S\d+_L\d+_))(R1|R2|r1|r2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(R1|R2|r1|r2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(1|2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_S\d+_L\d+)_)(R1|R2|r1|r2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(R1|R2|r1|r2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))"
    matched = re.match(pattern, file)

    if bool(matched) :
        matched_groups = matched.groups()
        a = []
        for val in matched_groups:
            if val != None :
                a.append(val)
        sample_name = matched.group(1)
        return sample_name, filename, file_extension, matched, sorted(set(a), key=a.index)

    else:
        logging.error(f"\033[;31;1m{file} has an unfamilier naming pattern.")
        return False

def check_extension(df):
    uniques = df['Extension'].unique()
    if len(uniques) > 1:
        logging.error(f"Y\033[;31;1mour input directory has multible fastq file extensions, please check directory.")
        exit(1)

    else:
        return uniques

def check_tail(df):
    uniques = df['tail'].unique()
    if len(uniques) > 1:
        logging.error(f"\033[;31;1mYour input directory has multible fastq file naming patterns, please check directory.")
        exit(1)

    else:
        return uniques

def check_PE(df):
    uniques = df['PE'].unique()
    if len(uniques) > 1:
        logging.error(f"\033[;31;1mYour input directory has both Paired and single files, please check directory.")
        exit(1)
    else:
        return uniques
def check_R(df):
    uniques = df['R'].unique()
    if len(uniques) > 1:
        logging.error(f"\033[;31;1mYour input directory has multible fastq file naming patterns, please check directory.")
        exit(1)
    else:
        return uniques


if args.threads is None:
    logging.warning(f"\033[;33;1mYou did not specify any number of threads, I will use all ({all_threads}). \033[;39;m")
    args.threads = all_threads
elif args.threads < 4:
    logging.error(f"\033[;31;1mGUAP 16s can't use less than 4 threads, would you please change the number of threads?")
    new_threads = input("Number of threads: (4 and above, enter any thing else to exit) ")
    try:
        int(new_threads)
        if new_threads < 4:
            logging.error(f"\033[;31;1mValue is smaller than 4. exiting...")
            exit(1)
        else:
            args.threads = new_threads
    except ValueError:
        logging.error(f"\033[;31;1mexiting...")
        exit(1)

outpath = os.path.abspath(args.output)

if not os.path.exists(outpath):
    print(f"\033[;33;1mcreating working directory @ {outpath} \033[;39;m")
    os.mkdir(outpath)
with open('config.yaml', 'w') as yaml_file:
    yaml.safe_dump(vars(args), yaml_file, default_flow_style=False, sort_keys=False)

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


samples= pd.DataFrame(samples).T
samples = samples.sort_values(by=['Sample ID'])

ext = str(check_extension(samples)[0])
tail = str(check_tail(samples)[0])
PE = bool(check_PE(samples))
R = str(check_R(samples)[0])
compressed = False

if ".gz" in ext:
    compressed = True

if PE is False:
    logging.error(f"\033[;31;1mSamples not in pairs... exiting...")
    exit(1)

if os.path.exists(outpath+"/"+"samples.tsv"):
    logging.warning(f"\033[;33;1mFound an exsiting sample.tsv file in output directory, will not override.\033[;39;m")
else:
    samples.to_csv(outpath+"/"+"samples.tsv",sep='\t')    

GUAP_DIR = (os.path.abspath(__file__)).replace("/bin/modules/modules.py","")

with open('config.yaml', 'a') as yaml_file:
    yaml_file.writelines(f"path: {path}\n")
    yaml_file.writelines(f"working_dir: {outpath}\n")
    yaml_file.writelines(f"ext: {ext}\n")
    yaml_file.writelines(f"tail: {tail}\n")
    yaml_file.writelines(f"R: {R}\n")
    yaml_file.writelines(f"compressed: {compressed}\n")
    yaml_file.writelines(f"total_mem: {all_mem}\n")
    yaml_file.writelines(f"GUAP_DIR: {GUAP_DIR}")
