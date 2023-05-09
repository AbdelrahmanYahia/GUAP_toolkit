"""This script parses the input dir and creates the sample table and config files"""

## need to add lane information recognition
## need to add validate input
import re 
import sys
import yaml 
import shutil
import pandas as pd
from collections import defaultdict
# GUAP modules import
from .globals import *

def check_extension(df): # takes pandas df and returns string
    # checks all files have same extension from pandas df, to use in generete sample table function
    uniques = df['ext'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has multible fastq file extensions, please check directory.")

    else:
        return uniques[0]


def check_PE(df): # takes pandas df and returns string
    """checks all files either single or paried ended from pandas df, to use in generete sample table function"""
    uniques = df['PE'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has both Paired and single files, please check directory.")
    else:
        return uniques[0]


def check_R(df): # takes pandas df and returns string
    """checks all files have same naming patterns from pandas df, to use in generete sample table function"""
    uniques = df['read_num'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has multible fastq file naming patterns, please check directory.")
    else:
        return uniques[0].replace("1","")


def check_pattern(df): # takes pandas df and returns a string 
    """checks all files have same naming patterns from pandas df, to use in generete sample table function"""
    uniques = df['matched_pattern'].unique()
    if len(uniques) > 1:
        glogger.prnt_fatel(f"Your input directory has multible fastq file naming patterns, please check directory.")
    else:
        return uniques[0]
    

def recogize_pattern(file_name): # takes string of fastq file name and returns dict with read info and id
    """ using re to recognize the naming pattern of samples (illumina, srr and general naming patten)"""
    # naming pattern for re 
    patterns = { # ! fix (_|\.) group for R pattern in dict config !
        "illumina": "(((.+)_(S\d+)_(L00\d))_(R1|R2|r1|r2|read1|read2)_(00\d)\.(fastq\.gz|fastq|fq\.gz|fq))",
        "SRR": "(((SRR)(\d+))(_|\.)(1|2|R1|R2|r1|r2|read1|read2)\.(fastq\.gz|fastq|fq\.gz|fq))",
        "general": "(((.+))(_|\.)(1|2|R1|R2|r1|r2|read1|read2)\.(fastq\.gz|fastq|fq\.gz|fq))"
    }

    matched_pattern = None
    ## loop on pattern to and checks whichs one matches 
    ## starting with illumina because general would match any ways
    ## breaks once successful 
    for ptrn_name, pattern in patterns.items():
        try:
            matched = re.match(pattern, file_name) 
        except:
            continue
        
        if bool(matched) :
            matched_pattern = ptrn_name
            break
            
        else:
            continue

    # Extracts sample information
    if matched_pattern == "illumina":
        file_name, sample_name, sample_id, sample_number, read_num, lane, tail, ext = matched.groups()[0], matched.groups()[1], matched.groups()[2], matched.groups()[3], matched.groups()[5], matched.groups()[4], matched.groups()[6], matched.groups()[7]

    elif matched_pattern == "SRR":
        file_name, sample_name, sample_id, read_num, lane, tail, ext = matched.groups()[0], matched.groups()[1], matched.groups()[3], matched.groups()[5], "", "", matched.groups()[6]

    elif matched_pattern == "general":
        file_name, sample_name, sample_id, read_num, lane, tail, ext = matched.groups()[0], matched.groups()[1], matched.groups()[1], matched.groups()[4], "", "", matched.groups()[5]

    else:
        file_name = sample_name = sample_id = read_num = lane = tail = ext = None

    # Returns a dictionary of sample information
    return {
        "file_name": file_name,
        "sample_name": sample_name,
        "sample_id": sample_id,
        "sample_number": sample_number,
        "read_num": read_num,
        "lane": lane,
        "tail": tail,
        "ext": ext,
        "matched_pattern": ptrn_name
    }


def parse_samples(inpath): # takes path return contains fastq files, returns df contains sample information
    ## takes input path
    ## gets the file names containg fastq and fq
    ## performs the recogize_pattern function to 
    ## capture sample information and stores it in 
    ## pandas df
    # input path to absolute path
    path = os.path.abspath(inpath)
    # list all files
    all_files = os.listdir(path)
    samples = defaultdict(dict)
    # takes fastq files only
    for file_name in all_files:
        if os.path.isfile(path + "/" + file_name) and ("fastq" in file_name or "fq" in file_name):
            # Captures the file path and name
            filename, file_extension = os.path.splitext(file_name)
            if "fastq" in filename or "fq" in filename:
                filename, new_ext = os.path.splitext(filename)
                file_extension = new_ext + file_extension
            # recogize_pattern function returns a dictitionary with sample names, id, and read information
            sample_info = recogize_pattern(file_name)

            # get only forward reads and replace the read number to get R2
            # appends sample information to a dict of dicts
            if "1" in sample_info["read_num"]:
                read_2 = sample_info["read_num"].replace("1","2")
                if sample_info["matched_pattern"] == "illumina":
                    read_1 = f"{sample_info['read_num']}_{sample_info['tail']}.f"
                    read_2 = f"{read_2}_{sample_info['tail']}.f"
                else:
                    read_1 = f"{sample_info['read_num']}.f"
                    read_2 = f"{read_2}.f"

                f2 = file_name.replace(read_1, read_2)
                if f2 in all_files:
                    sample_info["file2"] = f2
                    sample_info["PE"] = True
                    samples[sample_info["sample_id"]] = sample_info

                else:
                    sample_info["file2"] = ""
                    sample_info["PE"] = False

    # converts the dict to pandas df and returns the df
    m_samples = samples
    samples= pd.DataFrame(samples).T
    samples = samples.sort_values(by=['sample_id'])

    return samples

# sample table info
### file_name sample_name sample_id read_num lane tail ext matched_pattern file2 PE
### uniques = df['PE'].unique()

def parse_input_args(args): # takes args (object) returns dict of args information
    global glogger
    # set vars
    global verbose
    try:
        verbose = args.verbose
    except:
        verbose = False
    
    glogger.create_console_handler(verbose=verbose)
    glogger.prnt_info("Logger is set to INFO")
    outpath = os.path.abspath(args.output)
    path = os.path.abspath(args.input)

    # create output path if doesn't exsit 
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        glogger.create_file_handler(f"{outpath}/main_log.txt")
        glogger.prnt_info(f"Created Working Dir @ {outpath}")

    else:        
        try:
            if args.overwrite:
                # Remove directory and all its contents
                shutil.rmtree(outpath)
                os.mkdir(outpath)
                glogger.create_file_handler(f"{outpath}/main_log.txt")
                glogger.prnt_warning(f"Overwrited working dir {outpath}")
            else:
                glogger.create_file_handler(f"{outpath}/main_log.txt")

        except:
            glogger.create_file_handler(f"{outpath}/main_log.txt")


    # creates a dict with input args 
    all_args = vars(args)

    if all_mem < 7:
        glogger.prnt_error("Your System doesn't have enough memory to run the analyis")

    # validate samples
    samples = parse_samples(args.input)

    # check n of threads and use all threads if not supplied 
    if args.threads is None:
        glogger.prnt_warning(f"You did not specify any number of threads, I will use all ({all_threads}). ")
        args.threads = all_threads

    elif args.threads < 4:
        glogger.prnt_error(f"GUAP can't use less than 4 threads, would you please change the number of threads?  ")
        new_threads = int(input("Number of threads: (4 and above, enter any thing else to exit) "))
        try:
            int(new_threads)
            if int(new_threads) < 4:
                glogger.prnt_fatel(f"Value is smaller than 4. exiting...")
            else:
                args.threads = new_threads
        except ValueError:
            glogger.prnt_fatel(f"exiting...")


    ext = str(check_extension(samples))
    PE = bool(check_PE(samples))
    R = str(check_R(samples))
    compressed = False
    EXT = ext
    pattern = str(check_pattern(samples))
    tail = "001" if pattern == "illumina" else ""

    # to perform gunzipping 
    if ".gz" in ext:
        compressed = True
        EXT = ext.replace(".gz","")
    
    # check if analysis run before and created sample table 
    if os.path.exists(outpath+"/"+"samples.tsv"):
        glogger.prnt_warning(f"Found an exsiting sample.tsv file in output directory, will not override.")
    else:
        samples.to_csv(outpath+"/"+"samples.tsv",sep='\t')  
    
    extra_info = {
        "path": path,
        "working_dir": outpath,
        "ext": ext,
        "tail": tail,
        "R": R,
        "R1_pattern": f"_{R}1_{tail}.{EXT}",
        "R2_pattern": f"_{R}2_{tail}.{EXT}",
        "compressed" : compressed,
        "total_mem": all_mem,
        "GUAP_DIR": GUAP_DIR,
        "common_rules": f"{GUAP_DIR}/workflows/common/rules/"
    }
    if "decompress" not in all_args:
        all_args.update({"decompress":False})
    all_args.update(extra_info)
    # create config file 
    with open(f'{outpath}/config.yaml', 'w') as yaml_file:
        yaml.safe_dump(all_args, yaml_file, default_flow_style=False, sort_keys=False)
    
    # Check if a command was given as an argument
    if len(sys.argv) > 1:
        # Get the first argument (excluding the script name)
        command = " ".join(sys.argv[1:])
        
        # Write the command to a file
        with open(f"{outpath}/command.txt", "w") as f:
            f.write(command)
        # Write the command to a lastcommand file
        with open(f"{GUAP_DIR}/.last_run.txt", "w") as f:
            f.write(command)
    return all_args
