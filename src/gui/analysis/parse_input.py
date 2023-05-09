import os
import yaml 
import pandas as pd
from collections import defaultdict
import logging
import re
from .. import env

def full_match(file):
    filename, file_extension = os.path.splitext(file)
    if "fastq" in filename or "fq" in filename:
        filename , new_ext = os.path.splitext(filename)
        file_extension = new_ext + file_extension
    pattern = "(((\w+)(_S\d+_L\d+_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(R1|R2|r1|r2|read1|read2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(1|2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_S\d+_L\d+)_)(R1|R2|r1|r2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_S\d+_L\d+_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(R1|R2|r1|r2|read1|read2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(1|2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_S\d+_L\d+)_)(R1|R2|r1|r2|read1|read2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))"
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
        exit(1)
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


def full_match(file):
    filename, file_extension = os.path.splitext(file)
    if "fastq" in filename or "fq" in filename:
        filename , new_ext = os.path.splitext(filename)
        file_extension = new_ext + file_extension
    pattern = "(((\w+)(_S\d+_L\d+_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(R1|R2|r1|r2|read1|read2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(1|2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_S\d+_L\d+)_)(R1|R2|r1|r2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+)(_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_S\d+_L\d+_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(R1|R2|r1|r2|read1|read2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(1|2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_S\d+_L\d+)_)(R1|R2|r1|r2|read1|read2)\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))|(((\w+-\w+)(_))(R1|R2|r1|r2|read1|read2)_\d+\.((fastq\.gz)|(fastq)|(fq\.gz)|(fq)))"
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
        exit(1)
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


def validate_sample_input(args):
    outpath = os.path.abspath(args["output"])
    if os.path.exists(outpath):
        print(f"{outpath} exisits!! {args['output']}!")
    # create output path if doesn't exsit 
    if not os.path.exists(outpath):
        print(f"\033[;33;1mcreating working directory @ {outpath} \033[;39;m")
        os.mkdir(outpath)
    

    # create sample table 
    samples = defaultdict(dict)
    sample_names = set()
    pattern = ""
    path = os.path.abspath(args["input"])
    if not os.path.exists(path):
        exit(1)
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



    # create config file 
    with open('config.yaml', 'w') as yaml_file:
        yaml.safe_dump(args, yaml_file, default_flow_style=False, sort_keys=False)

    with open('config.yaml', 'a') as yaml_file:
        yaml_file.writelines(f"path: {path}\n")
        yaml_file.writelines(f"working_dir: {outpath}\n")
        yaml_file.writelines(f"ext: {ext}\n")
        yaml_file.writelines(f"tail: {tail}\n")
        yaml_file.writelines(f"R: {R}\n")
        yaml_file.writelines(f"R1_pattern: _{R}1{tail}{EXT}\n")
        yaml_file.writelines(f"R2_pattern: _{R}2{tail}{EXT}\n")
        yaml_file.writelines(f"compressed: {compressed}\n")
        yaml_file.writelines(f"total_mem: {env.all_mem}\n")
        yaml_file.writelines(f"GUAP_DIR: {env.GUAP_DIR}")

