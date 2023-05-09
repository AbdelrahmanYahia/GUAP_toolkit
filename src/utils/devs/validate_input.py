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
from pathlib import Path

def validate_basic_args(args, parser):
    pass


def validate_16s_args(args, parser, samples, m_samples):
    def qiime_table_checker(df):
        if ("sample-id" or "id" or "sampleid" or "sample id") in df.columns:
            if "#q2:types" in (df.iloc[:1]).values.tolist()[0]:
                return True
            else:
                return False
        else:
            return False
    
    # check if no classifier is selected 
    if args.train is False and args.classifier is None:
        parser.error("\033[;31;1m--classifier is required when --train is not set.")
        sys.exit()

    # check if classifier exsits 
    if not os.path.isfile(args.classifier):
        logging.error(f"\033[;31;1m{args.classifier} Doesn't exsit!")
        exit(1)

    # check if proper classifier is selected with the analysis
    if args.deblur is True and args.choose_classifier == "dada":
        parser.error("\033[;31;1m--choose-classifier dada is not currently supported with deblur.")
        sys.exit()

    if args.use_QIIME2 is True and args.choose_classifier == "dada":
        parser.error("\033[;31;1m--choose-classifier dada is not currently supported with QIIME2 DADA2.")
        sys.exit()

    # check metadata file 
    if args.metadata is None:
        create_meta = input(f"\033[;33;1mNo metadata file supplied,\033[;39;m do I create an empty one with sample IDs? (y/n) ")
        if create_meta == ( 'y' or 'Y'):
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



outpath = os.path.abspath(args.output)
# create output path if doesn't exsit 
if not os.path.exists(outpath):
    print(f"\033[;33;1mcreating working directory @ {outpath} \033[;39;m")
    os.mkdir(outpath)

# check if no arguments 
if not len(sys.argv) > 1:
    print(f"\033[;31;1mERROR: \033[;39;mNo Arguments supplied")
    print(help_message)
    exit(1)

# check for help
if args.help is True:
    print(help_message)
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

