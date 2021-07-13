## This script takes an input directory with counts data from 
## mapped to genome and mapped to mirbase, it adds the counts 
## data and export a table of counts for all samples.

# pandas package is required to be installed if not installed
import os
import glob
import argparse
import pandas as pd

# create parser for command line options 
parser = argparse.ArgumentParser(description='add file to run script on')
parser.add_argument('-i', '--input') # takes input files directory 
parser.add_argument('-o', '--output') # takes output file name 
args = parser.parse_args()

# change working directory to files directory 
os.chdir(args.input)

all_samples={} # empty dictionary for all samples
count_mirna_files = glob.glob("*.counts") # get all file names 
# print(count_mirna_files)
for sample in count_mirna_files: # create a key with file name and value of empty dictionary
    sample_name = sample.replace("_L001_R1_001.counts", "")
    all_samples[sample_name] = {}

for sample in all_samples: # get all counts data from mirbase step
    with open(sample+'_L001_R1_001.counts') as sfile: # opens the file and read lines 
        lines=sfile.readlines()
        for line in lines:
            line = line.strip().split("\t") # removes \n and make list first value is miRNA name and second value is count data
            all_samples[sample].update({str(line[0]):float(line[1])}) # adds key in nested dictionary as miRNA name and value is count data

# for sample in all_samples: # gets all count data from mapped to genome step 
#     with open(sample+'.human.counts') as rfile: # opens files and read lines 
#         lines=rfile.readlines()
#         for line in lines:
#             line = line.strip().split("\t")
#             if not line[0] in all_samples[sample]: # checks if miRNA name is abscent and adds it to the dictionary 
#                 all_samples[sample].update({str(line[0]):float(line[1])})
#             else:
#                 all_samples[sample][line[0]] +=  float(line[1]) # adds the values from mirbase file and mapped to genome files 

my_df = pd.DataFrame(all_samples) # converts nested dictionaries to dataframe using pandas 
my_df.to_csv(args.output) # exports csv file 