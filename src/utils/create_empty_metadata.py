import os
import argparse
import pandas as pd
import datetime 
now = datetime.datetime.now()
start_time = now.strftime("%Y-%m-%d %H:%M:%S")
os.environ['start_time'] = start_time

from .parse_input import parse_samples

parser = argparse.ArgumentParser(description='create empty metadata')
parser.add_argument('-i', '--input', help= "Input directory path", metavar='path',type=str, required=True)  
parser.add_argument('-o', '--output', help= "Output directory path", metavar='path',type=str, required=True)  
parser.add_argument('--half', action='store_true', help="create have the samples  number cond1 and the second half cond2")
args = parser.parse_args()

# validate samples
samples = parse_samples(args.input)
# Create the new DataFrame
samples_data = pd.DataFrame()

# Add the "sample_id" column
samples_data['sample_id'] = samples['sample_id']

if args.half:
    # Add the "condition" column with values 'cond1' and 'cond2'
    num_rows = len(samples)
    half_rows = num_rows // 2
    samples_data['condition'] = ['cond1'] * half_rows + ['cond2'] * (num_rows - half_rows)
else:
    # Add the "condition" column with "NA" values
    samples_data['condition'] = 'NA'

samples_data.to_csv(args.output,sep=',', index=False)