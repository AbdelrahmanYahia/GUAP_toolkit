import os
import sys
import argparse
from pathlib import Path
import yaml 

# create parser for command line options 
parser = argparse.ArgumentParser(description='16s analysis pipeline',
                                 prog='GUAP 16s',
                                 usage='%(prog)s -i [input] -o [output] -m [metadata] -c [classifier]',
                                 epilog='''This workflow is still under development
                                 thank you for testing for any bugs please contact:
                                 Abdelrhman @ abdelrhman.yahia@57357.org
                                 ''')
parser.version = '0.1'
parser.add_argument('-i', '--input', help= "Input directory path", required=True, metavar='path',type=str)  
parser.add_argument('-o', '--output', help= "Output directory path", required=True, metavar='path',type=str)  
parser.add_argument('-t', '--threads', help= "Number of threads", type=int, default=8)
parser.add_argument('-s', '--sample-table', help= "QIIME2 sample table")  
parser.add_argument('-m', '--metadata', help= "metadata file full path", required=True, metavar='path',type=str)  
parser.add_argument('--version', action='version', help="Print version number")
parser.add_argument("--deblur", action="store_true", help="Use deblur deniosing pipeline")
parser.add_argument('-L',"--deblur-trim-length", default=300, help='Deblur trim length, use with --deblur flag',required='--deblur' in sys.argv) 
parser.add_argument('-tf',"--trunc-f", default=0, help='DADA2 trunclength Forward', type=int) 
parser.add_argument('-tr',"--trunc-r", default=0, help='DADA2 trunclength Reverse', type=int) 
parser.add_argument('-l',"--trim-l", default=10, help='DADA2 trunclength Reverse', type=int) 
parser.add_argument('-r',"--trim-r", default=10, help='DADA2 trunclength Reverse', type=int) 
parser.add_argument('-ef',"--maxee-f", default=2, help='DADA2 maxEE Forward', type=int) 
parser.add_argument('-er',"--maxee-r", default=2, help='DADA2 maxEE Reverse', type=int) 
parser.add_argument("--train", action="store_true")
parser.add_argument('-fp',"--forward-primer", help='Train set forward primmer',required='--train' in sys.argv) 
parser.add_argument('-fr',"--reverse-primer", help='Train set reverse primmer',required='--train' in sys.argv)
parser.add_argument('-xl',"--max-length", default=550, help='Train set max length', type=int,required='--train' in sys.argv)
parser.add_argument('-ml',"--min-length", default=150, help='Train set min length', type=int,required='--train' in sys.argv)
parser.add_argument("--t-i-seqs", help='Train set input sequences (qza)',required='--train' in sys.argv)
parser.add_argument("--t-i-taxa", help='Train set input taxa (qza)',required='--train' in sys.argv)
parser.add_argument("-c", "--classifier", help='Classifier path', metavar='path',type=str)
parser.add_argument("-I", "--interactive", help='Interactive User input', action="store_true", default=False)
parser.add_argument( "--trimmomatic", help='Use Trimmomatic', action="store_true", default=False)
parser.add_argument('-n',"--name", default="/GUAP-16s-dada2", help='Name of files') 

args = parser.parse_args()
if args.train is False and args.classifier is None:
    parser.error("--classifier is required when --train is not set.")
    sys.exit()


with open('result.yml', 'w') as yaml_file:
    yaml.dump(vars(args), yaml_file, default_flow_style=False)