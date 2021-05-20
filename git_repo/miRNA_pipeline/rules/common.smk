import os
import pandas as pd

sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file,index_col=0)
SAMPLES = list(SampleTable.index)