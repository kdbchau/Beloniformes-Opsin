# merge_seqs.py

"""This script merges two multiple sequence alignment (msa) files. msa1 is the first round msa, and msa2
is the more up-to-date/makeshift MSA (after 2 or 3 rounds). if msa2 contains sequences with gaps for a species, but that same position contains
a base in msa1, replace the gap in msa2 with the base from msa1. And if bases exist in both, if msa1 has base1 but msa2 has base2,
this script will leave msa2 bases alone and so base2 is shown. Therefore, keep the most recent MSA (or most complete) as msa2."""

import os
import pandas as pd
from Bio import SeqIO
import sys
import numpy as np


mydir = os.getcwd()

def read_msa(msa):
    """Convert Oryzias latipes msa (original msa) into dataframe."""

    with open(msa, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            name = record.id # species name
            data = record.seq # DNA sequence

            if not 'df' in locals(): # If no dataframe 'df1' created, make an empty one
                df = pd.DataFrame(columns = range(int(len(data) / 1))) # 1 base per column
                df = df.rename(columns = lambda i: i+1) # Start column indexing at 1 instead of 0
            df.loc[name] = [str(data[i:i+1]) for i in range(0, len(data), 1)]

    return df.replace(regex = True, inplace = False, to_replace = "-", value = np.NaN)

input_msa1 = sys.argv[1]
input_msa2 = sys.argv[2]
df_msa1 = read_msa(input_msa1)
df_msa2 = read_msa(input_msa2)
#print(df_msa1)
#print(df_msa2)


df3 = pd.DataFrame() # create empty dataframe to hold result

def merge_msas(df1, df2):
   """Merge df1 and df2 so that between both, a new dataframe (df3 and thus, msa3) has the least amount of gaps."""

   df3 = df2.combine_first(df1)
   return df3.fillna("-")

merged_msa = merge_msas(df_msa1, df_msa2)
#print(merged_msa)

def merged_msafile(df, filename):
    """Write df3 to a new file. Manually check results too!"""

    with open(filename, 'w') as f:
        for name, seq in df.iterrows():
            bases = ''
            for position in seq: bases += position
            f.write(">" + name + "\n")
            f.write(bases + "\n")

outputfile = sys.argv[3]
merged_msafile(merged_msa, outputfile)
