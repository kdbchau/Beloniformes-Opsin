#ideal_species.py
"""After aligning each species to medaka reference, check to see which species has the most complete sequence alignment for realignment."""

import pandas as pd
import numpy as np
from Bio import SeqIO
import sys
import datetime

def read_codon_data(filename):
    """Convert a fasta MSA into a dataframe for editing."""
    with open(filename, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            name = record.id # species name
            data = record.seq # DNA sequence
            # If sequence length not divisible by 3, make it divisible by 3.
            # Must make sure the MSA is already in proper reading frame prior to this!
            if len(data) % 3 == 2: # i.e. if it ends with "AGT,TAA,NN"
                data = data + "-"  # add a hypen to make it a codon "AGT,TAA,NN-"
            elif len(data) % 3 == 1: # i.e. if it ends with "AGT,TAA,N"
                data = data + "--"    # add 2 hypens to make it a codon "AGT,TAA,N--"

            if not 'df' in locals(): # If no dataframe 'df' created, make an empty one
                df = pd.DataFrame(columns = range(int(len(data) / 3))) # 3 bases per column
                df = df.rename(columns = lambda i: i+1) # Start column indexing at 1 instead of 0
            df.loc[name] = [str(data[i:i+3]) for i in range(0, len(data), 3)]

    return df
	
input_msa = sys.argv[1]
df1 = read_codon_data(input_msa)

df_track = (~df1.apply(lambda x: x.str.contains('-'))).sum(1)	
print(df_track)
