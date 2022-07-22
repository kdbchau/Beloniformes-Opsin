#alignment_editor.py

"""This script edits multiple sequence alignment (MSA) for PAML (codeml).
Based on certain %, it removes gap columns, removes stop columns,
and edits remaining stop codons (potential errors) into gaps. It
also creates a new fasta MSA and corresponding phylip (.phy) file, and output stats.

Run this script in terminal after cd into the directory with the input_msa.fa file:
> python alignment_editor.py input_msa.fa stopfrac1 gapfrac2 output_msa.fa

This script was created with the help of Dr. Erik Spence from the SciNet HPC Consortium."""

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
number_of_species = len(df1)
original_alignment_bases = len(df1.columns)*3
original_alignment_columns = len(df1.columns)
############################################################################################################

#df_record_check = pd.DataFrame()
df_change = pd.DataFrame() #Create empty dataframe to collect stops changed to gaps
def change_stops(df, codon, frac1):
	"""This function changes stop codons to gaps if they are below a certain percentage (frac1).
    i.e. If frac1 = 0.1, this means that if less than 10% of the column are stop codons, change them
    to a gap instead. Stop codons are TAG, TGA, or TAA."""
	global df_record_check
	global df_change
	
    # First, get a record of which columns are going to have stop-to-gap conversion
	#df_record_check = df.loc[:, (df.isin(codon).sum()/len(df) <= frac1) & (df.isin(codon).sum() > 0)]
	#print(df_record_check)
	
	num_rows = len(df)	
	changes = [col for col in df.columns if 0 < sum((df[col].eq(codon[0]).sum(), df[col].eq(codon[1]).sum(), df[col].eq(codon[2]).sum())) / float(num_rows) < frac1]
	df_change = df[changes]
	#print(df_change)
	
	#For the columns where stop codons appear less than frac1, change only the stop codon into "---" in first dataframe (df1)
	for col in df.columns:
		if 0 < sum((df[col].eq(codon[0]).sum(), df[col].eq(codon[1]).sum(), df[col].eq(codon[2]).sum())) / float(num_rows) < frac1:
			df[col][df[col] == codon[0]] = "---"
			df[col][df[col] == codon[1]] = "---"
			df[col][df[col] == codon[2]] = "---"
	
	return df
	#print(df)
			
frac1 = float(sys.argv[2]) # fraction of stop codons present in column
codon = ["TGA", "TAA", "TAG"] # list of stop codons
df2 = change_stops(df1, codon, frac1)
change_alignment_columns = len(df_change.columns)
#print(df2)
# Alignment length should remain the same here, nothing is removed, only edited.
#############################################################################################################

df_stops = pd.DataFrame() # Create empty dataframe to collect removed stop columns

def remove_stops(df, codon, frac1):
	"""This function deletes columns that have stop codons above a certain fraction (frac1).
    i.e. If frac1 = 0.1, this means that if 10% of the data in a column has a stop codon
    (TAA, TAG, or TGA), delete that column."""
	global df_stops

    # Get a record of stop columns that will be removed from the dataframe (df2).
	num_rows = len(df)
	stops = [col for col in df.columns if sum((df[col].eq(codon[0]).sum(), df[col].eq(codon[1]).sum(), df[col].eq(codon[2]).sum())) / float(num_rows) > frac1]
	df_stops = df[stops]
	#print(df_stops)
	
	# If stop columns has over the % (frac1) of stop codons, delete that column in the dataframe (df2)
	for col in df.columns:
		if sum((df[col].eq(codon[0]).sum(), df[col].eq(codon[1]).sum(), df[col].eq(codon[2]).sum())) / float(num_rows) > frac1:
			df = df.drop(col, axis=1)
	#for cod in codon:
	#	for col in df.columns:
	#		if sum(df[col] == cod) / float(num_rows) > frac1:
	#			df = df.drop(col, axis = 1) # Delete the stop column

	return df

frac1 = float(sys.argv[2]) # fraction of stop codons present in column
codon = ["TGA", "TAA", "TAG"] # list of stop codons
df3 = remove_stops(df2, codon, frac1)
stop_alignment_bases = len(df3.columns)*3
stop_alignment_columns = len(df3.columns)
removed_stop_columns = len(df1.columns) - len(df3.columns)
#print(df3)
#############################################################################################################

df_gaps = pd.DataFrame() # Create empty dataframe to collect removed gap columns
max_codon = pd.DataFrame() # Create empty dataframe to collect most common codon per site

def remove_gaps(df, frac2):
    """This function deletes columns that have more than the fraction (frac2)
    percent of valid entries (hyphens). i.e. If frac2 = 0.2, if 20% of the data in a column
    has gaps, delete that column. This includes "---", "G--", and "GA-" situations."""

    global df_gaps
    global max_codon

    # Get a record of gap columns that will be removed from the dataframe (df3).
    num_rows = len(df)
    gaps = [col for col in df.columns if sum(df[col].str.contains("-")) / float(num_rows) > frac2]
    df_gaps = df[gaps]

    for col in df.columns:
        if sum(df[col].str.contains("-") / float(num_rows)) > frac2:
            df = df.drop(col, axis = 1) # Delete the gap column

    #print(df_gaps)
    return df

frac2 = float(sys.argv[3])
df4 = remove_gaps(df3, frac2) 
max_codon = df4.max()
gap_alignment_bases = len(df4.columns)*3
gap_alignment_columns = len(df4.columns)
removed_gap_columns = len(df3.columns) - len(df4.columns)
##############################################################################################################

def write_codon_data(df, filename):
    """This function writes out a new output MSA.
    Include the filename extension when writing the file name in terminal."""

    with open(filename, 'w') as f:
        for name, seq in df.iterrows():
            codons = ''
            for position in seq: codons += position
            f.write(">" + name + "\n")
            f.write(codons + "\n")

output_msa = sys.argv[4]
write_codon_data(df4, output_msa) # New edited msa output file created
##############################################################################################################

def write_phylip_data(df, filename):
    """This function writes a phylip (.phy) file from the final edited msa. This is always done automatically,
	no need to specify a new outputfile with .phy extension."""

    with open(filename +".phy", 'w') as f:
        # First line includes the number of species, and length of alignment in bases.
        f.write("{} {} \n".format(len(df), len(df.columns)*3))
        # Write in format of species name followed by its sequence on a new line (sequential format)
        for name, seq in df.iterrows():
            codons = ''
            for position in seq: codons += position
            f.write(name + "\n")
            f.write(codons + "\n")

output_phy = sys.argv[4].split(".fa")[0]
write_phylip_data(df4, output_phy)
#############################################################################################################

now = datetime.datetime.now()
today = now.strftime("%B %d, %Y at %H:%M")

total_removed_columns = len(df1.columns) - len(df4.columns)

def write_stats(filename):

	"""This function writes out all the stats of the msa edit."""
	with open(filename.split(".fa")[0]+"_STATS.txt", 'w') as f:
	
		f.write("BASIC GENE STATS FOR: %s. \n" % (filename.split(".fa")[0])) # Gene file name
		
		f.write("Date modified: %s \n" %today)
		
		f.write("\nLength of original alignment is: \n%d bases (%d columns) with %d species." '\n' \
                % (original_alignment_bases, original_alignment_columns, number_of_species))
		
		f.write("\nThe number of columns that had stop codons present in less than %d%% and changed to gaps occurred in %d columns." '\n' % (frac1*100, change_alignment_columns))
				
		f.write("\nRemoved columns with more than %d%% stop codons (TAA,TGA,TAG) present. Otherwise gaps used to replace stop codons." '\n' % (frac1*100))
		
		f.write("Alignment length is now: \n%d bases (%d columns). %d stop columns (%.1f%%) removed. \n" \
                % (stop_alignment_bases, stop_alignment_columns, removed_stop_columns, \
                   removed_stop_columns / float(original_alignment_columns) * 100))
				   
		f.write("\nDeleted columns with more than %d%% missing or incomplete codons (e.g. ---, G-, or GG-)." \
                % (frac2 * 100))
				
		f.write("\nAlignment length is now: \n%d bases (%d columns). %d gap columns (%.1f%%) removed. \n" \
                % (gap_alignment_bases, gap_alignment_columns, removed_gap_columns, \
                   removed_gap_columns / float(original_alignment_columns) * 100))
				   
		f.write("\nIn total, %d columns removed. %.1f%% of columns removed. \n" \
                % (total_removed_columns, (total_removed_columns / float(original_alignment_columns) * 100)))
				
	with open(filename.split(".fa")[0] + "_STATS.txt", 'a') as f:
	
		f.write("\nBelow are the most common codon per column in the final, edited alignment. \n")
		max_codon.to_csv(f, mode='a')

outputfile = sys.argv[4].split(".fa")[0]
write_stats(outputfile)

#########################################################################################################

def write_stat_files(gap_df, stop_df, change_df, filename):
    """This function creates excel and text records of removed column data and changed stop codons."""
    global df_gaps
    global df_stops
    global df_change

# Write an excel version for easy viewing
    writer = pd.ExcelWriter(filename.split(".fa")[0] + "_RemovedAttributes.xlsx")
    gap_df.to_excel(writer, sheet_name='Gap Columns')
    stop_df.to_excel(writer, sheet_name="Stop Columns")
    change_df.to_excel(writer, sheet_name="Change Columns")
    writer.save()

    # Write a txt file just in case, but does not look as pretty
    gap_df.to_csv(filename.split(".fa")[0]+"_RemovedGaps.txt", index=True, sep="\t")
    stop_df.to_csv(filename.split(".fa")[0] + "_RemovedStops.txt", index=True, sep="\t")
    change_df.to_csv(filename.split(".fa")[0] + "_RemovedChange.txt", index=True, sep="\t")

outputfile = sys.argv[4]
write_stat_files(df_gaps, df_stops, df_change, outputfile)
##########################################################################################################

# DONE!
## To turn in terminal, be in the directory with the fasta file and type. For example, to edit an ALIGNED and already in reading frame MSA 
## with parameters of 40% (frac1) stops and 30% (frac2) gaps the following is executed in terminal:
## python alignment_editor.py inputfilename.fa 0.1 0.3 outputfilename.fa
## Once completed, 7 new files are created
