#!/usr/local/bin/python3
#
# Program will go through the inverted mafft blocks and grab the block
# number, taxa names, start position, end position, and if it is
# reversed.
#
# to run: ~/CODE/inversions_block_info.py *.mafft
#
# file columns:
# block_name	taxa_name	block_start	block_end	rev_comp	inversion	
#
# NOTE: 1 = inversion, and 1 = reversed complement
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import Bio
import numpy as np
from Bio import SeqIO
#import pdb # FOR DEBUGGING ONLY import pdb
#import pandas as pd #importing pandas
#import numpy as np #importing Numpy
#import csv
import os #importing operating system commands
import glob #importing glob which helps read files
from Bio import AlignIO
#from functools import partial
#from collections import defaultdict
# FOR DEBUGGING ONLY pdb.set_trace()



#below will open the mafft file and get the aln info
align = AlignIO.read(sys.argv[1], "fasta")
#aln_len = len(align[ref_taxa].seq)
#ref_tax_test = align[ref_taxa].id
#print("REF TAX ID",ref_tax_test)

#getting the block name so we can print to 2 separate cod and non-cod
#files later
tmp_var = sys.argv[1]
sub_var = re.compile("\.txt")
tmp_var = sub_var.sub("",tmp_var)
block_name = tmp_var
sys.stderr.write(block_name)
sys.stderr.write("\n")

#edit the taxa IDs to get the start and end pos of each block
block_starts = []
block_ends = []
taxa_names = []
rev_comp = []
for i in range(0,len(align)):
    tmp_id = align[i].id
    #store 'a' which is the ':' character so we can sub it out later
    a = re.compile("\:")
    tmp_id = a.sub(" ",tmp_id)
    #store 'b' which is the '-' character so we can sub it out later
    b = re.compile("\-")
    tmp_id = b.sub(" ",tmp_id)
    #split the ID so we separate the actual ID from the start and end pos
    split_id = re.split(" ",tmp_id)
    #append start and stop and taxa name to list so we can use it later
    taxa_names.append(split_id[0])
    block_starts.append(split_id[1])
    block_ends.append(split_id[2])
    #storing info about if each seq is a reverse complement
    if (len(split_id) > 3):
        rev_comp.append(1)
        tmp_rev_comp = split_id[3]
    else:
        rev_comp.append(0)

block_starts = list(map(int,block_starts))
block_ends = list(map(int,block_ends))
num_taxa = len(block_ends)
print("taxa_names")
print(taxa_names)
print("rev_comp")
print(rev_comp)

sys.stderr.write("done reading in the files!\n")

#opening file to write to
out_file_name = block_name + "_info.txt"
out_file = open(out_file_name,"w+")
#out_file.write("Block\ttaxa\tstart\tend\trev_comp\tinversion\n")
for t in range(0,len(taxa_names)):
    out_file.write(block_name)
    out_file.write("\t")
    out_file.write(taxa_names[t])
    out_file.write("\t")
    out_file.write(str(block_starts[t]))
    out_file.write("\t")
    out_file.write(str(block_ends[t]))
    out_file.write("\t")
    out_file.write(str(rev_comp[t]))
    out_file.write("\t")
#    out_file.write("1\n")
    out_file.write("0\n")
