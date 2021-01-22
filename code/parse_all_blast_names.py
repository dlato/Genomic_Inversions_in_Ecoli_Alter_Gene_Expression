#!/usr/local/bin/python3
#
# will create a csv of all gene names/identifiers for blast proteome
# info (downloaded from UniProt)
# to run:parse_all_blast_names.py *.txt
#
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

#get file name
file_name = sys.argv[1]
sub_var = re.compile("blast_gene_names\.txt")
file_name = sub_var.sub("",file_name)

max_names = 0
all_names = []
#open blast proteome info file
f = open(sys.argv[1],"r")
for line in f:
    tmp_l = re.split("\t",line)
    tnames = tmp_l[4]
    names = re.split(" ",tnames)
    nlen = len(names)
    if (max_names < nlen):
        max_names = nlen
    all_names.append(names)

#opening file to write to
out_file_name = file_name + "_all_genes_info_blast.csv"
out_file = open(out_file_name,"w+")
for n in all_names:
    r = 0
    if (len(n) < max_names):
        r = max_names - len(n)
    for i in range(0,len(n)):
        out_file.write(str(n[i]))
        out_file.write(",")
    out_file.write(''.join(","*r))
    out_file.write("\n")


