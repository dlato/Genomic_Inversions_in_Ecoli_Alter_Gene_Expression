#!/usr/local/bin/python3
#
# Program will go through the genbank file and will obtain the locus
# ID for each gene as well as the start and end positions
# run as "get_gene_pos_locus.py *.gbk > *gene_locus_location" 
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import Bio
import numpy as np
from Bio import SeqIO
# FOR DEBUGGING ONLY pdb.set_trace()

#faa_filename = sys.argv[2]
#output_handle = open(faa_filename, "w")
#output_handle.write("gbk_start\tgbk_end\tgbk_midpoint\tgbk_gene_id\tgbk_locus_tag\tgbk_old_locus_tag\tgbk_strand\n")
#for gb_record in SeqIO.parse(open(sys.argv[1],"r"), "genbank") :
#    genes = gb_record.features
##    print(genes)
#    for seq_feature in gb_record.features :
##        if seq_feature.type=="RNA":
#        m = re.search("[a-z]RNA|CDS",seq_feature.type)
#        if m:
##        output_handle.write(str(type(seq_feature.type)))
##        match_object = re.search(r'*RNA', str(seq_feature.type))
##        output_handle.write(match_object)
##        if (seq_feature.type=="CDS" or match_object):
#            if (seq_feature.qualifiers):
##                if ('pseudo' in seq_feature.qualifiers):
##                    continue
##                #gene is not a pseudogene
##                else:
#                end = seq_feature.location.end
#                start = seq_feature.location.start
#                start = start + 1
#                tmp = end + start
#                midpoint = int(tmp/2)
#                if seq_feature.strand == 1:
#    		#leading strand = 0
#                    cod_val = 0
#                else:
#    		#lagging strand = 1
#                    cod_val = 1
#                #dealing with the fact that there might not be a
#                #gene id present
#                if ('gene' in seq_feature.qualifiers):
#                    gene_id = seq_feature.qualifiers['gene'][0]
#                #gene is is missing
#                else:
#                    gene_id = 'NA'
#                #dealing with the fact that there might not be a
#                #locus tag present
#                if ('locus_tag' in seq_feature.qualifiers):
#                    locus_tag = seq_feature.qualifiers['locus_tag'][0]
#                #gene is is missing
#                else:
#                    locus_tag = 'NA'
#                #dealing with the fact that there might not be a
#                #old locus tag present
#                if ('old_locus_tag' in seq_feature.qualifiers):
#                    old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
#                #gene is is missing
#                else:
#                    old_locus_tag = 'NA'
#                #print out all the info
#                output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
#                       start,
#                       seq_feature.location.end,
#                       midpoint,
#                       gene_id,
#                       locus_tag,
#                       old_locus_tag,
#                       cod_val))

#REMOVE PSEUDOGENES
faa_filename = sys.argv[2]
output_handle = open(faa_filename, "w")
output_handle.write("gbk_start\tgbk_end\tgbk_midpoint\tgbk_gene_id\tgbk_locus_tag\tgbk_old_locus_tag\tgbk_strand\n")
for gb_record in SeqIO.parse(open(sys.argv[1],"r"), "genbank") :
    genes = gb_record.features
#    print(genes)
    for seq_feature in gb_record.features :
#        if seq_feature.type=="RNA":
        m = re.search("[a-z]RNA|CDS",seq_feature.type)
        if m:
#        output_handle.write(str(type(seq_feature.type)))
#        match_object = re.search(r'*RNA', str(seq_feature.type))
#        output_handle.write(match_object)
#        if (seq_feature.type=="CDS" or match_object):
            if (seq_feature.qualifiers):
                if ('pseudo' in seq_feature.qualifiers):
                    continue
                #gene is not a pseudogene
                else:
                    end = seq_feature.location.end
                    start = seq_feature.location.start
                    start = start + 1
                    tmp = end + start
                    midpoint = int(tmp/2)
                    if seq_feature.strand == 1:
        		#leading strand = 0
                        cod_val = 0
                    else:
        		#lagging strand = 1
                        cod_val = 1
                    #dealing with the fact that there might not be a
                    #gene id present
                    if ('gene' in seq_feature.qualifiers):
                        gene_id = seq_feature.qualifiers['gene'][0]
                    #gene is is missing
                    else:
                        gene_id = 'NA'
                    #dealing with the fact that there might not be a
                    #locus tag present
                    if ('locus_tag' in seq_feature.qualifiers):
                        locus_tag = seq_feature.qualifiers['locus_tag'][0]
                    #gene is is missing
                    else:
                        locus_tag = 'NA'
                    #dealing with the fact that there might not be a
                    #old locus tag present
                    if ('old_locus_tag' in seq_feature.qualifiers):
                        old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
                    #gene is is missing
                    else:
                        old_locus_tag = 'NA'
                #print out all the info
                output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                       start,
                       seq_feature.location.end,
                       midpoint,
                       gene_id,
                       locus_tag,
                       old_locus_tag,
                       cod_val))












































































output_handle.write("0\t0\t0\t0\t0\t0\t0\n")

output_handle.close()
