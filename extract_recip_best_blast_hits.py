#!/usr/bin/python
# Extract Reciprocal Best Blast Hits
#######################################
# Input: tabular blast searches for AvsB and BvsA

from Bio import SeqIO
from optparse import OptionParser
import pandas as pd
import sys, re
import numpy as np

#properly read scientific notation
pd.set_option('display.float_format', lambda x: '%.6e' % x)

#flags and options
parser = OptionParser()
parser.add_option("-a", "--blast_output_A", dest="blast_result1", help="tabular blast output of AvsB", metavar="FILE")
parser.add_option("-b", "--blast_output_B", dest="blast_result2", help="tabular blast output of BvsA", metavar="FILE")
parser.add_option("-o", "--out_file", dest="out_handle", help="output of reciprocal best blast hits for A and B", metavar="FILE")
parser.add_option("-c", "--proteome_file_A", dest="proteome_fasta1", help="fasta file with proteins from A", metavar="FILE")
parser.add_option("-d", "--proteome_file_B", dest="proteome_fasta2", help="fasta file with proteins from B", metavar="FILE")
parser.add_option("-e", "--gbk_file_A", dest="gbk_file_A", help="gbk file for A", metavar="FILE")
parser.add_option("-f", "--gbk_file_B", dest="gbk_file_B", help="gbk file for B", metavar="FILE")
parser.add_option("-g", "--alt_names_A", dest="alt_names_A", help="tab separated file with alternate blast gene names for A", metavar="FILE")
parser.add_option("-i", "--alt_names_B", dest="alt_names_B", help="tab separated file with alternate blast gene names for B", metavar="FILE")
(options, args) = parser.parse_args()

#standard blast output column names:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

#open first blast output (AvsB)
AvsB_qseqid = []
AvsB_sseqid = []
AvsB_evalue = []
AvsB_bitscore = []
with open(options.blast_result1) as tabular:
    for line in tabular:
        splitting = line.split("\t")
        name = splitting[0]
        sid = splitting[1]
        evalue = splitting[10]
#        a = re.compile("e")
#        evalue = a.sub("E",evalue)
        btscore = splitting[11]
        a = re.compile("\n")
        btscore = a.sub("",btscore)
#        a = re.compile("e")
#        btscore = a.sub("E",btscore)
        AvsB_qseqid.append(name)
        AvsB_sseqid.append(sid)
        AvsB_evalue.append(float(evalue))
        AvsB_bitscore.append(float(btscore))

#open second blast output (BvsA)
BvsA_qseqid = []
BvsA_sseqid = []
BvsA_evalue = []
BvsA_bitscore = []
with open(options.blast_result2) as tabular:
    for line in tabular:
        splitting = line.split("\t")
        name = splitting[0]
        sid = splitting[1]
        evalue = splitting[10]
#        a = re.compile("e")
#        evalue = a.sub("E",evalue)
        btscore = splitting[11]
        a = re.compile("\n")
        btscore = a.sub("",btscore)
#        a = re.compile("e")
#        btscore = a.sub("E",btscore)
        BvsA_qseqid.append(name)
        BvsA_sseqid.append(sid)
        BvsA_evalue.append(float(evalue))
        BvsA_bitscore.append(float(btscore))
		#  identity = float(splitting[2])
              #  if name in 
                #        hits.append(name)
#find single best hit(s) in each blast output
#make df of just the columns we want
AvsB_df = {'qseqid':AvsB_qseqid, 'sseqid':AvsB_sseqid, 'evalue':AvsB_evalue, 'bitscore':AvsB_bitscore}
AvsB_df = pd.DataFrame(AvsB_df, columns=['qseqid','sseqid','evalue','bitscore'])
BvsA_df = {'qseqid':BvsA_qseqid, 'sseqid':BvsA_sseqid, 'evalue':BvsA_evalue, 'bitscore':BvsA_bitscore}
BvsA_df = pd.DataFrame(BvsA_df, columns=['qseqid','sseqid','evalue','bitscore'])

#all rows for each unique qseqid
uniq_AvsB_qseqid = AvsB_df.qseqid.unique()
uniq_BvsA_qseqid = BvsA_df.qseqid.unique()
AvsB_bh = pd.DataFrame()
BvsA_bh = pd.DataFrame()
for i in uniq_AvsB_qseqid:
    #rows that match each qseqid
    tmp_df = AvsB_df.loc[AvsB_df['qseqid'] == i]
    #all possible best hits
    AvsB_bh = AvsB_bh.append(tmp_df.loc[tmp_df['evalue'] == tmp_df['evalue'].min()])
for i in uniq_BvsA_qseqid:
    #rows that match each qseqid
    tmp_df = BvsA_df.loc[BvsA_df['qseqid'] == i]
    #all possible best hits
    BvsA_bh = BvsA_bh.append(tmp_df.loc[tmp_df['evalue'] == tmp_df['evalue'].min()])
#get reciprocal best hits
##create dataframes from just the columns we want
AvsB_rbh = AvsB_bh[['qseqid','sseqid']]
##make sure to switch the columns around
BvsA_rbh = BvsA_bh[['sseqid','qseqid']]
BvsA_rbh = BvsA_rbh.rename(columns={'sseqid': 'qseqid', 'qseqid': 'sseqid'})
#reciprocal best blast hits
intersected_df = pd.merge(AvsB_rbh, BvsA_rbh, how='inner')

#read in gbk files
A_gbk_gene_id = []
A_gbk_locus_tag = []
A_gbk_old_locus_tag = []
A_alt_gbk_gene_id = []
A_alt_gbk_locus_tag = []
A_alt_gbk_old_locus_tag = []
A_gbk_gene_syn = []
A_alt_gbk_gene_syn = []
A_gbk_start = []
A_gbk_end = []
for gb_record in SeqIO.parse(options.gbk_file_A, "genbank"):
    genes = gb_record.features
    for seq_feature in gb_record.features :
        m = re.search("[a-z]RNA|CDS",seq_feature.type)
        if m:
            if (seq_feature.qualifiers):
                #dealing with the fact that there might not be a
                #gene id present
                if ('gene' in seq_feature.qualifiers):
                    gene_id = seq_feature.qualifiers['gene'][0]
                    if ('/' in gene_id):
                       spl_id = gene_id.split("/")
                       gene_id = spl_id[0]
                       alt_gene_id = spl_id[1]
                    else:   
                       gene_id = seq_feature.qualifiers['gene'][0]
                       alt_gene_id = 'NA'
                #gene is is missing
                else:
                    gene_id = 'NA'
                    alt_gene_id = 'NA'
                #dealing with the fact that there might not be a
                #gene synonym present
                if ('gene_synonym' in seq_feature.qualifiers):
                    gene_syn = seq_feature.qualifiers['gene_synonym'][0]
                    if ('/' in gene_syn):
                       spl_id = gene_syn.split("/")
                       gene_syn = spl_id[0]
                       alt_gene_syn = spl_id[1]
                    else:   
                       gene_syn = seq_feature.qualifiers['gene_synonym'][0]
                       alt_gene_syn = 'NA'
                #gene is is missing
                else:
                    gene_syn = 'NA'
                    alt_gene_syn = 'NA'
                #dealing with the fact that there might not be a
                #locus tag present
                if ('locus_tag' in seq_feature.qualifiers):
                    locus_tag = seq_feature.qualifiers['locus_tag'][0]
                    if ('/' in locus_tag):
                       spl_id = locus_tag.split("/")
                       locus_tag = spl_id[0]
                       alt_locus_tag = spl_id[1]
                    else:   
                       locus_tag = seq_feature.qualifiers['locus_tag'][0]
                       alt_locus_tag = 'NA'
                #gene is is missing
                else:
                    locus_tag = 'NA'
                    alt_locus_tag = 'NA'
                #dealing with the fact that there might not be a
                #old locus tag present
                if ('old_locus_tag' in seq_feature.qualifiers):
                    old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
                    if ('/' in old_locus_tag):
                       spl_id = old_locus_tag.split("/")
                       old_locus_tag = spl_id[0]
                       alt_old_locus_tag = spl_id[1]
                    else:   
                       old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
                       alt_old_locus_tag = 'NA'
                #gene is is missing
                else:
                    old_locus_tag = 'NA'
                    alt_old_locus_tag = 'NA'
                start = seq_feature.location.start
                end = seq_feature.location.end
                #print out all the info
                A_gbk_gene_id.append(gene_id)
                A_gbk_gene_syn.append(gene_syn)
                A_alt_gbk_gene_syn.append(gene_syn)
                A_gbk_locus_tag.append(locus_tag)
                A_gbk_old_locus_tag.append(old_locus_tag)
                A_alt_gbk_gene_id.append(alt_gene_id)
                A_alt_gbk_locus_tag.append(alt_locus_tag)
                A_alt_gbk_old_locus_tag.append(alt_old_locus_tag)
                A_gbk_start.append(start)
                A_gbk_end.append(end)
#create gbk dfs
A_gbk_df = {'gene_id':A_gbk_gene_id, 'locus_tag':A_gbk_locus_tag, 'old_locus_tag':A_gbk_old_locus_tag,'gene_syn':A_gbk_gene_syn,'alt_gene_id':A_alt_gbk_gene_id, 'alt_locus_tag':A_alt_gbk_locus_tag, 'alt_old_locus_tag':A_alt_gbk_old_locus_tag,'alt_gene_syn':A_alt_gbk_gene_syn,'start':A_gbk_start, 'end':A_gbk_end}
A_gbk_df = pd.DataFrame(A_gbk_df,
columns=['gene_id','locus_tag','old_locus_tag','gene_syn','alt_gene_id','alt_locus_tag','alt_old_locus_tag','gene_syn','start','end'])

B_gbk_gene_id = []
B_gbk_locus_tag = []
B_gbk_old_locus_tag = []
B_alt_gbk_gene_id = []
B_alt_gbk_locus_tag = []
B_alt_gbk_old_locus_tag = []
B_gbk_gene_syn = []
B_alt_gbk_gene_syn = []
B_gbk_start = []
B_gbk_end = []
for gb_record in SeqIO.parse(options.gbk_file_B, "genbank"):
    genes = gb_record.features
    for seq_feature in gb_record.features :
        m = re.search("[a-z]RNA|CDS",seq_feature.type)
        if m:
            if (seq_feature.qualifiers):
                #dealing with the fact that there might not be a
                #gene id present
                if ('gene' in seq_feature.qualifiers):
                    gene_id = seq_feature.qualifiers['gene'][0]
                    if ('/' in gene_id):
                       spl_id = gene_id.split("/")
                       gene_id = spl_id[0]
                       alt_gene_id = spl_id[1]
                    else:   
                       gene_id = seq_feature.qualifiers['gene'][0]
                       alt_gene_id = 'NA'
                #gene is is missing
                else:
                    gene_id = 'NA'
                    alt_gene_id = 'NA'
                #dealing with the fact that there might not be a
                #gene synonym present
                if ('gene_synonym' in seq_feature.qualifiers):
                    gene_syn = seq_feature.qualifiers['gene_synonym'][0]
                    if ('/' in gene_syn):
                       spl_id = gene_syn.split("/")
                       gene_syn = spl_id[0]
                       alt_gene_syn = spl_id[1]
                    else:   
                       gene_syn = seq_feature.qualifiers['gene_synonym'][0]
                       alt_gene_syn = 'NA'
                #gene is is missing
                else:
                    gene_syn = 'NA'
                    alt_gene_syn = 'NA'
                #dealing with the fact that there might not be a
                #locus tag present
                if ('locus_tag' in seq_feature.qualifiers):
                    locus_tag = seq_feature.qualifiers['locus_tag'][0]
                    if ('/' in locus_tag):
                       spl_id = locus_tag.split("/")
                       locus_tag = spl_id[0]
                       alt_locus_tag = spl_id[1]
                    else:   
                       locus_tag = seq_feature.qualifiers['locus_tag'][0]
                       alt_locus_tag = 'NA'
                #gene is is missing
                else:
                    locus_tag = 'NA'
                    alt_locus_tag = 'NA'
                #dealing with the fact that there might not be a
                #old locus tag present
                if ('old_locus_tag' in seq_feature.qualifiers):
                    old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
                    if ('/' in old_locus_tag):
                       spl_id = old_locus_tag.split("/")
                       old_locus_tag = spl_id[0]
                       alt_old_locus_tag = spl_id[1]
                    else:   
                       old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
                       alt_old_locus_tag = 'NA'
                #gene is is missing
                else:
                    old_locus_tag = 'NA'
                    alt_old_locus_tag = 'NA'
                start = seq_feature.location.start
                end = seq_feature.location.end
                #print out all the info
                B_gbk_gene_id.append(gene_id)
                B_gbk_gene_syn.append(gene_syn)
                B_alt_gbk_gene_syn.append(gene_syn)
                B_gbk_locus_tag.append(locus_tag)
                B_gbk_old_locus_tag.append(old_locus_tag)
                B_alt_gbk_gene_id.append(alt_gene_id)
                B_alt_gbk_locus_tag.append(alt_locus_tag)
                B_alt_gbk_old_locus_tag.append(alt_old_locus_tag)
                B_gbk_start.append(start)
                B_gbk_end.append(end)
#create gbk dfs
B_gbk_df = {'gene_id':B_gbk_gene_id, 'locus_tag':B_gbk_locus_tag, 'old_locus_tag':B_gbk_old_locus_tag,'gene_syn':B_gbk_gene_syn,'alt_gene_id':B_alt_gbk_gene_id, 'alt_locus_tag':B_alt_gbk_locus_tag, 'alt_old_locus_tag':B_alt_gbk_old_locus_tag,'alt_gene_syn':B_alt_gbk_gene_syn,'start':B_gbk_start, 'end':B_gbk_end}
B_gbk_df = pd.DataFrame(B_gbk_df,
columns=['gene_id','locus_tag','old_locus_tag','gene_syn','alt_gene_id','alt_locus_tag','alt_old_locus_tag','gene_syn','start','end'])


#read in alt blast gene name files
#A
A_alt_names = pd.read_csv(options.alt_names_A, sep=" ", names=['b_name','b_locus'])
#make dictionary
A_alt_names_dic = A_alt_names.groupby('b_name')['b_locus'].apply(list).to_dict()
#B
B_alt_names = pd.read_csv(options.alt_names_B, sep=" ", names=['b_name','b_locus'])
#make dictionary
B_alt_names_dic = B_alt_names.groupby('b_name')['b_locus'].apply(list).to_dict()




##go through each unique value in qseqid column and find the highest value in evalue column
##Best Hit
#AvsB_bh = (AvsB_df.sort_values(['qseqid', 'evalue'], ascending=[True, True])
#             .drop_duplicates(['qseqid']).reset_index(drop=True)
#          )
#BvsA_bh = (BvsA_df.sort_values(['qseqid', 'evalue'], ascending=[True, True])
#             .drop_duplicates(['qseqid']).reset_index(drop=True)
#          )
###create dataframes from just the columns we want
#AvsB_rbh = AvsB_bh[['qseqid','sseqid']]
###make sure to switch the columns around
#BvsA_rbh = BvsA_bh[['sseqid','qseqid']]
##reciprocal best blast hits
#intersected_df = pd.merge(AvsB_rbh, BvsA_rbh, how='inner')

#getting additional info from fast files
#read in sequence data and store info
A_seqs = []
A_seq_names = []
A_blast_id = []
for record in SeqIO.parse(options.proteome_fasta1, "fasta"):
        rec_name = record.description
        splitting = rec_name.split("=")
        rec_name = splitting[3]
        a = re.compile(" PE")
        rec_name = a.sub("",rec_name)
        A_seq_names.append(rec_name)
        A_blast_id.append(record.id)
        A_seqs.append(record)
#put into a df
A_df = {'A_blast_id':A_blast_id, 'A_seq_name':A_seq_names}
A_df = pd.DataFrame(A_df, columns=['A_blast_id','A_seq_name'])

B_seqs = []
B_seq_names = []
B_blast_id = []
for record in SeqIO.parse(options.proteome_fasta2, "fasta"):
        rec_name = record.description
        splitting = rec_name.split("=")
        rec_name = splitting[3]
        a = re.compile(" PE")
        rec_name = a.sub("",rec_name)
        B_seq_names.append(rec_name)
        B_blast_id.append(record.id)
        B_seqs.append(record)
#put into a df
B_df = {'B_blast_id':B_blast_id, 'B_seq_name':B_seq_names}
B_df = pd.DataFrame(B_df, columns=['B_blast_id','B_seq_name'])

#add above blast info to df
#match blast to gene name
#make dictionary of A_df for faster searching
A_d = {}
for a, b in zip(A_df['A_blast_id'], A_df['A_seq_name']):
    if isinstance(a, list):
        for c in a:
            A_d[c] = b
    else:
        A_d[a] = b
intersected_df['A_gene'] = intersected_df['qseqid'].map(A_d)
#make dictionary of B_df for faster searching
B_d = {}
for a, b in zip(B_df['B_blast_id'], B_df['B_seq_name']):
    if isinstance(a, list):
        for c in a:
            B_d[c] = b
    else:
        B_d[a] = b
intersected_df['B_gene'] = intersected_df['sseqid'].map(B_d)

#search gbk file and add gene starts and ends to df
##make dictionary of A_df for faster searching
A_start_d = {}
for a, b in zip(A_gbk_df['gene_id'], A_gbk_df['start']):
    if isinstance(a, list):
        for c in a:
            A_start_d[c] = b
    else:
        A_start_d[a] = b
intersected_df['A_start'] = intersected_df['A_gene'].map(A_start_d)
##make dictionary of B_df for faster searching
B_start_d = {}
for a, b in zip(B_gbk_df['gene_id'], B_gbk_df['start']):
    if isinstance(a, list):
        for c in a:
            B_start_d[c] = b
    else:
        B_start_d[a] = b
intersected_df['B_start'] = intersected_df['B_gene'].map(B_start_d)
##make dictionary of A_df for faster searching
A_end_d = {}
for a, b in zip(A_gbk_df['gene_id'], A_gbk_df['end']):
    if isinstance(a, list):
        for c in a:
            A_end_d[c] = b
    else:
        A_end_d[a] = b
intersected_df['A_end'] = intersected_df['A_gene'].map(A_end_d)
##make dictionary of B_df for faster searching
B_end_d = {}
for a, b in zip(B_gbk_df['gene_id'], B_gbk_df['end']):
    if isinstance(a, list):
        for c in a:
            B_end_d[c] = b
    else:
        B_end_d[a] = b
intersected_df['B_end'] = intersected_df['B_gene'].map(B_end_d)

#use synteny to find single rbbh for each  duplicated qseqid
for i in intersected_df.qseqid.unique():
    tmp_df = intersected_df.loc[intersected_df['qseqid'] == i]
    num_rows = len(tmp_df.qseqid)
    if (num_rows > 1):
        #multiple rbbhs
#        print(tmp_df)
#        print(num_rows)
        #checking for na values (gene was not found in gbk file, mostly instertion seqs)   
        if (tmp_df['A_start'].isnull().values.any() or tmp_df['B_start'].isnull().values.any() or tmp_df['A_end'].isnull().values.any() or tmp_df['B_end'].isnull().values.any()):
#            print("ALL NA")
            intersected_df = intersected_df.drop(intersected_df.index[intersected_df['qseqid'] == i].tolist())
        else:
#            print("no NA")
            As_Bs_diff = []
            As_Be_diff = []
            Ae_Bs_diff = []
            Ae_Be_diff = []
            for x in range(0,num_rows):
                tmp_diff = abs(tmp_df.iloc[x]['A_start'] - tmp_df.iloc[x]['B_start'])
                As_Bs_diff.append(tmp_diff)  
                tmp_diff = abs(tmp_df.iloc[x]['A_start'] - tmp_df.iloc[x]['B_end'])
                As_Be_diff.append(tmp_diff)  
                tmp_diff = abs(tmp_df.iloc[x]['A_end'] - tmp_df.iloc[x]['B_start'])
                Ae_Bs_diff.append(tmp_diff)  
                tmp_diff = abs(tmp_df.iloc[x]['A_end'] - tmp_df.iloc[x]['B_end'])
                Ae_Be_diff.append(tmp_diff)  
            tmp_df['As_Bs_diff'] = np.array(As_Bs_diff)
            tmp_df['As_Be_diff'] = np.array(As_Be_diff)
            tmp_df['Ae_Bs_diff'] = np.array(Ae_Bs_diff)
            tmp_df['Ae_Be_diff'] = np.array(Ae_Be_diff)
            min_locs = [tmp_df.As_Bs_diff.idxmin(),tmp_df.As_Be_diff.idxmin(),tmp_df.Ae_Bs_diff.idxmin(),tmp_df.Ae_Be_diff.idxmin()]
            min_vals = [tmp_df.As_Bs_diff.min(),tmp_df.As_Be_diff.min(),tmp_df.Ae_Bs_diff.min(),tmp_df.Ae_Be_diff.min()]
            all_locs = list(tmp_df.index.values)
#            print('all_locs')
#            print(all_locs)
#            print("min")
#            print(min_vals)
#            print(min_locs)
            mm_loc = np.argmin(min_vals)
#            print(min_locs[mm_loc])
            all_locs.remove(min_locs[mm_loc])
#            print('all_locs')
#            print(all_locs)
            intersected_df = intersected_df.drop(all_locs)
            
#        print(tmp_df)
                
        
    
#add new cols with the alt blast names
intersected_df['A_alt_name'] = intersected_df['A_gene'].map(A_alt_names_dic)
intersected_df['B_alt_name'] = intersected_df['B_gene'].map(B_alt_names_dic)



out_file = open(options.out_handle, "w")
with open(options.out_handle, "w") as f:
    intersected_df.to_csv(f, sep='\t', index=False)
