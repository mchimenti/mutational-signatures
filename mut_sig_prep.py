#!/usr/bin/env python 
# -*- coding: utf-8 -*-

## Author: Michael Chimenti
## Date:   Oct 21, 2015

'''
NOTE: 
If you are using hg38 instead of hg19, just edit line 101 to use the hg38 dictionary.

Please edit line 46 (Entrez Email) to reflect your own email address.  

DESCRIPTION: 
A script to fetch SNPs from HG19 or HG38 reference genome and retrieve the trimer -1/+1 bases 
on either side of the SNP.  These are then sorted and tabulated. Outputs for creating
a MatLab .mat file from csv vectors are written to file as .csv files.

The script expects one column of your input dataframe/csv file to be called "sample"
containing sample names from which each SNP is originating and one to be called 
"VAR" containing the variant definitions in the common format:
    
chr1:159482:A>G

This code is meant to prepare data for input into the MatLab framework as part of 
the work described in Alexandrov et. al., Nature 500, 2013 and Cell Rep. 2013 Jan 31; 3(1): 246–259.   

In order for the MatLab code to work, you should process at least 10,000 SNPs. Ideally,
you'd have 100,000 - 300,000 SNPs across at least 10 samples.  This code needs enough
SNPs so that each of the 96 possible variations appears at least once in one sample. 

Output:

samples.csv : List of sample names
types.csv   : Vector of 96 substitution types (C>A, C>T, etc...matches the samples-by-counts matrix)
subtypes.csv: Vector of 96 substitution subtypes (ACA, ACT, ACG, etc...)
sample-by-counts.csv : (96 x # of samples) matrix with counts for all types/subtypes. 
'''

import pandas as pd
import sys
import csv

from Bio import Entrez, SeqIO       #Biopython version 1.65

#Global definitions
Entrez.email = "michael-chimenti@uiowa.edu"    


hg19_chrom = { 'chr1':'224589800', 'chr2':'224589811', 'chr3':'224589815',
               'chr4':'224589816', 'chr5':'224589817', 'chr6':'224589818', 
               'chr7':'224589819', 'chr8':'224589820', 'chr9':'224589821',
               'chr10':'224589801', 'chr11':'224589802', 'chr12':'224589803',
               'chr13':'224589804', 'chr14':'224589805', 'chr15':'224589806',
               'chr16':'224589807', 'chr17':'224589808', 'chr18':'224589809',
               'chr19':'224589810', 'chr20':'224589812', 'chr21':'224589813',
               'chr22':'224589814', 'chrX':'224589822', 'chrY':'224589823'}
               
#GI IDs from http://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.31/              
hg38_chrom = {'chr1':'568336023', 'chr2':'568336022', 'chr3':'568336021',
               'chr4':'568336020', 'chr5':'568336019', 'chr6':'568336018', 
               'chr7':'568336017', 'chr8':'568336016', 'chr9':'568336015',
               'chr10':'568336014', 'chr11':'568336013', 'chr12':'568336012',
               'chr13':'568336011', 'chr14':'568336010', 'chr15':'568336009',
               'chr16':'568336008', 'chr17':'568336007', 'chr18':'568336006',
               'chr19':'568336005', 'chr20':'568336004', 'chr21':'568336003',
               'chr22':'568336002', 'chrX':'568336001', 'chrY':'568336000'}

 
def sample_table(filename):
    '''read in a csv file with a header in the first row'''
    try:
        snp_table = pd.read_csv(filename, header=0)
    except IOError:
        print "File", filename, "cannot be found."
        raise
        
    return snp_table    
    

def group_table(table):
    '''group table on sample ID, return groups object'''
    try:
        grouped = table.groupby('sample')
    except KeyError:
        print "Could not find column called 'sample' in input"
        raise
    
    return grouped   


def get_snp_list(data):
    '''get SNPs as list from dataframe'''
    try:
        snplist = data.VAR.values.tolist()
    except ValueError:
        print "Please title the SNP column as 'var' in your input"
        raise
    
    return snplist
        
        
def split_vars(snp):  
    '''Split a variant field into chrom number, position, and variant (i.e., A>C)'''     
    snp_split = snp.split(':')
    chrom = snp_split[0]
    pos = snp_split[1]
    var = snp_split[2]  
    return chrom, pos, var      


def get_record(chrom, pos):               
    '''Fetch record, and minus/plus one base from Entrez'''
    try:
        handle = Entrez.efetch(db="nucleotide", 
                    id=hg19_chrom[chrom], 
                    rettype="fasta", 
                    strand=1, 
                    seq_start=int(pos) - 1, 
                    seq_stop=int(pos) + 1)
                       
        record = SeqIO.read(handle, "fasta")
        handle.close()
    except BaseException:
        return None
    
    return record
     
           
def format_trimer_table(table):
    '''Get the REF/ALT bases and create the trimer for downstream calculation'''
    trimer_tab = table.loc[:,["sample","VAR","minus_one","plus_one"]]
    trimer_tab["ref"] = trimer_tab.VAR.apply(lambda r: r.split(':')[2].split('>')[0])  
    trimer_tab["alt"] = trimer_tab.VAR.apply(lambda r: r.split(':')[2].split('>')[1])  
    trimer_tab["trimer"] = trimer_tab.minus_one + trimer_tab.ref + trimer_tab.plus_one 
    #filter the table
    trimer_tab = trimer_tab[trimer_tab.trimer.apply(lambda r: len(r) == 3)]  #drop longer indels
    trimer_tab = trimer_tab[trimer_tab.alt != '-']  # drop single nucleotide deletions
    #create code field
    trimer_tab["code"] = trimer_tab.VAR.apply(lambda r: r.split(':')[2]) + ':' + trimer_tab.trimer
    #MatLab framework only cares about C>X and T>X mutations
    trimer_tab = trimer_tab.loc[trimer_tab.code.apply(lambda r: r.split(':')[0] not in
                                                ['G>A', 'G>T', 'G>C', 'A>T','A>G','A>C'])]  
    return trimer_tab    
    
       
def write_data_for_matlab(df, colnames): 
    '''calculate the inputs for the matlab framework'''        
    df.columns = colnames
    df.fillna(0, inplace=True)

    codon_type = []
    codon_subtype = []
    for i in df.index:
        t = i.split(':')[0]
        s = i.split(':')[1]
        codon_type.append(t)
        codon_subtype.append(s)
    
    samples = colnames

    df.to_csv("samples_by_counts.csv", sep=",",header=False,index=False)

    with open('samples.csv', 'wb') as myfile:
        wr = csv.writer(myfile)
        wr.writerow(samples)

    with open('types.csv', 'wb') as myfile:
        wr = csv.writer(myfile)
        wr.writerow(codon_type)
    
    with open('subtypes.csv', 'wb') as myfile:
        wr = csv.writer(myfile)
        wr.writerow(codon_subtype)             

def main():
    if len(sys.argv) == 2:
        filename = str(sys.argv[1])
        
    table = sample_table(filename)
    grouped_table = group_table(table)
    
    colnames = []
    df = pd.DataFrame()
    for name, data in grouped_table:
        snp_list = get_snp_list(data)
        snp_flank_left = []
        snp_flank_right = []
        for snp in snp_list:
            chrom, pos, var = split_vars(snp)
            if len(var) == 3:
                record = get_record(chrom, pos)
                print "Got record:", record.seq
                if record.seq[1] == var.split('>')[0]: #compare to reference base
                    snp_flank_left.append(str(record.seq[0]))
                    snp_flank_right.append(str(record.seq[2]))
                else:
                    snp_flank_left.append("NA")
                    snp_flank_right.append("NA")
            
            else:
                snp_flank_left.append("NA")
                snp_flank_right.append("NA")
        
        data["minus_one"] = snp_flank_left
        data["plus_one"] = snp_flank_right
        #write out each snp table as they are completed
        data.to_csv(name + '.csv', index=False)
        
        trimer_tab = format_trimer_table(data)
        counts = trimer_tab.code.value_counts()  #should be 96 codes total
        colnames.append(name)
        df = pd.concat([df,counts], axis=1)

    write_data_for_matlab(df, colnames)

        
if __name__ == '__main__':
    main()
    
