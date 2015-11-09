#!/usr/bin/env python 

## Author: Michael Chimenti
## Date:   Oct 21, 2015

'''A script to fetch SNPs from HG19 reference genome and retrieve the -1/+1 bases 
on either side of the SNP'''

import pandas as pd
from Bio import Entrez, SeqIO       #"Bio" refers to the biopython version 1.65
import datetime


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


#Script 
colnames = ['Tumor','VAR']
snp_tab = pd.read_csv('~/iihg/schularick_mut_sig/consolidated-dedup.csv', usecols=[0,1], header=False)
snp_tab.columns = colnames

snp_tab_grouped = snp_tab.groupby('Tumor')

processed = 'BE3'
for name, data in snp_tab_grouped:
    if name in processed:
        print "Processing tumor group:", name
        snplist = data.VAR.values.tolist()

        start_time = datetime.datetime.now()
        snp_flank_left = []
        snp_flank_right = []
        
        for i, snp in enumerate(snplist):
            snp_split = snp.split(':')
            if  len(snp_split[2]) == 3:
                chrom = snp_split[0]
                pos = snp_split[1]
                
                try:
                    handle = Entrez.efetch(db="nucleotide", 
                       id=hg19_chrom[chrom], 
                       rettype="fasta", 
                       strand=1, 
                       seq_start=int(pos) - 1, 
                       seq_stop=int(pos) + 1)
                       
                    record = SeqIO.read(handle, "fasta")
                    handle.close()
                
                    if record.seq[1] == snp_split[2].split('>')[0]:        #compare to ref base
                        print "record: ", i, "chrom: ",chrom, "pos: ", pos, "trimer: ", record.seq
                        snp_flank_left.append(str(record.seq[0]))
                        snp_flank_right.append(str(record.seq[2]))
                    else:
                        print i, "NA"
                        snp_flank_left.append("NA")
                        snp_flank_right.append("NA")
                    
                except BaseException as e:   #catch any exception and move on
                    snp_flank_left.append("NA")
                    snp_flank_right.append("NA")
                    print e
                    continue
            else:
                print i, "NA"
        	snp_flank_left.append("NA")
                snp_flank_right.append("NA")
          

        data["minus_one"] = snp_flank_left
        data["plus_one"] = snp_flank_right


        time = datetime.datetime.now() - start_time
        print "Elapsed time for this group: ", time 
        print "Writing SNP table to CSV File."

        data.to_csv('/Users/mchimenti/iihg/schularick_mut_sig/snp_tab_' + name + '.csv', index=False)

