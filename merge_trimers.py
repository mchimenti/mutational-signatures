

import pandas as pd
import csv

trimer_tab = pd.read_csv('~/iihg/schularick_mut_sig/snp_tab_cat.csv', header=False)

trimer_tab["ref"] = trimer_tab.VAR.apply(lambda r: r.split(':')[2].split('>')[0])  #get ref base
trimer_tab["alt"] = trimer_tab.VAR.apply(lambda r: r.split(':')[2].split('>')[1])  #get the alt base
trimer_tab["trimer"] = trimer_tab.minus_one + trimer_tab.ref + trimer_tab.plus_one  #create the trimer

trimer_tab = trimer_tab[trimer_tab.trimer.isnull() == False]   #drop longer indels, retain only SNPs
trimer_tab = trimer_tab[trimer_tab.alt != '-']  # drop single nucleotide deletions

trimer_tab["code"] = trimer_tab.VAR.apply(lambda r: r.split(':')[2]) + ':' + trimer_tab.trimer
trimer_tab["mut"] = trimer_tab.VAR.apply(lambda r: r.split(':')[2])

#full_tab = pd.read_csv('~/iihg/schularick_mut_sig/consol_dedup_AY3.csv', header=False)
#merged = full_tab.merge(trimer_tab, on=['Tumor','VAR'])
#merged = merged.drop(merged.columns[2:70], axis =1)  
#merged.to_csv('file', sep=',')

#The MatLab code considers only C>A,G,T and T>A,C,G variants
trimer_tab = trimer_tab[trimer_tab.code.apply(lambda r: r.split(':')[0] not in
                                                ['G>A', 'G>T', 'G>C', 'A>T','A>G','A>C'])] 
#reorder the dataframe
cols = ['Tumor','VAR','minus_one','ref','alt','plus_one','mut','trimer','code']
trimer_tab = trimer_tab[cols]

trimer_tab_grouped = trimer_tab.groupby("Tumor")

#for each tumor sample, get value counts at "code"; should be 96 codes total
colnames = []
df = pd.DataFrame()
for name, data in trimer_tab_grouped:
    print name
    counts = data.code.value_counts()
    colnames.append(name)
    df = pd.concat([df,counts], axis=1)
    
df.columns = colnames

#need to replace NAs with zeros in df

df.fillna(0, inplace=True)

codon_type = []
codon_subtype = []
for i in df.index:
    t = i.split(':')[0]
    s = i.split(':')[1]
    codon_type.append(t)
    codon_subtype.append(s)
    
samples = df.columns

df.to_csv("samples_by_counts.csv", sep=",",header=False,index=False)

with open('samples.csv', 'wb') as myfile:
    wr = csv.writer(myfile) #quoting=csv.QUOTE_ALL)
    wr.writerow(samples)

with open('types.csv', 'wb') as myfile:
    wr = csv.writer(myfile)#, quoting=csv.QUOTE_ALL)
    wr.writerow(codon_type)
    
with open('subtypes.csv', 'wb') as myfile:
    wr = csv.writer(myfile)# quoting=csv.QUOTE_ALL)
    wr.writerow(codon_subtype)