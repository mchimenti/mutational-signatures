#!/usr/bin/python

"""A script to compare the cosine similarity of the computed signatures from the
raw data to the signatures in the COSMIC database.  Cosine similarity scores will
be ranked from lowest to highest for each signature.

ref_sigs = COSMIC database reference signatures
exp_sigs = Signatures calculated by MatLab framework from raw SNP data
"""

import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.style.use('ggplot')

ref_sigs = pd.read_csv("~/iihg/schularick_mut_sig/cosmic/signatures_probabilities.csv")

exp_sigs = pd.read_csv("~/iihg/schularick_mut_sig/tables/100iter_4_signatures_probabilities.csv")
exp_sigs.Trimer = exp_sigs.Trimer.apply(lambda r: r.strip("'"))
exp_sigs.Variant = exp_sigs.Variant.apply(lambda r: r.strip("'"))

#make sure to compare apples-to-apples
assert exp_sigs.Trimer.values.tolist() == ref_sigs.Trinucleotide.values.tolist()

ref_sigs_sub = ref_sigs.iloc[:,3:30]  
exp_sigs_sub = exp_sigs.iloc[:,2:6]

clist = []
rownames = []
colnames = []
df = pd.DataFrame()
for expsig in exp_sigs_sub:
    clist = []
    rownames = []
    for refsig in ref_sigs_sub:
        cos_sim = cosine_similarity(exp_sigs_sub[expsig], ref_sigs_sub[refsig])
        clist.append(cos_sim.item())
        rownames.append(refsig)
    
    colnames.append(expsig)
    df = pd.concat([df, pd.Series(clist)], axis=1)

colnames.append("COSMIC")

df = pd.concat([df, pd.Series(rownames)], axis=1)
df.columns = colnames

plt.figure()
df.plot(kind='bar',x = "COSMIC")
plt.ylabel("cosine similiarity (1 is identical)")   
plt.title("Cosine Similarity for Four Mutational Signatures to COSMIC Ref Signatures") 
plt.show()
    




