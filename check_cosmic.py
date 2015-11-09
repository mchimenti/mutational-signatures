#!/usr/bin/env python 

'''Read in Tumor SNP table, check if gene IDs are in COSMIC list of cancer
census genes.  Keep those rows that match.'''

import pandas as pd
import matplotlib.pyplot as plt

snp_tab = pd.read_csv('~/iihg/schularick_mut_sig/consolidated-dedup.csv', header=0)
cosmic = pd.read_csv('~/iihg/schularick_mut_sig/cosmic_grch38_cancer_gene_census.csv', header=0)

snp_tab_cosmic = snp_tab[snp_tab.GENE.isin(cosmic["Gene Symbol"])]
#snp_tab_cosmic.to_csv('snp_tab_gene_in_cosmic.csv', sep=',',index=False, header=True)

snp_tab_cosmic_counts = snp_tab_cosmic.GENE.value_counts()
snp_tab_cosmic_counts.plot(kind='barh', ylim=(0,30))

snp_tab_cosmic_grouped = snp_tab_cosmic.groupby("Tumor")


fig, axs = plt.subplots(6, 5)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.5, hspace=0.5)

for (i,j), ax in zip(snp_tab_cosmic_grouped, axs.flat):
    j.GENE.value_counts().plot(kind='barh',ax=ax,ylim=(0,10))

    
