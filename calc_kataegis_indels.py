#!/usr/bin/python

"""A script to calculate 'kataegis' (localized hypermutation) in cancer. Input
is list of SNPs with the format: chr1:159482:A>G.  Script assumes multiple samples
are present in the table and will group on the sample column.  The output is a 
plot called a 'rainfall' plot that will show regions of hypermutation."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import datetime
#import matplotlib.colors as colors
#import matplotlib.cm as cmx



matplotlib.style.use('ggplot')



snp_tab = pd.read_csv('~/iihg/schularick_mut_sig/raw_data/consolidated-dedup.csv', usecols=[0,1,59], header=0)

snp_tab["pos"] = snp_tab.VAR.apply(lambda r: r.split(':')[1])
snp_tab.pos = snp_tab.pos.astype('int64')
snp_tab["type"] = snp_tab.VAR.apply(lambda r: r.split(':')[2])

#drop long indels
#snp_tab = snp_tab.loc[snp_tab.type.apply(lambda r: len(r) == 3)]

#drop single indels
#snp_tab = snp_tab.loc[snp_tab.type.apply(lambda r: r.find("-") == -1)]

#just get T>x and C>x variants
snp_tab = snp_tab[snp_tab.type.apply(lambda r: r not in
                                                ['G>A', 'G>T', 'G>C', 'A>T','A>G','A>C'])] 
snp_tab["dist"] = snp_tab.pos.diff()

def get_colors(r):
    color_type = {"T>C":2,"T>A":3,"T>G":4,"C>T":5,"C>A":6,"C>G":7}
    if r in color_type.keys():
        return color_type[r]
    else:
       return 1
       
snp_tab["colors"] = snp_tab.type.apply(lambda r: get_colors(r))

snp_tab.fillna(value=0)
snp_tab_grouped = snp_tab.groupby("Tumor")

fig, axs = plt.subplots(6, 5)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.5, hspace=0.5)

#plot in one fig
#for (i,j), ax in zip(snp_tab_grouped, axs.flat):
#    j.reset_index(inplace=True)
#    j["number"] = j.index.tolist()
#    j.plot(kind='scatter',x='number',y="dist", ylim=(1,10000000), xlim=(0,4000),logy=True, title=i, ax=ax)

#plot sep figs
for (i,j), ax in zip(snp_tab_grouped, axs.flat):
    j.reset_index(inplace=True)
    j["number"] = j.index.tolist()
    j.plot(kind='scatter',x='number',y="dist", ylim=(1,10000000), xlim=(0,len(j)),
            logy=True, c = 'colors', cmap = 'spectral', title=i)
    
    
    plt.ylabel('Intermutational Distance (BP)', fontsize=20)
    j["time"] = j.index
    j.time = pd.to_datetime(j.time)
    j.set_index(j.time, inplace=True)
    pd.rolling_median(j.dist, 70, min_periods=5).plot(style='k')
   
    
    #plt.savefig(i+'.png')