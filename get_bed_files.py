
import pandas as pd
  
def get_bed(data):
    data["chrom"] = data.VAR.apply(lambda r: r.split(":")[0])
    data["chromStart"] = data.loc[:,"pos"]
    data["chromEnd"] = data.pos.apply(lambda r: r+1)
    bedfile = data.loc[:,("chrom","chromStart","chromEnd")]
    return bedfile
    

snp_tab = pd.read_csv('~/iihg/schularick_mut_sig/raw_data/consolidated-dedup.csv', usecols=[0,1,59], header=0)

snp_tab["pos"] = snp_tab.VAR.apply(lambda r: r.split(':')[1])
snp_tab.pos = snp_tab.pos.astype('int64')

snp_tab["type"] = snp_tab.VAR.apply(lambda r: r.split(':')[2])

#drop long indels
snp_tab = snp_tab.loc[snp_tab.type.apply(lambda r: len(r) == 3)]
#drop single indels
snp_tab = snp_tab.loc[snp_tab.type.apply(lambda r: r.find("-") == -1)]
#only take C>X and T>X SNPs
snp_tab = snp_tab[snp_tab.type.apply(lambda r: r not in
                                                ['G>A', 'G>T', 'G>C', 'A>T','A>G','A>C'])] 

snp_tab_grouped = snp_tab.groupby("Tumor")

for name, data in snp_tab_grouped:
    bedfile = get_bed(data)
    bedfile.to_csv(name + ".bed", sep = '\t', index=False, header=False)