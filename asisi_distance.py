#!/usr/bin/env python3
#this is run through the pipeline

import sys
import pandas as pd
import py2bit
import re

ID = sys.argv[1]
MOTIF = 'GCGATCGC'
filename = '/main/results/' + ID + '.formatted.txt'

# read file
with open(filename) as f:
    lines = f.readlines()

hg19 = py2bit.open('/main/hg19.2bit')
chroms_temp = hg19.chroms()

# get indices of motif in all chromosomes
keep = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
chroms = {}
Q=re.compile(MOTIF)
for chrom in chroms_temp.keys():
    if chrom in keep:
        print(chrom)
        seq = hg19.sequence(chrom).upper()
        chroms[chrom] = [item.start(0)+5 for item in Q.finditer(seq)]
hg19.close()

# # write bed file
# with open('/main/results/AsisI_sites.bed3', 'w') as file:
#     for chrom in chroms:
#         for site in chroms[chrom]:
#             file.write(chrom+'\t'+str(site)+'\t'+str(int(site)+8)+'\n')

# find distances
out = []
for line in lines:
    line = line.split('\t')
    if line[0] in keep and len(chroms[line[0]]) > 0:
        distance = min([abs(int(line[1])-i) for i in chroms[line[0]]])
        site = int(line[1])-distance
        if distance+int(line[1]) in chroms[line[0]]:
            site = distance+int(line[1])
        out.append([line[3].replace('\n',''),line[0],line[1],line[2],str(site),str(distance)])

# export to csv
df = pd.DataFrame(out, columns=['Seq ID','INS Chrom','INS Start','INS End','nearest_site','Distance'])
df.to_csv('/main/results/'+ID+'.AsiSIdistance.csv',index=False)