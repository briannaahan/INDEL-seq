#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import math

ID = sys.argv[1]
filename1 = '/main/results/'+ID+'.annotation.csv'
filename2 = '/main/results/'+ID+'.INS.fasta.out'

with open(filename1) as f:
    annotations = f.readlines()
with open(filename2) as f:
    repeats = f.readlines()

annotations = annotations[1:]
repeats = repeats[3:]

pie = {'Promoter': 0, 'Genebody': 0, 'Distal Intergenic': 0, 'Downstream': 0}
repeatpie = {'Promoter': 0, 'Genebody': 0, 'Distal Intergenic': 0, 'Downstream': 0}
noreppie = {'Promoter': 0, 'Genebody': 0, 'Distal Intergenic': 0, 'Downstream': 0}

ids = []
for line in repeats:
    ids.append(line.split()[4])

for line in annotations:
    if 'Promoter' in line:
        key = 'Promoter'
    elif 'UTR' in line:
        key = 'Genebody'
    elif 'Exon' in line:
        key = 'Genebody'
    elif 'Intron' in line:
        key = 'Genebody'
    elif 'Distal Intergenic' in line:
        key = 'Distal Intergenic'
    elif 'Downstream' in line:
        key = 'Downstream'
    else:
        print('brrrruuhhhhhhhhhhhh')
        print(line)
        key = None
    if line.split('\t')[5] in ids:
        if key: repeatpie[key] += 1
    else: 
        if key: noreppie[key] += 1
    if key: pie[key] += 1

total_samples = pie["Promoter"] + pie["Genebody"] + pie["Distal Intergenic"] + pie["Downstream"]
repeat_samples = repeatpie["Promoter"] + repeatpie["Genebody"] + repeatpie["Distal Intergenic"] + repeatpie["Downstream"]
nonrepeat_samples = noreppie["Promoter"] + noreppie["Genebody"] + noreppie["Distal Intergenic"] + noreppie["Downstream"]

fig, (ax1, ax2, ax3) = plt.subplots(1,3)
fig.set_size_inches(12, 5, forward=True)
ax1.pie(pie.values(), labels=pie.keys(), autopct='%1.1f%%')
ax1.set_title('All')
ax2.pie(noreppie.values(), labels=noreppie.keys(), autopct='%1.1f%%')
ax2.set_title('Non-Repeats (' + str(round(nonrepeat_samples/total_samples*100,2)) + '%)')
ax3.pie(repeatpie.values(), labels=repeatpie.keys(), autopct='%1.1f%%')
ax3.set_title('Repeats (' + str(round(repeat_samples/total_samples*100,2)) + '%)')
fig.suptitle('Gene Annotations', fontsize=16)
fig.tight_layout()

plt.savefig('/main/figures/'+ID+'pie.png', dpi=500)

all = []
non_repeats = []
repeats = []

for line in annotations:
    if line.split('\t')[5] in ids:
        repeats.append(int(line.split('\t')[14]))
    else:
        non_repeats.append(int(line.split('\t')[14]))
    all.append(int(line.split('\t')[14]))

#lim = max(abs(min(all)), abs(max(all)))
lim = 100000

BINS = 2000

fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.set_size_inches(7, 7, forward=True)
ax1.hist(all, bins=BINS)
ax1.set_title('All')
ax1.set_xlim([-lim, lim])
ax1.xaxis.set_visible(False)
#ax1.set_xticklabels([min(all), 0, max(all)], rotation=45, ha='right')
ax2.hist(non_repeats, bins=BINS)
ax2.set_title('Non-Repeats')
ax2.set_xlim([-lim, lim])
ax2.xaxis.set_visible(False)
ax3.hist(repeats, bins=BINS)
ax3.set_title('Repeats')
ax3.set_xlim([-lim, lim])
ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText=True)
ax3.set_xlabel('Distance to TSS (bp)', fontsize=16)

fig.tight_layout()

plt.savefig('/main/figures/'+ID+'tss.png', dpi=500)
