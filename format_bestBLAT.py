#!/usr/bin/env python3
#this is run through the pipeline

import sys
import pandas as pd

ID = sys.argv[1]

filename = '/main/results/' + ID + '.bestBLAT.csv'

# read file
with open(filename) as f:
    reads = f.readlines()

# pick the best
keep = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
filtered = []
current_id = ""
bestis = [0]
highscore = 0
for i in range(len(reads)):
    line = reads[i].split('\t')
    if line[0] in keep:
        id = ':'.join(line[3].split(':')[:7])
        if id == current_id:
            if int(line[4]) > highscore:
                bestis = [i]
                highscore = int(line[4])
            elif int(line[4]) == highscore:
                bestis.append(i)
        else:
            if current_id != "":
                for besti in bestis:
                    best = reads[besti]
                    best = best.split('\t')
                    best[3] = current_id
                    if float(best[5]) >= 95:
                        filtered.append('\t'.join(best))
                bestis = [i]
                highscore = int(line[4])
            current_id = id

# format
lines = []
for line in filtered:
    lines.append(line.split('\t')[:4])

# write file
with open('/main/results/'+ID+'.formatted.txt', 'w') as file:
    for line in lines:
        file.write('\t'.join(line) + '\n')

# no multiple matches
dupes = []
newlines = []
for i in range(len(lines)-1):
    id = lines[i][3]
    if not id in dupes:
        if id == lines[i+1][3]:
            dupes.append(id)
        else:
            newlines.append(lines[i])

with open('/main/results/'+ID+'.formatted_nomutliplematches.txt', 'w') as file:
    for line in newlines:
        file.write('\t'.join(line) + '\n')

with open('/main/results/'+ID+'.mutliplematchesid.txt', 'w') as file:
    for id in dupes:
        file.write(id + '\n')