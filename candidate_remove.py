#!/usr/bin/env python3

import sys
import pandas as pd

ID = sys.argv[1]

filename = '/main/results/' + ID + '.AsiSIdistance.csv'

with open(filename) as f:
    lines = f.readlines()
heading = lines[0].rstrip().split(',')
lines = lines[1:]

newlines = []

thisID = ""
thisHighscore = float('inf')
bestLine = ""
for line in lines:
    line = line.rstrip().split(',')
    id = line[0]
    distance = int(line[5])
    if id == thisID:
        if distance < thisHighscore:
            thisHighscore = distance
            bestLine = line
    else:
        if thisID != "":
            newlines.append(bestLine)
        thisID = id
        thisHighscore = distance
        bestLine = line

# overwrite old csv
df = pd.DataFrame(newlines, columns=heading)
# df.to_csv('/main/results/'+ID+'.AsiSIdistance.csv',index=False)

# format and filter for annotations
new = df[["INS Chrom","INS Start","INS End","Seq ID"]]
lines = new.values.tolist()
with open('/main/results/'+ID+'.formattedfiltered.txt', 'w') as file:
    for line in lines:
        file.write('\t'.join(line) + '\n')