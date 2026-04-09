#!/usr/bin/env python3
#this is run through the pipeline

import sys
import pandas as pd
import matplotlib.pyplot as plt
import math

ID = sys.argv[1]
filename1 = '/main/results/'+ID+'.INS.fasta.out' # '/Users/briannahan/Desktop/INDEL-seq-main/results/' + ID + '.INS.fasta.out'  
filename2 = '/main/results/'+ID+'.AsiSIdistance.csv' # '/Users/briannahan/Desktop/INDEL-seq-main/results/' + ID + '.AsiSIdistance.csv'   

with open(filename1) as f:
    repeats = f.readlines()

with open(filename2) as f:
    asisi = f.readlines()

repeat_heading = repeats[:3]
repeats = repeats[3:]

temp = []
ids = []

for line in repeats:
    temp.append(line.split()[4])


heading = asisi[0]
asisi = asisi[1:]
for line in asisi:
    if line.split(',')[0] in temp:
        ids.append(line.split(',')[0])

lines = []

for line in asisi:
    line = line.rstrip().split(',')
    if line[0] in ids:
        line.append('TRUE')
        line.append('blue')
    else:
        line.append('FALSE')
        line.append('red')
    lines.append(line)

heading = heading.rstrip().split(',')
heading.append('isRepeat')
heading.append('color')

df = pd.DataFrame(lines, columns=heading)
#df.to_csv('/main/results/'+ID+'.violin.csv',index=False) # '/Users/briannahan/Desktop/INDEL-seq-main/results/'+ID+'.violin.csv'
newlines = []
pie = {}
for line in repeats:
    line = line.split()
    if line[4] in ids:
        newlines.append(line)
        if line[10] in pie.keys():
            pie[line[10]] += 1
        else:
            pie[line[10]] = 1

retrotransposon = 0
total = 0
for label in pie.keys():
    total += pie[label]
    if "SINE" in label or "LINE" in label or "ERV" in label or "SVA" in label:
        retrotransposon += pie[label]

df = pd.DataFrame(newlines)
#df.to_csv('/main/results/'+ID+'.repeats.csv',index=False) #####

fig, ax = plt.subplots()
l = ax.pie(pie.values(), autopct='%1.1f%%', startangle=-90)
ax.set_title("Retrotransposons: " + str(round(retrotransposon/total*100,2)) + "%")

labels=pie.keys()

for label, t in zip(labels, l[1]):
    x, y = t.get_position()
    angle = int(math.degrees(math.atan2(y, x)))
    ha = "left"

    if x<0:
        angle -= 180
        ha = "right"

    plt.annotate(label, xy=(x,y), rotation=angle, ha=ha, va="center", rotation_mode="anchor", size=6)

plt.savefig('/main/figures/'+ID+'repeat_pie.png', dpi=500)