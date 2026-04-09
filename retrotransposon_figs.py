#!/usr/bin/env python3
#this is run through the pipeline

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

ID = sys.argv[1]
filename1 = '/main/results/'+ID+'.INS.fasta.out'
filename2 = '/main/results/'+ID+'.AsiSIdistance.csv' 

with open(filename1) as f:
    repeats = f.readlines()

with open(filename2) as f:
    asisi = f.readlines()

repeat_heading = repeats[:3]
repeats = repeats[3:]

temp = []
ids_distance = []

for line in repeats:
    temp.append(line.split()[4])


heading = asisi[0]
asisi = asisi[1:]
for line in asisi:
    if line.split(',')[0] in temp:
        ids_distance.append((line.split(',')[0],int(line.rstrip().split(',')[5])))
# ids_distance now contains the seq IDs in both AsisI and repeat files, along with AsisI distance. ex.: [seqID, distance]

class1 = {}
class2 = {}
class3 = {}

for line in repeats:
    line = line.split()
    for pair in ids_distance:
        if line[4] == pair[0]:
            if pair[1] < 2000:
                if line[10] in class1.keys():
                    class1[line[10]] += 1
                else:
                    class1[line[10]] = 1
            elif pair[1] < 100000:
                if line[10] in class2.keys():
                    class2[line[10]] += 1
                else:
                    class2[line[10]] = 1
            else:
                if line[10] in class3.keys():
                    class3[line[10]] += 1
                else:
                    class3[line[10]] = 1

tot1 = 0
tot2 = 0
tot3 = 0
for line in asisi:
    if int(line.rstrip().split(',')[5]) < 2000:
        tot1 +=1
    elif int(line.rstrip().split(',')[5]) < 100000:
        tot2 +=1
    else:
        tot3 += 1

retro1 = 0
retro2 = 0
retro3 = 0
for label in class1.keys():
    if "SINE" in label or "LINE" in label or "ERV" in label or "SVA" in label:
        retro1 += class1[label]
for label in class2.keys():
    if "SINE" in label or "LINE" in label or "ERV" in label or "SVA" in label:
        retro2 += class2[label]
for label in class3.keys():
    if "SINE" in label or "LINE" in label or "ERV" in label or "SVA" in label:
        retro3 += class3[label]

fig, ax = plt.subplots()
l = ax.pie(class1.values(), autopct='%1.1f%%', startangle=-90)
ax.set_title("Class I Repeats. Class I Retrotransposon Insertions: " + str(round(retro1/tot1*100,2)) + "%")

labels1=class1.keys()

for label, t in zip(labels1, l[1]):
    x, y = t.get_position()
    angle = int(math.degrees(math.atan2(y, x)))
    ha = "left"

    if x<0:
        angle -= 180
        ha = "right"

    plt.annotate(label, xy=(x,y), rotation=angle, ha=ha, va="center", rotation_mode="anchor", size=6)

plt.savefig('/main/figures/'+ID+'retro1.png', dpi=500)


fig, ax = plt.subplots()
l = ax.pie(class2.values(), autopct='%1.1f%%', startangle=-90)
ax.set_title("Class II Repeats. Class II Retrotransposon Insertions: " + str(round(retro2/tot2*100,2)) + "%")

labels2=class2.keys()

for label, t in zip(labels2, l[1]):
    x, y = t.get_position()
    angle = int(math.degrees(math.atan2(y, x)))
    ha = "left"

    if x<0:
        angle -= 180
        ha = "right"

    plt.annotate(label, xy=(x,y), rotation=angle, ha=ha, va="center", rotation_mode="anchor", size=6)

plt.savefig('/main/figures/'+ID+'retro2.png', dpi=500)


fig, ax = plt.subplots()
l = ax.pie(class3.values(), autopct='%1.1f%%', startangle=-90)
ax.set_title("Class III Repeats. Class III Retrotransposon Insertions: " + str(round(retro3/tot3*100,2)) + "%")

labels3=class3.keys()

for label, t in zip(labels3, l[1]):
    x, y = t.get_position()
    angle = int(math.degrees(math.atan2(y, x)))
    ha = "left"

    if x<0:
        angle -= 180
        ha = "right"

    plt.annotate(label, xy=(x,y), rotation=angle, ha=ha, va="center", rotation_mode="anchor", size=6)

plt.savefig('/main/figures/'+ID+'retro3.png', dpi=500)