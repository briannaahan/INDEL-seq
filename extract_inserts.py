#!/usr/bin/env python3
#this is run through the pipeline

import sys
import pandas as pd

ID = sys.argv[1]

filename = '/main/results/' + ID + '.sam'

# read file
with open(filename) as f:
    reads = f.readlines()

# remove duplicates, find insertions, and clean up
cigar = {}
def not_duplicate(read, i):
    if not(read[5] in cigar.keys()):
        return True
    for seq in cigar[read[5]]:
        difference = sum(seq[j] != read[9][int(i[0]):int(i[0])+int(i[1])][j] for j in range(len(seq)))
        if difference/len(seq) < 0.14:
            return False
    return True

def not_telomere(read, i):
    seq = read[9][int(i[0]):int(i[0])+int(i[1])]
    if 'TTAGGGTTAGGG' in seq or 'CCCTAACCCTAA' in seq:
        return False
    count1 = seq.count('TTAGGG')
    count2 = seq.count('CCCTAA')
    return count1*6 < 0.3*len(seq) and count2*6 < 0.3*len(seq)

insertions = []
reads = reads[2:]
for read in reads:
    read = read.split('\t')
    if read[5].count('I')==1 and read[5].count('M')==2 and read[5].count('D')==0:
        i = read[5].replace("I","M").split("M")
        if not_duplicate(read, i) and not_telomere(read, i):
            try:
                cigar[read[5]].append(read[9][int(i[0]):int(i[0])+int(i[1])])
            except:
                cigar[read[5]] = [read[9][int(i[0]):int(i[0])+int(i[1])]]
            insert = [read[0], read[5], read[9][int(i[0]):int(i[0])+int(i[1])]]
            insertions.append(insert)


# export to csv
df = pd.DataFrame(insertions, columns=['SEQ ID','CIGAR','INS'])
df.to_csv('/main/results/'+ID+'.INS.csv',index=False)

# export to fasta
with open('/main/results/'+ID+'.INS.fasta', 'w') as file:
    for insert in insertions:
        file.write('>'+insert[0]+'\n'+insert[2]+'\n')