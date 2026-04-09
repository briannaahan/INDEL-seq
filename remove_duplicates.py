#!/usr/bin/env python3
#this is run through the pipeline

import sys
import pandas as pd

ID = sys.argv[1]
filename = '/main/results/' + ID + '.sam'

# read file
with open(filename) as f:
    reads = f.readlines()


def insert(line):
    line[5].split

cigar = []
new = []
for line in reads:
    line = line.split('\t')
    if line[5] in cigar:
        asdfad
    else:
        cigar.append((line[5],insert(line)))
        new.append('\t'.join(line))
        