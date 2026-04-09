#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
from multiprocessing import Pool

fq1 = sys.argv[1]
fq2 = sys.argv[2]
barcodes = sys.argv[3] #cagatc

##preprocess
os.makedirs('/main/preprocess/',exist_ok=True)
#cut the first 6 bp NNNNNN
subprocess.call('cutadapt -u 6 -U 6 -j 0 -o /main/preprocess/' + os.path.basename(fq1).split('.')[0] + '.trimmed.fastq.gz -p /main/preprocess/' + os.path.basename(fq2).split('.')[0] + '.trimmed.fastq.gz ' + fq1 + ' ' + fq2,shell=True)

#demultx
subprocess.call('/opt/fastq-multx/fastq-multx -b -B ' + barcodes + ' /main/preprocess/' + os.path.basename(fq1).split('.')[0] + '.trimmed.fastq.gz preprocess/' + os.path.basename(fq2).split('.')[0] + '.trimmed.fastq.gz -o /main/preprocess/%_R1.fq.gz preprocess/%_R2.fq.gz',shell=True)

#################################

#merge forward and reverse paired-end reads
def merge(ID):
    subprocess.call('/opt/pear-0.9.11-linux-x86_64/bin/pear -f /main/preprocess/' + ID + '_R1.fq -r /main/preprocess/' + ID + '_R2.fq -o /main/preprocess/' + ID,shell=True)

N = os.cpu_count()
IDs = list(pd.read_csv(barcodes,sep="\t",header=None)[0].astype(str).values)
with Pool(N) as p:
    p.map(merge, IDs)

#alignment
os.makedirs('/main/results/',exist_ok=True)

def run_needle(ID):
    subprocess.call('needle /main/TRIM37.fasta /main/preprocess/' + ID + '.assembled.fastq -outfile /main/results/' + ID + '.sam -aformat sam -gapopen 10 -gapextend 0.5',shell=True)
print('starting needle')
with Pool(N) as p:
    p.map(run_needle, IDs)
print('needle completed')

# # keep this commented out
# # def run_gzip(ID):
# #     subprocess.call('gzip /main/results/' + ID + '.sam',shell=True)
# # with Pool(N) as p:
# #     p.map(run_gzip, IDs)




#get inserts
print('extracting inserts')
def extract_inserts(ID):
    subprocess.call('python3 /main/extract_inserts.py ' + ID, shell=True)
with Pool(N) as p:
    p.map(extract_inserts, IDs)
print('extracting inserts completed')

#################################

#BLAT
def blat(ID):
    subprocess.call('./blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 /main/hg19.2bit /main/results/' + ID + '.INS.fasta /main/results/' + ID + '.BLAToutput.psl && /main/BLAT/pslScore /main/results/'+ ID + '.BLAToutput.psl > /main/results/' + ID +'.bestBLAT.csv', shell=True)
    subprocess.call('python3 /main/format_bestBLAT.py ' + ID, shell=True)

print('starting BLAT')
with Pool(N) as p:
    p.map(blat, IDs)
print('BLAT completed')

#################################

#annotations
def annotate(ID):
    subprocess.call('Rscript /main/chipseq.R ' + ID, shell=True)
with Pool(N) as p:
    p.map(annotate, IDs)

#AsiSI distance
def distance(ID):
    subprocess.call('python3 /main/asisi_distance.py '+ID, shell=True)
with Pool(N) as p:
    p.map(distance, IDs)

#repeat annotation
def repeat_annotation(ID):
    subprocess.call('./RepeatMasker -e rmblast -species human /main/results/'+ ID +'.INS.fasta &&  ./ProcessRepeats /main/results/'+ID+'.INS.fasta.cat', shell=True)
with Pool(N) as p:
    p.map(repeat_annotation, IDs)

#figures
def figures(ID):
    subprocess.call('python3 /main/repeat_annotation.py '+ID+' && python3 /main/gene_anno_figs.py '+ID, shell=True)
with Pool(N) as p:
    p.map(figures, IDs)