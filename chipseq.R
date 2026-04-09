#BiocManager::install("ChIPSeeker")
#BiocManager::install("org.Hs.eg.db")

args <- commandArgs(trailingOnly = TRUE)
ID <- args[1]

library(readr)
library(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

path <- paste("/Users/briannahan/Desktop/INDEL-seq-main/results/", ID, ".formattedfiltered.txt", sep="")
peakAnno <- annotatePeak(path, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

df <- as.data.frame(peakAnno, row.names = NULL, optional = FALSE)
df1 <- subset(df, seqnames=='chr1' | seqnames=='chr10' | seqnames=='chr11' | seqnames=='chr12' | seqnames=='chr13' | seqnames=='chr14' | seqnames=='chr15' | seqnames=='chr16' | seqnames=='chr17' | seqnames=='chr18' | seqnames=='chr19' | seqnames=='chr2' | seqnames=='chr20' | seqnames=='chr21' | seqnames=='chr22' | seqnames=='chr23' | seqnames=='chr3' | seqnames=='chr4' | seqnames=='chr5' | seqnames=='chr6' | seqnames=='chr7' | seqnames=='chr8' | seqnames=='chr9' | seqnames=='chrX')

save_path <- paste("/Users/briannahan/Desktop/INDEL-seq-main/results/", ID, ".annotation.csv", sep="")
write_delim(df1, save_path, delim='\t')
