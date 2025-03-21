#!/bin/bash
### you should excute the below command before excuting cfDNA.smk
fa=/ChIP_seq_2/StemCells/data/genome/GRCh38.primary_assembly.genome.fa
indexdir=/ChIP_seq_2/Data/serum/cfDNA/data/index/GRCh38
mkdir /ChIP_seq_2/Data/serum/cfDNA/data/index/
bwa-mem2 index -p ${indexdir} ${fa}
####### make BSgenome.Hsapiens.Gencode.GRCh38
Rscript /ChIP_seq_2/Data/serum/cfDNA/scripts/BSgenome.R

########### make reference gene.bed#############
awk '$3 == "gene"' genes.gtf > genes_only.gtf