#!/bin/bash
### caculate the expected TE expression in genome
### (TE length)/(genome length)*(total reads)
#args 1 TE name
#args 2 gtf file
## illustrate: one line is seen as specific TE as long as it have the specific TE name in gend_id colum
function TEExpected(){
    TE=$1
    gtf=$2
    outfile=$3
    # cut is faster than awk in this Scenario
    grep "${TE}" ${gtf} | cut -d';' -f1 | cut -f4,5,9 >> ${outfile}

}
export -f TEExpected
gtf=/ChIP_seq_2/StemCells/data/genome/GRCh38_GENCODE_rmsk_TE.gtf
#eg: TE='gene_id "L1P5"'
outfile=a.txt
TE_list=$(cut -d';' -f1 ${gtf} | cut -f9 | sort -u)
cut -d';' -f1 ${gtf} | cut -f9 | sort -u | while read TE;do
    TEExpected "${TE}" ${gtf} ${outfile}
done

# 