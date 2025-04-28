#!/bin/bash
### NCBIdownload function introduction
##ID=$1 the sra run number
##outdir=$2 output dir
##log=$3 log file
function NCBIdownload(){
    ID=$1
    outdir=$2
    log=$3
    pre_donwload(){
        prefetch ${ID} -O ${outdir}
        mv ${outdir}/${ID}/${ID}.sra ${outdir}
        fastq-dump --split-3 --gzip ${outdir}/${ID}.sra
        #--split-files

    }
    if pre_donwload $ID -O /ChIP_seq_2/StemCells/data/fq;then
        echo "${ID} 下载成功" >> ${log}
    else
        echo "${ID} 下载失败" >> ${log}
    fi
    

}
export -f NCBIdownload
# outdir=/ChIP_seq_2/StemCells/data/fq
# sralist=/ChIP_seq_2/StemCells/log/download/ascp_download2.log
# nohup awk -F " " '$3 == "下载失败" {print $1,$2}' ${sralist} | parallel -j 3 --joblog /ChIP_seq_2/StemCells/log/download/parallel.log --colsep ' ' download {1} ${outdir} ${log} > /ChIP_seq_2/StemCells/log/download/ascp.log 2>&1 &