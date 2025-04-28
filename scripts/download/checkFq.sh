#!/bin/bash
## check whether fq has been download already or not


function check(){
    filePath=$1
    awk -F "\t" '{print $6,$7}' ${filePath} | while read -r SRR libraryOut;do
        dir=/ChIP_seq_2/StemCells/data/fq
        log=/ChIP_seq_2/StemCells/log/download
        if [ -e ${dir}/${SRR}_1.fastq.gz ] && [ -e ${dir}/${SRR}_2.fastq.gz ]; then
            echo -e "${SRR}\t${libraryOut}\tfq文件存在" >> ${log}/check.log
        else
            echo -e "${SRR}\t${libraryOut}\tfq文件不存在" >> ${log}/check.log
        fi
    done
}
:> /ChIP_seq_2/StemCells/log/download/check.log
file1=/ChIP_seq_2/StemCells/data/Human_bulk.csv
file2=/ChIP_seq_2/StemCells/data/Mouse_bulk.csv
file3=/ChIP_seq_2/StemCells/data/Mouse_ATAC.csv
check ${file1}
check ${file2}
check ${file3}