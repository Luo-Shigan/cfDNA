
#!/bin/bash
#conda环境：scRNA-seq、RNA-SNP（entrez-direct）
set -e
###PRJ2SRA function parameter introduction
## input=$1 the second column of infile muste be PRJ number (the serial number of BioProject)
## outdir=$2
function PRJ2SRA(){
    input=$1
    outdir=$2
    log=$3
    echo "PRJ2SRA begin" > ${log}
    arr_prj=$(awk '{print $2}' ${input})
    for i in ${arr_prj[@]};do
        echo "${i}开始处理"
        if esearch -db sra -query $i | efetch -format runinfo > "${outdir}/${i}-SRA.csv"; then
            echo "${i} 处理完成"
            # echo $?
        else
            echo "${i} 处理失败，跳过" >> ${log}
        fi
        if [ ! -s "${outdir}/${i}-SRA.csv" ]; then
            echo "${i} 没有信息" >> ${log}
            rm -f ${outdir}/${i}-SRA.csv
        fi
    done
    echo "完毕"
}
export -f PRJ2SRA
sra=/ChIP_seq_2/StemCells/data/download/GSE.txt
outdir=/ChIP_seq_2/StemCells/data/sra
log=/ChIP_seq_2/StemCells/log/download/PRJ2SRA-wrong.log
# PRJ2SRA ${sra} ${outdir} ${log}

###PRJ2SRASingle function parameter introduction
## PRJ=$1 the PRJ number
## outdir=$2
## log=$3 see if excute successfully
function PRJ2SRASingle(){
    PRJ=$1
    outdir=$2
    log=$3
    echo "${PRJ}开始处理"
    if esearch -db sra -query ${PRJ} | efetch -format runinfo > "${outdir}/${PRJ}-SRA.csv"; then
        echo "PRJ2SRASingle ${PRJ} 处理成功" >> ${log}
        # echo $?
    else
        echo "${PRJ} 处理失败，跳过" >> ${log}
    fi
    if [ ! -s "${outdir}/${PRJ}-SRA.csv" ]; then
        echo "${PRJ} 没有信息" >> ${log}
        rm -f ${outdir}/${PRJ}-SRA.csv
    fi
    echo "完毕"
}
export -f PRJ2SRASingle

###PRJ2SAMPLE function parameter introduction
## input=$1 the second column of infile muste be PRJ number (the serial number of BioProject)
## outdir=$2
function PRJ2SAMPLE(){
    input=$1
    outdir=$2
    log=$3
    echo "PRJ2SAMPLE begin" > ${log}
    arr_prj=$(awk '{print $2}' ${input})
    for i in ${arr_prj[@]};do
        echo "${i}开始处理"
        if esearch -db sra -query $i | efetch -format docsum > ${outdir}/${i}-SAMPLE.csv; then
            echo "${i} 处理完成"
            cat ${outdir}/${i}-SAMPLE.csv | xtract \
                -pattern DocumentSummary \
                -element Title,Bioproject,Biosample,Run@acc > tmp && mv tmp ${outdir}/${i}-SAMPLE.csv
        else
            echo "${i} 处理失败，跳过" >> ${log}
        fi

    done
    echo "完毕"
}
export -f PRJ2SAMPLE
# sra=/ChIP_seq_2/StemCells/data/download/GSE.txt
# outdir=/ChIP_seq_2/StemCells/data/sra
# log=/ChIP_seq_2/StemCells/log/download/PRJ2SRA-wrong.log
# PRJ2SAMPLE ${sra} ${outdir} ${log}
### RNAseq function parameter introduction
## input=$1 the second column of infile muste be PRJ number (the serial number of BioProject)
## outdir=$2

###PRJ2SAMPLESingle function parameter introduction
## PRJ=$1 the PRJ number
## outdir=$2
## log=$3 see if excute successfully
function PRJ2SAMPLESingle(){
    PRJ=$1
    outdir=$2
    log=$3
    echo "${PRJ}开始处理"
    if esearch -db sra -query ${PRJ} | efetch -format docsum > ${outdir}/${PRJ}-SAMPLE.csv; then
        echo "PRJ2SAMPLESingle ${PRJ} 处理成功" >> ${log}
        cat ${outdir}/${PRJ}-SAMPLE.csv | xtract \
            -pattern DocumentSummary \
            -element Title,Bioproject,Biosample,Run@acc > tmp && mv tmp ${outdir}/${PRJ}-SAMPLE.csv
    else
        echo "${PRJ} 处理失败，跳过" >> ${log}
    fi
    echo "完毕"
}
export -f PRJ2SAMPLESingle

###TargetSeq function parameter introduction

function TargetSeq(){
    input=$1
    outdir=$2
    type=$3
    arr_prj=$(awk '{print $2}' ${input})
    current_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    python_script="${current_dir}/extract.py"
    for i in ${arr_prj[@]};do
        infile=${outdir}/${i}-SRA.csv
        python ${python_script} --input $infile --sep , --colKey LibraryStrategy --colValue ${type} --output ${outdir}/${i}-SRA.ATAC-seq.csv
    done

}
export -f TargetSeq

### GetGSM function parameter introduction
## input=$1 the second column of infile muste be PRJ number (the serial number of BioProject)
## the GSM must be the NO. 30 col
function GetGSM(){
    input=$1
    outdir=$2
    type=$3
    : > ${outdir}/All-GSM.${type}.csv
    arr_prj=$(awk '{print $2}' ${input})
    for i in ${arr_prj[@]};do
        infile=${outdir}/${i}-SRA.${type}.csv
        #Run avgLength LibraryLayout ScientificName SampleName
        awk -F "\t" 'BEGIN {OFS = "\t"} NR > 1 {print $1,$7,$16,$29,$30}' ${infile} >>  ${outdir}/All-GSM.${type}.csv
    done
}
export -f GetGSM
### GetHtml function parameter introduction
## infile=$1 should inclue the GSM number that you want; for example the GSM libraryStrategy should be RNAseq which was decided by GetGSM
## it is no need to excute again this function as long as there exists 
function GetHtml(){
    infile=$1
    log=$2
    echo "GetHtml begin" > ${log}
    baseurl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
    awk '{print $1}' ${infile} | while read -r GSM;do
        echo $GSM
        url="${baseurl}${GSM}"
        if wget  -O "/ChIP_seq_2/StemCells/data/GSM/${GSM}.html" $url;then
            echo "Get ${url} successfully"
        else
            echo "Get ${url} unsuccessfully" >> ${log}
        fi
    done
}
export -f GetHtml
### GetHtml function parameter introduction
## the third col of infile must be GSM which was decided by GetGSM function
function JudgeRNAseq(){
    infile=$1
    outfile=$2
    cp ${infile} ${outfile}
    awk -F "\t" '{print $5}' ${infile} | while read -r GSM;do
        if grep -iqE "scRNA|Single-cell|10x" /ChIP_seq_2/StemCells/data/GSM/${GSM}.html;then
            sed -i "/${GSM}/ s/$/\tscRNA-seq/" ${outfile}
            # echo "${GSM} scRNA-seq" >> ${outfile}
        elif grep -iqE "bulk|VAHTS Universal V8 RNA-seq|Total RNA Extraction|magnetic beads with oligo" /ChIP_seq_2/StemCells/data/GSM/${GSM}.html;then
            sed -i "/${GSM}/ s/$/\tbulkRNA-seq/" ${outfile}
            # echo "${GSM} bulkRNA-seq" >> ${outfile}
        else
            sed -i "/${GSM}/ s/$/\totherRNA-seq/" ${outfile}
            # echo "${GSM} otherRNA-seq" >> ${outfile}
        fi
    done
    rm -f ${infile}

}
export -f JudgeRNAseq

### GSMRunInfo function parameter introduction
## GSM=$1
## outfile=$2
##根据GSM号获取SRA运行信息，当outfile不存在时，写入表头。不同outfile代表不同项目
# function GSMRunInfo(){
#     GSM=$1
#     outfile=$2
#     echo ${GSM}
#     esearch -db gds -query ${GSM} | elink -target sra | efetch -format runinfo > tmp
#     if [ ! -e ${outfile} ]; then
#         head -n 1 tmp > ${outfile}
#     fi
#     grep ${GSM} tmp >> ${outfile}
#     cat ${outfile}
#     rm -f tmp
# }
function GSMRunInfo(){
    local GSM="$1"
    local outfile="$2"
    if [[ -z "$GSM" || -z "$outfile" ]]; then
        echo "Usage: GSMRunInfo <GSM_ID> <output_file>"
        return 1
    fi

    echo "Querying GSM: ${GSM}..."
    
    esearch -db gds -query "${GSM}" | elink -target sra | efetch -format runinfo > tmp
    if [[ $? -ne 0 || ! -s tmp ]]; then
        echo "Error: Failed to fetch data for ${GSM}."
        rm -f tmp
        return 1
    fi

    if [[ ! -e "$outfile" ]]; then
        head -n 1 tmp > "$outfile"
    fi

    # 采用更严格的匹配方式，仅匹配完整的 GSM ID
    grep -w "${GSM}" tmp >> "$outfile"

    echo "Results saved to ${outfile}"

    rm -f tmp
}

export -f GSMRunInfo


# GSMRunInfo GSM5837959 /ChIP_seq_2/StemCells/RNASNP202503/data/GSMRunInfo.csv
# GSMRunInfo GSM5837962 /ChIP_seq_2/StemCells/RNASNP202503/data/GSMRunInfo.csv
# GSMRunInfo GSM5837965 /ChIP_seq_2/StemCells/RNASNP202503/data/GSMRunInfo.csv
# GSMRunInfo GSM5837968 /ChIP_seq_2/StemCells/RNASNP202503/data/GSMRunInfo.csv
# GSMRunInfo GSM5837971 /ChIP_seq_2/StemCells/RNASNP202503/data/GSMRunInfo.csv
# GSMRunInfo GSM5837974 /ChIP_seq_2/StemCells/RNASNP202503/data/GSMRunInfo.csv