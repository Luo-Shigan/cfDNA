#!/bin/bash
#### conda: RNA-SNP
### download function parameter introduction
## id=$1 should be run number,such as SRR16119550
## type=$2 should be LibraryLayout,such as PAIRED
## log=$3 should log path
## dest=$4 should be output dir
function ENAdownload(){
    #parllel子进程可能访问不到外部变量
    key=~/miniconda3/envs/RNA-SNP/etc/asperaweb_id_dsa.openssh
    id=$1
    type=$2
    log=$3
    dest=$4
    echo "${id} ${type}开始处理"
    charNumber=$(echo ${id} | wc -m)
    if [[ ${charNumber} -eq 12 ]];then
        x6=${id:0:6}
        x2=${id: -2}
        x2="0${x2}"
        echo ${x6}
        echo ${x2}
    elif [[ ${charNumber} -eq 11 ]];then
        x6=${id:0:6}
        x2=${id: -1}
        x2="00${x2}"
        echo ${x6}
        echo ${x2}
    else
        echo "sra Number既不是11位(老),也不是12位(新)"
        exit 1  # 终止整个脚本执行，并返回状态码 1

    fi
    file_path=(
    "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${x6}/${x2}/${id}/${id}_1.fastq.gz"
    "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${x6}/${x2}/${id}/${id}_2.fastq.gz"
    "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${x6}/${x2}/${id}/${id}.fastq.gz"
    )
    echo ${file_path[@]}
    function ascp_download(){
        ascp -k 1 -T -l 200m -P 33001 --file-checksum=md5 --overwrite=always -i ${key} $1 ${dest}
    }
    if [[ ${type} == "PAIRED" ]];then
        echo "双端测序"
        if ascp_download ${file_path[0]} && ascp_download ${file_path[1]};then
            echo "${id} ${type} 下载成功" >> ${log}
        elif ascp_download ${file_path[2]};then
            echo "${id} ${type}虽然为双端文件，但是可能进行了合并。更换下载地址后下载成功" >> ${log}
        else
            echo "${id} ${type} 下载失败" >> ${log}
        fi
    elif [[ ${type} == "SINGLE" ]];then
        echo "单端测序"
        if ascp_download ${file_path[2]};then
            echo "${id} ${type} 下载成功" >> ${log}
        elif ascp_download ${file_path[0]};then
            echo "${id} ${type} 虽然为单端测序文件，但可能命名习惯不一样。更换下载地址后下载成功" >> ${log}
            
        elif ascp_download ${file_path[1]};then
            echo "${id} ${type} 虽然为单端测序文件，但可能命名习惯不一样。更换下载地址后下载成功" >> ${log}
        else
            echo "${id} ${type} 下载失败" >> ${log}
        fi
    else
        echo "请输入正确的建库类型, PAIRED or SINGLE"
    fi
}
# 导出函数
export -f ENAdownload

