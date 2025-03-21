#!/bin/bash
bams=$@
function depth(){
    bam=$1
    sample_id=$(basename -s .bam ${bam})
    outfile=$2
        # 只在输出文件不存在时写入标题行
    if [ ! -f ${outfile} ]; then
        echo -e "sample\tgenome_depth\tcoverage_depth\tcoverage\tcoverage_ratio" > ${outfile}
    fi
    # the bam must be sorted
    samtools depth ${bam} | \
        awk -v sample=${sample_id} '{sum+=$3} END { print sample  "\t" sum/3000000000 "\t" sum/NR "\t" NR "\t" NR/3000000000 } ' \
        >> ${outfile}
}
export -f depth
depth /ChIP_seq_2/Data/serum/cfDNA/output/RG/V350318972.bam /ChIP_seq_2/Data/serum/cfDNA/output/RG/stats/depth_coverge.csv