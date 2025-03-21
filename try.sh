#!/bin/bash
bam=/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted/V350318972_L01_601.bam
outbam=V350318972_L01_601_dup.bam
outmetrics=V350318972_L01_601_dedup-metrics.txt
function dedup(){
    bam=$1
    outbam=$2
    outmetrics=$3
    gatk MarkDuplicates \
    -I ${bam} \
    -O ${outbam} \
    -M ${outmetrics} \
    --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate

}
function computeGCBias_(){
    threads=$1
    gn_size=$2
    png=$3
    gn_2bit_path=$4
    bam=$5
    freq=$6
    computeGCBias \
        --numberOfProcessors ${threads} \
        --effectiveGenomeSize ${gn_size} \
        --biasPlot ${png} \
        --plotFileFormat png \
        --genome ${gn_2bit_path} \
        --bamfile ${bam} \
        --GCbiasFrequenciesFile ${freq} 

}
threads=25
gn_size=2913022398
png=GC.png
gn_2bit_path=/ChIP_seq_2/Data/serum/cfDNA/output/gn/hg38.2bit
bam=/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_601.bam
freq=out.freq
# computeGCBias_ ${threads} ${gn_size} ${png} ${gn_2bit_path} ${bam} ${freq}
function correctGCBias_(){
    threads=$1
    gn_size=$2
    gn_2bit_path=$3
    bam=$4
    freq=$5
    outbam=$6
    correctGCBias \
    --numberOfProcessors ${threads} \
    --effectiveGenomeSize ${gn_size} \
    --genome ${gn_2bit_path} \
    --bamfile ${bam} \
    --GCbiasFrequenciesFile ${freq} \
    --correctedFile ${bam}
}
threads=25
gn_size=2913022398
gn_2bit_path=/ChIP_seq_2/Data/serum/cfDNA/output/gn/hg38.2bit
bam=/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_601.bam
freq=/ChIP_seq_2/Data/serum/cfDNA/output/GC_bias/V350318972_L01_604.txt
outbam=out.bam
# correctGCBias_  ${threads} ${gn_size} ${gn_2bit_path} ${bam} ${freq} ${outbam}
function wgs(){
    bam=$1
    out_sr=$2
    out_er=$3
    Rscript scripts/wgs_qc.R -i ${bam} -sr ${out_sr} -er ${out_er}
}
bam=/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_601.bam
out_sr=./temp/a.txt
out_er=./temp/b.txt
# wgs ${bam} ${out_sr} ${out_er}
function count_matrix_gene(){
    threads=$1
    region=$2
    gtf=$3
    outGene_matrix=$4
    shift 4
    bam=$@
    featureCounts -T ${threads} -O -t ${region} -g gene_id -M -p \
    -a ${gtf} \
    -o ${outGene_matrix} ${bam}

    /usr/bin/Rscript scripts/multiFeatureCounts2countMat.R \
            -i ${outGene_matrix}  \
            -o b.txt 
    /usr/bin/Rscript scripts/run-NormCountMat.R \
    -i b.txt \
    -o c.txt \
    -m TMM 

}
threads=25
region=gene
gtf=/ChIP_seq_2/StemCells/data/genome/gencode.v47.primary_assembly.annotation.gtf
outGene_matrix=a.txt
bam="/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_601.bam /ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_602.bam"
# count_matrix_gene ${threads} ${region} ${gtf} ${outGene_matrix} ${bam}
function FragmentHist(){
    out=$1
    shift 1
    files=$@
    python /ChIP_seq_2/Data/serum/cfDNA/scripts/getFragmentSize.py --input ${files} --out ${out}
}
files=$(find /ChIP_seq_2/Data/serum/cfDNA/output/stats -type f | tr '\n' ' ')
# echo ${files}
outfile=c.txt
# FragmentHist ${outfile} ${files}
function prepareMantaConfig(){
    bam=$1
    gn_fa_path=$2
    dir=$3
    configManta.py --bam ${bam} \
    --referenceFasta ${gn_fa_path} \
    --runDir ${dir}
}
bam="/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_601.bam"
gn_fa_path="/ChIP_seq_2/StemCells/data/genome/GRCh38.primary_assembly.genome.fa"
dir=temp
# prepareMantaConfig ${bam} ${gn_fa_path} ${dir}
function runMantan(){
    script=$1
    python ${script} -j 6
}
script=/ChIP_seq_2/Data/serum/cfDNA/temp/runWorkflow.py
# runMantan ${script}
function calFragSizeRatio(){
    shortMat=$1
    longMat=$2
    output=$3
    /usr/bin/Rscript scripts/frag/get-FragmentSizeMat.R \
    --short ${shortMat} \
    --long ${longMat} \
    --outfile ${output}
}
shortMat=/ChIP_seq_2/Data/serum/cfDNA/output/shortLong/matrix/CPM-TMM_matrix_gene_short.correctGC.txt
longMat=/ChIP_seq_2/Data/serum/cfDNA/output/shortLong/matrix/CPM-TMM_matrix_gene_long.correctGC.txt
output=a.txt
# calFragSizeRatio ${shortMat} ${longMat} ${output}
function count_matrix_TE(){
    region=$1
    gtf=$2
    gene_matrix=$3
    tmp=$4
    gene_CPM_TMM=$5
    shift 5
    bam=$@
    # featureCounts -T 16 -B -O -t ${region} -g gene_id -M -p \
    # -a ${gtf} \
    # -o ${gene_matrix} ${bam} 

    /usr/bin/Rscript scripts/multiFeatureCounts2countMat.R \
        -i ${gene_matrix}  \
        -o ${tmp} ;
    mv ${tmp} ${gene_matrix};

    /usr/bin/Rscript scripts/run-NormCountMat.R \
        -i ${gene_matrix} \
        -o ${gene_CPM_TMM} \
        -m TMM 
}
region="exon"
gtf="/ChIP_seq_2/StemCells/data/genome/GRCh38_GENCODE_rmsk_TE.gtf.gz"
gene_matrix=a.txt
tmp=b.txt
gene_CPM_TMM=TMM.txt
bam="/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_601.bam /ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_602.bam"
# count_matrix_TE $region $gtf ${gene_matrix} ${tmp} ${gene_CPM_TMM} ${bam}
function TEtranscripts(){
    geneGtf=$1
    teGtf=$2
    project=$3
    bam=$4
    TEcount -b ${bam} \
    --GTF ${geneGtf} \
    --TE ${teGtf} \
    --sortByPos --mode multi \
    --project ${project}
}
geneGtf=/ChIP_seq_2/StemCells/data/genome/gencode.v47.primary_assembly.annotation.gtf
teGtf=/ChIP_seq_2/StemCells/data/genome/GRCh38_GENCODE_rmsk_TE.gtf
project=bam
bam=/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972_L01_601.bam
# TEtranscripts ${geneGtf} ${teGtf} ${project} ${bam}
# TEtranscripts只适用于RNAseq数据,对于DNAseq的数据是否适用,仍然有很大的疑问

bam=/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/V350318972.bam
outfile=/ChIP_seq_2/Data/serum/cfDNA/output/bam-sorted-deduped/stats/depth_coverge.csv

# echo -e "sample""\t""genome_depth""\t""coverage_depth""\t""coverage""\t""coverage_ratio" >> ${outfile}
# time samtools depth ${bam} | \
#     awk '{sum+=$3} END { print "%"  "\t" sum/3000000000 "\t" sum/NR "\t" NR "\t" NR/3000000000 } ' \
#     >> ${outfile}