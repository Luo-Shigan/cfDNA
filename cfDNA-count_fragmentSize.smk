configfile: "config/cfDNA-count.yaml"
shell.prefix("set -x; set -e;")
###### control file
def request():
    output = dict()
    output['count_matrix_gene'] = [
        expand(outdir+"/matrix/gene/count_matrix_{region}.txt",region=gene_regions),
        expand(outdir+"/matrix/gene/count_matrix_{region}.txt.summary",region=gene_regions),
        expand(outdir+"/matrix/gene/CPM-TMM_matrix_{region}.txt",region=gene_regions)
    ]
    output['count_matrix_TE'] = [
        expand(outdir+"/matrix/TE/{level}/count_matrix_{region}.txt",level=levels,region=TE_regions),
        expand(outdir+"/matrix/TE/{level}/count_matrix_{region}.txt.summary",level=levels,region=TE_regions),
        expand(outdir+"/matrix/TE/{level}/CPM-TMM_matrix_{region}.txt",level=levels,region=TE_regions)
    ]
    output['bigwig'] = [ expand(outdir+"/wig/{sample_id}.bigwig",sample_id=samples)]
    output['fragment'] = [outdir+"/fragment/FragmentSize_ori.txt",outdir+"/fragment/FragmentSize_ori.png"]
    return list(output.values())

###### control variable

    # indir = "../../data/fq" receive from command line
    # outdir = "../../output" receive from command line
    index_dir = "/ChIP_seq_2/Data/index" # for convenience 
    samples = glob_wildcards("data/fq/{sample_id}_1.fq.gz").sample_id
    ########### Gene variable
    gene_regions = config['count_matrix_TE']['region']
    ###########TE variable
    TE_regions = config['count_matrix_TE']['region']
    levels = config['count_matrix_TE']['level']


rule all:
    input:request()

##############################################upstream###########################################################
rule trimming:
    input:
        fastq1=indir+"/{sample_id}_1.fq.gz",
        fastq2=indir+"/{sample_id}_2.fq.gz"
    output:
        fastq1=outdir+"/cutadapt/{sample_id}_1.fq.gz",
        fastq2=outdir+"/cutadapt/{sample_id}_2.fq.gz",
        report1=outdir+"/log/{sample_id}/trimming_statistics_1.txt",
        report2=outdir+"/log/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir=outdir+"/cutadapt",
        quality=30,
        trim_galore="/opt/TrimGalore-0.6.10/trim_galore"
    threads: 8
    shell:
        """
        {params.trim_galore} --phred33 --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fq.gz_trimming_report.txt {output.report2}
        """

# Quality control after adapter trimming
rule qc1:
    input:
        fastq1=outdir+"/cutadapt/{sample_id}_1.fq.gz",
        fastq2=outdir+"/cutadapt/{sample_id}_2.fq.gz"
    output:
        report1=outdir+"/cutadapt/qc1/{sample_id}/{sample_id}_1_fastqc.html",
        report2=outdir+"/cutadapt/qc1/{sample_id}/{sample_id}_2_fastqc.html"
    conda:
        config['conda']['cfDNA_base'] 
    log:
        outdir + "/log/{sample_id}/qc1.log"
    params:
        outdir=outdir+"/cutadapt/qc1/{sample_id}"
    shell:
        """
        fastqc -o {params.outdir} {input.fastq1} > {log} 2>&1
        fastqc -o {params.outdir} {input.fastq2} >> {log} 2>&1
        """
# make index for bwa
rule bwa_index:
    input:
        genome=config["genome"]
    output:
        infile = index_dir + "/bwa-mem2/GRCh38/GRCh38.0123"
    log:
        index_dir + "/bwa-mem2/GRCh38/bwa_index.log"
    params:
        index =  index_dir + "/bwa-mem2/GRCh38/GRCh38"
    conda:
        config['conda']['cfDNA_base']
    shell:
        """
        # bwa-mem2 index默认会使用所有可用的核
        bwa-mem2 index -p {params.index} {input.genome} > {log} 2>&1
        """

#Mapping with bwa
rule bwa_alignment:
    input:
        fastq1=outdir+"/cutadapt/{sample_id}_1.fq.gz",
        fastq2=outdir+"/cutadapt/{sample_id}_2.fq.gz",
        infile = index_dir + "/bwa-mem2/GRCh38/GRCh38.0123"
    output:
        bam=outdir+"/bam/{sample_id}.bam",
        sam=temp(outdir+"/bam/{sample_id}.sam")
    conda:
        config['conda']['cfDNA_base'] 
    params:
        index = index_dir + "/bwa-mem2/GRCh38/GRCh38"
    threads: 15
    log:
        outdir+"/log/{sample_id}/bwa-alignment.txt"
    shell:
        """
        bwa-mem2 mem -T 0 -t {threads} {params.index} {input.fastq1} {input.fastq2} -o {output.sam} > {log} 2>&1 
        samtools view -b {output.sam} > {output.bam}
        """
# Get unmapped reads from bwa"s bam file
rule getUnaligned:
    input:
        bam=outdir+"/bam/{sample_id}.bam"
    output:
        unmapped_1=outdir+"/bam/unmapped/{sample_id}_1.fq.gz",
        unmapped_2=outdir+"/bam/unmapped/{sample_id}_2.fq.gz"
    threads: 15
    conda:
        config['conda']['cfDNA_base']    
    log:
        outdir+"/log/{sample_id}/get-unaligned.txt"
    shell:
        """
        samtools fastq -@ {threads} -1 {output.unmapped_1} -2 {output.unmapped_2} -0 /dev/null -s /dev/null -f 13 {input.bam} > {log} 2>&1
        """
# Generate flag statistics in bam file
rule flagstat:
    input:
        bam=outdir+"/bam/{sample_id}.bam"
    output:
        flagstat=outdir+"/bam/flagstat/{sample_id}.txt"
    conda:
        config['conda']['cfDNA_base']
    log:
        outdir + "/log/{sample_id}/flagstat.log"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """ 
# sort bam files
rule sort:
    input:
        bam=outdir+"/bam/{sample_id}.bam"
    output:
        bam=temp(outdir+"/bam-sorted/{sample_id}.bam"),
        bai=temp(outdir+"/bam-sorted/{sample_id}.bam.bai")
    params:
        keep_proper_pair="-f 2" if config["sort"]["onlykeep_properpair"] else "",
    threads: 4
    conda:
        config['conda']['cfDNA_base']
    shell:
        """
        samtools view -h {params.keep_proper_pair} -F 0x4  {input.bam} | samtools sort  -@ {threads}  > {output.bam}
        samtools index -@ {threads} {output.bam}
        """
# get bam files with RG (gatk requires this Read Group)
rule addReadsGroup:
    input:
        bam=outdir+"/bam-sorted/{sample_id}.bam"
    output:
        bam=outdir+"/RG/{sample_id}.bam",
        bai=outdir+"/RG/{sample_id}.bam.bai"
    params:
        id="{sample_id}",
        java="--java-options -Xmx15G",
        tmp_dir=config["tmp_dir"]
    threads: 10
    conda:
        config['conda']['cfDNA_base']        
    log:
        log=outdir+"/log/{sample_id}/addReadsGroup.log"
    shell:
        """
        gatk AddOrReplaceReadGroups {params.java} \
            --TMP_DIR {params.tmp_dir} --INPUT {input.bam} --OUTPUT {output.bam} \
            -SO coordinate --RGLB cfDNA --RGPL BGI --RGPU DNBSEQ --RGSM {params.id} > {log.log} 2>&1
        samtools index -@ {threads} {output.bam} 2 >> {log.log}
        """
# remove duplicate reads
rule dedup:
    input:
        bam=outdir+"/RG/{sample_id}.bam"
    output:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam",
        bai=outdir+"/bam-sorted-deduped/{sample_id}.bam.bai",
        metrics=outdir+"/log/{sample_id}/dedup-metrics.txt"
    log:
        outdir+"/log/{sample_id}/MarkDuplicates.log"
    threads: 16
    conda:
        config['conda']['cfDNA_base']
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} -O {output.bam} -M {output.metrics} \
            --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate > {log} 2>&1
        samtools index -@ {threads} {output.bam}
        """
##############################################downstream###########################################################
# get count matrix
rule count_matrix_gene:
    input:
        bam=expand(outdir+"/bam-sorted-deduped/{sample_id}.bam",sample_id=samples)
    output:
        gene_matrix=outdir+"/matrix/gene/count_matrix_{region}.txt",
        gene_sum=outdir+"/matrix/gene/count_matrix_{region}.txt.summary",
        gene_CPM_TMM=outdir+"/matrix/gene/CPM-TMM_matrix_{region}.txt"
    log:
        log1=outdir+"/matrix/gene/log/count_matrix_{region}.log",
        CPM_TMM=outdir+"/matrix/gene/log/CPM-TMM_matrix_{region}.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 25
    params:
        tmp=outdir+"/matrix/gene/tmp_{region}",
        region="{region}",
        gtf=config['count_matrix_gene']['gtf'],
    shell:
        """
        #gene_name
        featureCounts -T {threads} -B -O -t {params.region} -g gene_id -M -p \
            -a {params.gtf} \
            -o {output.gene_matrix} {input.bam} \
            > {log.log1} 2>&1
        
        /usr/bin/Rscript scripts/quantification/multiFeatureCounts2countMat.R \
            -i {output.gene_matrix}  \
            -o {params.tmp} ;

        /usr/bin/Rscript scripts/quantification/run-NormCountMat.R \
            -i {params.tmp} \
            -o {output.gene_CPM_TMM} \
            -m TMM \
            > {log.CPM_TMM} 2>&1
        """ 
rule count_matrix_TE:
    input:
        bam=expand(outdir+"/bam-sorted-deduped/{sample_id}.bam",sample_id=samples)
    output:
        gene_matrix=outdir+"/matrix/TE/{level}/count_matrix_{region}.txt",
        gene_sum=outdir+"/matrix/TE/{level}/count_matrix_{region}.txt.summary",
        gene_CPM_TMM=outdir+"/matrix/TE/{level}/CPM-TMM_matrix_{region}.txt"
    log:
        log1=outdir+"/matrix/TE/{level}/log/count_matrix_{region}.log",
        CPM_TMM=outdir+"/matrix/TE/{level}/log/CPM-TMM_matrix_{region}.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 25
    params:
        tmp=outdir+"/matrix/TE/{level}/tmp_{region}",
        region="{region}",
        level="{level}",
        gtf=config['count_matrix_TE']['gtf'],
    shell:
        """
        featureCounts -T {threads} -B -O -t {params.region} -g {params.level} -M -p \
            -a {params.gtf} \
            -o {output.gene_matrix} {input.bam} \
            > {log.log1} 2>&1
        
        /usr/bin/Rscript scripts/quantification/multiFeatureCounts2countMat.R \
            -i {output.gene_matrix}  \
            -o {params.tmp} ;

        /usr/bin/Rscript scripts/quantification/run-NormCountMat.R \
            -i {params.tmp} \
            -o {output.gene_CPM_TMM} \
            -m TMM \
            > {log.CPM_TMM} 2>&1
        """ 

# get WIG files (for IGV visualization)
rule wig:
    input:
        bam = outdir+"/bam/{sample_id}.bam" if not config["remove_duplications"] else outdir+"/bam-sorted-deduped/{sample_id}.bam"
    output:
        bigwig = outdir+"/wig/{sample_id}.bigwig"
    log:
        log = outdir+"/log/{sample_id}/wig.log"
    conda:
        config['conda']['GCBias']
    threads: 15 
    params:
        bin=config['wig']['bin']
    shell:
        """
        bamCoverage --binSize {params.bin} --numberOfProcessors {threads} --extendReads --normalizeUsing CPM -b {input.bam} -o {output.bigwig} > {log.log} 2>&1 
        """

############ get (original) bam stat
rule getBamStatistics:
    input:
        bam = outdir+"/bam/{sample_id}.bam" if not config["getBamStatistics_duplications"] else outdir + "/bam-sorted-deduped/{sample_id}.bam"
    output:
    #please make sure the rule all is consistent with the config["getBamStatistics_duplications"]
        summary =  outdir + "/bam/stats/{sample_id}.txt" if not config["getBamStatistics_duplications"] else outdir + "/bam-sorted-deduped/stats/{sample_id}.txt"
    conda:
        config['conda']['cfDNA_base']
    log:
        outdir + "/log/{sample_id}/getBamStatistics.log"
    shell:
        """
        samtools stats {input.bam} > {output.summary} 2>{log}
        """
# Summarize histogram of fragment length from samtools stats output
rule getFragmentSize:
    input:
        stats = expand(outdir + "/bam/stats/{sample_id}.txt",sample_id=samples) if not config["getBamStatistics_duplications"] else expand(outdir + "/bam-sorted-deduped/stats/{sample_id}.txt",sample_id=samples)
    output:
        hist = outdir+"/fragment/FragmentSize_ori.txt" if not config["getBamStatistics_duplications"] else outdir+"/fragment/FragmentSize.txt",
        png = outdir+"/fragment/FragmentSize_ori.png" if not config["getBamStatistics_duplications"] else outdir+"/fragment/FragmentSize.png"
    conda:
        config['conda']['GCBias']
    log:
        outdir + "/log/getFragmentSize.log"
    shell:
        """
            python scripts/getFragmentSize.py --input {input.stats} --out {output.hist} > {log} 2>&1
            python scripts/plotFragmentSize.py --input {output.hist} --out {output.png} >> {log} 2>&1
        """
        

