configfile: "config/cfDNA-count.yaml"
shell.prefix("set -x; set -e;")
indir = config.get('indir', '../../data/fq')
outdir = config.get('outdir', '../../output')
onlykeep_properpair = config.get('onlykeep_properpair', True)
getBamStatistics_duplications = config.get('getBamStatistics_duplications', False)
remove_duplications = config.get('remove_duplications', False)



###### control variable
samples = glob_wildcards(indir + "/{sample_id}_1.fq.gz").sample_id
## Gene variable
gene_regions = config['count_matrix_TE']['region']
##TE variable
TE_regions = config['count_matrix_TE']['region']
levels = config['count_matrix_TE']['level']
genomes = ['mouse']

###### control file
def request():
    output = dict()
    output['count_matrix_gene'] = [
        expand(outdir + "/{genome}/matrix/gene/count_matrix_{region}.txt",region=gene_regions,genome=genomes),
        expand(outdir + "/{genome}/matrix/gene/count_matrix_{region}.txt.summary",region=gene_regions,genome=genomes),
        expand(outdir + "/{genome}/matrix/gene/CPM-TMM_matrix_{region}.txt",region=gene_regions,genome=genomes)
    ]
    output['count_matrix_TE'] = [
        expand(outdir + "/{genome}/matrix/TE/{level}/count_matrix_{region}.txt",level=levels,region=TE_regions,genome=genomes),
        expand(outdir + "/{genome}/matrix/TE/{level}/count_matrix_{region}.txt.summary",level=levels,region=TE_regions,genome=genomes),
        expand(outdir + "/{genome}/matrix/TE/{level}/CPM-TMM_matrix_{region}.txt",level=levels,region=TE_regions,genome=genomes)
    ]
    output['bigwig'] = [ expand(outdir + "/{genome}/wig/{sample_id}.bigwig",sample_id=samples,genome=genomes)]
    output['fragment'] = [ expand(outdir + "/{genome}/fragment/FragmentSize.txt",genome=genomes), 
                            expand(outdir + "/{genome}/fragment/FragmentSize.png",genome=genomes)
                    ]
    return list(output.values())




rule all:
    input:request()

##############################################upstream###########################################################
rule trimming:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        fastq1 = outdir + "/{genome}/cutadapt/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{genome}/cutadapt/{sample_id}_2.fq.gz",
        report1 = outdir + "/{genome}/cutadapt/{sample_id}/trimming_statistics_1.txt",
        report2 = outdir + "/{genome}/cutadapt/{sample_id}/trimming_statistics_2.txt"
    log:
        outdir + "/log/{genome}/{sample_id}/trimming.log"
    params:
        outdir = outdir + "/{genome}/cutadapt",
        quality = 30,
        trim_galore = "/opt/TrimGalore-0.6.10/trim_galore"
    threads: 8
    shell:
        """
        {params.trim_galore} --phred33 --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} 2> {log}
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fq.gz_trimming_report.txt {output.report2}
        """

# Quality control after adapter trimming
rule qc1:
    input:
        fastq1 = outdir + "/{genome}/cutadapt/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{genome}/cutadapt/{sample_id}_2.fq.gz"
    output:
        report1 = outdir + "/{genome}/cutadapt/qc1/{sample_id}/{sample_id}_1_fastqc.html",
        report2 = outdir + "/{genome}/cutadapt/qc1/{sample_id}/{sample_id}_2_fastqc.html"
    conda:
        config['conda']['cfDNA_base'] 
    log:
        outdir + "/log/{genome}/{sample_id}/qc1.log"
    params:
        outdir = outdir + "/{genome}/cutadapt/qc1/{sample_id}"
    shell:
        """
        fastqc -o {params.outdir} {input.fastq1} > {log} 2>&1
        fastqc -o {params.outdir} {input.fastq2} >> {log} 2>&1
        """

#Mapping with bwa
rule bwaMem2_alignment:
    input:
        fastq1 = outdir + "/{genome}/cutadapt/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{genome}/cutadapt/{sample_id}_2.fq.gz",
    output:
        bam = outdir + "/{genome}/bam/{sample_id}.bam",
        sam = temp(outdir + "/{genome}/bam/{sample_id}.sam")
    conda:
        config['conda']['cfDNA_base'] 
    params:
        index = lambda wildcards: config['bwaMem2'][wildcards.genome],
    threads: 15
    log:
        outdir+"/log/{genome}/{sample_id}/bwa-alignment.txt"
    shell:
        """
        bwa-mem2 mem -T 0 -t {threads} {params.index} {input.fastq1} {input.fastq2} -o {output.sam} > {log} 2>&1 
        samtools view -b {output.sam} > {output.bam}
        """
# Get unmapped reads from bwa"s bam file
rule getUnaligned:
    input:
        bam = outdir + "/{genome}/bam/{sample_id}.bam"
    output:
        unmapped_1 = outdir + "/{genome}/bam/unmapped/{sample_id}_1.fq.gz",
        unmapped_2 = outdir + "/{genome}/bam/unmapped/{sample_id}_2.fq.gz"
    threads: 15
    conda:
        config['conda']['cfDNA_base']    
    log:
        outdir + "/log/{genome}/{sample_id}/get-unaligned.txt"
    shell:
        """
        samtools fastq -@ {threads} -1 {output.unmapped_1} -2 {output.unmapped_2} -0 /dev/null -s /dev/null -f 13 {input.bam} > {log} 2>&1
        """
# Generate flag statistics in bam file
rule flagstat:
    input:
        bam = outdir + "/{genome}/bam/{sample_id}.bam"
    output:
        flagstat = outdir + "/{genome}/bam/flagstat/{sample_id}.txt"
    conda:
        config['conda']['cfDNA_base']
    log:
        outdir + "/log/{genome}/{sample_id}/flagstat.log"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """ 
# sort bam files
rule sort:
    input:
        bam = outdir + "/{genome}/bam/{sample_id}.bam"
    output:
        bam = temp(outdir + "/{genome}/bam-sorted/{sample_id}.bam"),
        bai = temp(outdir + "/{genome}/bam-sorted/{sample_id}.bam.bai")
    params:
        keep_proper_pair = "-f 2" if onlykeep_properpair else "",
    log:
        outdir + "/log/{genome}/{sample_id}/sort.log"
    threads: 4
    conda:
        config['conda']['cfDNA_base']
    shell:
        """
        samtools view -h {params.keep_proper_pair} -F 0x4  {input.bam} | samtools sort  -@ {threads}  > {output.bam} 2>{log}
        samtools index -@ {threads} {output.bam} 2>{log}
        """
# get bam files with RG (gatk requires this Read Group)
rule addReadsGroup:
    input:
        bam = outdir + "/{genome}/bam-sorted/{sample_id}.bam"
    output:
        bam = outdir + "/{genome}/RG/{sample_id}.bam",
        bai = outdir + "/{genome}/RG/{sample_id}.bam.bai"
    params:
        id = "{sample_id}",
        java = "--java-options -Xmx15G",
        tmp_dir = config["tmp_dir"]
    threads: 10
    conda:
        config['conda']['cfDNA_base']        
    log:
        log = outdir + "/log/{genome}/{sample_id}/addReadsGroup.log"
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
        bam = outdir + "/{genome}/RG/{sample_id}.bam"
    output:
        bam = outdir + "/{genome}/bam-sorted-deduped/{sample_id}.bam",
        bai = outdir + "/{genome}/bam-sorted-deduped/{sample_id}.bam.bai",
        metrics = outdir + "/{genome}/bam-sorted-deduped/{sample_id}_dedup-metrics.txt"
    log:
        outdir + "/log/{genome}/{sample_id}/MarkDuplicates.log"
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
##############################################get count matrix###########################################################
# get count matrix
rule count_matrix_gene:
    input:
        bam = expand( outdir + "/{genome}/bam-sorted-deduped/{sample_id}.bam",sample_id=samples,genome=genomes)
    output:
        gene_matrix = outdir + "/{genome}/matrix/gene/count_matrix_{region}.txt",
        gene_sum = outdir + "/{genome}/matrix/gene/count_matrix_{region}.txt.summary",
        gene_CPM_TMM = outdir + "/{genome}/matrix/gene/CPM-TMM_matrix_{region}.txt"
    log:
        log1 = outdir + "/{genome}/matrix/gene/log/count_matrix_{region}.log",
        CPM_TMM = outdir + "/{genome}/matrix/gene/log/CPM-TMM_matrix_{region}.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 25
    params:
        tmp = outdir + "/{genome}/matrix/gene/tmp_{region}",
        region = "{region}",
        gtf = config['count_matrix_gene']['gtf'],
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
        bam = expand(outdir + "/{genome}/bam-sorted-deduped/{sample_id}.bam",sample_id=samples,genome=genomes)
    output:
        gene_matrix = outdir + "/{genome}/matrix/TE/{level}/count_matrix_{region}.txt",
        gene_sum = outdir + "/{genome}/matrix/TE/{level}/count_matrix_{region}.txt.summary",
        gene_CPM_TMM = outdir + "/{genome}/matrix/TE/{level}/CPM-TMM_matrix_{region}.txt"
    log:
        log1 = outdir + "/{genome}/matrix/TE/{level}/log/count_matrix_{region}.log",
        CPM_TMM = outdir + "/{genome}/matrix/TE/{level}/log/CPM-TMM_matrix_{region}.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 25
    params:
        tmp = outdir + "/{genome}/matrix/TE/{level}/tmp_{region}",
        region = "{region}",
        level = "{level}",
        gtf = config['count_matrix_TE']['gtf'],
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
################################################# bam statistics #####################################################
# get WIG files (for IGV visualization)
rule wig:
    input:
        bam = outdir + "/{genome}/bam/{sample_id}.bam" if not remove_duplications else outdir+"/{genome}/bam-sorted-deduped/{sample_id}.bam"
    output:
        bigwig = outdir + "/{genome}/wig/{sample_id}.bigwig"
    log:
        log = outdir + "/log/{genome}/{sample_id}/wig.log"
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
        bam = outdir + "/{genome}/bam/{sample_id}.bam" if not getBamStatistics_duplications else outdir + "/{genome}/bam-sorted-deduped/{sample_id}.bam"
    output:
        summary =  outdir + "/{genome}/stats/{sample_id}/{sample_id}.txt" 
    conda:
        config['conda']['cfDNA_base']
    log:
        outdir + "/log/{genome}/{sample_id}/getBamStatistics.log"
    shell:
        """
        samtools stats {input.bam} > {output.summary} 2>{log}
        """
# Summarize histogram of fragment length from samtools stats output
rule getFragmentSize:
    input:
        stats = expand(outdir + "/{genome}/stats/{sample_id}/{sample_id}.txt"  ,sample_id=samples,genome=genomes)
    output:
        hist = outdir + "/{genome}/fragment/FragmentSize.txt",
        png = outdir + "/{genome}/fragment/FragmentSize.png"
    conda:
        config['conda']['GCBias']
    log:
        outdir + "/log/{genome}/getFragmentSize.log"
    shell:
        """
            python scripts/getFragmentSize.py --input {input.stats} --out {output.hist} > {log} 2>&1
            python scripts/plotFragmentSize.py --input {output.hist} --out {output.png} >> {log} 2>&1
        """
        

