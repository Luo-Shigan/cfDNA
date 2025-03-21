shell.prefix("set -x; set -e;")
configfile: "config/cfDNA-FragSize.yaml"
samples = glob_wildcards("data/fq/{sample_id}_1.fq.gz").sample_id
outdir="output"
rule all:
    input:
        outdir+"/shortLong/matrix/DNA-FragRatio_matrix_gene.txt"
################################################################
# Frag Size Ratio (gene-based) (2019 Nature DELFI)
################################################################
### get short&long bam

rule prepareShortLongBamGene:
    input:
        bam= outdir + "/GC/{sample_id}.bam", 
        # reference="genome/fasta/hg38.fa"
    output:
        shortBam=outdir+"/shortLong/{sample_id}-short.bam",
        longBam=outdir+"/shortLong/{sample_id}-long.bam",
    # log: "{outdir}/log/{sample_id}/prepareShortLongBam.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        outdir=outdir+"/shortLong"
    shell:
        """
        bash scripts/frag/get-FragmentSizeBam.sh \
            {wildcards.sample_id} \
            {input.bam} \
            {params.outdir}
        """


### count short bam
rule count_short_matrix_gene:
    input:
        bam=expand(outdir + "/shortLong/{sample_id}-short.bam",sample_id=samples),
    output:
        gene_matrix=outdir+"/shortLong/matrix/count_matrix_gene_short.txt",
        gene_sum=outdir+"/shortLong/matrix/count_matrix_gene_short.txt.summary",
        gene_CPM_TMM=outdir+"/shortLong/matrix/CPM-TMM_matrix_gene_short.correctGC.txt"
    log:
        log1=outdir+"/shortLong/matrix/log/count_matrix_gene_short.log",
        CPM_TMM=outdir+"/shortLong/matrix/log/CPM-TMM_matrix_gene_short.log",
    conda:
        config['conda']['cfDNA_base']
    threads: 8
    params:
        tmp=temp(outdir+"/shortLong/matrix/tmp_gene_short"),
        region="gene",
        gtf=config['gtf']
    shell:
        """
        featureCounts -T 4 -O -t {params.region} -g gene_id -M -p \
            -a {params.gtf} \
            -o {output.gene_matrix} {input.bam} \
            > {log.log1} 2>&1
        
        /usr/bin/Rscript scripts/multiFeatureCounts2countMat.R \
            -i {output.gene_matrix}  \
            -o {params.tmp} ;
        mv {params.tmp} {output.gene_matrix};


        /usr/bin/Rscript scripts/run-NormCountMat.R \
            -i {output.gene_matrix} \
            -o {output.gene_CPM_TMM} \
            -m TMM \
            > {log.CPM_TMM} 2>&1
        """ 

### count long bam
rule count_long_matrix_gene:
    input:
        bam=expand(outdir+"/shortLong/{sample_id}-long.bam",sample_id=samples)
    output:
        gene_matrix=outdir+"/shortLong/matrix/count_matrix_gene_long.txt",
        gene_sum=outdir+"/shortLong/matrix/count_matrix_gene_long.txt.summary",
        gene_CPM_TMM=outdir+"/shortLong/matrix/CPM-TMM_matrix_gene_long.correctGC.txt"
    log:
        log1=outdir+"/shortLong/matrix/log/count_matrix_gene_long.log",
        CPM_TMM=outdir+"/shortLong/matrix/log/CPM-TMM_matrix_gene_long.log",
    conda:
        config['conda']['cfDNA_base']
    threads: 8
    params:
        tmp=temp(outdir+"/shortLong/matrix/tmp_gene_long"),
        region="gene",
        gtf1=config['gtf']
    shell:
        """
        featureCounts -T 4 -O -t {params.region} -g gene_id -M -p \
            -a {params.gtf1} \
            -o {output.gene_matrix} {input.bam} \
            > {log.log1} 2>&1
        
        /usr/bin/Rscript scripts/multiFeatureCounts2countMat.R \
            -i {output.gene_matrix}  \
            -o {params.tmp} ;
        mv {params.tmp} {output.gene_matrix};

        /usr/bin/Rscript scripts/run-NormCountMat.R \
            -i {output.gene_matrix} \
            -o {output.gene_CPM_TMM} \
            -m TMM \
            > {log.CPM_TMM} 2>&1
        """ 

### cal FragSize ratio from correctGC CPM matrix
rule calFragSizeRatio:
    input:
        shortMat=outdir+"/shortLong/matrix/CPM-TMM_matrix_gene_short.correctGC.txt",
        longMat=outdir+"/shortLong/matrix/CPM-TMM_matrix_gene_long.correctGC.txt"
    output:
        outdir+"/shortLong/matrix/DNA-FragRatio_matrix_gene.txt",
    log: outdir+"/shortLong/matrix/log/calFragSizeRatio.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 8
    shell:
        """
        /usr/bin/Rscript scripts/frag/get-FragSizeMat.R \
            --short {input.shortMat} \
            --long {input.longMat} \
            --outfile {output} > {log} 2>&1
        """
