shell.prefix("set -x; set -e;")
configfile: "config/cfDNA-SNV.yaml"
outdir="output"
samples = glob_wildcards("data/fq/{sample_id}_1.fq.gz").sample_id
rule all:
    input:
        expand(outdir+"/vcf-filtered/{sample_id}.vcf.gz",sample_id=samples)
################################################################
# SNV (option1: GATK germline mode using Haplotypecaller)
################################################################
rule HaplotypeCaller:
    input:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam",
    output:
        vcf=outdir+"/vcf/{sample_id}.vcf.gz"
    log:
        outdir+"/log/{sample_id}/haplotypeCaller.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        gn_fa_path=config['genome'],
        tmp="temp",
        java="--java-options -Xmx35G" # "--java-options -Xmx10G", # if omit this option, a default value will be used.
    threads: 25
    shell:
        """
        gatk HaplotypeCaller {params.java} \
            -R {params.gn_fa_path} -I {input.bam} -O {output.vcf} \
            --tmp-dir {params.tmp} > {log} 2>&1
        """

rule filterHaplotypeCallerVcf:
    input:
        vcf=outdir+"/vcf/{sample_id}.vcf.gz"
    output:
        vcf=outdir+"/vcf-filtered/{sample_id}.vcf.gz"
    log:
        outdir+"/log/{sample_id}/haplotypeCaller-filtering.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 25
    params:
        gn_fa_path=config['genome'],
        java="--java-options -Xmx35G", # "--java-options -Xmx15G",
        tmp_dir="temp"
    shell:
        """
        gatk VariantFiltration {params.java}\
            -R {params.gn_fa_path} -V {input.vcf} -O {output.vcf} \
            -window 35 -cluster 3 \
            --filter-name FS20 -filter "FS > 20.0" \
            --filter-name QD2 -filter "QD < 2.0" \
            --filter-name DP10 -filter "DP < 10.0" \
            --filter-name QUAL20 -filter "QUAL < 20.0" \
            --tmp-dir {params.tmp_dir} > {log} 2>&1
        """ 


################################################################
# SNV (option2: GATK somatic mode using Mutect2)
################################################################
rule somaticMutect2:
    input:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam",
        gn_fa_path=config['genome'],
        gnomad_path=outdir+"/snpRef/somatic-hg38_af-only-gnomad.hg38.vcf.gz",
        genome1k_path=outdir+"/snpRef/somatic-hg38_1000g_pon.hg38.vcf.gz"
    output:
        vcf=outdir+"/mutect2-vcf/{sample_id}.vcf.gz"
    log:
        outdir+"/log/{sample_id}/mutect2.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        # tmp="tmp",
        # java="" # "--java-options -Xmx10G", # if omit this option, a default value will be used.
    threads: 16
    shell:
        """
        gatk Mutect2 \
            -R {input.gn_fa_path} -I {input.bam} -O {output.vcf} \
            --germline-resource {input.gnomad_path} \
            --panel-of-normals {input.genome1k_path} \
            --native-pair-hmm-threads 10 \
            --tmp-dir {tmp_dir} > {log} 2>&1
        """

rule filterMutect2Vcf:
    input:
        vcf=outdir+"/mutect2-vcf/{sample_id}.vcf.gz",
        gn_fa_path=config['genome']
    output:
        vcf=outdir+"/mutect2-vcf-filtered/{sample_id}.vcf.gz",
    log:
        outdir+"/log/{sample_id}/mutect2-filtering.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 16
    params:
        # java="--java-options -Xmx15G",
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.gn_fa_path} -V {input.vcf} -O {output.vcf} \
            -f-score-beta 1 \
            --tmp-dir {tmp_dir} > {log} 2>&1
        """ 