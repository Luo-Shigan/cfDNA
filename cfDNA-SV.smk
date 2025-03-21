configfile: "config/cfDNA-SV.yaml"
################################################################
# Structural Variation (manta)
################################################################
## Config manta structural variation caller
samples = glob_wildcards("data/fq/{sample_id}_1.fq.gz").sample_id
outdir="output"
rule all:
    input:
        expand(outdir+"/struatural-variation/{sample_id}/results/variants/candidateSV.vcf.gz",sample_id=samples)
rule prepareMantaConfig:
    input:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam"
        # reference="genome/fasta/hg38.fa"
    output:
        config=outdir+"/struatural-variation/{sample_id}/runWorkflow.py.config.pickle",
        script=outdir+"/struatural-variation/{sample_id}/runWorkflow.py" 
    log: 
        outdir+"/log/{sample_id}/prepareMantaConfig.log"
    params:
        gn_fa_path=config['genome'],
        outdir=outdir
    conda:
        config['conda']['manta']
    shell:
        """
        configManta.py --bam {input.bam} --referenceFasta {params.gn_fa_path} --runDir "{params.outdir}/struatural-variation/{wildcards.sample_id}"  > {log} 2>&1
        """

## Call structure variation
rule runMantan:
    input:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam", 
        # reference="genome/fasta/hg38.fa",
        config=outdir+"/struatural-variation/{sample_id}/runWorkflow.py.config.pickle",
        script=outdir+"/struatural-variation/{sample_id}/runWorkflow.py"
    output:
        sv=outdir+"/struatural-variation/{sample_id}/results/variants/candidateSV.vcf.gz"
    threads: 16
    log: outdir+"/log/{sample_id}/runMantan.log"
    conda:
        config['conda']['manta']
    shell:
        """
        python {input.script} -j 6 > {log} 2>&1
        """