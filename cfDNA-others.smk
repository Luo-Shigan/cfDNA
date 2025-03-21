########################################################################################
SNP_dir="output/SNP"
BQSR_dir="output/BQSR"       
## get BQSR bam files (correct single base quality)
rule baseQualityRecalibration:
    input:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam",
        gn_fa_path=config["genome"],
        dbSNP_path=SNP_dir+"/All_20180418.vcf.gz"
    output: 
        bam=BQSR_dir+"/{sample_id}.bam",
        bai=BQSR_dir+"/{sample_id}.bam.bai",
        grpFile=BQSR_dir+"/BQSR-grp/{sample_id}.grp"
    conda:
        config['conda']['cfDNA_base']
    log:
        prepareGrp=outdir+"/log/{sample_id}/BaseRecalibrator.log",
        BQSR=outdir+"/log/{sample_id}/ApplyBQSR.log"
    threads: 16
    # params:
    #     tmp_dir=config['tmp_dir']
    shell:
        """
        gatk BaseRecalibrator \
            --tmp-dir {tmp_dir} --input {input.bam} --output {output.grpFile} \
            --known-sites {input.dbSNP_path} --reference {input.gn_fa_path} > {log.prepareGrp} 2>&1
        gatk ApplyBQSR \
            --tmp-dir {tmp_dir} -R {input.gn_fa_path} -I {input.bam}  \
            --bqsr-recal-file {output.grpFile} -O {output.bam} > {log.BQSR} 2>&1
        samtool index -@ {threads} {output.bam}
        """
########################################################################################
# get correctGC bam files (correct GC bias for CNV )
## computeGC
rule fa2bit:
    input:
        fasta=config["genome"]
    output:
        gn_2bit_path = index_dir +"/fa2bit/GRCh38.2bit"
    conda:
        config['conda']['GCBias']
    log:
        index_dir + "/fa2bit/GRCh38.log"
    shell:
        """
        faToTwoBit {input.fasta} {output.gn_2bit_path} > {log} 2>&1
        """
rule computeGCBias:
    input:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam",
        gn_2bit_path=index_dir +"/fa2bit/GRCh38.2bit"
    output: 
        freq=outdir+"/GC_bias/{sample_id}.txt",
        png=outdir+"/GC_bias/{sample_id}.png"
    log:
        outdir+"/GC_bias/log/{sample_id}.log"
    conda:
        config['conda']['GCBias']
    threads: 16
    params:
        gn_size = config['GC']['genome_size'][1]['GRCh38']
    shell:
        """
        computeGCBias \
            --numberOfProcessors {threads} \
            --effectiveGenomeSize {params.gn_size} \
            --biasPlot {output.png} \
            --plotFileFormat png \
            --genome {input.gn_2bit_path} \
            --bamfile {input.bam} \
            --GCbiasFrequenciesFile {output.freq} \
            > {log} 2>&1
        """
## correctGC
rule correctGCBias:
    input:
        bam=outdir+"/bam-sorted-deduped/{sample_id}.bam",
        freq=outdir+"/GC_bias/{sample_id}.txt",
        gn_2bit_path=outdir+"/gn/hg38.2bit"
    output: 
        bam=outdir+"/GC/{sample_id}.bam",
        bai=outdir+"/GC/{sample_id}.bam.bai"
    log:
        outdir+"/log/{sample_id}.log"
    threads: 16
    conda:
        config['conda']['GCBias']
    params:
        gn_size = config['GC']['genome_size'][1]['GRCh38']
    shell:
        """
        correctGCBias \
            --numberOfProcessors {threads} \
            --effectiveGenomeSize {params.gn_size} \
            --genome {input.gn_2bit_path} \
            --bamfile {input.bam} \
            --GCbiasFrequenciesFile {input.freq} \
            --correctedFile {output.bam} \
            > {log} 2>&1
        """
########################################################################################
# sum QC, this rule didn't work
rule wgs_qc:
    input:
        bam=outdir+"/bam/{sample_id}.bam" if not config["remove_duplications"] else outdir+"/bam-sorted-deduped/{sample_id}.bam"
    output:
        saturation_qc=outdir+"/wgs_qc/saturation/{sample_id}.txt.300",
        saturations_30k=outdir+"/wgs_qc/saturation/{sample_id}.txt.3000",
        saturation_300k=outdir+"/wgs_qc/saturation/{sample_id}.txt.30000",
        enrich_qc=outdir+"/wgs_qc/enrich/{sample_id}.txt",
    log:
        outdir+"/wgs_qc/log/{sample_id}_wgsqc.log"
    params:
        outdir+"/wgs_qc/saturation/{sample_id}.txt"
    shell:
        """
        Rscript scripts/wgs_qc.R -i {input.bam} -sr {params} -er {output.enrich_qc} \
            > {log} 2>&1
        """
# this rule didn't work
rule summary_wgs_qc:
    input:
        enrichments=expand("{outdir}/wgs_qc/enrich/{sample_id}.txt",outdir=outdir,sample_id=samples),
        saturations=expand("{outdir}/wgs_qc/saturation/{sample_id}.txt.300",outdir=outdir,sample_id=samples),
        saturations_30k=expand("{outdir}/wgs_qc/saturation/{sample_id}.txt.3000",outdir=outdir,sample_id=samples),
        saturation_3M=expand("{outdir}/wgs_qc/saturation/{sample_id}.txt.30000",outdir=outdir,sample_id=samples),
    output:
        summary="{outdir}/wgs_qc/quality-control.txt"
    run:
        import pandas as pd
        sample_ids=[path.split("/")[-1].split(".")[0] for path in input.enrichments]
        records=[]
        for sample_id in sample_ids:
            enrichment_path=wildcards.outdir + "/wgs_qc/enrich/{}.txt".format(sample_id)
            with open(enrichment_path) as f:
                for line in f:
                    key,value=line.strip().split("\t")
                    if key == "enrichment.score.relH":
                        relH=value
                    elif key == "enrichment.score.GoGe":
                        GoGe=value
            saturation_path=wildcards.outdir + "/wgs_qc/saturation/{}.txt.300".format(sample_id)
            sat_df=pd.read_csv(saturation_path,sep="\t")
            es_sat_df=sat_df[sat_df["data"]=="estimated"]
            estimated_saturation=es_sat_df.sort_values(by="subset")["correlation"].iloc[-1]
            ob_sat_df=sat_df[sat_df["data"]=="observed"]
            observed_saturation=ob_sat_df.sort_values(by="subset")["correlation"].iloc[-1]

            saturation_path1=wildcards.outdir + "/wgs_qc/saturation/{}.txt.3000".format(sample_id)
            sat_df1=pd.read_csv(saturation_path1,sep="\t")
            es_sat_df1=sat_df1[sat_df1["data"]=="estimated"]
            estimated_saturation1=es_sat_df1.sort_values(by="subset")["correlation"].iloc[-1]
            ob_sat_df1=sat_df1[sat_df1["data"]=="observed"]
            observed_saturation1=ob_sat_df1.sort_values(by="subset")["correlation"].iloc[-1]

            saturation_path2=wildcards.outdir + "/wgs_qc/saturation/{}.txt.30000".format(sample_id)
            sat_df2=pd.read_csv(saturation_path2,sep="\t")
            es_sat_df2=sat_df2[sat_df2["data"]=="estimated"]
            estimated_saturation2=es_sat_df2.sort_values(by="subset")["correlation"].iloc[-1]
            ob_sat_df2=sat_df2[sat_df2["data"]=="observed"]
            observed_saturation2=ob_sat_df2.sort_values(by="subset")["correlation"].iloc[-1]

            records.append((sample_id,relH,GoGe,estimated_saturation,observed_saturation,estimated_saturation1,observed_saturation1,estimated_saturation2,observed_saturation2))
        table=pd.DataFrame.from_records(records)
        table.columns=["sample_id","enrichment.score.relH","enrichment.score.GoGe","estimated.max.saturation","observed.max.saturation","estimated.max.saturation.3k","observed.max.saturation.3k","estimated.max.saturation.30k","observed.max.saturation.30k"]
        table.to_csv(output.summary,sep="\t",index=False)
# get (original) bam AlignmentSummaryMetrics, this rule wasn't test
rule AlignmentSummaryMetrics:
    input:
        bam=outdir+"/bam-sorted/{sample_id}.bam",
        gn_fa_path=config['genome']
    output:
        outdir+"/firstAlignmentSummary/{sample_id}.AlignmentSummaryMetrics.txt"
    log:
        outdir+"/firstAlignmentSummary/log/{sample_id}.AlignmentSummaryMetrics.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        java="--java-options -Xmx10G"
    shell:
        """
        gatk {params.java} CollectAlignmentSummaryMetrics \
            -R {input.gn_fa_path} \
            -I {input.bam} \
            -O {output} \
            > {log} 2>&1
        """