# run.py
import argparse
import subprocess
import os

def help(pipleline):
    if pipleline == "cfDNA-count_fragmentSize.smk":
        print("""
            可用配置项说明：
            ----------------
            indir          : 输入数据目录，默认 ../../data/fq
            outdir         : 输出结果目录，默认 ../../output
            onlykeep_properpair : 是否仅保留proper pair的reads, 默认 true
            getBamStatistics_duplications : 是否获取去重后的bam统计信息, 默认 false
            remove_duplications: 在可视化bam文件时, 是否采用去除PCR重复的bam文件, 默认 false

            实际调用示例：
            snakemake --cores [核心数] \
            --directroy [snakefile目录] \
            --config indir = [fq.gz目录] \
            outdir = [输出目录] \
            ……

            注意: 
            - 命令行传递的目录要么为绝对路径, 要么是相对于snakfile所在目录的相对路径, 而不是相对于工作目录的相对路径
            - 第一次运行建议加上--dry-run参数, 以检查是否有错误
            
            """)
        
    exit(0)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run Snakemake pipeline of cfDNA")
    # pipeline
    parser.add_argument("--pipeline", type=str,default="cfDNA-count_fragmentSize.smk",required=True, help="Pipeline path or name")
    parser.add_argument("--cores", type=str, default="30", help="Number of cores to use")
    # help
    parser.add_argument("--pipeHelp",action="store_true",help="Show help for the pipeline")
    parser.add_argument("--dry-run", action="store_true", help="Dry run the pipeline")
    # input and output directories
    parser.add_argument("--indir", type=str, default="../../data/fq", help="fq.gz Input directory")
    parser.add_argument("--outdir", type=str, default="../../output", help="Output directory")
    # cfDNA-count_fragmentSize.smk
    parser.add_argument("--onlykeep_properpair", type=bool, default=True, help="Whether to keep only proper pair reads")
    parser.add_argument("--getBamStatistics_duplications", type=bool, default=False, help="Whether to get bam statistics after deduplication")
    parser.add_argument("--remove_duplications", type=bool, default=False, help="Whether to use deduplicated bam file for visualization")
    args = parser.parse_args()
    pipeline = args.pipeline
    directory = os.path.dirname(pipeline)
    cores = args.cores

    pipelineName = os.path.basename(pipeline)
    if pipelineName == "cfDNA-count_fragmentSize.smk":
        if args.pipeHelp:
            help(pipelineName)

        indir = args.indir
        outdir = args.outdir
        onlykeep_properpair = args.onlykeep_properpair
        getBamStatistics_duplications = args.getBamStatistics_duplications
        remove_duplications = args.remove_duplications
        indir = os.path.relpath(indir, start=directory) #路径末尾不会携带/，不用额外处理
        outdir = os.path.relpath(outdir, start=directory)
    # 构造 snakemake 命令
    if args.dry_run:
        cmd = [
            "snakemake",
            "--snakefile", pipeline,
            "--cores", cores,
            "--directory",directory,
            "--config",
            f"indir={indir}",
            f"outdir={outdir}",
            f"onlykeep_properpair={onlykeep_properpair}",
            f"getBamStatistics_duplications={getBamStatistics_duplications}",
            f"remove_duplications={remove_duplications}",
            "--use-conda",
            "--dry-run"
        ]
    else:
        cmd = [
            "snakemake",
            "--snakefile", pipeline,
            "--cores", cores,
            "--directory",directory,
            "--config",
            f"indir={indir}",
            f"outdir={outdir}",
            f"onlykeep_properpair={onlykeep_properpair}",
            f"getBamStatistics_duplications={getBamStatistics_duplications}",
            f"remove_duplications={remove_duplications}",
            "--use-conda"
        ]

    # 执行
    print("Executing command:", " ".join(cmd))
    # subprocess.run(cmd)
