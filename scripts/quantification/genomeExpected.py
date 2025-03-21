import pandas as pd

def calculate_relative_expression(gtf_file, genome_size, total_reads, output_csv):
    # 使用 Pandas 读取 GTF 文件，跳过注释行
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

    # 给列命名
    gtf_df.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

    # 提取 gene_id
    gtf_df['gene_id'] = gtf_df['attribute'].str.extract('gene_id "([^"]+)"')

    # 计算每个片段的大小
    gtf_df['size'] = gtf_df['end'] - gtf_df['start']

    # 按 gene_id 聚合，计算每个 gene_id 的片段大小总和
    gene_sizes = gtf_df.groupby('gene_id')['size'].sum().reset_index()

    # 计算相对表达量
    gene_sizes['relative_expression'] = (gene_sizes['size'] / genome_size) * total_reads

    # 计算所有基因的总片段大小
    total_size = gene_sizes['size'].sum()

    # CPM 标准化
    gene_sizes['CPM'] = (gene_sizes['size'] / total_size) * 1e6

    # 输出到 CSV 文件
    gene_sizes.to_csv(output_csv, index=False)

    print(f"Result saved to {output_csv}")

# 使用示例
gtf_file = "/ChIP_seq_2/StemCells/data/genome/GRCh38_GENCODE_rmsk_TE.gtf"  # GTF 文件路径
genome_size = 3200000000  # 基因组大小，单位为碱基对
total_reads = 998230554  # 测序总读数
output_csv = "relative_expression_with_CPM.csv"  # 输出 CSV 文件路径

calculate_relative_expression(gtf_file, genome_size, total_reads, output_csv)
