### Bioconductor上没有Gencode的FASTA PRI版本的，需要自己构建
library(BSgenomeForge)
seed_file="/ChIP_seq_2/Data/serum/cfDNA/scripts/BSgenome.seed"
#seqs_srcdir;destdir 序列文件所在以及输出的位置
forgeBSgenomeDataPkg(seed_file, seqs_srcdir=getwd(), destdir=getwd(), verbose=TRUE)
#forgeBSgenomeDataPkg(seed_file, verbose=TRUE)
#unlink参数表示是否overwrite已有的目录，seqs_srcdir是twoBit的目录，destdir为生成包的目录，这里需要一定的时候。