import pandas as pd
import argparse
### 从文件中提取指定列的值
# infile: 输入文件
# sep: 输入文件的分隔符
# colKey: 指定列的列名
# colValue: 指定列的值
# outfile: 输出文件
def keepcolmun(infile: str,sep: str,colKey: str,colValue: str,outfile: str):
    df = pd.read_table(infile,sep=sep)
    df = df[df[colKey] == colValue]
    df.to_csv(outfile,sep = "\t",index = False)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # 定义命令行参数 
    parser.add_argument('--input', type=str, required=True, help='Path to input file')
    parser.add_argument('--sep', type=str, required=True, help='the sep of your file is what ?')
    parser.add_argument('--colKey', type=str, required=True, help='the column name which you want to filter')
    parser.add_argument('--colValue', type=str, required=True, help='the column Value which you want to filter')
    parser.add_argument('--output', type=str, required=True, help='Directory to save output files')
    # 解析命令行参数 
    args = parser.parse_args()
    # 调用主函数，传递命令行参数 
    keepcolmun(infile=args.input, sep=args.sep, colKey= args.colKey, colValue = args.colValue ,outfile=args.output)
    
