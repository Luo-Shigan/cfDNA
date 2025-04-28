import pandas as pd
import argparse
def download_diff(all:str,downloaded:str,outfile:str):
    """
    find which run ID hasn't been downloaded
    """
    dfAll = pd.read_csv(all,sep="\t",header=None)
    dfDownloaded = pd.read_csv(downloaded,sep="\t",header=None)
    set_all = set(dfAll[0])
    set_downloaded = set(dfDownloaded[0])
    undownloaded = set_all - set_downloaded
    with open(outfile, 'w') as f:
        for elem in undownloaded:
            f.write(str(elem) + '\n')
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to find which run ID hasn't been downloaded.")
    # 定义命令行参数 
    parser.add_argument('--all', type=str, required=True, help='Path to input file contain all runId')
    parser.add_argument('--downloaded', type=str, required=True, help='Path to input file contain downloaded runId')
    parser.add_argument('--output', type=str, required=True, help='Path to output file')
    args = parser.parse_args()
    all = args.all
    downloaded = args.downloaded
    output = args.output
    download_diff(all,downloaded,output)




