import re
import pandas as pd
import argparse
def FragementHist(inlist:list,outfile:str):
    records=[]
    for path in inlist:
        sample_id=re.sub(".txt$","",path.split("/")[-1])
        with open(path) as f:
            for line in f:
                if not line.startswith("IS"):
                    continue
                #IS, insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs
                fields=line.strip().split("\t")
                IS, N=fields[1],fields[3]
                records.append((sample_id,int(IS),int(N)))
    fsHistTable=pd.DataFrame.from_records(records)
    fsHistTable.columns=["sample_id","insertion-size","count"]
    fsHist=fsHistTable.pivot(index="insertion-size",columns="sample_id",values="count")
    fsHist.to_csv(outfile,sep="\t")
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to read scRNA-seq data.")
    parser.add_argument('--input', type=str, nargs='+', required=True, help='List of input file paths')
    parser.add_argument('--out', type=str, required=True, help='Path to out file')
    args = parser.parse_args()
    FragementHist(args.input,args.out)