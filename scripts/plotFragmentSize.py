import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
def line(x:str = 'insertion-size',y:str=None,ax: plt.Axes=None):
    """
    在传入的 Axes 对象上绘制数据。
    :param ax: matplotlib 的 Axes 对象
    :param df: pandas DataFrame，包含数据
    """
    ax.plot(df[x], df[y], label=y)
    return ax
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to plot fragement size.")
    parser.add_argument('--input', type=str, required=True, help='Path to input file')
    parser.add_argument('--output', type=str, required=True, help='Path to output file')
    args = parser.parse_args()
    df = pd.read_table(args.input)
    fig, ax = plt.subplots(figsize=(25, 10), dpi=300)
    for i in df.columns:
        if i == 'insertion-size':
            continue
        else:
            sum = df[i].sum()
            df[i] = df[i]/sum
            line(y=i,ax=ax)
    # 获取当前 x 轴的刻度位置和标签
    current_xticks = ax.get_xticks()
    current_xticklabels = ax.get_xticklabels()

    # 定义要添加的新刻度位置和标签
    new_tick_position = 167
    new_tick_label = '167'

    # 将新刻度位置和标签添加到现有列表中
    new_xticks = np.append(current_xticks, new_tick_position)
    new_xticklabels = [label.get_text() for label in current_xticklabels]
    new_xticklabels.append(new_tick_label)

    # 重新设置 x 轴的刻度位置和标签
    ax.set_xticks(new_xticks)
    ax.set_xticklabels(new_xticklabels)
    ax.axvline(x=new_tick_position, color='r', linestyle='--', alpha=0.7)
    ax.set_xlabel('Fragment size (bp)')
    ax.set_ylabel('Frequency')
    # ax.set_title('insertion-size line plot')
    ax.legend()
    fig.savefig(args.output)