#

索引建立不包括在规则内，需要额外建立

## 开发日志

20250425:

- cfDNA-count_fragmentSize：新增小鼠支持

## 使用方法

运行：

```python
# 实际运行
python /path/to/run.py --pipeline /path/to/snakefile --indir /path/to/*fq.gz --outdir /path/to/outputDir
# 尝试运行
python /path/to/run.py --pipeline /path/to/snakefile --indir /path/to/*fq.gz --outdir /path/to/outputDir --dry-run
```

帮助：

```python
# run.py参数解析
python /path/to/run.py -h
# 具体模块参数解析
python /path/to/run.py --pipeline /path/to/snakefile --pipeHelp
```
