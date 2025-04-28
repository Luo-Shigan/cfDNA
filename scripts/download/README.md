# Description

all is done

## 开发日志

ascp.sh 加上 sleep 10s;10s对于下载时间的影响忽略不计，确能避免ascp连续下载突然降速或者下载

2025.04.25:
    - 增加gzipTest.sh 判断ascp下载文件是否压缩正确
    - 增加diff.py 判断那些文件没有下载
