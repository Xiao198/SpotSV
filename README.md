# SpotSV：a SV Validation Tool
### 运行环境
Linux

### 所需依赖包
python3
pysam
bcftools
matplotlib
os
math
random
time
signal
shutil

### 进入src文件夹
cd SpotSV/src

### 运行提示命令
python spotsv.py valid -h

### 根据提示命令输入文件即可运行
usage: spotsv.py valid [-h] -s SV_INPUT -o OUT_PATH -r REF_PATH -b BAM_INPUT [-t THREADS] [-p] --data_type {CCS,ONT,CLR}

optional arguments:
  -h, --help            show this help message and exit
  -s SV_INPUT           path to input file in bed or vcf format
  -o OUT_PATH           path to output
  -r REF_PATH           path to ref
  -b BAM_INPUT          path to bam
  -t THREADS            thread number
  -p                    whether to output processing images
  --data_type {CCS,ONT,CLR}
                        data type, choose from CCS,ONT and CLR
