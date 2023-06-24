# SpotSV：a SV Validation Tool
### To Start with:
cd SpotSV/src

### to get the help tips
python spotsv.py valid -h

### to valid SV
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
