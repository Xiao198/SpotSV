import sys
import logging
import argparse
import os
from valid.multiprocess_input import vcf_run

def run():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='consists of two major steps')
    subparsers = parser.add_subparsers(title='modules', dest='command', )  # two submodules

    parser_vcf = subparsers.add_parser('valid', help='if input structural variants are encoded in vcf format')
    parser_vcf.add_argument('-s', dest="sv_input", help='path to input file in bed or vcf format', type=os.path.abspath,
                            required=True)
    parser_vcf.add_argument("-o", dest="out_path", help="path to output", type=os.path.abspath, required=True)
    parser_vcf.add_argument("-r", dest="ref_path", help="path to ref", type=os.path.abspath, required=True)
    parser_vcf.add_argument('-b', dest="bam_input", help='path to bam', type=os.path.abspath, required=True)
    parser_vcf.add_argument('-t', dest="threads", type=int, help='thread number')
    parser_vcf.add_argument('-p', dest="process_images", help="whether to output processing images", default=False,
                            action='store_true')
    parser_vcf.add_argument('--data_type', help="data type, choose from CCS,ONT and CLR", choices=['CCS', 'ONT', 'CLR'],
                            required=True)
    args = parser.parse_args()

    if args.command == 'valid':
        print(args)
        if args.data_type == 'CCS':
            kmer_length = 31
        elif args.data_type in ['ONT', 'CLR']:
            kmer_length = 17
        else:
            kmer_length = 31
        vcf_run(args.threads, args.sv_input, args.out_path, args.ref_path, args.bam_input, args.process_images,
                kmer_length)


script_name = sys.argv[0]
if len(sys.argv) < 2:
    print('=======================================================')
    print('SpotSV: Validate SVs through vcf format files\n')
    print('Usage:')
    print('python spotSV.py [commands] <parameters>\n')
    print('commands:')
    print('valid: Validate SVs through vcf format files')
    print("=======================================================")
else:
    run()
