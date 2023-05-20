import os
import math
import shutil
import pysam
import random
import time
import signal
import multiprocessing
import matplotlib.pyplot as plt
from valid.sv_validation import run_validation
from valid.sequence import Sequence
from valid.alignment import Alignment
from valid.segment import Segment

def handler(signum, frame):
    raise AssertionError

def set_dir(filepath):
    '''
    如果文件夹不存在就创建，如果文件存在就清空！
    :param filepath:需要创建的文件夹路径
    :return:
    '''
    if not os.path.exists(filepath):
        os.mkdir(filepath)
    else:
        shutil.rmtree(filepath)
        os.mkdir(filepath)

def set_flank_length():
    return 1000

def analyze_sv_type(sv_str):
    '''
    analyze sv type though a sv type string.
    :param sv_str:
    :return:
    '''
    sv_type_list = []
    sv_str = sv_str.upper()
    # print(sv_str)
    if 'DEL' in sv_str:
        sv_type_list.append('DEL')
    if 'INS' in sv_str:
        sv_type_list.append('INS')
    if 'INV' in sv_str:
        sv_type_list.append('INV')
    if 'TDUP' in sv_str:
        sv_type_list.append('tDUP')
    if 'IDDUP' in sv_str:
        sv_type_list.append('idDUP')
    elif 'DDUP' in sv_str:
        sv_type_list.append('dDUP')
    return sv_type_list

def random_seq(length):
    '''
    generate random sequence
    :param length:
    :return:
    '''
    return ''.join(random.choice('ATCG') for x in range(length))

def sub_plot(seq1, seq2, segments, fig_pos, sub_title):
    '''
    subplot of ref vs ref, ref vs read, read vs read
    :param segments: a segment object list
    :param sub_num: such as (221)
    :param sub_title: title
    :return:
    '''
    # 遍历segment[]，获取每个segment对象，调用to_string方法
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    font = {
        'weight': 'normal',
        'size': 10,
    }
    for i in range(len(segments)):
        # x:read  y:ref
        read_coord = segments[i].get_x_coord()
        ref_coord = segments[i].get_y_coord()
        plt.subplot(fig_pos)
        plt.title(sub_title, font, y=-0.2)
        if segments[i].get_forward():
            plt.plot(eval(ref_coord), eval(read_coord), color='b', linewidth='2')
        else:
            plt.plot(eval(ref_coord), eval(read_coord), color='r', linewidth='2')
        plt.axis('square')
        ax = plt.gca()

        ax.spines['right'].set_visible(False)  # 右边框隐藏
        ax.spines['bottom'].set_visible(False)  # 下边框隐藏
        if fig_pos == 231:
            plt.axis([0, len_seq1, 0, len_seq1])
            ax.set_ylabel('REF', font)
            ax.set_xlabel('REF', font)
            ax.xaxis.set_label_coords(0.5, 1.05)
        elif fig_pos == 232:
            plt.axis([0, len_seq1, 0, len_seq2])
            ax.set_ylabel('READ', font)
            ax.set_xlabel('REF', font)
            ax.xaxis.set_label_coords(0.5, 1.05)
        elif fig_pos == 233:
            plt.axis([0, len_seq1, 0, len_seq2])
            ax.set_ylabel('READ', font)
            ax.set_xlabel('REF', font)
            ax.xaxis.set_label_coords(0.5, 1.05)

        elif fig_pos == 234:
            plt.axis([0, len_seq1, 0, len_seq1])
            ax.set_ylabel('ALT', font)
            ax.set_xlabel('ALT', font)
            ax.xaxis.set_label_coords(0.5, 1.05)
        elif fig_pos == 235:
            plt.axis([0, len_seq1, 0, len_seq2])
            ax.set_ylabel('READ', font)
            ax.set_xlabel('ALT', font)
            ax.xaxis.set_label_coords(0.5, 1.05)
        elif fig_pos == 236:
            plt.axis([0, len_seq1, 0, len_seq2])
            ax.set_ylabel('READ', font)
            ax.set_xlabel('ALT', font)
            ax.xaxis.set_label_coords(0.5, 1.05)
        elif fig_pos == 111:
            plt.axis([0, len_seq1, 0, len_seq2])
            ax.set_ylabel('READ', font)
            ax.set_xlabel('REF', font)
            # ax.xaxis.set_label_coords(0.5,1.05)
        ax.invert_yaxis()

def segment_plot(ref, read, alt, segments_ref_ref, segments_ref_read, filtered_ref_read_sorted, segments_alt_alt, segments_alt_read, filtered_alt_read_sorted, process_images):
    plt.clf()
    if process_images:
        sub_plot(ref, ref, segments_ref_ref, 231, 'REF vs REF')
        sub_plot(ref, read, segments_ref_read, 232, 'REF vs READ')
        sub_plot(ref, read, filtered_ref_read_sorted, 233, 'FILTERED REF vs READ')
        sub_plot(alt, alt, segments_alt_alt, 234, 'ALT vs ALT')
        sub_plot(alt, read, segments_alt_read, 235, 'ALT vs READ')
        sub_plot(alt, read, filtered_alt_read_sorted, 236, 'FILTERED ALT vs READ')
    else:
        sub_plot(ref, read, filtered_ref_read_sorted, 111, 'FILTERED REF vs READ')

def sig_vcf_process(ref_path, bam_input, chrom, start_coord, end_coord, sv_str, sv_length, ins_coord, vcf_sv, out_fig, process_images, vcf_record, kmer_length):
    print(os.getpid(), ' : ', start_coord, '-', end_coord)
    sv_type_list = analyze_sv_type(sv_str)
    flank_length = set_flank_length()
    reads_score_list = []
    bkps = []  # 存放每一个read的Breakpoint，[991, 2038]，如果是DEL,INS第一个最大，第二个最小
    plot_index_dict = {}  #{final_score: j,...}
    plot_tuple_dict = {}   #{j:plot_tuple,...}

    # read bam file
    bam = bam_input
    bam_reads, left_flank_length, right_flank_length = get_bam_reads(bam, chrom, int(start_coord) - flank_length, int(end_coord) + flank_length, kmer_length)
    # read fasta file
    fasta = pysam.FastaFile(ref_path)
    if len(bam_reads) > 0:
        if process_images:
            plt.figure(1, figsize=(13, 8))
        else:
            plt.figure(1, figsize=(5, 5))
        for j in range(len(bam_reads)):
            # if bam_reads[j][1] in ['S1_111878']:
            #     print(bam_reads[j][1])
            ref_seq = fasta.fetch(chrom, int(start_coord) - left_flank_length[j] + 1, int(end_coord) + right_flank_length[j])
            read_seq = bam_reads[j][0]
            read_id = bam_reads[j][-1]
            time_out_flag = 0

            if sv_type_list in [['DEL']]:
                alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[-right_flank_length[j]:]
            elif sv_type_list in [['INS']]:
                alt_seq = ref_seq[:left_flank_length[j]] + random_seq(int(sv_length)) + ref_seq[-right_flank_length[j]:]
            elif sv_type_list in [['INV']]:
                alt_seq = ref_seq[:left_flank_length[j]] + Sequence(ref_seq[left_flank_length[j]:(-right_flank_length[j])]).get_reverse_complement() + ref_seq[-right_flank_length[j]:]
            elif sv_type_list in [['tDUP']]:
                alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[left_flank_length[j]:(-right_flank_length[j])] + ref_seq[left_flank_length[j]:(-right_flank_length[j])] + ref_seq[-right_flank_length[j]:]
            elif sv_type_list in [['dDUP']]:
                if ins_coord == start_coord:
                    alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[-right_flank_length[j] - sv_length:-right_flank_length[j]] + ref_seq[left_flank_length[j]:(-right_flank_length[j])] + ref_seq[-right_flank_length[j]:]
                elif ins_coord == end_coord - sv_length:
                    alt_seq = ref_seq[:(-right_flank_length[j] - sv_length)] + ref_seq[left_flank_length[j]: left_flank_length[j] + sv_length] + ref_seq[-right_flank_length[j] - sv_length:]
                else:
                    alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[left_flank_length[j]:left_flank_length[j] + ins_coord - start_coord] + ref_seq[left_flank_length[j]:(-right_flank_length[j])] + ref_seq[left_flank_length[j] + ins_coord - start_coord: -right_flank_length[j]] + ref_seq[-right_flank_length[j]:]
            elif sv_type_list in [['idDUP']]:
                if ins_coord == start_coord:
                    alt_seq = ref_seq[:left_flank_length[j]] + Sequence(ref_seq[-right_flank_length[j] - sv_length:-right_flank_length[j]]).get_reverse_complement() + ref_seq[left_flank_length[j]:(-right_flank_length[j])] + ref_seq[-right_flank_length[j]:]
                elif ins_coord == end_coord - sv_length:
                    alt_seq = ref_seq[:(-right_flank_length[j] - sv_length)] + Sequence(ref_seq[left_flank_length[j]: left_flank_length[j] + sv_length]).get_reverse_complement() + ref_seq[-right_flank_length[j] - sv_length:]
                else:
                    alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[left_flank_length[j]:left_flank_length[j] + ins_coord - start_coord] + Sequence(ref_seq[left_flank_length[j]:(-right_flank_length[j])]).get_reverse_complement() + ref_seq[left_flank_length[j] + ins_coord - start_coord: -right_flank_length[j]] + ref_seq[-right_flank_length[j]:]

            elif sv_type_list in [['DEL', 'dDUP']]:
                if ins_coord[2] == start_coord:
                    alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[left_flank_length[j]:(-right_flank_length[j] + ins_coord[2] - ins_coord[3] + ins_coord[0] - ins_coord[1])] + ref_seq[left_flank_length[j]:left_flank_length[j] + ins_coord[3] - ins_coord[2]] + ref_seq[(-right_flank_length[j] + ins_coord[2] - ins_coord[3]):]
                elif ins_coord[0] == start_coord:
                    alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[(-right_flank_length[j] + ins_coord[2] - ins_coord[3]): -right_flank_length[j]] + ref_seq[left_flank_length[j] + ins_coord[1] - ins_coord[0]:]

            elif sv_type_list in [['DEL', 'idDUP']]:
                if ins_coord[2] == start_coord:
                    alt_seq = ref_seq[:left_flank_length[j]] + ref_seq[left_flank_length[j]:(-right_flank_length[j] + ins_coord[2] - ins_coord[3] + ins_coord[0] - ins_coord[1])] + Sequence(ref_seq[left_flank_length[j]:left_flank_length[j] + ins_coord[3] - ins_coord[2]]).get_reverse_complement() + ref_seq[(-right_flank_length[j] + ins_coord[2] - ins_coord[3]):]
                elif ins_coord[0] == start_coord:
                    alt_seq = ref_seq[:left_flank_length[j]] + Sequence(ref_seq[(-right_flank_length[j] + ins_coord[2] - ins_coord[3]): -right_flank_length[j]]).get_reverse_complement() + ref_seq[left_flank_length[j] + ins_coord[1] - ins_coord[0]:]
                else:
                    pass

            elif sv_type_list in [['DEL', 'INV']]:
                # del-inv
                if ins_coord[0] == start_coord:
                    alt_seq = ref_seq[:left_flank_length[j]] + Sequence(ref_seq[left_flank_length[j] + ins_coord[2] - start_coord : -right_flank_length[j]]).get_reverse_complement() + ref_seq[-right_flank_length[j]:]
                else:
                # inv-del
                    alt_seq = ref_seq[:left_flank_length[j]] + Sequence(ref_seq[left_flank_length[j]:left_flank_length[j] + ins_coord[1] - start_coord]) + ref_seq[-right_flank_length[j]:]

            try:
                signal.signal(signal.SIGALRM, handler)
                signal.alarm(600)
                score_ref_read, score_alt_read, sig_bkps, plot_tuple = run_validation(ref_seq, read_seq, alt_seq, sv_type_list, int(start_coord) - left_flank_length[j] + 1, sv_length, left_flank_length[j],right_flank_length[j], process_images, kmer_length)
                if sig_bkps != []:
                    bkps.append(sig_bkps)
            except AssertionError:
                # 超时的话，len(score_alt_read_list) == 0
                print(read_id, "time out")
                time_out_flag = 1
                break
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, signal.SIG_DFL)

            try:
                sv_score = 1 - score_alt_read / score_ref_read
            # 后续修改，sv score负值归为零
            except ZeroDivisionError:
                sv_score = 'NA'

            if sv_score == 'NA':
                final_score = 0
            else:
                final_score = float(sv_score)
                # print('final score:', final_score)
            if time_out_flag == 0:
                plot_index_dict[j] = final_score
                plot_tuple_dict[j] = plot_tuple
                if j == len(bam_reads) - 1:
                    max_plot_index = max(plot_index_dict, key=plot_index_dict.get)
                    segment_plot(plot_tuple_dict[max_plot_index][0], plot_tuple_dict[max_plot_index][1], plot_tuple_dict[max_plot_index][2],plot_tuple_dict[max_plot_index][3],plot_tuple_dict[max_plot_index][4],plot_tuple_dict[max_plot_index][5], plot_tuple_dict[max_plot_index][6],plot_tuple_dict[max_plot_index][7],plot_tuple_dict[max_plot_index][8],plot_tuple_dict[max_plot_index][9])
                    output_figname = os.path.join(out_fig, chrom + '_' + str(int(start_coord) - left_flank_length[max_plot_index] - 1) + '-' + str(int(end_coord) + right_flank_length[max_plot_index]) + '_' + sv_str + str(
                        plot_index_dict[max_plot_index]) + '.jpg')
                    plt.suptitle(chrom + ':' + str(int(start_coord) - left_flank_length[max_plot_index] - 1) + '-' + str(int(end_coord) + right_flank_length[max_plot_index]) + '_' + bam_reads[max_plot_index][-1] + '_' + sv_str)
                    plt.savefig(output_figname, bbox_inches='tight', pad_inches=0)

            reads_score_list.append(final_score)
    select_bkps_min = {}
    select_bkps_max = {}
    if bkps != []:
        for i in range(len(bkps)):
            if bkps[i][0] not in select_bkps_min.keys():
                select_bkps_min[bkps[i][0]] = 1
            else:
                select_bkps_min[bkps[i][0]] = select_bkps_min[bkps[i][0]] + 1

            if bkps[i][1] not in select_bkps_max.keys():
                select_bkps_max[bkps[i][1]] = 1
            else:
                select_bkps_max[bkps[i][1]] = select_bkps_max[bkps[i][1]] + 1

        final_bkps_min = max(select_bkps_min, key=select_bkps_min.get)
        final_bkps_max = max(select_bkps_max, key=select_bkps_max.get)
        final_bkps = str([final_bkps_min, final_bkps_max])

    else:
        final_bkps = '[]'
    back_reads_score_list = []
    back_reads_score_list.extend(reads_score_list)

    for k in range(len(reads_score_list)):
        if reads_score_list[k] == 'NA':
            reads_score_list[k] = -1
    support_reads_num = sum(i > 0 for i in reads_score_list)
    if len(reads_score_list) > 0:
        highest_score = max(reads_score_list)
        support_porpotion = "%.2f%%" % (support_reads_num * 100 / len(reads_score_list))
        if support_reads_num * 100 / len(reads_score_list) > 0.8:
            genotype_valid = '1/1'
        elif support_reads_num * 100 / len(reads_score_list) > 0.2:
            genotype_valid = '0/1'
        else:
            genotype_valid = '0/0'

    else:
        highest_score = 0
        support_porpotion = 0
        genotype_valid = 'None'

    # vcf record存入list
    chrom = vcf_sv.split('\t')[0]
    pos = vcf_sv.split('\t')[1]
    id = vcf_sv.split('\t')[2]
    ref = vcf_sv.split('\t')[3]
    alt = vcf_sv.split('\t')[4]
    qual = vcf_sv.split('\t')[5]
    filt = vcf_sv.split('\t')[6]
    info = vcf_sv.split('\t')[7].replace('\n', '')

    try:
        form = vcf_sv.split('\t')[8]
        simulate = vcf_sv.split('\t')[9]
        vcf_str = chrom + '\t' + pos + '\t' + id + '\t' + ref + '\t' + alt + '\t' + qual + '\t' + filt + '\t' + info + ';HIGHEST_SCORE=' + str(highest_score) + ';SUPPORT_PROPORPTION=' + str(
            support_porpotion) + ';GENOTYPE_VALID=' + genotype_valid + ';SCORE_LIST=[' + ','.join([str(x) for x in back_reads_score_list]) + ']' +';SV_BKPS='+ final_bkps + '\t' + form + '\t' + simulate + '\n'
        vcf_record.append(vcf_str)

    except IndexError:
        vcf_str = chrom + '\t' + pos + '\t' + id + '\t' + ref + '\t' + alt + '\t' + qual + '\t' + filt + '\t' + info + ';HIGHEST_SCORE=' + str(highest_score) + ';SUPPORT_PROPORPTION=' + str(
            support_porpotion) + ';GENOTYPE_VALID=' + genotype_valid + ';SCORE_LIST=[' + ','.join([str(x) for x in back_reads_score_list]) + ']' + ';SV_BKPS='+ final_bkps + '\n'
        vcf_record.append(vcf_str)
    return vcf_record

def read_vcf(vcf_file):
    '''

    :param vcf_file:
    :return:
    '''

    start_coord = []
    end_coord = []
    ins_coord = []
    chrom = []
    sv_type = []
    sv_length = []

    vcf_in = pysam.VariantFile(vcf_file)
    for rec in vcf_in:
        # ins_coord = 0
        chrom.append(rec.chrom)
        sv_type.append(rec.info['SVTYPE'])
        sv_length.append(rec.info['SVLEN'])
        if rec.info['SVTYPE'] in ['dDUP', 'idDUP']:
            dup_coord = []
            BKPS_str = str(rec.info['BKPS']).split(":")
            dup_coord.extend([int(rec.start + 1), int(rec.stop), int(BKPS_str[2]), int(BKPS_str[2]) + int(rec.info['SVLEN'])])
            start_coord.append(sorted(dup_coord)[0])
            end_coord.append(sorted(dup_coord)[-1])
            ins_coord.append(int(BKPS_str[2]))
        elif rec.info['SVTYPE'] in ['DEL+dDUP', 'DEL+idDUP']:
            BKPS_str = str(rec.info['BKPS']).split(":")
            ins_coord.append([int(BKPS_str[1].split("'")[0].split("-")[0]), int(BKPS_str[1].split("'")[0].split("-")[1]), int(BKPS_str[2].split("'")[0].split("-")[0]), int(BKPS_str[2].split("'")[0].split("-")[1]), int(BKPS_str[2].split("'")[0].split("-")[2])])
            start_coord.append(sorted([int(BKPS_str[2].split("'")[0].split("-")[0]),
                                       int(BKPS_str[2].split("'")[0].split("-")[1]),
                                       int(BKPS_str[2].split("'")[0].split("-")[2]),
                                       int(BKPS_str[2].split("'")[0].split("-")[2]) + int(BKPS_str[2].split("'")[0].split("-")[1])- int(BKPS_str[2].split("'")[0].split("-")[2]),
                                       int(BKPS_str[1].split("'")[0].split("-")[0]),
                                       int(BKPS_str[1].split("'")[0].split("-")[1])])[0])

            end_coord.append(sorted([int(BKPS_str[2].split("'")[0].split("-")[0]),
                                       int(BKPS_str[2].split("'")[0].split("-")[1]),
                                       int(BKPS_str[2].split("'")[0].split("-")[2]),
                                       int(BKPS_str[2].split("'")[0].split("-")[2]) + int(BKPS_str[2].split("'")[0].split("-")[1]) - int(BKPS_str[2].split("'")[0].split("-")[0]),
                                       int(BKPS_str[1].split("'")[0].split("-")[0]),
                                       int(BKPS_str[1].split("'")[0].split("-")[1])])[-1])

        elif rec.info['SVTYPE'] in ['DEL+INV']:
            BKPS_str = str(rec.info['BKPS']).split(":")
            start_coord.append(rec.start + 1)
            end_coord.append(rec.stop)
            ins_coord.append([int(BKPS_str[1].split("'")[0].split("-")[0]), int(BKPS_str[1].split("'")[0].split("-")[1]), int(BKPS_str[2].split("'")[0].split("-")[0]), int(BKPS_str[2].split("'")[0].split("-")[1])])
        else:
            start_coord.append(rec.start + 1)
            end_coord.append(rec.stop)
            ins_coord.append(0)
    return chrom, start_coord, end_coord, sv_type, sv_length, ins_coord


def get_seq_by_cigar(read, start_coord, end_coord):
    '''

    :param read:
    :param start_coord: ref的起始坐标
    :param end_coord: ref的结束坐标
    :return:
    '''
    ref_start = read.reference_start
    ref_end = read.reference_end
    # step 1：correct read start & read end
    cigar_tuple = read.cigartuples
    read_start = 0
    read_end = read.infer_read_length()
    for action, length in cigar_tuple:
        if action in [0, 7]:
            break
        if action in [4, 5]:
            read_start += length
    for action, length in cigar_tuple[::-1]:
        if action in [0, 7]:
            break
        if action in [4, 5]:
            read_end -= length
    if read.reference_start > start_coord and read.reference_end < end_coord:
        pointer_start, pointer_end = 0, read.infer_read_length()
        # alignment位于指定区间内，全部保留
        for action, length in cigar_tuple:
            if action in [0, 7]:
                break
            if action == 4:
                pointer_start += length

        for action, length in cigar_tuple[::-1]:
            if action in [0, 7]:
                break
            if action == 4:
                pointer_end -= length
    elif read.reference_start < start_coord + 1 and read.reference_end > end_coord + 1:
        # alignment start在区间左侧，alignment end在区间右侧，需要计算ref start， ref end
        # 相对位置
        read_rec = 0
        align_rec = ref_start
        start_count = 0
        pointer_start, start_miss_bp, pointer_end, end_miss_bp = 0, 0, read.infer_read_length(), 0
        for action, length in cigar_tuple:
            if action == 0:  # match
                read_rec += length
                align_rec += length
            if action == 1:  # insertion
                read_rec += length
            if action == 2:  # deletion
                align_rec += length
            if action == 4:  # soft clip or ref skip
                read_rec += length
            if action == 7:  # =, match
                read_rec += length
                align_rec += length
            if action == 8:  # mismatch(X)
                read_rec += length
                align_rec += length
            cigar_rec = (action, length)

            if align_rec > start_coord - 1:
                start_count += 1
                if start_count == 1:
                    start_dis = align_rec - start_coord
                    if cigar_rec[0] in [0, 7, 8]:
                        new_read_rec = read_rec - start_dis
                        new_start_dis = 0
                        pointer_start, start_miss_bp = new_read_rec, new_start_dis
                    else:
                        pointer_start, start_miss_bp = read_rec, start_dis

            if align_rec > end_coord - 1:
                end_dis = align_rec - end_coord
                if cigar_rec[0] in [0, 7, 8]:
                    new_end_rec = read_rec - end_dis
                    new_end_dis = 0
                    pointer_end = new_end_rec
                    end_miss_bp = new_end_dis
                else:
                    pointer_end = read_rec
                    end_miss_bp = end_dis
                break

    elif read.reference_start < start_coord + 1 and read.reference_end < end_coord + 1:
        # alignment start在区间左侧，需要计算ref start
        # 相对位置
        read_rec = 0
        align_rec = ref_start
        start_count = 0
        pointer_start, start_miss_bp, pointer_end, end_miss_bp = 0, 0, read.infer_read_length(), 0

        for action, length in cigar_tuple:
            if action == 0:  # match
                read_rec += length
                align_rec += length
            if action == 1:  # insertion
                read_rec += length
            if action == 2:  # deletion
                align_rec += length
            if action == 4:  # soft clip or ref skip
                read_rec += length
            if action == 7:  # =, match
                read_rec += length
                align_rec += length
            if action == 8:  # mismatch(X)
                read_rec += length
                align_rec += length
            cigar_rec = (action, length)

            if align_rec > start_coord - 1:
                start_count += 1
                if start_count == 1:
                    start_dis = align_rec - start_coord
                    if cigar_rec[0] in [0, 7, 8]:
                        new_read_rec = read_rec - start_dis
                        new_start_dis = 0
                        pointer_start, start_miss_bp = new_read_rec, new_start_dis
                    else:
                        pointer_start, start_miss_bp = read_rec, start_dis

        for action, length in cigar_tuple[::-1]:
            if action in [0, 7]:
                break
            if action == 4:
                pointer_end -= length
    elif read.reference_start > start_coord and read.reference_end > end_coord:
        # alignment end在区间右侧，需要计算ref end
        # 相对位置
        read_rec = 0
        align_rec = ref_start
        start_count = 0
        pointer_start, start_miss_bp, pointer_end, end_miss_bp = 0, 0, read.infer_read_length(), 0

        for action, length in cigar_tuple:
            if action == 0:  # match
                read_rec += length
                align_rec += length
            if action == 1:  # insertion
                read_rec += length
            if action == 2:  # deletion
                align_rec += length
            if action == 4:  # soft clip or ref skip
                read_rec += length
            if action == 7:  # =, match
                read_rec += length
                align_rec += length
            if action == 8:  # mismatch(X)
                read_rec += length
                align_rec += length
            cigar_rec = (action, length)

            if align_rec > end_coord - 1:
                end_dis = align_rec - end_coord
                if cigar_rec[0] in [0, 7, 8]:
                    new_end_rec = read_rec - end_dis
                    new_end_dis = 0
                    pointer_end, end_miss_bp = new_end_rec, new_end_dis
                else:
                    pointer_end, end_miss_bp = read_rec, end_dis
                break

        for action, length in cigar_tuple:
            if action in [0, 7]:
                break
            if action == 4:
                pointer_start += length
    else:
        pointer_start, pointer_end = 0, read.infer_read_length()
    align_seq = read.query_sequence[pointer_start:pointer_end]
    if read.is_reverse:
        # reverse
        infer_read_length = read.infer_read_length()
        read_start_reverse = infer_read_length - read_end
        read_end_reverse = infer_read_length - read_start
        sig_alignment = Alignment(ref_start, ref_end, read_start_reverse, read_end_reverse, read.is_reverse,
                                  Sequence(align_seq).get_reverse_complement(),'',read.qname)
    else:
        sig_alignment = Alignment(ref_start, ref_end, read_start, read_end, read.is_reverse, align_seq,'', read.qname)
    return sig_alignment


def get_bam_reads(bam_input, chrom, start_coord, end_coord, kmer_length):
    bam_file = pysam.AlignmentFile(bam_input, 'rb')
    all_alignments = bam_file.fetch(chrom, start_coord, end_coord)
    bam_reads = {}
    flank_length = set_flank_length()
    align_coord = {}  # {align_id:[[align1,align2][start1,end1,start2,end2,...]]}
    while True:
        try:
            sig_alignment = next(all_alignments)
            # 设置mapping quality > 0
            if sig_alignment.mapq > 0:
                align = get_seq_by_cigar(sig_alignment, start_coord, end_coord)
                align.read_seq = sig_alignment.query_sequence
                if align.align_id not in align_coord.keys():
                    align_coord[align.align_id] = []
                    align_coord[align.align_id].append([])    #存align对象
                    align_coord[align.align_id].append([])    #存align对象的reference坐标，排序用
                    # align_coord[align.align_id].append([])    #存align对象的read坐标，查漏unmapped seq用
                align_coord[align.align_id][0].append(align)
                align_coord[align.align_id][1].extend([align.ref_start, align.ref_end])

        except OSError:
            print('oserror')
            break
        except StopIteration:
            break
        except KeyboardInterrupt:
            print("Execution interrupted by KeyboardInterrupt!")
            break

    min_flank_length = []
    left_flank_length = []
    right_flank_length = []
    for key in align_coord.keys():
        if (min(align_coord[key][1]) < start_coord + flank_length - kmer_length) and (max(align_coord[key][1]) > end_coord - flank_length + kmer_length):
            # 16:length of kmer
            # 满足read start < sv start，但是可能不满足read start < section start (section start < sv start)
            # min_flank_length = [500,325,500,...]
            if min(align_coord[key][1]) - start_coord <= 0:
                left_flank_length.append(flank_length)
            elif min(align_coord[key][1]) - start_coord > 0:
                left_flank_length.append(start_coord + flank_length - min(align_coord[key][1]))

            if max(align_coord[key][1]) - end_coord >= 0:
                right_flank_length.append(flank_length)
            elif max(align_coord[key][1]) - end_coord < 0:
                right_flank_length.append(max(align_coord[key][1]) - end_coord + flank_length)
            if key not in bam_reads.keys():
                bam_reads[key] = []
            bam_reads[key].extend(align_coord[key][0])

    merge_reads = []
    for i in bam_reads.keys():
        # print(bam_reads[i])
        bam_reads_sorted = sorted(bam_reads[i], key=lambda aln: (aln.read_start, aln.read_end))
        read_str = ''
        for j in range(len(bam_reads_sorted)):
            if j != len(bam_reads_sorted) - 1 and (bam_reads_sorted[j + 1].read_start - bam_reads_sorted[j].read_end) > 50:
                if bam_reads_sorted[j].forward == False:
                    # 正向
                    read_str = read_str + bam_reads_sorted[j].seq
                    read_str = read_str + bam_reads_sorted[j].read_seq[bam_reads_sorted[j].read_end:bam_reads_sorted[j + 1].read_start]
                else:
                    # 反向
                    read_str = read_str + bam_reads_sorted[j].seq
                    read_str = read_str + Sequence(bam_reads_sorted[j].read_seq[len(bam_reads_sorted[j].read_seq) - bam_reads_sorted[j + 1].read_start: len(bam_reads_sorted[j].read_seq) - bam_reads_sorted[j].read_end]).get_reverse_complement()
            else:
                read_str = read_str + bam_reads_sorted[j].seq
        merge_reads.append([read_str, i])

    return merge_reads, left_flank_length, right_flank_length

def vcf_run(threads, sv_input, out_path,ref_path,bam_input,process_images, kmer_length):
    wall_time_start = time.time()
    cpu_time_start = time.process_time()

    #read vcf()
    chrom, start_coord, end_coord, sv_type, sv_length, ins_coord = read_vcf(sv_input)
    vcf_head = []
    vcf_title = []
    vcf_sv = []

    with open(sv_input, 'r') as ori_file:
        ori_file_list = ori_file.readlines()

    for i in range(len(ori_file_list)):
        if ori_file_list[i][0:2] == '##':
            vcf_head.append(ori_file_list[i])
        elif ori_file_list[i][0] == '#':
            vcf_title.append(ori_file_list[i])
        else:
            vcf_sv.append(ori_file_list[i])

    vcf_head.append('##INFO=<ID=HIGHEST_SCORE,Number=1,Type=String,Description="Hightet Validation Score">\n')
    vcf_head.append('##INFO=<ID=SUPPORT_PROPORPTION,Number=1,Type=String,Description="Proportion of Supporting SV reads">\n')
    vcf_head.append('##INFO=<ID=GENOTYPE_VALID,Number=1,Type=String,Description="Genotype of the Proposed SV as Assessed">\n')
    vcf_head.append('##INFO=<ID=SCORE_LIST,Number=1,Type=String,Description="All Validation Score">\n')
    vcf_head.append('##INFO=<ID=SV_BKPS,Number=1,Type=String,Description="Breakpoints of this SV section">\n')

    # mkdir outpath
    set_dir(out_path)
    out_fig = os.path.join(out_path, 'out_fig')
    out_result = os.path.join(out_path, 'out_result')
    set_dir(out_fig)
    set_dir(out_result)

    with open(os.path.join(out_result, 'Validation.vcf'),
              'a') as val_result_file:
        for item in vcf_head:
            val_result_file.write(item)
        for item in vcf_title:
            val_result_file.write(item)

    vcf_record = []
    process_pool = multiprocessing.Pool(threads)
    for i in range(len(chrom)):
        with open(os.path.join(out_result, 'output.log'), 'a') as log_file:
            if end_coord[i] - start_coord[i] <= 1000000:
                log_file.write('Processing record:' + str(chrom[i]) +'\t'+ str(start_coord[i]) +'\t'+ str(end_coord[i]) +'\t'+ str(sv_type[i]) + '\n')
                vcf_record = process_pool.apply_async(sig_vcf_process, (ref_path, bam_input, chrom[i], start_coord[i], end_coord[i], sv_type[i], sv_length[i], ins_coord[i], vcf_sv[i], out_fig, process_images, vcf_record, kmer_length)).get()
    process_pool.close()
    process_pool.join()

    with open(os.path.join(out_result, 'Validation.vcf'), 'a') as vcf_file:
        for r in vcf_record:
            vcf_file.write(r)

    os.system('bcftools sort ' + os.path.join(out_result, 'Validation.vcf ') + '-o ' + os.path.join(out_result, os.path.basename(sv_input).split('.')[0] + 'Validation.vcf'))
    wall_time_end = time.time()
    cpu_time_end = time.process_time()

    print("wall time：{}s".format(wall_time_end - wall_time_start))
    print("cpu time：{}s".format(cpu_time_end - cpu_time_start))

    with open(os.path.join(out_result, 'output.log'), 'a') as log_file:
        log_file.write("wall time：{}s\n".format(wall_time_end - wall_time_start))
        log_file.write("cpu time：{}s".format(cpu_time_end - cpu_time_start))