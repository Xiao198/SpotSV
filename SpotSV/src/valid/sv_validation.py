import signal
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from valid.segment import Segment
from valid.sequence import Sequence
from valid.hash_aligner import HashAligner

def run_aligner(ref, seq, kmer_length):
    k = kmer_length  # length of kmer
    window_size = kmer_length  # length extended to both sides
    mismatch_num = 0
    ref = Sequence(ref)  # generate ref sequence
    read = Sequence(seq)  # generate read sequence
    aligner_ref = HashAligner(k, window_size, mismatch_num)
    aligner_ref.run(read, ref)
    segments_ref = aligner_ref.getMergeSegments()
    return segments_ref

def get_segment_list(ref, read, alt, kmer_length):
    segments_ref_ref = run_aligner(ref, ref, kmer_length)
    segments_ref_read = run_aligner(ref, read, kmer_length)
    segments_alt_alt = run_aligner(alt, alt, kmer_length)
    segments_alt_read = run_aligner(alt, read, kmer_length)
    return segments_ref_ref, segments_ref_read, segments_alt_alt, segments_alt_read

def sv_distance_calculate(filtered_segments_sorted):
    sig_ref_read_distance = []
    for item in filtered_segments_sorted:
        read_start_coord = item.get_xstart_coord()
        read_mid_coord = item.get_xmid_coord()
        read_end_coord = item.get_xend_coord()

        ref_start_coord = item.get_ystart_coord()
        ref_mid_coord = item.get_ymid_coord()
        ref_end_coord = item.get_yend_coord()
        sig_ref_read_distance.append((abs(ref_start_coord - read_start_coord) + abs(ref_end_coord - read_end_coord)) / 2)
    sum_read_ref = 0
    for item in sig_ref_read_distance:
        sum_read_ref += item
    try:
        score_read_ref = abs(sum_read_ref) / len(sig_ref_read_distance)
    except ZeroDivisionError:
        score_read_ref = 0
    return score_read_ref


def calDiffBetTow(i, j):
    if abs(float(i.y_start - j.y_start)) == 0:
        return 5
    return abs(float(i.x_start - j.x_start)) / abs(float(i.y_start - j.y_start))

def linearOrNot(i, j):
    # merge condition1: different strand, pass directly
    if i.seg_forward != j.seg_forward:
        return False
    # merge condition2: colinear, not then false
    DIFF = calDiffBetTow(i, j)
    if DIFF > 1.1 or DIFF < 0.9:
        return False
    # merge condition3: colinear, yes, and dis
    DIS_X = abs(i.x_end - j.x_start)
    DIS_Y = abs(i.y_end - j.y_start)
    # maxDis = (i.length() + j.length()) * 1.5
    maxDis = 50
    if DIS_X > maxDis and DIS_Y > maxDis:
        return False
    # merge condition4: merged line's k is not -1 or 1
    tmp = float(j.x_end - i.x_start)
    if tmp == 0:
        tmp = 0.0001
    k = float(j.y_end - i.y_start) / tmp
    if abs((abs(k) - 1)) > 0.2:
        return False
    return True

def del_repeat_segment(diagonal_threshold, segments_ref_read, segments_ref_ref, sv_length, sv_type, left_flank_length, right_flank_length, read_length):
    # emputy list to store deleted repeat segments
    filtered_segments = []
    window_size = 20  # length extended to both sides
    remain_seg = []  #要保留的主对角线seg（SV左侧，SV右侧）

    for i in range(len(segments_ref_read)):
        find_min_length_dict = dict()
        find_left_length = []
        find_right_length = []
        leftset_dict = dict()
        rightset_dict = dict()

        # get the start coordinate of ref in the ref vs read graph
        readseg_ref_start_coord = segments_ref_read[i].get_ystart_coord()
        readseg_ref_end_coord = segments_ref_read[i].get_yend_coord()
        # ref section of segment in the ref vs read graph
        readseg_ref_range = set(range(readseg_ref_start_coord, readseg_ref_end_coord + 1))

        # get the start coordinate of read in the ref vs read graph
        readseg_read_start_coord = segments_ref_read[i].get_xstart_coord()
        readseg_read_end_coord = segments_ref_read[i].get_xend_coord()

        # Main diagonal segment preservation - add diagonal threshold
        # SV左侧， 全部SV主对角线（ins, del, inv, tdup, ddup）
        if readseg_read_start_coord <= left_flank_length and abs(readseg_ref_start_coord - readseg_read_start_coord) <= diagonal_threshold:
            # filtered_segments.append(segments_ref_read[i])
            remain_seg.append(segments_ref_read[i])
            continue

        # SV右侧，主对角线截距
        elif readseg_read_start_coord >= read_length - right_flank_length:
            if sv_type in [['INS'],['DUP'],['tDUP'],['itDUP'],['dDUP'],['idDUP']] and (readseg_read_start_coord - readseg_ref_start_coord) >= sv_length - diagonal_threshold and readseg_read_start_coord - readseg_ref_start_coord <= sv_length + diagonal_threshold:
                # filtered_segments.append(segments_ref_read[i])
                remain_seg.append(segments_ref_read[i])
                continue
            elif sv_type in [['DEL'],] and (readseg_ref_start_coord - readseg_read_start_coord) >= sv_length - diagonal_threshold and (readseg_ref_start_coord - readseg_read_start_coord) <= sv_length + diagonal_threshold:
                # filtered_segments.append(segments_ref_read[i])
                remain_seg.append(segments_ref_read[i])
                continue
            elif sv_type in [['INV'],] and abs(readseg_ref_start_coord - readseg_read_start_coord) <= diagonal_threshold:
                # filtered_segments.append(segments_ref_read[i])
                remain_seg.append(segments_ref_read[i])
                continue

        # else:
        not_repeat_count = 0
        for j in range(len(segments_ref_ref)):
            refseg_ref_start_coord = segments_ref_ref[j].get_ystart_coord()
            refseg_ref_end_coord = segments_ref_ref[j].get_yend_coord()
            refseg_ref_range = set(range(refseg_ref_start_coord, refseg_ref_end_coord + 1))

            # get the start coordinate of segments in ref vs ref graph
            refseg_read_start_coord = segments_ref_ref[j].get_xstart_coord()
            refseg_read_end_coord = segments_ref_ref[j].get_xend_coord()
            # ignore diagonal
            # Main diagonal segment preservation - add diagonal threshold
            if abs(refseg_ref_start_coord - refseg_read_start_coord) <= diagonal_threshold and abs(
                    refseg_ref_end_coord - refseg_read_end_coord) <= diagonal_threshold:
                not_repeat_count = not_repeat_count + 1
                if not_repeat_count == len(segments_ref_ref):
                    filtered_segments.append(segments_ref_read[i])
                if j != len(segments_ref_ref) - 1:
                    continue
            else:
                # 做差集
                difference_range = readseg_ref_range.difference(refseg_ref_range)
                # 差集不为空
                if difference_range != set():
                    # 差集等于ref range，没有overlap
                    if difference_range == readseg_ref_range:
                        not_repeat_count = not_repeat_count + 1
                        if not_repeat_count == len(segments_ref_ref):
                            filtered_segments.append(segments_ref_read[i])
                        if j != len(segments_ref_ref) - 1:
                            continue
                    # 差集有overlap
                    else:
                        min_value = min(difference_range)
                        max_value = max(difference_range)
                        leftset = set()
                        rightset = set()

                        for k in range(len(difference_range)):
                            if min_value + k in difference_range:
                                leftset.add(min_value + k)
                            else:
                                break

                        for l in range(len(difference_range)):
                            if max_value - l in difference_range:
                                rightset.add(max_value - l)
                            else:
                                break
                        leftset_dict[j] = leftset
                        rightset_dict[j] = rightset

                        if leftset == rightset:
                            update_length = max(leftset) - min(leftset) + 1
                            find_min_length_dict[j] = update_length
                        else:
                            left_update_length = max(leftset) - min(leftset) + 1
                            right_update_length = max(rightset) - min(rightset) + 1
                            find_left_length.append(left_update_length)
                            find_right_length.append(right_update_length)
                else:
                    # 如果有一条read seg被ref seg包围，直接下一条read
                    break

            if j == len(segments_ref_ref) - 1:
                if find_min_length_dict != dict() and find_left_length == [] and find_right_length == []:
                    update_length_index = min(find_min_length_dict, key=lambda x: find_min_length_dict[x])
                    update_length = find_min_length_dict[update_length_index]

                    if min(leftset_dict[update_length_index]) == readseg_ref_start_coord and update_length > window_size:
                        new_seg = Segment(readseg_read_start_coord, readseg_ref_start_coord, update_length,
                                          segments_ref_read[i].get_forward(), segments_ref_read[i].get_id())
                        filtered_segments.append(new_seg)
                    elif max(leftset_dict[update_length_index]) == readseg_ref_end_coord and update_length > window_size:
                        if segments_ref_read[i].get_forward():
                            new_seg = Segment(readseg_read_end_coord - update_length,
                                              readseg_ref_end_coord - update_length, update_length,
                                              segments_ref_read[i].get_forward(), segments_ref_read[i].get_id())
                            filtered_segments.append(new_seg)
                        else:
                            new_seg = Segment(readseg_read_end_coord + update_length,
                                              readseg_ref_end_coord - update_length, update_length,
                                              segments_ref_read[i].get_forward(), segments_ref_read[i].get_id())
                            filtered_segments.append(new_seg)

                    else:
                        pass
                elif find_min_length_dict == dict() and find_left_length != [] and find_right_length != []:

                    # 小于window size的left_length或right_length置零
                    for m in range(len(find_left_length)):
                        if find_left_length[m] <= window_size:
                            find_left_length[m] = 0

                    for n in range(len(find_right_length)):
                        if find_right_length[n] <= window_size:
                            find_right_length[n] = 0
                    sum_length = find_left_length[0] + find_right_length[0]

                    for p in range(len(find_left_length)):
                        if find_left_length[p] + find_right_length[p] < sum_length:
                            left_update_length = find_left_length[p]
                            right_update_length = find_right_length[p]
                            sum_length = left_update_length + right_update_length

                    if left_update_length > window_size:
                        left_new_seg = Segment(readseg_read_start_coord, readseg_ref_start_coord,
                                               left_update_length, segments_ref_read[i].get_forward(),
                                               str(segments_ref_read[i].get_id()) + '-1')
                        filtered_segments.append(left_new_seg)
                    if right_update_length > window_size:
                        if segments_ref_read[i].get_forward():
                            right_new_seg = Segment(readseg_read_end_coord - right_update_length,
                                                    readseg_ref_end_coord - right_update_length, right_update_length,
                                                    segments_ref_read[i].get_forward(),
                                                    str(segments_ref_read[i].get_id()) + '-2')
                        else:
                            right_new_seg = Segment(readseg_read_end_coord + right_update_length,
                                                    readseg_ref_end_coord - right_update_length, right_update_length,
                                                    segments_ref_read[i].get_forward(),
                                                    str(segments_ref_read[i].get_id()) + '-2')
                        filtered_segments.append(right_new_seg)

                elif find_min_length_dict != dict() and find_left_length != [] and find_right_length != []:
                    update_length_index = min(find_min_length_dict, key=lambda x: find_min_length_dict[x])
                    update_length = find_min_length_dict[update_length_index]

                    # 小于window size的left_length或right_length置零
                    for m in range(len(find_left_length)):
                        if find_left_length[m] <= window_size:
                            find_left_length[m] = 0

                    for n in range(len(find_right_length)):
                        if find_right_length[n] <= window_size:
                            find_right_length[n] = 0
                    sum_length = find_left_length[0] + find_right_length[0]

                    for p in range(len(find_left_length)):
                        if find_left_length[p] + find_right_length[p] < sum_length:
                            left_update_length = find_left_length[p]
                            right_update_length = find_right_length[p]
                            sum_length = left_update_length + right_update_length

                    if update_length <= sum_length:
                        if min(leftset_dict[update_length_index]) == readseg_ref_start_coord and update_length > window_size:
                            new_seg = Segment(readseg_read_start_coord, readseg_ref_start_coord, update_length,
                                              segments_ref_read[i].get_forward(), segments_ref_read[i].get_id())
                            filtered_segments.append(new_seg)
                        elif max(leftset_dict[update_length_index]) == readseg_ref_end_coord and update_length > window_size:
                            if segments_ref_read[i].get_forward():
                                new_seg = Segment(readseg_read_end_coord - update_length,
                                                  readseg_ref_end_coord - update_length, update_length,
                                                  segments_ref_read[i].get_forward(), segments_ref_read[i].get_id())
                                filtered_segments.append(new_seg)
                            else:
                                new_seg = Segment(readseg_read_end_coord + update_length,
                                                  readseg_ref_end_coord - update_length, update_length,
                                                  segments_ref_read[i].get_forward(), segments_ref_read[i].get_id())
                                filtered_segments.append(new_seg)

                        else:
                            pass
                    else:
                        if left_update_length > window_size:
                            left_new_seg = Segment(readseg_read_start_coord, readseg_ref_start_coord,
                                                   left_update_length, segments_ref_read[i].get_forward(),
                                                   str(segments_ref_read[i].get_id()) + '- 1')
                            filtered_segments.append(left_new_seg)
                        if right_update_length > window_size:
                            if segments_ref_read[i].get_forward():
                                right_new_seg = Segment(readseg_read_end_coord - right_update_length,
                                                        readseg_ref_end_coord - right_update_length,
                                                        right_update_length,
                                                        segments_ref_read[i].get_forward(),
                                                        str(segments_ref_read[i].get_id()) + '- 2')
                            else:
                                right_new_seg = Segment(readseg_read_end_coord + right_update_length,
                                                        readseg_ref_end_coord - right_update_length,
                                                        right_update_length,
                                                        segments_ref_read[i].get_forward(),
                                                        str(segments_ref_read[i].get_id()) + '- 2')
                            filtered_segments.append(right_new_seg)

    # 主对角线过滤
    curSegNum = 1

    while curSegNum < len(remain_seg):
        flag = 0
        curSeg = remain_seg[curSegNum]

        for i in range(curSegNum):
            candiSeg = remain_seg[i]

            if linearOrNot(candiSeg, curSeg):
                if curSeg.seg_forward == True:
                    candiSeg.set_length(max(abs(curSeg.x_end - candiSeg.x_start), abs(curSeg.y_end - candiSeg.y_start)))
                    candiSeg.set_x_end(candiSeg.x_start + candiSeg.seg_length)
                elif curSeg.seg_forward == False:
                    candiSeg.set_length(abs(candiSeg.seg_length) + max(abs(curSeg.x_end - candiSeg.x_end), abs(curSeg.y_end - candiSeg.y_end)))
                    candiSeg.set_x_end(candiSeg.x_start - candiSeg.seg_length)
                candiSeg.set_y_end(candiSeg.y_start + candiSeg.seg_length)

                # remove curSeg
                remain_seg.remove(curSeg)
                flag = 1
                break
        if flag == 0:
            curSegNum += 1
    after_filte_segs = []
    for seg in remain_seg:
        if (seg.y_end - seg.y_start) >= 20:
            # print(seg.yEnd() - seg.yStart())
            after_filte_segs.append(seg)
        remain_seg = after_filte_segs
    remain_seg_del = []
    for i in range(len(remain_seg) - 1):
        for j in range(i + 1, len(remain_seg)):
            if remain_seg[i].x_start <= remain_seg[j].x_start and remain_seg[i].x_end >= remain_seg[j].x_end:
                if remain_seg[j] not in remain_seg_del:
                    remain_seg_del.append(remain_seg[j])
    if remain_seg_del != []:
        for item in remain_seg_del:
            remain_seg.remove(item)
    filtered_segments.extend(remain_seg)

    # 根据read segment的read起始坐标进行排序
    filtered_segments_sorted = sorted(filtered_segments, key=lambda aln: (aln.x_start, aln.x_end))

    return filtered_segments_sorted


def judge_align_reverse(read, segments):
    '''
    判断是否反向比对
    :return: 
    '''
    reverse_flag = False
    ref_segments_sorted = sorted(segments, key=lambda aln: min(aln.x_start, aln.x_end))
    read_forward_list = []
    for item in ref_segments_sorted:
        read_forward_list.append(item.get_forward())

    if False not in read_forward_list:
        # 全是True，不反向
        reverse_flag = False
        pass
    elif True not in read_forward_list:
        # 全是False，反向
        reverse_flag = True
        pass
    elif read_forward_list[0] == False:
        # 第一条seg反向
        reverse_flag = True
        pass
    else:
        pass
    return reverse_flag

def get_reverse_segment(read, segments):
    for item in segments:
        # print(item.get_xstart_coord())
        item.set_x_start(len(read) - 1 - item.get_xstart_coord())
        # print(item.get_xstart_coord())
        item.set_x_end(len(read) - 1 - item.get_xend_coord())
        item.set_forward(bool(1 - item.get_forward()))

def search_bkps(segments, diagonal_threshold, ref_start_coord):
    # 根据read segment的read起始坐标进行排序
    segments_sorted = sorted(segments, key=lambda aln: (min(aln.x_start, aln.x_end)))
    bkps = []

    for i in range(len(segments_sorted) - 1):
        # print((i+1) * diagonal_threshold / len(segments_sorted),'???')
        if segments_sorted[i + 1].get_forward() != segments_sorted[i].get_forward():
            if segments_sorted[i].get_forward() == True:
                # 当前segment为True，下一条False
                bkps.append(segments_sorted[i].get_yend_coord() + ref_start_coord)
                bkps.append(segments_sorted[i + 1].get_yend_coord() + ref_start_coord)
            else:
                # 当前segment为False，下一条True
                bkps.append(segments_sorted[i].get_ystart_coord() + ref_start_coord)
                bkps.append(segments_sorted[i + 1].get_ystart_coord() + ref_start_coord)
        elif abs(abs(segments_sorted[i+1].get_ystart_coord() - segments_sorted[i].get_yend_coord()) - abs(segments_sorted[i+1].get_xstart_coord() - segments_sorted[i].get_xend_coord())) > diagonal_threshold:
            bkps.append(segments_sorted[i].get_yend_coord()+ ref_start_coord)
            bkps.append(segments_sorted[i+1].get_ystart_coord()+ ref_start_coord)
    bkps = sorted(bkps)
    return bkps

def run_validation(ref, read, alt, sv_type, ref_start_coord, sv_length, left_flank_length, right_flank_length, process_images, kmer_length):
    diagonal_threshold = int(0.015 * max(len(ref), len(read)))
    # step1: 比对
    segments_ref_ref, segments_ref_read, segments_alt_alt, segments_alt_read = get_segment_list(ref, read, alt, kmer_length)
    # 在这里判断是否反向比对
    reverse_flag1 = judge_align_reverse(read, segments_ref_read)
    reverse_flag2 = judge_align_reverse(read, segments_alt_read)
    if reverse_flag1:
        get_reverse_segment(read, segments_ref_read)
    if reverse_flag2:
        get_reverse_segment(read, segments_alt_read)
    # step2: 去重
    filtered_ref_read_sorted = del_repeat_segment(diagonal_threshold, segments_ref_read, segments_ref_ref, sv_length, sv_type, left_flank_length, right_flank_length, len(read))
    filtered_alt_read_sorted = del_repeat_segment(diagonal_threshold, segments_alt_read, segments_alt_alt, 0, sv_type, left_flank_length, right_flank_length, len(read))
    bkps = search_bkps(filtered_ref_read_sorted, diagonal_threshold, ref_start_coord)
    final_bkps = []
    if bkps != []:
        final_bkps.extend([min(bkps), max(bkps)])
    else:
        pass
    plot_tuple = (ref, read, alt, segments_ref_ref, segments_ref_read, filtered_ref_read_sorted, segments_alt_alt, segments_alt_read, filtered_alt_read_sorted, process_images)
    # setp4: sv score
    score_ref_read = sv_distance_calculate(filtered_ref_read_sorted)
    score_alt_read = sv_distance_calculate(filtered_alt_read_sorted)
    return score_ref_read, score_alt_read, final_bkps, plot_tuple