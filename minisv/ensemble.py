import os
import numpy as np
from .regex import re_info


def insilico_truth(msv_union, file_handler=None):

    def _emit(res_list, file_ids):
        sorted_lines = sorted(res_list, key=lambda x: (int(x[1]), int(x[4])))
        med_line = '\t'.join(map(str, sorted_lines[len(sorted_lines) >> 1]))
        print(med_line + '\t' + ','.join(file_ids), file=file_handler)

    group_ids = []
    file_ids = []
    res_list = []
    with open(msv_union) as inf:
        for line in inf:
            elements = line.strip().split('\t')
            # start
            elements[1] = int(elements[1])
            # end
            elements[4] = int(elements[4])
            group_id = ''
            file_id = ''
            for (key, val) in re_info.findall(elements[5]):
                 if key == 'group_id':
                     group_id = val
                 if key == 'file_id':
                     file_id  = val

            assert group_id != ''
            assert file_id != ''

            if len(group_ids) >= 1:
                if group_id == group_ids[-1]:
                    res_list.append(elements)
                    file_ids.append(file_id)
                else:
                    # output the completed group
                    _emit(res_list, file_ids)
                    group_ids.append(group_id)
                    file_ids = [file_id]
                    res_list = [elements]
            else:
                file_ids.append(file_id)
                group_ids.append(group_id)
                res_list.append(elements)

    # flush the final group: the loop only emits on a group_id transition,
    # so the last accumulated group is otherwise never output.
    if res_list:
        _emit(res_list, file_ids)


def double_strand_break(collapsed_msv_union):
    """ Identify four cases of double strand breaks
    """
    breakpts = []
    with open(collapsed_msv_union) as inf:
        for line in inf:
            elements = line.strip().split('\t')
            svtype = [ val for (key, val) in re_info.findall(elements[-1]) if key == 'SVTYPE' ][0]

            #if svtype in ['INS', 'DUP']:
            #    breakpts.append([elements[0], int(elements[1]), elements[2][0], elements[-1]])
            #    breakpts.append([elements[3], int(elements[4]), elements[2][1], elements[-1]])
            #    continue
            #if svtype in ['DEL']:
            #    breakpts.append([elements[0], int(elements[1]), elements[2][0], elements[-1], ">"])
            #    breakpts.append([elements[3], int(elements[4]), elements[2][1], elements[-1], "<"])
            #    continue

            # not sure, may inspect reads in asm
            breakpts.append([elements[0], int(elements[1]), elements[2][0], elements[-1]])
            breakpts.append([elements[3], int(elements[4]), elements[2][1], elements[-1]])


    breakpts.sort(key = lambda x: (x[0], int(x[1])))
    dtype_list = []

    for index, pt0 in enumerate(breakpts):
        for pt1 in breakpts[(index+1):]:
            dsbtype = ""
 
            # skip breakpoints on diff contigs
            if pt0[0] != pt1[0]:
                continue
 
            # skip very far away breakpoints
            if pt0[0] == pt1[0] and abs(pt1[1] - pt0[1]) >= 25e3:
                continue

            svtype0 = [ val for (key, val) in re_info.findall(pt0[-1]) if key == 'SVTYPE' ][0]
            svtype1 = [ val for (key, val) in re_info.findall(pt1[-1]) if key == 'SVTYPE' ][0]

            asmrid0 = [ val for (key, val) in re_info.findall(pt0[-1]) if key == 'asm_read_ids' ][0].split(',')
            asmrid1 = [ val for (key, val) in re_info.findall(pt1[-1]) if key == 'asm_read_ids' ][0].split(',')

            # ['chr10', 51935447, '<', '1,2,3', 'SVTYPE=BND;SVLEN=0;count=13;group_id=41;file_id=2', 'left(+)'] ['chr10', 51935447, '<', '0', 'SVTYPE=BND;SVLEN=0;count=13;group_id=935;file_id=0', 'right(-)'] BND BND
            # skip the same breakpoint coordinate from ensembling
            if pt0[0] == pt1[0] and pt0[1] == pt1[1]:
                continue

            # simple insertion skip
            if abs(pt1[1] - pt0[1]) <= 49:
                continue

            # NOTE: insertion first filter, then overlap second filter, last parallel, and case 1

            # from same read or different reads
            # insertion case 3
            if (svtype1 == "BND" or svtype0 == "BND") and (len(list(set(asmrid0) & set(asmrid1))) >= 1) and (len(list(set(asmrid0) & set(asmrid1))) == len(asmrid0)): # not work due to we sort the orientation and pt0[2] == "<" and pt1[2] == ">":
                dsbtype = "Insertions_case3"
                dtype_list.append(dsbtype)
                print('\t'.join(map(str, pt0)), '\t'.join(map(str, pt1)), svtype0, svtype1, dsbtype)
                continue

            # overlap reads
            # case 2
            if svtype0 == svtype1 and svtype0 == "DUP":
                dsbtype = "case2"
                dtype_list.append(dsbtype)
                print('\t'.join(map(str, pt0)), '\t'.join(map(str, pt1)), svtype0, svtype1, dsbtype)
                continue

            # inversion
            # case 4 from the same read
            if ((svtype0 == svtype1 and svtype0 == "INV") or ((svtype0 == "BND" or svtype1 == "BND") and ((pt0[2] == "<" and pt1[2] == "<") or (pt0[2] == ">" and pt1[2] == ">")))) and (len(list(set(asmrid0) & set(asmrid1))) >= 1):
                dsbtype = "Parallel_Breakpoints_case4"
                dtype_list.append(dsbtype)
                print('\t'.join(map(str, pt0)), '\t'.join(map(str, pt1)), svtype0, svtype1, dsbtype)
                continue

            # case 1 from the same read
            if svtype0 == svtype1 and svtype0 == "DEL" and (len(list(set(asmrid0) & set(asmrid1))) >= 1):
                dsbtype = "case1"
                dtype_list.append(dsbtype)
                print('\t'.join(map(str, pt0)), '\t'.join(map(str, pt1)), svtype0, svtype1, dsbtype)
                continue

            # case 1 from the diff read
            if pt0[-1] == '>' and pt1[-1] == '<' and (svtype0 == "BND" or svtype1 == "BND"):
                dsbtype = "case1"
                dtype_list.append(dsbtype)
                print('\t'.join(map(str, pt0)), '\t'.join(map(str, pt1)), svtype0, svtype1, dsbtype)
                continue


    val, cnt = np.unique(dtype_list, return_counts=True)
    print('#stat')
    print('\t'.join(list(val)))
    print('\t'.join(list(map(str, cnt))))

