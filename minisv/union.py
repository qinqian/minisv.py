from .eval import gc_eval_merge_sv, gc_parse_sv, gc_cmp_sv, iit_overlap, gc_cmp_same_sv1, gc_read_bed
import numpy as np
from operator import attrgetter
from .filtercaller import parse_readids, parse_svid, parse_msvasm, parse_msvasm_withtype
from .type import simple_type, get_type


def union_sv(msvs, read_min_len, opt, file_handler=None, min_file_count=None):
    if opt.bed is not None and not isinstance(opt.bed, dict):
        opt.bed = gc_read_bed(opt.bed)
    sv = []

    # cat all the svs
    for i, m in enumerate(msvs):
        v =  gc_parse_sv(m, read_min_len, opt.read_min_count, False, False)
        for j, vj in enumerate(v):
            t = vj
            t.file_id = i
            t.in_bed = True
            if opt.bed is not None:
                if t.ctg not in opt.bed or t.ctg2 not in opt.bed:
                    t.in_bed = False
                elif len(iit_overlap(opt.bed[t.ctg], t.pos, t.pos+1)) == 0:
                    t.in_bed = False
                elif len(iit_overlap(opt.bed[t.ctg2], t.pos2, t.pos2+1)) == 0:
                    t.in_bed = False
            sv.append(t)
   
    sv.sort(key=attrgetter('ctg', 'pos'))

    group = []
    n_ambi1 = 0
    n_ambi2 = 0

    # this follows the merge operation
    for i, v in enumerate(sv):
        merge_to = []
        is_ambi1 = False
        for j in range(len(group)-1, -1, -1):
            g = group[j]
            if g[0].ctg != v.ctg:
                 break
            # NOTE: maybe g[-1], to compare the closet one
            if v.pos - g[0].pos > opt.win_size:
                 break
            n_same = 0
            for k, gv in enumerate(g):
                if gc_cmp_same_sv1(opt.win_size, opt.min_len_ratio, g[k], sv[i]):
                    n_same += 1
            if n_same > 0:
                merge_to.append(j)
                # NOTE: save sv number not equal to group size
                if n_same != len(g):
                    is_ambi1 = True

        if is_ambi1:
            n_ambi1 += 1

        if len(merge_to) == 0:
            group.append([v])
        else:
            # NOTE: merge to the closet group
            group[merge_to[0]].append(v)
            # NOTE: two possible merge group
            if len(merge_to) > 1:
                n_ambi2 += 1

    # summarize the results
    cnt = []
    for i in range(0, 1<<len(msvs)):
        cnt.append(0)

    for j in range(len(group)):
        g = group[j]
        in_bed = False
        has_bnd = False
        max_len = 0
        max_cnt = 0
        x = 0
        for k in range(len(g)):
            if g[k].in_bed:
                in_bed = True

            max_cnt = max(max_cnt, g[k].count)
            max_len = max(max_len, abs(g[k].SVLEN))
            if g[k].SVTYPE == "BND": 
                has_bnd = True
            x |= 1<<g[k].file_id
            g[k].group_id = j
        # NOTE: different sv tools has different metrics of counts
        #       computed the max of each group to the threshold
        if max_cnt < opt.group_min_count:
            continue
        # non-tranlocation length condition
        if has_bnd == False and max_len < opt.min_len:
            continue
        if not in_bed:
            continue
        if min_file_count is not None and bin(x).count("1") < min_file_count:
            continue
        cnt[x] += 1
        if opt.print_sv:
            for k in range(len(g)):
                print("\t".join(map(str, gc_sv2array(g[k]))), file=file_handler)

    if not opt.print_sv:
        for x in range(1, len(cnt)):
            label = []
            for i in range(len(msvs)):
                label.append(x>>i&1)
            print("".join(map(str, label)), cnt[x], sep="\t", file=file_handler)


def gc_sv2array(s):
    attr = [f"SVTYPE={s.SVTYPE}", f"SVLEN={s.SVLEN}", f"count={s.count}"]
    if s.group_id is not None:
        attr.append(f"group_id={s.group_id}")
    if s.file_id is not None:
        attr.append(f"file_id={s.file_id}")
    if s.readids is not None:
        attr.append(f"orig_read_ids={','.join(s.readids)}")
    if s.asmreadids is not None:
        attr.append(f"asm_read_ids={','.join(s.asmreadids)}")
    if s.tr_info is not None:
        attr.append(f"tr_info={'#'.join(s.tr_info)}")
    if s.segdup_info is not None:
        attr.append(f"segdup_info={'#'.join(s.segdup_info)}")
    return [s.ctg, s.pos, s.ori, s.ctg2, s.pos2, ";".join(attr)]


def advunion_sv(msvs, readids, msvasm, read_min_len, opt):
    rids = parse_msvasm(msvasm)
    rids_withtype = parse_msvasm_withtype(msvasm)

    if opt.bed is not None:
        opt.bed = gc_read_bed(opt.bed)

    sv = []
    # cat all the svs
    for i, (m, r) in enumerate(zip(msvs, readids)):
        parsed_id_dict = parse_readids(r)
        v =  gc_parse_sv(m, read_min_len, opt.read_min_count, False, False)

        for j, vj in enumerate(v):
            t = vj
            t.file_id = i
            t.in_bed = True
            if opt.bed is not None:
                if t.ctg not in opt.bed or t.ctg2 not in opt.bed:
                    t.in_bed = False
                elif len(iit_overlap(opt.bed[t.ctg], t.pos, t.pos+1)) == 0:
                    t.in_bed = False
                elif len(iit_overlap(opt.bed[t.ctg2], t.pos2, t.pos2+1)) == 0:
                    t.in_bed = False

            query_svid = str(parse_svid(t.svid))
            t.clean_svid = query_svid
            assert query_svid in parsed_id_dict, f"{t.svid} {r}"
            t.readids = parsed_id_dict[query_svid]

            asm_ol_readn = set(t.readids) & rids
            read_num_wt_type = 0
            t_type_flag = simple_type(t.SVTYPE, t.ctg, t.ctg2)
            asm_read_ids = []
            for read_i in asm_ol_readn:
                assert read_i in rids_withtype

                # one read with multiple sv type
                # read_ids_withtype: asm sv results
                for sv_j in rids_withtype[read_i]:
                    # translocation
                    if t_type_flag == 8 and sv_j.flag == 8:
                        read_num_wt_type += 1
                        asm_read_ids.append(read_i)
                        break

                    # 1. >100k(cutoff) DUP/INS/DEL/INV can be translocation in self-assembly
                    if t_type_flag in [1, 2, 4] and abs(t.SVLEN) >= 1e5 and sv_j.flag == 8: 
                        read_num_wt_type += 1
                        asm_read_ids.append(read_i)
                        break

                    len_check = (
                        abs(sv_j.len) >= abs(t.SVLEN) * opt.min_len_ratio
                        and abs(t.SVLEN) >= abs(sv_j.len) * opt.min_len_ratio
                    ) or abs(abs(sv_j.len) - abs(t.SVLEN)) <= 1000

                    # strict sv type and length check
                    if (t_type_flag & sv_j.flag) and len_check:
                        read_num_wt_type += 1
                        asm_read_ids.append(read_i)
                        break

                    # patch complex BND such as template insertion
                    # or self-assembly intra-chrom BND but nanomonsv reports INV
                    if (t_type_flag == 0 or sv_j.flag == 0) and len_check:
                         read_num_wt_type += 1
                         asm_read_ids.append(read_i)
                         break

                    # 2. <1K INDEL could change type INS->DEL, DEL->INS
                    #    partly caused by VNTR
                    if t_type_flag in [1, 2] and sv_j.flag in [1, 2]:
                        if (t_type_flag == 1 and sv_j.flag == 2) or (t_type_flag == 2 and sv_j.flag == 1):
                            len_check1 = abs(t.SVLEN) + abs(sv_j.len) <= 1000
                            #ignore the SVTYPE for the filtering in this step
                            if len_check1:
                                read_num_wt_type += 1
                                asm_read_ids.append(read_i)
                                break
                            else:
                                # duplicated check
                                len_check2 = (
                                    abs(sv_j.len) >= abs(t.SVLEN) * opt.min_len_ratio
                                    and abs(t.SVLEN) >= abs(sv_j.len) * opt.min_len_ratio
                                )
                                if len_check2 and (t_type_flag & sv_j.flag):
                                    read_num_wt_type += 1
                                    asm_read_ids.append(read_i)
                                    break
            t.asmreadids = asm_read_ids
            sv.append(t)

    sv.sort(key=attrgetter('ctg', 'pos'))
    group = []
    n_ambi1 = 0
    n_ambi2 = 0
    # this follows the merge operation
    for i, v in enumerate(sv):
        merge_to = []
        is_ambi1 = False
        for j in range(len(group)-1, -1, -1):
            g = group[j]
            if g[0].ctg != v.ctg:
                 break
            # NOTE: maybe g[-1], to compare the closet one
            if v.pos - g[0].pos > opt.win_size:
                 break
            n_same = 0
            for k, gv in enumerate(g):
                if gc_cmp_same_sv1(opt.win_size, opt.min_len_ratio, g[k], sv[i]):
                    n_same += 1
            if n_same > 0:
                merge_to.append(j)
                # NOTE: save sv number not equal to group size
                if n_same != len(g):
                    is_ambi1 = True

        if is_ambi1:
            n_ambi1 += 1

        if len(merge_to) == 0:
            group.append([v])
        else:
            # NOTE: merge to the closet group
            group[merge_to[0]].append(v)
            # NOTE: two possible merge group
            if len(merge_to) > 1:
                n_ambi2 += 1
    # summarize the results
    cnt = []
    for i in range(0, 1<<len(msvs)):
        cnt.append(0)
    for j in range(len(group)):
        g = group[j]
        in_bed = False
        has_bnd = False
        max_len = 0
        max_cnt = 0
        x = 0
        for k in range(len(g)):
            if g[k].in_bed:
                in_bed = True

            max_cnt = max(max_cnt, g[k].count)
            max_len = max(max_len, abs(g[k].SVLEN))
            if g[k].SVTYPE == "BND": 
                has_bnd = True
            x |= 1<<g[k].file_id
            g[k].group_id = j
        # NOTE: different sv tools has different metrics of counts
        #       computed the max of each group to the threshold
        if max_cnt < opt.group_min_count:
            continue
        # non-tranlocation length condition
        if has_bnd == False and max_len < opt.min_len:
            continue
        if not in_bed:
            continue
        cnt[x] += 1
        if opt.print_sv:
            if opt.collapsed:
                max_k = np.argmax([ len(g[k].asmreadids) for k in range(len(g)) ])
                if len(g[max_k].asmreadids) >= opt.group_min_count:
                    print("\t".join(map(str, gc_sv2array(g[max_k]))))
            else:
                for k in range(len(g)):
                    if len(g[k].asmreadids) >= opt.group_min_count:
                        print("\t".join(map(str, gc_sv2array(g[k]))))


def parse_tr(msv):
    from intervaltree import Interval, IntervalTree
    tree_dict = {}
    with open(msv) as fin:
        for line in fin:
            line = line.strip().split('\t')
            if line[0] not in tree_dict:
                tree_dict[line[0]] = IntervalTree() 
            if line[1] != line[4]:
                tree_dict[line[0]].addi(int(line[1]), int(line[4]), line)
            else:
                tree_dict[line[0]].addi(int(line[1]), int(line[4])+1, line)
    return tree_dict


def parse_segdup(msv):
    from intervaltree import Interval, IntervalTree
    tree_dict = {}
    with open(msv) as fin:
        for line in fin:
            line = line.strip().split('\t')
            if line[0] not in tree_dict:
                tree_dict[line[0]] = IntervalTree() 
            tree_dict[line[0]].addi(int(line[1]), int(line[2])+1, line)
    return tree_dict


def union_sv_with_tr(msvs, read_min_len, opt):
    if opt.bed is not None:
        opt.bed = gc_read_bed(opt.bed)

    segdup_dict = parse_segdup(msvs[-1])

    tr_dict = parse_tr(msvs[-2])

    sv = []

    # cat all the svs
    for i, m in enumerate(msvs[:-2]):
        v =  gc_parse_sv(m, read_min_len, opt.read_min_count, False, False)
        for j, vj in enumerate(v):
            t = vj
            t.file_id = i
            t.in_bed = True
            if opt.bed is not None:
                if t.ctg not in opt.bed or t.ctg2 not in opt.bed:
                    t.in_bed = False
                elif len(iit_overlap(opt.bed[t.ctg], t.pos, t.pos+1)) == 0:
                    t.in_bed = False
                elif len(iit_overlap(opt.bed[t.ctg2], t.pos2, t.pos2+1)) == 0:
                    t.in_bed = False
            sv.append(t)
   
    sv.sort(key=attrgetter('ctg', 'pos'))

    group = []
    n_ambi1 = 0
    n_ambi2 = 0

    # this follows the merge operation
    for i, v in enumerate(sv):
        merge_to = []
        is_ambi1 = False
        for j in range(len(group)-1, -1, -1):
            g = group[j]
            if g[0].ctg != v.ctg:
                 break
            # NOTE: maybe g[-1], to compare the closet one
            if v.pos - g[0].pos > opt.win_size:
                 break
            n_same = 0
            for k, gv in enumerate(g):
                if gc_cmp_same_sv1(opt.win_size, opt.min_len_ratio, g[k], sv[i]):
                    n_same += 1
            if n_same > 0:
                merge_to.append(j)
                # NOTE: save sv number not equal to group size
                if n_same != len(g):
                    is_ambi1 = True

        if is_ambi1:
            n_ambi1 += 1

        if len(merge_to) == 0:
            group.append([v])
        else:
            # NOTE: merge to the closet group
            group[merge_to[0]].append(v)
            # NOTE: two possible merge group
            if len(merge_to) > 1:
                n_ambi2 += 1

    # summarize the results
    cnt = []
    for i in range(0, 1<<len(msvs[:-2])):
        cnt.append(0)

    read_cnt = []
    sv_lens = []
    for i in range(0, len(msvs[:-2])):
        read_cnt.append('NA')
        sv_lens.append('NA')

    within_tr = 0
    total_sv = 0

    for j in range(len(group)):
        g = group[j]
        in_bed = False
        has_bnd = False
        max_len = 0
        max_cnt = 0
        x = 0
        test1, test2, test3 = [], [], []
        for k in range(len(g)):
            if g[k].in_bed:
                in_bed = True

            max_cnt = max(max_cnt, g[k].count)
            max_len = max(max_len, abs(g[k].SVLEN))
            if g[k].SVTYPE == "BND": 
                has_bnd = True
            x |= 1<<g[k].file_id
            g[k].group_id = j

            # left breakpoint test
            if g[k].ctg in tr_dict:
               test1 = tr_dict[g[k].ctg].overlap(g[k].pos-25, g[k].pos+25)
            # right breakpoint test
            if g[k].ctg2 in tr_dict:
               test2 = tr_dict[g[k].ctg2].overlap(g[k].pos2-25, g[k].pos2+25)
            if g[k].ctg2 is None:
               # indel test
               test3 = tr_dict[g[k].ctg].overlap(g[k].st-25, g[k].en+25)
            if len(test1) > 0 or len(test2) > 0 or len(test3) > 0:
               test_trs = list(test1)+list(test2)+list(test3)
               g[k].tr_info = test_trs[len(test_trs)>>1].data

            if g[k].ctg in segdup_dict:
               test1 = segdup_dict[g[k].ctg].overlap(g[k].pos-25, g[k].pos+25)
            # right breakpoint test
            if g[k].ctg2 in segdup_dict:
               test2 = segdup_dict[g[k].ctg2].overlap(g[k].pos2-25, g[k].pos2+25)
            if g[k].ctg2 is None:
               # indel test
               test3 = segdup_dict[g[k].ctg].overlap(g[k].st-25, g[k].en+25)
            if len(test1) > 0 or len(test2) > 0 or len(test3) > 0:
               test_trs = list(test1)+list(test2)+list(test3)
               g[k].segdup_info = test_trs[len(test_trs)>>1].data

        # NOTE: different sv tools has different metrics of counts
        #       computed the max of each group to the threshold
        if max_cnt < opt.group_min_count:
            continue
        # non-tranlocation length condition
        if has_bnd == False and max_len < opt.min_len:
            continue
        if not in_bed:
            continue
        cnt[x] += 1
        if opt.print_sv:
            total_sv += 1
            tr_info_value = ""
            for k in range(len(g)):
                if g[k].tr_info != "":
                    tr_info_value = g[k].tr_info
                ##print("\t".join(map(str, gc_sv2array(g[k]))), tr_info_value, sep='\t')
                read_cnt[int(g[k].file_id)] = g[k].count

                if g[k].ctg2 is not None and g[k].ctg != g[k].ctg2:
                    sv_len_ind = np.inf
                else:
                    sv_len_ind = g[k].SVLEN

                sv_lens[int(g[k].file_id)] = sv_len_ind

            max_k = np.argmax([ g[k].count for k in range(len(g)) ])
            if tr_info_value != "":
                within_tr += 1

            print("\t".join(map(str, gc_sv2array(g[max_k]))), ','.join(map(str, sv_lens)), ','.join(map(str, read_cnt)), sep='\t')

#    if not opt.print_sv:
#        for x in range(1, len(cnt)):
#            label = []
#            for i in range(len(msvs)):
#                label.append(x>>i&1)
#            print("".join(map(str, label)), cnt[x], sep="\t")
