import gzip
import pandas as pd
import csv
import time
import os
import sys
import mappy as mp
from .regex import re_info
import math
import warnings
import re
from .util import is_gzipped, is_vcf
from .eval import (
    gc_parse_sv,
    iit_overlap,
    gc_read_bed,
    iit_sort_copy,
    iit_index,
    eval1,
)
from .merge import parse_sv, svinfo

from .type import get_type, simple_type
from pathlib import Path
import subprocess
from minisv.read_parser import load_reads
from .ensemble import insilico_truth


# The categories to group by
categories = ["l+s", "l+g+s", "l+g"]
callers = ["severus", "savana", "nanomonsv"]


def parse_readids(file_path):
    tsv = file_path
    read_id_dict = {}
    caller = 'snf' if is_vcf(tsv) or is_gzipped(file_path) else ''

    if is_gzipped(file_path):
        f = gzip.open(file_path, 'rt')
    else:
        f = open(file_path)

    line_no = 0
    for line in f:
        line_no += 1
        if caller == 'snf':
            if line.startswith("#"):
                continue
        if line.startswith("#SV_ID"):
            caller = 'severus'
            continue
        if line.startswith('VARIANT_ID'):
            caller = 'savana'
            continue
        if line.startswith('chr') and line.strip()[-1] in ['+', '-']:
            caller = 'nanomonsv'
        if 'avg_mapq' in line:
            caller = 'msv'

        if caller == 'severus':
            line = line.strip().split(',')
        elif caller == 'savana':
            line = line.strip().split('\t')
        elif caller == 'nanomonsv':
            line = line.strip().split('\t')
        elif caller == 'snf':
            line = line.strip().split('\t')
        elif caller == 'msv':
            line = line.strip().split('\t')
        else:
            raise Exception("not support yet")

        # only query the tumor SVs supported reads from the read id tsv file
        if caller == 'severus':
            read_id_dict[line[0]] = line[1].replace("\"", "").split('; ')
        elif caller == 'savana':
            read_id_dict[line[0]] = line[1].replace("\"", "").split(',')
        elif caller == 'nanomonsv':
            read_id_dict[line[7]] = read_id_dict.get(line[7], []) + [line[8]]
        elif caller == 'snf':
            info = re_info.findall(line[7])
            snf_reads = []
            for info_field, info_val in info:
                 if info_field.startswith('RNAMES'):
                     snf_reads = info_val.split(',')
                     break
            read_id_dict[line[2]] = snf_reads
        elif caller == 'msv':
            is_bp = False
            if re.match(r"[><]", line[2]):
                is_bp = True
            col_info = 8 if is_bp else 6
            info = line[col_info]
            info = re_info.findall(info)
            msv_reads = []
            for info_field, info_val in info:
                 if info_field.startswith('reads'):
                     msv_reads = info_val.split(',')
                     break
            # msv use line number as sv ids
            read_id_dict[str(line_no)] = msv_reads
        else:
            raise Exception("not support yet")
    f.close()
    return read_id_dict


def parse_msvasm(msvasm):
    with gzip.open(msvasm) as gzip_file:
        read_ids = []
        for line in gzip_file:
            line = line.strip().split()
            is_bp = False
            if re.match(r"[><]", line[2].decode('utf-8')):
                is_bp = True
            col_info = 8 if is_bp else 6
            read_ids.append(line[col_info-3].decode('utf-8'))

    read_ids = set(read_ids)
    return read_ids


def parse_msvasm_withtype(msvasm):

    h = {}
    with gzip.open(msvasm, 'rt') as gzip_file:
        read_ids = []
        for line in gzip_file:
            t = line.strip().split('\t')
            is_bp = False
            if re.match(r"[><]", t[2]):
                is_bp = True

            col_info = 8 if is_bp else 6
            name = t[col_info-3]
            if name not in h:
                h[name] = []
            h[name].append(get_type(t, col_info))
    return h


def classify_sv_len(svlen=0):
    if svlen == 0:
        sv_len_tag = 'translocation'
    elif svlen >= 1e6:
        sv_len_tag = '>1M'
    elif svlen < 1e6 and svlen >= 1e5:
        sv_len_tag = '100kb-1M'
    elif svlen < 1e5 and svlen >= 20e3:
        sv_len_tag = '20kb-100kb'
    elif svlen < 20e3 and svlen >= 100:
        sv_len_tag = '100bp-20kb'
    else:
        sv_len_tag = '<100bp'
    return sv_len_tag


def parse_svid(svid):
    svid = str(svid)
    # remove _1 and _2 for paired SV event
    # severus
    if svid.startswith('severus'):
        query_svid = re.sub("_[12]", "", svid)
    # savana
    elif svid.startswith("ID_"):
        query_svid = re.sub(r"(ID_\d+)_([1,2])", "\\1", svid)
    # nanomonsv
    elif svid.startswith("r_") or svid.startswith("i_") or svid.startswith("d_"):
        query_svid = re.sub(r"(r_\d+)_([0,1])", "\\1", svid)
    # sniffles2
    elif svid.startswith("Sniffles2"):
        query_svid = svid
    else:
        query_svid = svid
    return query_svid
        #raise Exception("not support yet")

def gnomad_filter(vcf_file, gnomad_bed, opt, both_ends=False, pad=0, out=None):
    """Drop caller VCF calls whose breakpoint overlaps a gnomAD-SV BED region.

    Type-agnostic breakpoint overlap (no allele-frequency parsing). By default a
    call is dropped if *either* breakpoint overlaps the BED; both_ends=True
    requires both. The original VCF is re-emitted to ``out`` minus the dropped
    calls; lines that gc_parse_sv never parses pass through unchanged.
    """
    if out is None:
        out = sys.stdout

    bed = gc_read_bed(gnomad_bed, is_gnomad=True)
    # parse permissively (min_len=0, min_count=0) so gnomAD overlap is the only filter
    # svs = gc_parse_sv(vcf_file, 0, 0, opt.ignore_flt, False)

    min_read_len = math.floor(opt.min_len * opt.read_len_ratio + 0.499)
    svs = gc_parse_sv(vcf_file, min_read_len, opt.min_count, opt.ignore_flt, opt.check_gt)

    drop = set()
    included = set()
    for t in svs:
        # Not long enough for non-BND type
        # NOTE: here use 100bp as the default length instead of 80
        if t.SVTYPE != "BND" and abs(t.SVLEN) < opt.min_len:
            continue

        if t.SVTYPE == "BND" and t.ctg == t.ctg2 and abs(t.SVLEN) < opt.min_len:
            continue

        if t.count > 0 and t.count < opt.min_count:
            continue

        # filter by VAF
        if t.vaf is not None and t.vaf < opt.min_vaf:
            continue

        if opt.bed is not None:
            if t.ctg not in opt.bed or t.ctg2 not in opt.bed:
                continue
            if len(iit_overlap(opt.bed[t.ctg], t.pos, t.pos + 1)) == 0:
                continue
            if len(iit_overlap(opt.bed[t.ctg2], t.pos2, t.pos2 + 1)) == 0:
                continue

        # all svs passed QC
        included.add(parse_svid(t.svid))
        hit1 = t.ctg in bed and len(iit_overlap(bed[t.ctg], t.pos - pad, t.pos + 1 + pad)) > 0
        hit2 = t.ctg2 in bed and len(iit_overlap(bed[t.ctg2], t.pos2 - pad, t.pos2 + 1 + pad)) > 0
        if (hit1 and hit2) if both_ends else (hit1 or hit2):
           # gnomad filtered svs
            drop.add(parse_svid(t.svid))

    if is_gzipped(vcf_file):
        f = gzip.open(vcf_file, 'rt')
    else:
        f = open(vcf_file)

    n = 0
    for line in f:
        if line[0] == "#":
            print(line.strip(), file=out)
            continue
        t = line.strip().split("\t")
        if len(t) >= 3 and parse_svid(t[2]) in drop:
            n += 1
            continue
        if parse_svid(t[2]) in included:
            print(line.strip(), file=out)

    print(len(included), len(drop), file=sys.stderr)
    f.close()


def _parse_union_msv(path):
    """Parse the 6-column .msv format that union_sv emits (via gc_sv2array).

    Format: ctg \\t pos \\t ori \\t ctg2 \\t pos2 \\t info_string
    where info_string carries SVTYPE/SVLEN/count/group_id/file_id/... .

    gc_parse_sv cannot read this format — it requires 9 columns (col_info=8)
    for breakpoint records. This is a tiny local parser that yields the
    (svinfo, raw_line) pairs needed by _consensus_lost_by_filter.
    """
    pairs = []
    if not os.path.exists(path):
        return pairs
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 6:
                continue
            info = t[5]
            svtype_m = re.search(r"SVTYPE=([^\s;]+)", info)
            svlen_m = re.search(r"SVLEN=(-?\d+)", info)
            svtype = svtype_m.group(1) if svtype_m else None
            svlen = int(svlen_m.group(1)) if svlen_m else 0
            s = svinfo(
                ctg=t[0], pos=int(t[1]), ctg2=t[3], pos2=int(t[4]),
                ori=t[2], SVTYPE=svtype, SVLEN=svlen, count=1, vaf=1,
            )
            pairs.append((s, line.rstrip("\n")))
    return pairs


def _consensus_lost_by_filter(consensus_msv, target_msv, opt, file_handler=None):
    """Emit lines from consensus_msv whose SV overlaps zero rows in target_msv.

    Pure set-difference. Matching reuses the IIT lookup pattern from
    gc_cmp_sv: a target's pos and pos2 are both indexed, so a consensus
    rep matches if EITHER of its breakpoints falls within opt.win_size of
    a target breakpoint AND gc_cmp_same_sv1 accepts the pair. This
    correctly handles BNDs whose second breakpoint lies on a different
    contig.

    Inputs are union_sv's 6-col .msv format (parsed via _parse_union_msv,
    not gc_parse_sv — see that helper's docstring). Eval-style call-quality
    filters (min_len/min_count/min_vaf/opt.bed) are intentionally NOT
    applied here: the consensus is already filtered upstream by
    union_sv + insilico_truth.
    """
    target_pairs = _parse_union_msv(target_msv)
    h = {}
    for s, _raw in target_pairs:
        if s.ctg not in h:
            h[s.ctg] = []
        if s.ctg2 not in h:
            h[s.ctg2] = []
        h[s.ctg].append({"st": s.pos, "en": s.pos + 1, "data": s})
        h[s.ctg2].append({"st": s.pos2, "en": s.pos2 + 1, "data": s})
    for ctg in h:
        h[ctg] = iit_sort_copy(h[ctg])
        iit_index(h[ctg])

    for r, raw in _parse_union_msv(consensus_msv):
        n = eval1(opt, h, r.ctg, r.pos, r) + eval1(opt, h, r.ctg2, r.pos2, r)
        if n == 0:
            print(raw, file=file_handler)


def othercaller_filterasm(vcf_file, opt, readidtsv, msvasm, outstat, consensus_sv_ids, out_filtered_vcf = None):
    """
    vcf: input vcf, support Severus, Sniffles, nanomonsv, and savana
    consensus_sv_ids: consensus SV ids from minisv annot
    """
    # NOTE: why some germline SV disappear and no somatic SV
    parsed_id_dict = parse_readids(readidtsv)
    read_ids = parse_msvasm(msvasm)
    read_ids_withtype = parse_msvasm_withtype(msvasm)

    min_read_len = math.floor(opt.min_len * opt.read_len_ratio + 0.499)
    vcf = gc_parse_sv(vcf_file, min_read_len, opt.min_count, opt.ignore_flt, opt.check_gt)

    consensus_ids = []
    if consensus_sv_ids != "":
        with open(consensus_sv_ids) as inf:
            for line in inf:
                consensus_ids.append(re.sub("_[12]", "", line.strip()))

    is_consensus = True
    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tsvid\tis_consensus\ttype\tcontig1\tpos1\tori\tcontig2\tpos2\tasm_support\tasm_support_onlyreadname\tread_name_number\n")

    # first vcf is to be annotated
    # read id, typing and length
    all_ids = []
    # only read name overlap
    all_ids_onlyname = []
    # msv sv id record
    for i in range(len(vcf)):
 
        t = vcf[i]

        # Not long enough for non-BND type
        # NOTE: here use 100bp as the default length instead of 80
        if t.SVTYPE != "BND" and abs(t.SVLEN) < opt.min_len:
            continue

        if t.SVTYPE == "BND" and t.ctg == t.ctg2 and abs(t.SVLEN) < opt.min_len:
            continue

        if t.count > 0 and t.count < opt.min_count:
            continue

        # filter by VAF
        if t.vaf is not None and t.vaf < opt.min_vaf:
            continue

        if opt.bed is not None:
            if t.ctg not in opt.bed or t.ctg2 not in opt.bed:
                continue
            if len(iit_overlap(opt.bed[t.ctg], t.pos, t.pos + 1)) == 0:
                continue
            if len(iit_overlap(opt.bed[t.ctg2], t.pos2, t.pos2 + 1)) == 0:
                continue

        query_svid = str(parse_svid(t.svid))
        
        if len(consensus_ids) > 0:
            is_consensus = query_svid in consensus_ids
        assert query_svid in parsed_id_dict, f"{query_svid}, {vcf_file}"

        ol_readn = set(parsed_id_dict[query_svid]) & read_ids

        read_num_wt_type = 0
        filt_cnt = len(ol_readn)
        t_type_flag = simple_type(t.SVTYPE, t.ctg, t.ctg2)
        for read_i in ol_readn:
            assert read_i in read_ids_withtype

            # one read with multiple sv type
            # read_ids_withtype: asm sv results
            for sv_j in read_ids_withtype[read_i]:
                # translocation
                if t_type_flag == 8 and sv_j.flag == 8:
                    read_num_wt_type += 1
                    break

                # 1. >100k(cutoff) DUP/INS/DEL/INV can be translocation in self-assembly
                if t_type_flag in [1, 2, 4] and abs(t.SVLEN) >= 1e5 and sv_j.flag == 8: 
                    read_num_wt_type += 1
                    break

                len_check = (abs(abs(sv_j.len) - abs(t.SVLEN)) <= opt.max_diff_len) or (
                    abs(sv_j.len) >= abs(t.SVLEN) * opt.min_len_ratio
                    and abs(t.SVLEN) >= abs(sv_j.len) * opt.min_len_ratio
                ) 
                # strict sv type and length check
                if (t_type_flag & sv_j.flag) and len_check:
                    read_num_wt_type += 1
                    break

                # patch complex BND such as template insertion
                # or self-assembly intra-chrom BND but nanomonsv reports INV
                if (t_type_flag == 0 or sv_j.flag == 0) and len_check:
                     read_num_wt_type += 1
                     break

                # 2. <1K INDEL could change type INS->DEL, DEL->INS
                #    partly caused by VNTR
                if t_type_flag in [1, 2] and sv_j.flag in [1, 2]:
                    if (t_type_flag == 1 and sv_j.flag == 2) or (t_type_flag == 2 and sv_j.flag == 1):
                        len_check1 = abs(t.SVLEN) + abs(sv_j.len) <= opt.max_diff_len
                        #ignore the SVTYPE for the filtering in this step
                        if len_check1:
                            read_num_wt_type += 1
                            break
                        else:
                            # duplicated check
                            len_check2 = (
                                abs(sv_j.len) >= abs(t.SVLEN) * opt.min_len_ratio
                                and abs(t.SVLEN) >= abs(sv_j.len) * opt.min_len_ratio
                            )
                            if len_check2 and (t_type_flag & sv_j.flag):
                                read_num_wt_type += 1
                                break

            assert read_num_wt_type <= filt_cnt
            if t.svid == opt.svid and opt.svid != "":
                print(t.svid, read_i, sv_j.flag, t_type_flag, sv_j.len, t.SVLEN, sv_j.type, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio, read_num_wt_type)

        if opt.print_all:
            outf.write('\t'.join(map(str, [t.SVLEN, classify_sv_len(abs(t.SVLEN)), t.svid, is_consensus, t.SVTYPE, t.ctg, t.pos, t.ori, t.ctg2, t.pos2, 
                       read_num_wt_type, filt_cnt, t.count]))+'\n')

         # output sv ids in vcf or msv
        if read_num_wt_type >= opt.min_count and filt_cnt >= opt.min_count and t.count >= opt.min_count:
            all_ids.append(query_svid)
        if opt.only_readname:
            if filt_cnt >= opt.min_count and t.count >= opt.min_count:
                all_ids_onlyname.append(query_svid)

    outf.close()

    if opt.svid != "":
        return

    assert set(all_ids).issubset(set(parsed_id_dict.keys()))
    assert set(all_ids_onlyname).issubset(set(parsed_id_dict.keys()))
    ## output vcf format with sv ids above
    if is_gzipped(vcf_file):
        f = gzip.open(vcf_file, 'rt')
    else:
        f = open(vcf_file)
    line_no = 0
    for line in f:
        line_no += 1
        if line[0] == "#":
            print(line.strip(), file=out_filtered_vcf)
            continue
        if 'avg_mapq' in line: # MSV
            query_svid = str(parse_svid(line_no))
        else: # other caller
            t = line.strip().split("\t")
            query_svid = str(parse_svid(t[2]))
        if opt.only_readname:
            if query_svid in all_ids_onlyname:
                print(line.strip())
        else:
            if query_svid in all_ids:
                print(line.strip(), file=out_filtered_vcf)
    f.close()


def call_filterseverus(severusvcf, readidtsv, msvasm, outstat, consensus_sv_ids, asm_count_cutoff=2):
    parsed_id_dict = parse_readids(readidtsv)
    read_ids = parse_msvasm(msvasm)
    consensus_ids = []
    with open(consensus_sv_ids) as inf:
        for line in inf:
            consensus_ids.append(re.sub("_[12]", "", line.strip()))

    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tsvid\tis_consensus\tDV\tread_name_number\tasm_support\n")
    vcf_svids = []
    inconsistent_read_num = 0
    total = 0
    with open(severusvcf) as inf:
        for line in inf:
            svlen = 0
            if line.startswith("#"):
                print(line.strip())
                continue

            total += 1
            elements = line.strip().split()
            info = elements[-3].split(';')
            for info_field in info:
                 if info_field.startswith('SVLEN='):
                     svlen = int(info_field.replace('SVLEN=', ''))

            sv_len_tag = classify_sv_len(abs(svlen))
            
            ## Severus FORMAT column: GT:GQ:VAF:hVAF:DR:DV
            svid = elements[2]
            svid = re.sub("_[12]", "", svid)
            vcf_svids.append(svid)

            # if there is no read information, do not output the SVs
            assert svid in parsed_id_dict
            form = elements[-1]
            form_list = form.split(':')
            assert int(form_list[-1]) >= len(set(parsed_id_dict[svid])), f"{svid} {form_list[-1]} {len(set(parsed_id_dict[svid]))} {form_list}"
            inconsistent_read_num += int(form_list[-1]) != len(set(parsed_id_dict[svid]))

            if int(form_list[-1]) != len(set(parsed_id_dict[svid])):
                 w_mess = f"{svid} {int(form_list[-1])} {len(set(parsed_id_dict[svid]))} unequal reads between output read names and vcf DV"
                 warnings.warn(w_mess, UserWarning)

            # if there is more than 2 SVs with self-assembly support, output the SVs
            if len(set(parsed_id_dict[svid]) & read_ids) >= asm_count_cutoff:
                print(line.strip())

            if elements[2].endswith('_2'): # BND_2, for paired events such as translocation, output one end for read filtering
                continue

            ## SVLEN, SVID, DV total tumor SV reads, read TSV file tumor SV read num, asm supported SV
            if sv_len_tag != '<100bp':
                outf.write('\t'.join([str(svlen), sv_len_tag, svid, str(svid in consensus_ids), form_list[-1], str(len(set(parsed_id_dict[svid]))), str(len(set(parsed_id_dict[svid]) & read_ids))])+'\n')

    assert set(vcf_svids).issubset(set(parsed_id_dict.keys()))
    outf.close()


def call_filtersnf(snfvcfgz, msvasm, outstat, asm_count_cutoff=2):
    # Sniffle2 suport, DV, read names issues:
    #   https://github.com/fritzsedlazeck/Sniffles/issues/496
    #   https://github.com/fritzsedlazeck/Sniffles/issues/382
    read_ids = parse_msvasm(msvasm)

    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tsupport\tsvid\tread_name_number\tDV\tasm_support\n")
    with gzip.open(snfvcfgz) as inf:
        for line in inf:
            line = line.decode('utf-8')
            svlen = 0
            support = 0
            snf_reads = []

            if line.startswith("#"):
                print(line.strip())
                continue
            elements = line.strip().split()
            svid = elements[2]
            # -2: tumor, -1: normal
            # GT:GQ:DR:DV:ID
            fmt = elements[8].split(':')
            tumor_form  = elements[-2]
            tumor_form_list = tumor_form.split(':')
            normal_form = elements[-1]
            normal_form_list = normal_form.split(':')

            info = elements[-4].split(';')
            for info_field in info:
                 if info_field.startswith('SVLEN='):
                     svlen = int(info_field.replace('SVLEN=', ''))
            sv_len_tag = classify_sv_len(abs(svlen))

            for info_field in info:
                 if info_field.startswith('RNAMES'):
                     snf_reads = info_field.replace('RNAMES=', '').split(',')
                 if info_field.startswith('SUPPORT'):
                     support = info_field.replace('SUPPORT=', '')

            for i in range(len(fmt)):
                if fmt[i] == 'DV':
                   ##assert len(set(snf_reads)) == int(tumor_form_list[-3]) + int(tumor_form_list[-2]) + int(normal_form_list[-2]) + int(normal_form_list[-3]) 
                   if int(normal_form_list[i]) == 0 and int(tumor_form_list[i]) >= asm_count_cutoff:
                       # SVID, support, svid, read name number, DV, asm supported reads
                       if sv_len_tag != '<100bp':
                           outf.write('\t'.join([str(svlen), sv_len_tag, support, svid, str(len(set(snf_reads))), tumor_form_list[-2], str(len(read_ids & set(snf_reads)))])+'\n')
                       # if there is more than 2 SVs with self-assembly support, output the SVs
                       if len(read_ids & set(snf_reads)) >= asm_count_cutoff:
                           print(line.strip())
    outf.close()


def call_filtermsv(msvtg, msvasm, outstat, asm_count_cutoff=2):

    read_ids = parse_msvasm(msvasm)

    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tread_name_number\tasm_support\n")
    msv_all_reads = []
    with open(msvtg) as inf:
        for line in inf:
            elements = line.strip().split()
            info = elements[-1].split(';')
            svlen = 0
            msv_reads = []

            for info_field in info:
                 if info_field.startswith('reads='):
                     msv_reads = info_field.replace('reads=', '').split(',')
                 if info_field.startswith('SVLEN='):
                     svlen = int(info_field.replace('SVLEN=', ''))
            msv_all_reads += msv_reads

            sv_len_tag = classify_sv_len(abs(svlen))
            if sv_len_tag != '<100bp':
                 # SVID, support, svid, read name number, DV, asm supported reads
                 outf.write('\t'.join([str(svlen), sv_len_tag, str(len(set(msv_reads))), str(len(read_ids & set(msv_reads)))])+'\n')
            if len(read_ids & set(msv_reads)) >= asm_count_cutoff:
                 print(line.strip())
    ##assert set(msv_all_reads).issubset(read_ids)
    outf.close()



def hit_to_paf(read_name, read_seq, hit):
    """
    Convert a mappy.Alignment object into a PAF line.
    """

    # Required PAF fields
    qname = read_name
    qlen = len(read_seq)
    qstart = hit.q_st
    qend = hit.q_en

    rname = hit.ctg
    rlen = hit.ctg_len
    rstart = hit.r_st
    rend = hit.r_en

    strand = "+" if hit.strand == 1 else "-"

    # Optional: minimap2-style tags
    tp = "P" if hit.is_primary else "S"       # primary or secondary
    nm = hit.NM if hit.NM is not None else 0
    ms = hit.mlen                            # number of matching bases
    bl = hit.blen                            # alignment block length
    mapq = hit.mapq
    cigar = hit.cigar_str or "*"
    md = hit.MD or ""
    cs = hit.cs or ""
    ds = hit.ds or ""

    # Build PAF line
    fields = [
        qname, qlen, qstart, qend,
        strand,
        rname, rlen, rstart, rend,
        ms, bl, mapq,
        f"NM:i:{nm}",
        f"ms:i:{ms}",
        f"AS:i:{ms - nm}",          # crude score estimate
        f"nn:i:0",
        f"tp:A:{tp}",
        f"cg:Z:{cigar}",
    ]

    if md:
        fields.append(f"MD:Z:{md}")
    if cs:
        fields.append(f"cs:Z:{cs}")
    if ds:
        fields.append(f"ds:Z:{ds}")

    return "\t".join(map(str, fields))


import time
import csv
import psutil
import os
import gc
from pathlib import Path

class Timer:
    """Context manager to track time and max memory usage per step"""
    def __init__(self, name, timings_list):
        self.name = name
        self.timings_list = timings_list
        self.process = psutil.Process(os.getpid())

    def __enter__(self):
        gc.collect()  # Clean up before measuring
        self.start_time = time.time()
        self.start_mem = self.process.memory_info().rss / 1024 / 1024  # MB
        self.max_mem = self.start_mem
        return self

    def __exit__(self, *args):
        gc.collect()
        elapsed = time.time() - self.start_time
        current_mem = self.process.memory_info().rss / 1024 / 1024  # MB

        # Update max memory if current is higher
        self.max_mem = max(self.max_mem, current_mem)

        record = {
            "step": self.name,
            "time_sec": round(elapsed, 3),
            "max_memory_mb": round(self.max_mem, 2)
        }
        self.timings_list.append(record)

        print(f"[TIMING] {self.name}: {elapsed:.2f}s | Max Memory: {self.max_mem:.1f} MB")


def isec(w, gsvs, file_handler=None):
    """Usage: minisv isec base.gsv alt.gsv [...]"""
    from .type import get_type

    g = []
    # file as filter
    for gsv in gsvs[1:]:
        h = {}
        with gzip.open(gsv, 'rt') as gsv_file:
            for line in gsv_file:
                t = line.strip().split("\t")
                if re.match(r"[><]", t[2]):
                    col_info = 8
                else:
                    col_info = 6
                name = t[col_info - 3]
                if name not in h:
                    h[name] = []
                h[name].append(get_type(t, col_info))
        g.append(h)
    # one read may contain somatic and germline sv
    # g: [{readname:[sv_dict]}]
    out = file_handler if file_handler is not None else None

    # file to be filtered
    with gzip.open(gsvs[0], 'rt') as base_gsv:
        for line in base_gsv:
            t = line.strip().split("\t")
            if re.match(r"[><]", t[2]):
                col_info = 8
            else:
                col_info = 6
            name = t[col_info - 3]
            x = get_type(t, col_info)
            n_found = 0

            for h in g:
                if name not in h:
                    break

                a = h[name]
                found = False
                for j in range(len(a)):
                    if x.st - w < a[j].en and a[j].st < x.en + w:
                        # NOTE: why a[i].flag & 8?? what if base sv is not translocation
                        if x.flag == 0 or (x.flag & a[j].flag) or (a[j].flag & 8):
                            found = True
                if not found:
                    break
                n_found += 1

            if n_found == len(g):
                print(line.strip(), file=out)


class MinisvReads:
    def __init__(self, som_vcfs, readid_tsvs, bam_path, ref, hap1_denovo_ref_path, hap2_denovo_ref_path, 
                 graph_ref_path = "/hlilab/hli/minigraph/HPRC-r2/CHM13-464.gfa.gz", work_dir="minisv_work", filtered_readcount_cutoff=2, platform='hifi', mm2="minimap2", mg="minigraph"):
        self.som_vcfs = som_vcfs
        self.readid_tsvs = readid_tsvs
        self.mm2 = mm2
        self.mg = mg

        self.bam_path = Path(bam_path)

        self.ref = Path(ref)
        self.hap1_denovo_ref_path = Path(hap1_denovo_ref_path)
        self.hap2_denovo_ref_path = Path(hap2_denovo_ref_path)
        self.graph_ref_path = Path(graph_ref_path)

        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.filtered_readcount_cutoff = filtered_readcount_cutoff
        self.platform = platform
        self.timings = []

    def extract_read_ids(self, bed, min_len, min_read_len, min_count, ignore_flt, check_gt):
        """ only extract somatic SV read ids """
        with Timer("extract_read_ids", self.timings):
            if bed is not None and not isinstance(bed, dict):
                bed = gc_read_bed(bed)

            self.read_ids_file = self.work_dir / "readid.names"
            som_read_ids = set()

            self.read_id_dict = {}
            for (readid_tsv, vcf) in zip(self.readid_tsvs, self.som_vcfs):
                assert os.path.exists(readid_tsv), f"{readid_tsv} not exists"
                assert os.path.exists(vcf), f"{vcf} not exists"

                self.read_id_dict |= parse_readids(readid_tsv)
                vcf = gc_parse_sv(vcf, min_read_len, min_count, ignore_flt, check_gt)

                for i in range(len(vcf)):
                    t = vcf[i]
                    # Not long enough for non-BND type
                    # NOTE: here use 100bp as the default length instead of 80
                    if t.SVTYPE != "BND" and abs(t.SVLEN) < min_len:
                        continue

                    if t.SVTYPE == "BND" and t.ctg == t.ctg2 and abs(t.SVLEN) < min_len:
                        continue

                    if t.count > 0 and t.count < min_count:
                        continue

                    if bed is not None:
                        if t.ctg not in bed or t.ctg2 not in bed:
                            continue
                        if len(iit_overlap(bed[t.ctg], t.pos, t.pos + 1)) == 0:
                            continue
                        if len(iit_overlap(bed[t.ctg2], t.pos2, t.pos2 + 1)) == 0:
                            continue

                    query_svid = str(parse_svid(t.svid))
                    assert query_svid in self.read_id_dict, 'somatic sv not in read id table'
                    som_read_ids.add(query_svid)

            read_names = set()
            for s in list(som_read_ids):
                read_names |= set(self.read_id_dict[s])
       
            with open(self.read_ids_file, "w") as fin:
                for i in sorted(list(read_names)):
                    fin.write(f"{i}\n")
   
    def extract_reads(self):
        """samtools fastq aln.bam | seqtk subseq - read-names.txt > reads.fq"""
        with Timer("extract_reads", self.timings):
            self.fastq_out = self.work_dir / "som_reads.fq.gz"
            if os.path.exists(self.fastq_out):
                return

            cmd = f"samtools fastq {self.bam_path} | seqtk subseq - {self.read_ids_file} | gzip  > {self.fastq_out}"
            subprocess.run(cmd, shell=True, check=True)
            print(f"Reads extracted to {self.fastq_out}")
            return

    def build_pooled_reference(self, out_fa="pooled_reference.fa"):
        with Timer("build_pooled_reference", self.timings):
            self.pooled_ref = self.work_dir / out_fa
            with open(self.pooled_ref, "w") as out:
                for ref in [self.hap1_denovo_ref_path, self.hap2_denovo_ref_path]:
                    with gzip.open(ref, "rt") as f:
                        out.write(f.read())
            return mp.Aligner(str(self.pooled_ref)) # preset??

    def align_reads_to_self(self, paf='denovo_aligned.paf.gz'):
        with Timer("align_reads_to_denovo", self.timings):
            self.paf_out = self.work_dir / paf
            if os.path.exists(self.paf_out):
                return

            if self.platform == 'hifi':
                cmd = f"{self.mm2} --ds -t 4 -cx map-hifi -s50 <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            else:
                cmd = f"{self.mm2} --ds -t 4 -cx lr:hq <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            #cmd = f"{self.mm2} --ds -t 4 -cx lr:hq <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to self assembly")

    def build_reference(self):
        return mp.Aligner(str(self.ref))

    def align_reads_to_grch38(self, paf='grch38_aligned.paf.gz'):
        with Timer("align_reads_to_grch38", self.timings):
            self.grch38_paf_out = self.work_dir / paf
            if os.path.exists(self.grch38_paf_out):
                return

            if self.platform == 'hifi':
                cmd = f"{self.mm2} --ds -t 4 -cx map-hifi -s50 {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            else:
                cmd = f"{self.mm2} --ds -t 4 -cx lr:hq {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to grch38")

    def align_reads_to_graph(self, gaf='denovo_aligned.gaf.gz'):
        with Timer("align_reads_to_graph", self.timings):
            self.gaf_out = self.work_dir / gaf
            if os.path.exists(self.gaf_out):
                return
            cmd = f"{self.mg} -cxlr -t 4 {self.graph_ref_path} {self.fastq_out} | gzip - > {self.gaf_out}"
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"Reads aligned to {self.graph_ref_path} to generate {self.gaf_out}")
            return

    def parse_raw_sv_grch38(self, opt, gsv='grch38.gsv.gz',
                            filtered_grch38_gsv='grch38_filtered.gsv.gz'):
        with Timer("parse_raw_sv_grch38", self.timings):
            self.grch38_gsv_out = self.work_dir / gsv
            self.filtered_grch38_gsv_out = self.work_dir / filtered_grch38_gsv
            print("**********************")
            print(self.filtered_grch38_gsv_out)
            f = gzip.open(self.grch38_gsv_out, 'wt')
            load_reads(str(self.grch38_paf_out), opt, f)
            f.close()

            ##opt.bed = gc_read_bed(opt.bed) if b is not None else None
            print("**************")
            print(opt.cen)
            print(opt.bed)
            print("**************")

            if opt.bed is not None and not isinstance(opt.bed, dict):
                opt.bed = gc_read_bed(opt.bed)

            if opt.bed is not None or opt.cen is not None:
                f = gzip.open(self.filtered_grch38_gsv_out, 'wt')
                with gzip.open(self.grch38_gsv_out, "rt") as fin:
                    for line in fin:
                        if opt.cen is not None:
                            t = line.strip().split("\t")
                            v = parse_sv(t)

                            #NOTE: test if we do not use centromere filtering
                            ## Step 1. post filter grch38-based sv signals by centromere
                            if v.cen_overlap is not None and v.cen_overlap > 0:
                                continue
                            if v.cen_dist is not None and v.cen_dist <= 5e5:
                                continue

                            ### Step 2. filter by overlapping with high conf region
                            #if opt.bed is not None:
                            #    if v.ctg not in opt.bed or v.ctg2 not in opt.bed:
                            #        continue
                            #    if len(iit_overlap(opt.bed[v.ctg], v.pos, v.pos + 1)) == 0:
                            #        continue
                            #    if len(iit_overlap(opt.bed[v.ctg2], v.pos2, v.pos2 + 1)) == 0:
                            #        continue
                            # TODO Step 3: filtered by SV typing consistency between minisv/severus/nanomonsv...
                            print(line.strip(), file=f)
                f.close()
            return

    def parse_raw_sv_self(self, opt, gsv='denovo.gsv.gz'):
        with Timer("parse_raw_sv_denovo", self.timings):
            self.denovo_gsv_out = self.work_dir / gsv

            f = gzip.open(self.denovo_gsv_out, 'wt')
            load_reads(str(self.paf_out), opt, f)
            f.close()
            return

    def parse_raw_sv_graph(self, opt, gsv='graph.gsv.gz'):
        with Timer("parse_raw_sv_graph", self.timings):
            self.graph_gsv_out = self.work_dir / gsv

            f = gzip.open(self.graph_gsv_out, 'wt')
            load_reads(str(self.gaf_out), opt, f)
            f.close()
            return

    def isec_g(self, opt, w, msv='l+g.gsv.gz'):
        # l + g
        with Timer("isec_graph", self.timings):
            self.isec_graph_out = self.work_dir / msv

            f = gzip.open(self.isec_graph_out, 'wt')
            if opt.bed is not None or opt.cen is not None:
                isec(w, [self.filtered_grch38_gsv_out, self.graph_gsv_out], f)
            else:
                isec(w, [self.grch38_gsv_out, self.graph_gsv_out], f)
            f.close()
            return

    def isec_s(self, opt, w, msv='l+s.gsv.gz'):
        # l+s
        with Timer("isec_graph", self.timings):
            self.isec_denovo_out = self.work_dir / msv

            f = gzip.open(self.isec_denovo_out, 'wt')
            if opt.bed is not None or opt.cen is not None:
                isec(w, [self.filtered_grch38_gsv_out, self.denovo_gsv_out], f)
            else:
                isec(w, [self.grch38_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def isec_gs(self, opt, w, msv='gs.gsv.gz'):
        # l+s
        with Timer("isec_gs", self.timings):
            self.isec_gs_out = self.work_dir / msv

            f = gzip.open(self.isec_gs_out, 'wt')
            #if opt.bed is not None or opt.cen is not None:
            #    isec(w, [self.filtered_grch38_gsv_out, self.graph_gsv_out, self.denovo_gsv_out], f)
            #else:
            #    isec(w, [self.grch38_gsv_out, self.graph_gsv_out, self.denovo_gsv_out], f)
            isec(w, [self.graph_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def parse_ids_from_gsv(self, msv):
        with gzip.open(self.work_dir / msv, 'rt') as fin:
            for line in fin:
                t = line.strip().split('\t')

                is_bp = False
                if re.match(r"[><]", t[2]):
                    is_bp = True

                col_info = 8 if is_bp else 6
                name = t[col_info-3]
                yield name

    def othercaller_filterasm(self, opt):
        """ the same interface as the former filterasm 
        """

        # filter each caller separately instead of jointly 
        # joint filtering is error-prone
        # l_only = self.work_dir / Path("grch38_filtered.gsv.gz") # filtered by centromere dist
        # ls = self.work_dir / Path("l+s.gsv.gz")
        # lgs = self.work_dir / Path("l+g+s.gsv.gz")
        # lg = self.work_dir / Path("l+g.gsv.gz")

        ls = self.work_dir / Path("denovo.gsv.gz")
        lgs = self.work_dir / Path("gs.gsv.gz")
        lg = self.work_dir / Path("graph.gsv.gz")

        for cat in categories:
            filtered_vcfs = [ str(self.work_dir / Path(f"{caller}_{cat}_{opt.min_count}_filtered.vcf")) for caller in callers + ['sniffles2'] ]
            filtered_vcf_stats = [ str(self.work_dir / Path(f"{caller}_{cat}_{opt.min_count}_filtered.stat")) for caller in callers + ['sniffles2'] ]

            if cat == 'l+s':
                msvasm = ls
            if cat == 'l+g':
                msvasm = lg
            if cat == 'l+g+s':
                msvasm = lgs
            for (readid_tsv, vcf, filtered_vcf, outstat) in zip(self.readid_tsvs, self.som_vcfs, filtered_vcfs, filtered_vcf_stats):
                filtered_vcf = open(filtered_vcf, 'w')
                othercaller_filterasm(vcf, opt, readid_tsv, msvasm, outstat, consensus_sv_ids="", out_filtered_vcf=filtered_vcf)
                filtered_vcf.close()
        return

    def union_filtered_vcf(self, read_min_len, opt):
        from .union import union_sv

        opt.print_sv = True
        f = open(self.work_dir / Path(f"l_only_union.msv"), 'w')
        union_sv(self.som_vcfs[:3], read_min_len, opt, file_handler=f)
        f.close()
        opt.print_sv = False
        f = open(self.work_dir / Path(f"l_only_union_stat.msv"), 'w')
        union_sv(self.som_vcfs[:3], read_min_len, opt, file_handler=f)
        f.close()

        for cat in categories:
            filtered_vcfs = [ str(self.work_dir / Path(f"{caller}_{cat}_{opt.read_min_count}_filtered.vcf")) for caller in callers ]

            opt.print_sv = True
            f = open(self.work_dir / Path(f"{cat}_union.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            opt.print_sv = False
            f = open(self.work_dir / Path(f"{cat}_union_stat.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            f = open(self.work_dir / Path(f"{cat}_union_dedup.msv"), 'w')
            insilico_truth(str(self.work_dir / Path(f"{cat}_union.msv")), f)
            f.close()


    def output_consensus(self, read_min_len, opt):
        """Export the 3-caller raw consensus and the consensus subsets dropped
        by the self-assembly (l+s) and pangenome (l+g) filters.

        Consumes self.som_vcfs[:3] (severus/savana/nanomonsv) and the
        {cat}_union.msv files written by self.union_filtered_vcf. Must be
        called AFTER union_filtered_vcf.
        """
        from .union import union_sv

        # 2a. Raw consensus: all-three-caller groups, then dedup to one rep per group.
        consensus_union = self.work_dir / "consensus_union.msv"
        opt.print_sv = True
        with open(consensus_union, "w") as fh:
            union_sv(self.som_vcfs[:3], read_min_len, opt,
                     file_handler=fh, min_file_count=3)

        consensus_dedup = self.work_dir / "consensus_3caller_dedup.msv"
        with open(consensus_dedup, "w") as fh:
            insilico_truth(str(consensus_union), fh)

        # 2b. Per-filter lost subsets: consensus reps overlapping zero rows in
        # the filter's union. l+s_union.msv and l+g_union.msv already exist on
        # disk from union_filtered_vcf. Scope is intentionally l+s and l+g
        # only (NOT l+g+s) per the locked design.
        for cat in ["l+s", "l+g"]:
            target = self.work_dir / f"{cat}_union.msv"
            out_path = self.work_dir / f"consensus_lost_by_{cat}.msv"
            with open(out_path, "w") as fh:
                _consensus_lost_by_filter(str(consensus_dedup), str(target),
                                          opt, file_handler=fh)


    def save_timings(self, tsv_path=None):
        """Save collected timings to a TSV file"""
        if not self.timings:
            print("[WARNING] No timing data collected.")
            return

        if tsv_path is None:
            tsv_path = self.work_dir / "minisv_timings.tsv"

        fieldnames = ["step", "time_sec", "max_memory_mb"]

        with open(tsv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for record in self.timings:
                writer.writerow(record)

        total_time = sum(r["time_sec"] for r in self.timings)
        overall_max_mem = max(r["max_memory_mb"] for r in self.timings)

        print(f"\n[SUMMARY] Timings saved to: {tsv_path}")
        print(f"   Total time: {total_time:.2f} seconds")
        print(f"   Highest memory used: {overall_max_mem:.1f} MB")


class MSVTestCutOff:
    def __init__(self, som_vcfs, readid_tsvs, work_dir="minisv_work", new_work_dir="new_work_dir", filtered_readcount_cutoff=2):
        self.som_vcfs = som_vcfs
        self.readid_tsvs = readid_tsvs

        self.work_dir = Path(work_dir)
        self.new_work_dir = Path(new_work_dir)
        self.new_work_dir.mkdir(parents=True, exist_ok=True)

        self.filtered_readcount_cutoff = filtered_readcount_cutoff
        self.timings = []

    def extract_read_ids(self, bed, min_len, min_read_len, min_count, ignore_flt, check_gt):
        """ only extract somatic SV read ids """
        with Timer("extract_read_ids", self.timings):
            if bed is not None and not isinstance(bed, dict):
                bed = gc_read_bed(bed)

            self.read_ids_file = self.work_dir / "readid.names"
            som_read_ids = set()

            self.read_id_dict = {}
            for (readid_tsv, vcf) in zip(self.readid_tsvs, self.som_vcfs):
                assert os.path.exists(readid_tsv), f"{readid_tsv} not exists"
                assert os.path.exists(vcf), f"{vcf} not exists"

                self.read_id_dict |= parse_readids(readid_tsv)
                vcf = gc_parse_sv(vcf, min_read_len, min_count, ignore_flt, check_gt)

                for i in range(len(vcf)):
                    t = vcf[i]
                    # Not long enough for non-BND type
                    # NOTE: here use 100bp as the default length instead of 80
                    if t.SVTYPE != "BND" and abs(t.SVLEN) < min_len:
                        continue

                    if t.SVTYPE == "BND" and t.ctg == t.ctg2 and abs(t.SVLEN) < min_len:
                        continue

                    if t.count > 0 and t.count < min_count:
                        continue

                    if bed is not None:
                        if t.ctg not in bed or t.ctg2 not in bed:
                            continue
                        if len(iit_overlap(bed[t.ctg], t.pos, t.pos + 1)) == 0:
                            continue
                        if len(iit_overlap(bed[t.ctg2], t.pos2, t.pos2 + 1)) == 0:
                            continue

                    query_svid = str(parse_svid(t.svid))
                    assert query_svid in self.read_id_dict, 'somatic sv not in read id table'
                    som_read_ids.add(query_svid)

            read_names = set()
            for s in list(som_read_ids):
                read_names |= set(self.read_id_dict[s])
       
            with open(self.read_ids_file, "w") as fin:
                for i in sorted(list(read_names)):
                    fin.write(f"{i}\n")
   
    def extract_reads(self):
        """samtools fastq aln.bam | seqtk subseq - read-names.txt > reads.fq"""
        with Timer("extract_reads", self.timings):
            self.fastq_out = self.work_dir / "som_reads.fq.gz"
            if os.path.exists(self.fastq_out):
                return

            cmd = f"samtools fastq {self.bam_path} | seqtk subseq - {self.read_ids_file} | gzip  > {self.fastq_out}"
            subprocess.run(cmd, shell=True, check=True)
            print(f"Reads extracted to {self.fastq_out}")
            return

    def build_pooled_reference(self, out_fa="pooled_reference.fa"):
        with Timer("build_pooled_reference", self.timings):
            self.pooled_ref = self.work_dir / out_fa
            with open(self.pooled_ref, "w") as out:
                for ref in [self.hap1_denovo_ref_path, self.hap2_denovo_ref_path]:
                    with gzip.open(ref, "rt") as f:
                        out.write(f.read())
            return mp.Aligner(str(self.pooled_ref)) # preset??

    def align_reads_to_self(self, paf='denovo_aligned.paf.gz'):
        with Timer("align_reads_to_denovo", self.timings):
            self.paf_out = self.work_dir / paf
            if os.path.exists(self.paf_out):
                return

            if self.platform == 'hifi':
                cmd = f"{self.mm2} --ds -t 4 -cx map-hifi -s50 <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            else:
                cmd = f"{self.mm2} --ds -t 4 -cx lr:hq <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            #cmd = f"{self.mm2} --ds -t 4 -cx lr:hq <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to self assembly")

    def build_reference(self):
        return mp.Aligner(str(self.ref))

    def align_reads_to_grch38(self, paf='grch38_aligned.paf.gz'):
        with Timer("align_reads_to_grch38", self.timings):
            self.grch38_paf_out = self.work_dir / paf
            if os.path.exists(self.grch38_paf_out):
                return

            if self.platform == 'hifi':
                cmd = f"{self.mm2} --ds -t 4 -cx map-hifi -s50 {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            else:
                cmd = f"{self.mm2} --ds -t 4 -cx lr:hq {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to grch38")

    def align_reads_to_graph(self, gaf='denovo_aligned.gaf.gz'):
        with Timer("align_reads_to_graph", self.timings):
            self.gaf_out = self.work_dir / gaf
            if os.path.exists(self.gaf_out):
                return
            cmd = f"{self.mg} -cxlr -t 4 {self.graph_ref_path} {self.fastq_out} | gzip - > {self.gaf_out}"
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"Reads aligned to {self.graph_ref_path} to generate {self.gaf_out}")
            return

    def parse_raw_sv_grch38(self, opt, gsv='grch38.gsv.gz',
                            filtered_grch38_gsv='grch38_filtered.gsv.gz'):
        with Timer("parse_raw_sv_grch38", self.timings):
            self.grch38_gsv_out = self.work_dir / gsv
            self.filtered_grch38_gsv_out = self.work_dir / filtered_grch38_gsv
            print("**********************")
            print(self.filtered_grch38_gsv_out)
            f = gzip.open(self.grch38_gsv_out, 'wt')
            load_reads(str(self.grch38_paf_out), opt, f)
            f.close()

            ##opt.bed = gc_read_bed(opt.bed) if b is not None else None
            print("**************")
            print(opt.cen)
            print(opt.bed)
            print("**************")

            if opt.bed is not None and not isinstance(opt.bed, dict):
                opt.bed = gc_read_bed(opt.bed)

            if opt.bed is not None or opt.cen is not None:
                f = gzip.open(self.filtered_grch38_gsv_out, 'wt')
                with gzip.open(self.grch38_gsv_out, "rt") as fin:
                    for line in fin:
                        if opt.cen is not None:
                            t = line.strip().split("\t")
                            v = parse_sv(t)

                            #NOTE: test if we do not use centromere filtering
                            ## Step 1. post filter grch38-based sv signals by centromere
                            if v.cen_overlap is not None and v.cen_overlap > 0:
                                continue
                            if v.cen_dist is not None and v.cen_dist <= 5e5:
                                continue

                            ### Step 2. filter by overlapping with high conf region
                            #if opt.bed is not None:
                            #    if v.ctg not in opt.bed or v.ctg2 not in opt.bed:
                            #        continue
                            #    if len(iit_overlap(opt.bed[v.ctg], v.pos, v.pos + 1)) == 0:
                            #        continue
                            #    if len(iit_overlap(opt.bed[v.ctg2], v.pos2, v.pos2 + 1)) == 0:
                            #        continue
                            # TODO Step 3: filtered by SV typing consistency between minisv/severus/nanomonsv...
                            print(line.strip(), file=f)
                f.close()
            return

    def parse_raw_sv_self(self, opt, gsv='denovo.gsv.gz'):
        with Timer("parse_raw_sv_denovo", self.timings):
            self.denovo_gsv_out = self.work_dir / gsv

            f = gzip.open(self.denovo_gsv_out, 'wt')
            load_reads(str(self.paf_out), opt, f)
            f.close()
            return

    def parse_raw_sv_graph(self, opt, gsv='graph.gsv.gz'):
        with Timer("parse_raw_sv_graph", self.timings):
            self.graph_gsv_out = self.work_dir / gsv

            f = gzip.open(self.graph_gsv_out, 'wt')
            load_reads(str(self.gaf_out), opt, f)
            f.close()
            return

    def isec_g(self, opt, w, msv='l+g.gsv.gz'):
        # l + g
        with Timer("isec_graph", self.timings):
            self.isec_graph_out = self.work_dir / msv

            f = gzip.open(self.isec_graph_out, 'wt')
            if opt.bed is not None or opt.cen is not None:
                isec(w, [self.filtered_grch38_gsv_out, self.graph_gsv_out], f)
            else:
                isec(w, [self.grch38_gsv_out, self.graph_gsv_out], f)
            f.close()
            return

    def isec_s(self, opt, w, msv='l+s.gsv.gz'):
        # l+s
        with Timer("isec_graph", self.timings):
            self.isec_denovo_out = self.work_dir / msv

            f = gzip.open(self.isec_denovo_out, 'wt')
            if opt.bed is not None or opt.cen is not None:
                isec(w, [self.filtered_grch38_gsv_out, self.denovo_gsv_out], f)
            else:
                isec(w, [self.grch38_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def isec_gs(self, opt, w, msv='gs.gsv.gz'):
        # l+s
        with Timer("isec_gs", self.timings):
            self.isec_gs_out = self.work_dir / msv

            f = gzip.open(self.isec_gs_out, 'wt')
            #if opt.bed is not None or opt.cen is not None:
            #    isec(w, [self.filtered_grch38_gsv_out, self.graph_gsv_out, self.denovo_gsv_out], f)
            #else:
            #    isec(w, [self.grch38_gsv_out, self.graph_gsv_out, self.denovo_gsv_out], f)
            isec(w, [self.graph_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def parse_ids_from_gsv(self, msv):
        with gzip.open(self.work_dir / msv, 'rt') as fin:
            for line in fin:
                t = line.strip().split('\t')

                is_bp = False
                if re.match(r"[><]", t[2]):
                    is_bp = True

                col_info = 8 if is_bp else 6
                name = t[col_info-3]
                yield name

    def othercaller_filterasm(self, opt):
        """ the same interface as the former filterasm 
        """
        ls = self.work_dir / Path("denovo.gsv.gz")
        lgs = self.work_dir / Path("gs.gsv.gz")
        lg = self.work_dir / Path("graph.gsv.gz")

        for cat in categories:
            filtered_vcfs = [ str(self.new_work_dir / Path(f"{caller}_{cat}_{opt.min_count}_filtered.vcf")) for caller in callers + ['sniffles2'] ]
            filtered_vcf_stats = [ str(self.new_work_dir / Path(f"{caller}_{cat}_{opt.min_count}_filtered.stat")) for caller in callers + ['sniffles2'] ]

            if cat == 'l+s':
                msvasm = ls
            if cat == 'l+g':
                msvasm = lg
            if cat == 'l+g+s':
                msvasm = lgs
            for (readid_tsv, vcf, filtered_vcf, outstat) in zip(self.readid_tsvs, self.som_vcfs, filtered_vcfs, filtered_vcf_stats):
                filtered_vcf = open(filtered_vcf, 'w')
                othercaller_filterasm(vcf, opt, readid_tsv, msvasm, outstat, consensus_sv_ids="", out_filtered_vcf=filtered_vcf)
                filtered_vcf.close()
        return

    def union_filtered_vcf(self, read_min_len, opt):
        from .union import union_sv

        opt.print_sv = True
        f = open(self.new_work_dir / Path(f"l_only_union.msv"), 'w')
        union_sv(self.som_vcfs[:3], read_min_len, opt, file_handler=f)
        f.close()
        opt.print_sv = False
        f = open(self.new_work_dir / Path(f"l_only_union_stat.msv"), 'w')
        union_sv(self.som_vcfs[:3], read_min_len, opt, file_handler=f)
        f.close()

        for cat in categories:
            filtered_vcfs = [ str(self.new_work_dir / Path(f"{caller}_{cat}_{opt.read_min_count}_filtered.vcf")) for caller in callers ]

            opt.print_sv = True
            f = open(self.new_work_dir / Path(f"{cat}_union.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            opt.print_sv = False
            f = open(self.new_work_dir / Path(f"{cat}_union_stat.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            f = open(self.new_work_dir / Path(f"{cat}_union_dedup.msv"), 'w')
            insilico_truth(str(self.new_work_dir / Path(f"{cat}_union.msv")), f)
            f.close()


    def save_timings(self, tsv_path=None):
        """Save collected timings to a TSV file"""
        if not self.timings:
            print("[WARNING] No timing data collected.")
            return

        if tsv_path is None:
            tsv_path = self.new_work_dir / "minisv_timings.tsv"

        fieldnames = ["step", "time_sec", "max_memory_mb"]

        with open(tsv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for record in self.timings:
                writer.writerow(record)

        total_time = sum(r["time_sec"] for r in self.timings)
        overall_max_mem = max(r["max_memory_mb"] for r in self.timings)

        print(f"\n[SUMMARY] Timings saved to: {tsv_path}")
        print(f"   Total time: {total_time:.2f} seconds")
        print(f"   Highest memory used: {overall_max_mem:.1f} MB")


class MinisvReadsTrio:
    def __init__(self, som_vcfs, readid_tsvs, bam_path, ref, dad_hap1_denovo_ref_path, dad_hap2_denovo_ref_path, mom_hap1_denovo_ref_path, mom_hap2_denovo_ref_path, 
                 work_dir="minisv_work", filtered_readcount_cutoff=2, platform='hifi', mm2="minimap2", mg="minigraph"):
        self.som_vcfs = [som_vcfs]
        self.readid_tsvs = [readid_tsvs]
        self.mm2 = mm2
        self.mg = mg

        self.bam_path = Path(bam_path)

        self.ref = Path(ref)
        self.dad_hap1_denovo_ref_path = Path(dad_hap1_denovo_ref_path)
        self.dad_hap2_denovo_ref_path = Path(dad_hap2_denovo_ref_path)

        self.mom_hap1_denovo_ref_path = Path(mom_hap1_denovo_ref_path)
        self.mom_hap2_denovo_ref_path = Path(mom_hap2_denovo_ref_path)

        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.filtered_readcount_cutoff = filtered_readcount_cutoff
        self.platform = platform
        self.timings = []

    def extract_read_ids(self, bed, min_len, min_read_len, min_count, ignore_flt, check_gt):
        """ only extract somatic SV read ids """
        with Timer("extract_read_ids", self.timings):
            if bed is not None and not isinstance(bed, dict):
                bed = gc_read_bed(bed)

            self.read_ids_file = self.work_dir / "readid.names"
            som_read_ids = set()

            self.read_id_dict = {}
            for (readid_tsv, vcf) in zip(self.readid_tsvs, self.som_vcfs):
                assert os.path.exists(readid_tsv), f"{readid_tsv} not exists"
                assert os.path.exists(vcf), f"{vcf} not exists"

                self.read_id_dict |= parse_readids(readid_tsv)
                vcf = gc_parse_sv(vcf, min_read_len, min_count, ignore_flt, check_gt)

                for i in range(len(vcf)):
                    t = vcf[i]
                    # Not long enough for non-BND type
                    # NOTE: here use 100bp as the default length instead of 80
                    if t.SVTYPE != "BND" and abs(t.SVLEN) < min_len:
                        continue

                    if t.SVTYPE == "BND" and t.ctg == t.ctg2 and abs(t.SVLEN) < min_len:
                        continue

                    if t.count > 0 and t.count < min_count:
                        continue

                    if bed is not None:
                        if t.ctg not in bed or t.ctg2 not in bed:
                            continue
                        if len(iit_overlap(bed[t.ctg], t.pos, t.pos + 1)) == 0:
                            continue
                        if len(iit_overlap(bed[t.ctg2], t.pos2, t.pos2 + 1)) == 0:
                            continue

                    query_svid = str(parse_svid(t.svid))
                    assert query_svid in self.read_id_dict, 'somatic sv not in read id table'
                    som_read_ids.add(query_svid)

            read_names = set()
            for s in list(som_read_ids):
                read_names |= set(self.read_id_dict[s])
       
            with open(self.read_ids_file, "w") as fin:
                for i in sorted(list(read_names)):
                    fin.write(f"{i}\n")
   
    def extract_reads(self):
        """samtools fastq aln.bam | seqtk subseq - read-names.txt > reads.fq"""
        with Timer("extract_reads", self.timings):
            self.fastq_out = self.work_dir / "som_reads.fq.gz"
            if os.path.exists(self.fastq_out):
                return

            cmd = f"samtools fastq {self.bam_path} | seqtk subseq - {self.read_ids_file} | gzip  > {self.fastq_out}"
            subprocess.run(cmd, shell=True, check=True)
            print(f"Reads extracted to {self.fastq_out}")
            return

    def build_pooled_reference(self, out_fa="pooled_reference.fa"):
        with Timer("build_pooled_reference", self.timings):
            self.pooled_ref = self.work_dir / out_fa
            with open(self.pooled_ref, "w") as out:
                for ref in [self.mom_hap1_denovo_ref_path, self.mom_hap2_denovo_ref_path, self.mom_hap1_denovo_ref_path, self.mom_hap2_denovo_ref_path]:
                    with gzip.open(ref, "rt") as f:
                        out.write(f.read())
            return mp.Aligner(str(self.pooled_ref)) # preset??

    def align_reads_to_self(self, paf='denovo_aligned.paf.gz'):
        with Timer("align_reads_to_denovo", self.timings):
            self.paf_out = self.work_dir / paf
            if os.path.exists(self.paf_out):
                return

            if self.platform == 'hifi':
                cmd = f"{self.mm2} --ds -t 4 -cx map-hifi -s50 <(cat {self.dad_hap1_denovo_ref_path} {self.dad_hap2_denovo_ref_path} {self.mom_hap1_denovo_ref_path} {self.mom_hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            else:
                cmd = f"{self.mm2} --ds -t 4 -cx lr:hq <(cat {self.dad_hap1_denovo_ref_path} {self.dad_hap2_denovo_ref_path} {self.mom_hap1_denovo_ref_path} {self.mom_hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            #cmd = f"{self.mm2} --ds -t 4 -cx lr:hq <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to self assembly")

    def build_reference(self):
        return mp.Aligner(str(self.ref))

    def align_reads_to_grch38(self, paf='grch38_aligned.paf.gz'):
        with Timer("align_reads_to_grch38", self.timings):
            self.grch38_paf_out = self.work_dir / paf
            if os.path.exists(self.grch38_paf_out):
                return

            if self.platform == 'hifi':
                cmd = f"{self.mm2} --ds -t 4 -cx map-hifi -s50 {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            else:
                cmd = f"{self.mm2} --ds -t 4 -cx lr:hq {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to grch38")

    def align_reads_to_graph(self, gaf='denovo_aligned.gaf.gz'):
        with Timer("align_reads_to_graph", self.timings):
            self.gaf_out = self.work_dir / gaf
            if os.path.exists(self.gaf_out):
                return
            cmd = f"{self.mg} -cxlr -t 4 {self.graph_ref_path} {self.fastq_out} | gzip - > {self.gaf_out}"
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"Reads aligned to {self.graph_ref_path} to generate {self.gaf_out}")
            return

    def parse_raw_sv_grch38(self, opt, gsv='grch38.gsv.gz',
                            filtered_grch38_gsv='grch38_filtered.gsv.gz'):
        with Timer("parse_raw_sv_grch38", self.timings):
            self.grch38_gsv_out = self.work_dir / gsv
            self.filtered_grch38_gsv_out = self.work_dir / filtered_grch38_gsv
            print("**********************")
            print(self.filtered_grch38_gsv_out)
            f = gzip.open(self.grch38_gsv_out, 'wt')
            load_reads(str(self.grch38_paf_out), opt, f)
            f.close()

            ##opt.bed = gc_read_bed(opt.bed) if b is not None else None
            print("**************")
            print(opt.cen)
            print(opt.bed)
            print("**************")

            if opt.bed is not None and not isinstance(opt.bed, dict):
                opt.bed = gc_read_bed(opt.bed)

            if opt.bed is not None or opt.cen is not None:
                f = gzip.open(self.filtered_grch38_gsv_out, 'wt')
                with gzip.open(self.grch38_gsv_out, "rt") as fin:
                    for line in fin:
                        if opt.cen is not None:
                            t = line.strip().split("\t")
                            v = parse_sv(t)

                            #NOTE: test if we do not use centromere filtering
                            ## Step 1. post filter grch38-based sv signals by centromere
                            if v.cen_overlap is not None and v.cen_overlap > 0:
                                continue
                            if v.cen_dist is not None and v.cen_dist <= 5e5:
                                continue

                            ### Step 2. filter by overlapping with high conf region
                            #if opt.bed is not None:
                            #    if v.ctg not in opt.bed or v.ctg2 not in opt.bed:
                            #        continue
                            #    if len(iit_overlap(opt.bed[v.ctg], v.pos, v.pos + 1)) == 0:
                            #        continue
                            #    if len(iit_overlap(opt.bed[v.ctg2], v.pos2, v.pos2 + 1)) == 0:
                            #        continue
                            # TODO Step 3: filtered by SV typing consistency between minisv/severus/nanomonsv...
                            print(line.strip(), file=f)
                f.close()
            return

    def parse_raw_sv_self(self, opt, gsv='denovo.gsv.gz'):
        with Timer("parse_raw_sv_denovo", self.timings):
            self.denovo_gsv_out = self.work_dir / gsv

            f = gzip.open(self.denovo_gsv_out, 'wt')
            load_reads(str(self.paf_out), opt, f)
            f.close()
            return

    def parse_raw_sv_graph(self, opt, gsv='graph.gsv.gz'):
        with Timer("parse_raw_sv_graph", self.timings):
            self.graph_gsv_out = self.work_dir / gsv

            f = gzip.open(self.graph_gsv_out, 'wt')
            load_reads(str(self.gaf_out), opt, f)
            f.close()
            return

    def isec_g(self, opt, w, msv='l+g.gsv.gz'):
        # l + g
        with Timer("isec_graph", self.timings):
            self.isec_graph_out = self.work_dir / msv

            f = gzip.open(self.isec_graph_out, 'wt')
            if opt.bed is not None or opt.cen is not None:
                isec(w, [self.filtered_grch38_gsv_out, self.graph_gsv_out], f)
            else:
                isec(w, [self.grch38_gsv_out, self.graph_gsv_out], f)
            f.close()
            return

    def isec_s(self, opt, w, msv='l+s.gsv.gz'):
        # l+s
        with Timer("isec_graph", self.timings):
            self.isec_denovo_out = self.work_dir / msv

            f = gzip.open(self.isec_denovo_out, 'wt')
            if opt.bed is not None or opt.cen is not None:
                isec(w, [self.filtered_grch38_gsv_out, self.denovo_gsv_out], f)
            else:
                isec(w, [self.grch38_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def isec_gs(self, opt, w, msv='gs.gsv.gz'):
        # l+s
        with Timer("isec_gs", self.timings):
            self.isec_gs_out = self.work_dir / msv

            f = gzip.open(self.isec_gs_out, 'wt')
            #if opt.bed is not None or opt.cen is not None:
            #    isec(w, [self.filtered_grch38_gsv_out, self.graph_gsv_out, self.denovo_gsv_out], f)
            #else:
            #    isec(w, [self.grch38_gsv_out, self.graph_gsv_out, self.denovo_gsv_out], f)
            isec(w, [self.graph_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def parse_ids_from_gsv(self, msv):
        with gzip.open(self.work_dir / msv, 'rt') as fin:
            for line in fin:
                t = line.strip().split('\t')

                is_bp = False
                if re.match(r"[><]", t[2]):
                    is_bp = True

                col_info = 8 if is_bp else 6
                name = t[col_info-3]
                yield name

    def othercaller_filterasm(self, opt):
        """ the same interface as the former filterasm 
        """

        # filter each caller separately instead of jointly 
        # joint filtering is error-prone
        ls = self.work_dir / Path("denovo.gsv.gz")

        for cat in ['l+s']:
            filtered_vcfs = [ str(self.work_dir / Path(f"{caller}_{cat}_{opt.min_count}_filtered.vcf")) for caller in ['sniffles2'] ]
            filtered_vcf_stats = [ str(self.work_dir / Path(f"{caller}_{cat}_{opt.min_count}_filtered.stat")) for caller in ['sniffles2'] ]

            if cat == 'l+s':
                msvasm = ls
            if cat == 'l+g':
                msvasm = lg
            if cat == 'l+g+s':
                msvasm = lgs
            for (readid_tsv, vcf, filtered_vcf, outstat) in zip(self.readid_tsvs, self.som_vcfs, filtered_vcfs, filtered_vcf_stats):
                filtered_vcf = open(filtered_vcf, 'w')
                othercaller_filterasm(vcf, opt, readid_tsv, msvasm, outstat, consensus_sv_ids="", out_filtered_vcf=filtered_vcf)
                filtered_vcf.close()
        return

    def union_filtered_vcf(self, read_min_len, opt):
        from .union import union_sv

        opt.print_sv = True
        f = open(self.work_dir / Path(f"l_only_union.msv"), 'w')
        union_sv(self.som_vcfs[:3], read_min_len, opt, file_handler=f)
        f.close()
        opt.print_sv = False
        f = open(self.work_dir / Path(f"l_only_union_stat.msv"), 'w')
        union_sv(self.som_vcfs[:3], read_min_len, opt, file_handler=f)
        f.close()

        for cat in categories:
            filtered_vcfs = [ str(self.work_dir / Path(f"{caller}_{cat}_{opt.read_min_count}_filtered.vcf")) for caller in callers ]

            opt.print_sv = True
            f = open(self.work_dir / Path(f"{cat}_union.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            opt.print_sv = False
            f = open(self.work_dir / Path(f"{cat}_union_stat.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            f = open(self.work_dir / Path(f"{cat}_union_dedup.msv"), 'w')
            insilico_truth(str(self.work_dir / Path(f"{cat}_union.msv")), f)
            f.close()


    def save_timings(self, tsv_path=None):
        """Save collected timings to a TSV file"""
        if not self.timings:
            print("[WARNING] No timing data collected.")
            return

        if tsv_path is None:
            tsv_path = self.work_dir / "minisv_timings.tsv"

        fieldnames = ["step", "time_sec", "max_memory_mb"]

        with open(tsv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for record in self.timings:
                writer.writerow(record)

        total_time = sum(r["time_sec"] for r in self.timings)
        overall_max_mem = max(r["max_memory_mb"] for r in self.timings)

        print(f"\n[SUMMARY] Timings saved to: {tsv_path}")
        print(f"   Total time: {total_time:.2f} seconds")
        print(f"   Highest memory used: {overall_max_mem:.1f} MB")

## new SV filter from bam -> fastq reads -> gsv
