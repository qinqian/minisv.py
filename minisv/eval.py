import math
from collections import defaultdict
import re
from dataclasses import dataclass
from typing import Optional

from .merge import svinfo
from .regex import re_info
from operator import itemgetter, attrgetter
from .util import is_gzipped
import gzip



SIZE_BINS = [
    (50, 100, "50-100"),
    (100, 500, "100-500"),
    (500, 1000, "500-1000"),
    (1000, 10000, "1k-10k"),
    (10000, 20000, "10k-20k"),
    (20000, 100000, "20k-100k"),
    (100000, float('inf'), ">100k"),
]



def eval(inputfiles, opt, size):
    """Evaluate the performance of the test SVs against the base SVs
    Args:
        base: list of base SVs
        test: list of test SVs
        opt: EvalConfig object
    Returns:
        None
    """

    min_read_len = math.floor(opt.min_len * opt.read_len_ratio + 0.499)

    if len(inputfiles) == 2:
        # base: truthset
        base = gc_parse_sv(inputfiles[0], min_read_len, opt.min_count, opt.ignore_flt, opt.check_gt)
        # compared vcf
        test = gc_parse_sv(inputfiles[1], min_read_len, opt.min_count, opt.ignore_flt, opt.check_gt)
        if not size:
            tot_fn, fn = gc_cmp_sv(opt, test, base, "FN")
            tot_fp, fp = gc_cmp_sv(opt, base, test, "FP")
        else:
            tot_fn, fn = gc_cmp_sv_wt_sizes(opt, test, base, "FN")
            tot_fp, fp = gc_cmp_sv_wt_sizes(opt, base, test, "FP")
        print("RN", tot_fn, fn, round(fn / tot_fn, 4), inputfiles[0], sep="\t")
        print("RP", tot_fp, fp, round(fp / tot_fp, 4), inputfiles[1], sep="\t")
        print("\n# SV size bin breakdown")
        print("SizeBin\tFN_tot\tFN\tRecall\tFP_tot\tFP\tPrecision")
    
        if size:
            for _, _, bin_name in SIZE_BINS:
                fn_stats = opt.size_bin_stats.get("FN", {})
                fp_stats = opt.size_bin_stats.get("FP", {})
    
                fn_t = fn_stats['tot'].get(bin_name, 0)
                fn_e = fn_stats['err'].get(bin_name, 0)
                fp_t = fp_stats['tot'].get(bin_name, 0)
                fp_e = fp_stats['err'].get(bin_name, 0)
    
                recall = 1 - (fn_e / fn_t) if fn_t > 0 else 1.0
                prec   = 1 - (fp_e / fp_t) if fp_t > 0 else 1.0
    
                print(f"{bin_name}\t{fn_t}\t{fn_e}\t{recall:.4f}\t{fp_t}\t{fp_e}\t{prec:.4f}")
    else:
        # multi-sample mode
        vcf = []
        for i in range(len(inputfiles)):
            vcf.append(
                gc_parse_sv(inputfiles[i], min_read_len, opt.min_count, opt.ignore_flt, opt.check_gt)
            )

        if opt.merge:
            for i in range(len(inputfiles)):
                other = []
                for j in range(len(inputfiles)):
                    if i != j:
                        other.append(vcf[j])

                merge = gc_eval_merge_sv(opt.win_size, opt.min_len_ratio, other)
                tot_fp, fp = gc_cmp_sv(opt, merge, vcf[i], 'merge')
                merge2 = []
                for k in range(len(merge)):
                    if merge[k].merge >= 2:
                        merge2.append(merge[k])
                tot_fn, fn = gc_cmp_sv(opt, vcf[i], merge2, 'merge')

                fn_val = round(fn / tot_fn, 4)
                fn_val = f"{fn_val:.4f}"
                fp_val = round(fp / tot_fp, 4)
                fp_val = f"{fp_val:.4f}"
                print("RN", tot_fn, fn, fn_val, inputfiles[i], sep='\t')
                print("RP", tot_fp, fp, fp_val, inputfiles[i], sep='\t')
        else:
            for i in range(len(inputfiles)):
                a = ["SN"]
                for j in range(len(inputfiles)):
                    cnt, err = gc_cmp_sv(opt, vcf[i], vcf[j], "XX")
                    if i != j:
                        val = round(1 - err / cnt, 4)
                        a.append(f"{val:.4f}")
                    else:
                        a.append(str(cnt))
                print("\t".join(a), inputfiles[i], sep='\t')



def gc_eval_merge_sv(win_size, min_len_ratio, all):
    """merge SV set"""
    vcf = []

    for i in range(len(all)):
        for j in range(len(all[i])):
            vcf.append(all[i][j])

    vcf.sort(key=attrgetter('ctg', 'pos'))
    out = []

    for i in range(len(vcf)):
        merge_to = -1
        vcf[i].merge = 0
        j = len(out) - 1
        while j >= 0 and out[j].ctg == vcf[i].ctg and vcf[i].pos - out[j].pos <= win_size:
            if gc_cmp_same_sv1(win_size, min_len_ratio, out[j], vcf[i]):
                 merge_to = j
                 break
            j -= 1

        if merge_to < 0:
            merge_to = len(out)
            out.append(vcf[i])
        out[merge_to].merge += 1
    
    return out


def gc_get_count_vcf(t):
    """
    Parse read counts from VCF info field
    """
    for m in re_info.findall(t[7]):
        if m[0] == "SUPPORT": # Sniffles2
            return int(m[1])
        elif m[0] == "TUMOUR_SUPPORT": # SAVANA 1.0.5
            return int(m[1])
        elif m[0] == "TUMOUR_READ_SUPPORT": # SAVANA 1.2
            return int(m[1])

    # DV: Severus & Sniffles2 & Svision
    # VR: nanomonsv
    if len(t) >= 10 and re.search("DV|VR", t[8]):
        fmt = t[8].split(":")
        fmt_i = -1
        n_fmt = 0
        for i in range(len(fmt)):
            if fmt[i] in ["DV", "VR"]:
                fmt_i = i
                n_fmt += 1
        if n_fmt == 1 and fmt_i >= 0:
            cnt = 0
            for i in range(9, len(t)):
                cnt += int(t[i].split(':')[fmt_i])
            return cnt
    return -1


def gc_get_count_msv(count_info):
    s = count_info.split("|") # NOTE: maybe support multiple samples?
    cnt = [0, 0]

    for i in range(len(s)):
        m = re.findall(r"([^\s:]+):(\d+),(\d+)", s[i])
        if len(m) > 0:
            cnt[0] += int(m[0][1])
            cnt[1] += int(m[0][2])
    return cnt
    

# for all SV callers
def gc_parse_sv(file_path, min_read_len, min_count: int, ignore_flt: bool, check_gt: bool):
    sv = []
    ignore_id = {}
    if is_gzipped(file_path):
        f = gzip.open(file_path, 'rt')
    else:
        f = open(file_path)

    line_no = 0
    for line in f:
        line_no += 1
        # skip vcf header lines
        if line[0] == "#":
            continue

        t = line.strip().split("\t")
        # POS must be number
        if re.match(r"^\d+", t[1]) is None:
            continue

        # start pos
        t[1] = int(t[1])

        type = 0
        info = None
        inv  = False

        # different format use different columns for SV type
        if re.match(r"^[><][><]$", t[2]):
            col_info = 8
            # breakpoint type
            type = 3
            info = t[col_info]
        # VCF format
        elif len(t) >= 8 and re.search(r";", t[7]):
            # other caller VCF
            type = 1
            info = t[7]
        elif re.match(r"^\d+$", t[2]) and re.search(r";", t[6]):
            # INDEL type
            col_info = 6
            type = 2
            info = t[col_info]

        if type == 0:
            raise Exception("No type available...")
            continue

        svtype = None
        svlen = 0
        cnt_tot = -1
        
        m = re.findall(r"\bSVTYPE=([^\s;]+)", info)
        if len(m) > 0:
            svtype = m[0]

        if svtype == "INV":
            inv = True
            
        m = re.findall(r"\bSVLEN=([^\s;]+)", info)
        if len(m) > 0:
            svlen = int(m[0])

        # parse read count
        # MSV
        if type in [2, 3]:
            regex = re.compile(r"\bcount=([^\s;]+)")
            m = regex.findall(info)
            if len(m) > 0:
               cf, cr = gc_get_count_msv(m[0])
               cnt_tot = cf + cr
        # OTHER VCFs
        elif type == 1:
            cnt_tot = gc_get_count_vcf(t)

        # COLO829 truth set do not have count info
        # if cnt_tot <= 0:
        #     continue

        if cnt_tot > 0 and cnt_tot < min_count:
            continue

        # Parse sv coordinate
        # MSV BED-like line
        # for INDEL
        if type == 2:
            t[2] = int(t[2])
            if t[1] > t[2]:
                raise Exception("incorrect BED format")

            # NOTE: min_read_len is 80 instead of 100 here???
            if abs(svlen) < min_read_len:
                continue
            # print(svtype, type, svlen, cnt_tot, min_read_len)
            # raise Exception("Test...")
            sv.append(
                svinfo(
                    ctg=t[0],
                    pos=t[1],
                    ctg2=t[0],
                    pos2=t[2],
                    ori=">>",
                    SVTYPE=svtype,
                    SVLEN=svlen,
                    inv=inv,
                    count=cnt_tot,
                    vaf=1,  # MSV do not have VAF
                    svid=str(line_no)
                )
            )
        elif type == 3:
            # breakpoint line
            t[4] = int(t[4])

            if t[0] == t[3] and abs(svlen) < min_read_len:
                continue

            if t[0] == t[3] and (t[2] == "><" or t[2] == "<>"):
                inv = True

            # print(svtype, type, svlen, cnt_tot, min_read_len)
            # raise Exception("Test...")

            sv.append(
                svinfo(
                    ctg=t[0],
                    pos=t[1],
                    ctg2=t[3],
                    pos2=t[4],
                    ori=t[2],
                    SVTYPE=svtype,
                    SVLEN=svlen,
                    inv=inv,
                    count=cnt_tot,
                    vaf=1, # MSV do not have VAF
                    svid=str(line_no),
                )
            )
        elif type == 1:
            # VCF line
            # VCF filter
            # ignore filtered calls
            if (not ignore_flt) and t[6] != "PASS" and t[6] != ".":
                continue

            # not a variant
            if check_gt and len(t) >= 9 and re.match(r"^0[\/\|]0", t[9]):
                continue

            # reference allele
            rlen = len(t[3])
            en = t[1] + rlen - 1

            # initialize a sv info object
            s = svinfo(ctg=t[0], pos=t[1] - 1, ctg2=t[0], pos2=en, ori=">>", 
                       inv=inv, count=cnt_tot,
                       svid=t[2], vaf=1)

            # need to improve the VAF parser
            # parse_vaf()
            # m = re.findall(r"\dVAF=([^\s;]+)", info)                
            # if len(m) > 0:
            #     s.vaf = float(m[0])
            # raise Exception("test vcf")

            # adjust sv length by allele sequences for indel
            if re.match(r"^[A-Z,\*]+$", t[4]) and t[4] != "SV" and t[4] != "CSV":

                # assume full allele sequence;
                # override SVTYPE/SVLEN even if present
                # multiple alleles
                alt = t[4].split(",")
                for i in range(len(alt)):
                    a = alt[i]
                    # alt allele - reference allele
                    length = len(a) - rlen

                    if abs(length) < min_read_len:
                        continue

                    s.ori = ">>"
                    s.SVLEN = length
                    if length < 0:
                        s.SVTYPE = "DEL"
                    else:
                        s.SVTYPE = "INS"
                    sv.append(s)
            # other SV type encoding
            else:
                if t[2] != ".":
                    # ignore previously visited ID
                    if t[2] in ignore_id:
                        continue
                    ignore_id[t[2]] = 1

                # Severus/nanomonsv/Savana has MATEID
                # NOTE: \b is not correct,
                #       need to be ;
                m = re.findall(r"\b(MATE_ID|MATEID)=([^\s;]+)", info)
                if len(m) > 0:
                    ignore_id[m[0][1]] = 1

                if svtype is None:
                    # we don't infer SVTYPE from breakpoint for existing data
                    raise Exception(f"cannot determine SVTYPE {t}")

                s.SVTYPE = svtype
                # patch nanomonsv insertion size
                if svtype == "INS":
                    m = re.findall(r"SVINSLEN=([^\s;]+)", info)
                    if len(m) > 0:
                        svlen = int(m[0])

                if svtype != "BND" and abs(svlen) < min_read_len:
                    continue  # too short

                # correct Deletion SVLEN consistently
                # severus has DEL >0 length
                if svtype == "DEL" and svlen > 0:
                    svlen = -svlen
                s.SVLEN = svlen

                # NOTE: \b not right, use ;
                m = re.findall(r"\bEND=(\d+)", info)
                if len(m) > 0:
                    s.pos2 = int(m[0])
                elif rlen == 1:
                    # ignore one-sided breakpoint
                    # NOTE: alt allele at least > 7 characters
                    #       e.g., ]chr7:152222658]C
                    if svtype == "BND" and len(t[4]) < 6:
                        continue
                    # NOTE: correct breakpoint end position?
                    if svtype == "DEL" or svtype == "DUP" or svtype == "INV":
                        s.pos2 = s.pos + abs(svlen)

                # parse VCF alt allele
                # https://samtools.github.io/hts-specs/VCFv4.2.pdf
                # t[p[
                m = re.findall(r"^[A-Z]+\[([^\s:]+):(\d+)\[$", t[4])
                # t]p]
                m1 = re.findall(r"^\]([^\s:]+):(\d+)\][A-Z]+$", t[4])
                # [p[t: reverse comp piece extending right of p is joined before t
                m2 = re.findall(r"^\[([^\s:]+):(\d+)\[[A-Z]+$", t[4])
                # t]p]: reverse comp piece extending left of p is joined after t
                m3 = re.findall(r"^[A-Z]+\]([^\s:]+):(\d+)\]$", t[4])
                if len(m) > 0:
                    m = m[0]
                    s.ctg2 = m[0]
                    s.pos2 = int(m[1])
                    s.ori = ">>"
                elif len(m1) > 0:
                    m1 = m1[0]
                    s.ctg2 = m1[0]
                    s.pos2 = int(m1[1])
                    s.ori = "<<"
                elif len(m2) > 0:
                    m2 = m2[0]
                    s.ctg2 = m2[0]
                    s.pos2 = int(m2[1])
                    s.ori = "<>"
                elif len(m3) > 0:
                    m3 = m3[0]
                    s.ctg2 = m3[0]
                    s.pos2 = int(m3[1])
                    s.ori = "><"

                if s.ctg == s.ctg2 and (s.ori == "><" or s.ori == "<>"):
                    s.inv = True

                if s.inv and s.SVTYPE == "BND":
                    s.SVTYPE = "INV"
                    svlen = abs(s.pos2 - s.pos)
                    s.SVLEN = svlen

                if s.SVTYPE == "DUP" and s.ori == ">>":
                    s.ori = "<<"

                if s.ctg == s.ctg2 and s.ori != "><" and s.ori != "<>" and s.SVTYPE == "INV":
                    s.ori = "><"
                    s.inv = True

                if s.ori == ">>" and s.SVTYPE == "BND" and s.ctg == s.ctg2:
                    m = re.findall(r"DETAILED_TYPE=([^\s;]+)", info)
                    if len(m) > 0:
                        if m[0] == 'Templated_ins':
                            # NOTE: this is not tandem duplication
                            s.SVTYPE = "BND"
                    else:
                        s.SVTYPE = "DEL"
                elif s.ori == "<<" and s.SVTYPE == "BND" and s.ctg == s.ctg2:
                    s.SVTYPE = "DUP"

                if s.SVTYPE == "DEL" and s.SVLEN > 0:
                    svlen = -svlen
                    s.SVLEN = svlen

                if s.SVTYPE != "BND" and s.ctg != s.ctg2:
                    # not possible for indel, dup and inversion
                    raise Exception("different contigs for non-BND type")

                if (
                    s.SVTYPE == "BND"
                    and s.ctg == s.ctg2
                ):
                    if s.SVLEN == 0 and abs(s.pos2 - s.pos) < min_read_len:
                        continue
                    if s.SVLEN != 0 and abs(s.SVLEN) < min_read_len:
                        continue

                if s.ctg == s.ctg2 and s.pos > s.pos2:
                    tmp = s.pos
                    s.pos = s.pos2
                    s.pos2 = tmp
                sv.append(s)
    f.close()
    return sv


def gc_cmp_sv(opt, base, test, label):
    """Compare two lists of SVs
    Args:
        base: list of base SVs
        test: list of test SVs
        opt: EvalConfig object
    Returns:
        tuple of (total_fn, fn)
    """
    h = {}
    for i in range(len(base)):
        s = base[i]
        if s.ctg not in h:
            h[s.ctg] = []
        if s.ctg2 not in h:
            h[s.ctg2] = []
        h[s.ctg].append({"st": s.pos, "en": s.pos + 1, "data": s})
        h[s.ctg2].append({"st": s.pos2, "en": s.pos2 + 1, "data": s})

    for ctg in h:
        h[ctg] = iit_sort_copy(h[ctg])
        iit_index(h[ctg])

    tot = error = 0
    for j in range(len(test)):
        t = test[j]

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

        tot += 1
        n = eval1(opt, h, t.ctg, t.pos, t) + eval1(opt, h, t.ctg2, t.pos2, t)

        if n == 0:
            error += 1
            if opt.print_err:
                #print(
                #    label,
                #    t.ctg,
                #    t.pos,
                #    t.ori,
                #    t.ctg2,
                #    t.pos2,
                #    t.SVTYPE,
                #    t.svid,
                #    t.SVLEN,
                #    n,
                #    t.merge,
                #    t.count,
                #    sep="\t",
                #)
                print(
                    label,
                    t.ctg,
                    t.pos,
                    t.ori,
                    t.ctg2,
                    t.pos2,
                    t.SVTYPE,
                    t.svid,
                    t.SVLEN,
                    n,
                    sep="\t",
                )
        if opt.print_all:
            ##print(t.ctg, t.pos, t.ori, t.ctg2, t.pos2, t.SVTYPE, t.SVLEN, t.svid, n, t.merge, t.count, sep="\t")
            print(t.ctg, t.pos, t.ori, t.ctg2, t.pos2, t.SVTYPE, t.SVLEN, t.svid, n, sep="\t")

    return [tot, error]


def gc_cmp_sv_wt_sizes(opt, base, test, label):
    """Compare two lists of SVs
    Args:
        base: list of base SVs
        test: list of test SVs
        opt: EvalConfig object
    Returns:
        tuple of (total_fn, fn)
    """
    h = {}
    for i in range(len(base)):
        s = base[i]
        if s.ctg not in h:
            h[s.ctg] = []
        if s.ctg2 not in h:
            h[s.ctg2] = []
        h[s.ctg].append({"st": s.pos, "en": s.pos + 1, "data": s})
        h[s.ctg2].append({"st": s.pos2, "en": s.pos2 + 1, "data": s})

    for ctg in h:
        h[ctg] = iit_sort_copy(h[ctg])
        iit_index(h[ctg])

    tot = error = 0

    bin_tot = {name: 0 for _, _, name in SIZE_BINS}
    bin_err = {name: 0 for _, _, name in SIZE_BINS}
    if not hasattr(opt, f'size_bin_stats'):
        opt.size_bin_stats = {}
        opt.size_bin_stats[label] = {}
    else:
        opt.size_bin_stats[label] = {}

    opt.size_bin_stats[label]['tot'] = defaultdict(int)
    opt.size_bin_stats[label]['err'] = defaultdict(int)

    for j in range(len(test)):
        t = test[j]

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

        # === NEW: Assign to size bin ===
        bin_name = None
        for min_l, max_l, name in SIZE_BINS:
            if min_l <= abs(t.SVLEN) <= max_l:
                bin_name = name
                break

        # patch translocation
        if t.ctg != t.ctg2:
            bin_name = ">100k"

        if bin_name is not None:
            bin_tot[bin_name] += 1
            opt.size_bin_stats[label]['tot'][bin_name] += 1
        else:
            print('out of case')
            print(t)

        tot += 1

        n = eval1(opt, h, t.ctg, t.pos, t) + eval1(opt, h, t.ctg2, t.pos2, t)
        if n == 0:
            error += 1
            if bin_name is not None:
                bin_err[bin_name] += 1
                opt.size_bin_stats[label]['err'][bin_name] += 1
            if opt.print_err:
                print(
                    label,
                    t.ctg,
                    t.pos,
                    t.ori,
                    t.ctg2,
                    t.pos2,
                    t.SVTYPE,
                    t.svid,
                    t.SVLEN,
                    n,
                    t.merge,
                    t.count,
                    sep="\t",
                )
        if opt.print_all:
            print(t.ctg, t.pos, t.ori, t.ctg2, t.pos2, t.SVTYPE, t.SVLEN, t.svid, n, t.merge, t.count, sep="\t")

    return [tot, error]


@dataclass
class iit_obj:
    st: int
    en: int
    max: Optional[int] = 0
    data: Optional[svinfo] = None


# interval query
def iit_sort_copy(a):
    """Sort a list of SVs by start position and return a copy
    """
    a.sort(key=lambda x: x["st"])
    b = []
    for i in range(len(a)):
        b.append(iit_obj(st=a[i]["st"], en=a[i]["en"], max=0, data=a[i]["data"]))
    return b


def iit_index(a):
    """
    reference: https://github.com/attractivechaos/plb2/blob/master/src/python/bedcov.py
    """
    if len(a) == 0:
        return -1

    for i in range(0, len(a), 2):
        last = a[i].en
        a[i].max = a[i].en
        last_i = i

    # NOTE: what is this loop for?
    k = 1
    while 1 << k <= len(a):
        i0 = (1 << k) - 1
        step = 1 << (k + 1)
        x = 1 << (k - 1)
        for i in range(i0, len(a), step):
            # NOTE: why exchange position
            a[i].max = a[i].en
            if a[i].max < a[i - x].max:
                a[i].max = a[i - x].max
            e = a[i + x].max if i + x < len(a) else last
            if a[i].max < e:
                a[i].max = e
        last_i = last_i - x if (last_i >> k) & 1 > 0 else last_i + x
        if last_i < len(a):
            last = last if last > a[last_i].max else a[last_i].max
        k += 1
    return k - 1


def eval1(opt, h, ctg, pos, t):
    if ctg not in h:
        return False

    st = pos - opt.win_size if pos > opt.win_size else 0
    en = pos + opt.win_size
    if opt.dbg:
        print(ctg, st, en, opt.win_size, sep="\t")
    a = iit_overlap(h[ctg], st, en)
    n = 0
    for i in range(len(a)):
        if gc_cmp_same_sv1(opt.win_size, opt.min_len_ratio, a[i].data, t):
            n += 1
    return n


def iit_overlap(a, st, en):
    """ inexplicit interval tree """
    h = 0
    stack = []
    b = []

    while 1 << h <= len(a):
        h += 1

    h -= 1
    stack.append([(1 << h) - 1, h, 0])

    while len(stack) > 0:
        t = stack.pop()
        # NOTE: what does x, h, w represent?
        x = t[0]
        h = t[1]
        w = t[2]
        # NOTE: why this should be 3
        if h <= 3:
            # NOTE: why right shift then left shift bit
            i0 = x >> h << h
            i1 = i0 + (1 << (h + 1)) - 1
            if i1 >= len(a):
                i1 = len(a)
            i = i0
            while i < i1 and a[i].st < en:
                if st < a[i].en:
                    b.append(a[i])
                i += 1
        # NOTE: why this should be 0
        elif w == 0:
            stack.append([x, h, 1])
            y = x - (1 << (h - 1))
            if y >= len(a) or a[y].max > st:
                stack.append([y, h - 1, 0])
        elif x < len(a) and a[x].st < en:
            if st < a[x].en:
                b.append(a[x])
            stack.append([x + (1 << (h - 1)), h - 1, 0])
    return b


def gc_cmp_same_sv1(win_size, min_len_ratio, b, t):
    # print(win_size, min_len_ratio)
    # check type
    if b.SVTYPE != t.SVTYPE:  # type mismatch
        if (
            (not (b.SVTYPE == "DUP" and t.SVTYPE == "INS"))
            and (not (b.SVTYPE == "INS" and t.SVTYPE == "DUP"))
            and b.SVTYPE != "BND"  # NOTE: because some tools cannot type correctly?
            and t.SVTYPE != "BND"
        ):  # special case for INS vs DUP
            return False

    # check length
    # NOTE: length should be similar to each other
    len_check = (abs(abs(b.SVLEN) - abs(t.SVLEN)) <= 1000) or (
        abs(b.SVLEN) >= abs(t.SVLEN) * min_len_ratio
        and abs(t.SVLEN) >= abs(b.SVLEN) * min_len_ratio
    )

    if b.SVTYPE != "BND" and t.SVTYPE != "BND" and (not len_check):
        return False

    match1 = match2 = 0
    # check the coordinates of end points within window
    # NOTE: match1 and match2 denote two breakpoints
    if (
        t.ctg == b.ctg
        and t.pos >= b.pos - win_size
        and t.pos <= b.pos + win_size
    ):
        match1 |= 1

    if (
        t.ctg == b.ctg2
        and t.pos >= b.pos2 - win_size
        and t.pos <= b.pos2 + win_size
    ):
        match1 |= 2

    if (
        t.ctg2 == b.ctg
        and t.pos2 >= b.pos - win_size
        and t.pos2 <= b.pos + win_size
    ):
        match2 |= 1
    if (
        t.ctg2 == b.ctg2
        and t.pos2 >= b.pos2 - win_size
        and t.pos2 <= b.pos2 + win_size
    ):
        match2 |= 2

    # NOTE: duplicated condition? seems to be a bug
    if b.SVTYPE == "DUP" and t.SVTYPE == "INS":
        return (match1 & 1) != 0
    elif b.SVTYPE == "INS" and t.SVTYPE == "DUP":
        return (match1 & 1) != 0
    elif b.SVTYPE == "BND" or t.SVTYPE == "BND":
        return ((match1 & 1) != 0 and (match2 & 2) != 0) or (
            (match1 & 2) != 0 and (match2 & 1) != 0
        )
    else:
        return (match1 & 1) != 0 and (match2 & 2) != 0


def gc_read_bed(fn, is_gnomad=False, gnomad_af=0.01):
    if is_gzipped(fn):
        f = gzip.open(fn, 'rt')
    else:
        f = open(fn)

    if is_gnomad:
        h_pt, h_sv = {}, {}
        n_skip = 0
        for line in f:
            t = line.strip().split("\t")
            if t[0].startswith("#"):
                continue
            if len(t) < 49:
                n_skip += 1
                continue
            af = 0 if t[48].strip() == 'NA' else float(t[48])
            if af < gnomad_af:
                continue
            svtype = t[4]
            ctg2_raw = t[9]
            if ctg2_raw != 'NA' and ctg2_raw != t[0]:
                # cross-chrom: index on both ends for SV-identity matching.
                # SVLEN=0: cross-chrom rows are matched as BND regardless of SVTYPE
                # so gc_cmp_same_sv1 skips the length check (eval.py:834).
                if t[13] == 'NA':
                    continue
                pos2 = int(t[13])
                s = svinfo(
                    ctg=t[0], pos=int(t[1]),
                    ctg2=ctg2_raw, pos2=pos2,
                    ori=">>",
                    SVTYPE=svtype,
                    SVLEN=0,
                    inv=(svtype == 'INV'),
                    count=0, vaf=1,
                    svid=t[3],
                )
                h_sv.setdefault(t[0], []).append(
                    {"st": int(t[1]), "en": int(t[1]) + 1, "data": s})
                h_sv.setdefault(ctg2_raw, []).append(
                    {"st": pos2, "en": pos2 + 1, "data": s})
            else:
                # same-chrom (or CHR2 == NA): point-overlap tree (current behavior)
                h_pt.setdefault(t[0], []).append(
                    {"st": int(t[1]), "en": int(t[2]), "data": None})
        f.close()
        if n_skip:
            import sys as _sys

            print(
                f"[gc_read_bed] WARNING: {n_skip} gnomAD row(s) skipped — fewer than 49 columns",
                file=_sys.stderr,
            )
        for tree in (h_pt, h_sv):
            for ctg in tree:
                tree[ctg] = iit_sort_copy(tree[ctg])
                iit_index(tree[ctg])
        return h_pt, h_sv
    else:
        h = {}
        for line in f:
            t = line.strip().split("\t")
            if t[0].startswith("#"):
                continue
            if len(t) < 3:
                continue
            h.setdefault(t[0], []).append(
                {"st": int(t[1]), "en": int(t[2]), "data": None})
        f.close()
        for ctg in h:
            h[ctg] = iit_sort_copy(h[ctg])
            iit_index(h[ctg])
        return h
