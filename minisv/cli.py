#!/usr/bin/env python
import math

import re
import gzip
import sys
from dataclasses import dataclass
from typing import Optional

import rich_click as click

from .eval import gc_read_bed
from .phase import extract_phase_HP, annotate_HP

##from .cygafcall import GafParser as cy_GafParser
from .minisv import GafParser
from .filtercaller import call_filterseverus, call_filtersnf, call_filtermsv
from .io import gc_cmd_view, merge_indel_breakpoints, parseNum, write_vcf
from .ensemble import insilico_truth, double_strand_break
from .union import union_sv, advunion_sv, union_sv_with_tr

__version__ = "0.1.2"


@dataclass
class opt:
    cen: dict
    bed: Optional[str] = None
    min_mapq: int = 5
    min_mapq_end: int = 30
    min_frac: float = 0.7
    min_len: int = 100
    min_aln_len_end: int = 2000
    min_aln_len_mid: int = 50
    max_cnt_10k: int = 5
    dbg: bool = False
    polyA_pen: int = 5
    polyA_drop: int = 100
    name: str = "foo"


@dataclass
class mergeopt:
    min_cnt: int = 3
    min_cnt_strand: int = 2
    min_cnt_rt: int = 1
    min_rt_len: int = 10
    win_size: int = 100
    max_diff: float = 0.05
    min_cen_dist: int = 5e5
    max_allele: int = 100
    max_check: int = 500


@dataclass
class unionopt:
    bed: Optional[str] = None
    min_len: int = 100
    read_min_count: int = 2
    group_min_count: int = 5
    read_len_ratio: int = 0.8
    win_size: int = 500
    min_len_ratio: float = 0.6
    print_sv: bool = False
    collapsed: bool = False

@dataclass
class EvalOpt:
    min_len: int = 100
    min_count: int = 2
    read_len_ratio: float = 0.8
    win_size: int = 500
    min_len_ratio: float = 0.6
    dbg: bool = False
    print_err: bool = False
    print_all: bool = False
    bed: Optional[str] = None
    min_vaf: float = 0
    ignore_flt: bool = False
    check_gt: bool = False
    merge: bool = False
    svid: str = ""
    only_readname: bool = False


@click.group(help="minisv tool commands")
@click.version_option(__version__)
@click.pass_context
def cli(ctx):
    pass


#@cli.command()
#@click.option(
#    "-i",
#    "--input",
#    required=True,
#    type=click.Path(exists=True),
#    help="tumor or normal sample input gaf/paf file, this is required input, if paired normal sample provided, this should be input tumor sample",
#)
#@click.option(
#    "--vntr",
#    required=False,
#    type=click.Path(exists=True),
#    help="optional vntr sites for flagging potential false discovery",
#)
#@click.option(
#    "--cent",
#    required=False,
#    type=click.Path(exists=True),
#    help="optional centromere sites for flagging potential false discovery",
#)
#@click.option(
#    "--l1",
#    required=False,
#    type=click.Path(exists=True),
#    help="optional L1 sequence fasta for filtering L1 elements analysis",
#)
#@click.option(
#    "-r", "--support_read", required=True, type=int, help="supported read threshold"
#)
#@click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
#@click.option(
#    "-c", "--cpu", required=True, type=int, default=1, help="read mapping quality"
#)
#@click.option("-l", "--mlen", required=True, type=int, help="min indel length")
#@click.option("-a", "--maplen", required=True, type=int, help="min mapping length")
#@click.option(
#    "-n",
#    "--normal",
#    required=False,
#    type=str,
#    help="paired normal sample input gaf/paf file if tumor-only and normal-only mode, skip this option",
#)
#@click.option(
#    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
#)
#@click.option(
#    "-s",
#    "--assembly",
#    required=True,
#    type=str,
#    help="assembly version, e.g., chm13graph,chm13linear,grch37graph,grch37linear,grch38graph,grch38linear",
#)
#@click.option("-d", "--ds", is_flag=True, help="support ds tag or not")
#@click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
#def getsv(
#    input: str,
#    vntr: str,
#    cent: str,
#    l1: str,
#    support_read: int,
#    mapq: int,
#    cpu: int,
#    mlen: int,
#    maplen: int,
#    normal: str,
#    prefix: str,
#    ds: bool,
#    assembly: str,
#    verbose: bool,
#):
#    """Get indel reads and merge the microhomology reads into large indels"""
#    print("get indel...")
#    command = " ".join(sys.argv)
#    print(command)
#
#    input_samples = [input, normal] if normal is not None else [input]
#    gafs = GafParser(input_samples, prefix, assembly, vntr, cent, l1)
#    all_breakpts = gafs.parse_sv_on_group_reads(
#        mapq, mlen, maplen, verbose, n_cpus=cpu, ds=ds
#    )
#
#    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
#    merged_breakpt_vcf_list = gafs.merge_breakpts(all_breakpts)
#
#    # output bedpe format
#    merge_indel_breakpoints(prefix, gafs.breakpt_file, gafs.indel_file)
#
#    # output vcf format
#    gafs.bed2vcf(command=command, merged_breakpt_vcf_list=merged_breakpt_vcf_list)


@cli.command()
@click.option("-n", required=False, default="foo", type=str, help="sample name")
@click.option("-q", required=False, default=5, type=int, help="minimum mapping quality")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="centromere bed file"
)
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option(
    "-f", required=False, default=0.7, type=float, help="min fraction of reads"
)
@click.option("-c", required=False, default=3, help="maximum count of sv per 10k")
@click.option("-a", required=False, type=int, default=5, help="polyA penalty")
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option(
    "-x",
    required=False,
    default=30,
    type=int,
    help="minimum mapping quality in the end",
)
@click.option(
    "-e",
    required=False,
    default=2000,
    type=int,
    help="minimum aligned length in the end",
)
@click.option(
    "-m",
    required=False,
    default=50,
    type=int,
    help="minimum aligned length in the middle",
)
@click.argument("filename", nargs=1)
def getsv(
    q: int,
    x: int,
    svlen: str,
    d: bool,
    f: float,
    c: int,
    a: int,
    e: int,
    m: int,
    n: str,
    b: str,
    filename: tuple,
):
    """Extract raw INDEL and breakpoints"""
    from .read_parser import load_reads

    options = opt(
        min_mapq=q,
        min_mapq_end=x,
        min_len=parseNum(svlen),
        dbg=d,
        min_frac=f,
        max_cnt_10k=c,
        polyA_pen=a,
        min_aln_len_end=e,
        min_aln_len_mid=m,
        name=n,
        cen={},
    )

    if b is not None:
        parse_centromere(b, options.cen)
    load_reads(filename, options)


@cli.command()
@click.option("-n", required=False, default="foo", type=str, help="sample name")
@click.option("-q", required=False, default=5, type=int, help="minimum mapping quality")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="centromere bed file"
)
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option(
    "-f", required=False, default=0.7, type=float, help="min fraction of reads"
)
@click.option("-c", required=False, default=3, help="maximum count of sv per 10k")
@click.option("-a", required=False, type=int, default=5, help="polyA penalty")
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option(
    "-x",
    required=False,
    default=30,
    type=int,
    help="minimum mapping quality in the end",
)
@click.option(
    "-e",
    required=False,
    default=2000,
    type=int,
    help="minimum aligned length in the end",
)
@click.option(
    "-m",
    required=False,
    default=50,
    type=int,
    help="minimum aligned length in the middle",
)
@click.argument("filename", nargs=1)
def extract(
    q: int,
    x: int,
    svlen: str,
    d: bool,
    f: float,
    c: int,
    a: int,
    e: int,
    m: int,
    n: str,
    b: str,
    filename: tuple,
):
    """Extract raw INDEL and breakpoints"""
    from .read_parser import load_reads

    options = opt(
        min_mapq=q,
        min_mapq_end=x,
        min_len=parseNum(svlen),
        dbg=d,
        min_frac=f,
        max_cnt_10k=c,
        polyA_pen=a,
        min_aln_len_end=e,
        min_aln_len_mid=m,
        name=n,
        cen={},
    )

    if b is not None:
        parse_centromere(b, options.cen)
    load_reads(filename, options)


def parse_centromere(b, cen):
    with open(b) as centromere_file:
        for line in centromere_file:
            t = line.strip().split("\t")
            if t[0] not in cen:
                cen[t[0]] = []
            cen[t[0]].append([int(t[1]), int(t[2])])
        for ctg in cen:
            cen[ctg].sort(key=lambda x: x[0])


@cli.command()
@click.option("-w", required=False, default=100, type=int, help="window size")
@click.option(
    "-d",
    required=False,
    default=0.05,
    type=float,
    help="maximum allele length different ratio",
)
@click.option("-c", required=False, default=4, type=int, help="minimum sv counts")
@click.option(
    "-s", required=False, default=2, type=int, help="minimum count per strand"
)
@click.option(
    "-r",
    required=False,
    default=0,
    type=int,
    help="min min(TSD_len,polyA_len) to tag a candidate RT",
)
@click.option("-R", "--rt", required=False, default=1, type=int, help="minimum count for RT")
@click.option("-a", required=False, default=100, type=int, help="maximum allele number")
@click.option(
    "-C",
    "--allelecount",
    required=False,
    default=500,
    type=int,
    help="compare up to INT reads per allele",
)
@click.option(
    "-e",
    required=False,
    default="500k",
    type=str,
    help="minimum distance to centromere",
)
@click.argument("input", type=click.File("r"), required=True, nargs=1)
def merge(
    w: int, d: float, c: int, s: int, r: int, rt: int, a: int, allelecount: int, e: str, input
):
    """Usage: sort -k 1,1 -k2,2n sv.bed | gafcall merge [options] -"""
    from .merge import merge_sv

    options = mergeopt(
        win_size=w,
        max_diff=d,
        min_cnt=c,
        min_cnt_rt=rt,
        min_cen_dist=parseNum(e),
        min_rt_len=r,
        min_cnt_strand=s,
        max_allele=a,
        max_check=allelecount,
    )
    merge_sv(options, input)


@cli.command()
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option("-b", type=click.Path(exists=True), help="bed to restrict the comparison")
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option(
    "-s", "--size", is_flag=True, help="size stratified evaluation"
)
@click.option("-c", required=False, default=0, type=int, help="minimum sv counts")
@click.option("-w", required=False, default="500", type=str, help="window size")
@click.option(
    "-r",
    required=False,
    default=0.8,  # NOTE: in js this might be a bug, these ratio name looks confusing
    type=float,
    help="read SVs longer than svlen*FLOAT",
)
@click.option(
    "-m",
    required=False,
    default=0.6,
    type=float,
    help="minimum length ratio similarity between read SV and truthset SV",
)
@click.option(
    "-M",
    "--merge",
    is_flag=True,
    help="merge SV as consensus truth"
)
@click.option(
    "-v", required=False, default=0, type=float, help="ignore VAF below FLOAT"
)
@click.option("-F", "--ignoreflt", is_flag=True, help="ignore VCF filter")
@click.option("-G", "--gt", is_flag=True, help="check GT")
@click.option("-e", is_flag=True, help="print errors")
@click.option("-a", is_flag=True, help="print all")
@click.argument("filename", nargs=-1)
def eval(
    w: str,
    svlen: str,
    size: bool,
    c: int,
    r: float,
    m: float,
    d: bool,
    e: bool,
    a: bool,
    v: float,
    b,
    gt,
    ignoreflt,
    merge,
    filename,
):
    """Evaluation of SV calls"""
    from .eval import eval
    print(svlen)
    print(merge)
    print(size)
    options = EvalOpt(
        min_len=parseNum(svlen),
        min_count=int(c),
        win_size=parseNum(w),
        read_len_ratio=r,
        min_len_ratio=m,  # NOTE: the option are not input
        dbg=d,
        bed=gc_read_bed(b) if b is not None else None,
        print_err=e,
        print_all=a,
        min_vaf=v,
        check_gt=gt,
        merge=merge,
        ignore_flt=ignoreflt,
    )
    eval(filename, options, size)


@cli.command()
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option("-b", type=click.Path(exists=True), help="bed to restrict the comparison")
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option("-c", required=False, default=0, type=int, help="minimum sv counts")
@click.option("-w", required=False, default="500", type=str, help="window size")
@click.option(
    "-r",
    required=False,
    default=0.8,  # NOTE: in js this might be a bug, these ratio name looks confusing
    type=float,
    help="read SVs longer than svlen*FLOAT",
)
@click.option(
    "-m",
    required=False,
    default=0.6,
    type=float,
    help="minimum length ratio similarity between read SV and truthset SV",
)
@click.option(
    "-M",
    "--merge",
    is_flag=True,
    help="merge SV as consensus truth"
)
@click.option(
    "-v", required=False, default=0, type=float, help="ignore VAF below FLOAT"
)
@click.option("-F", "--ignoreflt", is_flag=True, help="ignore VCF filter")
@click.option("-i", "--consensusid", required=False, type=str, default="", help="consensus id file")
@click.option("-s", "--svid", required=False, type=str, default="", help="debug one svid")
@click.option("-G", "--gt", is_flag=True, help="check GT")
@click.option("-r", "--onlyreadname", is_flag=True, help="only filter by read name for vcf output")
@click.option("-e", is_flag=True, help="print errors")
@click.option("-a", is_flag=True, help="print all")
@click.argument("readidtsv", type=str, nargs=1)
@click.argument("msvasm", type=str, nargs=1)
@click.argument("outstat", type=str, nargs=1)
@click.argument("vcffile", nargs=1)
def filterasm(
    w: str,
    svlen: str,
    c: int,
    r: float,
    m: float,
    d: bool,
    e: bool,
    a: bool,
    v: float,
    consensusid: str,
    svid: str,
    b,
    gt,
    ignoreflt,
    merge,
    readidtsv,
    msvasm,
    outstat,
    vcffile,
    onlyreadname
):
    """Evaluation of SV calls"""
    from .filtercaller import othercaller_filterasm
    options = EvalOpt(
        only_readname=onlyreadname,
        min_len=parseNum(svlen),
        min_count=int(c),
        win_size=parseNum(w),
        read_len_ratio=r,
        min_len_ratio=m,  # NOTE: the option are not input
        dbg=d,
        bed=gc_read_bed(b) if b is not None else None,
        print_err=e,
        print_all=a,
        min_vaf=v,
        check_gt=gt,
        merge=merge,
        ignore_flt=ignoreflt,
        svid=svid
    )
    othercaller_filterasm(vcffile, options, readidtsv, msvasm, outstat, consensusid)


@cli.command()
@click.option(
    "-n", "--name", required=False, default="test", type=str, help="test"
)
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option(
   "-p", "--platform", required=False, default="hifi", type=str, help="sequencing platform (default:hifi)"
)
@click.option(
    "-r",
    "--ratio",
    required=False,
    default=0.8,  # NOTE: in js this might be a bug, these ratio name looks confusing
    type=float,
    help="read SVs longer than svlen*FLOAT",
)
@click.option("-c", required=False, default=3, type=int, help="minimum sv counts")
@click.option("-g", required=False, default=5, type=int, help="min group read count")
@click.option("--uc", required=False, default=5, type=int, help="union sv counts")
@click.option("-F", "--ignoreflt", is_flag=True, help="ignore VCF filter")
@click.option("-G", "--gt", is_flag=True, help="check GT")
@click.option("-w", required=False, default=500, type=int, help="window size")
@click.option("-m", required=False, default=0.6, type=float, help="min sv length ratio")
@click.option("-b", required=False, default=None, type=str, help="union and evaluated within the high confident bed file")
@click.option("--maskb", required=False, default=None, type=str, help="centromere bed file to mask the grch38-based sv signal")
@click.option("--mm2", required=True, default=None, type=str, help="minimap2 path")
@click.option("--mg", required=True, default=None, type=str, help="minigraph path")
@click.option("--vcf", type=str, nargs=4, required=True, help="Exactly 4 VCFs following the order of severus savana nanomonsv sniffles2")
@click.option("--readid_tsv", type=str, nargs=4, required=True, help="Exactly 4 read-id TSV files corresponding to the VCFs")
@click.argument("bamfile", type=str, nargs=1)
@click.argument("ref", type=str, nargs=1)
@click.argument("hap1_denovo_ref", type=str, nargs=1)
@click.argument("hap2_denovo_ref", type=str, nargs=1)
@click.argument("graph_ref", type=str, nargs=1)
@click.argument("workdir", nargs=1)
def sv_cross_ref_filter(
    name, svlen, platform, ratio, c, uc, g,
    ignoreflt, gt, w, m, b, maskb,
    mm2, mg,
    vcf,
    readid_tsv,
    bamfile,
    ref,
    hap1_denovo_ref,
    hap2_denovo_ref,
    graph_ref,
    workdir
):
    """use samtools + seqtk + mappy to filter somatic SV with real-time alignment breakpoint evidences to denovo assembly """
    print(vcf)
    print(readid_tsv)

    assert len(vcf) == 4
    assert len(vcf) == len(readid_tsv)

    #vcf = vcf[0]
    #readid_tsv = readid_tsv[0]
    from .filtercaller import MinisvReads
    asm_read_cutoff = c

    minisv_reads = MinisvReads(vcf, readid_tsv, bamfile, ref, hap1_denovo_ref, hap2_denovo_ref, graph_ref, workdir, asm_read_cutoff, platform, mm2, mg)
    min_len = parseNum(svlen)

    min_read_len = math.floor(min_len * ratio + 0.499)
    minisv_reads.extract_read_ids(b, min_len, min_read_len, c, ignoreflt, gt)
    minisv_reads.extract_reads()

    #NOTE:svcall+l
    #minisv_reads.align_reads_to_grch38()
  
    #use svcall+s
    minisv_reads.align_reads_to_self()
    minisv_reads.align_reads_to_graph()

    # getsv options
    options = opt(
        min_mapq=0,
        min_mapq_end=0,
        min_len=min_len,
        min_frac=0.7,
        max_cnt_10k=5,
        polyA_pen=5,
        min_aln_len_end=2000,
        min_aln_len_mid=50,
        name=name,
        cen={}
    )
    minisv_reads.parse_raw_sv_self(options)
    if maskb is not None:
        parse_centromere(maskb, options.cen)
    if b is not None:
        options.bed = b

    options.min_mapq = 5
    options.min_mapq_end = 30

    #minisv_reads.parse_raw_sv_grch38(options)
    minisv_reads.parse_raw_sv_graph(options)

    # only svcall + g
    #      svcall + s
    #      svcall + gs
    #minisv_reads.isec_g(options, w)
    #minisv_reads.isec_s(options, w)
    minisv_reads.isec_gs(options, w)

    print(w)
    for cutoff in range(int(c), int(uc)+1):
        options = EvalOpt(
            only_readname=False,
            min_len=parseNum(svlen),
            min_count=cutoff,
            win_size=parseNum(str(w)),
            read_len_ratio=ratio,
            min_len_ratio=m,  # NOTE: the option are not input
            dbg=False,
            bed=gc_read_bed(b) if b is not None else None,
            print_err=False,
            print_all=True,
            min_vaf=0,
            check_gt=False,
            merge=False,
            ignore_flt=False,
            svid=""
        )
        minisv_reads.othercaller_filterasm(options)

    options = unionopt(
        bed=b,
        min_len=min_len,
        read_min_count=uc, # 
        group_min_count=g,
        read_len_ratio=ratio,
        win_size=w,
        min_len_ratio=m,
        print_sv=True
    )
    minisv_reads.union_filtered_vcf(min_read_len, options)
    minisv_reads.save_timings()


@cli.command()
@click.option(
    "-n", "--name", required=False, default="test", type=str, help="test"
)
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option(
   "-p", "--platform", required=False, default="hifi", type=str, help="sequencing platform (default:hifi)"
)
@click.option(
    "-r",
    "--ratio",
    required=False,
    default=0.8,  # NOTE: in js this might be a bug, these ratio name looks confusing
    type=float,
    help="read SVs longer than svlen*FLOAT",
)
@click.option("-c", required=False, default=3, type=int, help="minimum sv counts")
@click.option("-g", required=False, default=5, type=int, help="min group read count")
@click.option("--uc", required=False, default=5, type=int, help="union sv counts")
@click.option("-F", "--ignoreflt", is_flag=True, help="ignore VCF filter")
@click.option("-G", "--gt", is_flag=True, help="check GT")
@click.option("-w", required=False, default=500, type=int, help="window size")
@click.option("-m", required=False, default=0.6, type=float, help="min sv length ratio")
@click.option("-b", required=False, default=None, type=str, help="union and evaluated within the high confident bed file")
@click.option("--maskb", required=False, default=None, type=str, help="centromere bed file to mask the grch38-based sv signal")
@click.option("--mm2", required=True, default=None, type=str, help="minimap2 path")
@click.option("--mg", required=True, default=None, type=str, help="minigraph path")
@click.option("--vcf", type=str, nargs=1, required=True, help="Customized sv calling from sniffles2 that map SV to read names")
@click.option("--readid_tsv", type=str, nargs=1, required=True, help="the same as --vcf")
@click.argument("bamfile", type=str, nargs=1)
@click.argument("ref", type=str, nargs=1)
@click.argument("dad_hap1_denovo_ref", type=str, nargs=1)
@click.argument("dad_hap2_denovo_ref", type=str, nargs=1)
@click.argument("mom_hap1_denovo_ref", type=str, nargs=1)
@click.argument("mom_hap2_denovo_ref", type=str, nargs=1)
@click.argument("workdir", nargs=1)
def sv_trio_filter(
    name, svlen, platform, ratio, c, uc, g,
    ignoreflt, gt, w, m, b, maskb,
    mm2, mg,
    vcf,
    readid_tsv,
    bamfile,
    ref,
    dad_hap1_denovo_ref,
    dad_hap2_denovo_ref,
    mom_hap1_denovo_ref,
    mom_hap2_denovo_ref,
    workdir
):
    """use samtools + seqtk + mappy to filter somatic SV with real-time alignment breakpoint evidences to denovo assembly """
    print(vcf)
    print(readid_tsv)

    #vcf = vcf[0]
    #readid_tsv = readid_tsv[0]
    from .filtercaller import MinisvReadsTrio
    asm_read_cutoff = c

    minisv_reads = MinisvReadsTrio(vcf, readid_tsv, bamfile, ref, dad_hap1_denovo_ref, dad_hap2_denovo_ref, mom_hap1_denovo_ref, mom_hap2_denovo_ref, workdir, asm_read_cutoff, platform, mm2, mg)
    min_len = parseNum(svlen)

    min_read_len = math.floor(min_len * ratio + 0.499)
    minisv_reads.extract_read_ids(b, min_len, min_read_len, c, ignoreflt, gt)
    minisv_reads.extract_reads()

    #NOTE:svcall+l
    #minisv_reads.align_reads_to_grch38()
  
    #use svcall+s
    minisv_reads.align_reads_to_self()
    #minisv_reads.align_reads_to_graph()

    # getsv options
    options = opt(
        min_mapq=0,
        min_mapq_end=0,
        min_len=min_len,
        min_frac=0.7,
        max_cnt_10k=5,
        polyA_pen=5,
        min_aln_len_end=2000,
        min_aln_len_mid=50,
        name=name,
        cen={}
    )
    minisv_reads.parse_raw_sv_self(options)

    options.min_mapq = 5
    options.min_mapq_end = 30

    #minisv_reads.parse_raw_sv_grch38(options)
    #minisv_reads.parse_raw_sv_graph(options)

    # only svcall + g
    #      svcall + s
    #      svcall + gs
    #minisv_reads.isec_g(options, w)
    #minisv_reads.isec_s(options, w)
    #minisv_reads.isec_gs(options, w)

    print(w)
    for cutoff in range(int(c), int(uc)+1):
        options = EvalOpt(
            only_readname=False,
            min_len=parseNum(svlen),
            min_count=cutoff,
            win_size=parseNum(str(w)),
            read_len_ratio=ratio,
            min_len_ratio=m,  # NOTE: the option are not input
            dbg=False,
            bed=gc_read_bed(b) if b is not None else None,
            print_err=False,
            print_all=True,
            min_vaf=0,
            check_gt=False,
            merge=False,
            ignore_flt=False,
            svid=""
        )
        minisv_reads.othercaller_filterasm(options)

    ####options = unionopt(
    ####    bed=b,
    ####    min_len=min_len,
    ####    read_min_count=uc, # 
    ####    group_min_count=g,
    ####    read_len_ratio=ratio,
    ####    win_size=w,
    ####    min_len_ratio=m,
    ####    print_sv=True
    ####)
    ####minisv_reads.union_filtered_vcf(min_read_len, options)
    minisv_reads.save_timings()


@cli.command()
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option("-b", type=click.Path(exists=True), help="bed to restrict the comparison")
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option("-c", required=False, default=0, type=int, help="minimum sv counts")
@click.option("-w", required=False, default="500", type=str, help="window size")
@click.option(
    "-r",
    required=False,
    default=0.8,  # NOTE: in js this might be a bug, these ratio name looks confusing
    type=float,
    help="read SVs longer than svlen*FLOAT",
)
@click.option(
    "-m",
    required=False,
    default=0.6,
    type=float,
    help="minimum length ratio similarity between read SV and truthset SV",
)
@click.option(
    "-M",
    "--merge",
    is_flag=True,
    help="merge SV as consensus truth"
)
@click.option(
    "-v", required=False, default=0, type=float, help="ignore VAF below FLOAT"
)
@click.option("-F", "--ignoreflt", is_flag=True, help="ignore VCF filter")
@click.option("-G", "--gt", is_flag=True, help="check GT")
@click.option("-e", is_flag=True, help="print errors")
@click.option("-a", is_flag=True, help="print all")
@click.argument("filename", nargs=-1)
def annot(
    w: str,
    svlen: str,
    c: int,
    r: float,
    m: float,
    d: bool,
    e: bool,
    a: bool,
    v: float,
    b,
    gt,
    ignoreflt,
    merge,
    filename,
):
    """Evaluation of SV calls"""
    from .annot import annot

    options = EvalOpt(
        min_len=parseNum(svlen),
        min_count=int(c),
        win_size=parseNum(w),
        read_len_ratio=r,
        min_len_ratio=m,  # NOTE: the option are not input
        bed=None,
        print_err=e,
        print_all=a,
        min_vaf=v,
        check_gt=gt,
        merge=merge,
        ignore_flt=ignoreflt,
    )
    annot(filename, options)


@dataclass
class viewopt:
    min_read_len: int
    min_count: int
    ignore_flt: bool
    check_gt: bool
    count_long: bool
    bed: Optional[str]


@cli.command()
@click.option(
    "-minlen", required=False, default=100, type=int, help="minimum sv length"
)
@click.option("-ignoreflt", is_flag=True, help="ignore FILTER field in VCF")
@click.option("-gt", is_flag=True, help="check GT in VCF")
@click.option("-c", "--count", default=0, type=int, help="minimum count")
@click.option("-C", "--countlong", is_flag=True, help="count 20kb,100kb,1Mb and translocations")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="restricted bed file"
)
@click.argument("input", nargs=-1)
def view(minlen: int, ignoreflt: bool, gt: bool, countlong: bool, count: int, b: str, input):
    opt = viewopt(
        min_read_len=minlen, ignore_flt=ignoreflt, check_gt=gt, count_long=countlong, bed=b,
        min_count=count
    )
    gc_cmd_view(opt, input)


@cli.command(context_settings={"ignore_unknown_options": True})
@click.option(
    "-w", required=False, default=1000, type=int, help="window size"
)
@click.argument("gsvs", type=click.Path(), nargs=-1, required=True)
def isec(w, gsvs):
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
                print(line.strip())
                    


@cli.command()
@click.argument("input", type=click.File("r"), default=sys.stdin)
def vcf(input):
    for line in input:
        if line.startswith("#CHROM"):
            print(
                """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotyping quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">"""
            )
            print(line.strip())
        elif line.startswith("##"):
            print(line.strip())
        else:
            line = line.split()
            line[6] = "PASS"
            line = "\t".join(line)
            # print(line.strip(), "GT:GQ:VAF:DR:DV", "0/1:.:.:.:.", sep="\t")
            print(line.strip(), sep="\t")


@cli.command()
@click.argument("input", type=click.File("r"), default=sys.stdin)
def genvcf(input):
    options = None
    write_vcf(options, input)


@cli.command()
@click.argument("input", type=click.File("r"), default=sys.stdin)
def extracthp(input):
    options = None
    extract_phase_HP(input)


@cli.command()
@click.argument("hptagtsv", type=click.File("r"), default=sys.stdin, nargs=1)
@click.argument("msv", type=str, nargs=1)
def annotatehp(hptagtsv, msv):
    annotate_HP(hptagtsv, msv)


@cli.command()
@click.argument("severusvcf", type=str, nargs=1)
@click.argument("readidtsv", type=str, nargs=1)
@click.argument("msvasm", type=str, nargs=1)
@click.argument("outstat", type=str, nargs=1)
@click.argument("consensus_ids", type=str, nargs=1)
def filterseverus(severusvcf, readidtsv, msvasm, outstat, consensus_ids):
    """
    filter Severus results based on read ids overlap with graph alignment/self alignment
    """
    call_filterseverus(severusvcf, readidtsv, msvasm, outstat, consensus_ids)


@cli.command()
@click.argument("snfvcf", type=str, nargs=1)
@click.argument("msvasm", type=str, nargs=1)
@click.argument("outstat", type=str, nargs=1)
def filtersnf(snfvcf, msvasm, outstat):
    """
    filter Sniffles2 results based on read ids overlap with graph alignment/self alignment
    """
    call_filtersnf(snfvcf, msvasm, outstat)


@cli.command()
@click.option(
    "-t",
    required=False,
    default=1,
    type=int,
    help="tumor vcf sample index"
)
@click.option(
    "-n",
    required=False,
    default=2,
    type=int,
    help="normal vcf sample index"
)
@click.argument("snfvcf", type=str, nargs=1)
def snfpair(snfvcf, t, n):
    t_index = t
    n_index = n
    with gzip.open(snfvcf, 'rt') as inf:
        for line in inf:
            if line[0].startswith("#"):
                print(line.strip())
                continue
            t = line.strip().split("\t")
            sn = t[8+n_index].split(':')
            st = t[8+t_index].split(':')
            fmt = t[8].split(':')
            for i in range(len(fmt)):
                if fmt[i] == "DV" and int(st[i]) > 0 and (sn[i] == "NA" or int(sn[i]) == 0):
                    print(line.strip())
                    continue
         


@cli.command()
@click.argument("msvtg", type=str, nargs=1)
@click.argument("msvasm", type=str, nargs=1)
@click.argument("outstat", type=str, nargs=1)
def filtermsv(msvtg, msvasm, outstat):
    """
    filter Sniffles2 results based on read ids overlap with graph alignment/self alignment
    """
    call_filtermsv(msvtg, msvasm, outstat)


@cli.command()
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="tumor or normal sample input gaf/paf file, this is required input, if paired normal sample provided, this should be input tumor sample",
)
@click.option(
    "--vntr",
    required=False,
    type=click.Path(exists=True),
    help="optional vntr sites for flagging potential false discovery",
)
@click.option(
    "--cent",
    required=False,
    type=click.Path(exists=True),
    help="optional centromere sites for flagging potential false discovery",
)
@click.option(
    "--l1",
    required=False,
    type=click.Path(exists=True),
    help="optional L1 sequence fasta for filtering L1 elements analysis",
)
@click.option(
    "-r", "--support_read", required=True, type=int, help="supported read threshold"
)
@click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
@click.option(
    "-c", "--cpu", required=True, type=int, default=1, help="read mapping quality"
)
@click.option("-l", "--mlen", required=True, type=int, help="min indel length")
@click.option("-a", "--maplen", required=True, type=int, help="min mapping length")
@click.option(
    "-n",
    "--normal",
    required=False,
    type=str,
    help="paired normal sample input gaf/paf file if tumor-only and normal-only mode, skip this option",
)
@click.option(
    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
)
@click.option(
    "-s",
    "--assembly",
    required=True,
    type=str,
    help="assembly version, e.g., chm13graph,chm13linear,grch37graph,grch37linear,grch38graph,grch38linear",
)
@click.option("-d", "--ds", is_flag=True, help="support ds tag or not")
@click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
def getindel(
    input: str,
    vntr: str,
    cent: str,
    l1: str,
    support_read: int,
    mapq: int,
    cpu: int,
    mlen: int,
    maplen: int,
    normal: str,
    prefix: str,
    ds: bool,
    assembly: str,
    verbose: bool,
):
    """Get indel reads and merge the microhomology reads into large indels"""
    print("get indel...")
    command = " ".join(sys.argv)
    print(command)

    input_samples = [input, normal] if normal is not None else [input]
    gafs = GafParser(input_samples, prefix, assembly, vntr, cent, l1)
    gafs.parse_indel(mapq, mlen, verbose, n_cpus=cpu, ds=ds)
    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)


@cli.command()
@click.argument("msvunion", type=str, nargs=1)
def ensembleunion(msvunion):
    """
    ensemble and collapse the minisv.js union results
    """
    insilico_truth(msvunion)



@cli.command()
@click.option("-b", required=False, default=None, type=str, help="evaluated within the bed file")
@click.option("-l", required=False, default=100, type=int, help="minimum length")
@click.option("-c", required=False, default=5, type=int, help="min read count")
@click.option("-g", required=False, default=5, type=int, help="min group read count")
@click.option("-r", required=False, default=0.8, type=float, help="read length ratio")
@click.option("-w", required=False, default=500, type=int, help="window size")
@click.option("-m", required=False, default=0.6, type=float, help="min sv length ratio")
@click.option("-p", is_flag=True, help="print union of sv calls or not")
@click.argument("filename", nargs=-1, required=True)
def union(
    b, l, c, g, r, w, m, p, filename
):

    options = unionopt(
        bed=b,
        min_len=l,
        read_min_count=c,
        group_min_count=g,
        read_len_ratio=r,
        win_size=w,
        min_len_ratio=m,
        print_sv=p
    )

    read_min_len = math.floor(options.min_len * options.read_len_ratio + 0.499)
    union_sv(filename, read_min_len, options)


@cli.command()
@click.option("-b", required=False, default=None, type=str, help="evaluated within the bed file")
@click.option("-l", required=False, default=100, type=int, help="minimum length")
@click.option("-c", required=False, default=2, type=int, help="min read count")
@click.option("-g", required=False, default=5, type=int, help="min group read count")
@click.option("-r", required=False, default=0.8, type=float, help="read length ratio")
@click.option("-w", required=False, default=500, type=int, help="window size")
@click.option("-m", required=False, default=0.6, type=float, help="min sv length ratio")
@click.option("-p", is_flag=True, help="print union of sv calls or not")
@click.option("-u", is_flag=True, help="take only one sv call for one group of svs that has maximum read count from assembly")
@click.option('--input1', '-i1', type=click.Path(exists=True), multiple=True, help='List of SV call files')
@click.option('--input2', '-i2', type=click.Path(exists=True), multiple=True, help='List of read ids for the SV call files')
@click.argument("asmgsv", nargs=1, required=True)
def advunion(
    b, l, c, g, r, w, m, p, u,
    input1,
    input2, 
    asmgsv
):
    options = unionopt(
        bed=b,
        min_len=l,
        read_min_count=c,
        group_min_count=g,
        read_len_ratio=r,
        win_size=w,
        min_len_ratio=m,
        print_sv=p,
        collapsed=u
    )
    read_min_len = math.floor(options.min_len * options.read_len_ratio + 0.499)
    advunion_sv(input2, input1, asmgsv, read_min_len, options)


@cli.command()
@click.option("-b", required=False, default=None, type=str, help="evaluated within the bed file")
@click.option("-l", required=False, default=100, type=int, help="minimum length")
@click.option("-c", required=False, default=2, type=int, help="min read count")
@click.option("-g", required=False, default=5, type=int, help="min group read count")
@click.option("-r", required=False, default=0.8, type=float, help="read length ratio")
@click.option("-w", required=False, default=500, type=int, help="window size")
@click.option("-m", required=False, default=0.6, type=float, help="min sv length ratio")
@click.option("-p", is_flag=True, help="print union of sv calls or not")
@click.argument("filename", nargs=-1)
def union_tr(
    b, l, c, g, r, w, m, p, filename
):

    options = unionopt(
        bed=b,
        min_len=l,
        read_min_count=c,
        group_min_count=g,
        read_len_ratio=r,
        win_size=w,
        min_len_ratio=m,
        print_sv=p
    )

    read_min_len = math.floor(options.min_len * options.read_len_ratio + 0.499)
    union_sv_with_tr(filename, read_min_len, options)


@cli.command()
@click.argument("collapsedmsvunion", type=str, nargs=1)
def doublestrandbreak(collapsedmsvunion):
    """
    filter Sniffles2 results based on read ids overlap with graph alignment/self alignment
    """
    double_strand_break(collapsedmsvunion)


# @cli.command()
# @click.option(
#    "-i",
#    "--input",
#    required=True,
#    type=click.Path(exists=True),
#    help="input gaf/paf file",
# )
# @click.option(
#    "-r", "--support_read", required=True, type=int, help="supported read threshold"
# )
# @click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
# @click.option(
#    "-c", "--cpu", required=True, type=int, default=1, help="read mapping quality"
# )
# @click.option("-l", "--mlen", required=True, type=int, help="min indel length")
# @click.option(
#    "-n",
#    "--normal",
#    required=False,
#    type=str,
#    help="paired normal sample input gaf/paf file if tumor-only and normal-only mode, skip this option",
# )
# @click.option(
#    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
# )
# @click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
# def getindel_cython(
#    input: str,
#    support_read: int,
#    mapq: int,
#    cpu: int,
#    mlen: int,
#    normal: str,
#    prefix: str,
#    verbose: bool,
# ):
#    """Get indel reads and merge the microhomology reads into large indels"""
#    print("get cython indel...")
#    command = " ".join(sys.argv)
#    input_samples = [input, normal] if normal is not None else [input]
#    gafs = cy_GafParser(input_samples, prefix)
#    gafs.parse_indel(mapq, mlen, verbose, n_cpus=cpu)
#    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
#    gafs.bed2vcf(command=command)
