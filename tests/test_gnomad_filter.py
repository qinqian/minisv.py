"""Unit tests for iit_overlap (eval.py) and gnomad_filter (filtercaller.py)."""
import io
from types import SimpleNamespace

import pytest

from minisv.eval import gc_read_bed, iit_index, iit_overlap, iit_sort_copy
from minisv.filtercaller import gnomad_filter


def build_tree(intervals):
    """Build an indexed implicit interval tree the same way gc_read_bed does."""
    a = [{"st": st, "en": en, "data": None} for st, en in intervals]
    a = iit_sort_copy(a)
    iit_index(a)
    return a


def test_iit_overlap_hit():
    tree = build_tree([(100, 200), (300, 400)])
    assert len(iit_overlap(tree, 150, 151)) == 1
    assert len(iit_overlap(tree, 350, 360)) == 1


def test_iit_overlap_miss():
    tree = build_tree([(100, 200), (300, 400)])
    assert len(iit_overlap(tree, 250, 260)) == 0


def test_iit_overlap_halfopen_boundary():
    # interval is half-open [100, 200): position 199 overlaps, 200 does not
    tree = build_tree([(100, 200)])
    assert len(iit_overlap(tree, 199, 200)) == 1
    assert len(iit_overlap(tree, 200, 201)) == 0


def test_iit_overlap_multiple():
    tree = build_tree([(100, 200), (150, 250)])
    assert len(iit_overlap(tree, 160, 170)) == 2


# A 10 bp REF deletion: gc_parse_sv puts end1 at POS-1 and end2 at POS-1+rlen.
# POS=1000 -> end1=999, end2=1009;  POS=5000 -> end1=4999, end2=5009.
VCF_LINES = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1000\tdel_in\tATCGATCGAT\tA\t60\tPASS\tSVTYPE=DEL;SVLEN=-9",
    "chr1\t5000\tdel_out\tATCGATCGAT\tA\t60\tPASS\tSVTYPE=DEL;SVLEN=-9",
]

# gnomAD-shaped BED rows have 49 columns; CHR2==chrom routes to bed_pt (same-chrom).
def write_files(tmp_path, bed_regions):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text("\n".join(VCF_LINES) + "\n")
    bed = tmp_path / "gnomad.bed"
    lines = []
    for c, s, e in bed_regions:
        cols = ["NA"] * 49
        cols[0] = c
        cols[1] = str(s)
        cols[2] = str(e)
        cols[4] = "DEL"
        cols[9] = c        # CHR2 == chrom -> same-chrom -> bed_pt
        cols[43] = "9"
        cols[48] = "0.5"   # AF above default threshold
        lines.append("\t".join(cols))
    bed.write_text("\n".join(lines) + "\n")
    return str(vcf), str(bed)


_OPT = SimpleNamespace(
    ignore_flt=True, min_len=0, min_count=0, min_vaf=0,
    read_len_ratio=0.8, win_size=500, min_len_ratio=0.6,
    dbg=False, print_err=False, print_all=False,
    bed=None, gnomad_af=0.01, check_gt=False,
)


def run_filter(tmp_path, bed_regions, **kw):
    vcf, bed = write_files(tmp_path, bed_regions)
    out = io.StringIO()
    gnomad_filter(vcf, bed, _OPT, out=out, **kw)
    return out.getvalue()


def test_gnomad_filter_drops_overlap_keeps_outside(tmp_path):
    # region covers both ends of del_in (999, 1009), nothing near del_out
    out = run_filter(tmp_path, [("chr1", 990, 1010)])
    assert "del_in" not in out
    assert "del_out" in out
    # headers preserved
    assert "##fileformat=VCFv4.2" in out
    assert out.startswith("##fileformat")


def test_gnomad_filter_both_ends(tmp_path):
    # region [995,1005) covers only end1 (999), not end2 (1009) of del_in
    either = run_filter(tmp_path, [("chr1", 995, 1005)])
    assert "del_in" not in either  # default: one end is enough to drop

    both = run_filter(tmp_path, [("chr1", 995, 1005)], both_ends=True)
    assert "del_in" in both  # both_ends: a single overlapping end is not enough


def test_gnomad_filter_pad(tmp_path):
    # region [5015,5025) does not touch del_out (end2=5009) without padding
    no_pad = run_filter(tmp_path, [("chr1", 5015, 5025)])
    assert "del_out" in no_pad

    padded = run_filter(tmp_path, [("chr1", 5015, 5025)], pad=20)
    assert "del_out" not in padded  # pad expands the query to reach 5009


# ---------------------------------------------------------------------------
# Helper for gnomAD-shaped BED (49 cols, cross-chrom or same-chrom rows)
# ---------------------------------------------------------------------------

def _make_gnomad_bed(tmp_path, rows):
    """Write a minimal gnomAD-shaped BED (49 cols). rows = list of dicts."""
    lines = []
    for r in rows:
        cols = ["NA"] * 49
        cols[0]  = r["chrom"]
        cols[1]  = str(r["start"])
        cols[2]  = str(r["end"])
        cols[3]  = r.get("name", "gnomAD_test")
        cols[4]  = r.get("svtype", "BND")
        cols[9]  = r.get("chr2", "NA")
        cols[13] = str(r["end2"]) if "end2" in r else "NA"
        cols[43] = str(r.get("svlen", "NA"))
        cols[48] = str(r.get("af", 0.5))
        lines.append("\t".join(cols))
    p = tmp_path / "gnomad.bed"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


# ---------------------------------------------------------------------------
# gc_read_bed unit tests (Task 2)
# ---------------------------------------------------------------------------

def test_gc_read_bed_gnomad_same_chrom_goes_to_pt(tmp_path):
    """Same-chrom rows land in bed_pt; bed_sv is empty."""
    bed = _make_gnomad_bed(tmp_path, [
        {"chrom": "chr1", "start": 1000, "end": 2000, "svtype": "DEL",
         "chr2": "chr1", "af": 0.5},
    ])
    h_pt, h_sv = gc_read_bed(bed, is_gnomad=True)
    assert "chr1" in h_pt
    assert "chr1" not in h_sv


def test_gc_read_bed_gnomad_cross_chrom_goes_to_sv(tmp_path):
    """Cross-chrom rows land in bed_sv on BOTH contigs; bed_pt is empty."""
    bed = _make_gnomad_bed(tmp_path, [
        {"chrom": "chr1", "start": 1000, "end": 1001, "svtype": "BND",
         "chr2": "chr2", "end2": 2000, "af": 0.5},
    ])
    h_pt, h_sv = gc_read_bed(bed, is_gnomad=True)
    assert "chr1" not in h_pt
    assert "chr1" in h_sv
    assert "chr2" in h_sv
    node = h_sv["chr1"][0]
    assert node.data.ctg2 == "chr2"
    assert node.data.pos2 == 2000


def test_gc_read_bed_gnomad_af_filter(tmp_path):
    """Row with AF < gnomad_af is skipped; row above threshold is kept."""
    bed = _make_gnomad_bed(tmp_path, [
        {"chrom": "chr1", "start": 1000, "end": 2000, "svtype": "DEL",
         "chr2": "chr1", "af": 0.005},
    ])
    h_pt, h_sv = gc_read_bed(bed, is_gnomad=True, gnomad_af=0.01)
    assert "chr1" not in h_pt  # filtered out

    h_pt2, _ = gc_read_bed(bed, is_gnomad=True, gnomad_af=0.001)
    assert "chr1" in h_pt2    # passes lower threshold


@pytest.mark.parametrize("raw_svlen", ["NA", "", "-1"])
def test_gc_read_bed_gnomad_missing_svlen_is_zero(tmp_path, raw_svlen):
    bed = _make_gnomad_bed(tmp_path, [
        {"chrom": "chr1", "start": 1000, "end": 2000, "svtype": "DUP",
         "chr2": "chr1", "svlen": raw_svlen, "af": 0.5},
    ])
    h_pt, _ = gc_read_bed(bed, is_gnomad=True)
    assert h_pt["chr1"][0].data.SVLEN == 0


def test_gc_read_bed_non_gnomad_unchanged(tmp_path):
    """Non-gnomAD call (is_gnomad=False) returns a single dict, not a tuple."""
    p = tmp_path / "plain.bed"
    p.write_text("chr1\t100\t200\n")
    result = gc_read_bed(str(p))
    assert isinstance(result, dict)
    assert "chr1" in result


# ---------------------------------------------------------------------------
# gnomad_filter integration tests (Task 3)
# ---------------------------------------------------------------------------

VCF_BND = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    # BND chr1:1000 -> chr2:2000
    "chr1\t1001\tbnd_match\tN\tN[chr2:2000[\t60\tPASS\tSVTYPE=BND;SVLEN=0",
    # BND chr2:2000 -> chr3:7000  (mate-mismatch case)
    "chr2\t2001\tbnd_mismatch\tN\tN[chr3:7000[\t60\tPASS\tSVTYPE=BND;SVLEN=0",
]


def _run_filter_gnomad(tmp_path, gnomad_rows, vcf_lines, gnomad_af=0.01, **kw):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text("\n".join(vcf_lines) + "\n")
    bed = _make_gnomad_bed(tmp_path, gnomad_rows)
    out = io.StringIO()
    opt = SimpleNamespace(
        ignore_flt=True, min_len=0, min_count=0, min_vaf=0,
        read_len_ratio=0.8, win_size=500, min_len_ratio=0.6,
        dbg=False, print_err=False, print_all=False,
        bed=None, gnomad_af=gnomad_af, check_gt=False,
    )
    gnomad_filter(str(vcf), bed, opt, out=out, **kw)
    return out.getvalue()


def test_cross_chrom_sv_match_drops(tmp_path):
    """Caller BND chr1:1000<->chr2:2000 matches gnomAD row -> dropped."""
    out = _run_filter_gnomad(tmp_path, [
        {"chrom": "chr1", "start": 1000, "end": 1001,
         "svtype": "BND", "chr2": "chr2", "end2": 2000, "af": 0.5},
    ], VCF_BND)
    assert "bnd_match" not in out


def test_cross_chrom_mate_mismatch_kept(tmp_path):
    """
    Caller BND chr2:2000<->chr3:7000 — no single gnomAD record connects
    those two contigs. Must NOT be dropped.
    gnomAD has chr1:1000<->chr2:2000 and chr1:5000<->chr3:7000.
    """
    out = _run_filter_gnomad(tmp_path, [
        {"chrom": "chr1", "start": 1000, "end": 1001,
         "svtype": "BND", "chr2": "chr2", "end2": 2000, "af": 0.5},
        {"chrom": "chr1", "start": 5000, "end": 5001,
         "svtype": "BND", "chr2": "chr3", "end2": 7000, "af": 0.5},
    ], VCF_BND)
    assert "bnd_mismatch" in out


def test_same_chrom_same_type_similar_length_drops(tmp_path):
    """Same-chrom DEL with matching type and similar length is dropped."""
    out = _run_filter_gnomad(tmp_path, [
        {"chrom": "chr1", "start": 990, "end": 1010,
         "svtype": "DEL", "svlen": 9, "chr2": "chr1", "af": 0.5},
    ], VCF_LINES)
    assert "del_in" not in out
    assert "del_out" in out


def test_same_chrom_different_type_kept(tmp_path):
    """Same-chrom overlap with a different SVTYPE is kept."""
    out = _run_filter_gnomad(tmp_path, [
        {"chrom": "chr1", "start": 990, "end": 1010,
         "svtype": "DUP", "svlen": 9, "chr2": "chr1", "af": 0.5},
    ], VCF_LINES)
    assert "del_in" in out
    assert "del_out" in out


def test_same_chrom_bnd_does_not_match_del(tmp_path):
    """Same-chrom BND does not act as a wildcard for non-BND calls."""
    out = _run_filter_gnomad(tmp_path, [
        {"chrom": "chr1", "start": 990, "end": 1010,
         "svtype": "BND", "svlen": "NA", "chr2": "chr1", "af": 0.5},
    ], VCF_LINES)
    assert "del_in" in out
    assert "del_out" in out


def test_same_chrom_same_type_different_length_kept(tmp_path):
    """Same-chrom overlap with incompatible SVLEN is kept."""
    out = _run_filter_gnomad(tmp_path, [
        {"chrom": "chr1", "start": 990, "end": 1010,
         "svtype": "DEL", "svlen": 5000, "chr2": "chr1", "af": 0.5},
    ], VCF_LINES)
    assert "del_in" in out
    assert "del_out" in out


def test_gnomadaf_flag_cross_chrom(tmp_path):
    """Cross-chrom row at AF=0.005 ignored by default; dropped at gnomad_af=0.001."""
    gnomad_rows = [
        {"chrom": "chr1", "start": 1000, "end": 1001,
         "svtype": "BND", "chr2": "chr2", "end2": 2000, "af": 0.005},
    ]
    out_default = _run_filter_gnomad(tmp_path, gnomad_rows, VCF_BND)
    assert "bnd_match" in out_default  # AF too low -> row not indexed -> not dropped

    tmp2 = tmp_path / "sub2"
    tmp2.mkdir()
    out_low = _run_filter_gnomad(tmp2, gnomad_rows, VCF_BND, gnomad_af=0.001)
    assert "bnd_match" not in out_low  # gnomad_af=0.001 -> row included -> dropped
