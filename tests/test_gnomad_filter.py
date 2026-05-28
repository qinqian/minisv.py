"""Unit tests for iit_overlap (eval.py) and gnomad_filter (filtercaller.py)."""
import io
from types import SimpleNamespace

from minisv.eval import iit_index, iit_overlap, iit_sort_copy
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


def write_files(tmp_path, bed_regions):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text("\n".join(VCF_LINES) + "\n")
    bed = tmp_path / "gnomad.bed"
    bed.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in bed_regions))
    return str(vcf), str(bed)


def run_filter(tmp_path, bed_regions, **kw):
    vcf, bed = write_files(tmp_path, bed_regions)
    out = io.StringIO()
    gnomad_filter(vcf, bed, SimpleNamespace(ignore_flt=True), out=out, **kw)
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
