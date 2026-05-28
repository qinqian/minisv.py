"""Unit tests for the cross-caller consensus export feature."""
import io
from types import SimpleNamespace

from minisv.union import union_sv


# Three single-row caller VCFs. The DEL at chr1:1000 is shared by all three callers
# (the consensus). The DEL at chr1:8000 appears in only caller A.
SHARED_VCF = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1000\tshared\tATCGATCGAT\tA\t60\tPASS\tSVTYPE=DEL;SVLEN=-9",
]
UNIQUE_VCF_A = SHARED_VCF + [
    "chr1\t8000\tonly_A\tATCGATCGAT\tA\t60\tPASS\tSVTYPE=DEL;SVLEN=-9",
]


def _opt(print_sv=True):
    # Fields union_sv reads. Zero thresholds + huge win_size so synthetic rows survive
    # all the in-loop filters except the new min_file_count one.
    return SimpleNamespace(
        bed=None,
        read_min_count=0,
        win_size=500,
        min_len_ratio=0.5,
        group_min_count=0,
        min_len=0,
        print_sv=print_sv,
    )


def _write_vcfs(tmp_path):
    paths = []
    # caller A has both rows; B and C only carry the shared row
    for name, lines in [("a", UNIQUE_VCF_A), ("b", SHARED_VCF), ("c", SHARED_VCF)]:
        p = tmp_path / f"{name}.vcf"
        p.write_text("\n".join(lines) + "\n")
        paths.append(str(p))
    return paths


def test_union_sv_min_file_count_3_keeps_only_all_three_groups(tmp_path):
    vcfs = _write_vcfs(tmp_path)
    out = io.StringIO()
    union_sv(vcfs, 0, _opt(), file_handler=out, min_file_count=3)
    text = out.getvalue()
    # The shared DEL at chr1:1000 (pos=999 in 0-based MSV output) appears in all 3 callers.
    # The only_A DEL at chr1:8000 (pos=7999) has file-mask 0b001 (bit-count 1), so
    # min_file_count=3 drops it — its position must not appear in output.
    assert "chr1\t999\t" in text
    assert "chr1\t7999\t" not in text


def test_union_sv_min_file_count_none_is_unchanged(tmp_path):
    vcfs = _write_vcfs(tmp_path)
    out = io.StringIO()
    union_sv(vcfs, 0, _opt(), file_handler=out)  # default min_file_count=None
    text = out.getvalue()
    # both groups print under the default behavior
    assert "chr1\t999\t" in text
    assert "chr1\t7999\t" in text
