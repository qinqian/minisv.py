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
        dbg=False,
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


from minisv.filtercaller import _consensus_lost_by_filter


# Two consensus reps. far_DEL is at chr1:1000; near_DEL is at chr2:5000.
# The "target" union file will contain only a chr2:5000 row — so near_DEL has an
# overlap in the target (and is "lost"-by-the-target, i.e. excluded from the
# helper's output) while far_DEL does not (and is emitted).
CONSENSUS_VCF = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1000\tfar_DEL\tATCGATCGAT\tA\t60\tPASS\tSVTYPE=DEL;SVLEN=-9",
    "chr2\t5000\tnear_DEL\tATCGATCGAT\tA\t60\tPASS\tSVTYPE=DEL;SVLEN=-9",
]
TARGET_VCF = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr2\t5000\ttarget_hit\tATCGATCGAT\tA\t60\tPASS\tSVTYPE=DEL;SVLEN=-9",
]


def _build_msv(tmp_path, name, vcf_lines):
    """Produce a real .msv file by running union_sv on a single synthetic VCF."""
    vcf = tmp_path / f"{name}.vcf"
    vcf.write_text("\n".join(vcf_lines) + "\n")
    msv = tmp_path / f"{name}.msv"
    with open(msv, "w") as fh:
        union_sv([str(vcf)], 0, _opt(print_sv=True), file_handler=fh)
    return str(msv)


def test_consensus_lost_by_filter_excludes_overlapping_emits_disjoint(tmp_path):
    consensus = _build_msv(tmp_path, "consensus", CONSENSUS_VCF)
    target = _build_msv(tmp_path, "target", TARGET_VCF)

    out = io.StringIO()
    _consensus_lost_by_filter(consensus, target, _opt(), file_handler=out)
    text = out.getvalue()

    # far_DEL (chr1:1000 → 0-based 999) has nothing in target → emitted.
    # near_DEL (chr2:5000 → 0-based 4999) overlaps target_hit → excluded.
    lines = [ln for ln in text.splitlines() if ln.strip()]
    assert len(lines) == 1, f"expected exactly one emitted line, got: {text!r}"
    assert lines[0].startswith("chr1\t999\t")


def test_consensus_lost_by_filter_empty_target_emits_all(tmp_path):
    consensus = _build_msv(tmp_path, "consensus2", CONSENSUS_VCF)
    target = tmp_path / "empty.msv"
    target.write_text("")

    out = io.StringIO()
    _consensus_lost_by_filter(consensus, str(target), _opt(), file_handler=out)
    lines = [ln for ln in out.getvalue().splitlines() if ln.strip()]
    assert len(lines) == 2
