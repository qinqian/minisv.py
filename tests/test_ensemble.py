"""Unit tests for insilico_truth (ensemble.py)."""
import io

from minisv.ensemble import insilico_truth


# Two groups. g1 has two rows (median = the second after sort-by-(start,end)),
# g2 has a single row. The buggy code emits only g1 (flushed on the g1->g2
# transition) and drops g2 because nothing flushes after the loop ends.
MSV_LINES = [
    "chr1\t1000\t>>\tchr1\t1200\tgroup_id=g1;file_id=0",
    "chr1\t1010\t>>\tchr1\t1210\tgroup_id=g1;file_id=1",
    "chr2\t5000\t>>\tchr2\t5200\tgroup_id=g2;file_id=0",
]


def test_insilico_truth_emits_last_group(tmp_path):
    msv = tmp_path / "union.msv"
    msv.write_text("\n".join(MSV_LINES) + "\n")

    out = io.StringIO()
    insilico_truth(str(msv), file_handler=out)
    text = out.getvalue()
    lines = [ln for ln in text.splitlines() if ln.strip()]

    # g1 median: sorted by (start,end) -> [(1000,1200),(1010,1210)];
    # index len>>1 == 1 -> the chr1:1010 row.
    assert any(ln.startswith("chr1\t1010\t") for ln in lines), text
    # g2 (the final group) must be emitted too.
    assert any(ln.startswith("chr2\t5000\t") for ln in lines), text
    # exactly one median line per group
    assert len(lines) == 2, text
