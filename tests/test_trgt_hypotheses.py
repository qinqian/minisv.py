from minisv.trio_alleles import allele_hypotheses

def test_unknown_origin_keeps_parent_hypotheses_and_30_boundary():
    row=dict(denovo_status='Y:?',index='0',child_AL='130',allele_origin='?',father_AL='100,101',mother_AL='100,130')
    result=allele_hypotheses(row,minimum=30)
    assert [(r['parent'],r['delta']) for r in result]==[('F',30),('M',30)]
    assert allele_hypotheses(row)==[]
