from minisv.trio_alleles import allele_hypotheses

def test_non_denovo_input_does_not_generate_hypotheses():
    assert allele_hypotheses({'denovo_status':'N'})==[]
