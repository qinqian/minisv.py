"""Merge prepared caller records and allele-specific TRGT hypotheses."""
import gzip
import json
from pathlib import Path

from .trio_alleles import build_candidates


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def prepare_candidates(child, config, output):
    minimum = config.get('minimum_variant_length', 50)
    paths = list(child.get('candidate_vcfs', []))
    if child.get('candidate_vcf'):
        paths.insert(0, child['candidate_vcf'])
    if child.get('trgt_signal'):
        hypotheses = output / 'trgt_allele_hypotheses.vcf'
        spanning = Path(config['trgt_directory']) / f'{child["sample"]}.spanning.sorted.bam'
        build_candidates(child['trgt_signal'], str(spanning), str(hypotheses),
                         str(output / 'allele_audit.jsonl'), minimum=minimum, indexed=True)
        paths.append(str(hypotheses))
    if not paths:
        raise ValueError('No candidate input files')
    headers = {}
    for path in paths:
        with open_text(path) as source:
            for line in source:
                if not line.startswith('#'):
                    break
                if line.startswith(('##INFO=', '##FILTER=', '##contig=')):
                    key = (line.split('=', 1)[0], line.split('ID=', 1)[1].split(',', 1)[0])
                    headers.setdefault(key, line.rstrip() + '\n')
    combined = output / 'all_candidates.vcf'
    counts = []
    with combined.open('w') as target:
        target.write('##fileformat=VCFv4.2\n')
        target.writelines(headers.values())
        target.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for index, path in enumerate(paths):
            count = 0
            with open_text(path) as source:
                for line in source:
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.rstrip().split('\t')
                    info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
                    names = set(info.get('RNAMES', '').split(',')) - {'', '.'}
                    if fields[6] not in ('PASS', '.') or abs(int(info.get('SVLEN', 0))) < minimum or len(names) < 2:
                        continue
                    fields = fields[:8]
                    fields[2] = f'SOURCE{index}.{count + 1}.{fields[2]}'
                    target.write('\t'.join(fields) + '\n')
                    count += 1
            counts.append(dict(path=str(path), records=count))
    (output / 'candidate_summary.json').write_text(json.dumps(dict(inputs=counts,
        total_records=sum(x['records'] for x in counts), minimum_variant_length=minimum), indent=2))
    return str(combined)
